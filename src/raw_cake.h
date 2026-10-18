// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// Bespoke joint solver for the "V3 Raw Cake" -- the V3 "Tone Bypass"
// board (a 1979 unit): a late-1970s V3 board with all FOUR transistors
// MPSA18 (very high-gain, low-noise) and a TONE BYPASS DPDT switch (SW1, like the
// V5/V6). The switch routes Q2's collector either THROUGH the scooping tone stack
// or AROUND it via R28 (150k) -- flat, louder, "raw". Exposed as a live toggle
// via set_bypass() (wired to the plugin's Tone-Bypass port).
//
// The switch is modelled as low/high-resistance contacts (1 ohm closed / 100 Meg
// open) so the sparse pattern is FIXED: set_bypass() only flips four contact
// conductances. The R28 bypass path and the tone stack are BOTH always present;
// whichever is unselected sits behind open (100 Meg) contacts and is harmless.
//
// Architecture (see spice/the_raw_cake.cir):
//   in -R2- C1 -> Q4 -C4-> Sustain -C5- R19 -> Q3 (clip D3/D4, Re=100)
//      -C13- R12 -> Q2 (clip D1/D2, Re=100) -C15-> [SW1A] -> {tone stack | R28}
//      -> [SW1B] -C16- C3 -> Q1 (output) -C2-> Volume -> out
//
// All four transistors MPSA18. Full Ebers-Moll (Br=1). 31-unknown joint Newton
// per sample; ANALYTIC Jacobian + fixed-order sparse LU (gen::solve_raw_cake).
//
// Unknowns x[0..30]:
//   0 n_in  1 q4b  2 q4c  3 q4e
//   4 sus_top  5 sus_w  6 sus_bot  7 c5n
//   8 q3b  9 q3c  10 q3aux  11 q3e
//   12 c13n
//   13 q2b  14 q2c  15 q2aux  16 q2e
//   17 swa  18 tone_in  19 tt  20 tw  21 tb
//   22 r28a  23 r28b  24 swb  25 c16c3
//   26 q1b  27 q1c  28 q1e
//   29 vol_top  30 out

#pragma once

#include "dk_circuits.h"
#include "engine_select.h"

#include "newton_arbiter.h"

#include "cake_variant.h"
#include "bjt_dc_bias.h"
#include "nr_solver.h"
#include "solvers_generated.h"   // fixed-order sparse LU (gen::solve_raw_cake)

#include <algorithm>
#include <array>
#include <cmath>
#include "finite_guard.h"

namespace tocake {

class RawCakeNode {
public:
    using T = double;
    static constexpr std::size_t N = 31;
    enum class Tap { C2Out, Out };

    explicit RawCakeNode(const CakeVariant& v,
                         double p_tone = 0.5,
                         double p_vol  = 0.5,
                         double p_sus  = 1.0) noexcept
        : v_(v), p_tone_(p_tone), p_vol_(p_vol), p_sus_(p_sus) {}

    void prepare(double sample_rate) {
        // Rest seed from the engine's OWN DC: at a vanishing rate every companion conductance is zero -- the caps
        // are open, i.e. DC -- so this engine's own Newton, fed silence, lands on its operating point with exactly
        // the equations it runs. The constants in reset() are only the starting point, so a device-law change
        // cannot stale them.
        setup_rate(kDcSeedRate);
        reset();
        std::array<double, N> x_prev = x_;
        bool settled = false;
        for (int k = 0; k < kDcSeedMaxPasses && !settled; ++k) {
            (void) run_sample(0.0f, Tap::Out);
            double step = 0.0;
            for (std::size_t n = 0; n < N; ++n) step = std::max(step, std::fabs(x_[n] - x_prev[n]));
            x_prev = x_;
            settled = tocake::all_finite(x_) && step < kDcSeedTol;
        }
        if (settled) seed_ = x_;
        dc_seed_settled_ = settled;
        setup_rate(sample_rate);
        reset();
    }
    bool dc_seed_settled() const noexcept { return dc_seed_settled_; }

    void setup_rate(double sample_rate) {
        fs_ = sample_rate;
        const auto& ib   = v_.input_booster;
        const auto& sus  = v_.sustain;
        const auto& cl1  = v_.clipper_1;
        const auto& cl2  = v_.clipper_2;
        const auto& tone = v_.tone;
        const auto& ob   = v_.output_booster;

        G_C1_  = 2.0 * ib.C_input_couple   * fs_;
        G_C10_ = 2.0 * ib.C_feedback_hf    * fs_;
        G_C4_  = 2.0 * ib.C_output_couple  * fs_;
        G_C5_  = 2.0 * sus.C_wiper_couple  * fs_;
        G_C12_ = 2.0 * cl1.C_feedback_hf   * fs_;
        G_C6_  = 2.0 * cl1.C_diode_couple  * fs_;
        G_C13_ = 2.0 * cl1.C_output_couple * fs_;
        G_C11_ = 2.0 * cl2.C_feedback_hf   * fs_;
        G_C7_  = 2.0 * cl2.C_diode_couple  * fs_;
        G_C9_  = 2.0 * tone.C_hp_series    * fs_;
        G_C8_  = 2.0 * tone.C_lp_shunt     * fs_;
        G_C3_  = 2.0 * ob.C_input_couple   * fs_;
        G_C2_  = 2.0 * ob.C_output_couple  * fs_;
        G_C15_ = 2.0 * sp_C15_ * fs_;
        G_C16_ = 2.0 * sp_C16_ * fs_;

        G_R2_  = 1.0 / ib.R_input_series;
        G_R14_ = 1.0 / ib.R_base_to_gnd;
        G_R9_  = 1.0 / ib.R_feedback;
        G_R13_ = 1.0 / ib.R_collector_load;
        G_R22_ = 1.0 / ib.R_emitter;
        G_R19_ = 1.0 / cl1.R_input_series;
        G_R23_ = 1.0 / sus.R_pot_bottom;
        G_R20_ = 1.0 / cl1.R_base_to_gnd;
        G_R17_ = 1.0 / cl1.R_feedback;
        G_R18_ = 1.0 / cl1.R_collector_load;
        G_R21_ = 1.0 / cl1.R_emitter;
        G_R12_ = 1.0 / cl2.R_input_series;
        G_R16_ = 1.0 / cl2.R_base_to_gnd;
        G_R15_ = 1.0 / cl2.R_feedback;
        G_R11_ = 1.0 / cl2.R_collector_load;
        G_R10_ = 1.0 / cl2.R_emitter;
        G_R5_  = 1.0 / tone.R_hp_shunt;
        G_R8_  = 1.0 / tone.R_lp_series;
        G_R28_ = 1.0 / sp_R28_;
        G_R7_  = 1.0 / ob.R_bias_top;
        G_R3_  = 1.0 / ob.R_base_to_gnd;
        G_R6_  = 1.0 / ob.R_collector_load;
        G_R4_  = 1.0 / ob.R_emitter;

        twoIs_d_ = 2.0 * cl1.diodes.Is;  Vt_d_ = cl1.diodes.N * V_T_300K;

        recompute_sustain_g();
        recompute_tone_g();
        recompute_volume_g();
        set_bypass(bypass_);
        reset();
    }

    void set_sustain(double frac) noexcept { p_sus_ = std::clamp(frac,0.0,1.0); recompute_sustain_g(); }
    void set_tone(double frac)    noexcept { p_tone_= std::clamp(frac,0.0,1.0); recompute_tone_g(); }
    void set_volume(double frac)  noexcept { p_vol_ = std::clamp(frac,0.0,1.0); recompute_volume_g(); }

    // Tone-bypass switch: false = tone stack in circuit, true = flat R28 bypass.
    void set_bypass(bool b) noexcept {
        bypass_ = b;
        const double cl = 1.0, op = 1e-8;        // 1 ohm closed / 100 Meg open
        G_swA_t_ = b ? op : cl;  G_swA_b_ = b ? cl : op;
        G_swB_t_ = b ? op : cl;  G_swB_b_ = b ? cl : op;
    }

    void reset() noexcept {
        // Warm-start = the settled DC (model's fixed point; c16c3 is a cap-only
        // node whose true DC is 0, pinned by the 100M leak like SPICE's GMIN).
        x_ = seed_;   // the rest point prepare() polished (the hand-kept constants until it runs)
        init_cap_states();
        dbg_iters_ = dbg_calls_ = dbg_nonconv_ = 0;
        dbg_max_resid_ = 0.0;
    }

    float process_sample(float vin_volts) noexcept { return run_sample(vin_volts, Tap::Out); }
    float process_sample(float vin_volts, Tap tap) noexcept { return run_sample(vin_volts, tap); }

    const std::array<double,N>& state() const noexcept { return x_; }
    long   dbg_iters()     const noexcept { return dbg_iters_;   }
    long   dbg_calls()     const noexcept { return dbg_calls_;   }
    long   dbg_nonconv()   const noexcept { return dbg_nonconv_; }
    double dbg_max_resid() const noexcept { return dbg_max_resid_; }
    double dbg_max_jac_rel() const noexcept { return dbg_max_jac_rel_; }

private:
    void recompute_sustain_g() noexcept {
        const double R = v_.sustain.R_pot;
        G_Rst_ = 1.0 / (R * (1.0 - p_sus_) + 1e-6);
        G_Rsb_ = 1.0 / (R * p_sus_         + 1e-6);
    }
    void recompute_tone_g() noexcept {
        const double R = v_.tone.R_pot;
        G_Rtt_ = 1.0 / (R * (1.0 - p_tone_) + 1e-6);
        G_Rtb_ = 1.0 / (R * p_tone_         + 1e-6);
    }
    void recompute_volume_g() noexcept {
        const double R = v_.volume.R_pot;
        G_Rvt_ = 1.0 / (R * (1.0 - p_vol_) + 1e-6);
        G_Rvb_ = 1.0 / (R * p_vol_         + 1e-6);
    }

    void init_cap_states() noexcept {
        s_C1_  = G_C1_  * (x_[0]  - x_[1]);
        s_C10_ = G_C10_ * (x_[1]  - x_[2]);
        s_C4_  = G_C4_  * (x_[2]  - x_[4]);
        s_C5_  = G_C5_  * (x_[5]  - x_[7]);
        s_C12_ = G_C12_ * (x_[8]  - x_[9]);
        s_C6_  = G_C6_  * (x_[8]  - x_[10]);
        s_C13_ = G_C13_ * (x_[9]  - x_[12]);
        s_C11_ = G_C11_ * (x_[13] - x_[14]);
        s_C7_  = G_C7_  * (x_[13] - x_[15]);
        s_C15_ = G_C15_ * (x_[14] - x_[17]);
        s_C9_  = G_C9_  * (x_[18] - x_[19]);
        s_C8_  = G_C8_  * (x_[21] - 0.0);
        s_C16_ = G_C16_ * (x_[24] - x_[25]);
        s_C3_  = G_C3_  * (x_[25] - x_[26]);
        s_C2_  = G_C2_  * (x_[27] - x_[29]);
    }

    void residual(const std::array<double,N>& x, std::array<double,N>& out) const noexcept {
        const double V_T = V_T_300K;
        const double vcc = v_.vcc;
        out.fill(0.0);

        auto bjt = [V_T](const bjt_params& q, double Vb, double Vc, double Ve,
                         double& IC, double& IB, double& IE) noexcept {
            constexpr double Br = 1.0;
            const double ef    = std::exp((Vb - Ve) / V_T);
            const double er    = std::exp((Vb - Vc) / V_T);
            const double early = 1.0 + (Vc - Vb) / q.Vaf;
            IC = q.Is * (ef - er) * early - (q.Is / Br) * (er - 1.0);
            IB = (q.Is / q.Bf) * (ef - 1.0) + (q.Is / Br) * (er - 1.0);
            IE = IC + IB;
        };
        auto cap = [&](double G, int p, int n, double s) noexcept {
            const double i = G * (x[p] - x[n]) - s; out[p] += i; out[n] -= i; };
        auto capg = [&](double G, int p, double s) noexcept { out[p] += G * x[p] - s; };
        auto res = [&](double G, int p, int n) noexcept {
            const double i = G * (x[p] - x[n]); out[p] += i; out[n] -= i; };
        auto resg = [&](double G, int p) noexcept { out[p] += G * x[p]; };
        auto resv = [&](double G, int p, double Vf) noexcept { out[p] += G*(x[p]-Vf); };
        auto diode = [&](int p, int n) noexcept {
            const double i = twoIs_d_ * std::sinh((x[p] - x[n]) / Vt_d_);
            out[p] += i; out[n] -= i; };

        const auto& q4 = v_.input_booster.q;
        const auto& q3 = v_.clipper_1.q;
        const auto& q2 = v_.clipper_2.q;
        const auto& q1 = v_.output_booster.q;

        // --- Resistors ---
        resv(G_R2_, 0, vin_);                  // R2: n_in - in
        resg(G_R14_, 1); res(G_R9_, 1, 2); resv(G_R13_, 2, vcc); resg(G_R22_, 3);
        res(G_R19_, 7, 8); resg(G_R23_, 6); resg(G_R20_, 8); res(G_R17_, 8, 9);
        resv(G_R18_, 9, vcc); resg(G_R21_, 11); res(G_R12_, 12, 13); resg(G_R16_, 13);
        res(G_R15_, 13, 14); resv(G_R11_, 14, vcc); resg(G_R10_, 16);
        resg(G_R5_, 19); res(G_R8_, 18, 21); res(G_R28_, 22, 23);
        resv(G_R7_, 26, vcc); resg(G_R3_, 26); resv(G_R6_, 27, vcc); resg(G_R4_, 28);
        // Switch contacts (SW1A: swa-tone_in / swa-r28a ; SW1B: swb-tw / swb-r28b)
        res(G_swA_t_, 17, 18); res(G_swA_b_, 17, 22);
        res(G_swB_t_, 24, 20); res(G_swB_b_, 24, 23);
        // Pots
        res(G_Rst_, 4, 5);   res(G_Rsb_, 5, 6);
        res(G_Rtt_, 19, 20); res(G_Rtb_, 20, 21);
        res(G_Rvt_, 29, 30); resg(G_Rvb_, 30);
        resg(1e-8, 25);   // 100M leak: pin the C16-C3 junction (cap-only node) like SPICE GMIN

        // --- Caps ---
        cap(G_C1_, 0, 1, s_C1_);    cap(G_C10_, 1, 2, s_C10_);  cap(G_C4_, 2, 4, s_C4_);
        cap(G_C5_, 5, 7, s_C5_);    cap(G_C12_, 8, 9, s_C12_);  cap(G_C6_, 8, 10, s_C6_);
        cap(G_C13_, 9, 12, s_C13_); cap(G_C11_, 13, 14, s_C11_);cap(G_C7_, 13, 15, s_C7_);
        cap(G_C15_, 14, 17, s_C15_);cap(G_C9_, 18, 19, s_C9_);  capg(G_C8_, 21, s_C8_);
        cap(G_C16_, 24, 25, s_C16_);cap(G_C3_, 25, 26, s_C3_);  cap(G_C2_, 27, 29, s_C2_);

        // --- Diode pairs ---
        diode(10, 9);    // D3/D4
        diode(15, 14);   // D1/D2

        // --- Transistors (full EM) ---
        double IC,IB,IE;
        bjt(q4, x[1], x[2], x[3], IC,IB,IE);    out[1]+=IB; out[2]+=IC; out[3]-=IE;
        bjt(q3, x[8], x[9], x[11], IC,IB,IE);   out[8]+=IB; out[9]+=IC; out[11]-=IE;
        bjt(q2, x[13],x[14],x[16], IC,IB,IE);   out[13]+=IB;out[14]+=IC;out[16]-=IE;
        bjt(q1, x[26],x[27],x[28], IC,IB,IE);   out[26]+=IB;out[27]+=IC;out[28]-=IE;
    }

    void jacobian(const std::array<double,N>& x,
                  std::array<std::array<double,N>,N>& out) const noexcept {
        const double V_T = V_T_300K;
        for (auto& row : out) row.fill(0.0);

        auto jr = [&](double G, int p, int n) noexcept {
            out[p][p]+=G; out[p][n]-=G; out[n][p]-=G; out[n][n]+=G; };
        auto jg = [&](double G, int p) noexcept { out[p][p]+=G; };
        auto jd = [&](double g, int p, int n) noexcept {
            out[p][p]+=g; out[p][n]-=g; out[n][p]-=g; out[n][n]+=g; };

        auto bjt = [V_T](const bjt_params& q, double Vb, double Vc, double Ve,
                         double& gBB,double& gBC,double& gBE,
                         double& ibB,double& ibC,double& ibE) noexcept {
            constexpr double Br = 1.0;
            const double Is = q.Is, Vaf = q.Vaf;
            const double ef = std::exp((Vb-Ve)/V_T), er = std::exp((Vb-Vc)/V_T);
            const double early = 1.0 + (Vc-Vb)/Vaf;
            gBB = Is*(ef-er)/V_T*early - Is*(ef-er)/Vaf - (Is/Br)*er/V_T;
            gBC = Is*(er/V_T)*early + Is*(ef-er)/Vaf + (Is/Br)*er/V_T;
            gBE = -Is*ef/V_T*early;
            ibB = (Is/q.Bf)*ef/V_T + (Is/Br)*er/V_T;
            ibC = -(Is/Br)*er/V_T;
            ibE = -(Is/q.Bf)*ef/V_T;
        };

        // --- Resistors ---
        jg(G_R2_,0);
        jg(G_R14_,1); jr(G_R9_,1,2); jg(G_R13_,2); jg(G_R22_,3);
        jr(G_R19_,7,8); jg(G_R23_,6); jg(G_R20_,8); jr(G_R17_,8,9); jg(G_R18_,9);
        jg(G_R21_,11); jr(G_R12_,12,13); jg(G_R16_,13); jr(G_R15_,13,14); jg(G_R11_,14);
        jg(G_R10_,16); jg(G_R5_,19); jr(G_R8_,18,21); jr(G_R28_,22,23);
        jg(G_R7_,26); jg(G_R3_,26); jg(G_R6_,27); jg(G_R4_,28);
        // Switch contacts
        jr(G_swA_t_,17,18); jr(G_swA_b_,17,22); jr(G_swB_t_,24,20); jr(G_swB_b_,24,23);
        // Pots
        jr(G_Rst_,4,5); jr(G_Rsb_,5,6); jr(G_Rtt_,19,20); jr(G_Rtb_,20,21);
        jr(G_Rvt_,29,30); jg(G_Rvb_,30);
        jg(1e-8, 25);   // C16-C3 junction leak (matches residual)

        // --- Caps ---
        jr(G_C1_,0,1); jr(G_C10_,1,2); jr(G_C4_,2,4); jr(G_C5_,5,7); jr(G_C12_,8,9);
        jr(G_C6_,8,10); jr(G_C13_,9,12); jr(G_C11_,13,14); jr(G_C7_,13,15);
        jr(G_C15_,14,17); jr(G_C9_,18,19); jg(G_C8_,21); jr(G_C16_,24,25);
        jr(G_C3_,25,26); jr(G_C2_,27,29);

        // --- Diode pairs ---
        jd(twoIs_d_ * std::cosh((x[10]-x[9])  / Vt_d_) / Vt_d_, 10, 9);
        jd(twoIs_d_ * std::cosh((x[15]-x[14]) / Vt_d_) / Vt_d_, 15, 14);

        // --- Transistors ---
        const auto& q4 = v_.input_booster.q;
        const auto& q3 = v_.clipper_1.q;
        const auto& q2 = v_.clipper_2.q;
        const auto& q1 = v_.output_booster.q;
        double gBB,gBC,gBE,ibB,ibC,ibE;
        auto stampQ = [&](int b,int c,int e){
            out[b][b]+=ibB; out[b][c]+=ibC; out[b][e]+=ibE;
            out[c][b]+=gBB; out[c][c]+=gBC; out[c][e]+=gBE;
            out[e][b]-=gBB+ibB; out[e][c]-=gBC+ibC; out[e][e]-=gBE+ibE; };
        bjt(q4, x[1],x[2],x[3],   gBB,gBC,gBE, ibB,ibC,ibE); stampQ(1,2,3);
        bjt(q3, x[8],x[9],x[11],  gBB,gBC,gBE, ibB,ibC,ibE); stampQ(8,9,11);
        bjt(q2, x[13],x[14],x[16],gBB,gBC,gBE, ibB,ibC,ibE); stampQ(13,14,16);
        bjt(q1, x[26],x[27],x[28],gBB,gBC,gBE, ibB,ibC,ibE); stampQ(26,27,28);
    }

    float run_sample(float vin_volts, Tap tap) noexcept {
        vin_ = static_cast<double>(vin_volts);

        std::array<double,N> Fv{};
        std::array<std::array<double,N>,N> J{};
        ++dbg_calls_;
        bool conv=false; double lastF=0.0; int it=0;
        for (; it < 32; ++it) {
            residual(x_, Fv);
            double maxF=0.0; for (std::size_t i=0;i<N;++i) maxF=std::max(maxF,std::abs(Fv[i]));
            lastF=maxF;
            if (maxF < 1e-7 && it >= 1) { conv=true; break; }
            jacobian(x_, J);
            std::array<double,N> b{};
            for (std::size_t i=0;i<N;++i) b[i] = -Fv[i];
            if (!gen::solve_raw_cake(J, b)) break;
            for (std::size_t i=0;i<N;++i) {
                double d=b[i];
                if (d> 0.15) d= 0.15;
                if (d<-0.15) d=-0.15;
                x_[i]+=d;
            }
        }
        dbg_iters_ += it;
        if (!conv) ++dbg_nonconv_;
        if (lastF > dbg_max_resid_) dbg_max_resid_ = lastF;
        // Every unknown, and recover to the DC rest point: one index cannot
        // stand in for the vector, and rewriting the unknowns alone leaves the
        // capacitor states poisoned (see src/finite_guard.h).
        if (!tocake::all_finite(x_)) { reset(); ++dbg_nonfinite_; }

        s_C1_  = 2.0*G_C1_ *(x_[0] -x_[1])  - s_C1_;
        s_C10_ = 2.0*G_C10_*(x_[1] -x_[2])  - s_C10_;
        s_C4_  = 2.0*G_C4_ *(x_[2] -x_[4])  - s_C4_;
        s_C5_  = 2.0*G_C5_ *(x_[5] -x_[7])  - s_C5_;
        s_C12_ = 2.0*G_C12_*(x_[8] -x_[9])  - s_C12_;
        s_C6_  = 2.0*G_C6_ *(x_[8] -x_[10]) - s_C6_;
        s_C13_ = 2.0*G_C13_*(x_[9] -x_[12]) - s_C13_;
        s_C11_ = 2.0*G_C11_*(x_[13]-x_[14]) - s_C11_;
        s_C7_  = 2.0*G_C7_ *(x_[13]-x_[15]) - s_C7_;
        s_C15_ = 2.0*G_C15_*(x_[14]-x_[17]) - s_C15_;
        s_C9_  = 2.0*G_C9_ *(x_[18]-x_[19]) - s_C9_;
        s_C8_  = 2.0*G_C8_ *(x_[21]-0.0)    - s_C8_;
        s_C16_ = 2.0*G_C16_*(x_[24]-x_[25]) - s_C16_;
        s_C3_  = 2.0*G_C3_ *(x_[25]-x_[26]) - s_C3_;
        s_C2_  = 2.0*G_C2_ *(x_[27]-x_[29]) - s_C2_;

        switch (tap) {
            case Tap::C2Out: return static_cast<float>(x_[14] - 4.50);
            case Tap::Out:   return static_cast<float>(x_[30]);
        }
        return static_cast<float>(x_[30]);
    }

    const CakeVariant& v_;
    double p_tone_, p_vol_, p_sus_;
    double fs_ = 48000.0;
    double vin_ = 0.0;
    bool   bypass_ = false;

    // Tone-bypass extras (not in CakeVariant), from the schematic.
    double sp_C15_ = 0.1e-6, sp_C16_ = 0.1e-6, sp_R28_ = 150e3;

    double G_C1_=0,G_C10_=0,G_C4_=0,G_C5_=0,G_C12_=0,G_C6_=0,G_C13_=0,G_C11_=0,G_C7_=0;
    double G_C15_=0,G_C9_=0,G_C8_=0,G_C16_=0,G_C3_=0,G_C2_=0;
    double G_R2_=0,G_R14_=0,G_R9_=0,G_R13_=0,G_R22_=0,G_R19_=0,G_R23_=0,G_R20_=0,G_R17_=0,G_R18_=0,G_R21_=0;
    double G_R12_=0,G_R16_=0,G_R15_=0,G_R11_=0,G_R10_=0,G_R5_=0,G_R8_=0,G_R28_=0;
    double G_R7_=0,G_R3_=0,G_R6_=0,G_R4_=0;
    double G_Rst_=0,G_Rsb_=0,G_Rtt_=0,G_Rtb_=0,G_Rvt_=0,G_Rvb_=0;
    double G_swA_t_=1.0,G_swA_b_=1e-8,G_swB_t_=1.0,G_swB_b_=1e-8;
    double twoIs_d_=0,Vt_d_=0;

    double s_C1_=0,s_C10_=0,s_C4_=0,s_C5_=0,s_C12_=0,s_C6_=0,s_C13_=0,s_C11_=0,s_C7_=0;
    double s_C15_=0,s_C9_=0,s_C8_=0,s_C16_=0,s_C3_=0,s_C2_=0;

    long dbg_iters_=0, dbg_calls_=0, dbg_nonconv_=0;
    double dbg_max_resid_=0.0;
    double dbg_max_jac_rel_=0.0;

    std::array<double,N> x_{};
    std::array<double,N> seed_ {
            8.39222242e-13,     0.649113533, 3.85196336, 0.0336711307,   // n_in q4b q4c q4e
            2.71753746e-11,     2.71753746e-11,     2.69063115e-13,     0.649113507,    // sus_top sus_w sus_bot c5n
            0.649113507, 3.8519655, 3.8519655, 0.0336711165,    // q3b q3c q3aux(=q3c, diode-balance) q3e
            0.649113503,                               // c13n
            0.649113503, 3.85196582, 3.85196582, 0.0336711144,    // q2b q2c q2aux(=q2c) q2e
            -2.81730011e-09,     -2.81727406e-09,     -1.52283198e-10,     -4.98464286e-10,  -1.80122791e-09,  // swa tone_in tt tw tb
            -1.65874155e-09,     -1.65700372e-09,                          // r28a r28b
            -4.98445155e-10,     -3.83043992e-06,                          // swb c16c3
            1.79550409, 3.69664023, 1.16787881,             // q1b q1c q1e
            3.3128908e-11,     1.6564454e-11                            // vol_top out
        };   // hand-kept start for the DC polish in prepare()
    static constexpr double kDcSeedRate = 1.0e-9;   // Hz: companion conductances ~0 (as coloured_goat.h)
    static constexpr int    kDcSeedMaxPasses = 400;
    static constexpr double kDcSeedTol = 1.0e-6;    // V: above the Newton's own floor
    bool dc_seed_settled_ = false;
    // Times the guard fired. Exposed so a harness can SEE it: a guard
    // nothing can read is an alarm with no wire.
    long dbg_nonfinite_ = 0;
public:
    long nonfinite_resets() const noexcept { return dbg_nonfinite_; }
private:
};

// ---- the engine that ships -------------------------------------------------
// RawCake is the DK description of this circuit (src/dk_circuits.h);
// RawCakeNode above is the reference it is checked against, forced to converge.
using RawCake = shipped_engine_t<DkRawCake>;   // the DK, or its nodal twin (engine_select.h)

} // namespace tocake
