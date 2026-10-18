// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// Bespoke joint solver for the "V2 Goat&S OD" -- a Japanese clone of
// the circuit (1st edition, ~1974). Essentially the four-stage fuzz (4 transistors,
// 2 diode-clipper stages) but with a MODIFIED TONE STACK: an active RC
// network wrapped around the output transistor Q1 that fixes the bass value
// across the tone-pot sweep. The standard skeleton (passive tone stack) cannot
// represent it, so -- like the goats/PM-3/op-amp -- it gets its own solver.
//
// Architecture (see spice/the_goat_s_od.cir, cross-checked vs the published
// schematic):
//   in -C1- R2 -> Q4 (booster) -C4-> Sustain -C5- R19 -> Q3 (clip D3/D4, Re=180)
//      -C13- R12 -> Q2 (clip D1/D2, Re=180) -C14-> tone-top (ta)
//   Tone pot ta..tb, wiper -C9- R8 -> y -C3-> Q1 base; Q1 (Re=1k || C17 10u)
//   collector -C2-> tb. A ladder R5/C15/R28 - R27/C16/R29 - R30 couples ta..tb
//   into y. Output: tb -R31- volume -> out.
//
// All four transistors have emitter nodes (Q3/Q2 Re=180, not grounded). Full
// Ebers-Moll (Br=1). 30-unknown joint Newton per sample; ANALYTIC Jacobian +
// fixed-order sparse LU (gen::solve_goat_s_od), same optimization as the goats.
//
// Transistor C828 (Si NPN): provisional Is=1.5e-14 Bf=270 Vaf=80 (R-grade typ;
// no grade marked on the schematic -- refine Bf from the datasheet if desired).
//
// Unknowns x[0..29]:
//   0 n_c1   1 q4b   2 q4c   3 q4e
//   4 sus_top  5 sus_w  6 sus_bot   7 c5n
//   8 q3b   9 q3c   10 q3aux  11 q3e
//   12 c13n
//   13 q2b  14 q2c  15 q2aux  16 q2e
//   17 ta   18 tw   19 c9n   20 y   21 tb
//   22 x1   23 x2   24 x3
//   25 q1b  26 q1c  27 q1e
//   28 vol_top  29 out

#pragma once

#include "dk_circuits.h"
#include "engine_select.h"

#include "newton_arbiter.h"

#include "cake_variant.h"
#include "bjt_dc_bias.h"
#include "nr_solver.h"
#include "solvers_generated.h"   // fixed-order sparse LU (gen::solve_goat_s_od)

#include <algorithm>
#include <array>
#include <cmath>
#include "finite_guard.h"

namespace tocake {

class GoatSODNode {
public:
    using T = double;
    static constexpr std::size_t N = 30;
    enum class Tap { C2Out, Out };

    explicit GoatSODNode(const CakeVariant& v,
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
        const auto& ob   = v_.output_booster;

        // Cap conductances (G = 2*C*fs).
        G_C1_  = 2.0 * ib.C_input_couple   * fs_;
        G_C10_ = 2.0 * ib.C_feedback_hf    * fs_;
        G_C4_  = 2.0 * ib.C_output_couple  * fs_;
        G_C5_  = 2.0 * sus.C_wiper_couple  * fs_;
        G_C12_ = 2.0 * cl1.C_feedback_hf   * fs_;
        G_C6_  = 2.0 * cl1.C_diode_couple  * fs_;
        G_C13_ = 2.0 * cl1.C_output_couple * fs_;
        G_C11_ = 2.0 * cl2.C_feedback_hf   * fs_;
        G_C7_  = 2.0 * cl2.C_diode_couple  * fs_;
        G_C14_ = 2.0 * cl2.C_output_couple * fs_;   // C14: clipper2 -> tone top
        G_C3_  = 2.0 * ob.C_input_couple   * fs_;
        G_C2_  = 2.0 * ob.C_output_couple  * fs_;
        // Tone-network extras (not in CakeVariant).
        G_C9_  = 2.0 * sp_C9_  * fs_;
        G_C15_ = 2.0 * sp_C15_ * fs_;
        G_C16_ = 2.0 * sp_C16_ * fs_;
        G_C17_ = 2.0 * sp_C17_ * fs_;               // Q1 emitter bypass

        // Resistor conductances.
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
        G_R7_  = 1.0 / ob.R_bias_top;
        G_R3_  = 1.0 / ob.R_base_to_gnd;
        G_R6_  = 1.0 / ob.R_collector_load;
        G_R4_  = 1.0 / ob.R_emitter;
        // Tone-network resistors.
        G_R8_  = 1.0 / sp_R8_;
        G_R5_  = 1.0 / sp_R5_;
        G_R28_ = 1.0 / sp_R28_;
        G_R27_ = 1.0 / sp_R27_;
        G_R29_ = 1.0 / sp_R29_;
        G_R30_ = 1.0 / sp_R30_;
        G_R31_ = 1.0 / sp_R31_;

        // Clipper diode pairs (D1-D4) = 1N914.
        twoIs_d_ = 2.0 * cl1.diodes.Is;  Vt_d_ = cl1.diodes.N * V_T_300K;

        recompute_sustain_g();
        recompute_tone_g();
        recompute_volume_g();
        reset();
    }

    void set_sustain(double frac) noexcept { p_sus_ = std::clamp(frac,0.0,1.0); recompute_sustain_g(); }
    void set_tone(double frac)    noexcept { p_tone_= std::clamp(frac,0.0,1.0); recompute_tone_g(); }
    void set_volume(double frac)  noexcept { p_vol_ = std::clamp(frac,0.0,1.0); recompute_volume_g(); }

    void reset() noexcept {
        // Warm-start = the settled DC operating point (the model's own converged
        // fixed point, seeded so the 1uF/10uF couplers don't charge for tens of
        // ms after reset; see jumbo_goat.h).
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
        const double R = sp_Rvol_;   // 50k per the schematic (the netlist's 100k is wrong)
        G_Rvt_ = 1.0 / (R * (1.0 - p_vol_) + 1e-6);
        G_Rvb_ = 1.0 / (R * p_vol_         + 1e-6);
    }

    void init_cap_states() noexcept {
        s_C1_  = G_C1_  * (x_[0]  - 0.0);    // C1: n_c1 - in (in=0 at DC)
        s_C10_ = G_C10_ * (x_[1]  - x_[2]);
        s_C4_  = G_C4_  * (x_[2]  - x_[4]);
        s_C5_  = G_C5_  * (x_[5]  - x_[7]);
        s_C12_ = G_C12_ * (x_[8]  - x_[9]);
        s_C6_  = G_C6_  * (x_[8]  - x_[10]);
        s_C13_ = G_C13_ * (x_[9]  - x_[12]);
        s_C11_ = G_C11_ * (x_[13] - x_[14]);
        s_C7_  = G_C7_  * (x_[13] - x_[15]);
        s_C14_ = G_C14_ * (x_[14] - x_[17]);
        s_C9_  = G_C9_  * (x_[18] - x_[19]);
        s_C15_ = G_C15_ * (x_[22] - x_[23]);
        s_C16_ = G_C16_ * (x_[23] - x_[24]);
        s_C3_  = G_C3_  * (x_[20] - x_[25]);
        s_C2_  = G_C2_  * (x_[26] - x_[21]);
        s_C17_ = G_C17_ * (x_[27] - 0.0);
    }

    // Residual F(x) -> out. Full Ebers-Moll (Br=1) per transistor.
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
        auto capv = [&](double G, int p, double Vf, double s) noexcept {  // cap to fixed node
            out[p] += G * (x[p] - Vf) - s; };
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
        res (G_R2_, 0, 1);                     // R2: n_c1 - q4b
        resg(G_R14_, 1);                       // R14: q4b - gnd
        res (G_R9_, 1, 2);                     // R9:  q4b - q4c
        resv(G_R13_, 2, vcc);                  // R13: q4c - vcc
        resg(G_R22_, 3);                       // R22: q4e - gnd
        res (G_R19_, 7, 8);                    // R19: c5n - q3b
        resg(G_R23_, 6);                       // R23: sus_bot - gnd
        resg(G_R20_, 8);                       // R20: q3b - gnd
        res (G_R17_, 8, 9);                    // R17: q3b - q3c
        resv(G_R18_, 9, vcc);                  // R18: q3c - vcc
        resg(G_R21_, 11);                      // R21: q3e - gnd
        res (G_R12_, 12, 13);                  // R12: c13n - q2b
        resg(G_R16_, 13);                      // R16: q2b - gnd
        res (G_R15_, 13, 14);                  // R15: q2b - q2c
        resv(G_R11_, 14, vcc);                 // R11: q2c - vcc
        resg(G_R10_, 16);                      // R10: q2e - gnd
        res (G_R8_, 19, 20);                   // R8:  c9n - y
        res (G_R5_, 17, 22);                   // R5:  ta - x1
        res (G_R28_, 22, 23);                  // R28: x1 - x2
        res (G_R27_, 20, 23);                  // R27: y - x2
        res (G_R29_, 23, 24);                  // R29: x2 - x3
        res (G_R30_, 24, 21);                  // R30: x3 - tb
        resv(G_R7_, 25, vcc);                  // R7:  q1b - vcc
        resg(G_R3_, 25);                       // R3:  q1b - gnd
        resv(G_R6_, 26, vcc);                  // R6:  q1c - vcc
        resg(G_R4_, 27);                       // R4:  q1e - gnd
        res (G_R31_, 21, 28);                  // R31: tb - vol_top
        // Pots
        res(G_Rst_, 4, 5);   res(G_Rsb_, 5, 6);   // sustain
        res(G_Rtt_, 17, 18); res(G_Rtb_, 18, 21); // tone (top ta, wiper tw, bottom tb)
        res(G_Rvt_, 28, 29); resg(G_Rvb_, 29);    // volume top + bottom-to-gnd

        // --- Caps ---
        capv(G_C1_, 0, vin_, s_C1_);                                  // C1: n_c1 - in
        cap(G_C10_, 1, 2, s_C10_);  cap(G_C4_, 2, 4, s_C4_);
        cap(G_C5_, 5, 7, s_C5_);    cap(G_C12_, 8, 9, s_C12_);
        cap(G_C6_, 8, 10, s_C6_);   cap(G_C13_, 9, 12, s_C13_);
        cap(G_C11_, 13, 14, s_C11_);cap(G_C7_, 13, 15, s_C7_);
        cap(G_C14_, 14, 17, s_C14_);cap(G_C9_, 18, 19, s_C9_);
        cap(G_C15_, 22, 23, s_C15_);cap(G_C16_, 23, 24, s_C16_);
        cap(G_C3_, 20, 25, s_C3_);  cap(G_C2_, 26, 21, s_C2_);
        capg(G_C17_, 27, s_C17_);                                     // C17: q1e bypass

        // --- Diode pairs ---
        diode(10, 9);    // D3/D4: q3aux <-> q3c
        diode(15, 14);   // D1/D2: q2aux <-> q2c

        // --- Transistors (full EM) ---
        double IC,IB,IE;
        bjt(q4, x[1], x[2], x[3], IC,IB,IE);    out[1]+=IB; out[2]+=IC; out[3]-=IE;   // Q4
        bjt(q3, x[8], x[9], x[11], IC,IB,IE);   out[8]+=IB; out[9]+=IC; out[11]-=IE;  // Q3
        bjt(q2, x[13],x[14],x[16], IC,IB,IE);   out[13]+=IB;out[14]+=IC;out[16]-=IE;  // Q2
        bjt(q1, x[25],x[26],x[27], IC,IB,IE);   out[25]+=IB;out[26]+=IC;out[27]-=IE;  // Q1
    }

    // Analytic Jacobian dF/dx -- mirrors residual() element-by-element.
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

        // --- Resistors (mirror residual; resv/resg/capv/capg/cap all add G on diag pair) ---
        jr(G_R2_,0,1); jg(G_R14_,1); jr(G_R9_,1,2); jg(G_R13_,2); jg(G_R22_,3);
        jr(G_R19_,7,8); jg(G_R23_,6); jg(G_R20_,8); jr(G_R17_,8,9); jg(G_R18_,9);
        jg(G_R21_,11); jr(G_R12_,12,13); jg(G_R16_,13); jr(G_R15_,13,14); jg(G_R11_,14);
        jg(G_R10_,16); jr(G_R8_,19,20); jr(G_R5_,17,22); jr(G_R28_,22,23);
        jr(G_R27_,20,23); jr(G_R29_,23,24); jr(G_R30_,24,21);
        jg(G_R7_,25); jg(G_R3_,25); jg(G_R6_,26); jg(G_R4_,27); jr(G_R31_,21,28);
        // Pots
        jr(G_Rst_,4,5); jr(G_Rsb_,5,6); jr(G_Rtt_,17,18); jr(G_Rtb_,18,21);
        jr(G_Rvt_,28,29); jg(G_Rvb_,29);

        // --- Caps ---
        jg(G_C1_,0); jr(G_C10_,1,2); jr(G_C4_,2,4); jr(G_C5_,5,7); jr(G_C12_,8,9);
        jr(G_C6_,8,10); jr(G_C13_,9,12); jr(G_C11_,13,14); jr(G_C7_,13,15);
        jr(G_C14_,14,17); jr(G_C9_,18,19); jr(G_C15_,22,23); jr(G_C16_,23,24);
        jr(G_C3_,20,25); jr(G_C2_,26,21); jg(G_C17_,27);

        // --- Diode pairs: g = 2*Is*cosh(dV/Vt)/Vt ---
        jd(twoIs_d_ * std::cosh((x[10]-x[9])  / Vt_d_) / Vt_d_, 10, 9);   // D3/D4
        jd(twoIs_d_ * std::cosh((x[15]-x[14]) / Vt_d_) / Vt_d_, 15, 14);  // D1/D2

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
        bjt(q1, x[25],x[26],x[27],gBB,gBC,gBE, ibB,ibC,ibE); stampQ(25,26,27);
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
            if (!gen::solve_goat_s_od(J, b)) break;   // fixed-order sparse LU
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

        // Update bilinear cap states: s = 2*G*(Vp-Vn) - s.
        s_C1_  = 2.0*G_C1_ *(x_[0] -vin_)   - s_C1_;
        s_C10_ = 2.0*G_C10_*(x_[1] -x_[2])  - s_C10_;
        s_C4_  = 2.0*G_C4_ *(x_[2] -x_[4])  - s_C4_;
        s_C5_  = 2.0*G_C5_ *(x_[5] -x_[7])  - s_C5_;
        s_C12_ = 2.0*G_C12_*(x_[8] -x_[9])  - s_C12_;
        s_C6_  = 2.0*G_C6_ *(x_[8] -x_[10]) - s_C6_;
        s_C13_ = 2.0*G_C13_*(x_[9] -x_[12]) - s_C13_;
        s_C11_ = 2.0*G_C11_*(x_[13]-x_[14]) - s_C11_;
        s_C7_  = 2.0*G_C7_ *(x_[13]-x_[15]) - s_C7_;
        s_C14_ = 2.0*G_C14_*(x_[14]-x_[17]) - s_C14_;
        s_C9_  = 2.0*G_C9_ *(x_[18]-x_[19]) - s_C9_;
        s_C15_ = 2.0*G_C15_*(x_[22]-x_[23]) - s_C15_;
        s_C16_ = 2.0*G_C16_*(x_[23]-x_[24]) - s_C16_;
        s_C3_  = 2.0*G_C3_ *(x_[20]-x_[25]) - s_C3_;
        s_C2_  = 2.0*G_C2_ *(x_[26]-x_[21]) - s_C2_;
        s_C17_ = 2.0*G_C17_*(x_[27]-0.0)    - s_C17_;

        switch (tap) {
            case Tap::C2Out: return static_cast<float>(x_[14] - 4.50);  // clipper-2 collector AC
            case Tap::Out:   return static_cast<float>(x_[29]);
        }
        return static_cast<float>(x_[29]);
    }

    const CakeVariant& v_;
    double p_tone_, p_vol_, p_sus_;
    double fs_ = 48000.0;
    double vin_ = 0.0;

    // Tone-network + output extras (not in CakeVariant), from the schematic.
    double sp_R5_ = 10e3, sp_R8_ = 22e3, sp_R27_ = 10e3, sp_R28_ = 47e3;
    double sp_R29_ = 47e3, sp_R30_ = 10e3, sp_R31_ = 47e3, sp_Rvol_ = 50e3;
    double sp_C9_ = 0.0022e-6, sp_C15_ = 0.047e-6, sp_C16_ = 0.047e-6, sp_C17_ = 10e-6;

    double G_C1_=0,G_C10_=0,G_C4_=0,G_C5_=0,G_C12_=0,G_C6_=0,G_C13_=0,G_C11_=0,G_C7_=0;
    double G_C14_=0,G_C9_=0,G_C15_=0,G_C16_=0,G_C3_=0,G_C2_=0,G_C17_=0;
    double G_R2_=0,G_R14_=0,G_R9_=0,G_R13_=0,G_R22_=0,G_R19_=0,G_R23_=0,G_R20_=0,G_R17_=0,G_R18_=0,G_R21_=0;
    double G_R12_=0,G_R16_=0,G_R15_=0,G_R11_=0,G_R10_=0,G_R8_=0,G_R5_=0,G_R27_=0,G_R28_=0,G_R29_=0,G_R30_=0,G_R31_=0;
    double G_R7_=0,G_R3_=0,G_R6_=0,G_R4_=0;
    double G_Rst_=0,G_Rsb_=0,G_Rtt_=0,G_Rtb_=0,G_Rvt_=0,G_Rvb_=0;
    double twoIs_d_=0,Vt_d_=0;

    double s_C1_=0,s_C10_=0,s_C4_=0,s_C5_=0,s_C12_=0,s_C6_=0,s_C13_=0,s_C11_=0,s_C7_=0;
    double s_C14_=0,s_C9_=0,s_C15_=0,s_C16_=0,s_C3_=0,s_C2_=0,s_C17_=0;

    long dbg_iters_=0, dbg_calls_=0, dbg_nonconv_=0;
    double dbg_max_resid_=0.0;
    double dbg_max_jac_rel_=0.0;

    std::array<double,N> x_{};
    std::array<double,N> seed_ {
            0.666200784, 0.666200784, 4.52656716, 0.0440681276,   // n_c1 q4b q4c q4e
            3.21453977e-11,     3.21453977e-11,     3.18271264e-13,     0.697682625,    // sus_top sus_w sus_bot c5n
            0.697682625, 4.67922077, 4.67922077, 0.0765181974,    // q3b q3c q3aux q3e
            0.697682625,                               // c13n
            0.697682625, 4.67922077, 4.67922077, 0.0765181974,    // q2b q2c q2aux q2e
            -1.57266856e-11,     -1.29806922e-11,     -1.30430073e-11,     -1.30430106e-11,  -1.02346914e-11,  // ta tw c9n y tb
            -1.52532565e-11,     -1.30268202e-11,     -1.0724311e-11,                 // x1 x2 x3
            1.12420277, 6.21407315, 0.499219265,             // q1b q1c q1e
            -5.27561412e-12,     -2.63780706e-12                            // vol_top out
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
// GoatSOD is the DK description of this circuit (src/dk_circuits.h);
// GoatSODNode above is the reference it is checked against, forced to converge.
using GoatSOD = shipped_engine_t<DkGoatSOD>;   // the DK, or its nodal twin (engine_select.h)

} // namespace tocake
