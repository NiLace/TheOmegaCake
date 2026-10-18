// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// Bespoke joint solver for the "V2 Jumbo Goat" -- the biggest departure in the
// lineup: a SIX-transistor, THREE-diode-pair circuit with two parallel paths
// that recombine (by current summing) at the tone-stack wiper, and a GERMANIUM
// sub-stage. The fixed four-stage skeleton (and FullPedalV3) cannot represent it, so
// -- like the PM-3 and op-amp -- it gets its own solver.
//
// Architecture (see spice/the_jumbo_goat.cir, validated node-by-node):
//   MAIN fuzz path:
//     in -R2- C1 -> Q4 (booster) -C4-> Sustain -C5- R19 -> Q3 (clip D3/D4)
//        -R12- C13 -> Q2 (clip D1/D2) -C14-> tone_in -> passive tone stack
//   SIDE path (branched off tone_in):
//     tone_in -D6-> [C15 4.7u || R27 to gnd] -R28- Q5 (Ge; D5 clamps tone_in
//     AC127) -> Q6 (2N5172) collector, which is wired to the tone WIPER.
//   The tone wiper = Q6 collector = C3 input: the passive tone output and Q6's
//   collector current SUM here, then -C3-> Q1 (output booster) -C2-> Volume.
//
// Key behaviours (from the SPICE .op): Q2/Q3/Q5/Q6 emitters are grounded
// DIRECTLY (no emitter R). Q6 has NO collector pull-up -- it is biased ON (via
// R29 through its base) and SINKS current into the wiper, clamping it near 0 V
// at idle; Q5 (Ge) idles nearly OFF (Vbe ~13 mV) and is woken by the rectified
// envelope from D5/D6. So the Ge path is a signal-/envelope-driven mixer.
//
// Full Ebers-Moll (Br=1, B-C junction) for every transistor so saturation and
// the Ge leakage are modelled (the clean/clamped stages here run in/near
// saturation). 27-unknown joint Newton per sample; ANALYTIC Jacobian + fixed-order
// sparse LU (gen::solve_jumbo_goat).
//
// Main-path component values + the Si/diode device models come from the shared
// CakeVariant (the "V2 Jumbo Goat" spreadsheet tab). The SIDE-path extras
// (Q5 Ge model, R27/R28/R29, C15/C16, and C14 = clipper_2.C_output_couple) are
// carried as members here, exactly as the PM-3 carries its C15/C14/R27/R28.
//
// Unknowns x[0..26]:
//   0 q4b   1 q4c   2 q4e    3 n_c1(R2-C1)
//   4 sus_top  5 sus_w  6 sus_bot   7 c5n(C5-R19)
//   8 q3b   9 q3c   10 q3aux
//   11 c13n(R12-C13)
//   12 q2b  13 q2c  14 q2aux
//   15 tone_in  16 tt  17 tb  18 tw
//   19 d5k  20 q5b  21 q5c
//   22 q1b  23 q1c  24 q1e
//   25 vol_top  26 out

#pragma once

#include "dk_circuits.h"
#include "engine_select.h"

#include "newton_arbiter.h"

#include "cake_variant.h"
#include "bjt_dc_bias.h"
#include "nr_solver.h"
#include "solvers_generated.h"   // fixed-order sparse LU (gen::solve_jumbo_goat)

#include <algorithm>
#include <array>
#include <cmath>
#include "finite_guard.h"

namespace tocake {

class JumboGoatNode {
public:
    using T = double;
    static constexpr std::size_t N = 27;
    enum class Tap { C2Out, Out };

    explicit JumboGoatNode(const CakeVariant& v,
                       double p_tone = 0.5,
                       double p_vol  = 0.5,
                       double p_sus  = 1.0) noexcept
        : v_(v), p_tone_(p_tone), p_vol_(p_vol), p_sus_(p_sus) {}

    void prepare(double sample_rate) {
        fs_ = sample_rate;
        const auto& ib   = v_.input_booster;
        const auto& sus  = v_.sustain;
        const auto& cl1  = v_.clipper_1;
        const auto& cl2  = v_.clipper_2;
        const auto& tone = v_.tone;
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
        G_C14_ = 2.0 * cl2.C_output_couple * fs_;   // C14: clipper2 -> tone_in
        G_C9_  = 2.0 * tone.C_hp_series    * fs_;
        G_C8_  = 2.0 * tone.C_lp_shunt     * fs_;
        G_C3_  = 2.0 * ob.C_input_couple   * fs_;
        G_C2_  = 2.0 * ob.C_output_couple  * fs_;
        G_C17_ = 2.0 * sp_C17_ * fs_;               // C17: Q1 output Miller cap (q1c<->q1b)
        G_C15_ = 2.0 * sp_C15_ * fs_;               // side-path electrolytic
        G_C16_ = 2.0 * sp_C16_ * fs_;               // Q5 collector filter

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
        G_R12_ = 1.0 / cl2.R_input_series;
        G_R16_ = 1.0 / cl2.R_base_to_gnd;
        G_R15_ = 1.0 / cl2.R_feedback;
        G_R11_ = 1.0 / cl2.R_collector_load;
        G_R8_  = 1.0 / tone.R_lp_series;
        G_R5_  = 1.0 / tone.R_hp_shunt;
        G_R27_ = 1.0 / sp_R27_;
        G_R28_ = 1.0 / sp_R28_;
        G_R29_ = 1.0 / sp_R29_;
        G_R7_  = 1.0 / ob.R_bias_top;
        G_R3_  = 1.0 / ob.R_base_to_gnd;
        G_R6_  = 1.0 / ob.R_collector_load;
        G_R4_  = 1.0 / ob.R_emitter;

        // Clipper diode pairs (D1-D4) = 1N4148.
        twoIs_d_ = 2.0 * cl1.diodes.Is;  Vt_d_ = cl1.diodes.N * V_T_300K;
        // Gate rectifier pair (D5/D6) -- separate so it can be germanium.
        twoIs_g_ = 2.0 * gate_d_.Is;     Vt_g_ = gate_d_.N * V_T_300K;
        Is_g_    = gate_d_.Is;

        recompute_sustain_g();
        recompute_tone_g();
        recompute_volume_g();
        reset();
    }

    void set_sustain(double frac) noexcept { p_sus_ = std::clamp(frac,0.0,1.0); recompute_sustain_g(); }
    void set_tone(double frac)    noexcept { p_tone_= std::clamp(frac,0.0,1.0); recompute_tone_g(); }
    void set_volume(double frac)  noexcept { p_vol_ = std::clamp(frac,0.0,1.0); recompute_volume_g(); }

    void reset() noexcept {
        // Warm-start = the engine's own settled DC fixed point (digital silence
        // after reset() stays put). This MUST be accurate on the
        // big-cap nodes -- C15 (4.7u, d5k), C13 (.1u, c13n) and C16 (.1u, q5c) --
        // because a wrong seed there leaves those caps charging toward equilibrium
        // for tens of ms after reset() (R27.C15 ~= 32 ms), a slow startup
        // transient that does NOT exist in SPICE (it starts from its .op). Do not
        // hand-round these values: a seed rounded on c13n, d5k and q5c together
        // drops the SPICE null over the window from about -43 dB to -9 dB.
        // The values are the engine's own settled state: 30 s of digital silence at 96 kHz, default knobs,
        // Newton forced to 32 steps, diode auxiliaries set to their partner (unchanged between 30 and 60 s).
        // Re-derive them the same way whenever a device card of this circuit changes.
        x_ = {
            0.584046097, 3.87498456, 0.0185182842, -1.41272291e-15,        // q4b q4c q4e | n_c1
            2.77551743e-14,      2.77551743e-14,      2.74803706e-16,       0.56604282,    // sus_top sus_w sus_bot | c5n
            0.56604282, 3.78167604, 3.78167604,              // q3b q3c q3aux
            3.78167604,                                  // c13n
            0.56604282, 3.78167604, 3.78167604,              // q2b q2c q2aux
            0.0104455394, 0.00713808713, 0.0134305807, 0.0179533707,    // tone_in tt tb tw
            0.00696246798, 0.00718176389, 0.597161187,              // d5k q5b q5c
            1.65205574, 4.26404161, 1.07185552,              // q1b q1c q1e
            -1.74432773e-14,      -8.72163865e-15                              // vol_top out
        };
        // Companion cap states at the DC seed: s = G*(Vp - Vn).
        s_C1_  = G_C1_  * (x_[3]  - x_[0]);
        s_C10_ = G_C10_ * (x_[0]  - x_[1]);
        s_C4_  = G_C4_  * (x_[1]  - x_[4]);
        s_C5_  = G_C5_  * (x_[5]  - x_[7]);
        s_C12_ = G_C12_ * (x_[8]  - x_[9]);
        s_C6_  = G_C6_  * (x_[8]  - x_[10]);
        s_C13_ = G_C13_ * (x_[11] - x_[12]);
        s_C11_ = G_C11_ * (x_[12] - x_[13]);
        s_C7_  = G_C7_  * (x_[12] - x_[14]);
        s_C14_ = G_C14_ * (x_[13] - x_[15]);
        s_C9_  = G_C9_  * (x_[15] - x_[16]);
        s_C8_  = G_C8_  * (x_[17] - 0.0);
        s_C15_ = G_C15_ * (x_[19] - 0.0);
        s_C16_ = G_C16_ * (x_[21] - 0.0);
        s_C3_  = G_C3_  * (x_[18] - x_[22]);
        s_C2_  = G_C2_  * (x_[23] - x_[25]);
        s_C17_ = G_C17_ * (x_[22] - x_[23]);   // Q1 base - collector (cap open at DC)

        dbg_iters_ = dbg_calls_ = dbg_nonconv_ = 0;
        dbg_max_resid_ = 0.0;
    }

    float process_sample(float vin_volts) noexcept { return run_sample(vin_volts, Tap::Out); }
    float process_sample(float vin_volts, Tap tap) noexcept { return run_sample(vin_volts, tap); }

    void set_ge(double Is, double Bf, double Vaf) noexcept { q5_ge_ = {Is, Bf, Vaf}; }
    void set_gate_diode(double Is, double N) noexcept {
        gate_d_ = {Is, N, 1, 1}; twoIs_g_ = 2.0*Is; Is_g_ = Is; Vt_g_ = N*V_T_300K;
    }
    void set_r28(double r) noexcept { sp_R28_ = r; }   // call before prepare()
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

    // Residual F(x) -> out. Full Ebers-Moll (Br=1) per transistor; element-by-
    // element accumulation (each element's current added to both its nodes).
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
        // Cap current p->n: i = G*(Vp-Vn) - s; add +i at p, -i at n.
        auto cap = [&](double G, int p, int n, double s) noexcept {
            const double i = G * (x[p] - x[n]) - s;
            out[p] += i; out[n] -= i;
        };
        auto capg = [&](double G, int p, double s) noexcept {       // cap to gnd
            out[p] += G * x[p] - s;
        };
        auto res = [&](double G, int p, int n) noexcept {           // resistor p-n
            const double i = G * (x[p] - x[n]); out[p] += i; out[n] -= i;
        };
        auto resg = [&](double G, int p) noexcept { out[p] += G * x[p]; };      // R to gnd
        auto resv = [&](double G, int p, double Vfix) noexcept { out[p] += G*(x[p]-Vfix); }; // R to fixed
        auto diode = [&](int p, int n) noexcept {                   // clipper pair p<->n
            const double i = twoIs_d_ * std::sinh((x[p] - x[n]) / Vt_d_);
            out[p] += i; out[n] -= i;
        };
        // The two germanium gate diodes are ONE-WAY, not an antiparallel pair:
        // i = Is (e^(v/NVt) - 1).
        auto diode_1 = [&](int a, int k) noexcept {                 // anode a -> cathode k
            const double i = Is_g_ * (std::exp((x[a] - x[k]) / Vt_g_) - 1.0);
            out[a] += i; out[k] -= i;
        };
        auto diode_1g = [&](int k) noexcept {                       // anode GROUND -> cathode k
            const double i = Is_g_ * (std::exp(-x[k] / Vt_g_) - 1.0);
            out[k] -= i;
        };

        const auto& q4 = v_.input_booster.q;
        const auto& q3 = v_.clipper_1.q;
        const auto& q2 = v_.clipper_2.q;
        const auto& q1 = v_.output_booster.q;
        const auto& q6 = v_.input_booster.q;   // Q6 = 2N5172 (same Si part)

        // --- Resistors ---
        resv(G_R2_, 3, vin_);                  // R2: n_c1 - in
        resg(G_R14_, 0);                       // R14: q4b - gnd
        res (G_R9_, 0, 1);                     // R9:  q4b - q4c
        resv(G_R13_, 1, vcc);                  // R13: q4c - vcc
        resg(G_R22_, 2);                       // R22: q4e - gnd
        res (G_R19_, 7, 8);                    // R19: c5n - q3b
        resg(G_R23_, 6);                       // R23: sus_bot - gnd
        resg(G_R20_, 8);                       // R20: q3b - gnd
        res (G_R17_, 8, 9);                    // R17: q3b - q3c
        resv(G_R18_, 9, vcc);                  // R18: q3c - vcc
        res (G_R12_, 11, 9);                   // R12: c13n - q3c
        resg(G_R16_, 12);                      // R16: q2b - gnd
        res (G_R15_, 12, 13);                  // R15: q2b - q2c
        resv(G_R11_, 13, vcc);                 // R11: q2c - vcc
        res (G_R8_, 15, 17);                   // R8:  tone_in - tb
        resg(G_R5_, 16);                       // R5:  tt - gnd
        resg(G_R27_, 19);                      // R27: d5k - gnd
        res (G_R28_, 19, 20);                  // R28: d5k - q5b
        resv(G_R29_, 21, vcc);                 // R29: q5c - vcc
        resv(G_R7_, 22, vcc);                  // R7:  q1b - vcc
        resg(G_R3_, 22);                       // R3:  q1b - gnd
        resv(G_R6_, 23, vcc);                  // R6:  q1c - vcc
        resg(G_R4_, 24);                       // R4:  q1e - gnd
        // Pots
        res(G_Rst_, 4, 5);  res(G_Rsb_, 5, 6); // sustain top/bottom
        res(G_Rtt_, 16, 18); res(G_Rtb_, 18, 17); // tone top/bottom (wiper=18)
        res(G_Rvt_, 25, 26); resg(G_Rvb_, 26); // volume top + bottom-to-gnd

        // --- Caps ---
        cap(G_C1_, 3, 0, s_C1_);    cap(G_C10_, 0, 1, s_C10_);  cap(G_C4_, 1, 4, s_C4_);
        cap(G_C5_, 5, 7, s_C5_);    cap(G_C12_, 8, 9, s_C12_);  cap(G_C6_, 8, 10, s_C6_);
        cap(G_C13_, 11, 12, s_C13_);cap(G_C11_, 12, 13, s_C11_);cap(G_C7_, 12, 14, s_C7_);
        cap(G_C14_, 13, 15, s_C14_);cap(G_C9_, 15, 16, s_C9_);  capg(G_C8_, 17, s_C8_);
        capg(G_C15_, 19, s_C15_);   capg(G_C16_, 21, s_C16_);
        cap(G_C3_, 18, 22, s_C3_);  cap(G_C2_, 23, 25, s_C2_);
        cap(G_C17_, 22, 23, s_C17_);   // C17 Q1 collector<->base Miller

        // --- Diode pairs ---
        diode(10, 9);    // D3/D4: q3aux <-> q3c
        diode(14, 13);   // D1/D2: q2aux <-> q2c
        diode_1(15, 19); // D6: tone_in -> d5k (rectifies into C15 || R27)
        diode_1g(15);    // D5: gnd -> tone_in (clamps the negative half)

        // --- Transistors (full EM) ---
        double IC,IB,IE;
        bjt(q4, x[0], x[1], x[2], IC,IB,IE);  out[0]+=IB; out[1]+=IC; out[2]-=IE;   // Q4
        bjt(q3, x[8], x[9], 0.0, IC,IB,IE);   out[8]+=IB; out[9]+=IC;               // Q3 (e=gnd)
        bjt(q2, x[12],x[13],0.0, IC,IB,IE);   out[12]+=IB;out[13]+=IC;              // Q2 (e=gnd)
        bjt(q5_ge_, x[20],x[21],0.0, IC,IB,IE);out[20]+=IB;out[21]+=IC;             // Q5 Ge (e=gnd)
        bjt(q6, x[21],x[18],0.0, IC,IB,IE);   out[21]+=IB;out[18]+=IC;              // Q6 (b=q5c,c=tw,e=gnd)
        bjt(q1, x[22],x[23],x[24],IC,IB,IE);  out[22]+=IB;out[23]+=IC;out[24]-=IE;  // Q1
    }

    // Analytic Jacobian dF/dx -- mirrors residual() element-by-element. Linear
    // elements (resistors, companion caps, pots) contribute CONSTANT conductances
    // (the cap history term s is constant within a Newton solve); diode pairs and
    // transistors contribute their forward-EM partials (same model as residual(),
    // partials identical to coloured_goat.h's, validated there). Verified
    // entry-wise vs a central-difference J, and by the SPICE null.
    void jacobian(const std::array<double,N>& x,
                  std::array<std::array<double,N>,N>& out) const noexcept {
        const double V_T = V_T_300K;
        for (auto& row : out) row.fill(0.0);

        auto jr = [&](double G, int p, int n) noexcept {          // resistor/cap p-n
            out[p][p]+=G; out[p][n]-=G; out[n][p]-=G; out[n][n]+=G; };
        auto jg = [&](double G, int p) noexcept { out[p][p]+=G; };  // to gnd or fixed rail
        auto jd = [&](double g, int p, int n) noexcept {          // diode-pair conductance
            out[p][p]+=g; out[p][n]-=g; out[n][p]-=g; out[n][n]+=g; };

        // Forward-EM partials (gB*=dIC/d*, ib*=dIB/d*); identical to coloured_goat.h.
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

        // --- Resistors (mirror residual) ---
        jg(G_R2_,3); jg(G_R14_,0); jr(G_R9_,0,1); jg(G_R13_,1); jg(G_R22_,2);
        jr(G_R19_,7,8); jg(G_R23_,6); jg(G_R20_,8); jr(G_R17_,8,9); jg(G_R18_,9);
        jr(G_R12_,11,9); jg(G_R16_,12); jr(G_R15_,12,13); jg(G_R11_,13);
        jr(G_R8_,15,17); jg(G_R5_,16); jg(G_R27_,19); jr(G_R28_,19,20); jg(G_R29_,21);
        jg(G_R7_,22); jg(G_R3_,22); jg(G_R6_,23); jg(G_R4_,24);
        // Pots
        jr(G_Rst_,4,5); jr(G_Rsb_,5,6); jr(G_Rtt_,16,18); jr(G_Rtb_,18,17);
        jr(G_Rvt_,25,26); jg(G_Rvb_,26);

        // --- Caps (companion conductance; history term constant in x) ---
        jr(G_C1_,3,0); jr(G_C10_,0,1); jr(G_C4_,1,4); jr(G_C5_,5,7);
        jr(G_C12_,8,9); jr(G_C6_,8,10); jr(G_C13_,11,12); jr(G_C11_,12,13);
        jr(G_C7_,12,14); jr(G_C14_,13,15); jr(G_C9_,15,16); jg(G_C8_,17);
        jg(G_C15_,19); jg(G_C16_,21); jr(G_C3_,18,22); jr(G_C2_,23,25);
        jr(G_C17_,22,23);

        // --- Diode pairs: g = 2*Is*cosh(dV/Vt)/Vt ---
        jd(twoIs_d_ * std::cosh((x[10]-x[9])  / Vt_d_) / Vt_d_, 10, 9);   // D3/D4
        jd(twoIs_d_ * std::cosh((x[14]-x[13]) / Vt_d_) / Vt_d_, 14, 13);  // D1/D2
        jd(Is_g_ * std::exp((x[15]-x[19]) / Vt_g_) / Vt_g_, 15, 19);      // D6 rectifier
        jg(Is_g_ * std::exp(-x[15] / Vt_g_) / Vt_g_, 15);                 // D5 clamp to ground

        // --- Transistors (same node wiring as residual) ---
        const auto& q4 = v_.input_booster.q;
        const auto& q3 = v_.clipper_1.q;
        const auto& q2 = v_.clipper_2.q;
        const auto& q1 = v_.output_booster.q;
        const auto& q6 = v_.input_booster.q;   // Q6 = 2N5172 (same Si part)
        double gBB,gBC,gBE,ibB,ibC,ibE;
        // Q4 (full: b0 c1 e2)
        bjt(q4, x[0],x[1],x[2], gBB,gBC,gBE, ibB,ibC,ibE);
        out[0][0]+=ibB; out[0][1]+=ibC; out[0][2]+=ibE;
        out[1][0]+=gBB; out[1][1]+=gBC; out[1][2]+=gBE;
        out[2][0]-=gBB+ibB; out[2][1]-=gBC+ibC; out[2][2]-=gBE+ibE;
        // Q3 (grounded e: b8 c9)
        bjt(q3, x[8],x[9],0.0, gBB,gBC,gBE, ibB,ibC,ibE);
        out[8][8]+=ibB; out[8][9]+=ibC; out[9][8]+=gBB; out[9][9]+=gBC;
        // Q2 (grounded e: b12 c13)
        bjt(q2, x[12],x[13],0.0, gBB,gBC,gBE, ibB,ibC,ibE);
        out[12][12]+=ibB; out[12][13]+=ibC; out[13][12]+=gBB; out[13][13]+=gBC;
        // Q5 Ge (grounded e: b20 c21)
        bjt(q5_ge_, x[20],x[21],0.0, gBB,gBC,gBE, ibB,ibC,ibE);
        out[20][20]+=ibB; out[20][21]+=ibC; out[21][20]+=gBB; out[21][21]+=gBC;
        // Q6 (grounded e: b21 c18) -- base node 21 also carries Q5's collector stamp
        bjt(q6, x[21],x[18],0.0, gBB,gBC,gBE, ibB,ibC,ibE);
        out[21][21]+=ibB; out[21][18]+=ibC; out[18][21]+=gBB; out[18][18]+=gBC;
        // Q1 (full: b22 c23 e24)
        bjt(q1, x[22],x[23],x[24], gBB,gBC,gBE, ibB,ibC,ibE);
        out[22][22]+=ibB; out[22][23]+=ibC; out[22][24]+=ibE;
        out[23][22]+=gBB; out[23][23]+=gBC; out[23][24]+=gBE;
        out[24][22]-=gBB+ibB; out[24][23]-=gBC+ibC; out[24][24]-=gBE+ibE;
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
            jacobian(x_, J);                       // analytic
            std::array<double,N> b{};
            for (std::size_t i=0;i<N;++i) b[i] = -Fv[i];
            if (!gen::solve_jumbo_goat(J, b)) break;   // fixed-order sparse LU
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
        s_C1_  = 2.0*G_C1_ *(x_[3] -x_[0])  - s_C1_;
        s_C10_ = 2.0*G_C10_*(x_[0] -x_[1])  - s_C10_;
        s_C4_  = 2.0*G_C4_ *(x_[1] -x_[4])  - s_C4_;
        s_C5_  = 2.0*G_C5_ *(x_[5] -x_[7])  - s_C5_;
        s_C12_ = 2.0*G_C12_*(x_[8] -x_[9])  - s_C12_;
        s_C6_  = 2.0*G_C6_ *(x_[8] -x_[10]) - s_C6_;
        s_C13_ = 2.0*G_C13_*(x_[11]-x_[12]) - s_C13_;
        s_C11_ = 2.0*G_C11_*(x_[12]-x_[13]) - s_C11_;
        s_C7_  = 2.0*G_C7_ *(x_[12]-x_[14]) - s_C7_;
        s_C14_ = 2.0*G_C14_*(x_[13]-x_[15]) - s_C14_;
        s_C9_  = 2.0*G_C9_ *(x_[15]-x_[16]) - s_C9_;
        s_C8_  = 2.0*G_C8_ *(x_[17]-0.0)    - s_C8_;
        s_C15_ = 2.0*G_C15_*(x_[19]-0.0)    - s_C15_;
        s_C16_ = 2.0*G_C16_*(x_[21]-0.0)    - s_C16_;
        s_C3_  = 2.0*G_C3_ *(x_[18]-x_[22]) - s_C3_;
        s_C2_  = 2.0*G_C2_ *(x_[23]-x_[25]) - s_C2_;
        s_C17_ = 2.0*G_C17_*(x_[22]-x_[23]) - s_C17_;

        switch (tap) {
            case Tap::C2Out: return static_cast<float>(x_[13] - 4.020);  // clipper-2 collector AC
            case Tap::Out:   return static_cast<float>(x_[26]);
        }
        return static_cast<float>(x_[26]);
    }

    const CakeVariant& v_;
    double p_tone_, p_vol_, p_sus_;
    double fs_ = 48000.0;
    double vin_ = 0.0;

    // Side-path component values (not in CakeVariant), from the netlist /
    // "V2 Jumbo Goat" tab + the published schematic.
    //   q5_ge_ : Q5 germanium AC127, Is 1e-6: a leaky part, set so the gate opens on playing and closes in
    //     silence on a real guitar signal; the Pro Electron sheet gives hFE >= 50 and ICBO <= 10 uA but no VBE,
    //     so it cannot settle Is. On the SPICE reference the gate passes full fuzz for any Is from 1e-8 to
    //     1e-4; only 1e-9 closes it.
    //   sp_R28_ = 220 ohm. Not 220k: a 1000x error there is a major cause of
    //     a silent gate.
    bjt_params q5_ge_ { 1.0e-6, 100.0, 50.0 };   // AC127 (Ge), as dk_circuits.h kQ5Ge
    double sp_R27_ = 6.8e3, sp_R28_ = 220.0, sp_R29_ = 22e3;
    double sp_C15_ = 4.7e-6, sp_C16_ = 0.1e-6;
    double sp_C17_ = 470e-12;   // Q1 output-stage collector->base Miller cap (schematic C17)

    double G_C1_=0,G_C10_=0,G_C4_=0,G_C5_=0,G_C12_=0,G_C6_=0,G_C13_=0,G_C11_=0,G_C7_=0;
    double G_C14_=0,G_C9_=0,G_C8_=0,G_C15_=0,G_C16_=0,G_C3_=0,G_C2_=0,G_C17_=0;
    double G_R2_=0,G_R14_=0,G_R9_=0,G_R13_=0,G_R22_=0,G_R19_=0,G_R23_=0,G_R20_=0,G_R17_=0,G_R18_=0;
    double G_R12_=0,G_R16_=0,G_R15_=0,G_R11_=0,G_R8_=0,G_R5_=0,G_R27_=0,G_R28_=0,G_R29_=0;
    double G_R7_=0,G_R3_=0,G_R6_=0,G_R4_=0;
    double G_Rst_=0,G_Rsb_=0,G_Rtt_=0,G_Rtb_=0,G_Rvt_=0,G_Rvb_=0;
    double twoIs_d_=0,Vt_d_=0;
    // Gate diodes (D5/D6) -- GERMANIUM, ONE-WAY each: D6 rectifies into the
    // envelope node, D5 clamps tone_in against ground. Vf ~0.25 V; modelled as a
    // 1N34A-class Ge diode (Is~2e-7, N~1.06). Exact part unknown; BY-EAR /
    // data-pending like the AC127.
    diode_params gate_d_ { 2.0e-7, 1.06, 1, 1 };
    double twoIs_g_=0, Is_g_=0, Vt_g_=0;

    double s_C1_=0,s_C10_=0,s_C4_=0,s_C5_=0,s_C12_=0,s_C6_=0,s_C13_=0,s_C11_=0,s_C7_=0;
    double s_C14_=0,s_C9_=0,s_C8_=0,s_C15_=0,s_C16_=0,s_C3_=0,s_C2_=0,s_C17_=0;

    long dbg_iters_=0, dbg_calls_=0, dbg_nonconv_=0;
    double dbg_max_resid_=0.0;
    double dbg_max_jac_rel_=0.0;

    std::array<double,N> x_{};
    // Times the guard fired. Exposed so a harness can SEE it: a guard
    // nothing can read is an alarm with no wire.
    long dbg_nonfinite_ = 0;
public:
    long nonfinite_resets() const noexcept { return dbg_nonfinite_; }
private:
};

// ---- the engine that ships -------------------------------------------------
// JumboGoat is the DK description of this circuit (src/dk_circuits.h);
// JumboGoatNode above is the reference it is checked against, forced to converge.
using JumboGoat = shipped_engine_t<DkJumboGoat>;   // the DK, or its nodal twin (engine_select.h)

} // namespace tocake
