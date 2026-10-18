// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// Bespoke joint solver for the Pyramid Master PM-3 (a 1970s Japanese fuzz,
// boost removed) -- a UNIQUE one-off topology, so it gets its own solver rather
// than complicating the shared FullPedalV3. Same 4-transistor core,
// but with five deltas vs the V3 family (see spice/the_pyramid_master.cir):
//   1. C15 in PARALLEL with the input series R2 (bright cap)  -> explicit in_c1
//   2. R27 from input-booster output (sus_top) to GND
//   3. C14 AC-couples the clipper-2 collector to the tone stack -> explicit tin
//   4. VOLUME pot relocated: tone wiper -> Vol top; Vol wiper -> C3 -> Q1 base
//   5. R28 output pulldown; C2 drives the jack directly (no post-booster pot)
//
// Per-sample solve: Newton-Raphson on the residual F with a hand-derived
// ANALYTIC Jacobian J (validated entry-wise vs central differences). The 10 linear nodes are eliminated by a SCHUR
// complement so each Newton step solves a 14x14 instead of 24x24 (the linear
// block A_LL is factored once per knob move). The elimination is exact, so the
// audio matches the full solve to rounding. F is nulled against a SPICE deck of
// the same circuit.
//
// Unknowns x[0..23]:
//   0: V_B(Q4)   1: V_C(Q4)   2: V_E(Q4)        input booster
//   3: V_B(Q3)   4: V_C(Q3)   5: V_E(Q3)        clipper 1
//   6: V_d1                                      clipper-1 diode node
//   7: V_B(Q2)   8: V_C(Q2)   9: V_E(Q2)        clipper 2
//  10: V_d2                                      clipper-2 diode node
//  11: V_tone_top  12: V_tone_bot  13: V_tone_w  passive tone stack
//  14: V_B(Q1)  15: V_C(Q1)  16: V_E(Q1)        output booster
//  17: V_sus_top  18: V_sus_w  19: V_sus_bot     sustain pot
//  20: V_in_c1                                   input node (R2||C15 -> C1)
//  21: V_tin                                     tone-stack input (post-C14)
//  22: V_vol_w                                   relocated volume wiper
//  23: V_out                                     output (post-C2, R28 to GND)
//
// Caps with bilinear (trapezoidal) state (15): C15 C1 C10 C4 C5 C12 C6 C13
//   C11 C7 C14 C9 C8 C3 C2.  Series C5->R19 and C13->R12 are still eliminated
//   internally (Norton form, rooted at sus_w and q3c) exactly as in FullPedalV3.

#pragma once

#include "dk_circuits.h"
#include "engine_select.h"

#include "newton_arbiter.h"

#include "cake_variant.h"
#include "bjt_dc_bias.h"
#include "nr_solver.h"
#include "solvers_generated.h"   // fixed-order sparse LU (gen::solve_pm3_schur)
#include <array>
#include <cmath>
#include <algorithm>
#include "finite_guard.h"

namespace tocake {

class PyramidMasterNode {
public:
    static constexpr std::size_t N = 24;
    enum class Tap { IbOut, C1Out, C2Out, Out };

    explicit PyramidMasterNode(const CakeVariant& v,
                           double p_tone = 0.5,
                           double p_vol  = 0.5,
                           double p_sus  = 1.0) noexcept
        : v_(v), p_tone_(p_tone), p_vol_(p_vol), p_sus_(p_sus) {}

    void prepare(double sample_rate) {
        // Rest seed from the WHOLE circuit, not from each stage on its own.
        // The per-stage DC solves ignore how the stages load each other, so
        // they land several millivolts away from the engine's rest, which
        // would start every run with a transient the circuit does not
        // have. At a vanishing sample rate every companion conductance 2*C*fs
        // is zero -- the capacitors are open, i.e. DC -- so this engine's own
        // Newton, fed silence, lands on the full-circuit operating point with
        // exactly the equations it runs. reset() stays cheap: it only reads
        // the adopted values.
        setup_rate(kDcSeedRate);
        reset();
        std::array<double, N> x_prev = x_;
        bool settled = false;
        for (int k = 0; k < kDcSeedMaxPasses && !settled; ++k) {
            (void) run_sample(0.0, Tap::Out);
            double step = 0.0;
            for (std::size_t n = 0; n < N; ++n) step = std::max(step, std::fabs(x_[n] - x_prev[n]));
            x_prev = x_;
            settled = tocake::all_finite(x_) && step < kDcSeedTol;
        }
        const std::array<double, N> x = x_;
        setup_rate(sample_rate);
        if (settled) {
        bias_ib_.V_B = x[0];  bias_ib_.V_C = x[1];  bias_ib_.V_E = x[2];
        bias_cl1_.V_B = x[3];  bias_cl1_.V_C = x[4];  bias_cl1_.V_E = x[5];
        bias_cl2_.V_B = x[7];  bias_cl2_.V_C = x[8];  bias_cl2_.V_E = x[9];
        bias_ob_.V_B = x[14];  bias_ob_.V_C = x[15];  bias_ob_.V_E = x[16];
        }
        dc_seed_settled_ = settled;
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

        // Per-stage DC bias is UNAFFECTED by the PM-3 deltas (R27/R28/C14/C15
        // are all AC-coupled or shunt the AC-grounded nodes), so the standard
        // per-stage bias networks apply -- confirmed against the SPICE .op.
        bias_ib_  = solve_dc_bias(ib.q,  bias_network_for(ib),  v_.vcc, bjt_law::GummelPoon);
        bias_cl1_ = solve_dc_bias(cl1.q, bias_network_for(cl1), v_.vcc, bjt_law::GummelPoon);
        bias_cl2_ = solve_dc_bias(cl2.q,
                                  bias_network_for(cl2, tone_stack_dc_load(tone)),
                                  v_.vcc, bjt_law::GummelPoon);
        bias_ob_  = solve_dc_bias(ob.q,  bias_network_for(ob),  v_.vcc, bjt_law::GummelPoon);

        // Cap conductances (G = 2*C*fs). C15/C14 are PM-3-only (pm3_ fields).
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
        // PM-3-only caps carried in dedicated fields on the variant.
        G_C15_ = 2.0 * pm3_C15_ * fs_;
        G_C14_ = 2.0 * pm3_C14_ * fs_;

        // Input network: in -[R2 || C15]- in_c1 -C1- q4b  (in_c1 explicit).
        G_R2_ = 1.0 / ib.R_input_series;

        // PM-3 shunt resistors.
        G_R27_ = 1.0 / pm3_R27_;
        G_R28_ = 1.0 / pm3_R28_;

        G_R23_ = 1.0 / sus.R_pot_bottom;
        G_R19_ = 1.0 / cl1.R_input_series;
        const double G_C5R19_denom = G_C5_ + G_R19_;
        G_eff_C5R19_ = G_C5_ * G_R19_ / G_C5R19_denom;
        s_C5_denom_  = G_C5R19_denom;

        G_R12_          = 1.0 / cl2.R_input_series;
        G_C13R12_denom_ = G_C13_ + G_R12_;
        G_eff_C13R12_   = G_C13_ * G_R12_ / G_C13R12_denom_;

        G_R5_  = 1.0 / tone.R_hp_shunt;
        G_R8_  = 1.0 / tone.R_lp_series;

        G_R7_ = 1.0 / ob.R_bias_top;
        G_R3_ = 1.0 / ob.R_base_to_gnd;
        G_R6_ = 1.0 / ob.R_collector_load;
        G_R4_ = 1.0 / ob.R_emitter;

        twoIs_d1_ = 2.0 * cl1.diodes.Is;  Vt_d1_ = cl1.diodes.N * V_T_300K;
        twoIs_d2_ = 2.0 * cl2.diodes.Is;  Vt_d2_ = cl2.diodes.N * V_T_300K;

        recompute_sustain_g();
        recompute_tone_g();
        recompute_volume_g();
    }

    void set_sustain(double frac) noexcept {
        p_sus_ = std::clamp(frac, 0.0, 1.0);  recompute_sustain_g();
    }
    void set_tone(double frac) noexcept {
        p_tone_ = std::clamp(frac, 0.0, 1.0); recompute_tone_g();
    }
    void set_volume(double frac) noexcept {
        p_vol_ = std::clamp(frac, 0.0, 1.0);  recompute_volume_g();
    }

    void reset() noexcept {
        // Tone stack DC ladder (same as FullPedalV3) -- but PM-3 AC-couples the
        // tone stack via C14, so the whole tone/volume/output section sits at
        // 0 V DC. The transistor biases are unchanged.
        V_tt0_ = 0.0; V_tb0_ = 0.0; V_tw0_ = 0.0;

        x_ = { bias_ib_.V_B,  bias_ib_.V_C,  bias_ib_.V_E,
               bias_cl1_.V_B, bias_cl1_.V_C, bias_cl1_.V_E, bias_cl1_.V_C,
               bias_cl2_.V_B, bias_cl2_.V_C, bias_cl2_.V_E, bias_cl2_.V_C,
               V_tt0_, V_tb0_, V_tw0_,
               bias_ob_.V_B,  bias_ob_.V_C,  bias_ob_.V_E,
               0.0, 0.0, 0.0,     // sustain: AC-coupled -> 0 V at DC
               0.0,               // in_c1: tied to in(0) via R2 at DC -> 0 V
               0.0,               // tin: C14-blocked -> 0 V at DC
               0.0, 0.0 };        // vol_w, out: AC-coupled -> 0 V at DC

        // Cap states (trapezoidal companion): s = G*(V_p - V_n) at DC.
        s_C15_ = G_C15_ * (0.0           - 0.0);            // in - in_c1, both 0
        s_C1_  = G_C1_  * (0.0           - bias_ib_.V_B);   // in_c1 - q4b
        s_C10_ = G_C10_ * (bias_ib_.V_B  - bias_ib_.V_C);
        s_C4_  = G_C4_  * (bias_ib_.V_C  - 0.0);            // q4c - sus_top
        s_C5_  = G_C5_  * (0.0           - bias_cl1_.V_B);
        s_C12_ = G_C12_ * (bias_cl1_.V_B - bias_cl1_.V_C);
        s_C6_  = G_C6_  * (bias_cl1_.V_B - bias_cl1_.V_C);
        s_C13_ = G_C13_ * (bias_cl1_.V_C - bias_cl2_.V_B);
        s_C11_ = G_C11_ * (bias_cl2_.V_B - bias_cl2_.V_C);
        s_C7_  = G_C7_  * (bias_cl2_.V_B - bias_cl2_.V_C);
        s_C14_ = G_C14_ * (bias_cl2_.V_C - 0.0);            // q2c - tin
        s_C9_  = G_C9_  * (0.0           - V_tt0_);         // tin - tone_top
        s_C8_  = G_C8_  * (V_tb0_        - 0.0);
        s_C3_  = G_C3_  * (0.0           - bias_ob_.V_B);   // vol_w - q1b
        s_C2_  = G_C2_  * (bias_ob_.V_C  - 0.0);            // q1c - out
    }

    float process_sample(float vin_volts) noexcept {
        return run_sample(static_cast<double>(vin_volts), Tap::Out);
    }
    float process_sample(float vin_volts, Tap tap) noexcept {
        return run_sample(static_cast<double>(vin_volts), tap);
    }

    const dc_bias_result& bias_ib()  const noexcept { return bias_ib_;  }
    const dc_bias_result& bias_cl1() const noexcept { return bias_cl1_; }
    const dc_bias_result& bias_cl2() const noexcept { return bias_cl2_; }
    const dc_bias_result& bias_ob()  const noexcept { return bias_ob_;  }
    const std::array<double, N>& state() const noexcept { return x_; }
    double dbg_max_resid() const noexcept { return dbg_max_resid_; }
    long   dbg_noconv()    const noexcept { return dbg_noconv_; }
    long   dbg_iters()     const noexcept { return dbg_iters_; }

private:
    void recompute_sustain_g() noexcept {
        const double R = v_.sustain.R_pot;
        G_Rst_ = 1.0 / (R * (1.0 - p_sus_) + 1e-6);
        G_Rsb_ = 1.0 / (R * p_sus_         + 1e-6);
        schur_dirty_ = true;
    }
    void recompute_tone_g() noexcept {
        const double R = v_.tone.R_pot;
        R_tt_ = R * (1.0 - p_tone_) + 1e-6;
        R_tb_ = R * p_tone_         + 1e-6;
        G_Rtt_ = 1.0 / R_tt_;
        G_Rtb_ = 1.0 / R_tb_;
        schur_dirty_ = true;
    }
    void recompute_volume_g() noexcept {
        const double R = v_.volume.R_pot;
        G_Rvt_ = 1.0 / (R * (1.0 - p_vol_) + 1e-6);
        G_Rvb_ = 1.0 / (R * p_vol_         + 1e-6);
        schur_dirty_ = true;
    }

    float run_sample(double v_in, Tap tap) noexcept {
        const double R14 = v_.input_booster.R_base_to_gnd;
        const double R9  = v_.input_booster.R_feedback;
        const double R13 = v_.input_booster.R_collector_load;
        const double R22 = v_.input_booster.R_emitter;
        const double R20 = v_.clipper_1.R_base_to_gnd;
        const double R17 = v_.clipper_1.R_feedback;
        const double Rc1 = v_.clipper_1.R_collector_load;
        const double R21 = v_.clipper_1.R_emitter;
        const double R16 = v_.clipper_2.R_base_to_gnd;
        const double R15 = v_.clipper_2.R_feedback;
        const double Rc2 = v_.clipper_2.R_collector_load;
        const double R10 = v_.clipper_2.R_emitter;
        const auto&  q4 = v_.input_booster.q;
        const auto&  q3 = v_.clipper_1.q;
        const auto&  q2 = v_.clipper_2.q;
        const auto&  q1 = v_.output_booster.q;
        const double V_T = V_T_300K;
        const double vcc = v_.vcc;

        // Effective companion states for the two eliminated series paths.
        const double s_C5_eff  = s_C5_  * G_R19_ / s_C5_denom_;
        const double s_C13_eff = s_C13_ * G_R12_ / G_C13R12_denom_;

        auto F = [&](const std::array<double, N>& x,
                     std::array<double, N>&       out) noexcept
        {
            const double V_B4 = x[0],  V_C4 = x[1],  V_E4 = x[2];
            const double V_B3 = x[3],  V_C3 = x[4],  V_E3 = x[5];
            const double V_d1 = x[6];
            const double V_B2 = x[7],  V_C2 = x[8],  V_E2 = x[9];
            const double V_d2 = x[10];
            const double V_tt = x[11], V_tb = x[12], V_tw = x[13];
            const double V_B1 = x[14], V_C1 = x[15], V_E1 = x[16];
            const double V_st = x[17], V_sw = x[18], V_sb = x[19];
            const double V_inc1 = x[20], V_tin = x[21];
            const double V_volw = x[22], V_out = x[23];

            // Full Ebers-Moll transport model WITH the base-collector junction
            // (reverse current, Br=1 = SPICE default) so SATURATION matches the
            // SPICE NPN model. PM-3's 150k-load input booster is driven hard
            // into saturation; a forward-only model would miss that clip
            // shape. Forward Early factor as in FullPedalV3.
            auto bjt = [V_T](const bjt_params& q, double Vb, double Vc, double Ve,
                             double& IC, double& IB, double& IE) noexcept {
                constexpr double Br = 1.0;
                const double ef    = std::exp((Vb - Ve) / V_T);   // fwd BE junction
                const double er    = std::exp((Vb - Vc) / V_T);   // rev BC junction
                const double early = 1.0 + (Vc - Vb) / q.Vaf;
                IC = q.Is * (ef - er) * early - (q.Is / Br) * (er - 1.0);
                IB = (q.Is / q.Bf) * (ef - 1.0) + (q.Is / Br) * (er - 1.0);
                IE = IC + IB;
            };
            double IC4, IB4, IE4; bjt(q4, V_B4, V_C4, V_E4, IC4, IB4, IE4);
            double IC3, IB3, IE3; bjt(q3, V_B3, V_C3, V_E3, IC3, IB3, IE3);
            double IC2, IB2, IE2; bjt(q2, V_B2, V_C2, V_E2, IC2, IB2, IE2);
            double IC1, IB1, IE1; bjt(q1, V_B1, V_C1, V_E1, IC1, IB1, IE1);

            const double i_d1 = twoIs_d1_ * std::sinh((V_d1 - V_C3) / Vt_d1_);
            const double i_d2 = twoIs_d2_ * std::sinh((V_d2 - V_C2) / Vt_d2_);

            // Eliminated series-path currents (into the next base).
            const double i_R19_inB3 = G_eff_C5R19_ * (V_sw - V_B3) - s_C5_eff;
            const double i_C13      = G_eff_C13R12_ * (V_C3 - V_B2) - s_C13_eff;

            // Cap currents (p -> n direction): i = G*(Vp - Vn) - s.
            const double i_C1  = G_C1_  * (V_inc1 - V_B4) - s_C1_;   // in_c1 -> q4b
            const double i_C15 = G_C15_ * (V_inc1 - v_in) - s_C15_;  // in_c1 -> in
            const double i_C4  = G_C4_  * (V_C4  - V_st)  - s_C4_;   // q4c -> sus_top
            const double i_C14 = G_C14_ * (V_C2  - V_tin) - s_C14_;  // q2c -> tin
            const double i_C9  = G_C9_  * (V_tin - V_tt)  - s_C9_;   // tin -> tone_top
            const double i_C3  = G_C3_  * (V_volw- V_B1)  - s_C3_;   // vol_w -> q1b
            const double i_C2  = G_C2_  * (V_C1  - V_out) - s_C2_;   // q1c -> out

            // ---- Input booster Q4 ----
            // q4b: C1(from in_c1), R14(GND), R9+C10(to q4c), base current.
            out[0] = -i_C1 + V_B4 / R14 + (V_B4 - V_C4) / R9
                   + G_C10_ * (V_B4 - V_C4) - s_C10_ + IB4;
            out[1] = (V_C4 - vcc) / R13 + (V_C4 - V_B4) / R9
                   + G_C10_ * (V_C4 - V_B4) + s_C10_ + IC4 + i_C4;
            out[2] = V_E4 / R22 - IE4;

            // ---- Input node in_c1: R2(to in), C15(to in), C1(to q4b) ----
            out[20] = (V_inc1 - v_in) * G_R2_ + i_C15 + i_C1;

            // ---- Clipper 1 Q3 + diode d1 ----
            out[3] = -i_R19_inB3 + V_B3 / R20 + (V_B3 - V_C3) / R17
                   + G_C12_ * (V_B3 - V_C3) - s_C12_
                   + G_C6_  * (V_B3 - V_d1) - s_C6_ + IB3;
            out[4] = (V_C3 - vcc) / Rc1 + (V_C3 - V_B3) / R17
                   + G_C12_ * (V_C3 - V_B3) + s_C12_ + IC3 - i_d1 + i_C13;
            out[5] = V_E3 / R21 - IE3;
            out[6] = G_C6_ * (V_d1 - V_B3) + s_C6_ + i_d1;

            // ---- Clipper 2 Q2 + diode d2 ----
            out[7] = -i_C13 + V_B2 / R16 + (V_B2 - V_C2) / R15
                   + G_C11_ * (V_B2 - V_C2) - s_C11_
                   + G_C7_  * (V_B2 - V_d2) - s_C7_ + IB2;
            // q2c: collector load, R15+C11 fb, diode, AND C14 to tin (delta 3).
            out[8] = (V_C2 - vcc) / Rc2 + (V_C2 - V_B2) / R15
                   + G_C11_ * (V_C2 - V_B2) + s_C11_ + IC2 - i_d2 + i_C14;
            out[9] = V_E2 / R10 - IE2;
            out[10] = G_C7_ * (V_d2 - V_B2) + s_C7_ + i_d2;

            // ---- tin node: C14(from q2c), C9(to tone_top), R8(to tone_bot) ---
            out[21] = -i_C14 + i_C9 + (V_tin - V_tb) * G_R8_;

            // ---- Tone stack (fed from tin) ----
            out[11] = -i_C9 + V_tt * G_R5_ + (V_tt - V_tw) * G_Rtt_;
            out[12] = (V_tb - V_tin) * G_R8_
                    + (G_C8_ * V_tb - s_C8_) + (V_tb - V_tw) * G_Rtb_;
            // tone_w: R_tt, R_tb, AND R_vol_top to vol_w (delta 4).
            out[13] = (V_tw - V_tt) * G_Rtt_ + (V_tw - V_tb) * G_Rtb_
                    + (V_tw - V_volw) * G_Rvt_;

            // ---- Volume divider (relocated): tone_w -Rvt- vol_w -Rvb- GND ----
            // vol_w: R_vol_top(to tone_w), R_vol_bot(to GND), C3(to Q1 base).
            out[22] = (V_volw - V_tw) * G_Rvt_ + V_volw * G_Rvb_ + i_C3;

            // ---- Output booster Q1 (input from vol_w) ----
            out[14] = -i_C3 + (V_B1 - vcc) * G_R7_ + V_B1 * G_R3_ + IB1;
            out[15] = (V_C1 - vcc) * G_R6_ + i_C2 + IC1;
            out[16] = V_E1 * G_R4_ - IE1;

            // ---- Output node: C2(from q1c), R28(to GND) ----
            out[23] = -i_C2 + V_out * G_R28_;

            // ---- Sustain pot: sus_top has R27 to GND (delta 2) ----
            // C4 current ENTERS sus_top from q4c, so -i_C4 here (cf. q4c: +i_C4).
            out[17] = -i_C4 + (V_st - V_sw) * G_Rst_ + V_st * G_R27_;
            out[18] = (V_sw - V_st) * G_Rst_ + (V_sw - V_sb) * G_Rsb_
                    + (G_eff_C5R19_ * (V_sw - V_B3) - s_C5_eff);
            out[19] = (V_sb - V_sw) * G_Rsb_ + V_sb * G_R23_;
        };

        // Resistor conductances (reciprocals) for the analytic Jacobian.
        const double iR14=1.0/R14, iR9=1.0/R9, iR13=1.0/R13, iR22=1.0/R22;
        const double iR20=1.0/R20, iR17=1.0/R17, iRc1=1.0/Rc1, iR21=1.0/R21;
        const double iR16=1.0/R16, iR15=1.0/R15, iRc2=1.0/Rc2, iR10=1.0/R10;

        // Analytic Jacobian. Mirrors F term-by-term. Transistors use the FULL
        // Ebers-Moll partials (forward + base-collector reverse junction, Br=1),
        // matching the saturation model in F. Diodes contribute gd=2Is*cosh/Vt.
        // Everything else is constant conductance. Validated entry-wise against
        // a central-difference Jacobian.
        auto J = [&](const std::array<double, N>& x,
                     std::array<std::array<double, N>, N>& out) noexcept
        {
            const double V_B4=x[0],V_C4=x[1],V_E4=x[2];
            const double V_B3=x[3],V_C3=x[4],V_E3=x[5], V_d1=x[6];
            const double V_B2=x[7],V_C2=x[8],V_E2=x[9], V_d2=x[10];
            const double V_B1=x[14],V_C1=x[15],V_E1=x[16];

            // Full-EM partials: gB*=dIC/d*, ib*=dIB/d*, ie*=dIE/d*=gB*+ib*.
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
            double gBB4,gBC4,gBE4,ibB4,ibC4,ibE4; bjt(q4,V_B4,V_C4,V_E4,gBB4,gBC4,gBE4,ibB4,ibC4,ibE4);
            double gBB3,gBC3,gBE3,ibB3,ibC3,ibE3; bjt(q3,V_B3,V_C3,V_E3,gBB3,gBC3,gBE3,ibB3,ibC3,ibE3);
            double gBB2,gBC2,gBE2,ibB2,ibC2,ibE2; bjt(q2,V_B2,V_C2,V_E2,gBB2,gBC2,gBE2,ibB2,ibC2,ibE2);
            double gBB1,gBC1,gBE1,ibB1,ibC1,ibE1; bjt(q1,V_B1,V_C1,V_E1,gBB1,gBC1,gBE1,ibB1,ibC1,ibE1);
            const double ieB4=gBB4+ibB4, ieC4=gBC4+ibC4, ieE4=gBE4+ibE4;
            const double ieB3=gBB3+ibB3, ieC3=gBC3+ibC3, ieE3=gBE3+ibE3;
            const double ieB2=gBB2+ibB2, ieC2=gBC2+ibC2, ieE2=gBE2+ibE2;
            const double ieB1=gBB1+ibB1, ieC1=gBC1+ibC1, ieE1=gBE1+ibE1;
            const double gd1 = twoIs_d1_ * std::cosh((V_d1-V_C3)/Vt_d1_) / Vt_d1_;
            const double gd2 = twoIs_d2_ * std::cosh((V_d2-V_C2)/Vt_d2_) / Vt_d2_;

            for (auto& row : out) row.fill(0.0);

            // Q4 input booster (+ in_c1, C4 to sus_top)
            out[0][0]=G_C1_+iR14+iR9+G_C10_+ibB4; out[0][1]=-iR9-G_C10_+ibC4; out[0][2]=ibE4; out[0][20]=-G_C1_;
            out[1][0]=-iR9-G_C10_+gBB4; out[1][1]=iR13+iR9+G_C10_+gBC4+G_C4_; out[1][2]=gBE4; out[1][17]=-G_C4_;
            out[2][0]=-ieB4; out[2][1]=-ieC4; out[2][2]=iR22-ieE4;
            // in_c1: R2(to in) + C15(to in) + C1(to q4b)
            out[20][20]=G_R2_+G_C15_+G_C1_; out[20][0]=-G_C1_;
            // Q3 clipper 1 + diode 1
            out[3][3]=G_eff_C5R19_+iR20+iR17+G_C12_+G_C6_+ibB3; out[3][4]=-iR17-G_C12_+ibC3;
            out[3][5]=ibE3; out[3][6]=-G_C6_; out[3][18]=-G_eff_C5R19_;
            out[4][3]=-iR17-G_C12_+gBB3; out[4][4]=iRc1+iR17+G_C12_+gBC3+gd1+G_eff_C13R12_;
            out[4][5]=gBE3; out[4][6]=-gd1; out[4][7]=-G_eff_C13R12_;
            out[5][3]=-ieB3; out[5][4]=-ieC3; out[5][5]=iR21-ieE3;
            out[6][3]=-G_C6_; out[6][4]=-gd1; out[6][6]=G_C6_+gd1;
            // Q2 clipper 2 + diode 2 (+ C14 to tin)
            out[7][4]=-G_eff_C13R12_; out[7][7]=G_eff_C13R12_+iR16+iR15+G_C11_+G_C7_+ibB2;
            out[7][8]=-iR15-G_C11_+ibC2; out[7][9]=ibE2; out[7][10]=-G_C7_;
            out[8][7]=-iR15-G_C11_+gBB2; out[8][8]=iRc2+iR15+G_C11_+gBC2+gd2+G_C14_;
            out[8][9]=gBE2; out[8][10]=-gd2; out[8][21]=-G_C14_;
            out[9][7]=-ieB2; out[9][8]=-ieC2; out[9][9]=iR10-ieE2;
            out[10][7]=-G_C7_; out[10][8]=-gd2; out[10][10]=G_C7_+gd2;
            // tin (post-C14): C14(from q2c) + C9(to tt) + R8(to tb)
            out[21][8]=-G_C14_; out[21][21]=G_C14_+G_C9_+G_R8_; out[21][11]=-G_C9_; out[21][12]=-G_R8_;
            // Tone stack (fed from tin)
            out[11][21]=-G_C9_; out[11][11]=G_C9_+G_R5_+G_Rtt_; out[11][13]=-G_Rtt_;
            out[12][12]=G_R8_+G_C8_+G_Rtb_; out[12][21]=-G_R8_; out[12][13]=-G_Rtb_;
            out[13][11]=-G_Rtt_; out[13][12]=-G_Rtb_; out[13][13]=G_Rtt_+G_Rtb_+G_Rvt_; out[13][22]=-G_Rvt_;
            // Volume divider (relocated): tone_w -Rvt- vol_w -Rvb- GND; C3 to Q1
            out[22][13]=-G_Rvt_; out[22][22]=G_Rvt_+G_Rvb_+G_C3_; out[22][14]=-G_C3_;
            // Q1 output booster (input from vol_w, output C2 to `out`)
            out[14][22]=-G_C3_; out[14][14]=G_C3_+G_R7_+G_R3_+ibB1; out[14][15]=ibC1; out[14][16]=ibE1;
            out[15][14]=gBB1; out[15][15]=G_R6_+G_C2_+gBC1; out[15][16]=gBE1; out[15][23]=-G_C2_;
            out[16][14]=-ieB1; out[16][15]=-ieC1; out[16][16]=G_R4_-ieE1;
            // out node: C2(from q1c) + R28(to GND)
            out[23][15]=-G_C2_; out[23][23]=G_C2_+G_R28_;
            // Sustain pot (sus_top has R27 to GND)
            out[17][1]=-G_C4_; out[17][17]=G_C4_+G_Rst_+G_R27_; out[17][18]=-G_Rst_;
            out[18][17]=-G_Rst_; out[18][18]=G_Rst_+G_Rsb_+G_eff_C5R19_; out[18][19]=-G_Rsb_; out[18][3]=-G_eff_C5R19_;
            out[19][18]=-G_Rsb_; out[19][19]=G_Rsb_+G_R23_;
        };

        // ---- Schur-complement partition (see file header) --------------------
        // The 24 unknowns split into a NONLINEAR set (4 BJTs * B/C/E + 2 diode
        // nodes = 14) and a LINEAR set (tone/sustain/volume + in_c1/tin/vol_w/out
        // = 10). The linear equations are affine with CONSTANT coefficients, so
        // their 10x10 sub-Jacobian A_LL and the coupling A_LN are fixed until a
        // pot moves. We factor A_LL once (on knob change), precompute the
        // coupling M = A_LL^{-1} A_LN, and Newton-iterate only the 14 nonlinear
        // unknowns. The elimination is EXACT (an algebraic rearrangement of the
        // same 24x24 Newton step), so the iterates -- hence the audio -- match
        // the full solve to rounding. Both A_LL and A_LN are read straight out
        // of the already-validated analytic Jacobian J, so no circuit equation
        // is hand-transcribed here.
        static constexpr std::array<std::size_t, 14> NI =
            { 0,1,2,3,4,5,6,7,8,9,10,14,15,16 };
        static constexpr std::array<std::size_t, 10> LI =
            { 11,12,13,17,18,19,20,21,22,23 };

        if (schur_dirty_) {
            std::array<std::array<double, N>, N> Jf{};
            J(x_, Jf);                                  // linear rows are constant
            std::array<std::array<double, 10>, 10> ALL{};
            for (std::size_t k = 0; k < 10; ++k)
                for (std::size_t l = 0; l < 10; ++l) ALL[k][l] = Jf[LI[k]][LI[l]];
            lu_ALL_ = lu_factor<10>(ALL);
            for (std::size_t i = 0; i < 14; ++i)        // J_NL (constant coupling)
                for (std::size_t k = 0; k < 10; ++k) JNL_[i][k] = Jf[NI[i]][LI[k]];
            for (std::size_t j = 0; j < 14; ++j) {      // M = A_LL^{-1} A_LN
                std::array<double, 10> col{};
                for (std::size_t k = 0; k < 10; ++k) col[k] = Jf[LI[k]][NI[j]];
                lu_solve<10>(lu_ALL_, col);
                for (std::size_t k = 0; k < 10; ++k) M_[k][j] = col[k];
            }
            for (std::size_t i = 0; i < 14; ++i)        // C = J_NL * M (constant)
                for (std::size_t j = 0; j < 14; ++j) {
                    double s = 0.0;
                    for (std::size_t k = 0; k < 10; ++k) s += JNL_[i][k] * M_[k][j];
                    C_[i][j] = s;
                }
            schur_dirty_ = false;
        }

        std::array<double, N> Fv{};
        std::array<std::array<double, N>, N> Jc{};
        double maxF = 0.0;
        for (int it = 0; it < 40; ++it) {
            ++dbg_iters_;
            F(x_, Fv);
            maxF = 0.0;
            for (std::size_t i = 0; i < N; ++i)
                maxF = std::max(maxF, std::abs(Fv[i]));
            if (maxF < nr_tol_ && it >= 1) break;
            J(x_, Jc);
            // g = A_LL^{-1} F_L ; reduced rhs = -F_N + J_NL g ;
            // J_red = J_NN - C_ (C_ = J_NL M precomputed, constant until a knob moves).
            std::array<double, 10> g{};
            for (std::size_t k = 0; k < 10; ++k) g[k] = Fv[LI[k]];
            lu_solve<10>(lu_ALL_, g);
            std::array<std::array<double, 14>, 14> Jr{};
            std::array<double, 14> rhs{};
            for (std::size_t i = 0; i < 14; ++i) {
                const auto& Ji = Jc[NI[i]];
                double jg = 0.0;
                for (std::size_t k = 0; k < 10; ++k) jg += JNL_[i][k] * g[k];
                rhs[i] = -Fv[NI[i]] + jg;
                for (std::size_t j = 0; j < 14; ++j) Jr[i][j] = Ji[NI[j]] - C_[i][j];
            }
            // Fixed-order sparse LU (generated from the declared PM-3 pattern;
            // no pivot search / row swaps). Same reduced system as the dense
            // solve; validated against SPICE.
            if (!gen::solve_pm3_schur(Jr, rhs)) break;            // rhs <- dx_N
            // dx_L = -g - M dx_N, using the UNCLAMPED dx_N so the combined step
            // is identical to the full 24x24 Newton step before clamping.
            for (std::size_t i = 0; i < 14; ++i) {
                double d = rhs[i];
                if (d >  nr_step_) d =  nr_step_;
                if (d < -nr_step_) d = -nr_step_;
                x_[NI[i]] += d;
            }
            for (std::size_t k = 0; k < 10; ++k) {
                double dl = -g[k];
                for (std::size_t j = 0; j < 14; ++j) dl -= M_[k][j] * rhs[j];
                if (dl >  nr_step_) dl =  nr_step_;
                if (dl < -nr_step_) dl = -nr_step_;
                x_[LI[k]] += dl;
            }
        }
        if (maxF > dbg_max_resid_) dbg_max_resid_ = maxF;
        if (maxF > 1e-6) ++dbg_noconv_;
        // Every unknown, and recover to the DC rest point: one index cannot
        // stand in for the vector, and rewriting the unknowns alone leaves the
        // capacitor states poisoned (see src/finite_guard.h).
        if (!tocake::all_finite(x_)) { reset(); ++dbg_nonfinite_; }

        const double V_B4 = x_[0],  V_C4 = x_[1];
        const double V_B3 = x_[3],  V_C3 = x_[4],  V_d1 = x_[6];
        const double V_B2 = x_[7],  V_C2 = x_[8],  V_d2 = x_[10];
        const double V_tt = x_[11], V_tb = x_[12];
        const double V_B1 = x_[14], V_C1 = x_[15];
        const double V_st = x_[17], V_sw = x_[18];
        const double V_inc1 = x_[20], V_tin = x_[21];
        const double V_volw = x_[22], V_out = x_[23];

        // Series-path internal nodes (for cap-state updates).
        const double V_C5_node = (G_C5_ * V_sw + G_R19_ * V_B3 - s_C5_)
                               / (G_C5_ + G_R19_);
        const double V_x = (G_C13_ * V_C3 + G_R12_ * V_B2 - s_C13_)
                         / G_C13R12_denom_;

        // Update bilinear cap states: s = 2*G*(Vp - Vn) - s.
        s_C15_ = 2.0 * G_C15_ * (V_inc1 - v_in)  - s_C15_;
        s_C1_  = 2.0 * G_C1_  * (V_inc1 - V_B4)   - s_C1_;
        s_C10_ = 2.0 * G_C10_ * (V_B4   - V_C4)   - s_C10_;
        s_C4_  = 2.0 * G_C4_  * (V_C4   - V_st)   - s_C4_;
        s_C5_  = 2.0 * G_C5_  * (V_sw   - V_C5_node) - s_C5_;
        s_C12_ = 2.0 * G_C12_ * (V_B3   - V_C3)   - s_C12_;
        s_C6_  = 2.0 * G_C6_  * (V_B3   - V_d1)   - s_C6_;
        s_C13_ = 2.0 * G_C13_ * (V_C3   - V_x)    - s_C13_;
        s_C11_ = 2.0 * G_C11_ * (V_B2   - V_C2)   - s_C11_;
        s_C7_  = 2.0 * G_C7_  * (V_B2   - V_d2)   - s_C7_;
        s_C14_ = 2.0 * G_C14_ * (V_C2   - V_tin)  - s_C14_;
        s_C9_  = 2.0 * G_C9_  * (V_tin  - V_tt)   - s_C9_;
        s_C8_  = 2.0 * G_C8_  * (V_tb   - 0.0)    - s_C8_;
        s_C3_  = 2.0 * G_C3_  * (V_volw - V_B1)   - s_C3_;
        s_C2_  = 2.0 * G_C2_  * (V_C1   - V_out)  - s_C2_;

        switch (tap) {
            case Tap::IbOut: return static_cast<float>(V_st);
            case Tap::C1Out: return static_cast<float>(V_x);
            case Tap::C2Out: return static_cast<float>(V_C2 - bias_cl2_.V_C);
            case Tap::Out:   return static_cast<float>(V_out);
        }
        return static_cast<float>(V_out);
    }

    void reset_iterate() noexcept {
        x_ = { bias_ib_.V_B,  bias_ib_.V_C,  bias_ib_.V_E,
               bias_cl1_.V_B, bias_cl1_.V_C, bias_cl1_.V_E, bias_cl1_.V_C,
               bias_cl2_.V_B, bias_cl2_.V_C, bias_cl2_.V_E, bias_cl2_.V_C,
               V_tt0_, V_tb0_, V_tw0_,
               bias_ob_.V_B,  bias_ob_.V_C,  bias_ob_.V_E,
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    }

    const CakeVariant& v_;
    double p_tone_, p_vol_, p_sus_;
    static constexpr double kDcSeedRate = 1.0e-9;   // Hz: companion conductances ~0
    static constexpr int    kDcSeedMaxPasses = 400;
    static constexpr double kDcSeedTol = 1.0e-6;    // V: above the Newton's own floor (a diode node
                                                    // flickers ~2e-10 V pass to pass for ever)
    bool dc_seed_settled_ = false;
    double fs_ = 48000.0;

    // PM-3-only component values (not in the shared struct). Defaults match the
    // OmegaCake.xlsx "V1 Pyramid Master PM-3" tab; set via set_pm3_parts().
    double pm3_C15_ = 3300e-12, pm3_C14_ = 1e-6, pm3_R27_ = 100e3, pm3_R28_ = 100e3;

public:
    void set_pm3_parts(double C15, double C14, double R27, double R28) noexcept {
        pm3_C15_ = C15; pm3_C14_ = C14; pm3_R27_ = R27; pm3_R28_ = R28;
    }

private:
    dc_bias_result bias_ib_ {}, bias_cl1_ {}, bias_cl2_ {}, bias_ob_ {};

    double G_C1_=0,G_C10_=0,G_C4_=0,G_C5_=0,G_C12_=0,G_C6_=0,G_C13_=0;
    double G_C11_=0,G_C7_=0,G_C9_=0,G_C8_=0,G_C3_=0,G_C2_=0,G_C15_=0,G_C14_=0;
    double G_R2_=0.0, G_R27_=0.0, G_R28_=0.0;
    double G_Rst_=0.0, G_Rsb_=0.0, G_R23_=0.0;
    double G_R19_=0.0, G_eff_C5R19_=0.0, s_C5_denom_=1.0;
    double G_R12_=0.0, G_eff_C13R12_=0.0, G_C13R12_denom_=1.0;
    double G_R5_=0.0, G_R8_=0.0, G_Rtt_=0.0, G_Rtb_=0.0, R_tt_=1.0, R_tb_=1.0;
    double G_R7_=0.0, G_R3_=0.0, G_R6_=0.0, G_R4_=0.0;
    double G_Rvt_=0.0, G_Rvb_=0.0;
    double twoIs_d1_=0.0, Vt_d1_=0.0, twoIs_d2_=0.0, Vt_d2_=0.0;
    double V_tt0_=0.0, V_tb0_=0.0, V_tw0_=0.0;

    double s_C1_=0,s_C10_=0,s_C4_=0,s_C5_=0,s_C12_=0,s_C6_=0,s_C13_=0;
    double s_C11_=0,s_C7_=0,s_C9_=0,s_C8_=0,s_C3_=0,s_C2_=0,s_C15_=0,s_C14_=0;

    std::array<double, N> x_ {};
    double dbg_max_resid_ = 0.0; long dbg_noconv_ = 0;

    // Schur-complement state, all rebuilt only when schur_dirty_ is set (a pot
    // changed the linear-block conductances):
    //   lu_ALL_ : LU of the constant 10x10 linear block A_LL.
    //   M_      : A_LL^{-1} A_LN          (10x14) -- recovers dx_L from dx_N.
    //   JNL_    : nonlinear->linear coupling A_NL (14x10), all CONSTANT conductances.
    //   C_      : J_NL * A_LL^{-1} A_LN   (14x14) -- the constant Schur correction,
    //             so each Newton step only subtracts C_ from the live J_NN block.
    bool                                   schur_dirty_ = true;
    LUFactors<10>                          lu_ALL_ {};
    std::array<std::array<double, 14>, 10> M_   {};
    std::array<std::array<double, 10>, 14> JNL_ {};
    std::array<std::array<double, 14>, 14> C_   {};
    // Newton convergence policy. tol 1e-7 matches the rest of the family
    // (FullPedalV3): ~1e-4 of the signal currents, ~140 dB below the output --
    // inaudible. A 1e-9 tol costs ~8.8 iterations/sample (against ~2.6) for
    // no audible benefit. Step clamp nr_step_. (Tunable via set_nr_policy
    // for experiments.)
    double nr_tol_  = 1e-7;
    double nr_step_ = 0.15;
public:
    void set_nr_policy(double tol, double step) noexcept { nr_tol_ = tol; nr_step_ = step; }
    long  dbg_iters_ = 0;            // total Newton iterations across all samples
    double dbg_max_jac_rel_ = 0.0;   // worst analytic-vs-numerical J rel error
    // Times the guard fired. Exposed so a harness can SEE it: a guard
    // nothing can read is an alarm with no wire.
    long dbg_nonfinite_ = 0;
public:
    long nonfinite_resets() const noexcept { return dbg_nonfinite_; }
private:
};

// ---- the engine that ships -------------------------------------------------
// PyramidMaster is the DK description of this circuit (src/dk_circuits.h);
// PyramidMasterNode above is the reference it is checked against, forced to converge.
using PyramidMaster = shipped_engine_t<DkPyramidMaster>;   // the DK, or its nodal twin (engine_select.h)

} // namespace tocake
