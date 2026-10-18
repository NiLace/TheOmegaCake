// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// Bespoke joint solver for the "V2 Coloured Goat" -- an early-1970s board of the circuit
// with TWO topology departures from the V3 family that the shared FullPedalV3
// cannot represent, so (like the PM-3 and the op-amp) it gets its own solver:
//
//   1. Clipper 1 (Q3) has NO clipping diodes -- it is a plain common-emitter
//      gain stage (R+C feedback only). Only clipper 2 (Q2) clips. Every other
//      transistor variant in the lineup has two diode-clipping stages.
//   2. C14 (0.1 uF) AC-couples Q2's collector to the tone stack, where the V3
//      family DC-couples it directly. This adds a real circuit node (tin) and
//      a ~48 Hz 1st-order high-pass into the tone control, and floats the whole
//      tone/volume section at 0 V DC.
//
// Everything else is the standard four-stage core (input booster Q4, sustain pot,
// clipper 2 Q2, passive tone stack, output booster Q1, volume divider) and is
// lifted verbatim from FullPedalV3. The C14/tin handling is lifted from the
// PM-3 (PyramidMaster), which already AC-couples its tone stack via a C14.
//
// The two departures are read from the SHARED CakeVariant struct, no extra
// fields: clipper_1.diodes with count_fwd == count_rev == 0 means "no diode"
// (Is is also 0), and clipper_2.C_output_couple holds C14 (the V3 sets it 0).
//
// Unknowns x[0..21] (N = 22, same size as FullPedalV3: we DROP clipper-1's
// diode node and ADD the tone-input node):
//   0: V_B(Q4)   1: V_C(Q4)   2: V_E(Q4)        input booster
//   3: V_B(Q3)   4: V_C(Q3)   5: V_E(Q3)        clipper 1 (CLEAN -- no diode)
//   6: V_B(Q2)   7: V_C(Q2)   8: V_E(Q2)        clipper 2
//   9: V_d2                                      clipper-2 diode node
//  10: V_tin                                     tone-stack input (post-C14)
//  11: V_tt  12: V_tb  13: V_tw                  passive tone stack
//  14: V_B(Q1)  15: V_C(Q1)  16: V_E(Q1)         output booster
//  17: V_vt  18: V_out                           volume divider (out = wiper)
//  19: V_st  20: V_sw  21: V_sb                  sustain pot
//
// Full Ebers-Moll with the base-collector junction (Br = 1), so saturation
// matches SPICE. Analytic Jacobian, validated
// entry-wise against central differences. The
// per-sample solve is a 22x22 Newton step over a fixed-order sparse LU (the
// linear block is not Schur-eliminated as in PM-3/FullPedalV3). It is nulled
// against a SPICE deck of the same circuit.
//
// Caps with bilinear (trapezoidal) state (13): C1 C10 C4 C5 C12 C13 C11 C7 C14
//   C9 C8 C3 C2. Series C5->R19 and C13->R12 are eliminated internally (Norton
//   form, rooted at sus_w and q3c) exactly as in FullPedalV3.

#pragma once

#include "dk_circuits.h"
#include "engine_select.h"

#include "newton_arbiter.h"

#include "cake_variant.h"
#include "bjt_dc_bias.h"
#include "nr_solver.h"
#include "solvers_generated.h"   // fixed-order sparse LU (gen::solve_coloured_goat)

#include <algorithm>
#include <array>
#include <cmath>
#include "finite_guard.h"

namespace tocake {

class ColouredGoatNode {
public:
    using T = double;
    static constexpr std::size_t N = 22;
    enum class Tap { C1Out, C2Out, Out };

    explicit ColouredGoatNode(const CakeVariant& v,
                          double p_tone = 0.5,
                          double p_vol  = 0.5,
                          double p_sus  = 1.0) noexcept
        : v_(v), p_tone_(p_tone), p_vol_(p_vol), p_sus_(p_sus) {}

    void prepare(double sample_rate) {
        // Rest seed from the WHOLE circuit, not from each stage on its own.
        // The per-stage DC solves ignore how the stages load each other, so a
        // seed built from them sits several millivolts away from the engine's
        // rest and every run begins with a start-up transient the circuit does
        // not have. At a vanishing sample rate every companion conductance 2*C*fs
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
        bias_cl2_.V_B = x[6];  bias_cl2_.V_C = x[7];  bias_cl2_.V_E = x[8];
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

        // Per-stage DC bias. Clipper 1's (absent) diodes do not affect its DC
        // operating point. Clipper 2's collector is AC-coupled to the tone stack
        // via C14, so at DC it carries NO tone-stack load (unlike the V3, which
        // is direct-coupled) -- R_collector_to_gnd = 0, matching the SPICE .op.
        bias_ib_  = solve_dc_bias(ib.q,  bias_network_for(ib),  v_.vcc, bjt_law::GummelPoon);
        bias_cl1_ = solve_dc_bias(cl1.q, bias_network_for(cl1), v_.vcc, bjt_law::GummelPoon);
        bias_cl2_ = solve_dc_bias(cl2.q, bias_network_for(cl2), v_.vcc, bjt_law::GummelPoon);
        bias_ob_  = solve_dc_bias(ob.q,  bias_network_for(ob),  v_.vcc, bjt_law::GummelPoon);

        G_C1_  = 2.0 * ib.C_input_couple   * fs_;
        G_C10_ = 2.0 * ib.C_feedback_hf    * fs_;
        G_C4_  = 2.0 * ib.C_output_couple  * fs_;
        G_C5_  = 2.0 * sus.C_wiper_couple  * fs_;
        G_C12_ = 2.0 * cl1.C_feedback_hf   * fs_;
        G_C13_ = 2.0 * cl1.C_output_couple * fs_;
        G_C11_ = 2.0 * cl2.C_feedback_hf   * fs_;
        G_C7_  = 2.0 * cl2.C_diode_couple  * fs_;
        G_C14_ = 2.0 * cl2.C_output_couple * fs_;   // C14: clipper-2 -> tone stack
        G_C9_  = 2.0 * tone.C_hp_series    * fs_;
        G_C8_  = 2.0 * tone.C_lp_shunt     * fs_;
        G_C3_  = 2.0 * ob.C_input_couple   * fs_;
        G_C2_  = 2.0 * ob.C_output_couple  * fs_;

        // Input network: in -[R2]- (C1) - q4b, C1 eliminated (Norton) as in V3.
        const double R2 = ib.R_input_series;
        k_in_     = 1.0 + R2 * G_C1_;
        G_eff_in_ = G_C1_ / k_in_;

        // Sustain helpers + the eliminated C5->R19 series into Q3's base.
        G_R23_ = 1.0 / sus.R_pot_bottom;
        G_R19_ = 1.0 / cl1.R_input_series;
        const double G_C5R19_denom = G_C5_ + G_R19_;
        G_eff_C5R19_ = G_C5_ * G_R19_ / G_C5R19_denom;
        s_C5_denom_  = G_C5R19_denom;

        // Eliminated C13->R12 series into Q2's base.
        G_R12_          = 1.0 / cl2.R_input_series;
        G_C13R12_denom_ = G_C13_ + G_R12_;
        G_eff_C13R12_   = G_C13_ * G_R12_ / G_C13R12_denom_;

        G_R5_  = 1.0 / tone.R_hp_shunt;
        G_R8_  = 1.0 / tone.R_lp_series;

        G_R7_ = 1.0 / ob.R_bias_top;
        G_R3_ = 1.0 / ob.R_base_to_gnd;
        G_R6_ = 1.0 / ob.R_collector_load;
        G_R4_ = 1.0 / ob.R_emitter;

        // Only clipper 2 has diodes here. (Clipper 1's diode params are zero.)
        twoIs_d2_ = 2.0 * cl2.diodes.Is;  Vt_d2_ = cl2.diodes.N * V_T_300K;

        iR14_ = 1.0 / ib.R_base_to_gnd;    iR9_  = 1.0 / ib.R_feedback;
        iR13_ = 1.0 / ib.R_collector_load; iR22_ = 1.0 / ib.R_emitter;
        iR20_ = 1.0 / cl1.R_base_to_gnd;   iR17_ = 1.0 / cl1.R_feedback;
        iRc1_ = 1.0 / cl1.R_collector_load; iR21_ = 1.0 / cl1.R_emitter;
        iR16_ = 1.0 / cl2.R_base_to_gnd;   iR15_ = 1.0 / cl2.R_feedback;
        iRc2_ = 1.0 / cl2.R_collector_load; iR10_ = 1.0 / cl2.R_emitter;
        iBf4_ = 1.0 / ib.q.Bf;  iBf3_ = 1.0 / cl1.q.Bf;
        iBf2_ = 1.0 / cl2.q.Bf; iBf1_ = 1.0 / ob.q.Bf;

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
        // C14 AC-couples the tone stack, so the whole tone/volume/sustain
        // section sits at 0 V DC (as in PM-3). Transistor biases are unchanged.
        V_tt0_ = 0.0; V_tb0_ = 0.0; V_tw0_ = 0.0;

        x_ = { bias_ib_.V_B,  bias_ib_.V_C,  bias_ib_.V_E,
               bias_cl1_.V_B, bias_cl1_.V_C, bias_cl1_.V_E,
               bias_cl2_.V_B, bias_cl2_.V_C, bias_cl2_.V_E, bias_cl2_.V_C,
               0.0,                                  // tin: C14-blocked -> 0 V
               V_tt0_, V_tb0_, V_tw0_,
               bias_ob_.V_B,  bias_ob_.V_C,  bias_ob_.V_E,
               0.0, 0.0,                             // vol_top, out
               0.0, 0.0, 0.0 };                      // sustain (AC-coupled)

        s_C1_  = G_C1_  * (0.0           - bias_ib_.V_B);
        s_C10_ = G_C10_ * (bias_ib_.V_B  - bias_ib_.V_C);
        s_C4_  = G_C4_  * (bias_ib_.V_C  - 0.0);          // sus_top = 0 at DC
        s_C5_  = G_C5_  * (0.0           - bias_cl1_.V_B); // sus_w = 0 at DC
        s_C12_ = G_C12_ * (bias_cl1_.V_B - bias_cl1_.V_C);
        s_C13_ = G_C13_ * (bias_cl1_.V_C - bias_cl2_.V_B);
        s_C11_ = G_C11_ * (bias_cl2_.V_B - bias_cl2_.V_C);
        s_C7_  = G_C7_  * (bias_cl2_.V_B - bias_cl2_.V_C);
        s_C14_ = G_C14_ * (bias_cl2_.V_C - 0.0);          // tin = 0 at DC
        s_C9_  = G_C9_  * (0.0           - V_tt0_);        // tin = 0 at DC
        s_C8_  = G_C8_  * (V_tb0_        - 0.0);
        s_C3_  = G_C3_  * (V_tw0_        - bias_ob_.V_B);  // tone_w = 0 at DC
        s_C2_  = G_C2_  * (bias_ob_.V_C  - 0.0);

        dbg_iters_ = dbg_calls_ = dbg_nonconv_ = 0;
        dbg_max_resid_ = 0.0;
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

    long   dbg_iters()     const noexcept { return dbg_iters_;     }
    long   dbg_calls()     const noexcept { return dbg_calls_;     }
    long   dbg_nonconv()   const noexcept { return dbg_nonconv_;   }
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

    float run_sample(double v_in, Tap tap) noexcept {
        const double s_in_eff  = s_C1_  / k_in_;
        const double s_C5_eff  = s_C5_  * G_R19_ / s_C5_denom_;
        const double s_C13_eff = s_C13_ * G_R12_ / G_C13R12_denom_;

        const auto&  q4 = v_.input_booster.q;
        const auto&  q3 = v_.clipper_1.q;
        const auto&  q2 = v_.clipper_2.q;
        const auto&  q1 = v_.output_booster.q;

        const double V_T = V_T_300K;
        const double vcc = v_.vcc;

        auto F = [&](const std::array<double, N>& x,
                     std::array<double, N>&       out) noexcept
        {
            const double V_B4 = x[0],  V_C4 = x[1],  V_E4 = x[2];
            const double V_B3 = x[3],  V_C3 = x[4],  V_E3 = x[5];
            const double V_B2 = x[6],  V_C2 = x[7],  V_E2 = x[8];
            const double V_d2 = x[9];
            const double V_tin = x[10];
            const double V_tt = x[11], V_tb = x[12], V_tw = x[13];
            const double V_B1 = x[14], V_C1 = x[15], V_E1 = x[16];
            const double V_vt = x[17], V_out = x[18];
            const double V_st = x[19], V_sw = x[20], V_sb = x[21];

            // Full Ebers-Moll transport model WITH the base-collector junction
            // (reverse current, Br=1 = SPICE default) so SATURATION matches the
            // SPICE NPN model. The clean high-gain clipper-1 (Q3, no diodes to
            // clamp it) is driven hard into saturation here; the forward-only
            // model (FullPedalV3) would miss that clip shape. Same approach as
            // the PM-3 solver. Forward Early factor as in FullPedalV3.
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

            // Only clipper 2 clips.
            const double i_d2 = twoIs_d2_ * std::sinh((V_d2 - V_C2) / Vt_d2_);

            // Eliminated series paths (Norton form).
            const double i_C4_out   = G_C4_ * (V_C4 - V_st) - s_C4_;
            const double i_R19_inB3 = G_eff_C5R19_ * (V_sw - V_B3) - s_C5_eff;
            const double i_C13      = G_eff_C13R12_ * (V_C3 - V_B2) - s_C13_eff;
            // C14 (clipper-2 collector -> tin) and C9 (tin -> tone_top).
            const double i_C14 = G_C14_ * (V_C2  - V_tin) - s_C14_;
            const double i_C9  = G_C9_  * (V_tin - V_tt)  - s_C9_;

            // ---- Input booster Q4 ----
            out[0] = G_eff_in_ * (V_B4 - v_in) + s_in_eff
                   + V_B4 * iR14_ + (V_B4 - V_C4) * iR9_
                   + G_C10_ * (V_B4 - V_C4) - s_C10_ + IB4;
            out[1] = (V_C4 - vcc) * iR13_ + (V_C4 - V_B4) * iR9_
                   + G_C10_ * (V_C4 - V_B4) + s_C10_ + IC4 + i_C4_out;
            out[2] = V_E4 * iR22_ - IE4;

            // ---- Clipper 1 Q3 (CLEAN -- no diode, no C6) ----
            out[3] = -i_R19_inB3 + V_B3 * iR20_ + (V_B3 - V_C3) * iR17_
                   + G_C12_ * (V_B3 - V_C3) - s_C12_ + IB3;
            out[4] = (V_C3 - vcc) * iRc1_ + (V_C3 - V_B3) * iR17_
                   + G_C12_ * (V_C3 - V_B3) + s_C12_ + IC3 + i_C13;
            out[5] = V_E3 * iR21_ - IE3;

            // ---- Clipper 2 Q2 + diode d2 ----
            out[6] = -i_C13 + V_B2 * iR16_ + (V_B2 - V_C2) * iR15_
                   + G_C11_ * (V_B2 - V_C2) - s_C11_
                   + G_C7_  * (V_B2 - V_d2) - s_C7_ + IB2;
            out[7] = (V_C2 - vcc) * iRc2_ + (V_C2 - V_B2) * iR15_
                   + G_C11_ * (V_C2 - V_B2) + s_C11_ + IC2 - i_d2 + i_C14;
            out[8] = V_E2 * iR10_ - IE2;
            out[9] = G_C7_ * (V_d2 - V_B2) + s_C7_ + i_d2;

            // ---- tin node: C14(from q2c), C9(to tone_top), R8(to tone_bot) ---
            out[10] = -i_C14 + i_C9 + (V_tin - V_tb) * G_R8_;

            // ---- Tone stack (fed from tin) ----
            out[11] = -i_C9 + V_tt * G_R5_ + (V_tt - V_tw) * G_Rtt_;
            out[12] = (V_tb - V_tin) * G_R8_
                    + (G_C8_ * V_tb - s_C8_) + (V_tb - V_tw) * G_Rtb_;
            out[13] = (V_tw - V_tt) * G_Rtt_ + (V_tw - V_tb) * G_Rtb_
                    + (G_C3_ * (V_tw - V_B1) - s_C3_);

            // ---- Output booster Q1 ----
            out[14] = (G_C3_ * (V_B1 - V_tw) + s_C3_)
                    + (V_B1 - vcc) * G_R7_ + V_B1 * G_R3_ + IB1;
            out[15] = (V_C1 - vcc) * G_R6_
                    + (G_C2_ * (V_C1 - V_vt) - s_C2_) + IC1;
            out[16] = V_E1 * G_R4_ - IE1;

            // ---- Volume divider ----
            out[17] = (G_C2_ * (V_vt - V_C1) + s_C2_)
                    + (V_vt - V_out) * G_Rvt_;
            out[18] = (V_out - V_vt) * G_Rvt_ + V_out * G_Rvb_;

            // ---- Sustain pot ----
            out[19] = (G_C4_ * (V_st - V_C4) + s_C4_)
                    + (V_st - V_sw) * G_Rst_;
            out[20] = (V_sw - V_st) * G_Rst_ + (V_sw - V_sb) * G_Rsb_
                    + (G_eff_C5R19_ * (V_sw - V_B3) - s_C5_eff);
            out[21] = (V_sb - V_sw) * G_Rsb_ + V_sb * G_R23_;
        };

        // Analytic Jacobian (forward-only EM partials, as in FullPedalV3).
        auto J = [&](const std::array<double, N>& x,
                     std::array<std::array<double, N>, N>& out) noexcept
        {
            const double V_B4=x[0],V_C4=x[1],V_E4=x[2];
            const double V_B3=x[3],V_C3=x[4],V_E3=x[5];
            const double V_B2=x[6],V_C2=x[7],V_E2=x[8], V_d2=x[9];
            const double V_B1=x[14],V_C1=x[15],V_E1=x[16];

            // Full-EM partials matching F: gB*=dIC/d*, ib*=dIB/d*, ie*=dIE/d*.
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
            const double gd2 = twoIs_d2_ * std::cosh((V_d2 - V_C2) / Vt_d2_) / Vt_d2_;

            for (auto& row : out) row.fill(0.0);

            // Q4 (input booster)
            out[0][0]=G_eff_in_+iR14_+iR9_+G_C10_+ibB4; out[0][1]=-iR9_-G_C10_+ibC4; out[0][2]=ibE4;
            out[1][0]=-iR9_-G_C10_+gBB4; out[1][1]=iR13_+iR9_+G_C10_+gBC4+G_C4_; out[1][2]=gBE4; out[1][19]=-G_C4_;
            out[2][0]=-ieB4; out[2][1]=-ieC4; out[2][2]=iR22_-ieE4;
            // Q3 (clipper 1 -- clean, no diode node)
            out[3][3]=G_eff_C5R19_+iR20_+iR17_+G_C12_+ibB3; out[3][4]=-iR17_-G_C12_+ibC3;
            out[3][5]=ibE3; out[3][20]=-G_eff_C5R19_;
            out[4][3]=-iR17_-G_C12_+gBB3; out[4][4]=iRc1_+iR17_+G_C12_+gBC3+G_eff_C13R12_;
            out[4][5]=gBE3; out[4][6]=-G_eff_C13R12_;
            out[5][3]=-ieB3; out[5][4]=-ieC3; out[5][5]=iR21_-ieE3;
            // Q2 (clipper 2) + diode 2
            out[6][4]=-G_eff_C13R12_; out[6][6]=G_eff_C13R12_+iR16_+iR15_+G_C11_+G_C7_+ibB2;
            out[6][7]=-iR15_-G_C11_+ibC2; out[6][8]=ibE2; out[6][9]=-G_C7_;
            out[7][6]=-iR15_-G_C11_+gBB2; out[7][7]=iRc2_+iR15_+G_C11_+gBC2+gd2+G_C14_;
            out[7][8]=gBE2; out[7][9]=-gd2; out[7][10]=-G_C14_;
            out[8][6]=-ieB2; out[8][7]=-ieC2; out[8][8]=iR10_-ieE2;
            out[9][6]=-G_C7_; out[9][7]=-gd2; out[9][9]=G_C7_+gd2;
            // tin (post-C14): C14(from q2c) + C9(to tt) + R8(to tb)
            out[10][7]=-G_C14_; out[10][10]=G_C14_+G_C9_+G_R8_; out[10][11]=-G_C9_; out[10][12]=-G_R8_;
            // Tone stack (fed from tin)
            out[11][10]=-G_C9_; out[11][11]=G_C9_+G_R5_+G_Rtt_; out[11][13]=-G_Rtt_;
            out[12][10]=-G_R8_; out[12][12]=G_R8_+G_C8_+G_Rtb_; out[12][13]=-G_Rtb_;
            out[13][11]=-G_Rtt_; out[13][12]=-G_Rtb_; out[13][13]=G_Rtt_+G_Rtb_+G_C3_; out[13][14]=-G_C3_;
            // Q1 (output booster)
            out[14][13]=-G_C3_; out[14][14]=G_C3_+G_R7_+G_R3_+ibB1; out[14][15]=ibC1; out[14][16]=ibE1;
            out[15][14]=gBB1; out[15][15]=G_R6_+G_C2_+gBC1; out[15][16]=gBE1; out[15][17]=-G_C2_;
            out[16][14]=-ieB1; out[16][15]=-ieC1; out[16][16]=G_R4_-ieE1;
            // Volume divider
            out[17][15]=-G_C2_; out[17][17]=G_C2_+G_Rvt_; out[17][18]=-G_Rvt_;
            out[18][17]=-G_Rvt_; out[18][18]=G_Rvt_+G_Rvb_;
            // Sustain pot
            out[19][1]=-G_C4_; out[19][19]=G_C4_+G_Rst_; out[19][20]=-G_Rst_;
            out[20][3]=-G_eff_C5R19_; out[20][19]=-G_Rst_; out[20][20]=G_Rst_+G_Rsb_+G_eff_C5R19_; out[20][21]=-G_Rsb_;
            out[21][20]=-G_Rsb_; out[21][21]=G_Rsb_+G_R23_;
        };

        // Newton over the generated fixed-order sparse LU. Same convergence policy as
        // FullPedalV3: tol 1e-7, step-clamp 0.15, iteration cap 32. Warm-started from
        // the previous sample; typical 3-5 iterations.
        std::array<double, N> Fv{}, bvec{};
        std::array<std::array<double, N>, N> Jc{};
        ++dbg_calls_;
        bool   conv  = false;
        double lastF = 0.0;
        int    it    = 0;
        for (; it < 32; ++it) {
            F(x_, Fv);
            double maxF = 0.0;
            for (std::size_t i = 0; i < N; ++i)
                maxF = std::max(maxF, std::abs(Fv[i]));
            lastF = maxF;
            if (maxF < 1e-7 && it >= 1) { conv = true; break; }
            J(x_, Jc);
            for (std::size_t i = 0; i < N; ++i) bvec[i] = -Fv[i];
            if (!gen::solve_coloured_goat(Jc, bvec)) break;   // fixed-order sparse LU
            for (std::size_t i = 0; i < N; ++i) {
                double d = bvec[i];
                if (d >  0.15) d =  0.15;
                if (d < -0.15) d = -0.15;
                x_[i] += d;
            }
        }
        dbg_iters_ += it;
        if (!conv) ++dbg_nonconv_;
        if (lastF > dbg_max_resid_) dbg_max_resid_ = lastF;
        // Every unknown, and recover to the DC rest point: one index cannot
        // stand in for the vector, and rewriting the unknowns alone leaves the
        // capacitor states poisoned (see src/finite_guard.h).
        if (!tocake::all_finite(x_)) { reset(); ++dbg_nonfinite_; }

        const double V_B3 = x_[3], V_C3 = x_[4];
        const double V_B2 = x_[6], V_C2 = x_[7];
        const double V_B4 = x_[0], V_C4 = x_[1];
        const double V_tin = x_[10];
        const double V_C1 = x_[15], V_vt = x_[17], V_out = x_[18];
        const double V_sw = x_[20];

        const double V_C5_node = (G_C5_ * V_sw + G_R19_ * V_B3 - s_C5_)
                               / (G_C5_ + G_R19_);
        const double V_x = (G_C13_ * V_C3 + G_R12_ * V_B2 - s_C13_)
                         / G_C13R12_denom_;

        const double R2 = v_.input_booster.R_input_series;
        s_C1_  = 2.0 * G_C1_  * ((v_in - V_B4 + R2 * s_C1_) / k_in_) - s_C1_;
        s_C10_ = 2.0 * G_C10_ * (V_B4 - V_C4)        - s_C10_;
        s_C4_  = 2.0 * G_C4_  * (V_C4 - x_[19])      - s_C4_;
        s_C5_  = 2.0 * G_C5_  * (V_sw - V_C5_node)   - s_C5_;
        s_C12_ = 2.0 * G_C12_ * (V_B3 - V_C3)        - s_C12_;
        s_C13_ = 2.0 * G_C13_ * (V_C3 - V_x)         - s_C13_;
        s_C11_ = 2.0 * G_C11_ * (V_B2 - V_C2)        - s_C11_;
        s_C7_  = 2.0 * G_C7_  * (V_B2 - x_[9])       - s_C7_;
        s_C14_ = 2.0 * G_C14_ * (V_C2 - V_tin)       - s_C14_;
        s_C9_  = 2.0 * G_C9_  * (V_tin - x_[11])     - s_C9_;
        s_C8_  = 2.0 * G_C8_  * (x_[12] - 0.0)       - s_C8_;
        s_C3_  = 2.0 * G_C3_  * (x_[13] - x_[14])    - s_C3_;
        s_C2_  = 2.0 * G_C2_  * (V_C1 - V_vt)        - s_C2_;

        switch (tap) {
            case Tap::C1Out: return static_cast<float>(V_x);
            case Tap::C2Out: return static_cast<float>(V_C2 - bias_cl2_.V_C);
            case Tap::Out:   return static_cast<float>(V_out);
        }
        return static_cast<float>(V_out);
    }

    void reset_iterate() noexcept {
        x_ = { bias_ib_.V_B,  bias_ib_.V_C,  bias_ib_.V_E,
               bias_cl1_.V_B, bias_cl1_.V_C, bias_cl1_.V_E,
               bias_cl2_.V_B, bias_cl2_.V_C, bias_cl2_.V_E, bias_cl2_.V_C,
               0.0, V_tt0_, V_tb0_, V_tw0_,
               bias_ob_.V_B,  bias_ob_.V_C,  bias_ob_.V_E,
               0.0, 0.0, 0.0, 0.0, 0.0 };
    }

    const CakeVariant& v_;
    double p_tone_, p_vol_, p_sus_;
    static constexpr double kDcSeedRate = 1.0e-9;   // Hz: companion conductances ~0
    static constexpr int    kDcSeedMaxPasses = 400;
    static constexpr double kDcSeedTol = 1.0e-6;    // V: above the Newton's own floor (a diode node
                                                    // flickers ~2e-10 V pass to pass for ever)
    bool dc_seed_settled_ = false;
    double fs_ = 48000.0;

    dc_bias_result bias_ib_ {}, bias_cl1_ {}, bias_cl2_ {}, bias_ob_ {};

    double G_C1_=0,G_C10_=0,G_C4_=0,G_C5_=0,G_C12_=0,G_C13_=0;
    double G_C11_=0,G_C7_=0,G_C14_=0,G_C9_=0,G_C8_=0,G_C3_=0,G_C2_=0;
    double k_in_=1.0, G_eff_in_=0.0;
    double G_Rst_=0.0, G_Rsb_=0.0, G_R23_=0.0;
    double G_R19_=0.0, G_eff_C5R19_=0.0, s_C5_denom_=1.0;
    double G_R12_=0.0, G_eff_C13R12_=0.0, G_C13R12_denom_=1.0;
    double G_R5_=0.0, G_R8_=0.0, G_Rtt_=0.0, G_Rtb_=0.0;
    double G_R7_=0.0, G_R3_=0.0, G_R6_=0.0, G_R4_=0.0;
    double G_Rvt_=0.0, G_Rvb_=0.0;
    double twoIs_d2_=0.0, Vt_d2_=0.0;
    double V_tt0_=0.0, V_tb0_=0.0, V_tw0_=0.0;

    double iR14_=0,iR9_=0,iR13_=0,iR22_=0;
    double iR20_=0,iR17_=0,iRc1_=0,iR21_=0;
    double iR16_=0,iR15_=0,iRc2_=0,iR10_=0;
    double iBf4_=0,iBf3_=0,iBf2_=0,iBf1_=0;

    double dbg_max_jac_rel_ = 0.0;
    long   dbg_iters_ = 0, dbg_calls_ = 0, dbg_nonconv_ = 0;
    double dbg_max_resid_ = 0.0;

    double s_C1_=0,s_C10_=0,s_C4_=0,s_C5_=0,s_C12_=0,s_C13_=0;
    double s_C11_=0,s_C7_=0,s_C14_=0,s_C9_=0,s_C8_=0,s_C3_=0,s_C2_=0;

    std::array<double, N> x_ {};
    // Times the guard fired. Exposed so a harness can SEE it: a guard
    // nothing can read is an alarm with no wire.
    long dbg_nonfinite_ = 0;
public:
    long nonfinite_resets() const noexcept { return dbg_nonfinite_; }
private:
};

// ---- the engine that ships -------------------------------------------------
// ColouredGoat is the DK description of this circuit (src/dk_circuits.h);
// ColouredGoatNode above is the reference it is checked against, forced to converge.
using ColouredGoat = shipped_engine_t<DkColouredGoat>;   // the DK, or its nodal twin (engine_select.h)

} // namespace tocake
