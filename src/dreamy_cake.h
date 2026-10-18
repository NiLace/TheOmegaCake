// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// Bespoke joint solver for the "Vm Dreamy Cake" -- a clone of the
// original circuit with ONE topology departure the shared FullPedalV3 cannot do, so
// (like the Jumbo Open Goat) it gets its own solver:
//
//   * The emitter resistors of the FIRST THREE stages are REMOVED -- the input
//     booster (Q4) and BOTH clipper emitters (Q3, Q2) jumper STRAIGHT to ground
//     (R22 = R21 = R10 = 0 ohm). FullPedalV3 divides by R_emitter, so a 0-ohm
//     emitter is a divide-by-zero there. Here those three emitters are pinned to
//     0 V (the device sees Vbe = Vb) and their emitter current flows to ground,
//     untracked -- the three emitter nodes become decoupled dummies. Only the
//     OUTPUT booster (Q1) keeps a real emitter resistor (R4) + node.
//
// Everything else is the standard direct-coupled circuit (V3-family): both
// clippers carry a 1N914 diode pair, the tone stack is DC-fed from Q2's
// collector (C9 treble cap + R8 bass resistor, no C14), and the output booster
// + volume divider are stock. Lifted verbatim from FullPedalV3, with the BJTs
// upgraded to the FULL Ebers-Moll model (BC junction, Br = 1) the grounded
// max-gain stages need to match SPICE saturation (same model the Coloured Goat /
// Distortion & Cake use). NOTE on devices: schematic-correct mix -- Q4/Q3/Q2 =
// 2N5089, Q1 (output) = 2N5088 (the schematic overrides the netlist's
// all-2N5089). D5 is reverse-polarity protection, reverse-biased
// in use, and is left out. R27 (100 R in series with the 9 V) and C14 (330 uF)
// set the rail: the four stages draw 1.51 mA through R27, so it rests at
// 8.8487 V. This node engine holds the rail STATIC at that value (the variant's
// vcc); the shipped DK description carries R27 + C14 as a rail node
// (DreamyCircuitRail in dk_circuits.h), and the SPICE deck carries both.
//
// The tone pot is "wired in reverse" on the real pedal; that is a KNOB-DIRECTION
// (UX) detail handled in lv2.cpp by inverting the tone fraction, NOT here -- this
// solver and its SPICE deck use the same straight p_tone mapping so the null is
// exact at every tone setting.
//
// Unknowns x[0..21] (N = 22, FullPedalV3 layout):
//   0: V_B(Q4)  1: V_C(Q4)  2: V_E(Q4)=0 (dummy)   input booster (grounded E)
//   3: V_B(Q3)  4: V_C(Q3)  5: V_E(Q3)=0 (dummy)   clipper 1 (grounded E)
//   6: V_d1                                         clipper-1 diode node
//   7: V_B(Q2)  8: V_C(Q2)  9: V_E(Q2)=0 (dummy)   clipper 2 (grounded E)
//  10: V_d2                                         clipper-2 diode node
//  11: V_tt  12: V_tb  13: V_tw                     passive tone stack
//  14: V_B(Q1) 15: V_C(Q1) 16: V_E(Q1)             output booster (real R4)
//  17: V_vt  18: V_out                              volume divider (out = wiper)
//  19: V_st  20: V_sw  21: V_sb                     sustain pot
//
// Analytic Jacobian, validated entry-wise vs central differences. Per sample, a
// 22-unknown Newton over the generated sparse LU (gen::solve_dreamy_cake). It is
// nulled against a SPICE deck of the same circuit.

#pragma once

#include "dk_circuits.h"
#include "engine_select.h"

#include "newton_arbiter.h"

#include "cake_variant.h"
#include "bjt_dc_bias.h"
#include "nr_solver.h"
#include "solvers_generated.h"   // fixed-order sparse LU (gen::solve_dreamy_cake)

#include <algorithm>
#include <array>
#include <cmath>
#include "finite_guard.h"

namespace tocake {

class DreamyCakeNode {
public:
    using T = double;
    static constexpr std::size_t N = 22;
    enum class Tap { C1Out, C2Out, Out };

    explicit DreamyCakeNode(const CakeVariant& v,
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
        bias_ib_.V_B = x[0];  bias_ib_.V_C = x[1];
        bias_cl1_.V_B = x[3];  bias_cl1_.V_C = x[4];
        bias_cl2_.V_B = x[7];  bias_cl2_.V_C = x[8];
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

        // Per-stage DC bias. The grounded emitters (R_emitter = 0) make the
        // bias solver pin V_E = 0 for Q4/Q3/Q2 (it special-cases R_E = 0). The
        // tone stack is DC-fed from Q2's collector (direct, R8), so clipper 2's
        // bias carries that load (same as FullPedalV3).
        bias_ib_  = solve_dc_bias(ib.q,  bias_network_for(ib),  v_.vcc, bjt_law::GummelPoon);
        bias_cl1_ = solve_dc_bias(cl1.q, bias_network_for(cl1), v_.vcc, bjt_law::GummelPoon);
        bias_cl2_ = solve_dc_bias(cl2.q,
                                  bias_network_for(cl2, tone_stack_dc_load(tone)),
                                  v_.vcc, bjt_law::GummelPoon);
        bias_ob_  = solve_dc_bias(ob.q,  bias_network_for(ob),  v_.vcc, bjt_law::GummelPoon);

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

        const double R2 = ib.R_input_series;
        k_in_     = 1.0 + R2 * G_C1_;
        G_eff_in_ = G_C1_ / k_in_;

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

        // Transistor-stage resistor reciprocals. The first three emitters are
        // grounded (no R), so only the collector-load / base / feedback ones.
        iR14_ = 1.0 / ib.R_base_to_gnd;    iR9_  = 1.0 / ib.R_feedback;
        iR13_ = 1.0 / ib.R_collector_load;
        iR20_ = 1.0 / cl1.R_base_to_gnd;   iR17_ = 1.0 / cl1.R_feedback;
        iRc1_ = 1.0 / cl1.R_collector_load;
        iR16_ = 1.0 / cl2.R_base_to_gnd;   iR15_ = 1.0 / cl2.R_feedback;
        iRc2_ = 1.0 / cl2.R_collector_load;

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
        const auto& tone = v_.tone;

        const double I_tone = bias_cl2_.V_C / tone_stack_dc_load(tone);
        V_tb0_ = bias_cl2_.V_C - I_tone * tone.R_lp_series;
        V_tw0_ = V_tb0_        - I_tone * R_tb_;
        V_tt0_ = V_tw0_        - I_tone * R_tt_;

        x_ = { bias_ib_.V_B,  bias_ib_.V_C,  0.0,
               bias_cl1_.V_B, bias_cl1_.V_C, 0.0, bias_cl1_.V_C,
               bias_cl2_.V_B, bias_cl2_.V_C, 0.0, bias_cl2_.V_C,
               V_tt0_, V_tb0_, V_tw0_,
               bias_ob_.V_B,  bias_ob_.V_C,  bias_ob_.V_E,
               0.0, 0.0,
               0.0, 0.0, 0.0 };

        s_C1_  = G_C1_  * (0.0           - bias_ib_.V_B);
        s_C10_ = G_C10_ * (bias_ib_.V_B  - bias_ib_.V_C);
        s_C4_  = G_C4_  * (bias_ib_.V_C  - 0.0);
        s_C5_  = G_C5_  * (0.0           - bias_cl1_.V_B);
        s_C12_ = G_C12_ * (bias_cl1_.V_B - bias_cl1_.V_C);
        s_C6_  = G_C6_  * (bias_cl1_.V_B - bias_cl1_.V_C);
        s_C13_ = G_C13_ * (bias_cl1_.V_C - bias_cl2_.V_B);
        s_C11_ = G_C11_ * (bias_cl2_.V_B - bias_cl2_.V_C);
        s_C7_  = G_C7_  * (bias_cl2_.V_B - bias_cl2_.V_C);
        s_C9_  = G_C9_  * (bias_cl2_.V_C - V_tt0_);
        s_C8_  = G_C8_  * (V_tb0_        - 0.0);
        s_C3_  = G_C3_  * (V_tw0_        - bias_ob_.V_B);
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
        R_tt_  = R * (1.0 - p_tone_) + 1e-6;
        R_tb_  = R * p_tone_         + 1e-6;
        G_Rtt_ = 1.0 / R_tt_;
        G_Rtb_ = 1.0 / R_tb_;
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
            const double V_d1 = x[6];
            const double V_B2 = x[7],  V_C2 = x[8],  V_E2 = x[9];
            const double V_d2 = x[10];
            const double V_tt = x[11], V_tb = x[12], V_tw = x[13];
            const double V_B1 = x[14], V_C1 = x[15], V_E1 = x[16];
            const double V_vt = x[17], V_out = x[18];
            const double V_st = x[19], V_sw = x[20], V_sb = x[21];

            // Full Ebers-Moll (Br = 1) -- BC junction included so the grounded,
            // max-gain stages saturate like SPICE.
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
            // Q4/Q3/Q2 emitters grounded -> Ve = 0; Q1 has a real emitter node.
            double IC4, IB4, IE4; bjt(q4, V_B4, V_C4, 0.0,  IC4, IB4, IE4);
            double IC3, IB3, IE3; bjt(q3, V_B3, V_C3, 0.0,  IC3, IB3, IE3);
            double IC2, IB2, IE2; bjt(q2, V_B2, V_C2, 0.0,  IC2, IB2, IE2);
            double IC1, IB1, IE1; bjt(q1, V_B1, V_C1, V_E1, IC1, IB1, IE1);
            (void)IE4; (void)IE3; (void)IE2;   // grounded emitter current -> GND

            const double i_d1 = twoIs_d1_ * std::sinh((V_d1 - V_C3) / Vt_d1_);
            const double i_d2 = twoIs_d2_ * std::sinh((V_d2 - V_C2) / Vt_d2_);

            const double i_C4_out   = G_C4_ * (V_C4 - V_st) - s_C4_;
            const double i_R19_inB3 = G_eff_C5R19_ * (V_sw - V_B3) - s_C5_eff;
            const double i_C13      = G_eff_C13R12_ * (V_C3 - V_B2) - s_C13_eff;

            // ---- Input booster Q4 (grounded emitter) ----
            out[0] = G_eff_in_ * (V_B4 - v_in) + s_in_eff
                   + V_B4 * iR14_ + (V_B4 - V_C4) * iR9_
                   + G_C10_ * (V_B4 - V_C4) - s_C10_ + IB4;
            out[1] = (V_C4 - vcc) * iR13_ + (V_C4 - V_B4) * iR9_
                   + G_C10_ * (V_C4 - V_B4) + s_C10_ + IC4 + i_C4_out;
            out[2] = V_E4;                                   // dummy: grounded

            // ---- Clipper 1 Q3 (grounded emitter) + diode d1 ----
            out[3] = -i_R19_inB3 + V_B3 * iR20_ + (V_B3 - V_C3) * iR17_
                   + G_C12_ * (V_B3 - V_C3) - s_C12_
                   + G_C6_  * (V_B3 - V_d1) - s_C6_ + IB3;
            out[4] = (V_C3 - vcc) * iRc1_ + (V_C3 - V_B3) * iR17_
                   + G_C12_ * (V_C3 - V_B3) + s_C12_ + IC3 - i_d1 + i_C13;
            out[5] = V_E3;                                   // dummy: grounded
            out[6] = G_C6_ * (V_d1 - V_B3) + s_C6_ + i_d1;

            // ---- Clipper 2 Q2 (grounded emitter) + diode d2 + tone feed ----
            out[7] = -i_C13 + V_B2 * iR16_ + (V_B2 - V_C2) * iR15_
                   + G_C11_ * (V_B2 - V_C2) - s_C11_
                   + G_C7_  * (V_B2 - V_d2) - s_C7_ + IB2;
            out[8] = (V_C2 - vcc) * iRc2_ + (V_C2 - V_B2) * iR15_
                   + G_C11_ * (V_C2 - V_B2) + s_C11_ + IC2 - i_d2
                   + (G_C9_ * (V_C2 - V_tt) - s_C9_)
                   + G_R8_ * (V_C2 - V_tb);
            out[9] = V_E2;                                   // dummy: grounded
            out[10] = G_C7_ * (V_d2 - V_B2) + s_C7_ + i_d2;

            // ---- Tone stack (direct-coupled from Q2 collector) ----
            out[11] = (G_C9_ * (V_tt - V_C2) + s_C9_)
                    + V_tt * G_R5_ + (V_tt - V_tw) * G_Rtt_;
            out[12] = (V_tb - V_C2) * G_R8_
                    + (G_C8_ * V_tb - s_C8_) + (V_tb - V_tw) * G_Rtb_;
            out[13] = (V_tw - V_tt) * G_Rtt_ + (V_tw - V_tb) * G_Rtb_
                    + (G_C3_ * (V_tw - V_B1) - s_C3_);

            // ---- Output booster Q1 (real emitter R4) ----
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

        auto J = [&](const std::array<double, N>& x,
                     std::array<std::array<double, N>, N>& out) noexcept
        {
            const double V_B4=x[0],V_C4=x[1];
            const double V_B3=x[3],V_C3=x[4], V_d1=x[6];
            const double V_B2=x[7],V_C2=x[8], V_d2=x[10];
            const double V_B1=x[14],V_C1=x[15],V_E1=x[16];

            // Full-EM partials. For the grounded stages we pass Ve = 0 (constant)
            // and use only the Vb/Vc derivatives (the emitter node is decoupled).
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
            double gBB4,gBC4,gBE4,ibB4,ibC4,ibE4; bjt(q4,V_B4,V_C4,0.0, gBB4,gBC4,gBE4,ibB4,ibC4,ibE4);
            double gBB3,gBC3,gBE3,ibB3,ibC3,ibE3; bjt(q3,V_B3,V_C3,0.0, gBB3,gBC3,gBE3,ibB3,ibC3,ibE3);
            double gBB2,gBC2,gBE2,ibB2,ibC2,ibE2; bjt(q2,V_B2,V_C2,0.0, gBB2,gBC2,gBE2,ibB2,ibC2,ibE2);
            double gBB1,gBC1,gBE1,ibB1,ibC1,ibE1; bjt(q1,V_B1,V_C1,V_E1, gBB1,gBC1,gBE1,ibB1,ibC1,ibE1);
            (void)gBE4;(void)ibE4;(void)gBE3;(void)ibE3;(void)gBE2;(void)ibE2;
            const double ieB1=gBB1+ibB1, ieC1=gBC1+ibC1, ieE1=gBE1+ibE1;
            const double gd1 = twoIs_d1_ * std::cosh((V_d1 - V_C3) / Vt_d1_) / Vt_d1_;
            const double gd2 = twoIs_d2_ * std::cosh((V_d2 - V_C2) / Vt_d2_) / Vt_d2_;

            for (auto& row : out) row.fill(0.0);

            // Q4 (grounded emitter)
            out[0][0]=G_eff_in_+iR14_+iR9_+G_C10_+ibB4; out[0][1]=-iR9_-G_C10_+ibC4;
            out[1][0]=-iR9_-G_C10_+gBB4; out[1][1]=iR13_+iR9_+G_C10_+gBC4+G_C4_; out[1][19]=-G_C4_;
            out[2][2]=1.0;
            // Q3 (grounded emitter) + diode 1
            out[3][3]=G_eff_C5R19_+iR20_+iR17_+G_C12_+G_C6_+ibB3; out[3][4]=-iR17_-G_C12_+ibC3;
            out[3][6]=-G_C6_; out[3][20]=-G_eff_C5R19_;
            out[4][3]=-iR17_-G_C12_+gBB3; out[4][4]=iRc1_+iR17_+G_C12_+gBC3+gd1+G_eff_C13R12_;
            out[4][6]=-gd1; out[4][7]=-G_eff_C13R12_;
            out[5][5]=1.0;
            out[6][3]=-G_C6_; out[6][4]=-gd1; out[6][6]=G_C6_+gd1;
            // Q2 (grounded emitter) + diode 2 + tone feed
            out[7][4]=-G_eff_C13R12_; out[7][7]=G_eff_C13R12_+iR16_+iR15_+G_C11_+G_C7_+ibB2;
            out[7][8]=-iR15_-G_C11_+ibC2; out[7][10]=-G_C7_;
            out[8][7]=-iR15_-G_C11_+gBB2; out[8][8]=iRc2_+iR15_+G_C11_+gBC2+gd2+G_C9_+G_R8_;
            out[8][10]=-gd2; out[8][11]=-G_C9_; out[8][12]=-G_R8_;
            out[9][9]=1.0;
            out[10][7]=-G_C7_; out[10][8]=-gd2; out[10][10]=G_C7_+gd2;
            // Tone stack
            out[11][8]=-G_C9_; out[11][11]=G_C9_+G_R5_+G_Rtt_; out[11][13]=-G_Rtt_;
            out[12][8]=-G_R8_; out[12][12]=G_R8_+G_C8_+G_Rtb_; out[12][13]=-G_Rtb_;
            out[13][11]=-G_Rtt_; out[13][12]=-G_Rtb_; out[13][13]=G_Rtt_+G_Rtb_+G_C3_; out[13][14]=-G_C3_;
            // Q1 (output booster, real emitter)
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
            if (!gen::solve_dreamy_cake(Jc, bvec)) break;   // fixed-order sparse LU over the J pattern
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

        const double V_B4 = x_[0], V_C4 = x_[1];
        const double V_B3 = x_[3], V_C3 = x_[4], V_d1 = x_[6];
        const double V_B2 = x_[7], V_C2 = x_[8], V_d2 = x_[10];
        const double V_tt = x_[11], V_tb = x_[12], V_tw = x_[13];
        const double V_B1 = x_[14], V_C1 = x_[15];
        const double V_vt = x_[17], V_out = x_[18];
        const double V_st = x_[19], V_sw = x_[20];

        const double V_C5_node = (G_C5_ * V_sw + G_R19_ * V_B3 - s_C5_)
                               / (G_C5_ + G_R19_);
        const double V_x = (G_C13_ * V_C3 + G_R12_ * V_B2 - s_C13_)
                         / G_C13R12_denom_;

        const double R2 = v_.input_booster.R_input_series;
        s_C1_  = 2.0 * G_C1_  * ((v_in - V_B4 + R2 * s_C1_) / k_in_) - s_C1_;
        s_C10_ = 2.0 * G_C10_ * (V_B4 - V_C4)        - s_C10_;
        s_C4_  = 2.0 * G_C4_  * (V_C4 - V_st)        - s_C4_;
        s_C5_  = 2.0 * G_C5_  * (V_sw - V_C5_node)   - s_C5_;
        s_C12_ = 2.0 * G_C12_ * (V_B3 - V_C3)        - s_C12_;
        s_C6_  = 2.0 * G_C6_  * (V_B3 - V_d1)        - s_C6_;
        s_C13_ = 2.0 * G_C13_ * (V_C3 - V_x)         - s_C13_;
        s_C11_ = 2.0 * G_C11_ * (V_B2 - V_C2)        - s_C11_;
        s_C7_  = 2.0 * G_C7_  * (V_B2 - V_d2)        - s_C7_;
        s_C9_  = 2.0 * G_C9_  * (V_C2 - V_tt)        - s_C9_;
        s_C8_  = 2.0 * G_C8_  * (V_tb - 0.0)         - s_C8_;
        s_C3_  = 2.0 * G_C3_  * (V_tw - V_B1)        - s_C3_;
        s_C2_  = 2.0 * G_C2_  * (V_C1 - V_vt)        - s_C2_;

        switch (tap) {
            case Tap::C1Out: return static_cast<float>(V_x);
            case Tap::C2Out: return static_cast<float>(V_C2 - bias_cl2_.V_C);
            case Tap::Out:   return static_cast<float>(V_out);
        }
        return static_cast<float>(V_out);
    }

    void reset_iterate() noexcept {
        x_ = { bias_ib_.V_B,  bias_ib_.V_C,  0.0,
               bias_cl1_.V_B, bias_cl1_.V_C, 0.0, bias_cl1_.V_C,
               bias_cl2_.V_B, bias_cl2_.V_C, 0.0, bias_cl2_.V_C,
               V_tt0_, V_tb0_, V_tw0_,
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

    double G_C1_=0,G_C10_=0,G_C4_=0,G_C5_=0,G_C12_=0,G_C6_=0,G_C13_=0;
    double G_C11_=0,G_C7_=0,G_C9_=0,G_C8_=0,G_C3_=0,G_C2_=0;
    double k_in_=1.0, G_eff_in_=0.0;
    double G_Rst_=0.0, G_Rsb_=0.0, G_R23_=0.0;
    double G_R19_=0.0, G_eff_C5R19_=0.0, s_C5_denom_=1.0;
    double G_R12_=0.0, G_eff_C13R12_=0.0, G_C13R12_denom_=1.0;
    double G_R5_=0.0, G_R8_=0.0, G_Rtt_=0.0, G_Rtb_=0.0, R_tt_=1.0, R_tb_=1.0;
    double G_R7_=0.0, G_R3_=0.0, G_R6_=0.0, G_R4_=0.0;
    double G_Rvt_=0.0, G_Rvb_=0.0;
    double twoIs_d1_=0.0, Vt_d1_=0.0, twoIs_d2_=0.0, Vt_d2_=0.0;
    double V_tt0_=0.0, V_tb0_=0.0, V_tw0_=0.0;

    double iR14_=0,iR9_=0,iR13_=0;
    double iR20_=0,iR17_=0,iRc1_=0;
    double iR16_=0,iR15_=0,iRc2_=0;

    double dbg_max_jac_rel_ = 0.0;
    long   dbg_iters_ = 0, dbg_calls_ = 0, dbg_nonconv_ = 0;
    double dbg_max_resid_ = 0.0;

    double s_C1_=0,s_C10_=0,s_C4_=0,s_C5_=0,s_C12_=0,s_C6_=0,s_C13_=0;
    double s_C11_=0,s_C7_=0,s_C9_=0,s_C8_=0,s_C3_=0,s_C2_=0;

    std::array<double, N> x_ {};
    // Times the guard fired. Exposed so a harness can SEE it: a guard
    // nothing can read is an alarm with no wire.
    long dbg_nonfinite_ = 0;
public:
    long nonfinite_resets() const noexcept { return dbg_nonfinite_; }
private:
};

// ---- the engine that ships -------------------------------------------------
// DreamyCake is the DK description of this circuit (src/dk_circuits.h);
// DreamyCakeNode above is the reference it is checked against, forced to converge.
using DreamyCake = shipped_engine_t<DkDreamyCake>;   // the DK, or its nodal twin (engine_select.h)

} // namespace tocake
