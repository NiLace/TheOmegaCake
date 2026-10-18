// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// DC bias-point solver for the common-emitter NPN stages used by
// the four-stage transistor fuzz (input booster, clippers, output booster).
//
// Solves the same Ebers-Moll-with-Early-effect equations SPICE uses,
// so the bias points match SPICE's .op output to numerical noise.
// Newton-Raphson with analytical Jacobian, 3D system per stage.
//
// Two bias topologies, expressed in one unified system via conductances:
//   - Feedback-biased  (input booster + clippers): R_base_to_collector
//   - Divider-biased   (output booster):           R_base_to_vcc
// Set the unused resistor to 0; the solver treats 1/0 as zero conductance
// (== no current path) and the equation collapses to the right form.
//
// All currents and voltages in SI units (A, V). Temperature is baked in
// via V_T = kT/q at 300.15 K (LTSpice default), the value DiodeClipper
// uses against SPICE -- DO NOT change this without re-deriving against a
// new SPICE reference.

#pragma once

#include "cake_variant.h"
#include "nr_solver.h"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace tocake {

// Thermal voltage at 300.15 K (27 C, the simulators' TNOM), as the REFERENCE simulator computes it -- measured,
// not recalled: two identical diodes (IS 1e-14) at 1 mA and 10 mA give V = 0.65511880518 / 0.714675060213 V, so
// Vt = 25.8649529 mV from the slope and from the 1 mA point alike. That is k/q = 8.617343e-5 (CODATA 1998): neither
// CODATA 2014 (25.8649170 mV) nor CODATA 2018 (25.8649258). The ONE definition: every engine reads it from here.
inline constexpr double V_T_300K = 25.8649529e-3;

struct bias_network {
    double R_collector_load;       // V_C <-> V+ (R13/R18/R11/R6)
    double R_emitter;              // V_E <-> GND (R22/R21/R10/R4)
    double R_base_to_gnd;          // V_B <-> GND (R14/R20/R16/R3)
    double R_base_to_collector;    // V_B <-> V_C  (R9/R17/R15; 0 for output booster)
    double R_base_to_vcc;          // V_B <-> V+   (R7; 0 for input booster + clippers)
    double R_collector_to_gnd;     // V_C <-> GND  -- extra DC load on collector.
                                   // 0 for stages whose output is fully AC-coupled
                                   // (caps block DC -> no current to GND).
                                   // Nonzero for V3's clipper 2: the tone stack
                                   // is direct-coupled to that collector, so at
                                   // DC the collector sees R8+R_pot+R5 to GND.
};

struct dc_bias_result {
    double V_B  = 0.0;     // base voltage  (V)
    double V_C  = 0.0;     // collector voltage  (V)
    double V_E  = 0.0;     // emitter voltage  (V)
    double V_BE = 0.0;     // = V_B - V_E
    double V_CE = 0.0;     // = V_C - V_E
    double I_C  = 0.0;     // collector current  (A)
    double I_B  = 0.0;     // base current  (A)
    int    iterations = 0;
    bool   converged  = false;
};

namespace detail {

// 3x3 solver is the one in nr_solver.h. Local helper:
inline double conductance(double R) noexcept {
    return (R > 0.0) ? (1.0 / R) : 0.0;
}

} // namespace detail

// Solve the DC operating point for one common-emitter NPN stage.
//
// Unknowns: x = [V_B, V_C, I_C].  V_E and V_BE/V_CE are derived.
//
// Residual equations (must be zero at the operating point):
//
//   F0 (KCL at base):
//       G_FB * (V_C - V_B) + G_top * (V+ - V_B) - G_BG * V_B - I_C/Bf = 0
//
//   F1 (KCL at collector):
//       G_C  * (V+ - V_C) - G_FB * (V_C - V_B) - I_C = 0
//
//   F2 (BJT equation, Ebers-Moll active region with Early effect):
//       Is * (exp(V_BE / V_T) - 1) * (1 + V_CE / V_AF) - I_C = 0
//       where V_BE = V_B - V_E  and  V_CE = V_C - V_E
//             V_E  = I_C * (1 + 1/Bf) * R_E
//
// Conductance form makes both bias topologies one equation set:
//   feedback bias  -> G_top = 0, G_FB > 0
//   divider bias   -> G_top > 0, G_FB = 0
//
// Which device law to solve. This is a PARAMETER and not a compile switch
// because this solver is shared by several engines and each must get the law
// its own residual uses: a rest point from the other law makes the engine
// spend every record drifting from one operating point to the other, and the
// residual piles into the lowest bins -- a settling transient reading as
// model error.
//
//   Legacy      IC = Is(e-1)(1 + Vce/Vaf),  IB = IC/Bf
//   GummelPoon  IC = Is(e-1)(1 - Vbc/Vaf),  IB = Is/Bf (e-1)
//
// The second is what SPICE solves for a card with Br=1, Nr=1, Ikf=inf,
// Ise=0, Var=inf -- which is every transistor card in this repository.
// Against LTspice's .op for cake_v3.cir the two read:
//
//   node   LTspice     Legacy          GummelPoon
//   q3c    4.437241    4.45435 (+17.1 mV)  4.43730 (+0.06 mV)
//   q2c    4.391893    4.40771 (+15.8 mV)  4.39195 (+0.06 mV)
enum class bjt_law { Legacy, GummelPoon };

inline dc_bias_result solve_dc_bias(const bjt_params& q,
                                     const bias_network& net,
                                     double vcc,
                                     bjt_law law = bjt_law::Legacy,
                                     double v_t = V_T_300K,
                                     int max_iter = 50,
                                     double tol = 1e-12) noexcept
{
    const bool gp = (law == bjt_law::GummelPoon);
    using detail::conductance;

    const double G_C    = conductance(net.R_collector_load);
    const double G_E    = conductance(net.R_emitter);
    const double G_BG   = conductance(net.R_base_to_gnd);
    const double G_FB   = conductance(net.R_base_to_collector);
    const double G_top  = conductance(net.R_base_to_vcc);
    const double G_C_gnd = conductance(net.R_collector_to_gnd);   // extra DC load

    // alpha relates V_E to I_C: V_E = alpha * I_C  (KCL at emitter,
    // I_E = I_C * (1 + 1/Bf), V_E = I_E * R_E). Captured as R_E_eff.
    const double R_E_eff = (G_E > 0.0)
        ? (1.0 + 1.0 / q.Bf) / G_E
        : 0.0;     // emitter direct to GND (no R_E)

    // Initial guess. Start V_B from the bias network itself -- the resistor
    // divider voltage for a divider-biased stage (G_top > 0), else ~0.7 V for
    // a feedback-biased stage -- then derive a CONSISTENT I_C from the emitter
    // resistor. Do not seed I_C = Is*exp(0.65/V_T) when there is an emitter
    // resistor: it scales with the device's Is, and a large-Is part (e.g.
    // 2N5087) feeding a large emitter resistor (the 3.3 k output booster)
    // starts with V_E past the rail and Newton diverges. The seed does not move
    // the root (Newton converges to the same F = 0 solution regardless of start).
    double V_B = (G_top > 0.0)
        ? vcc * G_top / (G_top + G_BG + G_FB)   // divider, transistor off
        : 0.7;                                  // feedback-biased
    double V_C = 0.5 * vcc;
    const double V_E_guess = (V_B > 0.65) ? (V_B - 0.65) : 0.05;
    double I_C = (R_E_eff > 0.0)
        ? (V_E_guess / R_E_eff)                 // consistent with the emitter R
        : q.Is * std::exp(0.65 / v_t);          // no R_E: fall back to V_BE form

    dc_bias_result out;
    out.iterations = 0;
    out.converged  = false;

    for (int it = 0; it < max_iter; ++it) {
        const double V_E  = R_E_eff * I_C;
        const double V_BE = V_B - V_E;
        const double V_CE = V_C - V_E;

        const double e   = std::exp(V_BE / v_t);
        const double em1 = e - 1.0;
        // SPICE's forward Early factor is 1/qb with qb = 1/(1 - Vbc/Vaf),
        // and the base current carries none of it: Ibe = Is/Bf*(exp(Vbe/Vt)-1).
        const double early = gp ? (1.0 - (V_B - V_C) / q.Vaf)
                                : (1.0 + V_CE / q.Vaf);
        const double I_B_m = gp ? (q.Is * em1 / q.Bf) : (I_C / q.Bf);
        const double bjt_I = q.Is * em1 * early;

        // Residuals
        const double F[3] = {
            // KCL at base
            G_FB * (V_C - V_B) + G_top * (vcc - V_B) - G_BG * V_B - I_B_m,
            // KCL at collector. G_C_gnd adds an extra current path V_C -> GND
            // (used by V3's clipper 2 where the tone stack is direct-coupled).
            G_C * (vcc - V_C) - G_FB * (V_C - V_B) - I_C - G_C_gnd * V_C,
            // BJT model
            bjt_I - I_C,
        };

        // Convergence: residuals small *and* the next step is small
        const double maxF = std::max({std::abs(F[0]), std::abs(F[1]), std::abs(F[2])});
        if (maxF < tol && it > 0) {
            out.converged  = true;
            out.V_B = V_B; out.V_C = V_C; out.V_E = V_E;
            out.V_BE = V_BE; out.V_CE = V_CE;
            out.I_C = I_C; out.I_B = I_B_m;
            out.iterations = it;
            return out;
        }

        // Jacobian (partials of F w.r.t. [V_B, V_C, I_C])
        //
        // For F0 (KCL_base):
        //   d/dV_B = -G_FB - G_top - G_BG
        //   d/dV_C = G_FB
        //   d/dI_C = -1/Bf
        //
        // For F1 (KCL_col):
        //   d/dV_B = G_FB
        //   d/dV_C = -G_C - G_FB
        //   d/dI_C = -1
        //
        // For F2 (BJT):  I_C_model = Is * (exp(V_BE/V_T) - 1) * (1 + V_CE/V_AF)
        //   V_BE = V_B - V_E,   V_CE = V_C - V_E,   V_E = R_E_eff * I_C
        //   so V_BE depends on V_B and I_C; V_CE depends on V_C and I_C.
        //   d(V_BE)/dV_B =  1,   d(V_BE)/dI_C = -R_E_eff
        //   d(V_CE)/dV_C =  1,   d(V_CE)/dI_C = -R_E_eff
        //   d(bjt_I)/dV_BE = (Is/V_T) * e * early
        //   d(bjt_I)/dV_CE = (Is * em1) / V_AF
        //   so:
        //   d/dV_B = d(bjt_I)/dV_BE * 1
        //   d/dV_C = d(bjt_I)/dV_CE * 1
        //   d/dI_C = d(bjt_I)/dV_BE * (-R_E_eff) + d(bjt_I)/dV_CE * (-R_E_eff) - 1
        const double dIb_dVbe = (q.Is / v_t) * e * early;
        const double dIb_dVce = (q.Is * em1) / q.Vaf;

        // In the Gummel-Poon form F0's base-current term is Is/Bf*(e-1), so
        // it has derivatives in V_B and in I_C (through V_E) and none in V_C;
        // and F2's Early factor depends on V_B and V_C directly rather than
        // on V_CE.
        // V_E is taken as R_E_eff*I_C with R_E_eff = (1+1/Bf)/G_E, i.e.
        // I_E approximated as I_C*(1+1/Bf). In the Gummel-Poon form I_B
        // differs from I_C/Bf by the Early factor, which with Bf = 540 puts
        // the error in V_E at about 1.3e-4 of V_E -- a few microvolts.
        const double dIBm_dVb = (q.Is / v_t) * e / q.Bf;
        const double J[3][3] = {
            { gp ? (-G_FB - G_top - G_BG - dIBm_dVb) : (-G_FB - G_top - G_BG),
              G_FB,
              gp ? (R_E_eff * dIBm_dVb) : (-1.0 / q.Bf) },
            {  G_FB,                  -G_C - G_FB - G_C_gnd,        -1.0        },
            {  gp ? (dIb_dVbe - dIb_dVce) : dIb_dVbe,
               dIb_dVce,
               gp ? (-R_E_eff * dIb_dVbe - 1.0)
                  : (-R_E_eff * (dIb_dVbe + dIb_dVce) - 1.0) },
        };

        // solve J * dx = -F
        double dx[3];
        const double rhs[3] = { -F[0], -F[1], -F[2] };
        if (!detail::solve_3x3(J, rhs, dx)) {
            // Singular Jacobian -- bail
            break;
        }

        // Damped Newton: limit V_BE step to ~5 * V_T to keep the
        // exponential from blowing up on a bad iterate
        constexpr double max_vbe_step = 5.0 * V_T_300K;
        const double dV_BE_proposed = dx[0] - R_E_eff * dx[2];
        if (std::abs(dV_BE_proposed) > max_vbe_step) {
            const double scale = max_vbe_step / std::abs(dV_BE_proposed);
            dx[0] *= scale; dx[1] *= scale; dx[2] *= scale;
        }

        V_B += dx[0];
        V_C += dx[1];
        I_C += dx[2];
        // Clamp I_C to physically-meaningful positive values
        if (I_C < 1e-30) I_C = 1e-30;
    }

    // Failed to converge -- still report the last iterate so debugging
    // is possible. converged stays false.
    const double V_E  = R_E_eff * I_C;
    out.V_B = V_B; out.V_C = V_C; out.V_E = V_E;
    out.V_BE = V_B - V_E; out.V_CE = V_C - V_E;
    out.I_C = I_C; out.I_B = I_C / q.Bf;
    out.iterations = max_iter;
    return out;
}

// Convenience builders that pull the bias-network resistor values from
// the per-stage variant structs. Keeps callers from having to know
// which variant field maps to which solver field.

inline bias_network bias_network_for(const booster_params& p,
                                      double R_collector_to_gnd = 0.0) noexcept
{
    return {
        .R_collector_load    = p.R_collector_load,
        .R_emitter           = p.R_emitter,
        .R_base_to_gnd       = p.R_base_to_gnd,
        .R_base_to_collector = p.R_feedback,
        .R_base_to_vcc       = p.R_bias_top,
        .R_collector_to_gnd  = R_collector_to_gnd,
    };
}

inline bias_network bias_network_for(const clipper_params& p,
                                      double R_collector_to_gnd = 0.0) noexcept
{
    return {
        .R_collector_load    = p.R_collector_load,
        .R_emitter           = p.R_emitter,
        .R_base_to_gnd       = p.R_base_to_gnd,
        .R_base_to_collector = p.R_feedback,
        .R_base_to_vcc       = 0.0,
        .R_collector_to_gnd  = R_collector_to_gnd,
    };
}

// DC resistance from the V3-style direct-coupled clipper-2 collector to
// GND through the tone stack. The tone stack's caps (C_hp_series,
// C_lp_shunt) are open at DC; the only DC path is:
//     collector -> R_lp_series -> pot_bot -> full pot -> pot_top
//                                          -> R_hp_shunt -> GND
inline double tone_stack_dc_load(const tone_params& t) noexcept {
    return t.R_lp_series + t.R_pot + t.R_hp_shunt;
}

// Small-signal hybrid-pi parameters derived from a Q-point.
struct hybrid_pi {
    double g_m;     // transconductance (S)
    double r_pi;    // base-emitter intrinsic resistance (ohm)
    double r_o;     // collector-emitter Early-effect resistance (ohm)
};

inline hybrid_pi hybrid_pi_from(const dc_bias_result& q,
                                 const bjt_params& bjt,
                                 double v_t = V_T_300K) noexcept
{
    const double g_m  = q.I_C / v_t;
    const double r_pi = bjt.Bf / g_m;
    // r_o = (V_AF + V_CE) / I_C; we use V_AF + V_CE for accuracy. At
    // V_CE ~ 4 V, V_AF = 100, this is ~ 4% more than the simpler V_AF/I_C.
    const double r_o  = (bjt.Vaf + q.V_CE) / q.I_C;
    return { g_m, r_pi, r_o };
}

} // namespace tocake
