// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// ===========================================================================
//  The Op Cake -- the op-amp ("IC") version of the circuit -- real-time DSP model
// ===========================================================================
//
//  WHAT THIS IS
//  ------------
//  A physically-modelled emulation of the late-1970s op-amp version of the
//  circuit. Unlike the four-transistor original (modelled separately by
//  OmegaCake / cake_*.h), this circuit is built from op-amps:
//
//      input --[C1]--> U1A  --> U1B  --> U2 (LM741) --> tone stack --> Volume
//                   inverting   Sustain    diode-clamp    mid-scoop     out
//                    booster     stage      clipper         tone
//
//      * U1A  (one half of a JRC4558): inverting amp, gain ~ -R3/R2 ~ -5.9.
//      * U1B  (the other half): non-inverting "Sustain" stage, gain ~ x10.
//               The Sustain pot sits INSIDE this stage's feedback leg and also
//               taps off to the next stage, so the two roles share one network.
//      * U2   (an LM741): the clipper. Its feedback loop holds the diode clamp
//               -- THREE 1N4448 in series each direction, anti-parallel (the
//               original arrangement). That high, soft ~2 V clamp is the
//               "open/vocal" voice of the op-amp version.
//      * Tone: the classic passive mid-scoop, plus a Tone-Bypass
//               switch that routes a flat post-clip tap around the scoop.
//      * Single 9 V supply; U1 biased at +4.5 V, U2 at +4.945 V.
//
//  Component values live in op_cake_variant.h (transcribed from the KiCad
//  netlist of the circuit -- the source of truth for this variant).
//
//  MODELLING PHILOSOPHY -- why the circuit splits into independent blocks
//  ---------------------------------------------------------------------
//  Every op-amp output is low-impedance and there is no feedback path that
//  wraps from one stage back into an earlier one. So each stage hands a clean
//  voltage to the next and can be solved on its own, in order. (The transistor
//  circuit is the opposite -- its passive RC coupling makes every stage load
//  its neighbours -- which is why THAT model needs one big joint solve and
//  this one does not.)
//
//  Each block is solved in ABSOLUTE volts (not just the AC component) so that
//  the op-amp bias references and the supply-rail clipping are exact. Output
//  taps subtract the known DC bias to recover the AC signal.
//
//  THREE SOLVE TECHNIQUES, AND WHY EACH IS USED
//  --------------------------------------------
//  1. LINEAR networks (the passive tone stack, and an op-amp running in its
//     normal high-gain region) -> a direct matrix solve. The system matrix
//     only changes when a knob moves, so it is LU-factorised once per knob
//     change (see refactor_*()) and only back-substituted per sample
//     (lu_solve). Capacitors enter as trapezoidal companion sources (CapTrap).
//
//  2. OP-AMP CLIPPING (U1A/U1B hitting a supply rail; U2's LM741 hitting its
//     slew-rate limit) -> REGION DETECTION, never a stiff Newton nonlinearity.
//     Solve assuming the op-amp is in its linear region; if the result runs
//     past the limit, the op-amp is clipping -- pin it to the limit and re-
//     solve. This is exact and converges in one extra solve. (A hard limiter
//     fed straight into Newton has ZERO gradient once saturated, so Newton
//     ping-pongs between rails and never finds the real root -- on U1B's rails
//     and on the LM741 slew alike.)
//
//  3. The DIODE CLAMP at U2 -> a small 1-D Newton solve (the diodes are a
//     genuine smooth nonlinearity). Three identical series diodes collapse
//     EXACTLY into one element i = 2*Is*sinh(v / (3*N*Vt)), so there is just
//     one unknown to iterate.
//
//  VALIDATION
//  ----------
//  Every block is null-tested against the matching tap of a SPICE deck of the
//  same circuit (idealised the same way on both sides, so the null is honest).
//  Representative nulls at realistic input (0.1-0.5 V): U1A ~-100 dB,
//  U1B ~-77 dB, U2 clip ~-50 dB. The full output is held to a -60 dB null
//  against the circuit.
//
//  KNOWN APPROXIMATIONS (small; documented so they are not mistaken for bugs)
//  -------------------------------------------------------------------------
//   * The U1B->U2 boundary is SOLVED, not approximated:
//     U1B is linear in U2's inverting input su, so U1B is solved at su = bias
//     plus its response to su, U2's Newton carries that slope, and U1B settles
//     at the su U2 finds. Do not pin su at the bias:
//     the drift is not negligible where the LM741 slews and the sustain pot
//     gives nC5 a source impedance -- about 20 dB of null against the deck at
//     sus 0.3 / tone 1.0 / 0.1 V. At full sustain the wiper is driven by
//     U1B's output directly and the coupling is ~3e-5.
//   * U1A/U1B are solved as IDEAL (infinite-gain) op-amps; the SPICE deck uses
//     a very high finite gain (1e6) and, for U1B's rails, a smooth tanh rather
//     than this model's hard clamp. The knee mismatch is sub-millivolt and the
//     nulls confirm it is immaterial.
//
//  REAL-TIME SAFETY
//  ----------------
//  No heap allocation, no locks, no system calls on the audio path. All solver
//  storage is fixed-size std::array on the stack; every iteration count is
//  bounded. prepare()/set_*() (which refactorise matrices) run at block rate,
//  never per sample.
// ===========================================================================

#pragma once

#include "isa_hot.h"

#include <cstdint>

#include "op_cake_variant.h"
#include "bjt_dc_bias.h"
#include "nr_solver.h"
#include "finite_guard.h"

#include <array>
#include <cmath>

// Iteration cap of the U2 clipper's 1-D Newton. 60 ships; the macro exists so
// a census can ask whether a solve that runs out converges with more room.
#ifndef TOCAKE_OP_SU_MAXIT
#define TOCAKE_OP_SU_MAXIT 60
#endif

namespace tocake {

// ---------------------------------------------------------------------------
// Trapezoidal companion model of a capacitor.
//
// A capacitor between nodes a and b is replaced, each sample, by a resistor of
// conductance Gc = 2C/Ts in parallel with a current source `ieq` that carries
// the cap's history. The current flowing a->b is then  i = Gc*(Va-Vb) - ieq,
// which is exactly the trapezoidal-rule discretisation of i = C dv/dt. After
// the node voltages are solved, commit() updates the history for next sample.
//
//   set_dc(v): seed the history at a DC operating point (v = steady Va-Vb) so
//   the model starts already biased -- no cold-start transient (this mirrors
//   SPICE beginning its .tran from the solved .op point).
// ---------------------------------------------------------------------------
struct CapTrap {
    double Gc = 0.0;             // companion conductance, 2C/Ts
    double vp = 0.0, ip = 0.0;   // previous voltage (Va-Vb) and previous a->b current

    void   prepare(double C, double Ts) noexcept { Gc = 2.0 * C / Ts; }
    void   set_dc(double v)             noexcept { vp = v; ip = 0.0; }
    double ieq()                  const noexcept { return Gc * vp + ip; }
    void   commit(double v_new)         noexcept {
        ip = Gc * v_new - ieq();   // a->b current this step (uses the OLD ieq)
        vp = v_new;
    }
};

// Output tap, matching the SPICE deck's per-stage .wave files. Each returns the
// AC voltage at that node (the deck's tap subtracts the same DC bias). Only
// `Out` is used in production; the others exist for stage-by-stage validation.
enum class OpTap { U1A, U1B, Clip, Out };

class OpCake {
public:
    explicit OpCake(const OpCakeParams& p) noexcept : p_(p) {}

    // ----- lifecycle -------------------------------------------------------

    void prepare(double sample_rate) {
        fs_ = sample_rate;
        const double Ts = 1.0 / fs_;

        // Companion conductances for every capacitor.
        c1_.prepare(p_.C1, Ts);
        c2_.prepare(p_.C2, Ts);  c3_.prepare(p_.C3, Ts);
        c4_.prepare(p_.C4, Ts);  c5_.prepare(p_.C5, Ts);
        c6_.prepare(p_.C6, Ts);
        c7_.prepare(p_.C7, Ts);  c8_.prepare(p_.C8, Ts);
        c9_.prepare(p_.C9, Ts);  c13_.prepare(p_.C13, Ts);

        // Fixed (pot-independent) resistor conductances.
        G2_  = 1.0 / p_.R2;  G3_  = 1.0 / p_.R3;
        G4_  = 1.0 / p_.R4;  G5_  = 1.0 / p_.R5;  G6_  = 1.0 / p_.R6;
        G7_  = 1.0 / p_.R7;  G8_  = 1.0 / p_.R8;  G9_  = 1.0 / p_.R9;
        G10_ = 1.0 / p_.R10;
        G11_ = 1.0 / p_.R11; G12_ = 1.0 / p_.R12; G13_ = 1.0 / p_.R13;
        G19_ = 1.0 / p_.R19; G20_ = 1.0 / p_.R20;

        // U2 diode clamp: three identical series diodes per direction collapse
        // to one element with emission coefficient (n * N). Vt at SPICE's 27 C: the reference simulator's, measured
        // (bjt_dc_bias.h).
        constexpr double VT = V_T_300K;
        Is_   = p_.diodes.Is;
        nvt3_ = p_.diodes.count_fwd * p_.diodes.N * VT;

        // U2 LM741 macromodel (finite GBW + slew): a transconductance input
        // stage gm, whose output current is LIMITED to +/-Imax (= the slew
        // mechanism), charges the compensation cap Cc with leak Rp (= DC gain).
        //   GBW = gm/(2*pi*Cc) ;  slew = Imax/Cc ;  Aol = gm*Rp.
        // These MUST equal the SPICE deck's .param gm741/imax741/rp741/cc741.
        Cc_   = p_.u2.comp_cap;
        gm_   = 2.0 * M_PI * p_.u2.gbw_hz * Cc_;
        Imax_ = (p_.u2.slew_v_per_us * 1.0e6) * Cc_;
        Rp_   = p_.u2.open_gain / gm_;
        CcTs_ = Cc_ * fs_;                       // = Cc / Ts

        // U1A's network has no pot in it, so factor it once here. The pot-
        // dependent blocks are factored by the set_*() calls below.
        refactor_u1a();
        factored_ = false;       // the setters below must run in full
        prepared_ = true;
        set_sustain(p_sus_);     // -> refactor_u1b()
        set_tone(p_tone_);       // -> refactor_tone()
        set_volume(p_vol_);      // -> refactor_tone()
        set_bypass(bypass_);     // -> refactor_tone()
        factored_ = true;
        reset();
    }

    void reset() noexcept {
        // Seed every state at the DC operating point so there is no startup
        // transient. With no signal: U1A/U1B outputs rest at +4.5 V, U2 at
        // +4.945 V, and the post-clip (AC-coupled) nodes rest at 0 V.
        const double b1 = p_.u1b.v_bias;     // 4.5  V (U1 bias)
        const double su = p_.u2.v_bias;      // 4.945 V (U2 bias / virtual gnd)

        // U1A: input cap holds (nB - Vin) = 4.5 V at rest.
        c1_.set_dc(b1);

        // U1B nodes all rest at 4.5 V; nC5 is tied through R9 to U2's 4.945 V.
        x_ = { b1, b1, b1, b1, b1, b1, b1, su, 0.0 };
        c2_.set_dc(b1);          // C2:  up_p - gnd
        c3_.set_dc(0.0);         // C3:  nX  - vy   (both 4.5)
        c4_.set_dc(b1);          // C4:  nZ  - gnd
        c5_.set_dc(b1 - su);     // C5:  nsw - nC5
        // The Guru mod grounds R8 (the sustain pot's bottom), so at rest a DC
        // current runs from U1B's output (4.5 V: no current in R6/R7 with the
        // capacitors open) through the sustain pot's two legs and R8 to ground,
        // and nsw/nsb sit on that divider -- where C5 charges. Seeding them at
        // 4.5 V instead thumps on every entry (0.2-0.4 V, ~0.4 s).
        if (p_.guru_mod && Gst_ > 0.0 && Gsb_ > 0.0 && G8_ > 0.0) {
            const double r_st = 1.0 / Gst_, r_sb = 1.0 / Gsb_, r_8 = 1.0 / G8_;
            const double io   = b1 / (r_st + r_sb + r_8);
            x_[NSW] = io * (r_sb + r_8);
            x_[NSB] = io * r_8;
            x_[IO]  = io;
            c5_.set_dc(x_[NSW] - su);
        }

        // U2 state: comp voltage and integrator history both zero at rest.
        su_ = su; vcomp_ = 0.0; rprev_ = 0.0; vcomp_pp_ = 0.0;
        c6_.set_dc(0.0);         // C6:  u2out - su (both 4.945)

        // Tone caps: post-clip nodes rest at 0 V, U2 out at su -> C8/C7/C13
        // hold -su; C9 sits at 0.
        c8_.set_dc(-su);         // C8:  HI   - u2out
        c7_.set_dc(-su);         // C7:  R11A - u2out
        c9_.set_dc(0.0);         // C9:  LO   - gnd
        c13_.set_dc(-su);        // C13: NC13 - u2out

        // The previous-sample inputs of each block, as a fresh instance holds
        // them. Left untouched, a non-finite one would survive both the guard
        // and activate().
        vin_ = 0.0; va_ = 0.0; tone_in_ = 0.0;

        // The factorisations are built AROUND the state (factor_kcl takes
        // differences at the current point): rebuilt at the rest, so a reset
        // engine is the engine a fresh prepare() gives whatever ran before, and a
        // non-finite state cannot leave a non-finite factorisation behind.
        if (prepared_) { refactor_u1a(); refactor_u1b(); refactor_tone(); }
    }

    // ----- controls (physical wiper fractions 0..1; taper applied above) ---

    // A setter that is handed its current value does nothing. The host calls
    // them on every block, and a refactor is not neutral: factor_kcl() builds
    // the matrix by differences around the CURRENT state, so re-factoring at
    // each block boundary would make the audio depend on the block size.
    void set_sustain(double frac) noexcept {
        if (factored_ && clamp01(frac) == p_sus_) return;
        p_sus_ = clamp01(frac);
        // Sustain pot (top = U1B output vy, bottom = nsb, wiper = nsw).
        // frac = 1 -> wiper at the top = maximum drive into the clipper. The
        // 1 uOhm floor keeps a fully-CW/CCW pot from being a singular short.
        Gst_ = 1.0 / (p_.R_sustain * (1.0 - p_sus_) + 1e-6);
        Gsb_ = 1.0 / (p_.R_sustain *        p_sus_  + 1e-6);
        if (prepared_) refactor_u1b();   // before prepare() the companions are zero: a singular system
    }
    void set_tone(double frac) noexcept {
        if (factored_ && clamp01(frac) == p_tone_) return;
        p_tone_ = clamp01(frac);
        // Tone pot (top = treble leg HI, bottom = bass leg LO, wiper = TW).
        Gtt_ = 1.0 / (p_.R_tone * (1.0 - p_tone_) + 1e-6);
        Gtb_ = 1.0 / (p_.R_tone *        p_tone_  + 1e-6);
        if (prepared_) refactor_tone();
    }
    void set_volume(double frac) noexcept {
        if (factored_ && clamp01(frac) == p_vol_) return;
        p_vol_ = clamp01(frac);
        // Volume pot (top = vol_top, wiper = out, bottom = gnd). frac=1 -> max.
        Gvt_ = 1.0 / (p_.R_volume * (1.0 - p_vol_) + 1e-6);
        Gvb_ = 1.0 / (p_.R_volume *        p_vol_  + 1e-6);
        if (prepared_) refactor_tone();
    }
    // Tone-Bypass switch: false = tone stack engaged (scoop wiper -> volume);
    // true = bypassed (flat post-clip tap -> volume, around the scoop).
    void set_bypass(bool on) noexcept {
        if (factored_ && on == bypass_) return;
        bypass_ = on;
        Gswe_ = on ? 0.0   : 1.0e6;   // engaged leg  (tone wiper -> volume top)
        Gswb_ = on ? 1.0e6 : 0.0;     // bypass  leg  (flat tap    -> volume top)
        if (prepared_) refactor_tone();
    }

    const char* name() const noexcept { return p_.name; }

    // ----- per-sample processing -------------------------------------------
    // `vin` is the instantaneous input-jack voltage (the model's native unit).
    // Returns the AC voltage at the requested tap (the offline driver applies
    // the same x0.1 scale the SPICE deck uses; in the plugin only `Out` runs).
    TOCAKE_SPAN_BODY float process_sample(float vin, OpTap tap = OpTap::Out) noexcept {
        const double va = solve_u1a(static_cast<double>(vin));     // U1A output (abs V)
        if (tap == OpTap::U1A) return guard(va - p_.u1a.v_bias);

        solve_u1b(va);                                             // fills x_[]
        // U1B and U2 are one solve: U2 finds su, then U1B settles at it.
        const double u2out = solve_u2(x_[NC5], u1b_dsu_[u1b_region_][NC5]);   // LM741 clipper
        settle_u1b(su_);
        if (tap == OpTap::U1B) return guard(x_[VY] - p_.u1b.v_bias);
        if (tap == OpTap::Clip) return guard(u2out - p_.u2.v_bias);

        return guard(solve_tone(u2out));                           // tone + volume (AC)
    }

    // The ISA dispatch's two entries (src/isa_dispatch.h): the same engine
    // compiled twice, once per instruction set, each kept out of line.
    // This engine takes the OTHER shape of the two (see isa_hot.h): its body is
    // inlined INTO the loop rather than called per sample, which pays here where
    // it costs the DK engines, and roughly doubles what the ISA buys.
    // It also needs every per-sample callee INSIDE the entry: an out-of-line
    // callee is compiled once, at the baseline ISA, whichever entry calls it.
    // That is why `solve_tone` is forced inline.
    void span_run(float* buf, uint32_t n) noexcept {
        for (uint32_t j = 0; j < n; ++j) buf[j] = process_sample(buf[j], OpTap::Out);
    }
    TOCAKE_SPAN_RUN_ISA void span_run_isa_v3(float* buf, uint32_t n) noexcept {
        for (uint32_t j = 0; j < n; ++j) buf[j] = process_sample(buf[j], OpTap::Out);
    }

    // Count of samples on which the guard fired. Reported by the harness; a
    // guard nobody reads is an alarm with no wire.
    long nonfinite_resets() const noexcept { return dbg_nonfinite_; }

private:
    // ----- non-finite guard -------------------------------------------------
    //
    // This engine ships (lv2.cpp includes it for the three op-amp slots) and
    // is the only one that picks between two pre-factored solves with an
    // `if` -- so it is the likeliest to hand a bad solve back.
    //
    // Two properties (src/finite_guard.h):
    //
    //  * COVERAGE. Every piece of dynamic state, not the output node. This
    //    engine keeps three separate kinds: the 9 U1B unknowns, the LM741
    //    macromodel's own state (su_, vcomp_, rprev_), and ten trapezoidal
    //    capacitor histories. A NaN that appears in a capacitor history is
    //    invisible at the output for one sample and poisons every sample
    //    after it.
    //
    //  * THE ANCHOR IS ABSOLUTE. Recovery is reset(), the solved DC rest
    //    point, never "the last good state": a fallback that restores what
    //    already failed re-fails immediately and leaves the engine dead for
    //    the rest of the session. reset() re-seeds all three kinds of state,
    //    which is why it is the whole repair rather than half of one.
    //
    // Every public return passes through here -- the belt protects the paths
    // that are exercised, not just the one that ships, because the offline
    // taps are where a bad state is first seen.
    bool state_finite() const noexcept {
        if (!tocake::all_finite(x_)) return false;
        if (!tocake::all_finite({su_, vcomp_, rprev_})) return false;
        if (!tocake::all_finite({vin_, va_, tone_in_})) return false;
        for (const CapTrap* c : {&c1_, &c2_, &c3_, &c4_, &c5_,
                                 &c6_, &c7_, &c8_, &c9_, &c13_})
            if (!tocake::all_finite({c->vp, c->ip})) return false;
        return true;
    }

    float guard(double v) noexcept {
        if (!state_finite() || !std::isfinite(v)) {
            reset();     // rebuilds the factorisations at the rest too (a non-finite one would mean silence)
            ++dbg_nonfinite_;
            return 0.0f;      // the AC rest value of every tap
        }
        return static_cast<float>(v);
    }

    long dbg_nonfinite_ = 0;

    static double clamp01(double v) noexcept {
        return v < 0.0 ? 0.0 : (v > 1.0 ? 1.0 : v);
    }

    // ===== generic linear-network helpers ==================================
    //
    // Each linear block's Kirchhoff residual f(v) is AFFINE: f(v) = A*v - rhs,
    // where A depends only on conductances (so it changes only when a knob
    // moves) and `rhs` carries the per-sample constants (capacitor histories
    // and input voltages). Two helpers exploit that:
    //
    //   build_matrix:  A[:,c] = f(e_c) - f(0)   (the per-sample constants
    //                  cancel in the difference, so it can run any time).
    //   build_rhs:     b = -f(0) = A*x - f(x) at x=0, i.e. the solution of
    //                  A*x = b is f(x) = 0.   (uses the live constants.)
    //
    // We build_matrix + lu_factor on a knob change, and build_rhs + lu_solve
    // per sample.

    // A linear block's solver: its LU and the inverse built from it at the same
    // knob move. The matrices only change on
    // a knob move, so the per-sample solve becomes a matrix-vector product:
    // no division and no chain of dependent multiply-subtracts, which is what
    // the substitutions cost (about 60 % of this engine, profiled).
    template <std::size_t N>
    struct LinSolve {
        LUFactors<N> lu;
        std::array<std::array<double, N>, N> inv {};
    };

    template <std::size_t N, class KCL>
    static LinSolve<N> factor_kcl(KCL&& kcl) noexcept {
        std::array<double, N> z {}, f0 {};
        kcl(z, f0);
        std::array<std::array<double, N>, N> A {};
        std::array<double, N> col {}, fc {};
        for (std::size_t c = 0; c < N; ++c) {
            col = z; col[c] = 1.0;
            kcl(col, fc);
            for (std::size_t r = 0; r < N; ++r) A[r][c] = fc[r] - f0[r];
        }
        LinSolve<N> ls;
        ls.lu = lu_factor<N>(A);
        for (std::size_t c = 0; c < N; ++c) {
            std::array<double, N> e {};
            e[c] = 1.0;
            lu_solve<N>(ls.lu, e);
            for (std::size_t r = 0; r < N; ++r) ls.inv[r][c] = e[r];
        }
        return ls;
    }
    template <std::size_t N, class KCL>
    static std::array<double, N> solve_kcl(const LinSolve<N>& ls, KCL&& kcl) noexcept {
        std::array<double, N> z {}, f0 {};
        kcl(z, f0);
        std::array<double, N> b {};
        for (std::size_t r = 0; r < N; ++r) b[r] = -f0[r];
        std::array<double, N> x {};
        for (std::size_t c = 0; c < N; ++c) {
            const double bc = b[c];
            for (std::size_t r = 0; r < N; ++r) x[r] += ls.inv[r][c] * bc;
        }
        return x;
    }

    // ===== U1A -- inverting booster ========================================
    // Nodes: nB (after input cap C1), sa (inverting input), va (output), plus
    // the op-amp output current Io. The op-amp drives va so its inputs match:
    // ideal -> sa = +4.5 V; if that needs va outside [0, 9] V the op-amp is
    // railed -> pin va to the rail and re-solve (region detection).
    enum { A_NB = 0, A_SA = 1, A_VA = 2, A_IO = 3 };
    void u1a_kcl(const std::array<double, 4>& v, std::array<double, 4>& f,
                 bool ideal, double rail) const noexcept {
        const double nB = v[A_NB], sa = v[A_SA], va = v[A_VA], Io = v[A_IO];
        const double iC1 = c1_.Gc * (nB - vin_) - c1_.ieq();   // C1: nB - Vin
        f[A_NB] = iC1 + (nB - sa) * G2_;
        f[A_SA] = (sa - nB) * G2_ + (sa - va) * G3_;           // op-amp input: no current
        f[A_VA] = (va - sa) * G3_ - Io;                        // op-amp sources Io into va
        f[A_IO] = ideal ? (sa - p_.u1a.v_bias) : (va - rail);  // V+ = V-   OR   va pinned
    }
    void refactor_u1a() noexcept {
        lu_u1a_ideal_ = factor_kcl<4>([this](auto& v, auto& f){ u1a_kcl(v, f, true,  0.0); });
        lu_u1a_sat_   = factor_kcl<4>([this](auto& v, auto& f){ u1a_kcl(v, f, false, 0.0); });
    }
    double solve_u1a(double vin) noexcept {
        vin_ = vin;
        auto x = solve_kcl<4>(lu_u1a_ideal_, [this](auto& v, auto& f){ u1a_kcl(v, f, true, 0.0); });
        if (x[A_VA] < 0.0 || x[A_VA] > p_.vcc) {               // op-amp railed
            const double rail = (x[A_VA] > p_.vcc) ? p_.vcc : 0.0;
            x = solve_kcl<4>(lu_u1a_sat_, [this, rail](auto& v, auto& f){ u1a_kcl(v, f, false, rail); });
        }
        c1_.commit(x[A_NB] - vin_);
        return x[A_VA];
    }

    // ===== U1B -- non-inverting Sustain stage ==============================
    // A coupled feedback network (R6/R7/C4/C3) whose gain leg also contains the
    // Sustain pot, which simultaneously taps the wiper toward the clipper. It
    // squares off against the rails at high Sustain, so it uses the same region
    // detection as U1A. Solved in absolute volts; x_[] holds the node voltages.
    // (nC5 feeds the clipper; the rest are internal.)
    enum { NX = 0, UP = 1, UM = 2, VY = 3, NZ = 4, NSB = 5, NSW = 6, NC5 = 7, IO = 8 };
    void u1b_kcl(const std::array<double, 9>& v, std::array<double, 9>& f,
                 bool ideal, double rail) const noexcept {
        const double nX = v[NX], up = v[UP], um = v[UM], vy = v[VY], nZ = v[NZ],
                     nsb = v[NSB], nsw = v[NSW], nC5 = v[NC5], Io = v[IO];
        const double su  = su_kcl_;                            // U2's inverting input (the bias unless coupled)
        const double iC2 = c2_.Gc * (up)        - c2_.ieq();   // C2: up_p - gnd
        const double iC3 = c3_.Gc * (nX - vy)   - c3_.ieq();   // C3: nX   - vy
        const double iC4 = c4_.Gc * (nZ)        - c4_.ieq();   // C4: nZ   - gnd
        const double iC5 = c5_.Gc * (nsw - nC5) - c5_.ieq();   // C5: nsw  - nC5
        f[NX]  = (nX - va_) * G4_ + (nX - up) * G5_ + iC3;
        f[UP]  = (up - nX) * G5_ + iC2;                        // + input: no current
        f[UM]  = (um - vy) * G6_ + (um - nZ) * G7_;            // - input: no current
        f[VY]  = (vy - um) * G6_ - iC3 + (vy - nsw) * Gst_ - Io;
        // Feedback mod: R8 (sustain-pot bottom) grounded directly instead of
        // tied to feedback node Z. Guarded so the stock path stays bit-identical.
        const bool guru = p_.guru_mod;
        f[NZ]  = (nZ - um) * G7_ + iC4 + (guru ? 0.0 : (nZ - nsb) * G8_);
        f[NSB] = (guru ? nsb * G8_ : (nsb - nZ) * G8_) + (nsb - nsw) * Gsb_;
        f[NSW] = (nsw - vy) * Gst_ + (nsw - nsb) * Gsb_ + iC5;
        f[NC5] = -iC5 + (nC5 - su) * G9_;                      // R9 to U2's virtual ground
        f[IO]  = ideal ? (up - um) : (vy - rail);
    }
    void refactor_u1b() noexcept {
        lu_u1b_ideal_ = factor_kcl<9>([this](auto& v, auto& f){ u1b_kcl(v, f, true,  0.0); });
        lu_u1b_sat_   = factor_kcl<9>([this](auto& v, auto& f){ u1b_kcl(v, f, false, 0.0); });
        // U1B is linear in su: its response to a unit move of su, per region.
        // A difference of two solves is exact for a linear system whatever the
        // state (the state enters both alike).
        for (int sat = 0; sat < 2; ++sat) {
            const LinSolve<9>& ls = sat ? lu_u1b_sat_ : lu_u1b_ideal_;
            su_kcl_ = p_.u2.v_bias;
            const auto x0 = solve_kcl<9>(ls, [this, sat](auto& v, auto& f){ u1b_kcl(v, f, !sat, 0.0); });
            su_kcl_ = p_.u2.v_bias + 1.0;
            const auto x1 = solve_kcl<9>(ls, [this, sat](auto& v, auto& f){ u1b_kcl(v, f, !sat, 0.0); });
            for (int k = 0; k < 9; ++k) u1b_dsu_[sat][k] = x1[k] - x0[k];
        }
        su_kcl_ = p_.u2.v_bias;
    }
    void solve_u1b(double va) noexcept {
        va_ = va;
        x_ = solve_kcl<9>(lu_u1b_ideal_, [this](auto& v, auto& f){ u1b_kcl(v, f, true, 0.0); });
        if (x_[VY] < 0.0 || x_[VY] > p_.vcc) {                 // op-amp railed
            const double rail = (x_[VY] > p_.vcc) ? p_.vcc : 0.0;
            x_ = solve_kcl<9>(lu_u1b_sat_, [this, rail](auto& v, auto& f){ u1b_kcl(v, f, false, rail); });
            u1b_region_ = 1;
        } else {
            u1b_region_ = 0;
        }
    }
    // x_ at the solved su (coupled build), then the U1B capacitor histories.
    void settle_u1b(double su) noexcept {
        const double d = su - p_.u2.v_bias;
        for (int k = 0; k < 9; ++k) x_[k] += u1b_dsu_[u1b_region_][k] * d;
        commit_u1b();
    }
    void commit_u1b() noexcept {
        c2_.commit(x_[UP]);
        c3_.commit(x_[NX]  - x_[VY]);
        c4_.commit(x_[NZ]);
        c5_.commit(x_[NSW] - x_[NC5]);
    }

    // ===== U2 -- LM741 clipper (the heart) =================================
    // U2's inverting input `su` is its virtual ground (pinned near +4.945 V by
    // feedback; U2 never rails because the diodes clamp its output first). The
    // output is not instantaneous: it carries the LM741's compensation-
    // cap state (finite GBW + slew). Two coupled relations:
    //
    //   (1) KCL at su:  (vnc5-su)*G9 + q*(G10+Gc6) - ieq6 + 2*Is*sinh(q/3NVt) = 0
    //       with q = Vout - su,  Vout = bias + vcomp.   (the diode clamp + R10/C6
    //       feedback + R9 input balance to zero current into the op-amp input.)
    //
    //   (2) the op-amp integrator (trapezoidal):  the input-stage current
    //       Istage (LIMITED to +/-Imax = the slew) charges the comp cap:
    //          (Cc/Ts)*(vcomp - vcomp_prev) = 0.5*(R + R_prev),  R = Istage - vcomp/Rp.
    //
    // We substitute (2) into (1): in the linear region Istage = gm*(bias-su),
    // giving vcomp as a linear function of su, so the whole thing is a 1-D
    // Newton solve in su for the (smooth) diode sinh. The slew limit is handled
    // by REGION DETECTION (assume linear; if |Istage| > Imax the stage is
    // slewing -> clamp Istage to +/-Imax, fix vcomp, re-solve) -- NOT a stiff
    // limiter inside Newton, which would stall with zero gradient when slewing.
    // vnc5 is U1B's output at su = bias and dvnc5 its slope in su.
    double solve_u2(double vnc5, double dvnc5) noexcept {
        const double bias = p_.u2.v_bias;
        const double GpC  = G10_ + c6_.Gc;
        const double ie6  = c6_.ieq();
        const double vcp  = vcomp_, rp = rprev_;
        // Integrator solved for vcomp: vcomp = A + B*(bias - su) in the linear
        // region (Istage = gm*(bias-su)); vcomp = const when the stage slews.
        const double Dintg = CcTs_ + 0.5 / Rp_;
        const double Aintg = (CcTs_ * vcp + 0.5 * rp) / Dintg;
        const double Bintg = 0.5 * gm_ / Dintg;

        // Stop test of the 1-D U2 clipper solve: a residual stop, |KCL current|
        // < 0.1 nA (sub-microvolt node error, ~120 dB below the U2 -50 dB null
        // floor). Not a step test on |dsu|: that chases rounding noise for the
        // last few iterations. (An asinh change of variable adds a
        // non-convergent tail here; keep Newton in su.)
        const double inv_nvt3 = 1.0 / nvt3_;
        auto solve_su = [&](bool slewing, double vc_fixed) noexcept {
            double su = su_;
            for (int it = 0; it < TOCAKE_OP_SU_MAXIT; ++it) {
                const double vc      = slewing ? vc_fixed : (Aintg + Bintg * (bias - su));
                const double dvc_dsu = slewing ? 0.0      : -Bintg;
                const double q       = bias + vc - su;
                const double dq_dsu  = dvc_dsu - 1.0;
                // The pair's sinh and cosh from ONE exp (the DK does the same for its pairs): one libm
                // call plus a division per iteration instead of two.
                const double eq = std::exp(q * inv_nvt3), ieq = 1.0 / eq;
                const double f  = (vnc5 + dvnc5 * (su - bias) - su) * G9_ + q * GpC - ie6
                                  + Is_ * (eq - ieq);
                if (it > 0 && std::abs(f) < 1.0e-10) break;        // residual stop
                const double df = (dvnc5 - 1.0) * G9_ + (GpC + (Is_ * inv_nvt3) * (eq + ieq)) * dq_dsu;
                double dsu = f / df;
                if (dsu >  0.5) dsu =  0.5;          // step clamp (sinh is stiff)
                if (dsu < -0.5) dsu = -0.5;
                su -= dsu;
            }
            return su;
        };

        double su = solve_su(false, 0.0);                      // assume linear input stage
        const double istage = gm_ * (bias - su);
        if (std::abs(istage) > Imax_) {                        // ... it is actually slewing
            const double isat = (istage > 0.0) ? Imax_ : -Imax_;
            const double vc   = (CcTs_ * vcp + 0.5 * (isat + rp)) / Dintg;
            su      = solve_su(true, vc);
            vcomp_  = vc;
            rprev_  = isat - vcomp_ / Rp_;
        } else {
            vcomp_  = Aintg + Bintg * (bias - su);
            rprev_  = istage - vcomp_ / Rp_;
        }
        su_ = su;
        vcomp_pp_ = vcp;                                       // BDF2 arm: the state two samples back
        const double Vout = bias + vcomp_;
        c6_.commit(Vout - su_);                                // C6: u2out - su
        return Vout;
    }

    // ===== tone stack + volume + Tone-Bypass switch ========================
    // Passive, driven by U2's stiff output. DC-blocked by C8/C7/C13 so the
    // whole network is AC-only (output centred at 0 V). The switch is two
    // conductances: one ~short, one open, depending on set_bypass(). Linear ->
    // one cached LU solve per sample.
    enum { HI = 0, R11A = 1, LO = 2, TW = 3, NC13 = 4, NBYP = 5, VT = 6, OUT = 7 };
    void tone_kcl(const std::array<double, 8>& v, std::array<double, 8>& f) const noexcept {
        const double hi = v[HI], r11a = v[R11A], lo = v[LO], tw = v[TW],
                     nc13 = v[NC13], nbyp = v[NBYP], vt = v[VT], out = v[OUT];
        const double iC8  = c8_.Gc  * (hi   - tone_in_) - c8_.ieq();    // treble leg cap
        const double iC7  = c7_.Gc  * (r11a - tone_in_) - c7_.ieq();    // bass leg cap
        const double iC9  = c9_.Gc  * (lo)              - c9_.ieq();    // bass leg shunt cap
        const double iC13 = c13_.Gc * (nc13 - tone_in_) - c13_.ieq();   // bypass path cap
        f[HI]   = iC8 + hi * G12_ + (hi - tw) * Gtt_;
        f[R11A] = iC7 + (r11a - lo) * G11_;
        f[LO]   = (lo - r11a) * G11_ + iC9 + (lo - tw) * Gtb_;
        f[TW]   = (tw - hi) * Gtt_ + (tw - lo) * Gtb_ + (tw - vt) * Gswe_;
        f[NC13] = iC13 + (nc13 - nbyp) * G19_;
        f[NBYP] = (nbyp - nc13) * G19_ + nbyp * G20_ + (nbyp - vt) * Gswb_;
        f[VT]   = (vt - tw) * Gswe_ + (vt - nbyp) * Gswb_ + (vt - out) * Gvt_;
        f[OUT]  = (out - vt) * Gvt_ + out * Gvb_ + out * G13_;
    }
    void refactor_tone() noexcept {
        lu_tone_ = factor_kcl<8>([this](auto& v, auto& f){ tone_kcl(v, f); });
    }
    // Per-sample, so forced inline for the same reason as process_sample: an
    // out-of-line copy is outside the ISA clone that calls it.
    TOCAKE_SPAN_BODY double solve_tone(double u2out) noexcept {
        tone_in_ = u2out;
        auto x = solve_kcl<8>(lu_tone_, [this](auto& v, auto& f){ tone_kcl(v, f); });
        c8_.commit(x[HI]   - u2out);
        c7_.commit(x[R11A] - u2out);
        c9_.commit(x[LO]);
        c13_.commit(x[NC13] - u2out);
        return x[OUT];
    }

    // ===== state ===========================================================
    const OpCakeParams& p_;
    double fs_ = 48000.0;

    // knob positions (physical wiper fractions)
    double p_sus_ = 1.0, p_tone_ = 0.5, p_vol_ = 0.5;
    bool   bypass_ = false;
    bool   factored_ = false;   // the pot-dependent blocks match the knobs above
    bool   prepared_ = false;   // prepare() has run: the setters may factor

    // per-sample block inputs (set just before each solve)
    double vin_ = 0.0;        // U1A input  (jack volts)
    double va_ = 0.0;         // U1B input  (= U1A output, abs V)
    double tone_in_ = 0.0;    // tone input (= U2 output, abs V)

    // fixed resistor conductances
    double G2_ = 0, G3_ = 0, G4_ = 0, G5_ = 0, G6_ = 0, G7_ = 0, G8_ = 0, G9_ = 0;
    double G10_ = 0, G11_ = 0, G12_ = 0, G13_ = 0, G19_ = 0, G20_ = 0;
    // pot conductances (recomputed on a knob change)
    double Gst_ = 0, Gsb_ = 0;             // Sustain
    double Gtt_ = 0, Gtb_ = 0;             // Tone
    double Gvt_ = 0, Gvb_ = 0;             // Volume
    double Gswe_ = 1e6, Gswb_ = 0;         // Tone-Bypass switch legs

    // U2 diode + LM741 macromodel
    double Is_ = 0, nvt3_ = 1.0;
    double gm_ = 0, Imax_ = 0, Rp_ = 1e9, Cc_ = 0, CcTs_ = 1.0;
    double su_ = 0, vcomp_ = 0, rprev_ = 0;     // dynamic state
    double su_kcl_ = 0;                          // su as U1B's KCL sees it (the bias)
    std::array<std::array<double, 9>, 2> u1b_dsu_ {};   // U1B's response to su, ideal / railed
    int u1b_region_ = 0;
    double vcomp_pp_ = 0;                       // BDF2 arm only

    // capacitor companion models
    CapTrap c1_, c2_, c3_, c4_, c5_, c6_, c7_, c8_, c9_, c13_;

    // node solution of the U1B network (kept for the tap + cap commits)
    std::array<double, 9> x_ {};

    // cached LU factorisations (refactorised only when a knob moves)
    LinSolve<4> lu_u1a_ideal_, lu_u1a_sat_;
    LinSolve<9> lu_u1b_ideal_, lu_u1b_sat_;
    LinSolve<8> lu_tone_;
};

} // namespace tocake
