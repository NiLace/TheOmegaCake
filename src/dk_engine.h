// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// A DK engine for any of the transistor circuits: a Newton on the nonlinear
// port voltages (vbe and vbc of every transistor, one per antiparallel diode
// pair), with the linear network and the capacitor states folded into matrices
// built from the variant's component values.
//
// The circuit is a description, not code: a `Circuit` type lists its nodes,
// its resistors (as a function of the knobs), its capacitors, its transistors
// and diode pairs, and the nodes its taps read. Everything else -- the matrix
// build, the DC rest point, the per-sample Newton, the junction limiting, the
// guards -- is here, once. src/dk_standard.h is the reference circuit.
//
// Formulation. With x the node voltages, s the capacitor history terms and
// i(v) the device currents at the port voltages v = N x:
//
//     G x = b_c + b_u u + S s - N^T i(v)
//
// so with H = G^-1 every quantity the sample needs is linear in (1, u, s, i):
//
//     v   = p - K i,    p = N H (b_c + b_u u + S s),    K = N H N^T
//     vc  = Cm x        (capacitor voltages)            s' = 2 g_c vc - s
//     tap = e_tap x
//
// Only G depends on the knobs, so a knob move rebuilds the matrices and never
// touches the state: the capacitor histories are per-capacitor quantities,
// independent of G.
//
// Per sample: a quadratic extrapolation of v seeds the Newton on
// F(v) = v - p + K i(v) = 0, whose Jacobian I + K D has D block diagonal (2x2
// per transistor, 1x1 per pair). The exponential ports are limited like SPICE's
// pnjlim -- every vbe, every pair, and vbc when the law carries the reverse
// junction; a vbc of the forward-only law is linear and is not limited at all.
// The final currents are updated linearly from the last step (i + D dv), and a
// pair's sinh/cosh come from one exp.
//
// Circuit interface (all static):
//   NN, NB, ND, S, NT          node, transistor, pair, capacitor and tap counts
//   FULL_EM                    true: SPICE Gummel-Poon with the reverse term
//                              (Br = 1); false: forward-only (no BC junction)
//   kBjts[NB]   {c, b, e}      pins; a pin may be kGnd (an emitter to ground)
//   kPairs[ND]  {a, k}         the pair's port is v(a) - v(k)
//   NS, kSingles[NS] {a, k}    OPTIONAL single diodes (anode a, cathode k),
//                              i = Is (e^(v/NVt) - 1): a rectifier or a clamp,
//                              which an antiparallel pair cannot represent.
//                              A circuit that declares neither has NS = 0 and
//                              every loop below is empty.
//   single(variant, k) -> diode_params   (only if NS > 0)
//   kCaps[S]    {a, b}
//   kTaps[NT]   {node, ac}     Out, C1Out, C2Out[, IbOut]; ac: the tap reads its
//                              node minus its rest value
//   kSeedHigh[] nodes seeded at 0.6 V for the DC Newton (the bases)
//   resistors(variant, knobs, R)    calls R(a, b, ohms) once per resistor
//   cap_values(variant)             -> std::array<double, S>
//   bjt(variant, k)  -> bjt_params  pair(variant, k) -> diode_params
#pragma once

#include "isa_hot.h"

#include <cstdint>

#include "cake_variant.h"
#include "bjt_dc_bias.h"      // V_T_300K
#include "finite_guard.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <type_traits>

// Step criterion of the port Newton, in volts: the default of every instance,
// and what ships. A tighter criterion does not move the audio; arbiters that
// need the converged solution call set_tolerance().
#ifndef TOCAKE_DK_TOL
#define TOCAKE_DK_TOL 1.0e-3
#endif
// Seed of the Newton: 2 quadratic extrapolation of the last three solutions,
// 1 linear, 0 the previous solution.
#ifndef TOCAKE_DK_PREDICTOR
#define TOCAKE_DK_PREDICTOR 2
#endif

// The exponential of the device model.
#ifndef TOCAKE_EXP
#define TOCAKE_EXP(x) std::exp(x)
#endif

#ifndef TOCAKE_DK_NEWTON_CAP
#define TOCAKE_DK_NEWTON_CAP 32
#endif
#define TOCAKE_DK_INLINE [[gnu::always_inline]] inline
#ifndef TOCAKE_DK_REFINE
#define TOCAKE_DK_REFINE 1
#endif

namespace tocake {

// Known voltages a pin can name besides the unknown nodes.
inline constexpr int kDkGnd = -1, kDkVcc = -2, kDkIn = -3;

struct DkBjtPins  { int c, b, e; };
struct DkPairPins { int a, k; };
using DkSinglePins = DkPairPins;   // one diode (anode, cathode): same pins, not a pair
struct DkCapPins  { int a, b; };
struct DkTap      { int node; bool ac; };

// Knob positions as PHYSICAL wiper fractions (taper applied upstream), plus the
// tone-bypass switch of the circuits that have one.
struct DkKnobs { double sus = 1.0, tone = 0.5, vol = 0.5; bool bypass = false; };

// A circuit may ask for a tighter Newton step criterion than the global one by
// declaring `static constexpr double kNewtonTol`; the tighter of the two wins, so
// a build-wide -DTOCAKE_DK_TOL still lowers every family.
template <class C, class = void>
struct dk_circuit_tol { static constexpr double value = TOCAKE_DK_TOL; };
template <class C>
struct dk_circuit_tol<C, std::void_t<decltype(C::kNewtonTol)>> {
    static constexpr double value = C::kNewtonTol < TOCAKE_DK_TOL ? C::kNewtonTol : TOCAKE_DK_TOL;
};

// A circuit that carries its own supply network (a series resistor from the
// battery into a rail NODE) names the battery voltage with
// `static constexpr double kSupply`; kDkVcc is then the battery, not the rail,
// and the variant's `vcc` (the rail's rest value, which the static-rail
// arbiters use) is not read.
template <class C, class = void>
struct dk_circuit_supply { static constexpr bool has = false; };
template <class C>
struct dk_circuit_supply<C, std::void_t<decltype(C::kSupply)>> { static constexpr bool has = true; };

// While a knob GLIDES (lv2.cpp, knob_tick) the matrices are rebuilt every few
// samples; set_coarse(true) drops the iterative-refinement pass from those
// rebuilds, and set_coarse(false) -- the tick that lands on the port value --
// forces one full rebuild, so a settled knob runs exactly the matrices it ran
// before. While gliding, a pot moves only its two legs: instead of refactoring
// G, the solution of the LAST EXACT rebuild is corrected by the Sherman-
// Morrison-Woodbury identity for the few branches whose conductance changed.
template <class E> auto dk_set_coarse(E& e, bool c, int) -> decltype(e.set_coarse(c), void()) { e.set_coarse(c); }
template <class E> void dk_set_coarse(E&, bool, long) {}

// A circuit whose accuracy hangs on a FEW junctions may name them: the ports in
// `kTightPorts` must step under `kTightTol`, every other port under the global
// criterion. Without the declaration the whole circuit has one criterion.
template <class C, class = void>
struct dk_circuit_tight { static constexpr bool has = false; };
template <class C>
struct dk_circuit_tight<C, std::void_t<decltype(C::kTightPorts), decltype(C::kTightTol)>> {
    static constexpr bool has = true;
};

template <class Circuit> class NodalEngine;   // src/nodal_engine.h: the same description, a nodal Newton

template <class Circuit>
class DkEngine {
    // The nodal twin reuses this engine's device law, junction limiting and DC
    // rest point rather than carrying copies of them (one device equation, once).
    friend class NodalEngine<Circuit>;
public:
    // Taps read kTaps[0] (Out), [1] (C1Out), [2] (C2Out) and, where the circuit
    // has a fourth, [3] (IbOut, the input booster's output).
    enum class Tap { C1Out, C2Out, Out, IbOut };

    static constexpr int NN = Circuit::NN;
    static constexpr int NB = Circuit::NB;
    static constexpr int ND = Circuit::ND;
    // Optional third port block: single diodes. `NS` is read off the circuit
    // when it declares one (the same detection the predictor uses).
    template <class C> static constexpr int ns_of(...) { return 0; }
    template <class C> static constexpr int ns_of(decltype(C::NS)*) { return C::NS; }
    static constexpr int NS = ns_of<Circuit>(nullptr);
    static constexpr int kNDev = 4 * NB + ND + NS;   // devices()'s derivatives: 4 per BJT, 1 per diode, 1 per shunt
    static_assert(kNDev > 0, "a circuit with no device has nothing to solve");
    static constexpr int P  = 2 * NB + ND + NS;   // ports: (vbe, vbc) per transistor, the pairs, the single diodes
    // The per-sample products read TRANSPOSED copies of Ps, K, Xs and Xi, so the inner loop walks contiguous
    // memory, on the circuits large enough for that to pay.
    static constexpr bool kTransposed = P >= 14;
    static constexpr int PS = 2 * NB + ND;        // first single-diode port
    static constexpr int S  = Circuit::S;
    static constexpr int NT = Circuit::NT;
    static constexpr bool FULL_EM = Circuit::FULL_EM;

    // Does this circuit declare a generated sparse elimination for its port
    // Jacobian? Detected, so adding one to a circuit is a one-line change and
    // a circuit without one cannot silently lose it.
    template <class C, class = void> struct has_gen_lu : std::false_type {};
    template <class C> struct has_gen_lu<C, std::void_t<typename C::DkLu>> : std::true_type {};
    static constexpr bool kHasGenLu = has_gen_lu<Circuit>::value;
    // The Newton's seed: TOCAKE_DK_PREDICTOR unless the circuit names its own
    // (a circuit's best order is a measured property of that circuit).
    template <class C> static constexpr int predictor_of(...) { return TOCAKE_DK_PREDICTOR; }
    template <class C> static constexpr int predictor_of(decltype(C::kPredictor)*) { return C::kPredictor; }
    static constexpr int kPredictor = predictor_of<Circuit>(nullptr);

    explicit DkEngine(const CakeVariant& v,
                      double p_tone = 0.5,
                      double p_vol  = 0.5,
                      double p_sus  = 1.0) noexcept
        : v_(v) { k_.tone = p_tone; k_.vol = p_vol; k_.sus = p_sus; }

    void prepare(double sample_rate) {
        fs_ = sample_rate;
        dirty_ = true;
        reset();
    }

    // Rebuild at the current knobs and start from the whole circuit's DC rest
    // point at those knobs. It runs a DC Newton: bounded (kDcMaxIters dense
    // solves) but not cheap. The plugin
    // calls it from run() only on a deliberate model or HQ change; the
    // non-finite guard uses recover_nonfinite() (the rest re-solved only if a
    // knob has moved).
    void reset() noexcept {
        if (coarse_) { coarse_ = false; dirty_ = true; }   // a reset mid-glide starts from the EXACT matrices
        if (dirty_) rebuild();
        dc_seed_settled_ = solve_dc();
        rest_stale_ = !dc_seed_settled_;
        restore_rest();
        dbg_iters_ = dbg_calls_ = dbg_nonconv_ = dbg_nonfinite_ = dbg_rebuilds_ = 0;
        dbg_max_iters_ = 0;
        dbg_max_step_ = 0.0;
    }
    bool dc_seed_settled() const noexcept { return dc_seed_settled_; }

    // Physical wiper fractions (taper applied upstream). A move only flags the
    // matrices; they are rebuilt once, before the next sample. A setter handed
    // its current value does nothing.
    void set_sustain(double frac) noexcept { set_pot(k_.sus,  frac); }
    void set_tone(double frac)    noexcept { set_pot(k_.tone, frac); }
    void set_volume(double frac)  noexcept { set_pot(k_.vol,  frac); }
    void set_bypass(bool on)      noexcept { if (on != k_.bypass) { k_.bypass = on; dirty_ = true; rest_stale_ = true; } }
    void set_coarse(bool on) noexcept {
        if (coarse_ && !on) dirty_ = true;    // landing: the exact matrices again
        coarse_ = on;
    }

    // Newton step criterion in volts (default TOCAKE_DK_TOL). For arbiters.
    void set_tolerance(double volts) noexcept { tol_ = volts; }
    double tolerance() const noexcept { return tol_; }

    float process_sample(float vin_volts) noexcept {
        return run_sample(static_cast<double>(vin_volts), Tap::Out);
    }

    // The ISA dispatch's two entries (src/isa_dispatch.h). Both are the same
    // engine; they differ only in the instruction set they are compiled for, and
    // `noinline` keeps each one out of line on purpose: inlining this into the
    // caller's loop is slower (see isa_hot.h).
    TOCAKE_SPAN_BASE float span_sample(float vin_volts) noexcept {
        return run_sample(static_cast<double>(vin_volts), Tap::Out);
    }
    TOCAKE_SPAN_ISA float span_sample_isa_v3(float vin_volts) noexcept {
        return run_sample(static_cast<double>(vin_volts), Tap::Out);
    }

    // The span the dispatch drives. This engine keeps the body OUT of the loop
    // (see isa_hot.h), so these two are plain
    // loops over the pair above and carry no instructions of their own.
    void span_run(float* buf, uint32_t n) noexcept {
        for (uint32_t j = 0; j < n; ++j) buf[j] = span_sample(buf[j]);
    }
    void span_run_isa_v3(float* buf, uint32_t n) noexcept {
        for (uint32_t j = 0; j < n; ++j) buf[j] = span_sample_isa_v3(buf[j]);
    }
    float process_sample(float vin_volts, Tap tap) noexcept {
        return run_sample(static_cast<double>(vin_volts), tap);
    }

    long   dbg_iters()     const noexcept { return dbg_iters_; }
    long   dbg_calls()     const noexcept { return dbg_calls_; }
    long   dbg_nonconv()   const noexcept { return dbg_nonconv_; }
    int    dbg_max_iters() const noexcept { return dbg_max_iters_; }
    // Worst final Newton step since reset(), in volts. This engine stops on the
    // step, so that is its convergence figure; dbg_max_resid() returns the same
    // number for the tools written against the node engines' residual.
    double dbg_max_step()  const noexcept { return dbg_max_step_; }
    double dbg_max_resid() const noexcept { return dbg_max_step_; }
    long   nonfinite_resets() const noexcept { return dbg_nonfinite_; }
    long   rebuilds()      const noexcept { return dbg_rebuilds_; }
    // True once any network of this engine held more conductances than
    // kMaxBranches (see stamp()). Must stay false for every shipped circuit.
    bool   branch_overflow() const noexcept { return branch_overflow_; }
    // True once a resistor field was <= 0 or non-finite (see resistors()). Must stay false for every shipped circuit.
    bool   bad_part() const noexcept { return bad_part_; }

    // The DC rest point adopted by the last reset(), node by node.
    const std::array<double, NN>& rest_nodes() const noexcept { return x_rest_; }
    // Port voltages AFTER the last converged sample (never a trial point).
    const std::array<double, P>& solved_ports() const noexcept { return vp_; }

    // Matrices, exposed so a gate can compare them against an independent build.
    const std::array<double, P * P>& K() const noexcept { return K_; }
    const std::array<double, P * S>& Ps() const noexcept { return Ps_; }
    const std::array<double, P>& Pc() const noexcept { return Pc_; }
    const std::array<double, P>& Pu() const noexcept { return Pu_; }

private:
    void set_pot(double& p, double frac) noexcept {
        frac = std::clamp(frac, 0.0, 1.0);
        if (frac != p) { p = frac; dirty_ = true; rest_stale_ = true; }
    }

    // ---- the linear network ------------------------------------------------

    // Every conductance is also kept as a branch, so the DC residual can be
    // summed as g*(va - vb): a 1e-6 ohm pot leg is a 1e6 S entry, and G x
    // would cancel two ~4e6 A products into a residual of a few nanoamps.
    struct Branch { int a, b; double g; };
#ifndef TOCAKE_DK_MAX_BRANCHES
#define TOCAKE_DK_MAX_BRANCHES 96   // conductances one network can list as branches
#endif
    static constexpr int kMaxBranches = TOCAKE_DK_MAX_BRANCHES;
    struct Network {
        std::array<double, NN * NN> G{};
        std::array<double, NN> bc{}, bu{};
        std::array<Branch, kMaxBranches> br{};
        int nbr = 0;
        bool overflow = false;   // more conductances than kMaxBranches: see stamp()
        bool bad_part = false;   // a resistor field <= 0 or non-finite: not stamped, see resistors()
    };

    // A conductance goes into G always; into the branch list (used only to sum
    // the refinement residual) while there is room. A circuit with more than
    // kMaxBranches conductances is still solved, less refined, and says so:
    // branch_overflow(), which must stay false for every shipped circuit.
    static void stamp(Network& n, int a, int b, double g, double vcc) noexcept {
        if (n.nbr < kMaxBranches) n.br[n.nbr++] = { a, b, g };
        else n.overflow = true;
        for (int side = 0; side < 2; ++side) {
            const int i = side ? b : a, j = side ? a : b;
            if (i < 0) continue;
            n.G[i * NN + i] += g;
            if (j >= 0) n.G[i * NN + j] -= g;
            else if (j == kDkVcc) n.bc[i] += g * vcc;
            else if (j == kDkIn) n.bu[i] += g;
        }
    }

    // Resistors only: the network at DC, and the part of G that the rate does not touch.
    Network resistors() const noexcept {
        Network n;
        const double vcc = supply();
        // A field <= 0 or non-finite (a generator's "absent part") would stamp an infinite conductance and silence
        // the slot: it is skipped and flagged (bad_part(), which must stay false).
        Circuit::resistors(v_, k_, [&](int a, int b, double r) {
            if (!(r > 0.0) || !std::isfinite(r)) { n.bad_part = true; return; }
            stamp(n, a, b, 1.0 / r, vcc);
        });
        return n;
    }

    // Port p is the voltage between two pins: (plus, minus).
    struct PortPins { int plus, minus; };
    static constexpr std::array<PortPins, P> make_ports() noexcept {
        std::array<PortPins, P> p{};
        for (int k = 0; k < NB; ++k) {
            p[2 * k]     = { Circuit::kBjts[k].b, Circuit::kBjts[k].e };   // vbe
            p[2 * k + 1] = { Circuit::kBjts[k].b, Circuit::kBjts[k].c };   // vbc
        }
        for (int k = 0; k < ND; ++k) p[2 * NB + k] = { Circuit::kPairs[k].a, Circuit::kPairs[k].k };
        if constexpr (NS > 0)
            for (int k = 0; k < NS; ++k) p[PS + k] = { Circuit::kSingles[k].a, Circuit::kSingles[k].k };
        return p;
    }
    static constexpr std::array<PortPins, P> kPorts = make_ports();

    static double pins_of(const double* x, int a, int b) noexcept {
        return (a >= 0 ? x[a] : 0.0) - (b >= 0 ? x[b] : 0.0);
    }
    static double port_of(const double* x, int p) noexcept {
        return pins_of(x, kPorts[p].plus, kPorts[p].minus);
    }

    // In-place LU with partial pivoting of an NN x NN matrix; perm records the rows.
    static bool lu_pivot(std::array<double, NN * NN>& A, std::array<int, NN>& perm) noexcept {
        for (int i = 0; i < NN; ++i) perm[i] = i;
        for (int k = 0; k < NN; ++k) {
            int piv = k; double best = std::fabs(A[k * NN + k]);
            for (int r = k + 1; r < NN; ++r)
                if (std::fabs(A[r * NN + k]) > best) { best = std::fabs(A[r * NN + k]); piv = r; }
            if (!(best > 0.0)) return false;
            if (piv != k) {
                for (int c = 0; c < NN; ++c) std::swap(A[k * NN + c], A[piv * NN + c]);
                std::swap(perm[k], perm[piv]);
            }
            const double inv = 1.0 / A[k * NN + k];
            for (int r = k + 1; r < NN; ++r) {
                const double f = A[r * NN + k] * inv;
                A[r * NN + k] = f;
                if (f != 0.0)
                    for (int c = k + 1; c < NN; ++c) A[r * NN + c] -= f * A[k * NN + c];
            }
        }
        return true;
    }
    static void lu_pivot_solve(const std::array<double, NN * NN>& A, const std::array<int, NN>& perm,
                               const double* b, double* x) noexcept {
        std::array<double, NN> y{};
        for (int r = 0; r < NN; ++r) {
            double acc = b[perm[r]];
            for (int c = 0; c < r; ++c) acc -= A[r * NN + c] * y[c];
            y[r] = acc;
        }
        for (int r = NN - 1; r >= 0; --r) {
            double acc = y[r];
            for (int c = r + 1; c < NN; ++c) acc -= A[r * NN + c] * x[c];
            x[r] = acc / A[r * NN + r];
        }
    }

    // The same substitutions for M right-hand sides at once. B and X are NN x M,
    // row-major, so the innermost loop runs along a row of M independent columns.
    // Column j sees the operations of lu_pivot_solve in the same order.
    template <int M>
    static void lu_pivot_solve_many(const std::array<double, NN * NN>& A, const std::array<int, NN>& perm,
                                    const double* B, double* X) noexcept {
        std::array<double, NN * M> Y;
        for (int r = 0; r < NN; ++r) {
            double* yr = &Y[r * M];
            const double* br = B + perm[r] * M;
            for (int j = 0; j < M; ++j) yr[j] = br[j];
            for (int c = 0; c < r; ++c) {
                const double a = A[r * NN + c];
                const double* yc = &Y[c * M];
                for (int j = 0; j < M; ++j) yr[j] -= a * yc[j];
            }
        }
        for (int r = NN - 1; r >= 0; --r) {
            double* xr = X + r * M;
            const double* yr = &Y[r * M];
            for (int j = 0; j < M; ++j) xr[j] = yr[j];
            for (int c = r + 1; c < NN; ++c) {
                const double a = A[r * NN + c];
                const double* xc = X + c * M;
                for (int j = 0; j < M; ++j) xr[j] -= a * xc[j];
            }
            const double d = A[r * NN + r];
            for (int j = 0; j < M; ++j) xr[j] /= d;
        }
    }

    // Build every matrix the sample needs, at the current rate and knobs.
    void rebuild() noexcept {
        Network n = resistors();
        const auto C = Circuit::cap_values(v_);
        for (int k = 0; k < S; ++k) {
            gc_[k] = 2.0 * C[k] * fs_;
            stamp(n, Circuit::kCaps[k].a, Circuit::kCaps[k].b, gc_[k], supply());
        }
        branch_overflow_ |= n.overflow;
        bad_part_ |= n.bad_part;
        constexpr int M = 2 + S + P;
        std::array<double, NN * M> Xup{};
        const bool updated = coarse_ && base_ok_ && fs_ == base_fs_ && rank_update(n, Xup);
        std::array<int, NN> perm{};
        std::array<double, NN * NN> LU;      // factorised only when the rank update has not solved it
        if (!updated) {
            LU = n.G;
            if (!lu_pivot(LU, perm)) { dirty_ = false; return; }   // unreachable for a connected network
        }

        std::array<double, NN> x{};
        auto rows = [&](const double* xx, double* p_out, double* vc_out, double* o_out) {
            for (int p = 0; p < P; ++p) p_out[p] = port_of(xx, p);
            for (int k = 0; k < S; ++k) vc_out[k] = pins_of(xx, Circuit::kCaps[k].a, Circuit::kCaps[k].b);
            for (int t = 0; t < NT; ++t) o_out[t] = pins_of(xx, Circuit::kTaps[t].node, kDkGnd);
        };
        double pp[P], vc[S], oo[NT];

        // Columns: 0 = bc, 1 = bu, 2..2+S-1 = the capacitor sources, then the ports.
        std::array<double, NN * M> Bm{}, Xm{};
        if (updated) {
            Xm = Xup;
        } else {
        for (int k = 0; k < NN; ++k) { Bm[k * M + 0] = n.bc[k]; Bm[k * M + 1] = n.bu[k]; }
        for (int j = 0; j < S; ++j) {
            if (Circuit::kCaps[j].a >= 0) Bm[Circuit::kCaps[j].a * M + 2 + j] += 1.0;
            if (Circuit::kCaps[j].b >= 0) Bm[Circuit::kCaps[j].b * M + 2 + j] -= 1.0;
        }
        for (int j = 0; j < P; ++j) {
            if (kPorts[j].plus >= 0)  Bm[kPorts[j].plus * M + 2 + S + j] += 1.0;
            if (kPorts[j].minus >= 0) Bm[kPorts[j].minus * M + 2 + S + j] -= 1.0;
        }
        lu_pivot_solve_many<M>(LU, perm, Bm.data(), Xm.data());
        for (int pass = 0; pass < (coarse_ ? 0 : kRefinePasses); ++pass) {
            std::array<double, NN * M> Rm = Bm, Dm;
            for (int k = 0; k < n.nbr; ++k) {
                const int a = n.br[k].a, b = n.br[k].b;
                const double g = n.br[k].g;
                for (int j = 0; j < M; ++j) {
                    const double cur = g * ((a >= 0 ? Xm[a * M + j] : 0.0) - (b >= 0 ? Xm[b * M + j] : 0.0));
                    if (a >= 0) Rm[a * M + j] -= cur;
                    if (b >= 0) Rm[b * M + j] += cur;
                }
            }
            lu_pivot_solve_many<M>(LU, perm, Rm.data(), Dm.data());
            for (int i = 0; i < NN * M; ++i) Xm[i] += Dm[i];
        }
        if (!coarse_) {                  // an EXACT rebuild becomes the base of the updates
            lu0_ = LU; perm0_ = perm; xm0_ = Xm; base_fs_ = fs_;
            nbr0_ = n.nbr;
            for (int k = 0; k < n.nbr; ++k) br0_[k] = n.br[k];
            base_ok_ = true;
        }
        }
        auto column = [&](int j) { for (int k = 0; k < NN; ++k) x[k] = Xm[k * M + j]; return x.data(); };
        rows(column(0), pp, vc, oo);
        for (int p = 0; p < P; ++p) Pc_[p] = pp[p];
        for (int k = 0; k < S; ++k) Xc_[k] = vc[k];
        for (int t = 0; t < NT; ++t) Oc_[t] = oo[t];
        rows(column(1), pp, vc, oo);
        for (int p = 0; p < P; ++p) Pu_[p] = pp[p];
        for (int k = 0; k < S; ++k) Xu_[k] = vc[k];
        for (int t = 0; t < NT; ++t) Ou_[t] = oo[t];
        for (int j = 0; j < S; ++j) {
            rows(column(2 + j), pp, vc, oo);
            for (int p = 0; p < P; ++p) Ps_[p * S + j] = pp[p];
            for (int k = 0; k < S; ++k) Xs_[k * S + j] = vc[k];
            for (int t = 0; t < NT; ++t) Os_[t * S + j] = oo[t];
        }
        for (int j = 0; j < P; ++j) {
            rows(column(2 + S + j), pp, vc, oo);
            for (int p = 0; p < P; ++p) K_[p * P + j] = pp[p];
            for (int k = 0; k < S; ++k) Xi_[k * P + j] = vc[k];
            for (int t = 0; t < NT; ++t) Oi_[t * P + j] = oo[t];
        }
        // A capacitor pin on the input jack or the rail: its voltage carries
        // that known potential, which the projection of x cannot see.
        for (int k = 0; k < S; ++k) {
            const int a = Circuit::kCaps[k].a, b = Circuit::kCaps[k].b;
            if (a == kDkIn)  Xu_[k] += 1.0;
            if (b == kDkIn)  Xu_[k] -= 1.0;
            if (a == kDkVcc) Xc_[k] += supply();
            if (b == kDkVcc) Xc_[k] -= supply();
        }
        for (int k = 0; k < S; ++k) g2_[k] = 2.0 * gc_[k];
        if constexpr (kTransposed) {
        // After mask_cascade(), which edits K.
        for (int r = 0; r < P; ++r)
            for (int k = 0; k < S; ++k) PsT_[k * P + r] = Ps_[r * S + k];
        for (int r = 0; r < P; ++r)
            for (int c = 0; c < P; ++c) KT_[c * P + r] = K_[r * P + c];
        for (int r = 0; r < S; ++r)
            for (int k = 0; k < S; ++k) XsT_[k * S + r] = Xs_[r * S + k];
        for (int r = 0; r < S; ++r)
            for (int c = 0; c < P; ++c) XiT_[c * S + r] = Xi_[r * P + c];
        }
        dirty_ = false;
        ++dbg_rebuilds_;
    }

    // ---- devices -------------------------------------------------------------

    // Port currents i and their Jacobian D: 2x2 per transistor (stored as
    // d[4q..4q+3] = [be,be], [be,bc], [bc,be], [bc,bc]) and 1x1 per pair.
    //
    // The transistor law is SPICE's Gummel-Poon for these cards (Ikf, Ise and
    // Var at their defaults): with ef = exp(vbe/Vt), er = exp(vbc/Vt),
    //     IC = Is (ef - er)(1 - vbc/Vaf) - Is/Br (er - 1)
    //     IB = Is/Bf (ef - 1) + Is/Br (er - 1),      Br = 1
    // and the forward-only form drops every er term (IC = Is (ef - 1)(1 - vbc/Vaf)).
    // The port currents are i_be = IB + IC (it enters the base and leaves the
    // emitter) and i_bc = -IC.
    TOCAKE_DK_INLINE void devices(const double* vv, double* i, double* d) const noexcept {
        for (int k = 0; k < NB; ++k) {
            const double Is = bjt_is_[k];
            const int be = 2 * k, bc = be + 1;
            const double ef = TOCAKE_EXP(vv[be] * inv_vt_);
            const double early = 1.0 - vv[bc] * inv_vaf_[k];
            if constexpr (!FULL_EM) {
                const double tr = Is * (ef - 1.0);
                const double ic = tr * early, ib = tr * inv_bf_[k];
                const double g = Is * ef * inv_vt_;
                const double dic_be = g * early, dic_bc = -tr * inv_vaf_[k];
                i[be] = ib + ic; i[bc] = -ic;
                d[4 * k + 0] = g * inv_bf_[k] + dic_be;
                d[4 * k + 1] = dic_bc;
                d[4 * k + 2] = -dic_be;
                d[4 * k + 3] = -dic_bc;
            } else {
                // A reverse-biased BC junction (vbc of -3 to -5 V here) makes
                // er < e^-36 = 2.3e-16: below one ulp of every sum it enters
                // (Is*(er - 1) rounds to -Is, ef - er to ef), so it is not
                // computed. The branch is taken on almost every sample.
                const double xr = vv[bc] * inv_vt_;
                const double er = (xr < -36.0) ? 0.0 : TOCAKE_EXP(xr);
                const double tr = Is * (ef - er);
                const double rev = Is * (er - 1.0);               // Is/Br (er - 1), Br = 1
                const double ic = tr * early - rev;
                const double ib = Is * (ef - 1.0) * inv_bf_[k] + rev;
                const double gf = Is * ef * inv_vt_, gr = Is * er * inv_vt_;
                const double dic_be = gf * early;
                const double dic_bc = -gr * early - tr * inv_vaf_[k] - gr;
                const double dib_be = gf * inv_bf_[k], dib_bc = gr;
                i[be] = ib + ic; i[bc] = -ic;
                d[4 * k + 0] = dib_be + dic_be;
                d[4 * k + 1] = dib_bc + dic_bc;
                d[4 * k + 2] = -dic_be;
                d[4 * k + 3] = -dic_bc;
            }
        }
        for (int k = 0; k < ND; ++k) {
            const int p = 2 * NB + k;
            const double e = TOCAKE_EXP(vv[p] * inv_nvt_[k]), ie = 1.0 / e;
            i[p] = is_d_[k] * (e - ie);
            d[4 * NB + k] = is_d_[k] * (e + ie) * inv_nvt_[k];
        }
        // Single diodes: i = Is (e^(v/NVt) - 1), the ONE-way law. A rectifier
        // into an envelope node and a clamp to ground are not the same element
        // as an antiparallel pair, and a pair cannot stand in for either.
        if constexpr (NS > 0) {
            for (int k = 0; k < NS; ++k) {
                const int p = PS + k;
                const double e = TOCAKE_EXP(vv[p] * inv_nvt_s_[k]);
                i[p] = is_s_[k] * (e - 1.0);
                d[4 * NB + ND + k] = is_s_[k] * e * inv_nvt_s_[k];
            }
        }
    }

    // SPICE's junction limiting (pnjlim): above the critical voltage a step
    // that would move the junction by more than two thermal voltages is
    // replaced by the step that moves its CURRENT by the same factor the linear
    // model predicted. Without it the port Newton overshoots into the
    // exponential and oscillates on hard, broadband drive (white noise into the
    // high-gain slots never converged on a third of its samples). A pair is
    // symmetric, so it is limited on its magnitude.
    static constexpr int kJunctions = NB + ND + NS + (FULL_EM ? NB : 0);
    // Junction j: the vbe of transistor j (j < NB), then the pairs, then --
    // with the reverse term -- the vbc of every transistor. The kind is known
    // at compile time in every loop below (limit_all), so no port reads a flag.
    struct JunctionPort { int port; bool pair; };   // pair: limited on its magnitude
    static constexpr std::array<JunctionPort, kJunctions> make_junctions() noexcept {
        std::array<JunctionPort, kJunctions> j{};
        int n = 0;
        for (int k = 0; k < NB; ++k) j[n++] = { 2 * k, false };
        for (int k = 0; k < ND; ++k) j[n++] = { 2 * NB + k, true };
        for (int k = 0; k < NS; ++k) j[n++] = { PS + k, false };
        if constexpr (FULL_EM)
            for (int k = 0; k < NB; ++k) j[n++] = { 2 * k + 1, false };
        return j;
    }
    static constexpr std::array<JunctionPort, kJunctions> kJunctionPort = make_junctions();

    template <bool Pair>
    TOCAKE_DK_INLINE double limit_junction(double vnew, double vold, int j) const noexcept {
        const double vt = jvt_[j], vcrit = jvcrit_[j];
        const double sg = (Pair && vnew < 0.0) ? -1.0 : 1.0;   // pairs: fold onto the positive side
        const double vn = sg * vnew, vo = sg * vold;
        if (vn > vcrit && std::fabs(vn - vo) > 2.0 * vt) {
            double lim;
            if (vo > 0.0) {
                const double arg = 1.0 + (vn - vo) / vt;
                lim = (arg > 0.0) ? vo + vt * std::log(arg) : vcrit;
            } else {
                lim = vt * std::log(vn / vt);
            }
            return sg * lim;
        }
        return vnew;
    }

    // Limit every exponential port of vnew against vold, in junction order.
    TOCAKE_DK_INLINE void limit_all(double* vnew, const double* vold) const noexcept {
        for (int k = 0; k < NB; ++k)
            vnew[2 * k] = limit_junction<false>(vnew[2 * k], vold[2 * k], k);
        for (int k = 0; k < ND; ++k)
            vnew[2 * NB + k] = limit_junction<true>(vnew[2 * NB + k], vold[2 * NB + k], NB + k);
        for (int k = 0; k < NS; ++k)
            vnew[PS + k] = limit_junction<false>(vnew[PS + k], vold[PS + k], NB + ND + k);
        if constexpr (FULL_EM)
            for (int k = 0; k < NB; ++k)
                vnew[2 * k + 1] = limit_junction<false>(vnew[2 * k + 1], vold[2 * k + 1], NB + ND + NS + k);
    }

    void load_device_constants() noexcept {
        inv_vt_ = 1.0 / V_T_300K;
        for (int k = 0; k < NB; ++k) {
            const bjt_params q = Circuit::bjt(v_, k);
            bjt_is_[k] = q.Is; inv_bf_[k] = 1.0 / q.Bf; inv_vaf_[k] = 1.0 / q.Vaf;
        }
        for (int k = 0; k < ND; ++k) {
            const diode_params dp = Circuit::pair(v_, k);
            is_d_[k] = dp.Is; inv_nvt_[k] = 1.0 / (dp.N * V_T_300K);
        }
        if constexpr (NS > 0) {
            for (int k = 0; k < NS; ++k) {
                const diode_params dp = Circuit::single(v_, k);
                is_s_[k] = dp.Is; inv_nvt_s_[k] = 1.0 / (dp.N * V_T_300K);
            }
        }
        for (int j = 0; j < kJunctions; ++j) {
            const int p = kJunctionPort[j].port;
            double vt, is;
            if (kJunctionPort[j].pair)  { const diode_params dp = Circuit::pair(v_, p - 2 * NB); vt = dp.N * V_T_300K; is = dp.Is; }
            else if (p >= PS && p < PS + NS) {
                if constexpr (NS > 0) { const diode_params dp = Circuit::single(v_, p - PS); vt = dp.N * V_T_300K; is = dp.Is; }
                else                  { vt = V_T_300K; is = 1.0; }          // unreachable
            }
            else                        { vt = V_T_300K; is = Circuit::bjt(v_, p / 2).Is; }
            jvt_[j] = vt;
            jvcrit_[j] = vt * std::log(vt / (std::sqrt(2.0) * is));
        }
    }

    // ---- DC rest point --------------------------------------------------------

    // Nodal Newton on the resistive network with the capacitors open, the
    // input grounded, and the device laws the sample runs.
    bool solve_dc() noexcept {
        load_device_constants();
        const Network n = resistors();
        std::array<double, NN> x{};
        for (int b : Circuit::kSeedHigh) x[b] = 0.6;
        bool settled = false;
        for (int it = 0; it < kDcMaxIters && !settled; ++it) {
            double vv[P], i[P], d[kNDev];
            for (int p = 0; p < P; ++p) vv[p] = port_of(x.data(), p);
            devices(vv, i, d);
            std::array<double, NN * NN> J = n.G;
            std::array<double, NN> F{};
            const double vcc = supply();
            auto node_v = [&](int k) { return k >= 0 ? x[k] : (k == kDkVcc ? vcc : 0.0); };
            for (int k = 0; k < n.nbr; ++k) {             // input grounded: kDkIn reads 0
                const double cur = n.br[k].g * (node_v(n.br[k].a) - node_v(n.br[k].b));
                if (n.br[k].a >= 0) F[n.br[k].a] += cur;
                if (n.br[k].b >= 0) F[n.br[k].b] -= cur;
            }
            auto add_port = [&](int p, double cur) {
                if (kPorts[p].plus >= 0)  F[kPorts[p].plus] += cur;
                if (kPorts[p].minus >= 0) F[kPorts[p].minus] -= cur;
            };
            auto add_jac = [&](int p, int q, double g) {   // d i_p / d v_q, v_q = x[q+] - x[q-]
                const int rp[2] = { kPorts[p].plus, kPorts[p].minus };
                const int cq[2] = { kPorts[q].plus, kPorts[q].minus };
                for (int a = 0; a < 2; ++a)
                    for (int b = 0; b < 2; ++b)
                        if (rp[a] >= 0 && cq[b] >= 0)
                            J[rp[a] * NN + cq[b]] += (a ? -1.0 : 1.0) * (b ? -1.0 : 1.0) * g;
            };
            for (int p = 0; p < P; ++p) add_port(p, i[p]);
            for (int k = 0; k < NB; ++k) {
                const int be = 2 * k, bc = be + 1;
                add_jac(be, be, d[4 * k + 0]); add_jac(be, bc, d[4 * k + 1]);
                add_jac(bc, be, d[4 * k + 2]); add_jac(bc, bc, d[4 * k + 3]);
            }
            for (int k = 0; k < ND; ++k) add_jac(2 * NB + k, 2 * NB + k, d[4 * NB + k]);
            for (int k = 0; k < NS; ++k) add_jac(PS + k, PS + k, d[4 * NB + ND + k]);
            std::array<int, NN> perm{};
            if (!lu_pivot(J, perm)) return false;
            for (int r = 0; r < NN; ++r) F[r] = -F[r];
            std::array<double, NN> dx{};
            lu_pivot_solve(J, perm, F.data(), dx.data());
            double step = 0.0;
            for (int r = 0; r < NN; ++r) {
                const double s = std::clamp(dx[r], -0.15, 0.15);
                x[r] += s;
                step = std::max(step, std::fabs(s));
            }
            if (!all_finite(x)) return false;
            settled = step < kDcTol;
        }
        if (settled) x_rest_ = x;
        return settled;
    }

    // A pin's voltage at rest: a node's rest value, the rail, or 0 (ground and
    // the silent input).
    double rest_pin(int n) const noexcept {
        return n >= 0 ? x_rest_[n] : (n == kDkVcc ? supply() : 0.0);
    }
    // The non-finite guard's path: back to the DC rest AT THE CURRENT KNOBS. x_rest_ is solved by reset(); after a
    // knob move it describes the old operating point, so it is re-solved here first (bounded: kDcMaxIters). Cold and
    // out of line: it never runs on a healthy sample.
    [[gnu::noinline, gnu::cold]] void recover_nonfinite() noexcept {
        if (rest_stale_ && solve_dc()) rest_stale_ = false;
        restore_rest();
        ++dbg_nonfinite_;
    }
    void restore_rest() noexcept {
        for (int k = 0; k < S; ++k) {
            const int a = Circuit::kCaps[k].a, b = Circuit::kCaps[k].b;
            s_[k] = gc_[k] * (rest_pin(a) - rest_pin(b));
        }
        for (int p = 0; p < P; ++p) vp_[p] = vp1_[p] = vp2_[p] = port_of(x_rest_.data(), p);
        for (int t = 0; t < NT; ++t)
            tap_rest_[t] = Circuit::kTaps[t].ac ? pins_of(x_rest_.data(), Circuit::kTaps[t].node, kDkGnd) : 0.0;
    }

    // ---- the sample -------------------------------------------------------------

    TOCAKE_SPAN_BODY float run_sample(double u, Tap tap) noexcept {
        if (dirty_) rebuild();

        double p[P];
        for (int r = 0; r < P; ++r) p[r] = Pc_[r] + Pu_[r] * u;
        for (int k = 0; k < S; ++k) {
            const double sk = s_[k];
            if constexpr (kTransposed) {
                const double* col = &PsT_[k * P];
                for (int r = 0; r < P; ++r) p[r] += col[r] * sk;
            } else {
                for (int r = 0; r < P; ++r) p[r] += Ps_[r * S + k] * sk;
            }
        }
        for (int k = 0; k < P; ++k) {                 // quadratic extrapolation
            const double cur = vp_[k];
            vp_[k] = (kPredictor == 2) ? 3.0 * cur - 3.0 * vp1_[k] + vp2_[k]
                   : (kPredictor == 1) ? 2.0 * cur - vp1_[k] : cur;
            vp2_[k] = vp1_[k]; vp1_[k] = cur;
        }
        // An extrapolated junction voltage is limited like a Newton step: on a
        // fast edge the extrapolation lands deep in the exponential.
        limit_all(vp_.data(), vp1_.data());

        double i[P], d[kNDev], J[P * P], ip[P], R[P], dv[P];
        int it = 0;
        bool conv = false;
        double step = 0.0;
        for (; it < TOCAKE_DK_NEWTON_CAP; ++it) {
            devices(vp_.data(), i, d);
            for (int r = 0; r < P; ++r) R[r] = p[r] - vp_[r];
            for (int c = 0; c < P; ++c) {
                const double ic = i[c];
                if constexpr (kTransposed) {
                    const double* col = &KT_[c * P];
                    for (int r = 0; r < P; ++r) R[r] -= col[r] * ic;
                } else {
                    for (int r = 0; r < P; ++r) R[r] -= K_[r * P + c] * ic;
                }
            }
            // J = I + K D, column pairs per transistor, one column per pair.
            for (int r = 0; r < P; ++r) {
                const double* kr = &K_[r * P];
                for (int q = 0; q < NB; ++q) {
                    const int be = 2 * q, bc = be + 1;
                    J[r * P + be] = kr[be] * d[4 * q + 0] + kr[bc] * d[4 * q + 2];
                    J[r * P + bc] = kr[be] * d[4 * q + 1] + kr[bc] * d[4 * q + 3];
                }
                for (int q = 0; q < ND; ++q) J[r * P + 2 * NB + q] = kr[2 * NB + q] * d[4 * NB + q];
                for (int q = 0; q < NS; ++q) J[r * P + PS + q] = kr[PS + q] * d[4 * NB + ND + q];
                J[r * P + r] += 1.0;
            }

            // LU without pivoting: I + K D is not guaranteed diagonally dominant, so this rests on measurement: it
            // holds on every shipped circuit's equations, and a zero pivot would end in the non-finite guard. The
            // reciprocal of each pivot is stored for the back substitution.
            //
            // A circuit whose port Jacobian has STRUCTURAL zeros declares a
            // generated elimination that skips exactly the operations whose
            // operand is one of them. That is bit-identical, not an
            // approximation: the dense form computes `J[r][c] -= f * 0.0`,
            // which is `J[r][c] - 0.0`, which is J[r][c]. Only the Jumbo Goat
            // qualifies; every other circuit here is structurally dense and
            // takes the generic loop below.
            if constexpr (kHasGenLu) {
                // Out of line ON PURPOSE, so it runs at the baseline ISA even
                // from the v3 entry. Do not force it inline into the entry: the
                // entry grows by the whole elimination and the Jumbo Goat runs
                // slower.
                Circuit::DkLu::factor(J, ip);
            } else {
#pragma GCC unroll 16
            for (int k = 0; k < P; ++k) {
                const double inv = 1.0 / J[k * P + k];
                ip[k] = inv;
#pragma GCC unroll 16
                for (int r = k + 1; r < P; ++r) {
                    const double f = J[r * P + k] * inv;
                    J[r * P + k] = f;
                    for (int c = k + 1; c < P; ++c) J[r * P + c] -= f * J[k * P + c];
                }
            }
            }
            // Substitutions by columns: each finished unknown is subtracted
            // from every row still open. The row form is a chain of dependent
            // multiply-subtracts per unknown; here the updates of one column are
            // independent.
#pragma GCC unroll 16
            for (int c = 0; c < P - 1; ++c) {
                const double rc = R[c];
                for (int r = c + 1; r < P; ++r) R[r] -= J[r * P + c] * rc;
            }
#pragma GCC unroll 16
            for (int k = P - 1; k >= 0; --k) {
                const double rk = R[k] * ip[k];
                R[k] = rk;
                for (int r = 0; r < k; ++r) R[r] -= J[r * P + k] * rk;
            }
            step = 0.0;
            [[maybe_unused]] double step_tight = 0.0;
            double vnew[P];
            // No fixed step clamp: a vbc port of the forward-only law is linear
            // (it only scales the Early factor), and a 0.15 V clamp walks an
            // extrapolated vbc back from +4.4 V one clamp at a time until the
            // iteration cap. The exponential ports are limited by pnjlim below.
            for (int r = 0; r < P; ++r) vnew[r] = vp_[r] + R[r];
            limit_all(vnew, vp_.data());
            for (int r = 0; r < P; ++r) {
                const double s = vnew[r] - vp_[r];
                dv[r] = s; vp_[r] = vnew[r];
                if constexpr (dk_circuit_tight<Circuit>::has) {
                    if (kTightMask[r]) { step_tight = std::max(step_tight, std::fabs(s)); continue; }
                }
                step = std::max(step, std::fabs(s));
            }
            if constexpr (dk_circuit_tight<Circuit>::has) {
                // NaN-safe like the test below: a non-finite step fails both.
                const double t = Circuit::kTightTol < tol_ ? Circuit::kTightTol : tol_;
                if (!(step_tight < t)) step = step + tol_ + step_tight;   // >= tol_, and NaN stays NaN
            }
            // NaN-safe: a non-finite step never passes this test.
            if (step < tol_) { ++it; conv = true; break; }
        }
        ++dbg_calls_;
        dbg_iters_ += it;
        if (it > dbg_max_iters_) dbg_max_iters_ = it;
        if (!conv) ++dbg_nonconv_;
        if (!(step <= dbg_max_step_)) dbg_max_step_ = step;   // a NaN step is kept, not hidden

        // Currents at the solution, linearised from the last evaluation.
        for (int q = 0; q < NB; ++q) {
            const int be = 2 * q, bc = be + 1;
            const double dbe = dv[be], dbc = dv[bc];
            i[be] += d[4 * q + 0] * dbe + d[4 * q + 1] * dbc;
            i[bc] += d[4 * q + 2] * dbe + d[4 * q + 3] * dbc;
        }
        for (int q = 0; q < ND; ++q) i[2 * NB + q] += d[4 * NB + q] * dv[2 * NB + q];
        for (int q = 0; q < NS; ++q) i[PS + q] += d[4 * NB + ND + q] * dv[PS + q];

        const int t = (tap == Tap::Out) ? 0 : (tap == Tap::C1Out) ? 1 : (tap == Tap::C2Out) ? 2
                    : (NT > 3 ? 3 : 0);
        double out = Oc_[t] + Ou_[t] * u;
        for (int k = 0; k < S; ++k) out += Os_[t * S + k] * s_[k];
        for (int c = 0; c < P; ++c) out -= Oi_[t * P + c] * i[c];
        out -= tap_rest_[t];

        double vc[S];
        for (int r = 0; r < S; ++r) vc[r] = Xc_[r] + Xu_[r] * u;
        for (int k = 0; k < S; ++k) {
            const double sk = s_[k];
            if constexpr (kTransposed) {
                const double* col = &XsT_[k * S];
                for (int r = 0; r < S; ++r) vc[r] += col[r] * sk;
            } else {
                for (int r = 0; r < S; ++r) vc[r] += Xs_[r * S + k] * sk;
            }
        }
        for (int c = 0; c < P; ++c) {
            const double ic = i[c];
            if constexpr (kTransposed) {
                const double* col = &XiT_[c * S];
                for (int r = 0; r < S; ++r) vc[r] -= col[r] * ic;
            } else {
                for (int r = 0; r < S; ++r) vc[r] -= Xi_[r * P + c] * ic;
            }
        }
        for (int k = 0; k < S; ++k) s_[k] = g2_[k] * vc[k] - s_[k];

        // Every port and every state, and back to the DC rest: restoring the
        // last state would re-poison the next sample (see finite_guard.h).
        if (!all_finite(vp_) || !all_finite(s_) || !std::isfinite(out)) {
            recover_nonfinite();
            out = 0.0;
        }
        return static_cast<float>(out);
    }

    // ---- state -----------------------------------------------------------------

    static constexpr int    kRefinePasses = TOCAKE_DK_REFINE;
    static constexpr int    kDcMaxIters = 400;
    static constexpr double kDcTol = 1.0e-12;   // V

    // The potential kDkVcc names: the battery when the circuit models its
    // supply network, else the variant's (static) rail.
    double supply() const noexcept {
        if constexpr (dk_circuit_supply<Circuit>::has) return Circuit::kSupply;
        else return v_.vcc;
    }
    const CakeVariant& v_;
    DkKnobs k_;
    double fs_ = 96000.0;
    double tol_ = dk_circuit_tol<Circuit>::value;
    // The base of the rank updates: the last EXACT rebuild's factorisation,
    // solution and branch list.
    static constexpr int kMR = 2 + S + P;
    static constexpr int kMaxRank = 6;           // three pots, two legs each
    static constexpr double kRankMaxG = 1.0e-2;  // S: a leg under 100 ohm is rebuilt
    std::array<double, NN * NN> lu0_{};
    std::array<int, NN> perm0_{};
    std::array<double, NN * kMR> xm0_{};
    std::array<Branch, kMaxBranches> br0_{};
    int nbr0_ = 0;
    double base_fs_ = 0.0;
    bool base_ok_ = false;

    // X = X0 - W (I + D U^T W)^-1 D U^T X0, W = G0^-1 U, over the branches whose
    // conductance differs from the base. False (and a full rebuild) when the
    // branch list changed shape, too many changed, or one touches a known
    // voltage (its right-hand side would change too).
    bool rank_update(const Network& n, std::array<double, NN * kMR>& X) const noexcept {
        if (n.nbr != nbr0_) return false;
        int idx[kMaxRank]; double dg[kMaxRank]; int r = 0;
        for (int k = 0; k < n.nbr; ++k) {
            if (n.br[k].a != br0_[k].a || n.br[k].b != br0_[k].b) return false;
            if (n.br[k].g == br0_[k].g) continue;
            if (n.br[k].a < 0 && n.br[k].a != kDkGnd) return false;
            if (n.br[k].b < 0 && n.br[k].b != kDkGnd) return false;
            if (r == kMaxRank) return false;
            // A leg near a pot END (under 100 ohm) makes the update a near short,
            // and it loses the precision the refined rebuild keeps: rebuild.
            if (n.br[k].g > kRankMaxG || br0_[k].g > kRankMaxG) return false;
            idx[r] = k; dg[r] = n.br[k].g - br0_[k].g; ++r;
        }
        X = xm0_;
        if (r == 0) return true;
        double W[kMaxRank][NN];
        for (int c = 0; c < r; ++c) {
            double u[NN] = {};
            const int a = n.br[idx[c]].a, b = n.br[idx[c]].b;
            if (a >= 0) u[a] += 1.0;
            if (b >= 0) u[b] -= 1.0;
            lu_pivot_solve(lu0_, perm0_, u, W[c]);
            // One refinement pass against the BASE network, branch by branch,
            // for the same reason the full rebuild refines: the pot ends.
            double res[NN], dw[NN];
            for (int k = 0; k < NN; ++k) res[k] = u[k];
            for (int k = 0; k < nbr0_; ++k) {
                const int p = br0_[k].a, q = br0_[k].b;
                const double cur = br0_[k].g * ((p >= 0 ? W[c][p] : 0.0) - (q >= 0 ? W[c][q] : 0.0));
                if (p >= 0) res[p] -= cur;
                if (q >= 0) res[q] += cur;
            }
            lu_pivot_solve(lu0_, perm0_, res, dw);
            for (int k = 0; k < NN; ++k) W[c][k] += dw[k];
        }
        auto ut = [&](int c, const double* v) {       // u_c^T v
            const int a = n.br[idx[c]].a, b = n.br[idx[c]].b;
            return (a >= 0 ? v[a] : 0.0) - (b >= 0 ? v[b] : 0.0);
        };
        // A = I + D U^T W (r x r), T = D U^T X0 (r x M); solve A Z = T in place.
        double A[kMaxRank][kMaxRank], T[kMaxRank][kMR];
        for (int i = 0; i < r; ++i) {
            for (int j = 0; j < r; ++j) A[i][j] = (i == j ? 1.0 : 0.0) + dg[i] * ut(i, W[j]);
            const int a = n.br[idx[i]].a, b = n.br[idx[i]].b;
            for (int m = 0; m < kMR; ++m)
                T[i][m] = dg[i] * ((a >= 0 ? xm0_[a * kMR + m] : 0.0) - (b >= 0 ? xm0_[b * kMR + m] : 0.0));
        }
        for (int c = 0; c < r; ++c) {                 // Gaussian elimination, partial pivoting
            int piv = c;
            for (int i = c + 1; i < r; ++i) if (std::fabs(A[i][c]) > std::fabs(A[piv][c])) piv = i;
            if (!(std::fabs(A[piv][c]) > 0.0)) return false;
            if (piv != c) { for (int j = 0; j < r; ++j) std::swap(A[c][j], A[piv][j]);
                            for (int m = 0; m < kMR; ++m) std::swap(T[c][m], T[piv][m]); }
            for (int i = c + 1; i < r; ++i) {
                const double f = A[i][c] / A[c][c];
                for (int j = c; j < r; ++j) A[i][j] -= f * A[c][j];
                for (int m = 0; m < kMR; ++m) T[i][m] -= f * T[c][m];
            }
        }
        for (int c = r - 1; c >= 0; --c) {
            for (int m = 0; m < kMR; ++m) {
                double acc = T[c][m];
                for (int j = c + 1; j < r; ++j) acc -= A[c][j] * T[j][m];
                T[c][m] = acc / A[c][c];
            }
        }
        for (int k = 0; k < NN; ++k)
            for (int c = 0; c < r; ++c) {
                const double w = W[c][k];
                for (int m = 0; m < kMR; ++m) X[k * kMR + m] -= w * T[c][m];
            }
        for (double v : X) if (!std::isfinite(v)) return false;
        return true;
    }
    static constexpr std::array<bool, P> make_tight_mask() noexcept {
        std::array<bool, P> m{};
        if constexpr (dk_circuit_tight<Circuit>::has) for (int p : Circuit::kTightPorts) m[p] = true;
        return m;
    }
    static constexpr std::array<bool, P> kTightMask = make_tight_mask();
    bool dirty_ = true;
    bool coarse_ = false;
    bool dc_seed_settled_ = false;

    double inv_vt_ = 0.0;
    std::array<double, NB> bjt_is_{}, inv_bf_{}, inv_vaf_{};
    std::array<double, ND> is_d_{}, inv_nvt_{};
    std::array<double, (NS > 0 ? NS : 1)> is_s_{}, inv_nvt_s_{};
    std::array<double, kJunctions> jvt_{}, jvcrit_{};

    std::array<double, S> gc_{}, g2_{};
    std::array<double, P * P> K_{};
    std::array<double, P> Pc_{}, Pu_{};
    std::array<double, P * S> Ps_{};
    std::array<double, S> Xc_{}, Xu_{};
    std::array<double, S * S> Xs_{};
    std::array<double, S * P> Xi_{};
    std::array<double, NT> Oc_{}, Ou_{};
    std::array<double, NT * S> Os_{};
    std::array<double, NT * P> Oi_{};
    // Filled only when kTransposed; a few hundred doubles otherwise idle.
    std::array<double, S * P> PsT_{};   // column k of Ps, contiguous
    std::array<double, P * P> KT_{};    // column c of K
    std::array<double, S * S> XsT_{};   // column k of Xs
    std::array<double, P * S> XiT_{};   // column c of Xi

    std::array<double, NN> x_rest_{};
    std::array<double, NT> tap_rest_{};
    std::array<double, S> s_{};
    std::array<double, P> vp_{}, vp1_{}, vp2_{};

    long dbg_iters_ = 0, dbg_calls_ = 0, dbg_nonconv_ = 0, dbg_nonfinite_ = 0, dbg_rebuilds_ = 0;
    bool branch_overflow_ = false;
    bool bad_part_ = false;
    bool rest_stale_ = false;   // a knob moved since the last DC rest: the non-finite guard re-solves it first
    int  dbg_max_iters_ = 0;
    double dbg_max_step_ = 0.0;
};

}  // namespace tocake
