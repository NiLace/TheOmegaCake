// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// The Lost Pumpkin -- the two-dual-op-amp version of the circuit, from a late-1970s workshop drawing. Engine and
// component values.
//
// SOURCE OF TRUTH: the traced production schematic, in two independent drawings that agree on every audio-path value:
//   * a KiCad export of the board (values confirmed against a production unit);
//   * the production build document, page 3 (schematic + BOM).
// They differ in one part of the audio path: the input bias resistor, 133K (KiCad) against 130K (BOM), which moves the
// input high-pass from 12.0 to 12.2 Hz. The KiCad value is used.
//
// TOPOLOGY. Input: C2 100n into U1A(+), biased by R1 133K from Vref (1M and 47p to ground at the jack). Three
// non-inverting clipping stages (U1A, U1B, U3B), each with feedback 470K || 470p || a 1N4148 anti-parallel pair and a
// ground leg R + C to ground (stage 1: 4K7 + the SUSTAIN pot as a rheostat + 4u7; stages 2 and 3: 47K + 4u7), and a
// 4K7 / 4n7 low-pass after each. Tone: from the third low-pass's node, C12 3n3 + R18 27K to ground (treble arm, pot
// pin 3) and R17 27K + C13 10n to ground then C14 100n in SERIES (bass arm, pot pin 1); the wiper through C15 100n to
// U3A(+), biased by R19 240K from Vref. U3A gain 1 + R21/R20 (10K/10K) with C16 47u in the gain leg, then C17 10u into
// the VOLUME pot (10K to ground), the wiper the output. Vref: a TL071 buffer of a 100K/100K divider of the 9 V supply.
//
// THE OP-AMPS ARE THE REAL 4558. Each of the four halves is a data-sheet model of the NJM4558
// (src/njm4558_generated.h): an NPN input pair with emitter degeneration and a tail current (the input bias current
// and its dynamics), a two-stage Miller-compensated gain (GBW 3 MHz, slew 1 V/us, 100 dB), the output clamp at the
// sheet's load-dependent headroom, and RO. Ideal op-amps (hard rails at the supply, no input current, infinite
// bandwidth) are audibly different: the gain-bandwidth at small levels, the swing at large.
//
// ENGINE. Four blocks, each one op-amp with everything linear around it, solved as a nodal network reduced onto its
// junction PORTS: the linear part is
// factorised when the rate or a knob changes, and each sample runs a Newton on the ports only (five per clipping
// block: the pair's two base-emitter junctions, the two clamp diodes, the feedback diode pair; four in the output
// block). The blocks meet at the op-amps' (+) inputs: block k computes the node that drives block k+1's (+), and block
// k+1's base current flows back into that node with ONE SAMPLE OF LAG (a current of tens of nA into a node of a few
// kohm; TOCAKE_LP_SWEEPS > 1 iterates it). Block 1 holds the input
// network, block 3 the tone network, block 4 the output stage, C17 and the volume pot. Capacitors are trapezoidal,
// except the op-amps' compensation capacitor C2 (TOCAKE_LP_C2_RULE): its closed-loop mode sits far above Nyquist once
// the feedback diodes conduct, where the trapezoid does not damp.
// The rest point is solved on the whole nodal network with source stepping, once, in the constructor (not real time):
// it does not depend on the rate (the capacitors are open) nor on the knobs (every pot sits behind a coupling
// capacitor), so it is solved with the pots at their construction positions and prepare() does not repeat it (it
// would cost about 1.4 ms, longer than an audio block, on a model change or an HQ toggle). reset() restores it.
//
// THE NEWTON'S FORM. Every source is affine in one function of its port (exp(v / V_T), or the feedback
// pair's e - 1/e), so the Newton runs on the ports' functions, not on the sources. A clamp diode
// more than 40 V_T into reverse drops out of the system, and while both clamps of a block are off (every sample at
// playing levels) the Newton is a fixed-size 3 x 3 (2 x 2 in the output block) solved by its adjugate. The stop is
// 1e-4 of the port voltage (TOCAKE_LP_TOL). Measured against
// the per-source form with a 1e-9 stop, on 80 sines (0.01 to 3 V, 100 Hz to 8 kHz, four knob settings, 96 kHz): worst
// peak residual -108.9 dB re the peak output.
//
// IDEALISATION, the same in the deck: Vref exactly half the supply; the reverse-polarity diode and the supply filters
// drop nothing; the diodes the repository's D1N914 card. The pair's base-collector junctions are held at their reverse
// limit (their current is the constant IS/BR, 8e-16 A): the bases stay more than a volt under the collectors (V+ less
// RC times 15 uA) for any input below about 4 V peak.

#pragma once

#include "bjt_dc_bias.h"   // V_T_300K: the reference simulator's thermal voltage at 27 C
#include "isa_hot.h"
#include "njm4558_generated.h"

#ifndef TOCAKE_LP_SWEEPS
#define TOCAKE_LP_SWEEPS 1       // passes over the four blocks per sample; > 1 iterates the lagged base currents
#endif
#ifndef TOCAKE_LP_C2_RULE
#define TOCAKE_LP_C2_RULE 2      // the compensation capacitor's rule: 0 trapezoid, 1 backward Euler, 2 BDF2
#endif
#ifndef TOCAKE_LP_TOL
#define TOCAKE_LP_TOL 1e-4       // Newton stop: every port step below TOCAKE_LP_TOL * (1e-3 + max|p|) volts
#endif

#include <cmath>
#include <cstdint>
#include <cstring>
#include <utility>

namespace tocake {

struct LostPumpkinParams {
    const char* name;
    double vcc;                              // supply (V); Vref = vcc / 2, the op-amps between vcc and ground
    double C_in, R_bias;                     // C2 100n, R1 133K (to Vref)
    double R_sus_fixed, R_sustain, C_g1;     // R3 4K7, SUSTAIN B100K (rheostat), C4 4u7
    double R_g23, C_g23;                     // R9 / R14 47K, C6 / C9 4u7
    double R_f, C_f;                         // R4 / R10 / R15 470K, C3 / C7 / C10 470p
    double R_lp, C_lp;                       // R7 / R11 / R16 4K7, C5 / C8 / C11 4n7
    double C_treble, R_treble;               // C12 3n3, R18 27K
    double R_bass, C_bass, C_bass_series;    // R17 27K, C13 10n, C14 100n
    double R_tone;                           // TONE B100K
    double C_out_couple, R_out_bias;         // C15 100n, R19 240K (to Vref)
    double R_gain_leg, C_gain_leg, R_gain_fb;// R20 10K, C16 47u, R21 10K
    double C_vol, R_volume;                  // C17 10u, VOLUME B10K
    double diode_Is, diode_N;                // 1N4148, the repository's D1N914 card
};

inline constexpr LostPumpkinParams TheLostPumpkin = {
    "Lost Pumpkin", 9.0,
    100e-9, 133e3,
    4.7e3, 100e3, 4.7e-6,
    47e3, 4.7e-6,
    470e3, 470e-12,
    4.7e3, 4.7e-9,
    3.3e-9, 27e3,
    27e3, 10e-9, 100e-9,
    100e3,
    100e-9, 240e3,
    10e3, 47e-6, 10e3,
    10e-6, 10e3,
    2.52e-9, 1.752,
};

enum class LostTap { Stage1, Stage2, Stage3, Out };

namespace lp_detail {

constexpr int MAXN = 18;   // unknown nodes of the largest block (the third: the op-amp, its stage and the tone network)
constexpr int MAXC = 8;    // capacitors of a block
constexpr int MAXD = 3 + MAXC;   // drivers: 1, the input voltage, the lagged current, each capacitor's history
constexpr int MAXS = 7;    // nonlinear current sources
constexpr int MAXP = 5;    // ports
constexpr int MAXW = MAXC + 2;   // read-outs: each capacitor's voltage, the op-amp's output, the voltage to the next block
constexpr int D_ONE = 0, D_UIN = 1, D_LAG = 2, D_CAP = 3;
// Fixed nodes (negative references): ground, the supply, Vref (= the op-amps' mid-supply), the clamps' two offset
// nodes, and the block's input voltage.
constexpr int GND = -1, VCC = -2, VREF = -3, CHI = -4, CLO = -5, UIN = -6;
enum SrcKind { BJT_B, BJT_C, DIODE, DPAIR };

// SPICE3's junction limiting (forward steps only).
inline double pnjlim(double vnew, double vold, double nvt, double vcrit) noexcept {
    if (vnew > vcrit && std::fabs(vnew - vold) > 2.0 * nvt) {
        if (vold > 0.0) {
            const double arg = 1.0 + (vnew - vold) / nvt;
            return arg > 0.0 ? vold + nvt * std::log(arg) : vcrit;
        }
        return nvt * std::log(vnew / nvt);
    }
    return vnew;
}

// In-place LU with partial pivoting of the n x n matrix a; false = singular. Not real time (factorisations only).
inline bool lu(double (&a)[MAXN][MAXN], int n, int (&piv)[MAXN]) noexcept {
    for (int k = 0; k < n; ++k) {
        int p = k;
        for (int i = k + 1; i < n; ++i) if (std::fabs(a[i][k]) > std::fabs(a[p][k])) p = i;
        if (a[p][k] == 0.0) return false;
        piv[k] = p;
        if (p != k) for (int j = 0; j < n; ++j) std::swap(a[k][j], a[p][j]);
        const double inv = 1.0 / a[k][k];
        for (int i = k + 1; i < n; ++i) {
            const double f = (a[i][k] *= inv);
            if (f != 0.0) for (int j = k + 1; j < n; ++j) a[i][j] -= f * a[k][j];
        }
    }
    return true;
}
// Every swap applied to b first, then the substitutions (lu() swaps whole rows, multipliers included).
inline void lu_solve(const double (&a)[MAXN][MAXN], int n, const int (&piv)[MAXN], double* b) noexcept {
    for (int k = 0; k < n; ++k) if (piv[k] != k) std::swap(b[k], b[piv[k]]);
    for (int k = 0; k < n; ++k) for (int i = k + 1; i < n; ++i) b[i] -= a[i][k] * b[k];
    for (int k = n - 1; k >= 0; --k) {
        for (int j = k + 1; j < n; ++j) b[k] -= a[k][j] * b[j];
        b[k] /= a[k][k];
    }
}

// One block: its stamps (kept, so a knob rebuilds without re-describing the circuit), its reduction onto the ports and
// read-outs, and its state. Read-outs 0 .. nc-1 are the capacitors' voltages, in the order they were added.
struct Block {
    // ---- the description ---------------------------------------------------------------------------------------
    int n = 0, nc = 0, ns = 0, np = 0, nw = 0;
    struct R { int a, b; double g; };
    struct G { int from, to, cp, cn; double gm; };
    struct C { int a, b; double c; int rule; };
    R res[32]; int nres = 0;
    G vccs[8]; int nvccs = 0;
    C cap[MAXC];
    int iee_node = -1; double iee = 0.0;      // a constant current leaving this node to ground
    int lag_node = -1;                        // the lagged current leaves this node
    int src_from[MAXS], src_to[MAXS], src_port[MAXS], src_kind[MAXS];
    int port_a[MAXP], port_b[MAXP], port_kind[MAXP];
    int wa[MAXW], wb[MAXW];                   // read-out k = v(wa) - v(wb)
    int out_w = -1, next_w = -1;              // which read-outs are the op-amp's output and the next block's drive
    int ib2_src = -1;                         // the (+) input's base current (for the block upstream)
    double fixed[5] = {};                     // GND, VCC, VREF, CHI, CLO (index -ref - 1); UIN is a driver

    // ---- the reduction (rebuilt by factor()) -----------------------------------------------------------------
    double E[MAXP][MAXD] = {}, K[MAXP][MAXS] = {}, Fw[MAXW][MAXD] = {}, Kw[MAXW][MAXS] = {};
    double gc[MAXC] = {}, hc[MAXC] = {};      // the companion conductance, and C / T for the histories
    // the same reduction onto one function per port (rebuilt with it): ports p = Ef u + A phi(p), read-outs
    // w = Fwf u + Aw phi(p). Ef and Fwf are stored together by driver (row d of EwT: Ef's column d, then Fwf's, the
    // unused entries zero), Aw by port, so each product is a run of whole fixed-length columns: one pass over the
    // drivers gives the Newton's start p0 and the read-outs' driver part together.
    static constexpr int NEW = 16;   // MAXP + MAXW, rounded up to eight vectors of two
    alignas(32) double EwT[MAXD][NEW] = {};   // [driver][port r, then read-out MAXP + r]
    alignas(32) double AwT[MAXP][NEW] = {};   // [port][read-out MAXP + r]; the first MAXP entries stay 0
    double A[MAXP][MAXP] = {};

    // ---- the state ---------------------------------------------------------------------------------------------
    double v1[MAXC] = {}, v2[MAXC] = {}, i1[MAXC] = {};   // each capacitor: last voltage, the one before, last current
    double p[MAXP] = {};                                   // the ports at the last sample (Newton's start)
    double w[MAXW] = {};                                   // the read-outs of the last solve
    double isrc[MAXS] = {};                                // the sources' currents of the last solve
    // the rest point, restored by reset()
    double r_v1[MAXC] = {}, r_p[MAXP] = {}, r_w[MAXW] = {}, r_isrc[MAXS] = {};

    // ---- description helpers -----------------------------------------------------------------------------------
    void add_r(int a, int b, double r) noexcept { res[nres++] = R{a, b, 1.0 / r}; }
    void add_g(int from, int to, int cp, int cn, double gm) noexcept { vccs[nvccs++] = G{from, to, cp, cn, gm}; }
    int add_c(int a, int b, double c, int rule = 0) noexcept { cap[nc] = C{a, b, c, rule}; return nc++; }
    int add_port(int a, int b, int kind) noexcept { port_a[np] = a; port_b[np] = b; port_kind[np] = kind; return np++; }
    int add_src(int from, int to, int port, int kind) noexcept {
        src_from[ns] = from; src_to[ns] = to; src_port[ns] = port; src_kind[ns] = kind; return ns++;
    }
    // the extra read-outs come after the capacitors' (every capacitor is added before the first extra one)
    int add_w(int a, int b) noexcept {
        if (nw < nc) { for (int k = 0; k < nc; ++k) { wa[k] = cap[k].a; wb[k] = cap[k].b; } nw = nc; }
        wa[nw] = a; wb[nw] = b; return nw++;
    }

    double fixed_v(int ref) const noexcept { return fixed[-ref - 1]; }   // GND..CLO; UIN handled by the caller

    // The network matrix G and the driver matrix B (rhs = B u), the sources' incidence N (rhs += N i), the ports' and
    // read-outs' maps from the nodes (Pv, Wv) and the drivers (Pu, Wu). dc: capacitors open, no histories. KCL rows:
    // the currents LEAVING a node through the elements equal the currents injected.
    void build(bool dc, double (&Gm)[MAXN][MAXN], double (&B)[MAXN][MAXD], double (&N)[MAXN][MAXS],
               double (&Pv)[MAXP][MAXN], double (&Pu)[MAXP][MAXD], double (&Wv)[MAXW][MAXN],
               double (&Wu)[MAXW][MAXD]) const noexcept {
        std::memset(Gm, 0, sizeof Gm); std::memset(B, 0, sizeof B); std::memset(N, 0, sizeof N);
        std::memset(Pv, 0, sizeof Pv); std::memset(Pu, 0, sizeof Pu);
        std::memset(Wv, 0, sizeof Wv); std::memset(Wu, 0, sizeof Wu);
        // a fixed node's value as (constant, input coefficient)
        auto fx = [&](int ref, double& c, double& u) { c = 0.0; u = 0.0; if (ref == UIN) u = 1.0; else c = fixed_v(ref); };
        // conductance g between a and b
        auto cond = [&](int a, int b, double g) {
            for (int s = 0; s < 2; ++s) {
                const int x = s ? b : a, y = s ? a : b;
                if (x < 0) continue;
                Gm[x][x] += g;
                if (y >= 0) Gm[x][y] -= g;
                else { double c, u; fx(y, c, u); B[x][D_ONE] += g * c; B[x][D_UIN] += g * u; }
            }
        };
        for (int k = 0; k < nres; ++k) cond(res[k].a, res[k].b, res[k].g);
        for (int k = 0; k < nc; ++k) {
            if (dc) continue;
            cond(cap[k].a, cap[k].b, gc[k]);
            if (cap[k].a >= 0) B[cap[k].a][D_CAP + k] += 1.0;      // i = gc v - J: J enters a, leaves b
            if (cap[k].b >= 0) B[cap[k].b][D_CAP + k] -= 1.0;
        }
        // a VCCS: gm (v(cp) - v(cn)) leaves `from` and enters `to`
        for (int k = 0; k < nvccs; ++k) {
            const G& e = vccs[k];
            for (int s = 0; s < 2; ++s) {
                const int row = s ? e.to : e.from;
                const double sg = s ? -e.gm : e.gm;
                if (row < 0) continue;
                for (int t = 0; t < 2; ++t) {
                    const int col = t ? e.cn : e.cp;
                    const double sc = t ? -sg : sg;
                    if (col >= 0) Gm[row][col] += sc;
                    else { double c, u; fx(col, c, u); B[row][D_ONE] -= sc * c; B[row][D_UIN] -= sc * u; }
                }
            }
        }
        if (iee_node >= 0) B[iee_node][D_ONE] -= iee;
        if (lag_node >= 0) B[lag_node][D_LAG] -= 1.0;
        for (int s = 0; s < ns; ++s) {
            if (src_from[s] >= 0) N[src_from[s]][s] -= 1.0;
            if (src_to[s] >= 0)   N[src_to[s]][s]   += 1.0;
        }
        auto diff = [&](int a, int b, double* rv, double* ru) {
            for (int t = 0; t < 2; ++t) {
                const int x = t ? b : a;
                const double sg = t ? -1.0 : 1.0;
                if (x >= 0) rv[x] += sg;
                else { double c, u; fx(x, c, u); ru[D_ONE] += sg * c; ru[D_UIN] += sg * u; }
            }
        };
        for (int k = 0; k < np; ++k) diff(port_a[k], port_b[k], Pv[k], Pu[k]);
        for (int k = 0; k < nw; ++k) diff(wa[k], wb[k], Wv[k], Wu[k]);
    }

    // The reduction at the sample rate: E = Pv G^-1 B + Pu, K = Pv G^-1 N, and the same for the read-outs.
    bool factor() noexcept {
        double Gm[MAXN][MAXN], B[MAXN][MAXD], N[MAXN][MAXS], Pv[MAXP][MAXN], Pu[MAXP][MAXD], Wv[MAXW][MAXN], Wu[MAXW][MAXD];
        build(false, Gm, B, N, Pv, Pu, Wv, Wu);
        int piv[MAXN];
        if (!lu(Gm, n, piv)) return false;
        const int nd = D_CAP + nc;
        double x[MAXN];
        for (int col = 0; col < nd + ns; ++col) {
            for (int i = 0; i < n; ++i) x[i] = col < nd ? B[i][col] : N[i][col - nd];
            lu_solve(Gm, n, piv, x);
            for (int r = 0; r < np; ++r) {
                double s = 0.0;
                for (int i = 0; i < n; ++i) s += Pv[r][i] * x[i];
                if (col < nd) E[r][col] = s + Pu[r][col]; else K[r][col - nd] = s;
            }
            for (int r = 0; r < nw; ++r) {
                double s = 0.0;
                for (int i = 0; i < n; ++i) s += Wv[r][i] * x[i];
                if (col < nd) Fw[r][col] = s + Wu[r][col]; else Kw[r][col - nd] = s;
            }
        }
        return true;
    }
};

}  // namespace lp_detail

class LostPumpkin {
    using Block = lp_detail::Block;
public:
    explicit LostPumpkin(const LostPumpkinParams& p) noexcept : p_(p) {
        describe();
        apply_sustain(); apply_tone();   // the pots at their construction positions (the rest does not depend on them)
        solve_rest();
    }

    void prepare(double sample_rate) noexcept {
        fs_ = sample_rate;
        const double T = 1.0 / fs_;
        for (Block& b : blk_)
            for (int k = 0; k < b.nc; ++k) {
                const double c = b.cap[k].c;
                b.hc[k] = c / T;
                b.gc[k] = b.cap[k].rule == 0 ? 2.0 * c / T : b.cap[k].rule == 1 ? c / T : 1.5 * c / T;
            }
        prepared_ = true;
        apply_sustain(); apply_tone(); apply_volume();
        for (Block& b : blk_) refactor(b);
        reset();
    }

    // Every state back to the rest point (solved in the constructor; it depends on neither the rate nor the knobs).
    void reset() noexcept {
        for (Block& b : blk_) {
            for (int k = 0; k < b.nc; ++k) { b.v1[k] = b.v2[k] = b.r_v1[k]; b.i1[k] = 0.0; }
            std::memcpy(b.p, b.r_p, sizeof b.p);
            std::memcpy(b.w, b.r_w, sizeof b.w);
            std::memcpy(b.isrc, b.r_isrc, sizeof b.isrc);
        }
        for (int k = 0; k < 3; ++k) lag_[k] = blk_[k + 1].r_isrc[blk_[k + 1].ib2_src];
    }

    // Sustain: the pot is a rheostat between pin 3 and the wiper (pin 1 open); 1 = the wiper at pin 3 = most gain.
    // A setter that receives its current value does nothing (the audio must not depend on the host's block size).
    void set_sustain(double frac) noexcept {
        const double s = clamp01(frac);
        if (s == sus_ && prepared_) return;
        sus_ = s;
        if (prepared_) { apply_sustain(); refactor(blk_[0]); }
    }
    // Tone: 1 = the wiper at pin 3, the treble arm.
    void set_tone(double frac) noexcept {
        const double t = clamp01(frac);
        if (t == tone_ && prepared_) return;
        tone_ = t;
        if (prepared_) { apply_tone(); refactor(blk_[2]); }
    }
    // Volume: 1 = the wiper at the top (most output). The pot's whole 10K loads C17 at any position, so the network
    // does not change: only the divider does.
    void set_volume(double frac) noexcept {
        vol_ = clamp01(frac);
        apply_volume();
    }

    TOCAKE_SPAN_BODY float process_sample(float vin, LostTap tap = LostTap::Out) noexcept {
        double u[4][lp_detail::MAXD];
        double drive = vin;
        for (int sweep = 0; sweep < TOCAKE_LP_SWEEPS; ++sweep) {
            drive = vin;
            for (int k = 0; k < 4; ++k) {
                Block& b = blk_[k];
                load_drivers(b, drive, k < 3 ? lag_[k] : 0.0, u[k]);
                solve_ports(b, u[k]);
                drive = b.w[b.next_w];
            }
            if (TOCAKE_LP_SWEEPS > 1)   // the lagged currents from this pass, for the next one
                for (int k = 0; k < 3; ++k) lag_[k] = blk_[k + 1].isrc[blk_[k + 1].ib2_src];
        }
        for (int k = 0; k < 4; ++k) commit(blk_[k]);
        for (int k = 0; k < 3; ++k) lag_[k] = blk_[k + 1].isrc[blk_[k + 1].ib2_src];
        const double vr = 0.5 * p_.vcc;
        switch (tap) {
        case LostTap::Stage1: return guard(blk_[0].w[blk_[0].out_w] - vr);
        case LostTap::Stage2: return guard(blk_[1].w[blk_[1].out_w] - vr);
        case LostTap::Stage3: return guard(blk_[2].w[blk_[2].out_w] - vr);
        default: break;
        }
        return guard(drive * vol_div_);
    }

    void span_run(float* buf, uint32_t n) noexcept {
        for (uint32_t j = 0; j < n; ++j) buf[j] = process_sample(buf[j], LostTap::Out);
    }
    TOCAKE_SPAN_RUN_ISA void span_run_isa_v3(float* buf, uint32_t n) noexcept {
        for (uint32_t j = 0; j < n; ++j) buf[j] = process_sample(buf[j], LostTap::Out);
    }

    long nonfinite_resets() const noexcept { return dbg_nonfinite_; }
    const char* name() const noexcept { return p_.name; }
    // The port Newton's counters (always kept: a few adds per block solve): block solves, iterations, and solves
    // that ran out of iterations.
    long dbg_solves() const noexcept { return dbg_solves_; }
    long dbg_iters()  const noexcept { return dbg_iters_; }
    long dbg_capped() const noexcept { return dbg_capped_; }
    // the largest last Newton step (V) of any solve: under the stop tolerance unless a solve was capped
    double dbg_worst_step() const noexcept { return dbg_worst_step_; }
    bool rest_ok() const noexcept { return rest_ok_; }
    // The rest point's op-amp outputs (V, absolute), for the record.
    double rest_out(int k) const noexcept { return blk_[k].r_w[blk_[k].out_w]; }

private:
    static constexpr int MAXIT = 60;
    static double clamp01(double f) noexcept { return f < 0.0 ? 0.0 : (f > 1.0 ? 1.0 : f); }

    // ---- the circuit -------------------------------------------------------------------------------------------
    // One op-amp into block b, the deck's subcircuit element for element: (+) at `pin`, (-) at `nin`, output at
    // `out`; its seven internal nodes are the block's next. Sources Ib1, Ic1, Ib2, Ic2, the two clamps; ports vbe1,
    // vbe2, the two clamps.
    void add_opamp(Block& b, int pin, int nin, int out) noexcept {
        namespace O = njm4558;
        using namespace lp_detail;
        const int c1 = b.n++, c2 = b.n++, e1 = b.n++, e2 = b.n++, t = b.n++, n6 = b.n++, n7 = b.n++;
        b.add_r(VCC, c1, O::RC); b.add_r(VCC, c2, O::RC);
        b.add_r(e1, t, O::RE);   b.add_r(e2, t, O::RE);
        b.iee_node = t; b.iee = O::IEE;           // IEE 10 4 (V- is ground)
        b.add_r(t, VREF, O::REE);
        b.add_g(n6, VREF, c1, c2, O::GA);         // GA 6 99 11 12
        b.add_g(n6, VREF, t, VREF, O::GCM);       // GCM 6 99 10 99
        b.add_r(n6, VREF, O::R2);
        b.add_g(n7, VREF, n6, VREF, O::GB);       // GB 7 99 6 99
        b.add_r(n7, VREF, O::RO2);
        b.add_c(n6, n7, O::C2, TOCAKE_LP_C2_RULE);
        // EOUT + RO: the buffer copies v7 through RO. As a Norton, (v7 - v_out) / RO enters `out` and nothing is drawn
        // from node 7: two VCCS, gRO (v7 - vref) and gRO (vref - v_out).
        b.add_g(VREF, out, n7, VREF, 1.0 / O::RO);
        b.add_g(VREF, out, VREF, out, 1.0 / O::RO);
        const int q1 = b.add_port(nin, e1, BJT_B), q2 = b.add_port(pin, e2, BJT_B);
        b.add_src(nin, e1, q1, BJT_B); b.add_src(c1, e1, q1, BJT_C);
        b.ib2_src = b.add_src(pin, e2, q2, BJT_B); b.add_src(c2, e2, q2, BJT_C);
        const int pc = b.add_port(n7, CHI, DIODE), pe = b.add_port(CLO, n7, DIODE);
        b.add_src(n7, CHI, pc, DIODE); b.add_src(CLO, n7, pe, DIODE);
    }

    void describe() noexcept {
        using namespace lp_detail;
        const double vcc = p_.vcc;
        for (Block& b : blk_) {
            b.fixed[0] = 0.0; b.fixed[1] = vcc; b.fixed[2] = 0.5 * vcc;
            b.fixed[3] = vcc - njm4558::VC; b.fixed[4] = 0.0 + njm4558::VE;
        }
        // Blocks 1-3: a clipping stage. Nodes: m (-), out, the ground leg's middle; the op-amp's seven; the next node.
        for (int k = 0; k < 3; ++k) {
            Block& b = blk_[k];
            int pin = UIN;
            if (k == 0) {                       // the input network: C2 from the jack (the input) to (+), R1 to Vref
                pin = b.n++;
                b.add_c(UIN, pin, p_.C_in);
                b.add_r(pin, VREF, p_.R_bias);
            }
            const int m = b.n++, out = b.n++, leg = b.n++;
            add_opamp(b, pin, m, out);
            b.add_r(m, out, p_.R_f);
            b.add_c(m, out, p_.C_f);
            const int pf = b.add_port(m, out, DPAIR);
            b.add_src(m, out, pf, DPAIR);
            leg_r_[k] = b.nres;                 // stage 1's is re-set by the sustain knob
            b.add_r(m, leg, k == 0 ? p_.R_sus_fixed : p_.R_g23);
            b.add_c(leg, GND, k == 0 ? p_.C_g1 : p_.C_g23);
            const int nx = b.n++;               // the low-pass node (stages 1, 2) or the tone network's input (3)
            b.add_r(out, nx, p_.R_lp);
            b.add_c(nx, GND, p_.C_lp);
            if (k < 2) {
                b.lag_node = nx;
                b.out_w = b.add_w(out, GND);
                b.next_w = b.add_w(nx, GND);
            } else {                            // the tone network: G = nx, H, I, K, W, JN
                const int H = b.n++, I = b.n++, Kn = b.n++, W = b.n++, JN = b.n++;
                b.add_c(nx, H, p_.C_treble);
                b.add_r(H, GND, p_.R_treble);
                b.add_r(nx, I, p_.R_bass);
                b.add_c(I, GND, p_.C_bass);
                b.add_c(I, Kn, p_.C_bass_series);
                tone_r_ = b.nres;
                b.add_r(H, W, 1.0);             // pin 3 to the wiper (set by the knob)
                b.add_r(W, Kn, 1.0);            // the wiper to pin 1
                b.add_c(W, JN, p_.C_out_couple);
                b.add_r(JN, VREF, p_.R_out_bias);
                b.lag_node = JN;
                b.out_w = b.add_w(out, GND);
                b.next_w = b.add_w(JN, GND);
            }
        }
        // Block 4: U3A (+) driven by JN; gain leg R20 + C16, feedback R21; C17 into the volume pot's whole 10K.
        {
            Block& b = blk_[3];
            const int m = b.n++, out = b.n++, leg = b.n++, vt = b.n++;
            add_opamp(b, UIN, m, out);
            b.add_r(m, out, p_.R_gain_fb);
            b.add_r(m, leg, p_.R_gain_leg);
            b.add_c(leg, GND, p_.C_gain_leg);
            b.add_c(out, vt, p_.C_vol);
            b.add_r(vt, GND, p_.R_volume + 2e-6);   // top + bottom, each with the deck's 1 micro-ohm
            b.out_w = b.add_w(out, GND);
            b.next_w = b.add_w(vt, GND);
        }
    }

    void apply_sustain() noexcept {
        blk_[0].res[leg_r_[0]].g = 1.0 / (p_.R_sus_fixed + p_.R_sustain * (1.0 - sus_) + 1e-6);
    }
    // Each leg carries a 1 ohm end resistance, in the deck too. Not 1 micro-ohm as the other pots: inside the block
    // matrix 1e6 S beside the pair's 1e-8 S costs the reduction its last digits, and the rest point solved on the
    // unreduced network then misses the sample solve's by 41 uV at tone 0 (measured: a start thump; 1.9 nV with 1 ohm).
    void apply_tone() noexcept {
        blk_[2].res[tone_r_].g     = 1.0 / (p_.R_tone * (1.0 - tone_) + 1.0);
        blk_[2].res[tone_r_ + 1].g = 1.0 / (p_.R_tone * tone_ + 1.0);
    }
    void apply_volume() noexcept {
        const double top = p_.R_volume * (1.0 - vol_) + 1e-6, bot = p_.R_volume * vol_ + 1e-6;
        vol_div_ = bot / (top + bot);
    }

    // ---- the sample ----------------------------------------------------------------------------------------------
    // Each capacitor's history current J (i = gc v - J): trapezoid gc v1 + i1, backward Euler (C/T) v1, BDF2
    // (C/T)(2 v1 - v2 / 2).
    TOCAKE_SPAN_BODY static void load_drivers(const Block& b, double uin, double lag, double* u) noexcept {
        u[lp_detail::D_ONE] = 1.0; u[lp_detail::D_UIN] = uin; u[lp_detail::D_LAG] = lag;
        for (int k = 0; k < b.nc; ++k) {
            const int r = b.cap[k].rule;
            u[lp_detail::D_CAP + k] = r == 0 ? b.gc[k] * b.v1[k] + b.i1[k]
                                    : r == 1 ? b.hc[k] * b.v1[k]
                                             : b.hc[k] * (2.0 * b.v1[k] - 0.5 * b.v2[k]);
        }
    }

    // A source's current and its derivative against its port.
    TOCAKE_SPAN_BODY void source(int kind, double v, double& i, double& g) const noexcept {
        using namespace lp_detail;
        if (kind == DPAIR) {
            const double x = v * inv_nvt_d_;
            const double e = std::exp(x < 80.0 ? (x > -80.0 ? x : -80.0) : 80.0), ie = 1.0 / e;
            i = p_.diode_Is * (e - ie);
            g = p_.diode_Is * inv_nvt_d_ * (e + ie);
            return;
        }
        const double x = v * inv_vt_;
        const double e = std::exp(x < 80.0 ? x : 80.0);
        if (kind == DIODE) { i = njm4558::IS_D * (e - 1.0); g = njm4558::IS_D * inv_vt_ * e; return; }
        // the BJT with its base-collector junction at the reverse limit: Ib = IS/BF (e - 1) - IS/BR, Ic = IS e + IS/BR
        // (BR = 1, the default)
        constexpr double isb = njm4558::IS_Q / njm4558::BF;
        if (kind == BJT_B) { i = isb * (e - 1.0) - njm4558::IS_Q; g = isb * inv_vt_ * e; }
        else               { i = njm4558::IS_Q * (e + 1.0);       g = njm4558::IS_Q * inv_vt_ * e; }
    }

    // The real-time form of source(), split in two so each PORT pays one exp() however many sources share it (a BJT's
    // base and collector share its base-emitter junction): port_exp() per port, source_e() per source. A diode
    // further than 40 V_T into reverse takes e = 0: e^-40 = 4e-18 is below half an ulp of 1, so IS (e - 1) is the
    // same double, and its conductance (1e-33 S) moves no Jacobian entry. Most samples the output clamp's two diodes
    // sit hundreds of V_T into reverse.
    TOCAKE_SPAN_BODY void port_exp(int kind, double v, double& e, double& ie) const noexcept {
        using namespace lp_detail;
        if (kind == DPAIR) {
            const double x = v * inv_nvt_d_;
            e = std::exp(x < 80.0 ? (x > -80.0 ? x : -80.0) : 80.0);
            ie = 1.0 / e;
            return;
        }
        const double x = v * inv_vt_;
        ie = 0.0;
        e = (kind == DIODE && x < -40.0) ? 0.0 : std::exp(x < 80.0 ? x : 80.0);
    }
    TOCAKE_SPAN_BODY void source_e(int kind, double e, double ie, double& i, double& g) const noexcept {
        using namespace lp_detail;
        if (kind == DPAIR) { i = p_.diode_Is * (e - ie); g = p_.diode_Is * inv_nvt_d_ * (e + ie); return; }
        if (kind == DIODE) { i = njm4558::IS_D * (e - 1.0); g = njm4558::IS_D * inv_vt_ * e; return; }
        constexpr double isb = njm4558::IS_Q / njm4558::BF;
        if (kind == BJT_B) { i = isb * (e - 1.0) - njm4558::IS_Q; g = isb * inv_vt_ * e; }
        else               { i = njm4558::IS_Q * (e + 1.0);       g = njm4558::IS_Q * inv_vt_ * e; }
    }

    TOCAKE_SPAN_BODY double limit(int kind, double vnew, double vold) const noexcept {
        using namespace lp_detail;
        if (kind == DPAIR) {
            if (vnew >= 0.0) return pnjlim(vnew, vold > 0.0 ? vold : 0.0, nvt_d_, vcrit_d_);
            return -pnjlim(-vnew, vold < 0.0 ? -vold : 0.0, nvt_d_, vcrit_d_);
        }
        return pnjlim(vnew, vold, V_T_300K, kind == DIODE ? vcrit_dio_ : vcrit_q_);
    }

    // The block's ports by Newton from the last sample's, then its read-outs. Nothing is committed here.
    TOCAKE_SPAN_BODY void solve(Block& b, const double* u) noexcept {
        using namespace lp_detail;
        const int np = b.np, ns = b.ns, nd = D_CAP + b.nc;
        double p0[MAXP], i[MAXS], g[MAXS];
        for (int r = 0; r < np; ++r) {
            double s = 0.0;
            for (int d = 0; d < nd; ++d) s += b.E[r][d] * u[d];
            p0[r] = s;
        }
        double* p = b.p;
        ++dbg_solves_;
        int it = 0;
        bool done = false;
        double dlast = 0.0;
        for (; it < MAXIT; ++it) {
            double ex[MAXP], ix[MAXP];
            for (int k = 0; k < np; ++k) port_exp(b.port_kind[k], p[k], ex[k], ix[k]);
            for (int s = 0; s < ns; ++s) source_e(b.src_kind[s], ex[b.src_port[s]], ix[b.src_port[s]], i[s], g[s]);
            double J[MAXP][MAXP + 1];
            for (int r = 0; r < np; ++r) {
                double f = p[r] - p0[r];
                for (int c = 0; c < np; ++c) J[r][c] = (r == c) ? 1.0 : 0.0;
                for (int s = 0; s < ns; ++s) {
                    f -= b.K[r][s] * i[s];
                    J[r][b.src_port[s]] -= b.K[r][s] * g[s];
                }
                J[r][np] = -f;
            }
            // Gaussian elimination with partial pivoting on the np x np system
            for (int k = 0; k < np; ++k) {
                int q = k;
                for (int r = k + 1; r < np; ++r) if (std::fabs(J[r][k]) > std::fabs(J[q][k])) q = r;
                if (q != k) for (int c = k; c <= np; ++c) std::swap(J[k][c], J[q][c]);
                const double inv = 1.0 / J[k][k];
                for (int r = k + 1; r < np; ++r) {
                    const double f = J[r][k] * inv;
                    for (int c = k; c <= np; ++c) J[r][c] -= f * J[k][c];
                }
                J[k][k] = inv;
            }
            double dp[MAXP];
            for (int k = np - 1; k >= 0; --k) {
                double s = J[k][np];
                for (int c = k + 1; c < np; ++c) s -= J[k][c] * dp[c];
                dp[k] = s * J[k][k];
            }
            double dmax = 0.0, pmax = 0.0;
            for (int k = 0; k < np; ++k) {
                const double pn = limit(b.port_kind[k], p[k] + dp[k], p[k]);
                const double d = pn - p[k];
                if (std::fabs(d) > dmax) dmax = std::fabs(d);
                // the sources at the new point to first order (the step is below the tolerance at convergence)
                for (int s = 0; s < ns; ++s) if (b.src_port[s] == k) i[s] += g[s] * d;
                p[k] = pn;
                if (std::fabs(pn) > pmax) pmax = std::fabs(pn);
            }
            dlast = dmax;
            if (dmax < TOCAKE_LP_TOL * (1e-3 + pmax)) { done = true; ++it; break; }
        }
        dbg_iters_ += it;
        if (!done) ++dbg_capped_;
        if (dlast > dbg_worst_step_) dbg_worst_step_ = dlast;
        for (int s = 0; s < ns; ++s) b.isrc[s] = i[s];
        for (int r = 0; r < b.nw; ++r) {
            double s = 0.0;
            for (int d = 0; d < nd; ++d) s += b.Fw[r][d] * u[d];
            for (int q = 0; q < ns; ++q) s += b.Kw[r][q] * i[q];
            b.w[r] = s;
        }
    }

    // ---- the same Newton on one function per port ------------------------------------------
    // Every source is affine in its port's function phi: a BJT's base IS/BF (e - 1) - IS and collector IS (e + 1), a
    // clamp diode IS (e - 1), with phi = e = exp(v / V_T); the feedback pair Is (e - 1/e), with phi = e - 1/e. So the
    // K i term is A phi + a constant, the constant folded into the drivers' first column: a system in the ports
    // whatever the number of sources. A clamp diode with e = 0 (port_exp's reverse shortcut, most samples) has phi = 0
    // and a zero Jacobian column: its row does not feed the others, so the Newton eliminates only the ON ports (3 of 5
    // when both clamps are off) and the OFF ports follow by substitution. The same iterates as the per-source form, to
    // rounding.
    static void src_affine(int kind, double is_d, double& a, double& c) noexcept {
        using namespace lp_detail;
        constexpr double isb = njm4558::IS_Q / njm4558::BF;
        if (kind == BJT_B)      { a = isb; c = -(isb + njm4558::IS_Q); }
        else if (kind == BJT_C) { a = njm4558::IS_Q; c = njm4558::IS_Q; }
        else if (kind == DIODE) { a = njm4558::IS_D; c = -njm4558::IS_D; }
        else                    { a = is_d; c = 0.0; }
    }
    void refactor(Block& b) noexcept {
        using namespace lp_detail;
        b.factor();
        const int nd = D_CAP + b.nc;
        for (int r = 0; r < b.np; ++r) {
            for (int d = 0; d < nd; ++d) b.EwT[d][r] = b.E[r][d];
            for (int k = 0; k < b.np; ++k) b.A[r][k] = 0.0;
            for (int s = 0; s < b.ns; ++s) {
                double a, c; src_affine(b.src_kind[s], p_.diode_Is, a, c);
                b.A[r][b.src_port[s]] += b.K[r][s] * a;
                b.EwT[D_ONE][r] += b.K[r][s] * c;
            }
        }
        for (int r = 0; r < b.nw; ++r) {
            for (int d = 0; d < nd; ++d) b.EwT[d][MAXP + r] = b.Fw[r][d];
            for (int k = 0; k < b.np; ++k) b.AwT[k][MAXP + r] = 0.0;
            for (int s = 0; s < b.ns; ++s) {
                double a, c; src_affine(b.src_kind[s], p_.diode_Is, a, c);
                b.AwT[b.src_port[s]][MAXP + r] += b.Kw[r][s] * a;
                b.EwT[D_ONE][MAXP + r] += b.Kw[r][s] * c;
            }
        }
    }
    // phi and its derivative; a clamp diode past the reverse shortcut gives exactly 0 and 0 (the port is OFF)
    TOCAKE_SPAN_BODY void port_phi(int kind, double v, double& ph, double& dph) const noexcept {
        double e, ie;
        port_exp(kind, v, e, ie);
        if (kind == lp_detail::DPAIR) { ph = e - ie; dph = (e + ie) * inv_nvt_d_; }
        else                          { ph = e;      dph = e * inv_vt_; }
    }

    // The M x M system [J | r] by its adjugate: one division, no pivot search.
    template <int M>
    TOCAKE_SPAN_BODY static void cramer(const double (&J)[M][M + 1], double (&x)[M]) noexcept {
        if constexpr (M == 2) {
            const double inv = 1.0 / (J[0][0] * J[1][1] - J[0][1] * J[1][0]);
            x[0] = (J[1][1] * J[0][2] - J[0][1] * J[1][2]) * inv;
            x[1] = (J[0][0] * J[1][2] - J[1][0] * J[0][2]) * inv;
        } else {
            const double c00 = J[1][1] * J[2][2] - J[1][2] * J[2][1];
            const double c01 = J[1][2] * J[2][0] - J[1][0] * J[2][2];
            const double c02 = J[1][0] * J[2][1] - J[1][1] * J[2][0];
            const double inv = 1.0 / (J[0][0] * c00 + J[0][1] * c01 + J[0][2] * c02);
            const double r0 = J[0][3], r1 = J[1][3], r2 = J[2][3];
            x[0] = (c00 * r0 + (J[0][2] * J[2][1] - J[0][1] * J[2][2]) * r1
                             + (J[0][1] * J[1][2] - J[0][2] * J[1][1]) * r2) * inv;
            x[1] = (c01 * r0 + (J[0][0] * J[2][2] - J[0][2] * J[2][0]) * r1
                             + (J[0][2] * J[1][0] - J[0][0] * J[1][2]) * r2) * inv;
            x[2] = (c02 * r0 + (J[0][1] * J[2][0] - J[0][0] * J[2][1]) * r1
                             + (J[0][0] * J[1][1] - J[0][1] * J[1][0]) * r2) * inv;
        }
    }

    // The clamps-off path: while both clamp diodes sit past the reverse shortcut (phi = 0), the Newton on the M ON
    // ports (the pair's two junctions, and the feedback pair in a clipping block) with fixed sizes and indices; the
    // clamps' ports are linear in the ON ports' phi and follow once at the end. false = a clamp is (or comes) on at
    // the start or at the solution: the general Newton takes over from where this one stopped.
    template <int M>
    TOCAKE_SPAN_BODY bool solve_clamps_off(Block& b, const double* p0, double* ph, int& it, double& dlast) noexcept {
        using namespace lp_detail;
        constexpr int on[3] = {0, 1, 4};
        constexpr int PC = 2, PE = 3;
        double* p = b.p;
        auto clamps_off = [&] {
            return p[PC] * inv_vt_ < -40.0 && p[PE] * inv_vt_ < -40.0;
        };
        if (!clamps_off()) return false;
        double a[M][M];
        for (int i = 0; i < M; ++i) for (int j = 0; j < M; ++j) a[i][j] = b.A[on[i]][on[j]];
        double q[M], dq[M];
        bool conv = false;
        for (; it < MAXIT; ++it) {
            for (int j = 0; j < M; ++j) {
                const int k = on[j];
                port_phi(b.port_kind[k], p[k], q[j], dq[j]);
            }
            double J[M][M + 1];
            for (int i = 0; i < M; ++i) {
                double f = p[on[i]] - p0[on[i]];
                for (int j = 0; j < M; ++j) { f -= a[i][j] * q[j]; J[i][j] = (i == j ? 1.0 : 0.0) - a[i][j] * dq[j]; }
                J[i][M] = -f;
            }
            double x[M];
            cramer<M>(J, x);
            double dmax = 0.0, pmax = 0.0;
            for (int j = 0; j < M; ++j) {
                const int k = on[j];
                const double pn = limit(b.port_kind[k], p[k] + x[j], p[k]);
                const double d = pn - p[k];
                if (std::fabs(d) > dmax) dmax = std::fabs(d);
                q[j] += dq[j] * d;
                p[k] = pn;
                if (std::fabs(pn) > pmax) pmax = std::fabs(pn);
            }
            dlast = dmax;
            if (dmax < TOCAKE_LP_TOL * (1e-3 + pmax)) { conv = true; ++it; break; }
        }
        for (int k : {PC, PE}) {
            double s = p0[k];
            for (int j = 0; j < M; ++j) s += b.A[k][on[j]] * q[j];
            p[k] = s;
        }
        if (!conv || !clamps_off()) return false;
        for (int k = 0; k < b.np; ++k) ph[k] = 0.0;
        for (int j = 0; j < M; ++j) ph[on[j]] = q[j];
        return true;
    }

    // acc += x col over sixteen doubles, as eight vectors of two: one SSE2 register each in the generic build (vectors
    // of four measured 1.6 times slower there, lowered piecewise; 4 % faster in the x86-64-v3 clone)
    typedef double v2d __attribute__((vector_size(16)));
    TOCAKE_SPAN_BODY static void axpy8(v2d (&acc)[8], const double* col, double x) noexcept {
        const v2d xv = {x, x};
        for (int i = 0; i < 8; ++i) {
            v2d c;
            std::memcpy(&c, col + 2 * i, sizeof c);
            acc[i] += c * xv;
        }
    }

    TOCAKE_SPAN_BODY void solve_ports(Block& b, const double* u) noexcept {
        using namespace lp_detail;
        const int np = b.np, nd = D_CAP + b.nc;
        double p0[MAXP], ph[MAXP] = {}, dph[MAXP];
        v2d ew[8] = {};
        for (int d = 0; d < nd; ++d) axpy8(ew, b.EwT[d], u[d]);
        double ews[Block::NEW];
        std::memcpy(ews, ew, sizeof ews);
        for (int r = 0; r < MAXP; ++r) p0[r] = ews[r];
        double* p = b.p;
        ++dbg_solves_;
        int it = 0;
        bool done = false;
        double dlast = 0.0;
        if (b.np == 5) done = solve_clamps_off<3>(b, p0, ph, it, dlast);
        else           done = solve_clamps_off<2>(b, p0, ph, it, dlast);
        if (!done)
        for (; it < MAXIT; ++it) {
            int on[MAXP], off[MAXP], m = 0, mo = 0;
            for (int k = 0; k < np; ++k) {
                port_phi(b.port_kind[k], p[k], ph[k], dph[k]);
                if (dph[k] != 0.0) on[m++] = k; else off[mo++] = k;
            }
            double f[MAXP];
            for (int r = 0; r < np; ++r) {
                double s = p[r] - p0[r];
                for (int c = 0; c < m; ++c) s -= b.A[r][on[c]] * ph[on[c]];
                f[r] = s;
            }
            double J[MAXP][MAXP + 1];
            for (int i = 0; i < m; ++i) {
                const int r = on[i];
                for (int c = 0; c < m; ++c) J[i][c] = (i == c ? 1.0 : 0.0) - b.A[r][on[c]] * dph[on[c]];
                J[i][m] = -f[r];
            }
            for (int k = 0; k < m; ++k) {
                int q = k;
                for (int r = k + 1; r < m; ++r) if (std::fabs(J[r][k]) > std::fabs(J[q][k])) q = r;
                if (q != k) for (int c = k; c <= m; ++c) std::swap(J[k][c], J[q][c]);
                const double inv = 1.0 / J[k][k];
                for (int r = k + 1; r < m; ++r) {
                    const double fr = J[r][k] * inv;
                    for (int c = k + 1; c <= m; ++c) J[r][c] -= fr * J[k][c];
                }
                J[k][k] = inv;
            }
            double dp[MAXP], x[MAXP];
            for (int k = m - 1; k >= 0; --k) {
                double s = J[k][m];
                for (int c = k + 1; c < m; ++c) s -= J[k][c] * x[c];
                x[k] = s * J[k][k];
                dp[on[k]] = x[k];
            }
            for (int i = 0; i < mo; ++i) {
                const int r = off[i];
                double s = -f[r];
                for (int c = 0; c < m; ++c) s += b.A[r][on[c]] * dph[on[c]] * x[c];
                dp[r] = s;
            }
            double dmax = 0.0, pmax = 0.0;
            for (int k = 0; k < np; ++k) {
                const double pn = limit(b.port_kind[k], p[k] + dp[k], p[k]);
                const double d = pn - p[k];
                if (std::fabs(d) > dmax) dmax = std::fabs(d);
                ph[k] += dph[k] * d;   // phi at the new point to first order (the step is below the stop at the end)
                p[k] = pn;
                if (std::fabs(pn) > pmax) pmax = std::fabs(pn);
            }
            dlast = dmax;
            if (dmax < TOCAKE_LP_TOL * (1e-3 + pmax)) { done = true; ++it; break; }
        }
        dbg_iters_ += it;
        if (!done) ++dbg_capped_;
        if (dlast > dbg_worst_step_) dbg_worst_step_ = dlast;
        {
            double a, c; src_affine(BJT_B, p_.diode_Is, a, c);
            b.isrc[b.ib2_src] = a * ph[b.src_port[b.ib2_src]] + c;
        }
        for (int k = 0; k < np; ++k) axpy8(ew, b.AwT[k], ph[k]);
        std::memcpy(ews, ew, sizeof ews);
        for (int r = 0; r < MAXW; ++r) b.w[r] = ews[MAXP + r];
    }

    TOCAKE_SPAN_BODY static void commit(Block& b) noexcept {
        for (int k = 0; k < b.nc; ++k) {
            const double v = b.w[k];
            if (b.cap[k].rule == 0) {
                const double J = b.gc[k] * b.v1[k] + b.i1[k];
                b.i1[k] = b.gc[k] * v - J;
            }
            b.v2[k] = b.v1[k];
            b.v1[k] = v;
        }
    }

    // ---- the rest point (not real time) ----------------------------------------------------------------------------
    // Each block's DC solution on its whole nodal network (capacitors open), Newton with source stepping and the
    // junctions limited: in port space the open pair is badly conditioned. The
    // blocks meet as in the sample: block k drives k+1 and k+1's base current flows back into k, iterated to a fixed
    // point.
    bool dc_block(Block& b, double uin, double lag, double* x) noexcept {
        using namespace lp_detail;
        double Gm[MAXN][MAXN], B[MAXN][MAXD], N[MAXN][MAXS], Pv[MAXP][MAXN], Pu[MAXP][MAXD], Wv[MAXW][MAXN], Wu[MAXW][MAXD];
        b.build(true, Gm, B, N, Pv, Pu, Wv, Wu);
        const int n = b.n, np = b.np, ns = b.ns;
        double plim[MAXP];
        for (int k = 0; k < np; ++k) plim[k] = b.port_kind[k] == BJT_B ? 0.5 : 0.0;
        for (int step = 1; step <= 10; ++step) {
            const double lam = 0.1 * step;
            const double u[3] = { lam, lam * uin, lam * lag };
            int it = 0;
            for (; it < 500; ++it) {
                double A[MAXN][MAXN], rhs[MAXN];
                std::memcpy(A, Gm, sizeof A);
                for (int r = 0; r < n; ++r) rhs[r] = B[r][D_ONE] * u[0] + B[r][D_UIN] * u[1] + B[r][D_LAG] * u[2];
                double pu[MAXP];
                for (int k = 0; k < np; ++k) {
                    double s = Pu[k][D_ONE] * u[0] + Pu[k][D_UIN] * u[1];
                    pu[k] = s;
                    for (int c = 0; c < n; ++c) s += Pv[k][c] * x[c];
                    plim[k] = limit(b.port_kind[k], s, plim[k]);
                }
                for (int s = 0; s < ns; ++s) {
                    const int k = b.src_port[s];
                    double i, g;
                    source(b.src_kind[s], plim[k], i, g);
                    // i(p) ~ i + g (Pv x + pu - plim): G -= N g Pv, rhs += N (i + g (pu - plim))
                    for (int r = 0; r < n; ++r) {
                        if (N[r][s] == 0.0) continue;
                        for (int c = 0; c < n; ++c) A[r][c] -= N[r][s] * g * Pv[k][c];
                        rhs[r] += N[r][s] * (i + g * (pu[k] - plim[k]));
                    }
                }
                int piv[MAXN];
                if (!lu(A, n, piv)) return false;
                lu_solve(A, n, piv, rhs);
                double d = 0.0, m = 0.0;
                for (int r = 0; r < n; ++r) {
                    d = std::fmax(d, std::fabs(rhs[r] - x[r]));
                    m = std::fmax(m, std::fabs(rhs[r]));
                    x[r] = rhs[r];
                }
                if (!std::isfinite(d)) return false;
                if (d < 1e-13 + 1e-11 * m) break;
            }
            if (it == 500) return false;
        }
        // the rest state from the nodal solution
        const double u[3] = { 1.0, uin, lag };
        auto val = [&](const double* rv, const double* ru) {
            double s = ru[D_ONE] * u[0] + ru[D_UIN] * u[1];
            for (int c = 0; c < n; ++c) s += rv[c] * x[c];
            return s;
        };
        for (int k = 0; k < np; ++k) b.r_p[k] = val(Pv[k], Pu[k]);
        for (int k = 0; k < b.nw; ++k) b.r_w[k] = val(Wv[k], Wu[k]);
        for (int k = 0; k < b.nc; ++k) b.r_v1[k] = b.r_w[k];
        for (int s = 0; s < ns; ++s) { double g; source(b.src_kind[s], b.r_p[b.src_port[s]], b.r_isrc[s], g); }
        return true;
    }

    void solve_rest() noexcept {
        double x[4][lp_detail::MAXN] = {};
        double lag[3] = {};
        rest_ok_ = true;
        for (int round = 0; round < 6; ++round) {
            double drive = 0.0;                 // the jack rests at ground
            for (int k = 0; k < 4; ++k) {
                Block& b = blk_[k];
                if (!dc_block(b, drive, k < 3 ? lag[k] : 0.0, x[k])) rest_ok_ = false;
                drive = b.r_w[b.next_w];
            }
            for (int k = 0; k < 3; ++k) lag[k] = blk_[k + 1].r_isrc[blk_[k + 1].ib2_src];
        }
    }

    // Every state, not one node: a history poisoned while its node stayed finite would otherwise survive the guard.
    bool state_finite() const noexcept {
        // one sum per block and per array, so the additions do not form one long chain (the unused entries are 0)
        double s[4][3];
        for (int i = 0; i < 4; ++i) {
            const Block& b = blk_[i];
            double a = 0.0, c = 0.0, e = 0.0;
            for (int k = 0; k < lp_detail::MAXC; ++k) { a += b.v1[k]; c += b.i1[k]; }
            for (int k = 0; k < lp_detail::MAXP; ++k) e += b.p[k];
            s[i][0] = a; s[i][1] = c; s[i][2] = e;
        }
        double t = lag_[0] + lag_[1] + lag_[2];
        for (int i = 0; i < 4; ++i) t += (s[i][0] + s[i][1]) + s[i][2];
        return std::isfinite(t);
    }
    // A non-finite output or state: back to the rest point.
    float guard(double y) noexcept {
        if (!std::isfinite(y) || !state_finite()) {
            reset();
            ++dbg_nonfinite_;
            return 0.0f;
        }
        return static_cast<float>(y);
    }

    const LostPumpkinParams& p_;
    bool prepared_ = false, rest_ok_ = false;
    double fs_ = 48000.0;
    double sus_ = 1.0, tone_ = 0.5, vol_ = 0.5, vol_div_ = 0.5;
    Block blk_[4];
    int leg_r_[3] = {}, tone_r_ = 0;
    double lag_[3] = {};
    const double inv_vt_ = 1.0 / V_T_300K;
    const double nvt_d_ = p_.diode_N * V_T_300K, inv_nvt_d_ = 1.0 / nvt_d_;
    const double vcrit_d_ = nvt_d_ * std::log(nvt_d_ / (1.4142135623730951 * p_.diode_Is));
    const double vcrit_q_ = V_T_300K * std::log(V_T_300K / (1.4142135623730951 * njm4558::IS_Q));
    const double vcrit_dio_ = V_T_300K * std::log(V_T_300K / (1.4142135623730951 * njm4558::IS_D));
    long dbg_nonfinite_ = 0, dbg_solves_ = 0, dbg_iters_ = 0, dbg_capped_ = 0;
    double dbg_worst_step_ = 0.0;
};

}  // namespace tocake
