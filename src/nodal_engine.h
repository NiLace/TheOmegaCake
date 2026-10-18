// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// A NODAL engine for the same circuit descriptions src/dk_engine.h solves: the
// Newton runs on the voltages of the nodes the devices touch, over the Schur
// complement of the linear network, instead of on the port voltages over the
// dense K = N H N^T.
//
// Why a second form of the same Newton. Per iteration the DK
// factors a dense P x P matrix I + K D (P = 10 on the four-stage fuzz, 15 on the Jumbo
// Goat): 483 / 1462 multiply-adds, about half of run_sample on the shipped
// profile. The stages of the four-stage fuzz are RC-coupled, so after the linear-only
// nodes are eliminated the Schur complement over the q device nodes stays
// SPARSE (q = 14 for P = 10, 52 nonzeros, no fill under a static Markowitz
// order): a fixed-pattern factor and solve cost 102 / 123 multiply-adds. Which
// form is faster once the maps, the Jacobian build and the knob path are counted
// is measured per circuit (engine_select.h).
//
// Formulation. With the DK's own stamps, G x = b_c + b_u u + S s - N^T i(v),
// v = N x. Split the unknowns into D (device nodes, q of them) and L (the rest):
//
//     x_L = W_c + W_u u + W_s s - W_D x_D,     W_* = G_LL^-1 [b_c b_u S G_LD]_L
//     F(x_D) = Sc x_D - (Q_c + Q_u u + Q_s s) + N_D^T i(N_D x_D) = 0
//     Sc = G_DD - G_DL W_D,  Q_* = [.]_D - G_DL W_*
//
// so the per-sample Newton sees only q unknowns and a Jacobian Sc + N_D^T D N_D
// with a pattern fixed by the circuit. The capacitor voltages and the taps are
// affine in (1, u, s, x_D) through the same maps. The device law, the junction
// limiting (SPICE's pnjlim) and the DC rest point are the DK's, called through
// its friend access: no device equation is written twice.
//
// Junction limiting in a nodal Newton is SPICE's: the limited junction voltage
// is where the device is EVALUATED, and the node equations see the device
// linearised about it (i + D (v - v_lim)). The DK limits the port voltages
// themselves; the two engines take different iterates and agree at the root.
//
// Unlike the DK, a rebuild while a knob glides re-eliminates the whole linear
// block (|L|^3/3 + a few hundred multiply-adds): there is no low-rank update.

#pragma once

#include "dk_engine.h"
#include "dk_standard.h"
#include "dk_circuits.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <type_traits>

// Refinement passes of the linear block's elimination (the DK refines once).
#ifndef TOCAKE_NODAL_REFINE
#define TOCAKE_NODAL_REFINE 1
#endif

namespace tocake {

template <class C> struct nodal_lu { using type = void; };

template <class Circuit>
class NodalEngine {
public:
    using Dk = DkEngine<Circuit>;
    using Tap = typename Dk::Tap;
    using Gen = typename nodal_lu<Circuit>::type;
    // A circuit with a generated elimination (src/nodal_solver_generated.h, written from this engine's own
    // schedule) runs it: the same operations in the same order with literal indices, so the compiler unrolls them.
    static constexpr bool kGen = !std::is_void_v<Gen>;

    static constexpr int NN = Circuit::NN;
    static constexpr int NB = Circuit::NB;
    static constexpr int ND = Circuit::ND;
    static constexpr int NS = Dk::NS;
    static constexpr int kNDev = 4 * NB + ND + NS;   // devices()'s derivatives: 4 per BJT, 1 per diode, 1 per shunt
    static_assert(kNDev > 0, "a circuit with no device has nothing to solve");
    static constexpr int P  = Dk::P;
    static constexpr int PS = Dk::PS;
    static constexpr int S  = Circuit::S;
    static constexpr int NT = Circuit::NT;

    // ---- the partition: device nodes D, linear-only nodes L -----------------------
    // D holds every node a port reads; L is eliminated once per rebuild. The
    // Jumbo Goat nulls the DK at -113 dB at tone 0 and -183 at tone 0.5 (its gate
    // diode hangs on the tone wiper, whose pot leg is 1e-6 ohm at the end). Do not
    // promote the leg's far node into D: that puts a 1e6 S entry inside the
    // unpivoted factor.
    static constexpr std::array<bool, NN> make_is_port_node() noexcept {
        std::array<bool, NN> d{};
        for (int p = 0; p < P; ++p) {
            if (Dk::kPorts[p].plus  >= 0) d[Dk::kPorts[p].plus]  = true;
            if (Dk::kPorts[p].minus >= 0) d[Dk::kPorts[p].minus] = true;
        }
        return d;
    }
    static constexpr std::array<bool, NN> kIsPortNode = make_is_port_node();
    static constexpr int count_port_nodes() noexcept { int n = 0; for (bool b : kIsPortNode) n += b; return n; }
    static constexpr int Q = count_port_nodes();
    static constexpr int QMAX = Q;
    static constexpr int M = 2 + S + QMAX;          // columns of the elimination: bc, bu, S, G_LD
    struct PortDev { int plus, minus; };            // D indices, -1 = ground
    static constexpr std::array<int, NN> make_dev_index() noexcept {
        std::array<int, NN> m{}; int n = 0;
        for (int k = 0; k < NN; ++k) m[k] = kIsPortNode[k] ? n++ : -1;
        return m;
    }
    static constexpr std::array<int, NN> kDev = make_dev_index();
    static constexpr std::array<PortDev, P> make_port_dev() noexcept {
        std::array<PortDev, P> m{};
        for (int p = 0; p < P; ++p)
            m[p] = { Dk::kPorts[p].plus >= 0 ? kDev[Dk::kPorts[p].plus] : -1,
                     Dk::kPorts[p].minus >= 0 ? kDev[Dk::kPorts[p].minus] : -1 };
        return m;
    }
    static constexpr std::array<PortDev, P> kPortDev = make_port_dev();

    explicit NodalEngine(const CakeVariant& v, double p_tone = 0.5, double p_vol = 0.5, double p_sus = 1.0) noexcept
        : dk_(v, p_tone, p_vol, p_sus), v_(v) { k_.tone = p_tone; k_.vol = p_vol; k_.sus = p_sus; }

    void prepare(double sample_rate) { fs_ = sample_rate; dirty_ = true; reset(); }

    // The DC rest is the DK's (its nodal DC Newton on the resistive network);
    // then this engine's own matrices at the current knobs and rate.
    void reset() noexcept {
        dk_.set_sustain(k_.sus); dk_.set_tone(k_.tone); dk_.set_volume(k_.vol); dk_.set_bypass(k_.bypass);
        dk_.set_tolerance(tol_);
        dk_.prepare(fs_);
        dc_seed_settled_ = dk_.dc_seed_settled();
        rest_stale_ = !dc_seed_settled_;
        x_rest_ = dk_.rest_nodes();
        partition();
        coarse_ = false; have_order_ = false;
        rebuild();
        restore_rest();
        dbg_iters_ = dbg_calls_ = dbg_nonconv_ = dbg_nonfinite_ = dbg_rebuilds_ = 0; dbg_max_iters_ = 0; dbg_max_step_ = 0.0;
    }
    bool dc_seed_settled() const noexcept { return dc_seed_settled_; }
    bool branch_overflow() const noexcept { return dk_.branch_overflow(); }
    bool bad_part() const noexcept { return bad_part_ || dk_.bad_part(); }

    void set_sustain(double f) noexcept { set_pot(k_.sus, f); }
    void set_tone(double f)    noexcept { set_pot(k_.tone, f); }
    void set_volume(double f)  noexcept { set_pot(k_.vol, f); }
    void set_bypass(bool on)   noexcept { if (on != k_.bypass) { k_.bypass = on; dirty_ = true; rest_stale_ = true; } }
    // While a knob glides the product rebuilds every few samples: the rebuilds skip the elimination's refinement
    // pass (the DK's coarse glide), and the tick that lands forces one exact rebuild, so a settled knob runs the
    // exact matrices.
    void set_coarse(bool on) noexcept {
        if (coarse_ && !on) dirty_ = true;
        coarse_ = on;
    }
    void set_tolerance(double volts) noexcept { tol_ = volts; }
    double tolerance() const noexcept { return tol_; }

    float process_sample(float vin_volts) noexcept { return run_sample(static_cast<double>(vin_volts), Tap::Out); }
    // The product's entry shape (src/isa_hot.h, src/isa_dispatch.h): the sample body out of line, in the two
    // instruction sets the DK ships, and the span loops the dispatch drives.
    TOCAKE_SPAN_BASE float span_sample(float vin_volts) noexcept { return run_sample(static_cast<double>(vin_volts), Tap::Out); }
    TOCAKE_SPAN_ISA float span_sample_isa_v3(float vin_volts) noexcept { return run_sample(static_cast<double>(vin_volts), Tap::Out); }
    void span_run_isa_v3(float* buf, uint32_t n) noexcept { for (uint32_t j = 0; j < n; ++j) buf[j] = span_sample_isa_v3(buf[j]); }
    void span_run(float* buf, uint32_t n) noexcept { for (uint32_t j = 0; j < n; ++j) buf[j] = span_sample(buf[j]); }
    float process_sample(float vin_volts, Tap tap) noexcept { return run_sample(static_cast<double>(vin_volts), tap); }

    long   dbg_iters()     const noexcept { return dbg_iters_; }
    long   dbg_calls()     const noexcept { return dbg_calls_; }
    long   dbg_nonconv()   const noexcept { return dbg_nonconv_; }
    int    dbg_max_iters() const noexcept { return dbg_max_iters_; }
    double dbg_max_step()  const noexcept { return dbg_max_step_; }
    double dbg_max_resid() const noexcept { return dbg_max_step_; }   // the DK's name for the same figure (tools written against the node engines)
    long   nonfinite_resets() const noexcept { return dbg_nonfinite_; }
    long   rebuilds()      const noexcept { return dbg_rebuilds_; }
    // The Schur block's size and pattern: nonzeros before and after the symbolic fill.
    int    block_size()    const noexcept { return q_; }
    int    pattern_nnz()   const noexcept { return nnz0_; }
    int    pattern_fill()  const noexcept { return fill_; }
    // Rebuilds whose live pattern was not the generated one (must stay 0).
    long   pattern_mismatches() const noexcept { return dbg_pat_mismatch_; }
    // The schedule, for scripts/gen_nodal_solver.py.
    void dump_schedule(std::FILE* f, const char* name) const {
        std::fprintf(f, "NODAL %s Q=%d NN=%d P=%d S=%d\n", name, q_, NN, P, S);
        std::fprintf(f, "DEV"); for (int d = 0; d < q_; ++d) std::fprintf(f, " %d", dev_node_[d]); std::fprintf(f, "\n");
        std::fprintf(f, "PAT"); for (int i = 0; i < q_ * q_; ++i) if (pat0_[i]) std::fprintf(f, " %d", i); std::fprintf(f, "\n");
        std::fprintf(f, "SC"); for (int k = 0; k < nsc_; ++k) std::fprintf(f, " %d", sc_r_[k] * q_ + sc_c_[k]); std::fprintf(f, "\n");
        std::fprintf(f, "FILL"); for (int i = 0; i < q_ * q_; ++i) if (patf_[i] && !pat0_[i]) std::fprintf(f, " %d", i); std::fprintf(f, "\n");
        std::fprintf(f, "ORDER"); for (int s = 0; s < q_; ++s) std::fprintf(f, " %d", order_[s]); std::fprintf(f, "\n");
        for (int s = 0; s < q_; ++s) {
            std::fprintf(f, "STEP %d rows", s); for (int a = 0; a < nrows_[s]; ++a) std::fprintf(f, " %d", rows_[s][a]);
            std::fprintf(f, " cols"); for (int a = 0; a < ncols_[s]; ++a) std::fprintf(f, " %d", cols_[s][a]);
            std::fprintf(f, " urows"); for (int a = 0; a < nurows_[s]; ++a) std::fprintf(f, " %d", urows_[s][a]);
            std::fprintf(f, "\n");
        }
        std::fprintf(f, "END\n");
    }
    // Port voltages after the last sample (the device's evaluation point).
    const std::array<double, P>& solved_ports() const noexcept { return veval_; }
    // The DC rest point adopted by the last reset(), node by node (the DK's, copied there).
    const std::array<double, NN>& rest_nodes() const noexcept { return x_rest_; }

private:
    void set_pot(double& p, double f) noexcept { f = std::clamp(f, 0.0, 1.0); if (f != p) { p = f; dirty_ = true; rest_stale_ = true; } }

    double supply() const noexcept { return dk_.supply(); }

    // ---- the linear network, stamped as the DK stamps it -------------------------
    struct Branch { int a, b; double g; };
    static constexpr int kMaxBranches = TOCAKE_DK_MAX_BRANCHES;
    struct Net { std::array<double, NN * NN> G{}; std::array<double, NN> bc{}, bu{}; std::array<Branch, kMaxBranches> br{}; int nbr = 0; };
    static void stamp(Net& n, int a, int b, double g, double vcc) noexcept {
        if (n.nbr < kMaxBranches) n.br[n.nbr++] = { a, b, g };
        for (int side = 0; side < 2; ++side) {
            const int i = side ? b : a, j = side ? a : b;
            if (i < 0) continue;
            n.G[i * NN + i] += g;
            if (j >= 0) n.G[i * NN + j] -= g;
            else if (j == kDkVcc) n.bc[i] += g * vcc;
            else if (j == kDkIn) n.bu[i] += g;
        }
    }

    // Dense LU with partial pivoting on the nl x nl linear block (row-major, leading dimension nl), m right-hand sides.
    // In double: long double does not move the Jumbo Goat's tone-end null and doubles the cost of a rebuild.
    using ld = double;
    static bool lu_l(std::array<ld, NN * NN>& A, std::array<int, NN>& perm, int nl) noexcept {
        for (int i = 0; i < nl; ++i) perm[i] = i;
        for (int k = 0; k < nl; ++k) {
            int piv = k; ld best = std::fabs(A[k * nl + k]);
            for (int r = k + 1; r < nl; ++r)
                if (std::fabs(A[r * nl + k]) > best) { best = std::fabs(A[r * nl + k]); piv = r; }
            if (!(best > 0.0)) return false;
            if (piv != k) {
                for (int c = 0; c < nl; ++c) std::swap(A[k * nl + c], A[piv * nl + c]);
                std::swap(perm[k], perm[piv]);
            }
            const ld inv = 1.0 / A[k * nl + k];
            for (int r = k + 1; r < nl; ++r) {
                const ld f = A[r * nl + k] * inv;
                A[r * nl + k] = f;
                if (f != 0.0) for (int c = k + 1; c < nl; ++c) A[r * nl + c] -= f * A[k * nl + c];
            }
        }
        return true;
    }
    static void solve_l(const std::array<ld, NN * NN>& A, const std::array<int, NN>& perm, int nl, int m,
                        const ld* B, ld* X, ld* Y) noexcept {   // B, X, Y: nl x m row-major
        for (int r = 0; r < nl; ++r) {
            ld* yr = &Y[r * m]; const ld* br = B + perm[r] * m;
            for (int j = 0; j < m; ++j) yr[j] = br[j];
            for (int c = 0; c < r; ++c) { const ld a = A[r * nl + c]; const ld* yc = &Y[c * m]; for (int j = 0; j < m; ++j) yr[j] -= a * yc[j]; }
        }
        for (int r = nl - 1; r >= 0; --r) {
            ld* xr = X + r * m; const ld* yr = &Y[r * m];
            for (int j = 0; j < m; ++j) xr[j] = yr[j];
            for (int c = r + 1; c < nl; ++c) { const ld a = A[r * nl + c]; const ld* xc = X + c * m; for (int j = 0; j < m; ++j) xr[j] -= a * xc[j]; }
            const ld d = A[r * nl + r];
            for (int j = 0; j < m; ++j) xr[j] /= d;
        }
    }

    // The partition, once per reset(): the port nodes, plus every unknown node a
    // pot leg joins to one of them (a branch whose conductance differs between
    // two knob settings).
    void partition() noexcept {
        q_ = 0;
        for (int k = 0; k < NN; ++k) { dev_[k] = kIsPortNode[k] ? q_ : -1; if (kIsPortNode[k]) dev_node_[q_++] = k; }
        nl_ = 0;
        for (int k = 0; k < NN; ++k) lin_[k] = kIsPortNode[k] ? -1 : nl_++;
        for (int p = 0; p < P; ++p)
            port_dev_[p] = { Dk::kPorts[p].plus >= 0 ? dev_[Dk::kPorts[p].plus] : -1,
                             Dk::kPorts[p].minus >= 0 ? dev_[Dk::kPorts[p].minus] : -1 };
    }

    // Build the Schur complement and every map, at the current rate and knobs.
    void rebuild() noexcept {
        const int q = q_, nl = nl_, m = 2 + S + q;
        Net n;
        const double vcc = supply();
        Circuit::resistors(v_, k_, [&](int a, int b, double r) {   // as the DK's resistors(): a bad field is flagged
            if (!(r > 0.0) || !std::isfinite(r)) { bad_part_ = true; return; }
            stamp(n, a, b, 1.0 / r, vcc);
        });
        const auto C = Circuit::cap_values(v_);
        for (int k = 0; k < S; ++k) { gc_[k] = 2.0 * C[k] * fs_; stamp(n, Circuit::kCaps[k].a, Circuit::kCaps[k].b, gc_[k], vcc); }
        for (int k = 0; k < S; ++k) g2_[k] = 2.0 * gc_[k];

        // Right-hand sides on L: [bc bu S-columns G_LD-columns], nl x m; G_LL (nl x nl).
        auto& GLL = wk_gll_; auto& B = wk_b_; auto& W = wk_w_;
        GLL.fill(0.0); B.fill(0.0); W.fill(0.0);
        for (int i = 0; i < NN; ++i) {
            const int li = lin_[i];
            if (li < 0) continue;
            for (int j = 0; j < NN; ++j) {
                const double g = n.G[i * NN + j];
                if (g == 0.0) continue;
                if (lin_[j] >= 0) GLL[li * nl + lin_[j]] = g;
                else              B[li * m + 2 + S + dev_[j]] = g;        // G_LD column
            }
            B[li * m + 0] = n.bc[i]; B[li * m + 1] = n.bu[i];
        }
        for (int j = 0; j < S; ++j) {
            const int a = Circuit::kCaps[j].a, b = Circuit::kCaps[j].b;
            if (a >= 0 && lin_[a] >= 0) B[lin_[a] * m + 2 + j] += 1.0;
            if (b >= 0 && lin_[b] >= 0) B[lin_[b] * m + 2 + j] -= 1.0;
        }
        if (nl > 0) {
            auto& LU = wk_lu_; LU = GLL;
            if (!lu_l(LU, perm_, nl)) { dirty_ = false; return; }       // unreachable for a connected network
            solve_l(LU, perm_, nl, m, B.data(), W.data(), wk_y_.data());
            // Refinement, the residual summed BRANCH BY BRANCH as g (xa - xb), which a
            // 1e6 S pot leg does not spoil (the DK's reason: G x cancels two ~1e6
            // products into nanoamps). A branch to a D node or a known voltage
            // contributes g xa alone: that side is not an unknown of this block.
            for (int pass = 0; pass < (coarse_ ? 0 : TOCAKE_NODAL_REFINE); ++pass) {
                auto& R = wk_r_; auto& D = wk_d_;
                R = B;
                for (int k = 0; k < n.nbr; ++k) {
                    const int a = n.br[k].a, b = n.br[k].b;
                    const int la = (a >= 0) ? lin_[a] : -1, lb = (b >= 0) ? lin_[b] : -1;
                    if (la < 0 && lb < 0) continue;
                    const ld g = n.br[k].g;
                    for (int j = 0; j < m; ++j) {
                        const ld cur = g * ((la >= 0 ? W[la * m + j] : 0.0) - (lb >= 0 ? W[lb * m + j] : 0.0));
                        if (la >= 0) R[la * m + j] -= cur;
                        if (lb >= 0) R[lb * m + j] += cur;
                    }
                }
                solve_l(LU, perm_, nl, m, R.data(), D.data(), wk_y_.data());
                for (int i = 0; i < nl * m; ++i) W[i] += D[i];
            }
        }
        // Sc = G_DD - G_DL W_D;  Q_* = [.]_D - G_DL W_*   (T = [Qc Qu Qs Sc] over D rows)
        auto& T = wk_t_; T.fill(0.0);
        for (int i = 0; i < NN; ++i) {
            const int di = dev_[i];
            if (di < 0) continue;
            T[di * m + 0] = n.bc[i]; T[di * m + 1] = n.bu[i];
            for (int j = 0; j < NN; ++j) {
                const double g = n.G[i * NN + j];
                if (g == 0.0) continue;
                if (dev_[j] >= 0) T[di * m + 2 + S + dev_[j]] += g;
                else { const ld* w = &W[lin_[j] * m]; for (int c = 0; c < m; ++c) T[di * m + c] -= static_cast<ld>(g) * w[c]; }
            }
        }
        for (int j = 0; j < S; ++j) {
            const int a = Circuit::kCaps[j].a, b = Circuit::kCaps[j].b;
            if (a >= 0 && dev_[a] >= 0) T[dev_[a] * m + 2 + j] += 1.0;
            if (b >= 0 && dev_[b] >= 0) T[dev_[b] * m + 2 + j] -= 1.0;
        }
        for (int d = 0; d < q; ++d) {                      // back to double for the sample
            Qc_[d] = static_cast<double>(T[d * m + 0]); Qu_[d] = static_cast<double>(T[d * m + 1]);
            for (int j = 0; j < S; ++j) Qs_[d * S + j] = static_cast<double>(T[d * m + 2 + j]);
            for (int e = 0; e < q; ++e) Sc_[d * q + e] = static_cast<double>(T[d * m + 2 + S + e]);
        }
        // A pin's voltage as an affine map in (1, u, s, x_D): the capacitor voltages and the taps.
        auto pin_map = [&](int node, double sign, double* c, double* u, double* sv, double* dv) {
            if (node >= 0 && dev_[node] >= 0) { dv[dev_[node]] += sign; return; }
            if (node >= 0) {                                    // an L node: its W row
                const ld* w = &W[lin_[node] * m];
                *c += static_cast<double>(sign * w[0]); *u += static_cast<double>(sign * w[1]);
                for (int j = 0; j < S; ++j) sv[j] += static_cast<double>(sign * w[2 + j]);
                for (int e = 0; e < q; ++e) dv[e] -= static_cast<double>(sign * w[2 + S + e]);
                return;
            }
            if (node == kDkIn)  *u += sign;
            if (node == kDkVcc) *c += sign * vcc;
        };
        Xc_.fill(0.0); Xu_.fill(0.0); Xs_.fill(0.0); Xd_.fill(0.0);
        for (int k = 0; k < S; ++k) {
            pin_map(Circuit::kCaps[k].a, +1.0, &Xc_[k], &Xu_[k], &Xs_[k * S], &Xd_[k * q]);
            pin_map(Circuit::kCaps[k].b, -1.0, &Xc_[k], &Xu_[k], &Xs_[k * S], &Xd_[k * q]);
        }
        Oc_.fill(0.0); Ou_.fill(0.0); Os_.fill(0.0); Od_.fill(0.0);
        for (int t = 0; t < NT; ++t)
            pin_map(Circuit::kTaps[t].node, +1.0, &Oc_[t], &Ou_[t], &Os_[t * S], &Od_[t * q]);
        // The elimination order depends on the PATTERN only, which the knobs do not move: recompute it when the
        // pattern changed (first rebuild, a bypass switch), not on every glide tick (Markowitz is O(q^3) per step).
        bool same = have_order_;
        if (same) for (int i = 0; i < q * q; ++i) if ((Sc_[i] != 0.0) != sc_pat_[i]) { same = false; break; }
        if (!same) { symbolic(); have_order_ = true; for (int i = 0; i < q * q; ++i) sc_pat_[i] = Sc_[i] != 0.0; }
        else if constexpr (!kGen) { nsc_ = 0; for (int i = 0; i < q * q; ++i) if (Sc_[i] != 0.0) { sc_r_[nsc_] = i / q; sc_c_[nsc_] = i % q; ++nsc_; } }
        dirty_ = false;
        ++dbg_rebuilds_;
    }

    // ---- the fixed-pattern elimination ----------------------------------------------
    // The pattern is Sc's nonzeros plus the device blocks; a static Markowitz
    // order is chosen on it and the fill it creates is recorded. Factor and
    // solve then touch only those cells. No pivoting: the linear part is a passive admittance, but the devices add
    // unsymmetric terms, so the fixed order rests on measurement: it holds on every shipped circuit's equations,
    // and a zero pivot would end in the non-finite guard.
    void symbolic() noexcept {
        const int q = q_;
        std::array<bool, QMAX * QMAX> pat{};
        nsc_ = 0;
        for (int i = 0; i < q * q; ++i) {
            pat[i] = Sc_[i] != 0.0;
            if (pat[i]) { sc_r_[nsc_] = i / q; sc_c_[nsc_] = i % q; ++nsc_; }
        }
        for (int p = 0; p < P; ++p) {                      // the device stamps over the port pins
            const int ids[2] = { port_dev_[p].plus, port_dev_[p].minus };
            for (int a : ids) for (int b : ids) if (a >= 0 && b >= 0) pat[a * q + b] = true;
        }
        for (int k = 0; k < NB; ++k) {                     // the 2x2 transistor block couples vbe and vbc
            const int ids[4] = { port_dev_[2 * k].plus, port_dev_[2 * k].minus, port_dev_[2 * k + 1].plus, port_dev_[2 * k + 1].minus };
            for (int a : ids) for (int b : ids) if (a >= 0 && b >= 0) pat[a * q + b] = true;
        }
        for (int d = 0; d < q; ++d) pat[d * q + d] = true;
        nnz0_ = 0; for (int i = 0; i < q * q; ++i) nnz0_ += pat[i];
        pat0_ = pat;
        if constexpr (kGen) {                              // the generated schedule's pattern must be the live one
            bool same = (Gen::Q == q);
            if (same) for (int i = 0; i < q * q; ++i) same &= (pat[i] == Gen::kPat[i]);
            if (!same) ++dbg_pat_mismatch_;
        }
        std::array<bool, QMAX> alive{}; for (int i = 0; i < q; ++i) alive[i] = true;
        fill_ = 0;
        for (int step = 0; step < q; ++step) {
            int best = -1; long cost = 0;
            for (int k = 0; k < q; ++k) {
                if (!alive[k]) continue;
                long r = 0, c = 0;
                for (int j = 0; j < q; ++j) if (alive[j] && j != k) { r += pat[k * q + j]; c += pat[j * q + k]; }
                if (best < 0 || r * c < cost) { best = k; cost = r * c; }
            }
            const int k = best;
            order_[step] = k;
            int nr = 0, nc = 0;
            for (int i = 0; i < q; ++i) if (alive[i] && i != k && pat[i * q + k]) rows_[step][nr++] = i;
            for (int j = 0; j < q; ++j) if (alive[j] && j != k && pat[k * q + j]) cols_[step][nc++] = j;
            nrows_[step] = nr; ncols_[step] = nc;
            for (int a = 0; a < nr; ++a)
                for (int b = 0; b < nc; ++b)
                    if (!pat[rows_[step][a] * q + cols_[step][b]]) { pat[rows_[step][a] * q + cols_[step][b]] = true; ++fill_; }
            alive[k] = false;
        }
        patf_ = pat;
        std::array<int, QMAX> pos{};
        for (int s = 0; s < q; ++s) pos[order_[s]] = s;
        for (int s = 0; s < q; ++s) {
            const int k = order_[s]; int nu = 0;
            for (int i = 0; i < q; ++i) if (pos[i] < s && pat[i * q + k]) urows_[s][nu++] = i;
            nurows_[s] = nu;
        }
    }
    // J in place (leading dimension q_): L below the pivots as multipliers, U above; ip = 1/pivot.
    void factor(double* J, double* ip) const noexcept {
        constexpr int q = Q;
        for (int s = 0; s < q; ++s) {
            const int k = order_[s];
            const double inv = 1.0 / J[k * q + k];
            ip[s] = inv;
            for (int a = 0; a < nrows_[s]; ++a) {
                const int r = rows_[s][a];
                const double f = J[r * q + k] * inv;
                J[r * q + k] = f;
                for (int b = 0; b < ncols_[s]; ++b) { const int c = cols_[s][b]; J[r * q + c] -= f * J[k * q + c]; }
            }
        }
    }
    void solve(const double* J, const double* ip, double* y) const noexcept {
        constexpr int q = Q;
        for (int s = 0; s < q; ++s) {
            const int k = order_[s]; const double yk = y[k];
            for (int a = 0; a < nrows_[s]; ++a) { const int r = rows_[s][a]; y[r] -= J[r * q + k] * yk; }
        }
        for (int s = q - 1; s >= 0; --s) {
            const int k = order_[s]; const double yk = y[k] * ip[s]; y[k] = yk;
            for (int a = 0; a < nurows_[s]; ++a) { const int r = urows_[s][a]; y[r] -= J[r * q + k] * yk; }
        }
    }

    // ---- rest and state -------------------------------------------------------------
    double rest_pin(int n) const noexcept { return n >= 0 ? x_rest_[n] : (n == kDkVcc ? supply() : 0.0); }
    // The non-finite guard's path: the DC rest AT THE CURRENT KNOBS (the DK twin re-solves it after a knob move; its
    // reset() result describes the old operating point). Cold and out of line.
    [[gnu::noinline, gnu::cold]] void recover_nonfinite() noexcept {
        if (rest_stale_) {
            dk_.set_sustain(k_.sus); dk_.set_tone(k_.tone); dk_.set_volume(k_.vol); dk_.set_bypass(k_.bypass);
            if (dk_.solve_dc()) { x_rest_ = dk_.x_rest_; rest_stale_ = false; }
        }
        restore_rest();
        ++dbg_nonfinite_;
    }
    void restore_rest() noexcept {
        for (int k = 0; k < S; ++k) {
            const int a = Circuit::kCaps[k].a, b = Circuit::kCaps[k].b;
            s_[k] = gc_[k] * (rest_pin(a) - rest_pin(b));
        }
        for (int d = 0; d < q_; ++d) { xd_[d] = x_rest_[dev_node_[d]]; xd1_[d] = xd2_[d] = xd_[d]; }
        ports_of(xd_.data(), veval_.data());
        for (int t = 0; t < NT; ++t)
            tap_rest_[t] = Circuit::kTaps[t].ac ? rest_pin(Circuit::kTaps[t].node) : 0.0;
    }
    static void ports_of(const double* xd, double* v) noexcept {
        for (int p = 0; p < P; ++p)
            v[p] = (kPortDev[p].plus >= 0 ? xd[kPortDev[p].plus] : 0.0) - (kPortDev[p].minus >= 0 ? xd[kPortDev[p].minus] : 0.0);
    }
    // N^T applied to a port vector: current into plus, out of minus.
    static void add_ports(double* f, const double* ip_) noexcept {
        for (int p = 0; p < P; ++p) {
            if (kPortDev[p].plus  >= 0) f[kPortDev[p].plus]  += ip_[p];
            if (kPortDev[p].minus >= 0) f[kPortDev[p].minus] -= ip_[p];
        }
    }
    // + g * (e_plus_p - e_minus_p)(e_plus_q - e_minus_q)^T into J (leading dimension q_).
    static void add_jac(double* J, int p, int qq, double g) noexcept {
        const int rp[2] = { kPortDev[p].plus, kPortDev[p].minus }, cq[2] = { kPortDev[qq].plus, kPortDev[qq].minus };
        for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b)
            if (rp[a] >= 0 && cq[b] >= 0) J[rp[a] * Q + cq[b]] += ((a ^ b) ? -g : g);
    }

    // ---- the sample --------------------------------------------------------------------
    TOCAKE_SPAN_BODY float run_sample(double u, Tap tap) noexcept {
        if (dirty_) rebuild();
        constexpr int q = Q;                               // q_ == Q by construction (partition())
        double bt[QMAX];
        for (int d = 0; d < q; ++d) bt[d] = Qc_[d] + Qu_[d] * u;
        for (int k = 0; k < S; ++k) { const double sk = s_[k]; for (int d = 0; d < q; ++d) bt[d] += Qs_[d * S + k] * sk; }
        for (int d = 0; d < q; ++d) {                      // the DK's extrapolation, on the nodes
            const double cur = xd_[d];
            xd_[d] = (Dk::kPredictor == 2) ? 3.0 * cur - 3.0 * xd1_[d] + xd2_[d]
                   : (Dk::kPredictor == 1) ? 2.0 * cur - xd1_[d] : cur;
            xd2_[d] = xd1_[d]; xd1_[d] = cur;
        }
        double vp[P], vlim[P], i[P], d[kNDev], J[QMAX * QMAX], ip[QMAX], F[QMAX];
        ports_of(xd_.data(), vp);
        dk_.limit_all(vp, veval_.data());                  // the extrapolated junctions, limited like a step
        for (int p = 0; p < P; ++p) vlim[p] = vp[p];
        int it = 0; bool conv = false; double step = 0.0;
        for (; it < TOCAKE_DK_NEWTON_CAP; ++it) {
            dk_.devices(vlim, i, d);
            // Residual: Sc x - b + N^T (i + D (v - v_lim)), v = N x the node-implied ports.
            ports_of(xd_.data(), vp);
            double ilin[P];
            for (int k = 0; k < NB; ++k) {
                const int be = 2 * k, bc = be + 1;
                const double dbe = vp[be] - vlim[be], dbc = vp[bc] - vlim[bc];
                ilin[be] = i[be] + d[4 * k + 0] * dbe + d[4 * k + 1] * dbc;
                ilin[bc] = i[bc] + d[4 * k + 2] * dbe + d[4 * k + 3] * dbc;
            }
            for (int k = 0; k < ND; ++k) ilin[2 * NB + k] = i[2 * NB + k] + d[4 * NB + k] * (vp[2 * NB + k] - vlim[2 * NB + k]);
            for (int k = 0; k < NS; ++k) ilin[PS + k] = i[PS + k] + d[4 * NB + ND + k] * (vp[PS + k] - vlim[PS + k]);
            for (int r = 0; r < q; ++r) F[r] = -bt[r];
            if constexpr (kGen) Gen::residual(F, Sc_.data(), xd_.data());
            else for (int k = 0; k < nsc_; ++k) F[sc_r_[k]] += Sc_[sc_r_[k] * q + sc_c_[k]] * xd_[sc_c_[k]];
            add_ports(F, ilin);
            if constexpr (kGen) Gen::assemble(J, Sc_.data());
            else for (int k = 0; k < q * q; ++k) J[k] = Sc_[k];
            for (int k = 0; k < NB; ++k) {
                const int be = 2 * k, bc = be + 1;
                add_jac(J, be, be, d[4 * k + 0]); add_jac(J, be, bc, d[4 * k + 1]);
                add_jac(J, bc, be, d[4 * k + 2]); add_jac(J, bc, bc, d[4 * k + 3]);
            }
            for (int k = 0; k < ND; ++k) add_jac(J, 2 * NB + k, 2 * NB + k, d[4 * NB + k]);
            for (int k = 0; k < NS; ++k) add_jac(J, PS + k, PS + k, d[4 * NB + ND + k]);
            if constexpr (kGen) Gen::factor(J, ip); else factor(J, ip);
            for (int r = 0; r < q; ++r) F[r] = -F[r];
            if constexpr (kGen) Gen::solve(J, ip, F); else solve(J, ip, F);
            for (int r = 0; r < q; ++r) xd_[r] += F[r];
            // The step is read on the junctions, limited against the last evaluation point.
            double vnew[P];
            ports_of(xd_.data(), vnew);
            dk_.limit_all(vnew, vlim);
            step = 0.0;
            [[maybe_unused]] double step_tight = 0.0;
            for (int p = 0; p < P; ++p) {
                const double sd = vnew[p] - vlim[p];
                vlim[p] = vnew[p];
                if constexpr (dk_circuit_tight<Circuit>::has) {
                    if (Dk::kTightMask[p]) { step_tight = std::max(step_tight, std::fabs(sd)); continue; }
                }
                step = std::max(step, std::fabs(sd));
            }
            if constexpr (dk_circuit_tight<Circuit>::has) {
                const double t = Circuit::kTightTol < tol_ ? Circuit::kTightTol : tol_;
                if (!(step_tight < t)) step = step + tol_ + step_tight;
            }
            if (step < tol_) { ++it; conv = true; break; }
        }
        ++dbg_calls_; dbg_iters_ += it;
        if (it > dbg_max_iters_) dbg_max_iters_ = it;
        if (!conv) ++dbg_nonconv_;
        if (!(step <= dbg_max_step_)) dbg_max_step_ = step;
        for (int p = 0; p < P; ++p) veval_[p] = vlim[p];

        const int t = (tap == Tap::Out) ? 0 : (tap == Tap::C1Out) ? 1 : (tap == Tap::C2Out) ? 2 : (NT > 3 ? 3 : 0);
        double out = Oc_[t] + Ou_[t] * u;
        for (int k = 0; k < S; ++k) out += Os_[t * S + k] * s_[k];
        for (int e = 0; e < q; ++e) out += Od_[t * q + e] * xd_[e];
        out -= tap_rest_[t];

        double vc[S];
        for (int r = 0; r < S; ++r) vc[r] = Xc_[r] + Xu_[r] * u;
        for (int k = 0; k < S; ++k) { const double sk = s_[k]; for (int r = 0; r < S; ++r) vc[r] += Xs_[r * S + k] * sk; }
        for (int e = 0; e < q; ++e) { const double xe = xd_[e]; for (int r = 0; r < S; ++r) vc[r] += Xd_[r * q + e] * xe; }
        for (int k = 0; k < S; ++k) s_[k] = g2_[k] * vc[k] - s_[k];

        if (!all_finite(xd_) || !all_finite(s_) || !std::isfinite(out)) {
            recover_nonfinite();
            out = 0.0;
        }
        return static_cast<float>(out);
    }

    Dk dk_;
    const CakeVariant& v_;
    DkKnobs k_;
    double fs_ = 96000.0;
    double tol_ = dk_circuit_tol<Circuit>::value;
    bool dirty_ = true, dc_seed_settled_ = false, coarse_ = false, have_order_ = false;
    bool bad_part_ = false, rest_stale_ = false;
    std::array<bool, QMAX * QMAX> sc_pat_{};

    // the partition (reset) and the working buffers of a rebuild (members, not stack: a rebuild runs on the audio thread)
    int q_ = 0, nl_ = NN;
    std::array<int, NN> dev_{}, lin_{}, dev_node_{};
    std::array<PortDev, P> port_dev_{};
    std::array<int, NN> perm_{};
    std::array<ld, NN * NN> wk_gll_{}, wk_lu_{};
    std::array<ld, NN * M> wk_b_{}, wk_w_{}, wk_r_{}, wk_d_{}, wk_y_{};
    std::array<ld, QMAX * M> wk_t_{};

    std::array<double, S> gc_{}, g2_{};
    std::array<double, QMAX * QMAX> Sc_{};
    std::array<double, QMAX> Qc_{}, Qu_{};
    std::array<double, QMAX * S> Qs_{};
    std::array<double, S> Xc_{}, Xu_{};
    std::array<double, S * S> Xs_{};
    std::array<double, S * QMAX> Xd_{};
    std::array<double, NT> Oc_{}, Ou_{};
    std::array<double, NT * S> Os_{};
    std::array<double, NT * QMAX> Od_{};

    std::array<int, QMAX> order_{}, nrows_{}, ncols_{}, nurows_{};
    std::array<std::array<int, QMAX>, QMAX> rows_{}, cols_{}, urows_{};
    std::array<int, QMAX * QMAX> sc_r_{}, sc_c_{}; int nsc_ = 0;   // Sc's nonzeros, listed
    int nnz0_ = 0, fill_ = 0;
    std::array<bool, QMAX * QMAX> pat0_{}, patf_{};
    long dbg_pat_mismatch_ = 0;

    std::array<double, NN> x_rest_{};
    std::array<double, NT> tap_rest_{};
    std::array<double, S> s_{};
    std::array<double, QMAX> xd_{}, xd1_{}, xd2_{};
    std::array<double, P> veval_{};

    long dbg_iters_ = 0, dbg_calls_ = 0, dbg_nonconv_ = 0, dbg_nonfinite_ = 0, dbg_rebuilds_ = 0;
    int dbg_max_iters_ = 0;
    double dbg_max_step_ = 0.0;
};

}  // namespace tocake

#if __has_include("nodal_solver_generated.h")
#include "nodal_solver_generated.h"
#endif
