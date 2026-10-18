// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// Small-system Newton-Raphson solver for per-sample non-linear circuit
// equations. Generic over dimension; specialised 3x3 and 4x4 linear
// solves for speed.
//
// Usage: provide a residual function F(x) returning array<N> and a
// Jacobian function J(x) returning array<array<N>, N>. The solver
// iterates x <- x - J^{-1} F until ||F|| < tol or max_iter is reached.
// Caller is expected to warm-start with the previous sample's solution
// (typical convergence: 2-4 iterations).

#pragma once

#include <array>
#include <cmath>

namespace tocake {

template <std::size_t N>
struct nr_result {
    std::array<double, N> x;
    int    iterations;
    bool   converged;
    double final_residual_inf;   // max(|F_i|) at the final iterate
};

namespace detail {

// 3x3 linear solve via Cramer's rule. Returns false if J is singular.
inline bool solve_3x3(const double J[3][3], const double rhs[3], double out[3]) noexcept
{
    auto det3 = [](double a,double b,double c,
                   double d,double e,double f,
                   double g,double h,double i) noexcept {
        return a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
    };
    const double D = det3(J[0][0], J[0][1], J[0][2],
                          J[1][0], J[1][1], J[1][2],
                          J[2][0], J[2][1], J[2][2]);
    if (std::abs(D) < 1e-30) return false;
    const double invD = 1.0 / D;
    out[0] = invD * det3(rhs[0], J[0][1], J[0][2],
                         rhs[1], J[1][1], J[1][2],
                         rhs[2], J[2][1], J[2][2]);
    out[1] = invD * det3(J[0][0], rhs[0], J[0][2],
                         J[1][0], rhs[1], J[1][2],
                         J[2][0], rhs[2], J[2][2]);
    out[2] = invD * det3(J[0][0], J[0][1], rhs[0],
                         J[1][0], J[1][1], rhs[1],
                         J[2][0], J[2][1], rhs[2]);
    return true;
}

// 4x4 linear solve via partial-pivot Gaussian elimination. Returns
// false if the matrix is singular (no usable pivot). Modifies M and b
// in place; result in b.
inline bool solve_4x4_inplace(double M[4][4], double b[4]) noexcept
{
    constexpr int N = 4;
    for (int col = 0; col < N; ++col) {
        int pivot = col;
        double maxAbs = std::abs(M[col][col]);
        for (int r = col + 1; r < N; ++r) {
            if (std::abs(M[r][col]) > maxAbs) {
                maxAbs = std::abs(M[r][col]);
                pivot = r;
            }
        }
        if (maxAbs < 1e-30) return false;
        if (pivot != col) {
            for (int c = 0; c < N; ++c) std::swap(M[col][c], M[pivot][c]);
            std::swap(b[col], b[pivot]);
        }
        const double inv = 1.0 / M[col][col];
        for (int r = col + 1; r < N; ++r) {
            const double f = M[r][col] * inv;
            for (int c = col; c < N; ++c) M[r][c] -= f * M[col][c];
            b[r] -= f * b[col];
        }
    }
    for (int r = N - 1; r >= 0; --r) {
        double s = b[r];
        for (int c = r + 1; c < N; ++c) s -= M[r][c] * b[c];
        b[r] = s / M[r][r];
    }
    return true;
}

// Generic NxN partial-pivot Gaussian elimination, templated on N. Same
// algorithm as solve_4x4_inplace but works for any compile-time N.
template <std::size_t N>
inline bool solve_NxN_inplace(std::array<std::array<double, N>, N>& M,
                               std::array<double, N>& b) noexcept
{
    for (std::size_t col = 0; col < N; ++col) {
        std::size_t pivot = col;
        double maxAbs = std::abs(M[col][col]);
        for (std::size_t r = col + 1; r < N; ++r) {
            if (std::abs(M[r][col]) > maxAbs) {
                maxAbs = std::abs(M[r][col]);
                pivot = r;
            }
        }
        if (maxAbs < 1e-30) return false;
        if (pivot != col) {
            for (std::size_t c = 0; c < N; ++c) std::swap(M[col][c], M[pivot][c]);
            std::swap(b[col], b[pivot]);
        }
        const double inv = 1.0 / M[col][col];
        for (std::size_t r = col + 1; r < N; ++r) {
            const double f = M[r][col] * inv;
            for (std::size_t c = col; c < N; ++c) M[r][c] -= f * M[col][c];
            b[r] -= f * b[col];
        }
    }
    for (std::size_t r = N; r-- > 0; ) {
        double s = b[r];
        for (std::size_t c = r + 1; c < N; ++c) s -= M[r][c] * b[c];
        b[r] = s / M[r][r];
    }
    return true;
}

} // namespace detail

// Generic Newton-Raphson template. F and J are callables.
//   F(const array<N>& x, array<N>& out_F)            -> void
//   J(const array<N>& x, array<array<N>, N>& out_J)  -> void
// max_step_per_var caps each step component's absolute value to prevent
// the exponential terms (BJT, diode) from blowing up the iterate. Set
// to a sensible value for the problem (e.g., 0.2 V for circuit nodes).
template <std::size_t N, class FFunc, class JFunc>
inline nr_result<N> newton_raphson(const std::array<double, N>& x0,
                                    FFunc&& F, JFunc&& J,
                                    int    max_iter        = 30,
                                    double tol             = 1e-10,
                                    double max_step_per_var = 0.2) noexcept
{
    nr_result<N> r {};
    r.x          = x0;
    r.converged  = false;
    r.iterations = 0;

    std::array<double, N>                Fvec {};
    std::array<std::array<double, N>, N> Jmat {};

    for (int it = 0; it < max_iter; ++it) {
        F(r.x, Fvec);

        double maxF = 0.0;
        for (std::size_t i = 0; i < N; ++i) {
            const double a = std::abs(Fvec[i]);
            if (a > maxF) maxF = a;
        }
        r.final_residual_inf = maxF;
        if (maxF < tol && it > 0) {
            r.converged  = true;
            r.iterations = it;
            return r;
        }

        J(r.x, Jmat);

        // RHS = -F. Solve J * dx = -F: Cramer's rule for N = 3,
        // partial-pivot Gaussian elimination for any other N.
        std::array<double, N> bvec {};
        for (std::size_t i = 0; i < N; ++i) bvec[i] = -Fvec[i];

        bool ok;
        if constexpr (N == 3) {
            double Jraw[3][3];
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) Jraw[i][j] = Jmat[i][j];
            double rhs[3] = { bvec[0], bvec[1], bvec[2] };
            double dx3[3];
            ok = detail::solve_3x3(Jraw, rhs, dx3);
            if (ok) for (int i = 0; i < 3; ++i) bvec[i] = dx3[i];
        } else {
            ok = detail::solve_NxN_inplace<N>(Jmat, bvec);
        }
        if (!ok) {
            r.iterations = it + 1;
            return r;
        }

        for (std::size_t i = 0; i < N; ++i) {
            double dx = bvec[i];
            if (dx >  max_step_per_var) dx =  max_step_per_var;
            if (dx < -max_step_per_var) dx = -max_step_per_var;
            r.x[i] += dx;
        }
        r.iterations = it + 1;
    }

    return r;
}

// ---------------------------------------------------------------------------
// Cached LU factorisation (factor once, solve many).
//
// `solve_NxN_inplace` above factorises AND solves in a single call, destroying
// the matrix -- ideal for the per-sample Newton step where the matrix changes
// every iteration. But our LINEAR blocks (the op-amp tone/feedback networks)
// have a system matrix that only changes when a POT moves, while the
// right-hand side changes every sample. For those, factorise the matrix once
// when a knob changes (lu_factor) and then do just the cheap forward/back
// substitution per sample (lu_solve). That turns an O(N^3) elimination per
// sample into an O(N^2) substitution -- the dominant cost saving in the model.
//
// Standard Doolittle LU with partial pivoting: the unit-lower-triangular L and
// upper-triangular U are packed into one matrix M (L below the diagonal, U on
// and above it); `piv` records the row swaps.
// ---------------------------------------------------------------------------

template <std::size_t N>
struct LUFactors {
    std::array<std::array<double, N>, N> M   {};   // packed L\U
    std::array<std::size_t, N>           piv {};   // row permutation
    bool                                 ok  = false;
};

// Factorise A (taken by value -- the caller's matrix is left untouched).
// Returns ok = false if A is singular (no usable pivot).
template <std::size_t N>
inline LUFactors<N> lu_factor(std::array<std::array<double, N>, N> A) noexcept
{
    LUFactors<N> lu;
    for (std::size_t i = 0; i < N; ++i) lu.piv[i] = i;
    for (std::size_t col = 0; col < N; ++col) {
        std::size_t pivot = col;
        double maxAbs = std::abs(A[col][col]);
        for (std::size_t r = col + 1; r < N; ++r) {
            if (std::abs(A[r][col]) > maxAbs) { maxAbs = std::abs(A[r][col]); pivot = r; }
        }
        if (maxAbs < 1e-30) { lu.ok = false; return lu; }
        if (pivot != col) { std::swap(A[col], A[pivot]); std::swap(lu.piv[col], lu.piv[pivot]); }
        const double inv = 1.0 / A[col][col];
        for (std::size_t r = col + 1; r < N; ++r) {
            const double f = A[r][col] * inv;
            A[r][col] = f;                                   // store L multiplier
            for (std::size_t c = col + 1; c < N; ++c) A[r][c] -= f * A[col][c];
        }
    }
    lu.M  = A;
    lu.ok = true;
    return lu;
}

// Solve A x = b using a prior factorisation. b is overwritten with x.
template <std::size_t N>
inline void lu_solve(const LUFactors<N>& lu, std::array<double, N>& b) noexcept
{
    std::array<double, N> x {};
    for (std::size_t i = 0; i < N; ++i) x[i] = b[lu.piv[i]];          // apply pivots
    for (std::size_t i = 0; i < N; ++i) {                            // forward (L, unit diag)
        double s = x[i];
        for (std::size_t j = 0; j < i; ++j) s -= lu.M[i][j] * x[j];
        x[i] = s;
    }
    for (std::size_t i = N; i-- > 0; ) {                             // back (U)
        double s = x[i];
        for (std::size_t j = i + 1; j < N; ++j) s -= lu.M[i][j] * x[j];
        x[i] = s / lu.M[i][i];
    }
    b = x;
}

} // namespace tocake
