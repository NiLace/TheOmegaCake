// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors

#pragma once
//
// Shared non-finite guard for the per-sample solvers.
//
// Two properties every engine's guard must have:
//
//   * COVERAGE. Test ALL N unknowns, not only the output node. A solve that
//     diverges rarely does so at the output first: the poison appears at an
//     internal node and only reaches the output a sample later, if at all.
//     Between those two moments the engine would run on a state that is
//     already garbage, and nothing would say so.
//
//   * THE ANCHOR. Recover to the full DC rest point -- an ABSOLUTE state --
//     not by rewriting the unknowns alone: that leaves the capacitor states
//     untouched, a non-finite capacitor state survives it and re-poisons the
//     very next sample, so the guard fires for ever and the engine never comes
//     back. For the same reason a fallback must never restore "the last good
//     state": restoring what already failed perpetuates the failure.
//
// The check is a read over N doubles once per sample, off the iteration loop,
// so it is not on the critical path of the Newton solve.
//
#include <array>
#include <cmath>
#include <cstddef>

namespace tocake {

// True when every unknown is finite. Checked in full: a single index cannot
// stand in for the vector.
template <std::size_t N>
inline bool all_finite(const std::array<double, N>& v) noexcept
{
    for (std::size_t i = 0; i < N; ++i)
        if (!std::isfinite(v[i])) return false;
    return true;
}

// For the stored states, which are plain members rather than an array.
inline bool all_finite(std::initializer_list<double> xs) noexcept
{
    for (double x : xs)
        if (!std::isfinite(x)) return false;
    return true;
}

}  // namespace tocake
