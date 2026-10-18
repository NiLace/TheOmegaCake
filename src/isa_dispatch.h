// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// ISA dispatch for the one hot loop: one block through the oversampler and the
// engine, compiled twice and chosen at load time by the CPU.
//
// WHY THIS SHAPE, and it is the whole design decision.
//
// The baseline of everything published here stays generic x86-64: a binary that
// emits AVX2 globally does not run slowly on a CPU without it, it does not START,
// and it takes the host down with it. A run-time dispatch is the other thing --
// the default clone is the generic one, and only a CPU that has the instructions
// selects the fast one.
//
// WHERE THE BOUNDARY GOES. A dispatch accelerates only what it clones, and it
// does not cover what the cloned function CALLS. So the boundary is put around
// the BLOCK, not the sample:
//
//   * everything hot is INSIDE. `Oversampler::process` and the engine's
//     `process_sample` are templates, so both inline whole into each clone and
//     the entire per-sample path is compiled at the clone's ISA.
//   * the indirection is paid once per BLOCK (about 375 times a second at 128
//     samples and 48 kHz), not once per sample. A dispatch on a per-sample
//     function also prevents inlining, which would eat most of the gain.
//
// Confined to the hot function, the gain matches a whole-build
// `-march=x86-64-v3` on the DK engines. Any function the hot path CALLS out of
// line stays at the default ISA and loses its share.
//
// NOT BIT-IDENTICAL, and that is by construction: GCC contracts `a*b+c` into an
// FMA where the target has one, so the v3 clone's rounding differs from the
// default clone's. The audio moves about -169 dB, so a bit-exact comparison has
// to run the DEFAULT clone (e.g. under `qemu-x86_64 -cpu Nehalem`).
//
// Verification is an instruction census, not the clock: post-Haswell
// instructions inside the v3 clone and ZERO outside it, plus a run on emulated
// CPUs with a canary that must fault -- a build whose clone never executes and
// a build with no clone at all time the same.
#pragma once

#include "isa_hot.h"
#include "oversampler.h"

#include <cstdint>

namespace tocake {

#if (defined(__x86_64__) || defined(_M_X64)) && defined(__GNUC__)

// The guard. Read ONCE at load time, not per block: `__builtin_cpu_supports`
// compiles to a load and a test, but putting it behind a constant lets the
// branch below predict perfectly and keeps the two loops from sharing code.
// `__builtin_cpu_init()` first: the feature words live in libgcc's
// `__cpu_model`, and reading them before it is filled in answers "no feature at
// all" -- which is safe but would ship the dispatch permanently switched off,
// and nothing would say so; only a run on a CPU that DOES have the ISA shows it.
// The guard asks for EXACTLY the level the clones are compiled for
// (`target("arch=x86-64-v3")`), which licenses BMI1/BMI2/LZCNT/MOVBE/F16C as
// well as AVX2 and FMA: a CPU (or a virtual one) with AVX2+FMA but without the
// rest must take the base path, or a future clone could fault on it.
inline bool isa_v3_detect() noexcept
{
    __builtin_cpu_init();
    return __builtin_cpu_supports("x86-64-v3");
}
inline const bool kIsaV3 = isa_v3_detect();
#endif

// Run `eng` over `n` oversampled samples in place, at the best instruction set
// this CPU has.
//
// The branch is at the SPAN, so it is taken once per chunk (once per host
// block in practice), not once per sample.
//
// The LOOP belongs to the engine, not to this function, and that is the point:
// whether the per-sample body should be inlined into the loop or called out of
// line is a per-engine MEASUREMENT, so each engine owns its own pair and this
// function only chooses between them. See isa_hot.h.
template <class Engine>
inline void run_engine_span(float* buf, uint32_t n, Engine* eng) noexcept
{
#if (defined(__x86_64__) || defined(_M_X64)) && defined(__GNUC__)
    if (kIsaV3) { eng->span_run_isa_v3(buf, n); return; }
#endif
    eng->span_run(buf, n);
}

// One block, in place, through the oversampler and `eng`.
template <class Engine>
inline void run_block(Oversampler& os, float* buf, uint32_t n, Engine* eng) noexcept
{
    os.process_span(buf, buf, n,
                    [eng](float* b, uint32_t c) noexcept { run_engine_span(b, c, eng); });
}

}  // namespace tocake
