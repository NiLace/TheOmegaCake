// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// The three attributes the ISA dispatch is built out of, in one place so the
// engines agree on them. The design and its measurements live in
// src/isa_dispatch.h; this header is only the spelling.
#pragma once

// A per-sample entry compiled for x86-64-v3 (AVX2 + FMA + BMI). Anything
// marked with this may only ever be CALLED behind a run-time guard: the
// instructions are not conditional, the guard is.
#if (defined(__x86_64__) || defined(_M_X64)) && defined(__GNUC__)
#  define TOCAKE_SPAN_ISA  __attribute__((target("arch=x86-64-v3"), noinline))
#else
#  define TOCAKE_SPAN_ISA
#endif

// Its baseline twin. `noinline` is deliberate and it is a measurement, not a
// habit: letting the engine inline into the caller's span loop is slower on the
// DK engines. A 25 kB body belongs in its own function.
#define TOCAKE_SPAN_BASE __attribute__((noinline))

// The body the two entries share. Forced inline so each entry gets its OWN copy
// compiled at its own instruction set -- that duplication IS the dispatch.
#define TOCAKE_SPAN_BODY __attribute__((always_inline)) inline

// A whole-SPAN entry compiled for x86-64-v3, for an engine whose per-sample work
// is better inlined into the loop than called. Same rule as above: only ever
// reached behind the run-time guard.
//
// Which of the two an engine wants is a MEASUREMENT, not a style. The DK
// engines lose time when their 25 kB body is inlined into the loop, so they
// take the per-sample pair; the op-amp GAINS from it and doubles its share of
// the ISA, so it takes the span pair.
#if (defined(__x86_64__) || defined(_M_X64)) && defined(__GNUC__)
#  define TOCAKE_SPAN_RUN_ISA __attribute__((target("arch=x86-64-v3")))
#else
#  define TOCAKE_SPAN_RUN_ISA
#endif
