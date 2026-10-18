// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// Flush-to-zero / denormals-are-zero guard.
//
// IIR filter and distortion tails decay toward zero and produce denormal
// floats; on x86 denormal arithmetic is 10-100x slower, which shows up as
// CPU spikes on note decays / silence. Flushing denormals to zero is
// standard audio-DSP practice and inaudible (denormals are < 1e-308, ~300 dB
// below any signal).
//
// THE REGISTER IS PER-THREAD, NOT PER-INSTANCE, and that is why setting it
// and walking away is a defect rather than a shortcut. MXCSR belongs to the
// audio thread, not to this plugin: leaving flush-to-zero on changes the
// floating-point behaviour of everything that runs after us on that thread --
// the other plugins in the chain, and the host's own summing -- which is an
// undeclared side effect on shared state. A plugin does not get to choose the
// arithmetic mode of its neighbours.
//
// So the mode is saved on entry and restored on exit. The cost is a handful of
// nanoseconds once per processing call, against the tens of microseconds the
// block itself takes.
//
// Restoring is only correct if the restore happens AFTER the audio work.
// RAII gives that here, but it is verified rather than assumed: a long silent
// tail must still produce no denormal output, matching a build that never
// restores.

#pragma once

#if defined(__x86_64__) || defined(__i386__)
#include <xmmintrin.h>   // _MM_SET_FLUSH_ZERO_MODE, _mm_getcsr
#include <pmmintrin.h>   // _MM_SET_DENORMALS_ZERO_MODE
#endif

namespace tocake {

// For callers that genuinely own the thread: offline drivers, benches, tests.
inline void disable_denormals() noexcept {
#if defined(__x86_64__) || defined(__i386__)
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#endif
}

// What an audio callback should use: flush denormals for the duration of the
// block, then hand the thread back exactly as it was found.
class DenormalGuard {
public:
    DenormalGuard() noexcept {
#if defined(__x86_64__) || defined(__i386__)
        saved_ = _mm_getcsr();
        disable_denormals();
#endif
    }
    ~DenormalGuard() noexcept {
#if defined(__x86_64__) || defined(__i386__)
        _mm_setcsr(saved_);
#endif
    }
    DenormalGuard(const DenormalGuard&)            = delete;
    DenormalGuard& operator=(const DenormalGuard&) = delete;

private:
#if defined(__x86_64__) || defined(__i386__)
    unsigned int saved_ = 0;
#endif
};

}  // namespace tocake
