// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// Integer oversampler for the nonlinear DSP (anti-aliasing).
//
// Hard clipping generates harmonics above Nyquist that fold back as aliasing
// at a 48 kHz host rate. This runs the model at N x the host rate so those
// harmonics land above the (oversampled) Nyquist, then a steep low-pass on
// the way back down removes them before decimation.
//
// Polyphase windowed-sinc FIR (Kaiser window), N = 1..MAX_FACTOR (the product
// uses 2, 3 and 4). Latency: latency_samples() host-rate samples (TPP - 1).
// Buffers are sized for the MAX factor at prepare() time, so changing factor
// at runtime only recomputes coefficients (no allocation on the audio thread).
//
// process(in, out, n, proc): upsample -> run `proc` per oversampled sample
// -> downsample. At N=1 it is a bit-transparent pass-through (proc per
// input sample): the offline tools' 1x path. The product runs 2x, 3x or 4x
// (src/os_policy.h).

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

// The filter: TPP taps per polyphase branch, cutoff at OS_FC x the HOST
// Nyquist, Kaiser beta OS_BETA. With the cutoff AT host Nyquist the transition
// band is symmetric about it: at 44.1 kHz the passband reaches 20 kHz and what
// folds back lands above 20 kHz. 32 taps make that transition narrow enough.
// Measured through the plugin at 44.1 kHz, the fuzz's harmonics at 18 / 20 kHz
// against the same engine unfiltered: -0.30 / -0.50 dB. A transition that
// starts below 20 kHz rolls off the top octave, and an ANMR against SPICE does
// not see such a roll-off. Worst folded-back gain into 0-20 kHz: -41 dB at
// 44.1, -58 at 48. Latency TPP - 1 = 31 host samples.
#ifndef TOCAKE_OS_TPP
#define TOCAKE_OS_TPP 32
#endif
#ifndef TOCAKE_OS_FC
#define TOCAKE_OS_FC 1.00      // cutoff as a fraction of the HOST Nyquist
#endif
#ifndef TOCAKE_OS_BETA
#define TOCAKE_OS_BETA 5.0
#endif

namespace tocake {

class Oversampler {
public:
    static constexpr int MAX_FACTOR = 4;
    static constexpr int TPP        = TOCAKE_OS_TPP;   // FIR taps per polyphase branch (a power of two)
    static_assert((TPP & (TPP - 1)) == 0, "TPP must be a power of two (the up-ring wraps with a mask)");
    static_assert(TPP % 16 == 0, "the dot products take 8 lanes, the symmetric one over half of TPP * N taps");
    // The downsampler's history ring: a power of two at least TPP * MAX_FACTOR
    // long, the same for every factor, so the wrap is a mask for 3x too.
    static constexpr int kDownRing  = [] { int r = 1; while (r < TPP * MAX_FACTOR) r <<= 1; return r; }();

    void prepare(int max_block) {
        max_chunk_ = (max_block > 0) ? max_block : 8192;
        os_buf_.assign(static_cast<std::size_t>(MAX_FACTOR) * max_chunk_, 0.0f);
        const int Lmax = TPP * MAX_FACTOR;
        up_coeffs_.assign(Lmax, 0.0f);
        down_coeffs_.assign(Lmax, 0.0f);
        up_branch_.assign(Lmax, 0.0f);
        up_ring_.assign(2 * TPP, 0.0f);          // doubled: see process_chunk
        up_head_ = 0;
        down_ring_.assign(2 * kDownRing, 0.0f);  // doubled: see process_chunk
        set_factor(N_);
    }

    // Recompute filters for a new factor (1..MAX_FACTOR). No allocation.
    void set_factor(int f) {
        if (f < 1) f = 1;
        if (f > MAX_FACTOR) f = MAX_FACTOR;
        N_ = f;
        L_ = TPP * N_;
        if (N_ > 1) design_filters();
        reset();
    }

    int factor()          const noexcept { return N_; }
    // Both filters count, not one: (L_-1)/(2*N_) would be only one of them.
    //
    // The chain is upsample THEN downsample, and each is a linear-phase FIR of
    // length L at the oversampled rate, so each contributes (L-1)/2 oversampled
    // samples. The downsampler emits its host sample at oversampled index
    // m*N + (N-1), which gives back (N-1)/N, so the host-rate total is
    //
    //     (L-1)/N - (N-1)/N  =  (L - N)/N  =  L/N - 1  =  TPP - 1
    //
    // -- independent of the factor. Checked by measurement as well, on a
    // 16-tap build (an impulse through process() with proc = identity, no
    // circuit in the path): peak delay and energy centroid both read 15.00 =
    // TPP - 1 at N = 2 and at N = 4, with sum(out) = 1.0. N = 1 is a bit-exact passthrough and really is 0.
    //
    // This is the one number a null can never check: every null here aligns
    // by cross-correlation, so it normalises away the quantity being declared.
    // A wrong value leaves a host doing delay compensation out by the
    // difference -- inaudible alone, a comb filter in any parallel or blended
    // path.
    int latency_samples() const noexcept { return (N_ > 1) ? (L_ / N_ - 1) : 0; }

    void reset() noexcept {
        std::fill(up_ring_.begin(),   up_ring_.end(),   0.0f);
        up_head_ = 0;
        std::fill(down_ring_.begin(), down_ring_.end(), 0.0f);
        down_head_ = 0;
    }

    // proc: float(float) -- processes one (oversampled) sample.
    template <class Proc>
    void process(const float* in, float* out, uint32_t n, Proc&& proc) {
        process_span(in, out, n, [&proc](float* b, uint32_t c) noexcept {
            for (uint32_t j = 0; j < c; ++j) b[j] = proc(b[j]);
        });
    }

    // span: void(float* buf, uint32_t count) -- processes `count` oversampled
    // samples IN PLACE.
    //
    // This form exists for the ISA dispatch (src/isa_dispatch.h), and the
    // reason is the difference between one indirect call per BLOCK and one per
    // SAMPLE. A dispatched function cannot be inlined into its caller, so if the
    // engine is reached through a per-sample callable the resolver's call lands
    // on the hot path 96 000 to 384 000 times a second and eats the gain it is
    // there to collect. Handing the engine a SPAN moves that call to the chunk
    // boundary, which is once per host block.
    //
    // The per-sample `process` above is written in terms of this one, so
    // there is a single implementation of the resampling rather than two that
    // can drift apart. Same operations in the same order ⇒ bit-identical.
    template <class Span>
    void process_span(const float* in, float* out, uint32_t n, Span&& span) {
        if (N_ == 1) {                       // bypass: exact, no resampling
            if (in != out) std::copy(in, in + n, out);
            span(out, n);
            return;
        }
        uint32_t off = 0;
        while (off < n) {
            const uint32_t c = std::min<uint32_t>(max_chunk_, n - off);
            process_chunk(in + off, out + off, c, span);
            off += c;
        }
    }

private:
    template <class Span>
    void process_chunk(const float* in, float* out, uint32_t c, Span&& span) {
        const int N = N_, L = L_;
        // Both histories are DOUBLED rings: every sample is written at `head`
        // and at `head + size`, so the last `taps` samples are always one
        // contiguous run ending at `head + size`, oldest first. A contiguous
        // dot product is what a compiler can turn into SIMD; a masked index,
        // `(head - k) & mask`, is not.
        //
        // Each dot product keeps kLanes independent partial sums, summed at the
        // end: a single accumulator makes every add wait for the previous one,
        // and without -ffast-math the compiler may not reorder them itself. The
        // order of the additions is therefore not the direct form's, and the
        // output differs from it in the last bits (~-140 dB), not in sound.
        //
        // Upsample: each input sample -> N oversampled samples (polyphase).
        // up_branch_ holds branch p's taps contiguously, oldest-first, so
        // branch p is one straight dot product against the history window.
        const int umask = TPP - 1;
        for (uint32_t m = 0; m < c; ++m) {
            up_head_ = (up_head_ + 1) & umask;
            up_ring_[up_head_] = in[m];
            up_ring_[up_head_ + TPP] = in[m];
            const float* x = &up_ring_[up_head_ + 1];          // TPP samples, oldest first
            for (int p = 0; p < N; ++p)
                os_buf_[m * N + p] = dot(&up_branch_[p * TPP], x, TPP);
        }
        // Run the nonlinear model at the oversampled rate.
        const uint32_t os_n = c * static_cast<uint32_t>(N);
        span(os_buf_.data(), os_n);
        // Downsample: anti-alias FIR, one output per N oversampled samples. The
        // filter is linear-phase (h[k] = h[L-1-k]): the two samples that share a
        // tap are added first and multiplied once, half the multiplies.
        const int mask = kDownRing - 1;
        for (uint32_t m = 0; m < c; ++m) {
            for (int p = 0; p < N; ++p) {
                down_head_ = (down_head_ + 1) & mask;
                const float v = os_buf_[m * N + p];
                down_ring_[down_head_] = v;
                down_ring_[down_head_ + kDownRing] = v;
            }
            const float* x = &down_ring_[down_head_ + kDownRing - L + 1];   // L samples, oldest first
            out[m] = dot_symmetric(down_coeffs_.data(), x, L);
        }
    }

    static constexpr int kLanes = 8;

    // sum c[k] * x[k], k < n (n a multiple of kLanes).
    static float dot(const float* __restrict c, const float* __restrict x, int n) noexcept {
        float acc[kLanes] = {};
        for (int k = 0; k < n; k += kLanes)
            for (int j = 0; j < kLanes; ++j) acc[j] += c[k + j] * x[k + j];
        return reduce(acc);
    }

    // sum h[k] * x[k] over k < n for a symmetric h (n even, n/2 a multiple of
    // kLanes): sum h[k] * (x[k] + x[n-1-k]) over k < n/2.
    static float dot_symmetric(const float* __restrict h, const float* __restrict x, int n) noexcept {
        float acc[kLanes] = {};
        const int half = n / 2;
        for (int k = 0; k < half; k += kLanes)
            for (int j = 0; j < kLanes; ++j) acc[j] += h[k + j] * (x[k + j] + x[n - 1 - k - j]);
        return reduce(acc);
    }

    static float reduce(const float (&a)[kLanes]) noexcept {
        return ((a[0] + a[4]) + (a[1] + a[5])) + ((a[2] + a[6]) + (a[3] + a[7]));
    }

    void design_filters() {
        const double fc   = 0.5 / N_ * TOCAKE_OS_FC;   // cutoff (norm. to OS rate)
        const double beta = TOCAKE_OS_BETA;            // Kaiser window shape
        const double i0b  = i0(beta);
        const double pi   = 3.14159265358979323846;
        double sum = 0.0;
        for (int k = 0; k < L_; ++k) {
            const double mc = k - (L_ - 1) / 2.0;
            const double s  = (std::abs(mc) < 1e-9)
                                  ? 2.0 * fc
                                  : std::sin(2.0 * pi * fc * mc) / (pi * mc);
            const double r  = 2.0 * k / (L_ - 1) - 1.0;
            const double w  = i0(beta * std::sqrt(std::max(0.0, 1.0 - r * r))) / i0b;
            const double h  = s * w;
            down_coeffs_[k] = static_cast<float>(h);
            sum += h;
        }
        const double inv = 1.0 / sum;          // DC -> unity
        for (int k = 0; k < L_; ++k) {
            down_coeffs_[k] = static_cast<float>(down_coeffs_[k] * inv);
            up_coeffs_[k]   = static_cast<float>(down_coeffs_[k] * N_); // DC gain N
        }
        // Branch p, tap k (k = 0 newest) is up_coeffs_[k*N + p]; stored oldest
        // first so it lines up with the history window in process_chunk.
        for (int p = 0; p < N_; ++p)
            for (int k = 0; k < TPP; ++k)
                up_branch_[p * TPP + (TPP - 1 - k)] = up_coeffs_[k * N_ + p];
    }

    static double i0(double x) {               // modified Bessel I0 (series)
        double sum = 1.0, t = 1.0;
        for (int k = 1; k < 40; ++k) {
            t *= (x * x) / (4.0 * k * k);
            sum += t;
            if (t < 1e-12 * sum) break;
        }
        return sum;
    }

    int N_ = 2, L_ = TPP * 2;
    int max_chunk_ = 8192;
    int down_head_ = 0;                    // ring-buffer write index (downsampler)
    std::vector<float> up_coeffs_, down_coeffs_, up_branch_, up_ring_, down_ring_, os_buf_;
    int up_head_ = 0;
};

} // namespace tocake
