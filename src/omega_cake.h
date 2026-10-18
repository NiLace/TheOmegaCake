// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// Top-level orchestrator of the four-stage transistor fuzz.
//
// Owns the full-pedal joint solver (StandardEngine) and presents a simple
// knob-based API: set_sustain/tone/volume take the user-facing knob angle
// (0..1), apply the variant's taper to get the physical wiper fraction,
// and forward that to the solver. Same code path for every board
// revision -- pass a different CakeVariant to model a different pedal.
//
// The DSP itself is the StandardEngine (standard_engine.h: the DK engine of
// dk_standard.h);
// this class is the thin orchestration + taper layer that the LV2 wrapper
// and the offline harness both drive.

#pragma once

#include "cake_variant.h"
#include "dk_engine.h"
#include "standard_engine.h"

#include <cstdint>

namespace tocake {

class OmegaCake {
public:
    using T = double;

    explicit OmegaCake(const CakeVariant& v) noexcept
        : variant_(v), pedal_(v) {}

    void prepare(double sample_rate) { pedal_.prepare(sample_rate); }
    void reset() noexcept            { pedal_.reset(); }

    // Knob positions: 0.0 = fully CCW, 1.0 = fully CW. Taper maps the knob
    // angle to the physical wiper fraction the circuit solves at.
    void set_sustain(double knob) noexcept {
        pedal_.set_sustain(apply_taper(knob, variant_.taper_sustain));
    }
    void set_coarse(bool on) noexcept { dk_set_coarse(pedal_, on, 0); }
    void set_tone(double knob) noexcept {
        pedal_.set_tone(apply_taper(knob, variant_.taper_tone));
    }
    void set_volume(double knob) noexcept {
        pedal_.set_volume(apply_taper(knob, variant_.taper_volume));
    }

    const char* variant_name() const noexcept { return variant_.name; }

    // Process one sample (input in volts at the jack, output in volts at
    // the output jack).
    float process_sample(float vin_volts) noexcept {
        return pedal_.process_sample(vin_volts);
    }

    // The ISA dispatch's entries, forwarded to the inner engine's (which are the
    // ones that carry the attribute). These forwarders vanish; the engine's own
    // pair stays out of line.
    __attribute__((always_inline)) inline float span_sample(float vin_volts) noexcept {
        return pedal_.span_sample(vin_volts);
    }
    __attribute__((always_inline)) inline float span_sample_isa_v3(float vin_volts) noexcept {
        return pedal_.span_sample_isa_v3(vin_volts);
    }

    void span_run(float* buf, uint32_t n) noexcept { pedal_.span_run(buf, n); }
    void span_run_isa_v3(float* buf, uint32_t n) noexcept { pedal_.span_run_isa_v3(buf, n); }

    // Convergence telemetry forwarded from the inner solver.
    long   dbg_iters()     const noexcept { return pedal_.dbg_iters();     }
    long   dbg_calls()     const noexcept { return pedal_.dbg_calls();     }
    long   dbg_nonconv()   const noexcept { return pedal_.dbg_nonconv();   }
    double dbg_max_resid() const noexcept { return pedal_.dbg_max_resid(); }

private:
    const CakeVariant& variant_;
    StandardEngine        pedal_;
};

} // namespace tocake
