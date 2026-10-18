// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// Variant data: every value that changes between board revisions of the circuit.
//
// The DSP code is shared across every board revision; only this
// struct's contents differ. Adding a new variant (another 1970s board
// revision, a Russian-made version, a modern reissue, ...) is a data-entry exercise -- write another
// constexpr literal of this type.
//
// Values for the V3 below mirror the verified KiCad netlist
// at netlists/TheCake.net. Device-model
// parameters match the idealised models in spice/cake_v3.cir, the
// SPICE reference the engines null against.

#pragma once

namespace tocake {

// -- Pot taper ---------------------------------------------------------------
//
// Taper maps the knob angle the user turns (0..1) to the PHYSICAL wiper
// fraction the circuit sees (0..1). It is a layer ON TOP of the validated
// circuit solve, never inside it: the DSP only ever cares about the wiper
// fraction. The SPICE references were generated at physical fractions, so
// Linear (knob == fraction) reproduces the validated behaviour exactly.
//
// The original pedals typically use audio/log tapers on Volume and Sustain and a
// linear taper on Tone, but this varies by revision -- set per-variant once
// confirmed against the schematic. Defaults are Linear so an un-annotated
// variant behaves exactly as validated.
enum class PotTaper { Linear, Log, AntiLog };

// Map knob angle (0..1) to physical wiper fraction (0..1). Log/AntiLog use
// a simple power-law approximation of an audio taper (refine the curve if a
// closer match to a specific pot is needed).
inline double apply_taper(double knob, PotTaper taper) noexcept {
    if (knob < 0.0) knob = 0.0;
    if (knob > 1.0) knob = 1.0;
    switch (taper) {
        case PotTaper::Log:     return knob * knob;                 // slow start
        case PotTaper::AntiLog: return 1.0 - (1.0 - knob) * (1.0 - knob);
        case PotTaper::Linear:  break;
    }
    return knob;
}

// -- Device model parameters ------------------------------------------------

// Ebers-Moll-level NPN BJT parameters. The bias-solver uses these to
// find the DC operating point; nothing derived from it is stored in the
// variant.
struct bjt_params {
    double Is;     // reverse saturation current (A)
    double Bf;     // forward common-emitter current gain
    double Vaf;    // forward Early voltage (V)  -- gives the r_o we need
};

// Shockley diode parameters, plus a count per direction. The transistor
// engines model ONE diode per direction (or none: 0/0 on the circuits without
// the pair), so a row may only say 1/1 or 0/0 -- anything else would be
// silently ignored. scripts/gen_variants.py only ever writes those two.
struct diode_params {
    double Is;          // reverse saturation current (A)
    double N;           // emission coefficient (unitless)
    int    count_fwd;   // diodes in the "forward" antiparallel branch
    int    count_rev;   // diodes in the "reverse" antiparallel branch
};

// -- Per-stage parameter groups --------------------------------------------

// Common-emitter NPN with collector-to-base feedback (R + C in
// parallel). No diode. Used as both the input booster and the
// output booster, with different component values.
struct booster_params {
    bjt_params q;

    // Resistor network
    double R_input_series;     // base series resistor (R2 / -)
    double R_base_to_gnd;      // base pulldown (R14 / R3)
    double R_feedback;         // collector-to-base feedback R (R9 / -)
    double R_collector_load;   // collector to V+ (R13 / R6)
    double R_emitter;          // emitter to GND (R22 / R4)

    // Optional second bias resistor (output booster uses it; input
    // booster uses a single feedback R instead and sets this to 0).
    double R_bias_top;         // base to V+ (- / R7)

    // Caps
    double C_input_couple;     // input DC block (C1 / C3)
    double C_feedback_hf;      // HF roll-off across R_feedback (C10 / -)
    double C_output_couple;    // output DC block (C4 / C2)

    // Optional bespoke output-stage parts (0 = absent on every standard
    // variant). Used only by the simplified-tone mod family (Vm Precision /
    // Vm German Cake) via src/simple_tone_cake.h; the shared FullPedalV3 and
    // every other solver leave these 0 and never read them.
    double R_emitter_bypass = 0.0;   // R28: series R of the emitter AC-bypass leg
    double C_emitter_bypass = 0.0;   // C15: series cap of that leg (to GND)
    double R_output_load    = 0.0;   // R29: shunt load on the output couple node
};

// Common-emitter NPN with collector-to-base feedback PLUS an
// antiparallel diode pair (optionally AC-coupled via a series cap)
// inside that feedback path. The clipper heart of the circuit.
struct clipper_params {
    bjt_params   q;
    diode_params diodes;

    double R_input_series;     // (R19 / R12)
    double R_base_to_gnd;      // (R20 / R16)
    double R_feedback;         // (R17 / R15)
    double R_collector_load;   // (R18 / R11)
    double R_emitter;          // (R21 / R10)

    double C_feedback_hf;      // (C12 / C11)
    double C_diode_couple;     // AC-coupling for the diode pair (C6 / C7).
                               // Set to 0 to DC-couple instead -- some
                               // Russian-made variants do this.
    double C_output_couple;    // to next stage (C13 on clipper_1; the
                               // V3's clipper_2 has 0 -- tone stack is
                               // wired directly).
};

// Inter-stage attenuator (post-input-booster, pre-clipper_1).
// 100k pot in series, with a 1k pulldown on the bottom end and a
// 100n coupling cap from the wiper to the next stage.
struct sustain_params {
    double R_pot;              // total pot resistance (Sustain1)
    double R_pot_bottom;       // pulldown from pot bottom to GND (R23)
    double C_wiper_couple;     // coupling cap on the wiper (C5)
};

// Passive Baxandall-scoop tone stack. Topology fixed across all
// known board revisions; only component values vary.
struct tone_params {
    double R_pot;              // total pot resistance (Tone1)

    // HP branch (series cap, shunt R)
    double C_hp_series;        // signal -> pot_top  (C9)
    double R_hp_shunt;         // pot_top -> GND     (R5)

    // LP branch (series R, shunt cap)
    double R_lp_series;        // signal -> pot_bot  (R8)
    double C_lp_shunt;         // pot_bot -> GND     (C8)
};

// Output volume pot. Just a divider; no extra components in V3 (and
// likely all variants).
struct volume_params {
    double R_pot;              // total pot resistance (Volume1)
};

// -- Top-level variant struct ----------------------------------------------

struct CakeVariant {
    const char* name;

    double vcc;                     // supply rail (V). 9.0 for every known
                                    // board; carried in the variant so
                                    // future oddballs (12V mods, sagged-
                                    // battery emulations) need no code change.

    booster_params  input_booster;
    sustain_params  sustain;
    clipper_params  clipper_1;
    clipper_params  clipper_2;
    tone_params     tone;
    booster_params  output_booster;
    volume_params   volume;

    // Knob->wiper-fraction taper per pot. Default Linear == validated
    // behaviour; set the real tapers per-variant once confirmed against
    // the schematic. See PotTaper / apply_taper above.
    PotTaper taper_sustain = PotTaper::Linear;
    PotTaper taper_tone    = PotTaper::Linear;
    PotTaper taper_volume  = PotTaper::Linear;
};

// -- V3 variant ---------------------------------------------------
//
// Values from netlists/TheCake.net, cross-referenced against
// spice/cake_v3.cir component-by-component.
//
// Transistor: 2N5088 (the V3's NPN). The .q params below are a 3-param
// idealised reduction of the Fairchild 2N5088 datasheet at this circuit's
// operating point (Ic ~ 0.16-0.43 mA): Bf=540 effective DC beta there,
// Is=5.9e-15 from the datasheet Vbe, Vaf=62 from Fairchild's printed SPICE
// model. All four stages are the same part, so they share these params --
// the per-stage operating-point differences come from the bias solver, not
// from the device. See the V9 literal below for the BC550C reduction.

inline constexpr CakeVariant V3 = {
    .name = "V3",
    .vcc  = 9.0,

    .input_booster = {
        .q = { .Is = 5.9e-15, .Bf = 540.0, .Vaf = 62.0 },
        .R_input_series   = 39.0e3,
        .R_base_to_gnd    = 47.0e3,
        .R_feedback       = 470.0e3,
        .R_collector_load = 10.0e3,
        .R_emitter        = 100.0,
        .R_bias_top       = 0.0,            // unused: feedback-biased
        .C_input_couple   = 1.0e-6,
        .C_feedback_hf    = 470.0e-12,
        .C_output_couple  = 1.0e-6,
    },

    .sustain = {
        .R_pot         = 100.0e3,
        .R_pot_bottom  = 1.0e3,
        .C_wiper_couple = 100.0e-9,
    },

    .clipper_1 = {
        .q      = { .Is = 5.9e-15, .Bf = 540.0, .Vaf = 62.0 },
        .diodes = { .Is = 2.52e-9, .N  = 1.752,
                    .count_fwd = 1, .count_rev = 1 },
        .R_input_series   = 10.0e3,
        .R_base_to_gnd    = 100.0e3,
        .R_feedback       = 470.0e3,
        .R_collector_load = 10.0e3,
        .R_emitter        = 150.0,
        .C_feedback_hf    = 470.0e-12,
        .C_diode_couple   = 1.0e-6,
        .C_output_couple  = 100.0e-9,       // C13 to next stage
    },

    .clipper_2 = {
        .q      = { .Is = 5.9e-15, .Bf = 540.0, .Vaf = 62.0 },
        .diodes = { .Is = 2.52e-9, .N  = 1.752,
                    .count_fwd = 1, .count_rev = 1 },
        .R_input_series   = 10.0e3,
        .R_base_to_gnd    = 100.0e3,
        .R_feedback       = 470.0e3,
        .R_collector_load = 10.0e3,
        .R_emitter        = 150.0,
        .C_feedback_hf    = 470.0e-12,
        .C_diode_couple   = 1.0e-6,
        .C_output_couple  = 0.0,            // V3 wires tone stack directly
    },

    .tone = {
        .R_pot       = 100.0e3,
        .C_hp_series = 4.0e-9,
        .R_hp_shunt  = 22.0e3,
        .R_lp_series = 39.0e3,
        .C_lp_shunt  = 10.0e-9,
    },

    .output_booster = {
        .q = { .Is = 5.9e-15, .Bf = 540.0, .Vaf = 62.0 },
        .R_input_series   = 0.0,            // unused: input is the tone wiper
        .R_base_to_gnd    = 100.0e3,        // R3
        .R_feedback       = 0.0,            // no feedback resistor on output booster
        .R_collector_load = 15.0e3,         // R6
        .R_emitter        = 3.3e3,          // R4
        .R_bias_top       = 430.0e3,        // R7  (the KiCad netlist's value; some V3 board revs show 470k)
        .C_input_couple   = 100.0e-9,       // C3
        .C_feedback_hf    = 0.0,            // unused
        .C_output_couple  = 100.0e-9,       // C2
    },

    .volume = {
        .R_pot = 100.0e3,
    },

    // Standard tapers of the circuit: audio/log on Volume and Sustain, linear on
    // Tone. The log curve here is a power-law approximation -- verify the
    // exact tapers against the V3 schematic and refine apply_taper() if a
    // closer match to the real pots is wanted. (Taper is knob->fraction
    // only; the SPICE-validated physics is unaffected.)
    .taper_sustain = PotTaper::Log,
    .taper_tone    = PotTaper::Linear,
    .taper_volume  = PotTaper::Log,
};

// -- V9 variant (NYC "Classic" Reissue, EC3003-C / 2008+) ---------------------
//
// The V9 is the SAME topology as the V3/V6 family ("same tonal family as V6";
// V6 = the Tone-Bypass V3 reworked). So this literal is the V3 with only the
// EC3003-C deltas applied -- everything not listed below is inherited from the
// V3 and should be cross-checked against a real EC3003-C schematic before this
// is treated as final.
//
// EC3003-C (2008+) deltas vs V3, from a public reference of the board:
//   * Transistor: BC550C (V3: 2N5088). 3-param idealised reduction from the
//     NXP/ON datasheet at our operating point (Ic ~ 0.2-0.45 mA):
//       Bf = 500     -- DC hFE, flat ~500 across 0.1-10 mA (NXP Fig.2; typ 520
//                       at 2 mA). Far below the 90 mA Ikf corner, so no rolloff.
//       Is = 1.5e-14 -- from datasheet Vbe(typ)=660 mV at Ic=2 mA, Vce=5 V.
//       Vaf = 53     -- from NXP's own BC547C SPICE card (same die, C-grade);
//                       datasheets state no Early voltage. ro >> Rc, weak effect.
//     vs the 2N5088: ~same beta, but ~2.5x larger Is -> idles ~24 mV lower Vbe.
//   * Tone HP shunt R5: 22k -> 27k.   Tone LP series R8: 39k -> 27k.
//   * Volume pot: 100k -> 500k (log).  (EC3003-C "thicker/bassier" voicing.)
//   * Diodes: 1N4148 (2008+) -- IDENTICAL to the V3 model, so unchanged here.
//     (The earlier 2000-07 EC3003 used 1N6263 Schottky -- a different variant.)
inline constexpr CakeVariant V9 = {
    .name = "V9",
    .vcc  = 9.0,

    .input_booster = {
        .q = { .Is = 1.5e-14, .Bf = 500.0, .Vaf = 53.0 },   // BC550C
        .R_input_series   = 39.0e3,
        .R_base_to_gnd    = 47.0e3,
        .R_feedback       = 470.0e3,
        .R_collector_load = 10.0e3,
        .R_emitter        = 100.0,
        .R_bias_top       = 0.0,
        .C_input_couple   = 1.0e-6,
        .C_feedback_hf    = 470.0e-12,
        .C_output_couple  = 1.0e-6,
    },

    .sustain = {
        .R_pot         = 100.0e3,
        .R_pot_bottom  = 1.0e3,
        .C_wiper_couple = 100.0e-9,
    },

    .clipper_1 = {
        .q      = { .Is = 1.5e-14, .Bf = 500.0, .Vaf = 53.0 },   // BC550C
        .diodes = { .Is = 2.52e-9, .N  = 1.752,
                    .count_fwd = 1, .count_rev = 1 },            // 1N4148
        .R_input_series   = 10.0e3,
        .R_base_to_gnd    = 100.0e3,
        .R_feedback       = 470.0e3,
        .R_collector_load = 10.0e3,
        .R_emitter        = 150.0,
        .C_feedback_hf    = 470.0e-12,
        .C_diode_couple   = 1.0e-6,
        .C_output_couple  = 100.0e-9,
    },

    .clipper_2 = {
        .q      = { .Is = 1.5e-14, .Bf = 500.0, .Vaf = 53.0 },   // BC550C
        .diodes = { .Is = 2.52e-9, .N  = 1.752,
                    .count_fwd = 1, .count_rev = 1 },            // 1N4148
        .R_input_series   = 10.0e3,
        .R_base_to_gnd    = 100.0e3,
        .R_feedback       = 470.0e3,
        .R_collector_load = 10.0e3,
        .R_emitter        = 150.0,
        .C_feedback_hf    = 470.0e-12,
        .C_diode_couple   = 1.0e-6,
        .C_output_couple  = 0.0,            // tone stack wired directly
    },

    .tone = {
        .R_pot       = 100.0e3,
        .C_hp_series = 4.0e-9,
        .R_hp_shunt  = 27.0e3,              // V9 delta: 22k -> 27k
        .R_lp_series = 27.0e3,              // V9 delta: 39k -> 27k
        .C_lp_shunt  = 10.0e-9,
    },

    .output_booster = {
        .q = { .Is = 1.5e-14, .Bf = 500.0, .Vaf = 53.0 },   // BC550C
        .R_input_series   = 0.0,
        .R_base_to_gnd    = 100.0e3,        // R3
        .R_feedback       = 0.0,
        .R_collector_load = 15.0e3,         // R6
        .R_emitter        = 3.3e3,          // R4
        .R_bias_top       = 430.0e3,        // R7  (inherited from V3; confirm)
        .C_input_couple   = 100.0e-9,       // C3
        .C_feedback_hf    = 0.0,
        .C_output_couple  = 100.0e-9,       // C2
    },

    .volume = {
        .R_pot = 500.0e3,                   // V9 delta: 100k -> 500k log
    },

    .taper_sustain = PotTaper::Log,
    .taper_tone    = PotTaper::Linear,
    .taper_volume  = PotTaper::Log,
};

// -- V1 & V2 variants -------------------------------------------------------
//
// The V1 (earliest) and V2 (early-1970s) families are generated from the
// authoritative OmegaCake.xlsx spreadsheet -- see variants_generated.h, pulled
// in at the bottom of this namespace. Its schematic-exact tabs cover
// V2_Violet_Bighorn, V1_66_5/6/8, The_Pyramid, V1_72_1/2, V1_Super_Omega and
// the V2_47/Classic/Pink/Blue_Bighorn variants. V3, V6 and V9 are hand-
// written here: they have no spreadsheet tab (modern high-gain family).

// -- V6 variant ("reverse logo", ~1979-84, board 3034) ---------------------------
//
// The reference: "one standardised schematic" (consistent unit-to-unit) = the Tone-
// Bypass V3 reworked onto a half-size board. High-confidence data. Reuses the
// V3's exact 2N5088 + 1N914 models.
//
// Deltas vs V3 (see v6-verification-checklist.md):
//   * Tone HP shunt R5: 22k -> 39k  (HP == LP == 39k: the LEAST-scooped board)
//   * All three pots 100k -> 150k, and LINEAR taper (the reference: "Pots 150k lin")
// Everything else (clip caps, 470 pF feedback, R7, diodes) inherited from V3.
inline constexpr CakeVariant V6 = {
    .name = "V6",
    .vcc  = 9.0,

    .input_booster = {
        .q = { .Is = 5.9e-15, .Bf = 540.0, .Vaf = 62.0 },   // 2N5088
        .R_input_series   = 39.0e3,
        .R_base_to_gnd    = 47.0e3,
        .R_feedback       = 470.0e3,
        .R_collector_load = 10.0e3,
        .R_emitter        = 100.0,
        .R_bias_top       = 0.0,
        .C_input_couple   = 1.0e-6,
        .C_feedback_hf    = 470.0e-12,
        .C_output_couple  = 1.0e-6,
    },

    .sustain = {
        .R_pot         = 150.0e3,           // V6 delta: 100k -> 150k
        .R_pot_bottom  = 1.0e3,
        .C_wiper_couple = 100.0e-9,
    },

    .clipper_1 = {
        .q      = { .Is = 5.9e-15, .Bf = 540.0, .Vaf = 62.0 },   // 2N5088
        .diodes = { .Is = 2.52e-9, .N  = 1.752,
                    .count_fwd = 1, .count_rev = 1 },
        .R_input_series   = 10.0e3,
        .R_base_to_gnd    = 100.0e3,
        .R_feedback       = 470.0e3,
        .R_collector_load = 10.0e3,
        .R_emitter        = 150.0,
        .C_feedback_hf    = 470.0e-12,
        .C_diode_couple   = 1.0e-6,
        .C_output_couple  = 100.0e-9,
    },

    .clipper_2 = {
        .q      = { .Is = 5.9e-15, .Bf = 540.0, .Vaf = 62.0 },   // 2N5088
        .diodes = { .Is = 2.52e-9, .N  = 1.752,
                    .count_fwd = 1, .count_rev = 1 },
        .R_input_series   = 10.0e3,
        .R_base_to_gnd    = 100.0e3,
        .R_feedback       = 470.0e3,
        .R_collector_load = 10.0e3,
        .R_emitter        = 150.0,
        .C_feedback_hf    = 470.0e-12,
        .C_diode_couple   = 1.0e-6,
        .C_output_couple  = 0.0,
    },

    .tone = {
        .R_pot       = 150.0e3,             // V6 delta: 100k -> 150k
        .C_hp_series = 4.0e-9,
        .R_hp_shunt  = 39.0e3,              // V6 delta: 22k -> 39k (least scoop)
        .R_lp_series = 39.0e3,
        .C_lp_shunt  = 10.0e-9,
    },

    .output_booster = {
        .q = { .Is = 5.9e-15, .Bf = 540.0, .Vaf = 62.0 },   // 2N5088
        .R_input_series   = 0.0,
        .R_base_to_gnd    = 100.0e3,
        .R_feedback       = 0.0,
        .R_collector_load = 15.0e3,
        .R_emitter        = 3.3e3,
        .R_bias_top       = 430.0e3,
        .C_input_couple   = 100.0e-9,
        .C_feedback_hf    = 0.0,
        .C_output_couple  = 100.0e-9,
    },

    .volume = {
        .R_pot = 150.0e3,                   // V6 delta: 100k -> 150k
    },

    // V6 pots are 150k LINEAR.
    .taper_sustain = PotTaper::Linear,
    .taper_tone    = PotTaper::Linear,
    .taper_volume  = PotTaper::Linear,
};

// Spreadsheet-generated variants (V1/V2 families). Included INSIDE the
// namespace so the literals see CakeVariant / PotTaper defined above.
#include "variants_generated.h"

} // namespace tocake
