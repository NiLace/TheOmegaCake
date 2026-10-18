// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// The op-amp version of the circuit ("The Op Cake") -- component values.
//
// SOURCE OF TRUTH for this variant is the KiCad netlist at
// netlists/TheOpCake.net (NOT the OmegaCake.xlsx spreadsheet -- the op-amp
// circuit is a different topology from the transistor one and has no
// spreadsheet tab). Every value below is transcribed directly from that
// netlist; the matching SPICE reference deck is spice/the_op_cake.cir.
//
// This is the late-1970s "IC" / op-amp version of the circuit: two
// JRC4558 op-amp gain stages (U1A inverting booster, U1B non-inverting
// Sustain stage) feeding a single LM741 (U2) clipping stage whose feedback
// loop holds the diode clamp -- THREE 1N4448 in series each direction,
// anti-parallel (6 diodes total), confirmed against the original
// schematic. A passive mid-scoop tone stack + a tone-bypass switch + volume
// finish the chain. Single 9 V supply; U1 biased at +4.5 V, U2 at +4.95 V.
//
// Unlike the transistor model (CakeVariant), the nonlinearity here is
// localised to ONE node: U2's summing junction. The op-amps idealise to
// infinite-gain (their outputs are low-impedance and isolate the blocks),
// so the DSP solves the diode clamp implicitly at that node and treats
// everything else as linear filters. See op_cake.h.

#pragma once

namespace tocake {

// Shockley diode parameters + per-direction series count. The op-amp version's
// clamp is 3 series diodes each way (count_fwd = count_rev = 3), which puts
// the clip threshold near ~2 V -- the high, soft "open" voice of the op-amp version.
// The engine (OpCake) reads count_fwd for BOTH directions: the two must be
// equal, and every model below is checked for it.
struct op_diode_params {
    double Is;          // reverse saturation current (A)
    double N;           // emission coefficient (unitless)
    int    count_fwd;   // diodes in series, forward antiparallel branch
    int    count_rev;   // diodes in series, reverse antiparallel branch
};

// Op-amp idealisation knobs.
//   * The fast 4558 stages (U1A/U1B) solve as ideal infinite-gain with a rail
//     clamp -- only open_gain / v_bias matter (gbw_hz = comp_cap = 0).
//   * The LM741 (U2) uses a DYNAMIC macromodel: a transconductance input stage
//     whose output current is LIMITED (= the slew mechanism) charges the
//     compensation cap (= dominant pole / finite GBW), buffered to the output.
//     Set open_gain (DC Aol), gbw_hz, slew_v_per_us, comp_cap; the derived
//     gm/Imax/Rp are computed in OpCake::prepare and MUST match the SPICE
//     deck's .param block in spice/the_op_cake.cir.
struct op_amp_params {
    double v_bias;          // non-inverting (+) input reference (V)
    double open_gain;       // open-loop DC gain Aol (V/V)
    double slew_v_per_us;   // output slew-rate limit (0 = unlimited / ideal)
    double gbw_hz;          // gain-bandwidth product (0 = ideal, no pole)
    double comp_cap;        // compensation cap (F) -- the macromodel scale
};

struct OpCakeParams {
    const char* name;
    double vcc;                 // supply rail (V)

    // -- Input network + U1A inverting booster ----------------------------
    double C1;                  // input DC block (.15 uF)
    double R1;                  // input pulldown (1 M) -- ignored for a
                                //   voltage-source drive; carried for record
    double R2;                  // inverting input series (56 k)
    double R3;                  // U1A feedback (330 k)  -> gain ~= -R3/R2
    op_amp_params u1a;          // +4.5 V, fast (4558)

    // -- Interstage shaping + U1B non-inverting Sustain stage -------------
    double R4;                  // U1A out -> interstage node X (10 k)
    double C3;                  // X <-> U1B out, 4700 pF (feedback HF shape)
    double R5;                  // X -> U1B(+) (47 k)
    double C2;                  // U1B(+) -> GND, .01 uF (LP into + input)
    double R6;                  // U1B feedback (560 k)
    double R7;                  // U1B(-) -> node Z (62 k)
    double C4;                  // Z -> GND, 10 uF (AC-grounds the gain leg)
    double R8;                  // Z -> Sustain bottom (47 R)
    double R_sustain;           // Sustain pot (10 k)
    double C5;                  // Sustain wiper coupling (4.7 uF)
    op_amp_params u1b;          // +4.5 V, fast (4558)

    // -- U2 LM741 clipper (the nonlinear heart) ---------------------------
    double R9;                  // U2 inverting input series (8.2 k)
    double R10;                 // U2 feedback (470 k) -> max gain ~ R10/R9
    double C6;                  // U2 feedback HF cap (150 pF) -- STATEFUL,
                                //   makes the clamp a dynamic nonlinearity
    op_diode_params diodes;     // 3x 1N4448 series, each direction
    op_amp_params u2;           // +4.95 V; LM741 (with slew)

    // -- Passive tone stack (the mid scoop) ------------------------------
    double C8;                  // treble leg series cap (.1 uF) -> tone top
    double R12;                 // treble leg shunt to GND (1.2 k)
    double C7;                  // bass leg series cap (1 uF)
    double R11;                 // bass leg series R (5.6 k)
    double C9;                  // bass leg shunt cap to GND (.12 uF)
    double R_tone;              // Tone pot (10 k)

    // -- Tone-bypass switch network (R19/R20/C13) -------------------------
    // SPDT picks the tone-stack wiper (engaged) or a flat post-clip tap.
    double C13;                 // .15 uF, flat-path series cap from U2 out
    double R19;                 // 47 k
    double R20;                 // 47 k

    // -- Output volume ----------------------------------------------------
    double R_volume;            // Volume pot (100 k)
    double R13;                 // output pulldown (100 k)

    // -- Feedback mod of the 1322 board (the "Guru Cake") -----------------
    // false = the stock op-amp circuit. true = the recommended mod: the sustain-pot
    // bottom (R8) is grounded DIRECTLY (J2 cut, J3 placed) instead of being
    // tied to the U1B feedback node Z -- this unloads the IC1b gain leg from
    // the sustain pot and gives the "slight crunch". The mod's C6 (10uF) into
    // the sustain pot is a sub-2 Hz DC block; with U1B modelled as an ideal
    // (load-independent) op-amp it is sonically transparent, so it is omitted.
    bool guru_mod = false;
};

// -- "The Op Cake" -- transcribed from netlists/TheOpCake.net ---------------
//
// Op-amp models (must stay in lock-step with the SPICE deck so the null is
// honest):
//   * U1A, U1B (JRC4558 halves): idealised as infinite-gain in the DSP, very-
//     high-finite-gain (1e6) in the deck. The exact figure is immaterial in
//     deep feedback; what matters is the supply-rail clipping, modelled on
//     both sides. No slew model (the 4558 is fast and never the limiter here).
//   * U2 (LM741): a full dynamic macromodel -- finite open-loop gain (Aol),
//     gain-bandwidth (gbw_hz) and slew rate (slew_v_per_us), realised via a
//     compensation cap (comp_cap). The values below are real-741 specs; the
//     derived gm/Imax/Rp in OpCake::prepare must equal the deck's .param block.
inline constexpr OpCakeParams TheOpCake = {
    .name = "OpCake",
    .vcc  = 9.0,

    .C1 = 0.15e-6,  .R1 = 1.0e6,   .R2 = 56.0e3,  .R3 = 330.0e3,
    .u1a = { .v_bias = 4.5,  .open_gain = 1.0e6, .slew_v_per_us = 0.0,
             .gbw_hz = 0.0, .comp_cap = 0.0 },   // ideal (4558, no slew model)

    .R4 = 10.0e3,   .C3 = 4700.0e-12, .R5 = 47.0e3, .C2 = 0.01e-6,
    .R6 = 560.0e3,  .R7 = 62.0e3,     .C4 = 10.0e-6, .R8 = 47.0,   // schematic R7 = 47R (47k would leave the Sustain knob inert)
    .R_sustain = 10.0e3, .C5 = 4.7e-6,
    .u1b = { .v_bias = 4.5,  .open_gain = 1.0e6, .slew_v_per_us = 0.0,
             .gbw_hz = 0.0, .comp_cap = 0.0 },   // ideal (4558, no slew model)

    .R9 = 8.2e3,    .R10 = 470.0e3,   .C6 = 150.0e-12,
    .diodes = { .Is = 2.52e-9, .N = 1.752, .count_fwd = 3, .count_rev = 3 },
    // LM741 macromodel: Aol = 2e5, GBW = 1 MHz, SR = 0.5 V/us, Cc = 30 pF
    // (real-741 specs). Derived gm/Imax/Rp in OpCake::prepare must equal the
    // .param block in spice/the_op_cake.cir.
    .u2 = { .v_bias = 4.945, .open_gain = 2.0e5, .slew_v_per_us = 0.5,
            .gbw_hz = 1.0e6, .comp_cap = 30.0e-12 },

    .C8 = 0.1e-6,   .R12 = 1.2e3,  .C7 = 1.0e-6,
    .R11 = 5.6e3,   .C9 = 0.12e-6, .R_tone = 10.0e3,

    .C13 = 0.15e-6, .R19 = 47.0e3, .R20 = 47.0e3,  // C13 150n per the 1322 board's schematic

    .R_volume = 100.0e3, .R13 = 100.0e3,
};

// -- "The Little Pumpkin" -- the small one-knob op-amp version (V4-5 era) -----
// Same op-amp core as TheOpCake (IC1 4558, IC2 741, 1N4148 3+3 clamp)
// with the small version's values from a published trace of it
// (2001): input cap .068uF (vs .15u), input series 68k (vs 56k), the sustain
// wiper coupling 1uF (vs 4.7u) and a 25k volume pot (vs 100k). The small
// version has NO Sustain or Tone pots -- the sustain pot is a fixed 10k with the
// output taken at its top (fuzz at max) and Tone is a 2-position switch at the
// tone-pot ends. Those are applied at the lv2 layer (sustain pinned, tone
// driven by the Tone Switch port, tone_switch); only Volume is a live knob.
// The 25k is not among the parts the trace marks as differing from the
// earlier drawing it builds on, and the 1322 board's sheet has 100k: it is the
// trace's reading, kept.
inline constexpr OpCakeParams TheLittlePumpkin = [] {
    OpCakeParams p = TheOpCake;
    p.C1 = 0.068e-6;   // input DC-block (TheOpCake .15uF)
    p.R2 = 68.0e3;     // inverting input series (TheOpCake 56k)
    p.C5 = 1.0e-6;     // sustain-wiper coupling (TheOpCake 4.7uF)
    p.R_volume = 25.0e3;  // volume pot (TheOpCake 100k)
    return p;
}();

// -- "The Guru Cake" -- the op-amp circuit with the 1322 board's feedback mod --
// Identical to TheOpCake but with the recommended negative-feedback mod
// (guru_mod = true): the sustain-pot bottom grounded directly, unloading the
// IC1b gain leg for a slight crunch. All three pots are live (like TheOpCake).
inline constexpr OpCakeParams TheGuruCake = [] {
    OpCakeParams p = TheOpCake;
    p.guru_mod = true;
    return p;
}();

// OpCake models a symmetric clamp (count_fwd for both directions).
static_assert(TheOpCake.diodes.count_fwd == TheOpCake.diodes.count_rev
              && TheLittlePumpkin.diodes.count_fwd == TheLittlePumpkin.diodes.count_rev
              && TheGuruCake.diodes.count_fwd == TheGuruCake.diodes.count_rev,
              "OpCake models a symmetric clamp: count_fwd must equal count_rev");

} // namespace tocake
