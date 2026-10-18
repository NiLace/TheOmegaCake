// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// The Standard circuit (27 slots) as a DK engine: the circuit description for
// src/dk_engine.h. Ten nonlinear ports (four transistors, two diode pairs), 25
// nodes, 13 capacitors.
//
// Same circuit and same device laws as FullPedalV3 -- the node set, the pot legs
// (R * fraction + 1e-6), the trapezoidal companion of every capacitor, SPICE's
// law for these transistor cards (with or without its reverse term, as
// src/standard_law.h selects: the shipped build keeps it, because Vbc > 0 does
// occur on some Standard slots at a hot input) and the antiparallel pair as
// 2*Is*sinh -- so the
// two engines solve the same equations with different algebra. What differs is
// what each does with them per sample: FullPedalV3 takes one Newton step on 22
// node voltages; this engine converges on 10 port voltages.
//
// Measured against FullPedalV3 FORCED to converge:
// every Standard slot, four knob points, two levels, -164.7 dB worst -- the
// same equations to the last bits a double carries through this network.

#pragma once

#include "dk_engine.h"
#include "standard_law.h"

namespace tocake {

struct StandardCircuit {
    // Node indices of the linear network (ground, the 9 V rail and the input
    // jack are known voltages and not unknowns).
    enum Node : int {
        IN_C1, Q4B, Q4C, Q4E, SUS_TOP, SUS_W, SUS_BOT, C5_R19,
        Q3B, Q3C, Q3E, Q3_AUX, C13_R12,
        Q2B, Q2C, Q2E, Q2_AUX,
        TONE_TOP, TONE_BOT, TONE_W,
        Q1B, Q1C, Q1E, VOL_TOP, OUT,
        NODES
    };
    static constexpr int NN = NODES, NB = 4, ND = 2, S = 13, NT = 3;
    static constexpr bool FULL_EM = true;

    static constexpr std::array<DkBjtPins, NB> kBjts = {{
        {Q4C, Q4B, Q4E}, {Q3C, Q3B, Q3E}, {Q2C, Q2B, Q2E}, {Q1C, Q1B, Q1E},
    }};
    static constexpr std::array<DkPairPins, ND> kPairs = {{
        {Q3_AUX, Q3C},              // clipper-1 pair (D3/D4)
        {Q2_AUX, Q2C},              // clipper-2 pair (D1/D2)
    }};
    static constexpr std::array<DkCapPins, S> kCaps = {{
        {IN_C1, Q4B},       // C1
        {Q4B, Q4C},         // C10
        {Q4C, SUS_TOP},     // C4
        {SUS_W, C5_R19},    // C5
        {Q3B, Q3C},         // C12
        {Q3B, Q3_AUX},      // C6
        {Q3C, C13_R12},     // C13
        {Q2B, Q2C},         // C11
        {Q2B, Q2_AUX},      // C7
        {Q2C, TONE_TOP},    // C9
        {TONE_BOT, kDkGnd}, // C8
        {TONE_W, Q1B},      // C3
        {Q1C, VOL_TOP},     // C2
    }};
    // Out: the volume wiper. C1Out: the node after C13 (as the deck's tap).
    // C2Out: Q2's collector about its rest (the deck subtracts its bias).
    static constexpr std::array<DkTap, NT> kTaps = {{ {OUT, false}, {C13_R12, false}, {Q2C, true} }};
    static constexpr int kSeedHigh[] = { Q4B, Q3B, Q2B, Q1B };

    static std::array<double, S> cap_values(const CakeVariant& v) noexcept {
        const auto& ib = v.input_booster; const auto& sus = v.sustain;
        const auto& c1 = v.clipper_1;     const auto& c2 = v.clipper_2;
        const auto& t  = v.tone;          const auto& ob = v.output_booster;
        return { ib.C_input_couple, ib.C_feedback_hf, ib.C_output_couple, sus.C_wiper_couple,
                 c1.C_feedback_hf, c1.C_diode_couple, c1.C_output_couple,
                 c2.C_feedback_hf, c2.C_diode_couple,
                 t.C_hp_series, t.C_lp_shunt, ob.C_input_couple, ob.C_output_couple };
    }

    template <class R>
    static void resistors(const CakeVariant& v, const DkKnobs& k, R&& r) noexcept {
        const auto& ib = v.input_booster; const auto& sus = v.sustain;
        const auto& c1 = v.clipper_1;     const auto& c2 = v.clipper_2;
        const auto& t  = v.tone;          const auto& ob = v.output_booster;
        const auto& vol = v.volume;
        r(kDkIn, IN_C1, ib.R_input_series);                     // R2
        r(Q4B, kDkGnd, ib.R_base_to_gnd);                       // R14
        r(Q4B, Q4C, ib.R_feedback);                             // R9
        r(kDkVcc, Q4C, ib.R_collector_load);                    // R13
        r(Q4E, kDkGnd, ib.R_emitter);                           // R22
        r(SUS_TOP, SUS_W, sus.R_pot * (1.0 - k.sus) + 1e-6);    // sustain, top leg
        r(SUS_W, SUS_BOT, sus.R_pot * k.sus + 1e-6);            // sustain, bottom leg
        r(SUS_BOT, kDkGnd, sus.R_pot_bottom);                   // R23
        r(C5_R19, Q3B, c1.R_input_series);                      // R19
        r(Q3B, kDkGnd, c1.R_base_to_gnd);                       // R20
        r(Q3B, Q3C, c1.R_feedback);                             // R17
        r(kDkVcc, Q3C, c1.R_collector_load);                    // R18
        r(Q3E, kDkGnd, c1.R_emitter);                           // R21
        r(C13_R12, Q2B, c2.R_input_series);                     // R12
        r(Q2B, kDkGnd, c2.R_base_to_gnd);                       // R16
        r(Q2B, Q2C, c2.R_feedback);                             // R15
        r(kDkVcc, Q2C, c2.R_collector_load);                    // R11
        r(Q2E, kDkGnd, c2.R_emitter);                           // R10
        r(TONE_TOP, kDkGnd, t.R_hp_shunt);                      // R5
        r(Q2C, TONE_BOT, t.R_lp_series);                        // R8
        r(TONE_TOP, TONE_W, t.R_pot * (1.0 - k.tone) + 1e-6);   // tone, top leg
        r(TONE_W, TONE_BOT, t.R_pot * k.tone + 1e-6);           // tone, bottom leg
        r(kDkVcc, Q1B, ob.R_bias_top);                          // R7
        r(Q1B, kDkGnd, ob.R_base_to_gnd);                       // R3
        r(kDkVcc, Q1C, ob.R_collector_load);                    // R6
        r(Q1E, kDkGnd, ob.R_emitter);                           // R4
        r(VOL_TOP, OUT, vol.R_pot * (1.0 - k.vol) + 1e-6);      // volume, top leg
        r(OUT, kDkGnd, vol.R_pot * k.vol + 1e-6);               // volume, bottom leg
    }

    static bjt_params bjt(const CakeVariant& v, int k) noexcept {
        const bjt_params* q[4] = { &v.input_booster.q, &v.clipper_1.q, &v.clipper_2.q, &v.output_booster.q };
        return *q[k];
    }
    static diode_params pair(const CakeVariant& v, int k) noexcept {
        return k == 0 ? v.clipper_1.diodes : v.clipper_2.diodes;
    }
};

using DkStandard = DkEngine<StandardCircuit>;

}  // namespace tocake
