// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// Circuit descriptions of the bespoke families for the generic DK engine
// (src/dk_engine.h). Each one is transcribed from the equations of its node
// engine -- the same nodes, branches and device laws -- and is
// checked against that engine FORCED to converge, then against its SPICE deck.
//
// The pot legs keep the node engines' convention, R * fraction + 1e-6.
#pragma once

#include "dk_engine.h"
#include "dk_solver_generated.h"

namespace tocake {

// Component values shared by every four-stage circuit that keeps the Standard
// set of thirteen capacitors in the Standard order.
struct StandardCircuitValues {
    static std::array<double, 13> caps(const CakeVariant& v) noexcept {
        const auto& ib = v.input_booster; const auto& sus = v.sustain;
        const auto& c1 = v.clipper_1;     const auto& c2 = v.clipper_2;
        const auto& t  = v.tone;          const auto& ob = v.output_booster;
        return { ib.C_input_couple, ib.C_feedback_hf, ib.C_output_couple, sus.C_wiper_couple,
                 c1.C_feedback_hf, c1.C_diode_couple, c1.C_output_couple,
                 c2.C_feedback_hf, c2.C_diode_couple,
                 t.C_hp_series, t.C_lp_shunt, ob.C_input_couple, ob.C_output_couple };
    }
    static bjt_params bjt(const CakeVariant& v, int k) noexcept {
        const bjt_params* q[4] = { &v.input_booster.q, &v.clipper_1.q, &v.clipper_2.q, &v.output_booster.q };
        return *q[k];
    }
    static diode_params pair(const CakeVariant& v, int k) noexcept {
        return k == 0 ? v.clipper_1.diodes : v.clipper_2.diodes;
    }
};

// ---- V2 Distortion & Cake, V2 Coloured Goat --------------------------------
//
// The Standard circuit with C14 AC-coupling Q2's collector to the tone stack
// (a node TIN feeding both C9 and R8), every transistor on the Gummel-Poon law
// with the reverse term. Distortion & Cake keeps both diode pairs; the Coloured
// Goat's clipper 1 is a clean gain stage (no pair, no C6).
template <bool Clip1Diodes>
struct C14CakeCircuit {
    enum Node : int {
        IN_C1, Q4B, Q4C, Q4E, SUS_TOP, SUS_W, SUS_BOT, C5_R19,
        Q3B, Q3C, Q3E, C13_R12,
        Q2B, Q2C, Q2E, Q2_AUX, TIN,
        TONE_TOP, TONE_BOT, TONE_W,
        Q1B, Q1C, Q1E, VOL_TOP, OUT,
        Q3_AUX,                     // only with clipper-1 diodes
        NODES_WITH_AUX
    };
    static constexpr int NN = Clip1Diodes ? NODES_WITH_AUX : Q3_AUX;
    static constexpr int NB = 4, ND = Clip1Diodes ? 2 : 1, S = Clip1Diodes ? 14 : 13, NT = 3;
    static constexpr bool FULL_EM = true;

    static constexpr std::array<DkBjtPins, NB> kBjts = {{
        {Q4C, Q4B, Q4E}, {Q3C, Q3B, Q3E}, {Q2C, Q2B, Q2E}, {Q1C, Q1B, Q1E},
    }};
    static constexpr std::array<DkPairPins, ND> make_pairs() noexcept {
        std::array<DkPairPins, ND> p{};
        if constexpr (Clip1Diodes) { p[0] = {Q3_AUX, Q3C}; p[1] = {Q2_AUX, Q2C}; }
        else                       { p[0] = {Q2_AUX, Q2C}; }
        return p;
    }
    static constexpr std::array<DkPairPins, ND> kPairs = make_pairs();
    static constexpr std::array<DkCapPins, S> make_caps() noexcept {
        std::array<DkCapPins, S> c{};
        int n = 0;
        c[n++] = {IN_C1, Q4B};      // C1
        c[n++] = {Q4B, Q4C};        // C10
        c[n++] = {Q4C, SUS_TOP};    // C4
        c[n++] = {SUS_W, C5_R19};   // C5
        c[n++] = {Q3B, Q3C};        // C12
        if constexpr (Clip1Diodes) c[n++] = {Q3B, Q3_AUX};   // C6
        c[n++] = {Q3C, C13_R12};    // C13
        c[n++] = {Q2B, Q2C};        // C11
        c[n++] = {Q2B, Q2_AUX};     // C7
        c[n++] = {Q2C, TIN};        // C14
        c[n++] = {TIN, TONE_TOP};   // C9
        c[n++] = {TONE_BOT, kDkGnd};// C8
        c[n++] = {TONE_W, Q1B};     // C3
        c[n++] = {Q1C, VOL_TOP};    // C2
        return c;
    }
    static constexpr std::array<DkCapPins, S> kCaps = make_caps();
    static constexpr std::array<DkTap, NT> kTaps = {{ {OUT, false}, {C13_R12, false}, {Q2C, true} }};
    static constexpr int kSeedHigh[] = { Q4B, Q3B, Q2B, Q1B };

    static std::array<double, S> cap_values(const CakeVariant& v) noexcept {
        const auto& ib = v.input_booster; const auto& sus = v.sustain;
        const auto& c1 = v.clipper_1;     const auto& c2 = v.clipper_2;
        const auto& t  = v.tone;          const auto& ob = v.output_booster;
        std::array<double, S> c{};
        int n = 0;
        c[n++] = ib.C_input_couple; c[n++] = ib.C_feedback_hf; c[n++] = ib.C_output_couple;
        c[n++] = sus.C_wiper_couple; c[n++] = c1.C_feedback_hf;
        if constexpr (Clip1Diodes) c[n++] = c1.C_diode_couple;
        c[n++] = c1.C_output_couple; c[n++] = c2.C_feedback_hf; c[n++] = c2.C_diode_couple;
        c[n++] = c2.C_output_couple;                    // C14
        c[n++] = t.C_hp_series; c[n++] = t.C_lp_shunt;
        c[n++] = ob.C_input_couple; c[n++] = ob.C_output_couple;
        return c;
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
        r(TIN, TONE_BOT, t.R_lp_series);                        // R8, from the C14 node
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
        if constexpr (Clip1Diodes) return k == 0 ? v.clipper_1.diodes : v.clipper_2.diodes;
        else return v.clipper_2.diodes;
    }
};

using DkDistortionCake = DkEngine<C14CakeCircuit<true>>;
using DkColouredGoat   = DkEngine<C14CakeCircuit<false>>;

// ---- Vm Precision Cake, Vm German Cake (simplified tone stack) --------------
//
// The Distortion & Cake core with the treble cap C9 and the bass resistor R8
// removed: C14 couples Q2's collector straight onto the tone-pot top (TIN, with
// a shunt R to ground), the pot bottom carries C8 to ground, the wiper feeds
// C3. The output booster's emitter has an AC bypass leg (C15 then R28 to
// ground) across R4, and the volume top a load R29 to ground. Per-stage diodes
// (the German Cake runs germanium on clipper 2).
struct SimpleToneCircuit {
    enum Node : int {
        IN_C1, Q4B, Q4C, Q4E, SUS_TOP, SUS_W, SUS_BOT, C5_R19,
        Q3B, Q3C, Q3E, Q3_AUX, C13_R12,
        Q2B, Q2C, Q2E, Q2_AUX, TIN,
        TONE_BOT, TONE_W,
        Q1B, Q1C, Q1E, E_BYP, VOL_TOP, OUT,
        NODES
    };
    static constexpr int NN = NODES, NB = 4, ND = 2, S = 14, NT = 3;
    static constexpr bool FULL_EM = true;

    static constexpr std::array<DkBjtPins, NB> kBjts = {{
        {Q4C, Q4B, Q4E}, {Q3C, Q3B, Q3E}, {Q2C, Q2B, Q2E}, {Q1C, Q1B, Q1E},
    }};
    static constexpr std::array<DkPairPins, ND> kPairs = {{ {Q3_AUX, Q3C}, {Q2_AUX, Q2C} }};
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
        {Q2C, TIN},         // C14 (the ~4.7 uF clip-2 coupling)
        {TONE_BOT, kDkGnd}, // C8
        {TONE_W, Q1B},      // C3
        {Q1C, VOL_TOP},     // C2
        {Q1E, E_BYP},       // C15, emitter bypass
    }};
    static constexpr std::array<DkTap, NT> kTaps = {{ {OUT, false}, {C13_R12, false}, {Q2C, true} }};
    static constexpr int kSeedHigh[] = { Q4B, Q3B, Q2B, Q1B };

    static std::array<double, S> cap_values(const CakeVariant& v) noexcept {
        const auto& ib = v.input_booster; const auto& sus = v.sustain;
        const auto& c1 = v.clipper_1;     const auto& c2 = v.clipper_2;
        const auto& t  = v.tone;          const auto& ob = v.output_booster;
        return { ib.C_input_couple, ib.C_feedback_hf, ib.C_output_couple, sus.C_wiper_couple,
                 c1.C_feedback_hf, c1.C_diode_couple, c1.C_output_couple,
                 c2.C_feedback_hf, c2.C_diode_couple, c2.C_output_couple,
                 t.C_lp_shunt, ob.C_input_couple, ob.C_output_couple, ob.C_emitter_bypass };
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
        r(TIN, kDkGnd, t.R_hp_shunt);                           // pot-top shunt (R8 / R27)
        r(TIN, TONE_W, t.R_pot * (1.0 - k.tone) + 1e-6);        // tone, top leg
        r(TONE_W, TONE_BOT, t.R_pot * k.tone + 1e-6);           // tone, bottom leg
        r(kDkVcc, Q1B, ob.R_bias_top);                          // R7
        r(Q1B, kDkGnd, ob.R_base_to_gnd);                       // R3
        r(kDkVcc, Q1C, ob.R_collector_load);                    // R6
        r(Q1E, kDkGnd, ob.R_emitter);                           // R4
        r(E_BYP, kDkGnd, ob.R_emitter_bypass);                  // R28, after C15
        r(VOL_TOP, kDkGnd, ob.R_output_load);                   // R29
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

using DkSimpleToneCake = DkEngine<SimpleToneCircuit>;

// ---- Vm Dreamy Cake ----------------------------------------------------------
//
// The Standard (direct-coupled) circuit with the emitter resistors of the first
// three stages removed: Q4, Q3 and Q2 have their emitters on ground, so those
// emitters are not nodes. Gummel-Poon with the reverse term (the grounded,
// maximum-gain stages saturate). The tone pot is wired in reverse on the real
// pedal; that is a knob direction and lives in lv2.cpp, not here.
//
// The supply network: the four stages draw their current from the 9 V through
// R27 (100 R) into a rail held by C14 (330 uF), so the rail rests below 9 V and
// moves with the signal (tau = R27 * C14 = 33 ms). With `Rail` the rail is a
// node and kDkVcc is the battery; without it the rail is the variant's static
// `vcc` (its rest value), which is what the node engine arbiter holds.

template <bool Rail>
struct DreamyCircuitT {
    enum Node : int {
        IN_C1, Q4B, Q4C, SUS_TOP, SUS_W, SUS_BOT, C5_R19,
        Q3B, Q3C, Q3_AUX, C13_R12,
        Q2B, Q2C, Q2_AUX,
        TONE_TOP, TONE_BOT, TONE_W,
        Q1B, Q1C, Q1E, VOL_TOP, OUT,
        RAIL,                       // only with the supply network
        NODES_WITH_RAIL
    };
    static constexpr int NN = Rail ? NODES_WITH_RAIL : RAIL, NB = 4, ND = 2, S = Rail ? 14 : 13, NT = 3;
    static constexpr bool FULL_EM = true;
    static constexpr int V = Rail ? RAIL : kDkVcc;       // where the collector loads return
    static constexpr double kR27 = 100.0, kC14 = 330e-6;

    static constexpr std::array<DkBjtPins, NB> kBjts = {{
        {Q4C, Q4B, kDkGnd}, {Q3C, Q3B, kDkGnd}, {Q2C, Q2B, kDkGnd}, {Q1C, Q1B, Q1E},
    }};
    static constexpr std::array<DkPairPins, ND> kPairs = {{ {Q3_AUX, Q3C}, {Q2_AUX, Q2C} }};
    static constexpr std::array<DkCapPins, S> make_caps() noexcept {
        std::array<DkCapPins, S> c{{
            {IN_C1, Q4B}, {Q4B, Q4C}, {Q4C, SUS_TOP}, {SUS_W, C5_R19},       // C1 C10 C4 C5
            {Q3B, Q3C}, {Q3B, Q3_AUX}, {Q3C, C13_R12},                        // C12 C6 C13
            {Q2B, Q2C}, {Q2B, Q2_AUX}, {Q2C, TONE_TOP}, {TONE_BOT, kDkGnd},   // C11 C7 C9 C8
            {TONE_W, Q1B}, {Q1C, VOL_TOP},                                    // C3 C2
        }};
        if constexpr (Rail) c[13] = {RAIL, kDkGnd};                           // C14
        return c;
    }
    static constexpr std::array<DkCapPins, S> kCaps = make_caps();
    static constexpr std::array<DkTap, NT> kTaps = {{ {OUT, false}, {C13_R12, false}, {Q2C, true} }};
    static constexpr int kSeedHigh[] = { Q4B, Q3B, Q2B, Q1B };

    static std::array<double, S> cap_values(const CakeVariant& v) noexcept {
        const auto c = StandardCircuitValues::caps(v);
        std::array<double, S> out{};
        for (int k = 0; k < 13; ++k) out[k] = c[k];
        if constexpr (Rail) out[13] = kC14;
        return out;
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
        r(V, Q4C, ib.R_collector_load);                         // R13
        r(SUS_TOP, SUS_W, sus.R_pot * (1.0 - k.sus) + 1e-6);    // sustain, top leg
        r(SUS_W, SUS_BOT, sus.R_pot * k.sus + 1e-6);            // sustain, bottom leg
        r(SUS_BOT, kDkGnd, sus.R_pot_bottom);                   // R23
        r(C5_R19, Q3B, c1.R_input_series);                      // R19
        r(Q3B, kDkGnd, c1.R_base_to_gnd);                       // R20
        r(Q3B, Q3C, c1.R_feedback);                             // R17
        r(V, Q3C, c1.R_collector_load);                         // R18
        r(C13_R12, Q2B, c2.R_input_series);                     // R12
        r(Q2B, kDkGnd, c2.R_base_to_gnd);                       // R16
        r(Q2B, Q2C, c2.R_feedback);                             // R15
        r(V, Q2C, c2.R_collector_load);                         // R11
        r(TONE_TOP, kDkGnd, t.R_hp_shunt);                      // R5
        r(Q2C, TONE_BOT, t.R_lp_series);                        // R8
        r(TONE_TOP, TONE_W, t.R_pot * (1.0 - k.tone) + 1e-6);   // tone, top leg
        r(TONE_W, TONE_BOT, t.R_pot * k.tone + 1e-6);           // tone, bottom leg
        r(V, Q1B, ob.R_bias_top);                               // R7
        r(Q1B, kDkGnd, ob.R_base_to_gnd);                       // R3
        r(V, Q1C, ob.R_collector_load);                         // R6
        r(Q1E, kDkGnd, ob.R_emitter);                           // R4
        r(VOL_TOP, OUT, vol.R_pot * (1.0 - k.vol) + 1e-6);      // volume, top leg
        r(OUT, kDkGnd, vol.R_pot * k.vol + 1e-6);               // volume, bottom leg
        if constexpr (Rail) r(kDkVcc, RAIL, kR27);             // R27
    }

    static bjt_params bjt(const CakeVariant& v, int k) noexcept { return StandardCircuitValues::bjt(v, k); }
    static diode_params pair(const CakeVariant& v, int k) noexcept { return StandardCircuitValues::pair(v, k); }
};

// The battery, read only when the rail is a node (dk_circuit_supply).
struct DreamyCircuitRail : DreamyCircuitT<true> { static constexpr double kSupply = 9.0; };
using DreamyCircuitStatic = DreamyCircuitT<false>;
using DreamyCircuit = DreamyCircuitRail;

using DkDreamyCake = DkEngine<DreamyCircuit>;
using DkDreamyCakeStaticRail = DkEngine<DreamyCircuitStatic>;

// ---- V2 Jumbo Open Goat ---------------------------------------------------------
//
// The Jumbo Open Goat main path: clipper emitters on ground (R10 = R21 = 0), the
// series R12 BEFORE C13 on the way to Q2's base, C14 coupling Q2's collector to
// the tone input, and a 470 pF Miller cap C17 across Q1's base and collector.
// Both clipper pairs take clipper 1's diode card (the node engine does the same).
struct JumboOpenCircuit {
    enum Node : int {
        Q4B, Q4C, Q4E, N_C1,
        SUS_TOP, SUS_W, SUS_BOT, C5N,
        Q3B, Q3C, Q3AUX, C13N,
        Q2B, Q2C, Q2AUX,
        TONE_IN, TT, TB, TW,
        Q1B, Q1C, Q1E,
        VOL_TOP, OUT,
        NODES
    };
    static constexpr int NN = NODES, NB = 4, ND = 2, S = 15, NT = 3;
    static constexpr bool FULL_EM = true;
    static constexpr double kC17 = 470e-12;   // Q1 collector-base Miller cap (schematic C17)

    static constexpr std::array<DkBjtPins, NB> kBjts = {{
        {Q4C, Q4B, Q4E}, {Q3C, Q3B, kDkGnd}, {Q2C, Q2B, kDkGnd}, {Q1C, Q1B, Q1E},
    }};
    static constexpr std::array<DkPairPins, ND> kPairs = {{ {Q3AUX, Q3C}, {Q2AUX, Q2C} }};
    static constexpr std::array<DkCapPins, S> kCaps = {{
        {N_C1, Q4B}, {Q4B, Q4C}, {Q4C, SUS_TOP}, {SUS_W, C5N},          // C1 C10 C4 C5
        {Q3B, Q3C}, {Q3B, Q3AUX}, {C13N, Q2B},                           // C12 C6 C13
        {Q2B, Q2C}, {Q2B, Q2AUX}, {Q2C, TONE_IN}, {TONE_IN, TT},         // C11 C7 C14 C9
        {TB, kDkGnd}, {TW, Q1B}, {Q1C, VOL_TOP}, {Q1B, Q1C},             // C8 C3 C2 C17
    }};
    // Out: the volume wiper. C1Out: the node between R12 and C13. C2Out: Q2's
    // collector about its rest.
    static constexpr std::array<DkTap, NT> kTaps = {{ {OUT, false}, {C13N, false}, {Q2C, true} }};
    static constexpr int kSeedHigh[] = { Q4B, Q3B, Q2B, Q1B };

    static std::array<double, S> cap_values(const CakeVariant& v) noexcept {
        const auto& ib = v.input_booster; const auto& sus = v.sustain;
        const auto& c1 = v.clipper_1;     const auto& c2 = v.clipper_2;
        const auto& t  = v.tone;          const auto& ob = v.output_booster;
        return { ib.C_input_couple, ib.C_feedback_hf, ib.C_output_couple, sus.C_wiper_couple,
                 c1.C_feedback_hf, c1.C_diode_couple, c1.C_output_couple,
                 c2.C_feedback_hf, c2.C_diode_couple, c2.C_output_couple, t.C_hp_series,
                 t.C_lp_shunt, ob.C_input_couple, ob.C_output_couple, kC17 };
    }

    template <class R>
    static void resistors(const CakeVariant& v, const DkKnobs& k, R&& r) noexcept {
        const auto& ib = v.input_booster; const auto& sus = v.sustain;
        const auto& c1 = v.clipper_1;     const auto& c2 = v.clipper_2;
        const auto& t  = v.tone;          const auto& ob = v.output_booster;
        const auto& vol = v.volume;
        r(kDkIn, N_C1, ib.R_input_series);                      // R2
        r(Q4B, kDkGnd, ib.R_base_to_gnd);                       // R14
        r(Q4B, Q4C, ib.R_feedback);                             // R9
        r(kDkVcc, Q4C, ib.R_collector_load);                    // R13
        r(Q4E, kDkGnd, ib.R_emitter);                           // R22
        r(C5N, Q3B, c1.R_input_series);                         // R19
        r(SUS_BOT, kDkGnd, sus.R_pot_bottom);                   // R23
        r(Q3B, kDkGnd, c1.R_base_to_gnd);                       // R20
        r(Q3B, Q3C, c1.R_feedback);                             // R17
        r(kDkVcc, Q3C, c1.R_collector_load);                    // R18
        r(C13N, Q3C, c2.R_input_series);                        // R12
        r(Q2B, kDkGnd, c2.R_base_to_gnd);                       // R16
        r(Q2B, Q2C, c2.R_feedback);                             // R15
        r(kDkVcc, Q2C, c2.R_collector_load);                    // R11
        r(TONE_IN, TB, t.R_lp_series);                          // R8
        r(TT, kDkGnd, t.R_hp_shunt);                            // R5
        r(kDkVcc, Q1B, ob.R_bias_top);                          // R7
        r(Q1B, kDkGnd, ob.R_base_to_gnd);                       // R3
        r(kDkVcc, Q1C, ob.R_collector_load);                    // R6
        r(Q1E, kDkGnd, ob.R_emitter);                           // R4
        r(SUS_TOP, SUS_W, sus.R_pot * (1.0 - k.sus) + 1e-6);    // sustain, top leg
        r(SUS_W, SUS_BOT, sus.R_pot * k.sus + 1e-6);            // sustain, bottom leg
        r(TT, TW, t.R_pot * (1.0 - k.tone) + 1e-6);             // tone, top leg
        r(TW, TB, t.R_pot * k.tone + 1e-6);                     // tone, bottom leg
        r(VOL_TOP, OUT, vol.R_pot * (1.0 - k.vol) + 1e-6);      // volume, top leg
        r(OUT, kDkGnd, vol.R_pot * k.vol + 1e-6);               // volume, bottom leg
    }

    static bjt_params bjt(const CakeVariant& v, int k) noexcept { return StandardCircuitValues::bjt(v, k); }
    static diode_params pair(const CakeVariant& v, int) noexcept { return v.clipper_1.diodes; }
};

using DkJumboOpenGoat = DkEngine<JumboOpenCircuit>;

// ---- V3 Jumbo Colour Cake ---------------------------------------------------------
//
// The Jumbo Colour Cake: three transistors (no output booster), a
// clean first gain stage (no pair), one clipper pair on Q2, C14 into the tone
// input TA, then TA -R5- TT -C8|gnd and TA -C9- TB -R8|gnd with the pot between
// TT and TB, and the wiper -C3-> the volume pot.
struct JumboColourCircuit {
    enum Node : int {
        N_IN, Q4B, Q4C, Q4E, SUS_TOP, SUS_W, SUS_BOT, C5N,
        Q3B, Q3C, Q3E, C13N,
        Q2B, Q2C, Q2AUX, Q2E,
        TA, TT, TW, TB,
        VOL_TOP, OUT,
        NODES
    };
    static constexpr int NN = NODES, NB = 3, ND = 1, S = 12, NT = 3;
    static constexpr bool FULL_EM = true;

    static constexpr std::array<DkBjtPins, NB> kBjts = {{ {Q4C, Q4B, Q4E}, {Q3C, Q3B, Q3E}, {Q2C, Q2B, Q2E} }};
    static constexpr std::array<DkPairPins, ND> kPairs = {{ {Q2AUX, Q2C} }};
    static constexpr std::array<DkCapPins, S> kCaps = {{
        {N_IN, Q4B}, {Q4B, Q4C}, {Q4C, SUS_TOP}, {SUS_W, C5N},          // C1 C10 C4 C5
        {Q3B, Q3C}, {Q3C, C13N}, {Q2B, Q2C}, {Q2B, Q2AUX},               // C12 C13 C11 C7
        {Q2C, TA}, {TA, TB}, {TT, kDkGnd}, {TW, VOL_TOP},                // C14 C9 C8 C3
    }};
    static constexpr std::array<DkTap, NT> kTaps = {{ {OUT, false}, {C13N, false}, {Q2C, true} }};
    static constexpr int kSeedHigh[] = { Q4B, Q3B, Q2B };

    static std::array<double, S> cap_values(const CakeVariant& v) noexcept {
        const auto& ib = v.input_booster; const auto& sus = v.sustain;
        const auto& c1 = v.clipper_1;     const auto& c2 = v.clipper_2;
        const auto& t  = v.tone;          const auto& ob = v.output_booster;
        return { ib.C_input_couple, ib.C_feedback_hf, ib.C_output_couple, sus.C_wiper_couple,
                 c1.C_feedback_hf, c1.C_output_couple, c2.C_feedback_hf, c2.C_diode_couple,
                 c2.C_output_couple, t.C_hp_series, t.C_lp_shunt, ob.C_input_couple };
    }

    template <class R>
    static void resistors(const CakeVariant& v, const DkKnobs& k, R&& r) noexcept {
        const auto& ib = v.input_booster; const auto& sus = v.sustain;
        const auto& c1 = v.clipper_1;     const auto& c2 = v.clipper_2;
        const auto& t  = v.tone;          const auto& vol = v.volume;
        r(kDkIn, N_IN, ib.R_input_series);                      // R2
        r(Q4B, kDkGnd, ib.R_base_to_gnd);                       // R14
        r(Q4B, Q4C, ib.R_feedback);                             // R9
        r(kDkVcc, Q4C, ib.R_collector_load);                    // R13
        r(Q4E, kDkGnd, ib.R_emitter);                           // R22
        r(C5N, Q3B, c1.R_input_series);                         // R19
        r(SUS_BOT, kDkGnd, sus.R_pot_bottom);                   // R23
        r(Q3B, kDkGnd, c1.R_base_to_gnd);                       // R20
        r(Q3B, Q3C, c1.R_feedback);                             // R17
        r(kDkVcc, Q3C, c1.R_collector_load);                    // R18
        r(Q3E, kDkGnd, c1.R_emitter);                           // R21
        r(C13N, Q2B, c2.R_input_series);                        // R12
        r(Q2B, kDkGnd, c2.R_base_to_gnd);                       // R16
        r(Q2B, Q2C, c2.R_feedback);                             // R15
        r(kDkVcc, Q2C, c2.R_collector_load);                    // R11
        r(Q2E, kDkGnd, c2.R_emitter);                           // R10
        r(SUS_TOP, SUS_W, sus.R_pot * (1.0 - k.sus) + 1e-6);    // sustain, top leg
        r(SUS_W, SUS_BOT, sus.R_pot * k.sus + 1e-6);            // sustain, bottom leg
        r(TA, TT, t.R_hp_shunt);                                // R5
        r(TB, kDkGnd, t.R_lp_series);                           // R8
        r(TT, TW, t.R_pot * (1.0 - k.tone) + 1e-6);             // tone, top leg
        r(TW, TB, t.R_pot * k.tone + 1e-6);                     // tone, bottom leg
        r(VOL_TOP, OUT, vol.R_pot * (1.0 - k.vol) + 1e-6);      // volume, top leg
        r(OUT, kDkGnd, vol.R_pot * k.vol + 1e-6);               // volume, bottom leg
    }

    static bjt_params bjt(const CakeVariant& v, int k) noexcept { return StandardCircuitValues::bjt(v, k); }
    static diode_params pair(const CakeVariant& v, int) noexcept { return v.clipper_2.diodes; }
};

using DkJumboColourCake = DkEngine<JumboColourCircuit>;

// ---- V2 Dirty Goat ------------------------------------------------------------------
//
// The Dirty Goat: a standard four-transistor circuit with R12 before
// C13 on the way to Q2's base, and an "additional cap" C14 (22 nF) in series on
// the tone pot's bottom leg (TB -C14- POT_BOT, the pot between TT and POT_BOT).
// Both pairs take clipper 1's diode card (as the node engine does).
struct DirtyGoatCircuit {
    enum Node : int {
        N_IN, Q4B, Q4C, Q4E, SUS_TOP, SUS_W, SUS_BOT, C5N,
        Q3B, Q3C, Q3AUX, Q3E, C13N,
        Q2B, Q2C, Q2AUX, Q2E,
        TT, TW, TB, POT_BOT,
        Q1B, Q1C, Q1E, VOL_TOP, OUT,
        NODES
    };
    static constexpr int NN = NODES, NB = 4, ND = 2, S = 14, NT = 3;
    static constexpr bool FULL_EM = true;
    static constexpr double kC14 = 0.022e-6;   // the extra cap on the pot bottom

    static constexpr std::array<DkBjtPins, NB> kBjts = {{
        {Q4C, Q4B, Q4E}, {Q3C, Q3B, Q3E}, {Q2C, Q2B, Q2E}, {Q1C, Q1B, Q1E},
    }};
    static constexpr std::array<DkPairPins, ND> kPairs = {{ {Q3AUX, Q3C}, {Q2AUX, Q2C} }};
    static constexpr std::array<DkCapPins, S> kCaps = {{
        {N_IN, Q4B}, {Q4B, Q4C}, {Q4C, SUS_TOP}, {SUS_W, C5N},          // C1 C10 C4 C5
        {Q3B, Q3C}, {Q3B, Q3AUX}, {C13N, Q2B},                           // C12 C6 C13
        {Q2B, Q2C}, {Q2B, Q2AUX}, {Q2C, TT}, {TB, kDkGnd},               // C11 C7 C9 C8
        {TB, POT_BOT}, {TW, Q1B}, {Q1C, VOL_TOP},                        // C14 C3 C2
    }};
    static constexpr std::array<DkTap, NT> kTaps = {{ {OUT, false}, {C13N, false}, {Q2C, true} }};
    static constexpr int kSeedHigh[] = { Q4B, Q3B, Q2B, Q1B };

    static std::array<double, S> cap_values(const CakeVariant& v) noexcept {
        const auto& ib = v.input_booster; const auto& sus = v.sustain;
        const auto& c1 = v.clipper_1;     const auto& c2 = v.clipper_2;
        const auto& t  = v.tone;          const auto& ob = v.output_booster;
        return { ib.C_input_couple, ib.C_feedback_hf, ib.C_output_couple, sus.C_wiper_couple,
                 c1.C_feedback_hf, c1.C_diode_couple, c1.C_output_couple,
                 c2.C_feedback_hf, c2.C_diode_couple, t.C_hp_series, t.C_lp_shunt,
                 kC14, ob.C_input_couple, ob.C_output_couple };
    }

    template <class R>
    static void resistors(const CakeVariant& v, const DkKnobs& k, R&& r) noexcept {
        const auto& ib = v.input_booster; const auto& sus = v.sustain;
        const auto& c1 = v.clipper_1;     const auto& c2 = v.clipper_2;
        const auto& t  = v.tone;          const auto& ob = v.output_booster;
        const auto& vol = v.volume;
        r(kDkIn, N_IN, ib.R_input_series);                      // R2
        r(Q4B, kDkGnd, ib.R_base_to_gnd);                       // R14
        r(Q4B, Q4C, ib.R_feedback);                             // R9
        r(kDkVcc, Q4C, ib.R_collector_load);                    // R13
        r(Q4E, kDkGnd, ib.R_emitter);                           // R22
        r(C5N, Q3B, c1.R_input_series);                         // R19
        r(SUS_BOT, kDkGnd, sus.R_pot_bottom);                   // R23
        r(Q3B, kDkGnd, c1.R_base_to_gnd);                       // R20
        r(Q3B, Q3C, c1.R_feedback);                             // R17
        r(kDkVcc, Q3C, c1.R_collector_load);                    // R18
        r(Q3E, kDkGnd, c1.R_emitter);                           // R21
        r(Q3C, C13N, c2.R_input_series);                        // R12
        r(Q2B, kDkGnd, c2.R_base_to_gnd);                       // R16
        r(Q2B, Q2C, c2.R_feedback);                             // R15
        r(kDkVcc, Q2C, c2.R_collector_load);                    // R11
        r(Q2E, kDkGnd, c2.R_emitter);                           // R10
        r(TT, kDkGnd, t.R_hp_shunt);                            // R5
        r(Q2C, TB, t.R_lp_series);                              // R8
        r(kDkVcc, Q1B, ob.R_bias_top);                          // R7
        r(Q1B, kDkGnd, ob.R_base_to_gnd);                       // R3
        r(kDkVcc, Q1C, ob.R_collector_load);                    // R6
        r(Q1E, kDkGnd, ob.R_emitter);                           // R4
        r(SUS_TOP, SUS_W, sus.R_pot * (1.0 - k.sus) + 1e-6);    // sustain, top leg
        r(SUS_W, SUS_BOT, sus.R_pot * k.sus + 1e-6);            // sustain, bottom leg
        r(TT, TW, t.R_pot * (1.0 - k.tone) + 1e-6);             // tone, top leg
        r(TW, POT_BOT, t.R_pot * k.tone + 1e-6);                // tone, bottom leg (to C14)
        r(VOL_TOP, OUT, vol.R_pot * (1.0 - k.vol) + 1e-6);      // volume, top leg
        r(OUT, kDkGnd, vol.R_pot * k.vol + 1e-6);               // volume, bottom leg
    }

    static bjt_params bjt(const CakeVariant& v, int k) noexcept { return StandardCircuitValues::bjt(v, k); }
    static diode_params pair(const CakeVariant& v, int) noexcept { return v.clipper_1.diodes; }
};

using DkDirtyGoat = DkEngine<DirtyGoatCircuit>;

// ---- V6 Reversed Cake, V3 Raw Cake ----------------------------------------------
//
// The tone-bypass boards: C13 before R12 into Q2, a coupling cap from Q2's
// collector to the switch common SWA, and a DPDT switch that routes SWA either
// through the tone stack (TONE_IN ... wiper TW) or around it through R28, and
// back to SWB. The switch is four contacts, 1 ohm closed and 100 Mohm open (as
// the node engines model them), so the network keeps one topology and the
// bypass knob only changes four conductances. Both pairs take clipper 1's diode
// card. The V6 couples SWB to Q1 through C3 alone; the Raw Cake has C16 and C3
// in series, and their junction carries a 100 Mohm leak to ground (SPICE's GMIN).
template <bool RawCake>
struct BypassCakeCircuit {
    enum Node : int {
        N_IN, Q4B, Q4C, Q4E, SUS_TOP, SUS_W, SUS_BOT, C5N,
        Q3B, Q3C, Q3AUX, Q3E, C13N,
        Q2B, Q2C, Q2AUX, Q2E,
        SWA, TONE_IN, TT, TW, TB, R28A, R28B, SWB,
        Q1B, Q1C, Q1E, VOL_TOP, OUT,
        C16C3,                      // only in the Raw Cake
        NODES_WITH_C16
    };
    static constexpr int NN = RawCake ? NODES_WITH_C16 : C16C3;
    static constexpr int NB = 4, ND = 2, S = RawCake ? 15 : 14, NT = 3;
    static constexpr bool FULL_EM = true;
    static constexpr double kCouple = 0.1e-6;                    // C14 (V6) / C15 (Raw Cake)
    static constexpr double kC16 = 0.1e-6;
    static constexpr double kR28 = RawCake ? 150e3 : 160e3;      // the bypass resistor
    static constexpr double kClosed = 1.0, kOpen = 1e8;          // contact resistances, ohms
    static constexpr double kLeak = 1e8;                         // the C16-C3 junction to ground

    static constexpr std::array<DkBjtPins, NB> kBjts = {{
        {Q4C, Q4B, Q4E}, {Q3C, Q3B, Q3E}, {Q2C, Q2B, Q2E}, {Q1C, Q1B, Q1E},
    }};
    static constexpr std::array<DkPairPins, ND> kPairs = {{ {Q3AUX, Q3C}, {Q2AUX, Q2C} }};
    static constexpr std::array<DkCapPins, S> make_caps() noexcept {
        std::array<DkCapPins, S> c{};
        int n = 0;
        c[n++] = {N_IN, Q4B}; c[n++] = {Q4B, Q4C}; c[n++] = {Q4C, SUS_TOP}; c[n++] = {SUS_W, C5N};   // C1 C10 C4 C5
        c[n++] = {Q3B, Q3C}; c[n++] = {Q3B, Q3AUX}; c[n++] = {Q3C, C13N};                          // C12 C6 C13
        c[n++] = {Q2B, Q2C}; c[n++] = {Q2B, Q2AUX}; c[n++] = {Q2C, SWA}; c[n++] = {TONE_IN, TT};   // C11 C7 C14|C15 C9
        c[n++] = {TB, kDkGnd};                                                                     // C8
        if constexpr (RawCake) { c[n++] = {SWB, C16C3}; c[n++] = {C16C3, Q1B}; }                  // C16 C3
        else                   { c[n++] = {SWB, Q1B}; }                                            // C3
        c[n++] = {Q1C, VOL_TOP};                                                                   // C2
        return c;
    }
    static constexpr std::array<DkCapPins, S> kCaps = make_caps();
    static constexpr std::array<DkTap, NT> kTaps = {{ {OUT, false}, {C13N, false}, {Q2C, true} }};
    static constexpr int kSeedHigh[] = { Q4B, Q3B, Q2B, Q1B };

    static std::array<double, S> cap_values(const CakeVariant& v) noexcept {
        const auto& ib = v.input_booster; const auto& sus = v.sustain;
        const auto& c1 = v.clipper_1;     const auto& c2 = v.clipper_2;
        const auto& t  = v.tone;          const auto& ob = v.output_booster;
        std::array<double, S> c{};
        int n = 0;
        c[n++] = ib.C_input_couple; c[n++] = ib.C_feedback_hf; c[n++] = ib.C_output_couple; c[n++] = sus.C_wiper_couple;
        c[n++] = c1.C_feedback_hf; c[n++] = c1.C_diode_couple; c[n++] = c1.C_output_couple;
        c[n++] = c2.C_feedback_hf; c[n++] = c2.C_diode_couple; c[n++] = kCouple; c[n++] = t.C_hp_series;
        c[n++] = t.C_lp_shunt;
        if constexpr (RawCake) c[n++] = kC16;
        c[n++] = ob.C_input_couple; c[n++] = ob.C_output_couple;
        return c;
    }

    template <class R>
    static void resistors(const CakeVariant& v, const DkKnobs& k, R&& r) noexcept {
        const auto& ib = v.input_booster; const auto& sus = v.sustain;
        const auto& c1 = v.clipper_1;     const auto& c2 = v.clipper_2;
        const auto& t  = v.tone;          const auto& ob = v.output_booster;
        const auto& vol = v.volume;
        r(kDkIn, N_IN, ib.R_input_series);                      // R2
        r(Q4B, kDkGnd, ib.R_base_to_gnd);                       // R14
        r(Q4B, Q4C, ib.R_feedback);                             // R9
        r(kDkVcc, Q4C, ib.R_collector_load);                    // R13
        r(Q4E, kDkGnd, ib.R_emitter);                           // R22
        r(C5N, Q3B, c1.R_input_series);                         // R19
        r(SUS_BOT, kDkGnd, sus.R_pot_bottom);                   // R23
        r(Q3B, kDkGnd, c1.R_base_to_gnd);                       // R20
        r(Q3B, Q3C, c1.R_feedback);                             // R17
        r(kDkVcc, Q3C, c1.R_collector_load);                    // R18
        r(Q3E, kDkGnd, c1.R_emitter);                           // R21
        r(C13N, Q2B, c2.R_input_series);                        // R12
        r(Q2B, kDkGnd, c2.R_base_to_gnd);                       // R16
        r(Q2B, Q2C, c2.R_feedback);                             // R15
        r(kDkVcc, Q2C, c2.R_collector_load);                    // R11
        r(Q2E, kDkGnd, c2.R_emitter);                           // R10
        r(TT, kDkGnd, t.R_hp_shunt);                            // R5
        r(TONE_IN, TB, t.R_lp_series);                          // R8
        r(R28A, R28B, kR28);                                    // R28, the bypass path
        r(kDkVcc, Q1B, ob.R_bias_top);                          // R7
        r(Q1B, kDkGnd, ob.R_base_to_gnd);                       // R3
        r(kDkVcc, Q1C, ob.R_collector_load);                    // R6
        r(Q1E, kDkGnd, ob.R_emitter);                           // R4
        const double through = k.bypass ? kOpen : kClosed, around = k.bypass ? kClosed : kOpen;
        r(SWA, TONE_IN, through);                               // SW1A to the tone stack
        r(SWA, R28A, around);                                   // SW1A to the bypass
        r(SWB, TW, through);                                    // SW1B from the tone wiper
        r(SWB, R28B, around);                                   // SW1B from the bypass
        r(SUS_TOP, SUS_W, sus.R_pot * (1.0 - k.sus) + 1e-6);    // sustain, top leg
        r(SUS_W, SUS_BOT, sus.R_pot * k.sus + 1e-6);            // sustain, bottom leg
        r(TT, TW, t.R_pot * (1.0 - k.tone) + 1e-6);             // tone, top leg
        r(TW, TB, t.R_pot * k.tone + 1e-6);                     // tone, bottom leg
        r(VOL_TOP, OUT, vol.R_pot * (1.0 - k.vol) + 1e-6);      // volume, top leg
        r(OUT, kDkGnd, vol.R_pot * k.vol + 1e-6);               // volume, bottom leg
        if constexpr (RawCake) r(C16C3, kDkGnd, kLeak);         // the cap-only junction's leak
    }

    static bjt_params bjt(const CakeVariant& v, int k) noexcept { return StandardCircuitValues::bjt(v, k); }
    static diode_params pair(const CakeVariant& v, int) noexcept { return v.clipper_1.diodes; }
};

using DkReversedCake = DkEngine<BypassCakeCircuit<false>>;
using DkRawCake  = DkEngine<BypassCakeCircuit<true>>;

// ---- V2 Goat&S OD -------------------------------------------------------------------
//
// The Goat&S OD: a four-transistor circuit whose input cap C1 sits
// on the jack side of R2, and whose tone stack is an active network around
// Q1: C14 from Q2's collector to the pot top TA; the wiper TW -C9- R8 -> Y -C3->
// Q1's base; a ladder TA -R5- X1 -C15- X2 -C16- X3 -R30- TB with R28 across C15,
// R29 across C16 and R27 from Y to X2; Q1's collector -C2-> TB; TB -R31-> the
// 50 k volume pot; Q1's emitter bypassed by C17. Both pairs take clipper 1's
// diode card. Every value the variant does not carry is the node engine's.
struct GoatSODCircuit {
    enum Node : int {
        N_C1, Q4B, Q4C, Q4E, SUS_TOP, SUS_W, SUS_BOT, C5N,
        Q3B, Q3C, Q3AUX, Q3E, C13N,
        Q2B, Q2C, Q2AUX, Q2E,
        TA, TW, C9N, Y, TB, X1, X2, X3,
        Q1B, Q1C, Q1E, VOL_TOP, OUT,
        NODES
    };
    static constexpr int NN = NODES, NB = 4, ND = 2, S = 16, NT = 3;
    static constexpr bool FULL_EM = true;
    static constexpr double kR5 = 10e3, kR8 = 22e3, kR27 = 10e3, kR28 = 47e3;
    static constexpr double kR29 = 47e3, kR30 = 10e3, kR31 = 47e3, kRvol = 50e3;
    static constexpr double kC9 = 0.0022e-6, kC15 = 0.047e-6, kC16 = 0.047e-6, kC17 = 10e-6;

    static constexpr std::array<DkBjtPins, NB> kBjts = {{
        {Q4C, Q4B, Q4E}, {Q3C, Q3B, Q3E}, {Q2C, Q2B, Q2E}, {Q1C, Q1B, Q1E},
    }};
    static constexpr std::array<DkPairPins, ND> kPairs = {{ {Q3AUX, Q3C}, {Q2AUX, Q2C} }};
    static constexpr std::array<DkCapPins, S> kCaps = {{
        {N_C1, kDkIn}, {Q4B, Q4C}, {Q4C, SUS_TOP}, {SUS_W, C5N},        // C1 C10 C4 C5
        {Q3B, Q3C}, {Q3B, Q3AUX}, {Q3C, C13N},                           // C12 C6 C13
        {Q2B, Q2C}, {Q2B, Q2AUX}, {Q2C, TA}, {TW, C9N},                  // C11 C7 C14 C9
        {X1, X2}, {X2, X3}, {Y, Q1B}, {Q1C, TB}, {Q1E, kDkGnd},          // C15 C16 C3 C2 C17
    }};
    static constexpr std::array<DkTap, NT> kTaps = {{ {OUT, false}, {C13N, false}, {Q2C, true} }};
    static constexpr int kSeedHigh[] = { Q4B, Q3B, Q2B, Q1B };

    static std::array<double, S> cap_values(const CakeVariant& v) noexcept {
        const auto& ib = v.input_booster; const auto& sus = v.sustain;
        const auto& c1 = v.clipper_1;     const auto& c2 = v.clipper_2;
        const auto& ob = v.output_booster;
        return { ib.C_input_couple, ib.C_feedback_hf, ib.C_output_couple, sus.C_wiper_couple,
                 c1.C_feedback_hf, c1.C_diode_couple, c1.C_output_couple,
                 c2.C_feedback_hf, c2.C_diode_couple, c2.C_output_couple, kC9,
                 kC15, kC16, ob.C_input_couple, ob.C_output_couple, kC17 };
    }

    template <class R>
    static void resistors(const CakeVariant& v, const DkKnobs& k, R&& r) noexcept {
        const auto& ib = v.input_booster; const auto& sus = v.sustain;
        const auto& c1 = v.clipper_1;     const auto& c2 = v.clipper_2;
        const auto& t  = v.tone;          const auto& ob = v.output_booster;
        r(N_C1, Q4B, ib.R_input_series);                        // R2
        r(Q4B, kDkGnd, ib.R_base_to_gnd);                       // R14
        r(Q4B, Q4C, ib.R_feedback);                             // R9
        r(kDkVcc, Q4C, ib.R_collector_load);                    // R13
        r(Q4E, kDkGnd, ib.R_emitter);                           // R22
        r(C5N, Q3B, c1.R_input_series);                         // R19
        r(SUS_BOT, kDkGnd, sus.R_pot_bottom);                   // R23
        r(Q3B, kDkGnd, c1.R_base_to_gnd);                       // R20
        r(Q3B, Q3C, c1.R_feedback);                             // R17
        r(kDkVcc, Q3C, c1.R_collector_load);                    // R18
        r(Q3E, kDkGnd, c1.R_emitter);                           // R21
        r(C13N, Q2B, c2.R_input_series);                        // R12
        r(Q2B, kDkGnd, c2.R_base_to_gnd);                       // R16
        r(Q2B, Q2C, c2.R_feedback);                             // R15
        r(kDkVcc, Q2C, c2.R_collector_load);                    // R11
        r(Q2E, kDkGnd, c2.R_emitter);                           // R10
        r(C9N, Y, kR8);                                         // R8
        r(TA, X1, kR5);                                         // R5
        r(X1, X2, kR28);                                        // R28
        r(Y, X2, kR27);                                         // R27
        r(X2, X3, kR29);                                        // R29
        r(X3, TB, kR30);                                        // R30
        r(kDkVcc, Q1B, ob.R_bias_top);                          // R7
        r(Q1B, kDkGnd, ob.R_base_to_gnd);                       // R3
        r(kDkVcc, Q1C, ob.R_collector_load);                    // R6
        r(Q1E, kDkGnd, ob.R_emitter);                           // R4
        r(TB, VOL_TOP, kR31);                                   // R31
        r(SUS_TOP, SUS_W, sus.R_pot * (1.0 - k.sus) + 1e-6);    // sustain, top leg
        r(SUS_W, SUS_BOT, sus.R_pot * k.sus + 1e-6);            // sustain, bottom leg
        r(TA, TW, t.R_pot * (1.0 - k.tone) + 1e-6);             // tone, top leg
        r(TW, TB, t.R_pot * k.tone + 1e-6);                     // tone, bottom leg
        r(VOL_TOP, OUT, kRvol * (1.0 - k.vol) + 1e-6);          // volume (50 k), top leg
        r(OUT, kDkGnd, kRvol * k.vol + 1e-6);                   // volume, bottom leg
    }

    static bjt_params bjt(const CakeVariant& v, int k) noexcept { return StandardCircuitValues::bjt(v, k); }
    static diode_params pair(const CakeVariant& v, int) noexcept { return v.clipper_1.diodes; }
};

using DkGoatSOD = DkEngine<GoatSODCircuit>;

// ---- V1 Master Pyramid (PM-3, a 1970s Japanese fuzz with its boost removed) --------
//
// The four-transistor core with five departures: C15 across the input resistor
// R2 (a bright cap, so the node IN_C1 between them and C1); R27 from the input
// booster's output to ground; C14 coupling Q2's collector to the tone stack
// (TIN, feeding C9 and R8); the volume pot moved between the tone wiper and C3;
// and R28 as the output pulldown after C2. C5 then R19 and C13 then R12, as in
// the Standard circuit. Parts the variant does not carry are the node engine's.
struct PyramidMasterCircuit {
    enum Node : int {
        IN_C1, Q4B, Q4C, Q4E, SUS_TOP, SUS_W, SUS_BOT, C5_R19,
        Q3B, Q3C, Q3E, Q3_AUX, C13_R12,
        Q2B, Q2C, Q2E, Q2_AUX, TIN,
        TONE_TOP, TONE_BOT, TONE_W, VOL_W,
        Q1B, Q1C, Q1E, OUT,
        NODES
    };
    static constexpr int NN = NODES, NB = 4, ND = 2, S = 15, NT = 4;
    static constexpr bool FULL_EM = true;
    static constexpr double kC15 = 3300e-12, kC14 = 1e-6, kR27 = 100e3, kR28 = 100e3;
    // Linear extrapolation seeds this circuit's Newton better than quadratic:
    // paired on the .so, 0.935-0.977 of order 2 (7/8 pairs), where every
    // other family measured got slower with it.
    static constexpr int kPredictor = 1;

    static constexpr std::array<DkBjtPins, NB> kBjts = {{
        {Q4C, Q4B, Q4E}, {Q3C, Q3B, Q3E}, {Q2C, Q2B, Q2E}, {Q1C, Q1B, Q1E},
    }};
    static constexpr std::array<DkPairPins, ND> kPairs = {{ {Q3_AUX, Q3C}, {Q2_AUX, Q2C} }};
    static constexpr std::array<DkCapPins, S> kCaps = {{
        {IN_C1, kDkIn}, {IN_C1, Q4B}, {Q4B, Q4C}, {Q4C, SUS_TOP}, {SUS_W, C5_R19},  // C15 C1 C10 C4 C5
        {Q3B, Q3C}, {Q3B, Q3_AUX}, {Q3C, C13_R12},                                   // C12 C6 C13
        {Q2B, Q2C}, {Q2B, Q2_AUX}, {Q2C, TIN}, {TIN, TONE_TOP}, {TONE_BOT, kDkGnd},  // C11 C7 C14 C9 C8
        {VOL_W, Q1B}, {Q1C, OUT},                                                    // C3 C2
    }};
    // Out: after C2 and R28. C1Out: the node after C13. C2Out: Q2's collector
    // about its rest. IbOut: the sustain-pot top (the input booster's output).
    static constexpr std::array<DkTap, NT> kTaps = {{ {OUT, false}, {C13_R12, false}, {Q2C, true}, {SUS_TOP, false} }};
    static constexpr int kSeedHigh[] = { Q4B, Q3B, Q2B, Q1B };

    static std::array<double, S> cap_values(const CakeVariant& v) noexcept {
        const auto& ib = v.input_booster; const auto& sus = v.sustain;
        const auto& c1 = v.clipper_1;     const auto& c2 = v.clipper_2;
        const auto& t  = v.tone;          const auto& ob = v.output_booster;
        return { kC15, ib.C_input_couple, ib.C_feedback_hf, ib.C_output_couple, sus.C_wiper_couple,
                 c1.C_feedback_hf, c1.C_diode_couple, c1.C_output_couple,
                 c2.C_feedback_hf, c2.C_diode_couple, kC14, t.C_hp_series, t.C_lp_shunt,
                 ob.C_input_couple, ob.C_output_couple };
    }

    template <class R>
    static void resistors(const CakeVariant& v, const DkKnobs& k, R&& r) noexcept {
        const auto& ib = v.input_booster; const auto& sus = v.sustain;
        const auto& c1 = v.clipper_1;     const auto& c2 = v.clipper_2;
        const auto& t  = v.tone;          const auto& ob = v.output_booster;
        const auto& vol = v.volume;
        r(kDkIn, IN_C1, ib.R_input_series);                     // R2 (C15 across it)
        r(Q4B, kDkGnd, ib.R_base_to_gnd);                       // R14
        r(Q4B, Q4C, ib.R_feedback);                             // R9
        r(kDkVcc, Q4C, ib.R_collector_load);                    // R13
        r(Q4E, kDkGnd, ib.R_emitter);                           // R22
        r(SUS_TOP, kDkGnd, kR27);                               // R27
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
        r(TIN, TONE_BOT, t.R_lp_series);                        // R8, from the C14 node
        r(TONE_TOP, TONE_W, t.R_pot * (1.0 - k.tone) + 1e-6);   // tone, top leg
        r(TONE_W, TONE_BOT, t.R_pot * k.tone + 1e-6);           // tone, bottom leg
        r(TONE_W, VOL_W, vol.R_pot * (1.0 - k.vol) + 1e-6);     // volume, top leg (from the tone wiper)
        r(VOL_W, kDkGnd, vol.R_pot * k.vol + 1e-6);             // volume, bottom leg
        r(kDkVcc, Q1B, ob.R_bias_top);                          // R7
        r(Q1B, kDkGnd, ob.R_base_to_gnd);                       // R3
        r(kDkVcc, Q1C, ob.R_collector_load);                    // R6
        r(Q1E, kDkGnd, ob.R_emitter);                           // R4
        r(OUT, kDkGnd, kR28);                                   // R28, output pulldown
    }

    static bjt_params bjt(const CakeVariant& v, int k) noexcept { return StandardCircuitValues::bjt(v, k); }
    static diode_params pair(const CakeVariant& v, int k) noexcept { return StandardCircuitValues::pair(v, k); }
};

using DkPyramidMaster = DkEngine<PyramidMasterCircuit>;

// ---- V2 Jumbo Goat ----------------------------------------------------------------------
//
// The Jumbo Goat: the Jumbo Open Goat's main path plus a germanium side
// path off the tone input.
// The two germanium gate diodes are NOT an antiparallel pair (read off the
// published schematic and the KiCad netlist): ONE of them rectifies the tone input
// into C15 || R27, the other sits between GROUND and the tone input and clamps
// the negative half. They are two SINGLE diodes (kSingles, the engine's
// one-way port).
// R28 feeds Q5 (AC127, germanium), whose collector (R29 to the rail, C16 to
// ground) drives Q6's base; Q6's collector IS the tone wiper, so its current sums
// with the passive tone output before C3. Q3, Q2, Q5 and Q6 have their emitters
// on ground. Six transistors, two diode pairs and two single diodes: 16 ports. The germanium models are
// the node engine's (Q5: see jumbo_goat.h).
struct JumboGoatCircuit {
    // The AC127 gate (Is 1e-6; the figures below were measured at 6e-8) sits on a low-current junction
    // where a 1e-3 V step criterion stops one iteration early: -58.1 dB against
    // the converged engine (ceiling -90).
    // Only the gate's junctions carry it: Q5 (ports 6-7), Q6 (8-9) and the two
    // gate diodes (14-15) step under 1e-4 V, the main path under the global
    // criterion (dk_circuit_tight). -87.3 dB against the converged engine for
    // 1.23 iterations per sample (all 16 ports tight: -96.0 dB, 1.31).
    static constexpr std::array<int, 6> kTightPorts = { 6, 7, 8, 9, 14, 15 };
    static constexpr double kTightTol = 1.0e-4;
    enum Node : int {
        Q4B, Q4C, Q4E, N_C1,
        SUS_TOP, SUS_W, SUS_BOT, C5N,
        Q3B, Q3C, Q3AUX, C13N,
        Q2B, Q2C, Q2AUX,
        TONE_IN, TT, TB, TW,
        D5K, Q5B, Q5C,
        Q1B, Q1C, Q1E,
        VOL_TOP, OUT,
        NODES
    };
    static constexpr int NN = NODES, NB = 6, ND = 2, S = 17, NT = 3;
    static constexpr int NS = 2;          // the two gate diodes, one way each
    static constexpr bool FULL_EM = true;
    // The only circuit here whose port Jacobian is not structurally dense: the
    // germanium gate and its envelope node (ports 6-8) barely couple to the
    // rest, so 58 of 256 entries are exactly zero at every rate, knob and
    // switch position. The generated elimination skips exactly those
    // operations and is bit-identical -- `x -= f * 0.0` is `x`.
    using DkLu = dkgen::DkLu_jumbogoat;
    static constexpr double kR27 = 6.8e3, kR28 = 220.0, kR29 = 22e3;
    static constexpr double kC15 = 4.7e-6, kC16 = 0.1e-6, kC17 = 470e-12;
    // AC127: Is 1e-6, a leaky germanium part, set so the gate opens on playing and closes in silence on a
    // real guitar signal (the class value from Vbe ~0.25 V at 1 mA, 6e-8, keeps it shut most of the time).
    // The Pro Electron sheet gives hFE >= 50 and ICBO <= 10 uA but no VBE, so it cannot settle Is. Bf 100
    // sits in that hFE band.
#ifndef TOCAKE_AC127_IS
#define TOCAKE_AC127_IS 1.0e-6
#endif
    static constexpr bjt_params kQ5Ge { TOCAKE_AC127_IS, 100.0, 50.0 };
    static constexpr diode_params kGate { 2.0e-7, 1.06, 1, 1 };   // D5/D6, germanium

    static constexpr std::array<DkBjtPins, NB> kBjts = {{
        {Q4C, Q4B, Q4E}, {Q3C, Q3B, kDkGnd}, {Q2C, Q2B, kDkGnd},
        {Q5C, Q5B, kDkGnd},         // Q5, germanium
        {TW, Q5C, kDkGnd},          // Q6: base on Q5's collector, collector on the tone wiper
        {Q1C, Q1B, Q1E},
    }};
    static constexpr std::array<DkPairPins, ND> kPairs = {{ {Q3AUX, Q3C}, {Q2AUX, Q2C} }};
    // D6 rectifies the tone input into the envelope node; D5 clamps the
    // negative half of the tone input against ground (anode, cathode).
    static constexpr std::array<DkSinglePins, NS> kSingles = {{ {TONE_IN, D5K}, {kDkGnd, TONE_IN} }};
    static constexpr std::array<DkCapPins, S> kCaps = {{
        {N_C1, Q4B}, {Q4B, Q4C}, {Q4C, SUS_TOP}, {SUS_W, C5N},          // C1 C10 C4 C5
        {Q3B, Q3C}, {Q3B, Q3AUX}, {C13N, Q2B},                           // C12 C6 C13
        {Q2B, Q2C}, {Q2B, Q2AUX}, {Q2C, TONE_IN}, {TONE_IN, TT},         // C11 C7 C14 C9
        {TB, kDkGnd}, {D5K, kDkGnd}, {Q5C, kDkGnd},                      // C8 C15 C16
        {TW, Q1B}, {Q1C, VOL_TOP}, {Q1B, Q1C},                           // C3 C2 C17
    }};
    static constexpr std::array<DkTap, NT> kTaps = {{ {OUT, false}, {C13N, false}, {Q2C, true} }};
    static constexpr int kSeedHigh[] = { Q4B, Q3B, Q2B, Q1B };

    static std::array<double, S> cap_values(const CakeVariant& v) noexcept {
        const auto& ib = v.input_booster; const auto& sus = v.sustain;
        const auto& c1 = v.clipper_1;     const auto& c2 = v.clipper_2;
        const auto& t  = v.tone;          const auto& ob = v.output_booster;
        return { ib.C_input_couple, ib.C_feedback_hf, ib.C_output_couple, sus.C_wiper_couple,
                 c1.C_feedback_hf, c1.C_diode_couple, c1.C_output_couple,
                 c2.C_feedback_hf, c2.C_diode_couple, c2.C_output_couple, t.C_hp_series,
                 t.C_lp_shunt, kC15, kC16, ob.C_input_couple, ob.C_output_couple, kC17 };
    }

    template <class R>
    static void resistors(const CakeVariant& v, const DkKnobs& k, R&& r) noexcept {
        const auto& ib = v.input_booster; const auto& sus = v.sustain;
        const auto& c1 = v.clipper_1;     const auto& c2 = v.clipper_2;
        const auto& t  = v.tone;          const auto& ob = v.output_booster;
        const auto& vol = v.volume;
        r(kDkIn, N_C1, ib.R_input_series);                      // R2
        r(Q4B, kDkGnd, ib.R_base_to_gnd);                       // R14
        r(Q4B, Q4C, ib.R_feedback);                             // R9
        r(kDkVcc, Q4C, ib.R_collector_load);                    // R13
        r(Q4E, kDkGnd, ib.R_emitter);                           // R22
        r(C5N, Q3B, c1.R_input_series);                         // R19
        r(SUS_BOT, kDkGnd, sus.R_pot_bottom);                   // R23
        r(Q3B, kDkGnd, c1.R_base_to_gnd);                       // R20
        r(Q3B, Q3C, c1.R_feedback);                             // R17
        r(kDkVcc, Q3C, c1.R_collector_load);                    // R18
        r(C13N, Q3C, c2.R_input_series);                        // R12
        r(Q2B, kDkGnd, c2.R_base_to_gnd);                       // R16
        r(Q2B, Q2C, c2.R_feedback);                             // R15
        r(kDkVcc, Q2C, c2.R_collector_load);                    // R11
        r(TONE_IN, TB, t.R_lp_series);                          // R8
        r(TT, kDkGnd, t.R_hp_shunt);                            // R5
        r(D5K, kDkGnd, kR27);                                   // R27
        r(D5K, Q5B, kR28);                                      // R28
        r(kDkVcc, Q5C, kR29);                                   // R29
        r(kDkVcc, Q1B, ob.R_bias_top);                          // R7
        r(Q1B, kDkGnd, ob.R_base_to_gnd);                       // R3
        r(kDkVcc, Q1C, ob.R_collector_load);                    // R6
        r(Q1E, kDkGnd, ob.R_emitter);                           // R4
        r(SUS_TOP, SUS_W, sus.R_pot * (1.0 - k.sus) + 1e-6);    // sustain, top leg
        r(SUS_W, SUS_BOT, sus.R_pot * k.sus + 1e-6);            // sustain, bottom leg
        r(TT, TW, t.R_pot * (1.0 - k.tone) + 1e-6);             // tone, top leg
        r(TW, TB, t.R_pot * k.tone + 1e-6);                     // tone, bottom leg
        r(VOL_TOP, OUT, vol.R_pot * (1.0 - k.vol) + 1e-6);      // volume, top leg
        r(OUT, kDkGnd, vol.R_pot * k.vol + 1e-6);               // volume, bottom leg
    }

    static bjt_params bjt(const CakeVariant& v, int k) noexcept {
        switch (k) {
            case 0: return v.input_booster.q;
            case 1: return v.clipper_1.q;
            case 2: return v.clipper_2.q;
            case 3: return kQ5Ge;
            case 4: return v.input_booster.q;       // Q6: the same Si part as Q4
            default: return v.output_booster.q;
        }
    }
    static diode_params pair(const CakeVariant& v, int k) noexcept {
        (void) k; return v.clipper_1.diodes;
    }
    static diode_params single(const CakeVariant& v, int k) noexcept {
        (void) v; (void) k; return kGate;
    }
};

using DkJumboGoat = DkEngine<JumboGoatCircuit>;

}  // namespace tocake
