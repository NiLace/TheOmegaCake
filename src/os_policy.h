// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// The oversampling factor each slot runs at, with the HQ switch off and on.
// The plugin (lv2.cpp) and the engine twin that renders the golden corpus
// (render_corpus.cpp) both read it from here, so they cannot drift apart.
//
// HQ off, the defaults, chosen by the worst-frame ANMR on real material
// (criterion: below 0 dB; reference: the same engine at 384 kHz, so it judges
// the rate, not the model; docs/anmr_real_*.tsv, gated by anmr_real_check.py):
//   * 2x for every transistor slot, the Jumbo Goat included: its germanium noise
//     gate (a half-conducting Q6 on the tone wiper) reads -3.19 dB at 2x on real
//     material, under the criterion; on the SPICE grid it reads -7.03 dB at 2x
//     against -13.52 at 3x, which HQ restores;
//   * 4x for the op-amp slots: their LM741 macromodel is numerically unstable
//     below 4x (it diverges to NaN at 2x under hard drive: 13 of 32 stress
//     settings at 2x, 0 of 32 at 4x).
// HQ on: 4x everywhere.
#pragma once

#include "variant_registry_generated.h"

namespace tocake {

inline int slot_factor(VariantKind kind, bool hq) noexcept
{
    if (hq) return 4;
    switch (kind) {
        case VariantKind::OpAmp:
        case VariantKind::OpAmpLittle:
        case VariantKind::OpAmpGuru: return 4;
        // The Lost Pumpkin: its own engine (data-sheet 4558s, their compensation integrated by BDF2, so no
        // stability floor). 2x, like the transistor slots.
        case VariantKind::LostPumpkin: return 2;
        default:                     return 2;
    }
}

} // namespace tocake
