// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// The transistor law of the Standard circuit (dk_standard.h, FULL_EM): SPICE's
// Gummel-Poon with the reverse term, Br = 1, as the reference decks run it. The
// forward-only law drops the reverse (Vbc) term, which is exact only while no
// transistor saturates; some Standard slots saturate at a hot input.
#pragma once

