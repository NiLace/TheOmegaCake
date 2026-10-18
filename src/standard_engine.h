// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// The engine that ships for the 27 Standard slots, named once.
//
// DkStandard: a converged Newton on 10 port voltages.
#pragma once

#include "dk_standard.h"
#include "engine_select.h"
namespace tocake { using StandardEngine = shipped_engine_t<DkStandard>; }   // the DK, or its nodal twin (engine_select.h)
#define TOCAKE_STANDARD_ENGINE_NAME "NodalStandard"
