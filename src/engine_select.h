// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// Which engine a family SHIPS: its DK description (src/dk_engine.h) or the
// nodal twin of the same description
// (src/nodal_engine.h). Every family header names its shipped
// engine through shipped_engine_t<Dk...>, so the switch cannot miss one.
#pragma once

#include "dk_engine.h"
#include "nodal_engine.h"

namespace tocake {

template <class E> struct shipped_engine { using type = E; };
template <class C> struct shipped_engine<DkEngine<C>> { using type = NodalEngine<C>; };
// Per circuit, the faster of the two engines: the two smallest port counts stay on the DK (Coloured Goat,
// P = 9; Jumbo Colour Cake, P = 7), where the nodal form is slower; every other DK family ships the nodal twin.
template <> struct shipped_engine<DkEngine<C14CakeCircuit<false>>> { using type = DkEngine<C14CakeCircuit<false>>; };
template <> struct shipped_engine<DkEngine<JumboColourCircuit>>    { using type = DkEngine<JumboColourCircuit>; };
template <class E> using shipped_engine_t = typename shipped_engine<E>::type;

}  // namespace tocake
