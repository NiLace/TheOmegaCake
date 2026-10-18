// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// The bespoke node engines may stop their Newton on the residual test from the
// first step on: the absolute 1e-7 A test is met after that step on 100-470k
// nodes, so they take about one step per sample on plucked notes.
#pragma once
