# The Omega Cake

A fuzz pedal plugin (LV2, Linux) with **45 models** of one family of 1970s-and-later fuzz circuits:
four-stage transistor fuzzes with a passive tone stack, their op-amp descendants, and the clones and
mods that grew around them. Each model is simulated from its schematic (transistors, diodes and
op-amps solved as a circuit, sample by sample) and checked against a SPICE simulation of the same
circuit.

- Mono in, mono out. One plugin, one model selector.
- Three knobs (**Sustain**, **Tone**, **Volume**), plus the switches the original circuits had.
- An **HQ** switch (4x oversampling on every model).
- A footswitch with a latency-matched, faded bypass.
- A photoreal X11 editor.

Full guide: [`docs/USER_GUIDE.md`](docs/USER_GUIDE.md) · every model with its year and engine:
[`docs/MODEL_LIST.md`](docs/MODEL_LIST.md) · release notes: [`CHANGELOG.md`](CHANGELOG.md).

## The models

The list in the editor groups them as below. The number is the `model` port value; the year is the
year of the circuit.

**Classics**: 0 The Pyramid (1971) · 14 Violet Goat (1973) · 29 Pumpkin Cake (1978) ·
34 Crimson Army (1993) · 36 Black Cake (1995) · 37 NYC Cake (2000)

| era | original maker | clones and mods |
|---|---|---|
| **1969-72** | 1 Gain Pyramid · 2 Scooped Pyramid · 3 Fat Pyramid · 4 Green Pyramid · 7 Bass Pyramid · 8 Inverted Pyramid | 5 Vintage Pyramid |
| **1973-75** | 11 Alpha Goat · 12 Beta Violet Goat · 13 Pink Goat · 15 Blue Goat · 16 White Goat | 6 Master Pyramid · 9 Super Pyramid · 20 Coloured Goat · 21 Jumbo Goat · 22 Jumbo Open Goat · 23 Goat&S OD · 24 Dirty Goat |
| **1976-79** | 17 Singer Goat · 18 Boxed Goat · 19 Beaver Goat · 43 Little Pumpkin · 28 Classic Cake · 30 Mini Pumpkin · 26 Raw Cake · 32 Reversed Cake | 27 Jumbo Colour Cake · 31 Guru Cake · 25 Distortion & Cake |
| **1990-now** | 33 Red Army · 35 Cake War · 38 Brooklyn Cake · 44 Lost Pumpkin | 39 Inflatable Johnson · 42 Dreamy Cake · 41 German Cake · 40 Precision Cake · 10 Super New Pyramid |

## Controls (ports)

| port | symbol | range (default) | what it does |
|---|---|---|---|
| Sustain | `sustain` | 0-1 (0.8) | the fuzz: gain into the clipping stages |
| Tone | `tone` | 0-1 (0.5) | the tone control: clockwise brighter (Dreamy Cake: darker, as on its circuit) |
| Volume | `volume` | 0-1 (0.6) | output level |
| HQ | `hq` | off/on (off) | on: every model at 4x oversampling; off: each model at its own default (2x, or 4x for the op-amp models) |
| Model | `model` | 0-44 (0, The Pyramid) | which circuit |
| Tone Bypass | `tone_bypass` | off/on (off) | takes the tone stack out: Raw Cake, Reversed Cake, Pumpkin Cake, Guru Cake only |
| Tone Switch (Little Pumpkin) | `tone_switch` | Dark/Bright (Dark) | the two-way tone switch of Mini Pumpkin and Little Pumpkin, which have no Tone or Sustain knob |
| Footswitch (Enabled) | `enabled` | off/on (on) | the plugin's bypass, designated `lv2:enabled` so a host's own bypass button drives it |
| Latency | `latency` | output | reported latency: 31 samples at the host rate, on every model and setting |

**Footswitch.** The bypassed signal is the input delayed by the same 31 samples, so the latency the
host compensates never changes. Engaging and bypassing crossfade over about 10.7 ms. While bypassed
the circuit stops computing (no CPU) and holds its state.

## Install

Build (Linux, x86-64; needs CMake 3.21+, a C++20 compiler, pkg-config, the LV2 headers, and the
development packages of `cairo`, `cairo-xlib`, `x11` and `freetype2`):

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

This assembles the bundle in `build/TheOmegaCake.lv2/`. Install it by copying the whole bundle:

```sh
rm -rf ~/.lv2/TheOmegaCake.lv2
cp -r build/TheOmegaCake.lv2 ~/.lv2/
```

Restart the host afterwards: hosts cache the plugin binary. The plugin URI is
`https://nylarea.com/plugins/omegacake`.

The default build targets generic x86-64. It carries a second, faster copy of the hot code for newer
CPUs and picks it at run time, so the same binary runs on older CPUs too.

## Third-party parts

- **Fonts** in the bundle (Doto, Space Mono) are under the SIL Open Font License 1.1; their licence
  texts travel beside them. They stay under the OFL and are not relicensed.
- **Interface images** (`gui/photo/`) are renders of our own 3D scene of the pedal, lit with the
  `studio_small_08` HDRI from [Poly Haven](https://polyhaven.com/a/studio_small_08) (CC0). The
  lettering in the renders is set in Holtwood One SC and Doto (both SIL Open Font License 1.1). The
  front drawing and the rest of the interface are covered by this project's licence.

## How this was made, and the licence

This plugin was built mainly with **Claude Code** (Anthropic's AI coding agent), working from my
direction on the DSP, the design and the testing. I have not reviewed every line by hand, so please
treat it accordingly:

- **No warranty.** It works in my own testing, but it carries no guarantees of any kind.
- **No support.** This is a personal project; I can't promise help, fixes or answers.
- **Help is welcome.** If you read code, spot bugs or want to improve it, review and contributions are
  genuinely appreciated.

Released under the **GNU Affero General Public License v3.0, that version only** (`AGPL-3.0-only`) —
see [`LICENSE`](LICENSE). You may use it, study it, modify it and sell it; what nobody may do is close
it: any version you distribute, or let others use over a network, has to come with its complete source
under this same licence. The third-party parts listed above keep their own licences. Contributions are
welcome — see [`CONTRIBUTING.md`](CONTRIBUTING.md).

Author: NiLace <nilace@nylarea.com> (shown in hosts as "NL Sounds").
