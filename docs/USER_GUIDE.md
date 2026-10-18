# The Omega Cake — user guide

This guide covers what the controls do on each model, the switches, latency, sample rates and CPU
cost. The model list with years and engines is in [`MODEL_LIST.md`](MODEL_LIST.md); install
instructions are in the [`README`](../README.md).

## 1. The signal path

The plugin is the pedal and nothing else: one mono input, one mono output. It takes its input as the
voltage at the pedal's input jack, so the level you feed it matters, as it does with the real circuit. A
hotter source (active pickups, a booster in front) drives the first stage harder. The plugin has no
input trim, input meter or pickup emulation; set the level earlier in the chain.

Each model is a circuit simulation: the transistors, diodes and op-amps are solved together, sample by
sample, at an oversampled rate (see §4). Clipping, bias shift and the interaction between the knobs come
from the circuit, not from a curve fitted to it.

**Level matching.** Each model has a fixed output gain after the circuit, measured on a guitar DI so that,
with Volume at full, every model has the same MEDIAN loudness over its Sustain × Tone travel (Sustain
0.25 / 0.5 / 1 × Tone 0 / 0.5 / 1), and none peaks above −2 dBTP at any knob or switch setting.
Switching models therefore does not jump in level on average, even though the original circuits differ
widely. At a given knob setting the models still differ, because each circuit's tone and gain controls
change its level in their own way: across the models, eight in ten sit within about 2 dB of each other at
Tone noon and within 3.5-6 dB at the Tone extremes. One model has more: the **Lost Pumpkin**'s
tone network loses far more level towards bright than the others' (about 4 dB more at Tone fully
clockwise), so it carries a second gain after the circuit that follows its Tone knob and keeps it on the
typical model's level at every Tone position. Its Tone therefore changes the colour as the original does,
but not the volume the way the original does. This gain sits after the circuit
and does not change how a model clips.

## 2. The knobs

All three knobs run from 0 (fully anticlockwise) to 1 (fully clockwise). Each one turns a
potentiometer in the simulated circuit, through the taper that model's schematic gives the pot:

- **Linear**: the wiper moves in proportion to the knob.
- **Log** (audio taper): the wiper position is the knob position squared, so at noon the wiper sits
  at 25 % of its travel. This curve approximates an audio-taper pot.

Knob changes are smoothed (a glide of about 20 ms), so automation and fast turns do not click.

### Knob laws per model

| models (port value) | Sustain | Tone | Volume |
|---|---|---|---|
| The Pyramid (0), Gain Pyramid (1), Scooped Pyramid (2), Fat Pyramid (3), Bass Pyramid (7), Inverted Pyramid (8), Alpha Goat (11), Beta Violet Goat (12), Pink Goat (13), Violet Goat (14), Blue Goat (15), White Goat (16), Singer Goat (17), Boxed Goat (18), Beaver Goat (19), Jumbo Goat (21), Jumbo Open Goat (22), Goat&S OD (23), Dirty Goat (24), Raw Cake (26), Jumbo Colour Cake (27), Classic Cake (28), Reversed Cake (32), Red Army (33), Crimson Army (34), Cake War (35), Black Cake (36), NYC Cake (37), Brooklyn Cake (38), Inflatable Johnson (39), Precision Cake (40), Dreamy Cake (42), Lost Pumpkin (44) | Linear | Linear | Linear |
| Green Pyramid (4), Vintage Pyramid (5), Super Pyramid (9), Super New Pyramid (10), Coloured Goat (20) | Log | Log | Log |
| Master Pyramid (6) | Log | Linear | Linear |
| Distortion & Cake (25) | Linear | Linear | Log |
| German Cake (41) | Log | Linear | Log |
| Pumpkin Cake (29), Guru Cake (31) | Linear | Linear | Log |
| Mini Pumpkin (30), Little Pumpkin (43) | no knob (fuzz fixed at maximum) | no knob (Tone Switch, §3) | Log |

Where a schematic does not state the taper (Jumbo Goat, Jumbo Open Goat, Dirty Goat, Jumbo Colour Cake,
Inflatable Johnson), the pots are linear.

### Tone direction

On every model with a Tone knob, clockwise is brighter, with one exception: on **Dreamy Cake** the
tone pot is wired the other way round on the original circuit, and the plugin keeps that: clockwise is
darker. (On Jumbo Colour Cake the plugin turns the pot round so that clockwise is brighter, as on the
original pedal.)

## 3. The switches

Each switch acts only on its own models. On every other model it does nothing: those models stay
bit-identical when the switch moves.

| switch | models | what it does |
|---|---|---|
| **Tone Bypass** (`tone_bypass`) | Raw Cake, Reversed Cake, Pumpkin Cake, Guru Cake | Takes the tone stack out, as the switch on those circuits does: the fuzz reaches the volume pot flat, and the Tone knob stops acting. Louder (+3 to +5 dB on a broadband test signal) and with more mids. |
| **Tone Switch** (`tone_switch`) | Mini Pumpkin (op-amp), Little Pumpkin (transistor) | These circuits have no tone knob, only a two-way switch at the two ends of the tone network: **Dark** or **Bright**. Bright against Dark: above 2 kHz +10 / +7.6 dB, below 500 Hz −1.7 / −6.9 dB. The Tone and Sustain knobs do not act on these two models; Volume does. |
| **HQ** (`hq`) | all | See §4. |
| **Footswitch** (`enabled`) | all | See §5. |

The editor shows only the controls the selected model has: three knobs; three knobs and the Tone
Bypass toggle; or one Volume knob and the Tone Switch toggle.

## 4. Oversampling, HQ, latency and sample rates

The circuits run at a multiple of the host rate. With **HQ off**, each model runs at its own default
factor:

| models | factor, HQ off |
|---|---|
| every transistor model, and the Lost Pumpkin | 2x |
| Pumpkin Cake, Mini Pumpkin, Guru Cake | 4x (their op-amp model needs it to stay stable) |

How the defaults were chosen: on a real guitar DI, each model at each factor against the same model at 384 kHz,
the factor where the worst frame of aliasing stays below the masking threshold (worst-frame ANMR below 0 dB). That measures
the rate, not the model; the model against its circuit is the SPICE check of §7.

The Jumbo Goat at its 2x default rolls the top octave off slightly: its 17th and 19th harmonics of a 1 kHz note
(18 and 20 kHz) come out up to 1.3 and 1.8 dB lower than in the circuit. HQ on brings them back.

With **HQ on**, every model runs at 4x. Switching HQ, like changing model, restarts the circuit from its
resting point. It is not crossfaded, so it can click if you switch while playing.

**Latency.** The plugin reports **31 samples** at the host rate (0.65 ms at 48 kHz) on its `latency`
port. The figure is the same for every model, every factor, HQ on or off, and with the footswitch
engaged or bypassed, so the host's delay compensation never has to change.

**Sample rates.** The plugin runs at the host's rate; no rate is fixed. At
44.1 and 48 kHz the resampling filter keeps the passband up to 20 kHz.

**Changing model at 96 or 192 kHz with very small buffers.** Selecting another model (or switching HQ)
settles the new circuit at its operating point before it plays, and that happens inside one audio block.
At 96 or 192 kHz with a 32-sample buffer, that block can take longer than the block lasts, and the host
may drop it: a short dropout when the model changes, which some hosts report as an xrun. Playback after
the change is unaffected, and so is the CPU use while playing. On an i7-8700K, 32-sample blocks,
HQ on: the model-change block takes 150-300 µs at 192 kHz (the block lasts 167 µs) and about 370 µs at
96 kHz (333 µs). At 48 kHz it fits (about 330 of 667 µs). At 44.1 or 48 kHz, or with buffers of 128
samples or more, it does not happen.

## 5. The footswitch (bypass)

The `enabled` port carries the LV2 bypass designation, so a host that has its own bypass button drives
this port.

- **Bypassed**, the output is the input **delayed by the same 31 samples** the plugin reports, so
  bypassing does not shift the track in time against others, and parallel chains do not comb-filter.
- **Engaging and bypassing crossfade** over about 10.7 ms instead of switching hard.
- **While bypassed, the circuit stops** (no CPU) and keeps its state, like a real pedal with its signal
  cut. On re-engage the circuit first runs for 10 ms with the output still dry, and then the crossfade
  starts, so the circuit's own settling transient is not heard.

## 6. CPU per model

Share of **one core** at 48 kHz, HQ off (each model at its default factor), on **one machine**: an
Intel i7-8700K. Each cell is the worse of two settings: the
default knobs (Sustain 0.8) and Sustain at 1.0. The Jumbo Goat row is at its default 2x.
Other CPUs give other figures, so read these as relative costs.

| port | model | % of one core |
|---:|---|---:|
| 0 | The Pyramid | 3.84 |
| 1 | Gain Pyramid | 3.77 |
| 2 | Scooped Pyramid | 3.70 |
| 3 | Fat Pyramid | 3.76 |
| 4 | Green Pyramid | 3.78 |
| 5 | Vintage Pyramid | 3.65 |
| 6 | Master Pyramid | 4.03 |
| 7 | Bass Pyramid | 3.76 |
| 8 | Inverted Pyramid | 3.81 |
| 9 | Super Pyramid | 3.74 |
| 10 | Super New Pyramid | 3.75 |
| 11 | Alpha Goat | 3.65 |
| 12 | Beta Violet Goat | 3.83 |
| 13 | Pink Goat | 3.77 |
| 14 | Violet Goat | 3.83 |
| 15 | Blue Goat | 3.80 |
| 16 | White Goat | 3.66 |
| 17 | Singer Goat | 3.65 |
| 18 | Boxed Goat | 3.74 |
| 19 | Beaver Goat | 3.67 |
| 20 | Coloured Goat | 3.87 |
| 21 | Jumbo Goat | 5.24 |
| 22 | Jumbo Open Goat | 3.67 |
| 23 | Goat&S OD | 4.02 |
| 24 | Dirty Goat | 3.83 |
| 25 | Distortion & Cake | 3.95 |
| 26 | Raw Cake | 3.80 |
| 27 | Jumbo Colour Cake | 2.74 |
| 28 | Classic Cake | 3.72 |
| 29 | Pumpkin Cake | 2.43 |
| 30 | Mini Pumpkin | 2.35 |
| 31 | Guru Cake | 2.43 |
| 32 | Reversed Cake | 3.87 |
| 33 | Red Army | 3.63 |
| 34 | Crimson Army | 3.61 |
| 35 | Cake War | 3.59 |
| 36 | Black Cake | 3.62 |
| 37 | NYC Cake | 3.69 |
| 38 | Brooklyn Cake | 3.80 |
| 39 | Inflatable Johnson | 3.74 |
| 40 | Precision Cake | 3.87 |
| 41 | German Cake | 3.62 |
| 42 | Dreamy Cake | 3.32 |
| 43 | Little Pumpkin | 3.71 |
| 44 | Lost Pumpkin | 4.04 |

HQ on runs every model at 4x and costs more. HQ-on figures are not published for this
release. Only the selected model computes; the other 44 cost nothing while unselected, and a bypassed
plugin costs nothing.

## 7. How close each model is to its circuit

Every model is checked against a SPICE simulation of its schematic, on plucked guitar notes over a grid
of knob settings and input levels. The score is ANMR, the noise-to-mask ratio of the difference (A-weighted)
in the worst frame. The project's criterion is below −10 dB. At the default factor (HQ off), some
settings miss it. Every one of them has Sustain or Tone fully clockwise. The worst cell per model:

| model | worst setting that misses (Sustain / Tone / input) | HQ off | with HQ on |
|---|---|---:|---:|
| Master Pyramid | 1.0 / 1.0 / 0.3 V | −1.23 dB | −12.91 dB |
| Distortion & Cake | 1.0 / 1.0 / 0.3 V | −3.98 dB | −14.20 dB |
| Dirty Goat | 1.0 / 1.0 / 0.3 V | −5.84 dB | −14.03 dB |
| Raw Cake | 1.0 / 1.0 / 0.3 V | −6.29 dB | −17.05 dB |
| Reversed Cake | 1.0 / 1.0 / 0.3 V | −6.36 dB | −14.34 dB |
| Precision Cake | 1.0 / 1.0 / 0.3 V | −7.99 dB | −15.29 dB |
| Dreamy Cake | 1.0 / 1.0 / 0.3 V | −8.17 dB | −19.90 dB |
| Jumbo Goat | 1.0 / 1.0 / 0.3 V | −7.03 dB | −15.83 dB |
| Jumbo Open Goat | 1.0 / 1.0 / 0.3 V | −9.29 dB | −18.00 dB |
| Lost Pumpkin | 1.0 / 1.0 / 0.1 V | −9.77 dB | −18.55 dB (whole grid) |
| shared four-transistor engine (measured on its reference circuit) | 1.0 / 1.0 / 0.3 V | −9.83 dB | passes |

If you play with Sustain and Tone near maximum and want the closest match to the circuit, turn
**HQ on**: every setting in the table then passes.

One case stays outside the criterion with HQ on as well: the **Jumbo Goat** with Sustain low (0.3-0.6),
Tone at noon and a quiet input (0.1 V), where its noise gate is almost shut and the output is close to
silence (−7.9 to −9.5 dB). The difference is about −133 dBV in absolute terms, and the score is high
only because it is measured relative to near-silence.
