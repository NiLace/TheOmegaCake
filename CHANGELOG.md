# Changelog

## 1.0.2 — 2026-10-19

LV2 version pair in `lv2/theomegacake.ttl`: `lv2:minorVersion 1`, `lv2:microVersion 2` (an odd micro marks a
development build in LV2, so there is no 1.0.1).

### Changed

- **Brooklyn Cake** (model 38): its Volume pot is now **logarithmic** (500k log), as the written sources give the
  board it models. With Volume at full nothing changes; at part travel it is quieter than in 1.0.0 (about 6 dB at
  noon), as a log pot is. Sustain and Tone stay linear. Sessions that use model 38 with Volume below full will sound
  quieter; raise Volume to match.

### Fixed

- **Lost Pumpkin** (model 44): selecting it, or switching **HQ** while it is selected, no longer takes about
  1.3 ms inside one audio block (longer than a 32-sample block at any rate, so a host with a small buffer could
  drop out). Its rest point does not depend on the rate or the knobs and is now worked out once, when the plugin
  is loaded: the switch now costs about 60-80 us on an i7-8700K, below every other model. The sound is unchanged
  (at most one least significant bit of a 24-bit sample).

- The editor now asks the host for the **HQ**, **Tone Bypass** and **Tone Switch** ports
  (`ui:portNotification` in `lv2/theomegacake_ui.ttl`). A host that sends the editor only the ports it lists
  could show those three switches in a stale position after automation or a preset change.

## 1.0.0 — first public version

LV2 version pair in `lv2/theomegacake.ttl`: `lv2:minorVersion 1`, `lv2:microVersion 0`.

### What it is

- An LV2 fuzz pedal plugin (Linux, x86-64) with **45 models** (`model` port 0-44) of one family of
  fuzz circuits: the four-stage transistor originals from 1969 on, their op-amp versions, and the
  clones and mods that followed, up to a production circuit of 2025. See `docs/MODEL_LIST.md`.
- Each model is a **circuit simulation**: transistors, diodes and op-amps solved together, sample by
  sample, with the component values of its schematic. Each one is checked against a SPICE simulation of
  the same circuit.
- **Sustain, Tone, Volume**, with each model's own pot tapers from its schematic (`docs/USER_GUIDE.md`
  §2), and a ~20 ms knob glide so automation does not click.
- The switches the circuits had: **Tone Bypass** (Raw Cake, Reversed Cake, Pumpkin Cake, Guru Cake) and a
  two-way **Tone Switch**, Dark/Bright (Mini Pumpkin, Little Pumpkin). Each acts only on its own models.
- **HQ switch**: every model at 4x oversampling. Off, each model runs at its own default (2x; 4x for
  the op-amp models).
- **Constant latency** of 31 samples, reported to the host, the same for every model and setting.
- **Footswitch** on the LV2 `enabled` designation: the bypassed signal is delayed by the reported latency,
  engaging and bypassing crossfade over ~10.7 ms, and a bypassed plugin uses no CPU.
- **Level-matched models**: with Volume at full, every model has the same median loudness over its
  Sustain × Tone travel on a guitar DI, and none peaks above −2 dBTP. The Lost Pumpkin also follows
  the typical model's level as its Tone turns (its own circuit drops about 4 dB more towards bright).
- **Real-time safe**: no memory allocation on the audio thread, including model changes
  (`lv2:hardRTCapable`).
- **Generic x86-64 baseline**, with a faster code path chosen at run time on CPUs that support it.
- **Photoreal X11 editor** that shows only the controls the selected model has.
- **CPU**: 2.4 to 4.0 % of one core per model at 48 kHz, HQ off, except the Jumbo Goat at 5.2 %
  (Intel i7-8700K at rest, one machine; full table in `docs/USER_GUIDE.md` §6).

### Known limitations

- **Not every setting is inaudibly close to the circuit at the default oversampling.** Against SPICE,
  some settings with Sustain or Tone fully clockwise miss the project's audibility criterion (ANMR
  below −10 dB) at HQ off: on nine models with an engine of their own, and on the shared four-transistor
  engine (measured on Classic Cake), the worst being the Master Pyramid (−1.2 dB). With HQ on, all of
  them pass. The list is in `docs/USER_GUIDE.md` §7.
- **Jumbo Goat, gate nearly closed**: with low Sustain and a quiet input, its noise gate leaves the output
  close to silence, and the difference from the circuit (about −133 dBV absolute) misses the criterion
  even with HQ on.
- **Model change and HQ change are not crossfaded**: they restart the circuit and can click if made while
  playing.
- **At 96/192 kHz with 32-sample buffers**, the block in which the model (or HQ) changes can overrun and
  be dropped by the host. Playback afterwards is unaffected.
- **Little Pumpkin**: the trace it is built from does not name the transistor type; the one used is an
  assumption.
- **The editor is X11 only** (it runs under XWayland).
- **CPU figures come from one machine** (i7-8700K). HQ-on costs are not published for this version.
