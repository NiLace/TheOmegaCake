# The Omega Cake — model list

One row per model. **Slot** is the value of the `model` port: a saved session stores this number, and
it never changes. **Year** is the year of the circuit. **Year source** summarises what the year rests on.
Unless a row says otherwise, that is
the dating printed on a published, dated trace or schematic of the circuit. **Values from** is where the
component values come from.

Engines:

- **four-transistor**: the shared engine for the classic circuit: input booster, two diode-clipping
  gain stages with the sustain pot between them, a passive tone stack and an output recovery stage.
- **own engine**: a circuit that departs from that topology and has a solver of its own; the
  departure is named in the row.
- **op-amp**: two op-amp gain stages and an op-amp clipper with a diode clamp in its feedback, then
  the passive tone stack.

Component values: "spreadsheet" is the project's component sheet, one tab per model, from which the
values in `src/variants_generated.h` are generated.

| slot | name | engine family | year | year source | values from |
|---:|---|---|---:|---|---|
| 0 | The Pyramid | four-transistor | 1971 | trace dated "circa 1971" | spreadsheet |
| 1 | Gain Pyramid | four-transistor | 1970 | trace dated "late 1969-early 1970" | spreadsheet |
| 2 | Scooped Pyramid | four-transistor | 1970 | trace dated "around 1970-71" | spreadsheet |
| 3 | Fat Pyramid | four-transistor | 1970 | trace dated "around 1970" | spreadsheet |
| 4 | Green Pyramid | four-transistor | 1971 | trace dated "around 1971" | spreadsheet |
| 5 | Vintage Pyramid | four-transistor | 1971 | an early clone, dated "around 1971" | spreadsheet |
| 6 | Master Pyramid | own engine (bright input cap, relocated volume pot, AC-coupled tone stack) | 1973 | maker's timeline, "circa 1973" | spreadsheet |
| 7 | Bass Pyramid | four-transistor | 1972 | trace dated "circa 1972-73", pots dated 1972 | spreadsheet |
| 8 | Inverted Pyramid | four-transistor | 1972 | same caption as the Bass Pyramid | spreadsheet |
| 9 | Super Pyramid | four-transistor | 1973 | trace dated "circa 1973" | spreadsheet |
| 10 | Super New Pyramid | four-transistor | 2022 | traces of January 2022 | spreadsheet |
| 11 | Alpha Goat | four-transistor | 1973 | trace dated "circa 1973-74" | spreadsheet |
| 12 | Beta Violet Goat | four-transistor | 1973 | trace dated "circa 1973-74" | spreadsheet |
| 13 | Pink Goat | four-transistor | 1973 | trace dated "circa 1973" | spreadsheet |
| 14 | Violet Goat | four-transistor | 1973 | trace dated "around 1973-1974" | spreadsheet |
| 15 | Blue Goat | four-transistor | 1974 | trace dated "circa 1974" | spreadsheet |
| 16 | White Goat | four-transistor | 1974 | trace dated "circa 1974-75", pots dated 1974 | spreadsheet |
| 17 | Singer Goat | four-transistor | 1976 | trace dated "circa 1976", pots dated 1975 | spreadsheet |
| 18 | Boxed Goat | four-transistor | 1976 | trace dated "circa 1976-77" | spreadsheet |
| 19 | Beaver Goat | four-transistor | 1976 | trace dated "around 1976" | spreadsheet |
| 20 | Coloured Goat | own engine (first gain stage without clipping diodes, AC-coupled tone stack) | 1973 | clone history, "around 1973" | spreadsheet |
| 21 | Jumbo Goat | own engine (six transistors, three diode pairs, germanium noise gate) | 1973 | maker history, "in 1973" | spreadsheet |
| 22 | Jumbo Open Goat | own engine (the Jumbo Goat without its noise gate) | 1973 | the same production as the Jumbo Goat | spreadsheet |
| 23 | Goat&S OD | own engine (active tone network around the output transistor) | 1974 | first edition, "circa 1974-1976" | spreadsheet |
| 24 | Dirty Goat | own engine (a capacitor in the tone pot's bottom leg) | 1975 | "around 1974 or '75", pots dated 1975-76 | spreadsheet |
| 25 | Distortion & Cake | own engine (AC-coupled tone stack) | 1979 | third edition, "circa 1979-80" | spreadsheet |
| 26 | Raw Cake | own engine (tone-bypass switch) | 1979 | trace dated "circa 1979" | spreadsheet |
| 27 | Jumbo Colour Cake | own engine (three transistors, one diode pair, no output stage) | 1976 | clone history, "around 1976" | spreadsheet |
| 28 | Classic Cake | four-transistor | 1977 | "not available until 1977" (chosen over the pots' 1976) | spreadsheet |
| 29 | Pumpkin Cake | op-amp | 1978 | dated "circa 1978-1980" | the circuit's KiCad netlist |
| 30 | Mini Pumpkin | op-amp (fuzz fixed at maximum, two-way tone switch) | 1978 | maker's timeline, 1978 | the Pumpkin Cake netlist with four values from a published trace (2001) |
| 31 | Guru Cake | op-amp (negative-feedback mod) | 1978 | year of the 1978 op-amp circuit; the mod itself is not dated | the Pumpkin Cake netlist with the mod |
| 32 | Reversed Cake | own engine (tone-bypass switch, one coupling cap per side) | 1979 | trace dated "circa 1979" | spreadsheet |
| 33 | Red Army | four-transistor | 1991 | first edition, "around 1991" | spreadsheet |
| 34 | Crimson Army | four-transistor | 1993 | second edition, "circa 1993-1994" | spreadsheet |
| 35 | Cake War | four-transistor | 1993 | no single production unit; the year of the schematic it follows, chosen | spreadsheet |
| 36 | Black Cake | four-transistor | 1995 | "from around 1995" | spreadsheet |
| 37 | NYC Cake | four-transistor | 2000 | first edition of the reissue, "in 2000": the first reissue board, with Schottky clipping diodes (its next revision differs only in the bypass wiring, which is not modelled) | spreadsheet; pot tapers measured on a unit of this board |
| 38 | Brooklyn Cake | four-transistor | 2008 | fourth edition, "in 2008" | spreadsheet; Volume 500k log as the written sources give this board; Sustain and Tone linear, as measured on the earlier board |
| 39 | Inflatable Johnson | four-transistor | 1995 | "got its start in 1995" | spreadsheet |
| 40 | Precision Cake | own engine (simplified tone stack without the mid notch, output-stage emitter bypass) | 2008 | "in the late 2000s" (2008 chosen) | spreadsheet |
| 41 | German Cake | own engine (as Precision Cake) | 2002 | "circa 2002" | spreadsheet |
| 42 | Dreamy Cake | own engine (no emitter resistors on the first three stages) | 1999 | "hit the market in 1999" | spreadsheet |
| 43 | Little Pumpkin | four-transistor (no sustain pot, two-way tone switch) | 1976 | maker's timeline, 1976 (chosen) | spreadsheet, from a published trace (2003); the transistor type is not named on the trace and is assumed |
| 44 | Lost Pumpkin | own engine (four op-amp halves: three diode-clipping stages and an output gain stage) | 2025 | production began December 2025 (press reports of 11 and 12 December 2025) | two drawings of the production board (2026), agreeing on every audio-path value but one input bias resistor |

"chosen" marks a year where the sources disagree and the project picked one.
