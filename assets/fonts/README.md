# Bundled fonts — The Omega Cake GUI

These TTFs are bundled into the LV2 UI so the panel renders identically on any machine, independent of
what fonts the host system has installed. They are loaded directly via FreeType into Cairo (see
`src/ui_fonts.h`), bypassing fontconfig.

The panel itself is a set of Blender renders (`gui/photo/`); its
printed labels are part of those images. The UI draws only one line of live text: the display's
readout and the model list.

| File                  | Family / weight                        | Used for                            |
|-----------------------|----------------------------------------|-------------------------------------|
| `Doto-Round-Bold.ttf` | Doto, instanced ROND 100 / wght 700    | the display readout and model list  |
| `SpaceMono-400.ttf`   | Space Mono Regular                     | fallback text face                  |
| `SpaceMono-700.ttf`   | Space Mono Bold                        | fallback text face                  |

## License — IMPORTANT for publishing

Both families are licensed under the **SIL Open Font License 1.1** (`OFL-Doto.txt`,
`OFL-SpaceMono.txt`). The OFL permits bundling them with software under any licence -- including this
plugin's AGPL-3.0-only -- but the font files remain under the OFL and are not relicensed. The `OFL-*.txt`
files are copied into the LV2 bundle next to them (OFL 1.1 §2 requires the licence with every copy).
A static instance of Doto is a Modified Version: it keeps the OFL and does not reuse a reserved name.

Provenance: Space Mono from Google Fonts; Doto from google/fonts
`ofl/doto`, instanced as above.
