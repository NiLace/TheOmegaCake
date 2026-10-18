// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// ui_fonts.h — load the bundled TTFs into Cairo font faces via FreeType.
//
// The photoreal panel draws only one line of live text, the display's readout and the model list,
// in Doto (a dot-matrix face, the one the plate is rendered with); Space Mono is the fallback. We do
// NOT rely on the host system having these installed (it usually won't); the
// .ttf files travel inside the LV2 bundle and are loaded here directly, so the
// panel renders identically everywhere. fontconfig is bypassed entirely.
//
// Each cairo_font_face_t owns a PRIVATE FT_Library holding its one FT_Face, hung on the face as user
// data whose destructor is FreeType's own FT_Done_FreeType: when cairo finally drops the face (from
// its process-wide scaled-font cache, possibly long after this UI is gone and the module unloaded),
// FreeType frees the library and the face with it, and nothing calls back into this module once it is
// unloaded.
//
// The destructor is FT_Done_FreeType (FT_Error(*)(FT_Library)) called as cairo's void(*)(void*): a
// conversion the C++ standard leaves undefined, and exact on the platform ABI this plugin targets
// (x86-64 System V: the pointer argument travels in the same register, the int result is discarded).
// The alternative -- a wrapper of the right type -- is code in THIS module, which a host that draws text after
// unloading the module would call into.

#pragma once

#include <cairo/cairo.h>
#include <cairo/cairo-ft.h>
#include <ft2build.h>
#include FT_FREETYPE_H

#include <string>

namespace tocake {
namespace ui {

// The five faces the panel needs. Any may be null if its file was missing;
// the draw code falls back to a Cairo toy font in that case.
struct Fonts {
    cairo_font_face_t* doto    = nullptr;  // the display's dot-matrix readout and the model list
    cairo_font_face_t* mono400 = nullptr;  // fallback text face
    cairo_font_face_t* mono700 = nullptr;
    bool ok = false;                        // true only if all three loaded
};

// Key used to hang the private FT_Library off its cairo font face.
inline cairo_user_data_key_t* ft_library_key() {
    static cairo_user_data_key_t key;
    return &key;
}

// One face from one bundled file, owning its own FreeType library (see the top of this file).
inline cairo_font_face_t* load_face(const std::string& path) {
    FT_Library lib = nullptr;
    if (FT_Init_FreeType(&lib) != 0) return nullptr;
    FT_Face face = nullptr;
    if (FT_New_Face(lib, path.c_str(), 0, &face) != 0) { FT_Done_FreeType(lib); return nullptr; }
    cairo_font_face_t* cf = cairo_ft_font_face_create_for_ft_face(face, 0);
    if (cairo_font_face_status(cf) != CAIRO_STATUS_SUCCESS) {
        cairo_font_face_destroy(cf);
        FT_Done_FreeType(lib);
        return nullptr;
    }
    if (cairo_font_face_set_user_data(cf, ft_library_key(), lib,
                                      reinterpret_cast<cairo_destroy_func_t>(
                                          reinterpret_cast<void (*)()>(&FT_Done_FreeType)))
        != CAIRO_STATUS_SUCCESS) {
        cairo_font_face_destroy(cf);
        FT_Done_FreeType(lib);
        return nullptr;
    }
    return cf;
}

// Load all five faces from `dir` (the directory holding the .ttf files — the
// LV2 bundle path at runtime, or assets/fonts for the standalone preview).
inline Fonts load_fonts(const std::string& dir) {
    Fonts f;
    const std::string base = dir.empty() ? "." : dir;
    f.doto    = load_face(base + "/Doto-Round-Bold.ttf");
    f.mono400 = load_face(base + "/SpaceMono-400.ttf");
    f.mono700 = load_face(base + "/SpaceMono-700.ttf");
    f.ok = f.doto && f.mono400 && f.mono700;
    return f;
}

// Drops this UI's references; each face's FreeType library goes when cairo drops the face.
inline void free_fonts(Fonts& f) {
    if (f.doto)    cairo_font_face_destroy(f.doto);
    if (f.mono400) cairo_font_face_destroy(f.mono400);
    if (f.mono700) cairo_font_face_destroy(f.mono700);
    f = Fonts{};
}

} // namespace ui
} // namespace tocake
