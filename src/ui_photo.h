// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// The photoreal panel: pre-rendered Blender images composited with Cairo.
//
// No X in here, so the preview tool draws exactly what the host shows. The images are rendered at
// twice the logical size (HiDPI) and scaled once per scale factor into a cache, never per paint.
//
// Two knob layouts (cfg): THREE (Sustain, Tone, Volume) and ONE (one bigger Volume, the Little
// Pumpkin models). Two states, active and bypass: the print itself changes look with the footswitch
// (livelier colours, the cherries lit, the labels lit), so EVERY layer exists in both states and the
// change crossfades between them.
//
// Decoding the PNGs is most of the cost of opening the editor. Decoding and scaling run on a few
// threads, one image per job: each image is decoded or scaled by exactly the code of the serial
// path, so the pixels do not depend on how the jobs fall on the threads.
//
// Layers, back to front:
//   plate      the whole top for the cfg and state: knobs at rest, NO toggle, HQ off, display empty
//   knobs      one frame of each knob's strip, through a soft circular mask over its base circle
//   toggle     the toggle with its two legends (only the selected one lit), on the models that have
//              one, through a feathered rectangle
//   hq         the lit HQ button, pasted only while HQ is on (and only lit with the pedal on)
//   fsw        one frame of the 3PDT strip (rest -> pressed)
//   readout    the model's name and the caret, drawn live in the dot-matrix face
#pragma once

#include "ui_photo_layout.h"

#include <cairo.h>
#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <string>
#include <thread>
#include <vector>

namespace tocake {
namespace ui {
namespace photo {

inline constexpr double LOGICAL_W = ASSET_W / 2.0;   // the window at scale 1
inline constexpr double LOGICAL_H = ASSET_H / 2.0;
inline constexpr int CFG_THREE = 0, CFG_ONE = 1;
inline constexpr const char* CFG_NAME[2] = { "three", "one" };
inline constexpr int SET_BYPASS = 0, SET_TONE_SWITCH = 1;        // the toggle's legends: BYPASS/TONE, BRIGHT/DARK
inline constexpr const char* SET_NAME[2] = { "bypass", "tone_switch" };
inline constexpr int ST_ACTIVE = 0, ST_BYPASS = 1;
inline constexpr const char* ST_NAME[2] = { "active", "bypass" };
inline constexpr unsigned BACKDROP_HEX = 0x161618;    // around the pedal: the colour its renders fade to
// The lit readout as the scene renders it (median of its lit pixels on the plate render).
inline constexpr unsigned READOUT_HEX = 0xD45EFF;
inline constexpr double KNOB_MIN_DEG = -135.0, KNOB_MAX_DEG = 135.0;
inline constexpr double FADE_S = 0.120;               // crossfade on the footswitch
inline constexpr double FSW_S = 0.080;                // plunger travel, down or up
inline constexpr double FEATHER = 3.0;                // asset px of soft edge on a pasted region

// Which layout and which toggle legends a model takes, from its controls code
// (VARIANT_CONTROLS: 0 three knobs, 1 three knobs + Tone Bypass, 2 Volume + Tone Switch).
inline int cfg_of(int controls) { return controls == 2 ? CFG_ONE : CFG_THREE; }
inline int set_of(int controls) { return controls == 2 ? SET_TONE_SWITCH : SET_BYPASS; }
inline int n_knobs(int cfg) { return cfg == CFG_ONE ? 1 : 3; }
inline const KnobGeo& knob_geo(int cfg, int k) { return cfg == CFG_ONE ? KNOBS_ONE[0] : KNOBS_THREE[k]; }
inline const IBox& toggle_box(int set) { return set == SET_TONE_SWITCH ? TOGGLE_TONE_SWITCH : TOGGLE_BYPASS; }

// --- a few threads over a list of independent jobs ----------------------------------

inline constexpr int MAX_THREADS = 6;

// Runs fn(0..n-1), each job once, on up to MAX_THREADS threads (the caller's among them), and
// returns when all are done. The jobs are taken in list order, so a list puts its longest first.
// A thread that cannot be started leaves its share to the others.
template <class Fn>
inline void run_jobs(int n, Fn&& fn)
{
    std::atomic<int> next{ 0 };
    auto worker = [&] { for (int i; (i = next.fetch_add(1)) < n;) fn(i); };
    const unsigned hw = std::thread::hardware_concurrency();
    const int nt = std::min({ int(hw ? hw : 1), MAX_THREADS, n });
    std::vector<std::thread> t;
    for (int i = 0; i < nt - 1; ++i) {
        try { t.emplace_back(worker); } catch (...) { break; }
    }
    worker();
    for (auto& th : t) th.join();
}

// --- the images, as loaded from the bundle --------------------------------------

// The images of one cfg. Only one cfg is ever on screen; the other is read when the model changes
// layout. The images both cfgs show (toggles, footswitch, HQ) may be shared between two sets: each
// set holds its own reference to them (cairo counts its surfaces' references).
struct Assets {
    int cfg = -1;
    cairo_surface_t* plate[2] = {};        // [state]
    cairo_surface_t* knob[3][2] = {};      // [knob][state], frames stacked vertically
    cairo_surface_t* toggle[2][2] = {};    // [set][state], two frames: down, up
    cairo_surface_t* hq_on = nullptr;
    cairo_surface_t* fsw[2] = {};          // [state]
    bool ok = false;
};

// PNG decodes since the process started: tells a decode from a cache hit.
inline std::atomic<unsigned> g_png_decodes{ 0 };

inline cairo_surface_t* load_png(const std::string& path)
{
    g_png_decodes.fetch_add(1);
    cairo_surface_t* s = cairo_image_surface_create_from_png(path.c_str());
    if (cairo_surface_status(s) != CAIRO_STATUS_SUCCESS) { cairo_surface_destroy(s); return nullptr; }
    return s;
}

inline void free_assets(Assets& a)
{
    auto d = [](cairo_surface_t*& s) { if (s) cairo_surface_destroy(s); s = nullptr; };
    for (auto& s : a.plate) d(s);
    for (auto& k : a.knob) for (auto& s : k) d(s);
    for (auto& t : a.toggle) for (auto& s : t) d(s);
    d(a.hq_on);
    for (auto& s : a.fsw) d(s);
    a = Assets{};
}

inline bool sized(cairo_surface_t* s, int w, int h)
{
    return s && cairo_image_surface_get_width(s) == w && cairo_image_surface_get_height(s) == h;
}

// One image to decode: its file, where it goes, and the size it must have.
struct LoadJob {
    std::string name;
    cairo_surface_t** dst;
    int w, h;
    bool both_cfgs;     // shown by both cfgs: may be taken from the other cfg's set
};

// Every image of the cfg, into `a`'s slots. The knob strips are most of the bytes, so they come
// first, then the plates, then the small parts.
inline std::vector<LoadJob> load_list(Assets& a, int cfg)
{
    std::vector<LoadJob> l;
    const std::string c = CFG_NAME[cfg];
    for (int st = 0; st < 2; ++st)
        for (int k = 0; k < n_knobs(cfg); ++k) {
            const KnobGeo& kg = knob_geo(cfg, k);
            l.push_back({ "knob_" + c + "_" + ST_NAME[st] + "_" + kg.name + ".png", &a.knob[k][st],
                          kg.x1 - kg.x0, (kg.y1 - kg.y0) * kg.frames, false });
        }
    for (int st = 0; st < 2; ++st)
        l.push_back({ "plate_" + c + "_" + ST_NAME[st] + ".png", &a.plate[st], int(ASSET_W), int(ASSET_H), false });
    for (int st = 0; st < 2; ++st)
        l.push_back({ std::string("fsw_") + ST_NAME[st] + ".png", &a.fsw[st], FSW.x1 - FSW.x0,
                      (FSW.y1 - FSW.y0) * FSW_FRAMES, true });
    for (int st = 0; st < 2; ++st)
        for (int t = 0; t < 2; ++t) {
            const IBox& b = toggle_box(t);
            l.push_back({ std::string("toggle_") + SET_NAME[t] + "_" + ST_NAME[st] + ".png", &a.toggle[t][st],
                          b.x1 - b.x0, (b.y1 - b.y0) * 2, true });
        }
    l.push_back({ "hq_on.png", &a.hq_on, HQ_LIT.x1 - HQ_LIT.x0, HQ_LIT.y1 - HQ_LIT.y0, true });
    return l;
}

// The name an image is looked up by in the other cfg's set.
inline std::string share_key(const LoadJob& j, int /* cfg: the same name in both */)
{
    return j.name;
}

// `dir` is the bundle's `photo/` directory. A missing or unreadable image leaves `ok` false and the
// UI draws a plain placeholder: a panel with a hole in it would look like the plugin.
// `other`, when given and complete, is the other cfg's set: the images both cfgs show are taken
// from it instead of decoded again. `threaded` false decodes on the calling thread alone.
inline Assets load_assets(const std::string& dir, int cfg, const Assets* other = nullptr, bool threaded = true)
{
    Assets a;
    a.cfg = cfg;
    std::vector<LoadJob> l = load_list(a, cfg);
    std::vector<const LoadJob*> todo;
    std::vector<LoadJob> theirs;
    Assets other_slots;                     // the other set's pointers, looked up through its own list
    if (other && other->ok && other->cfg != cfg) {
        other_slots = *other;
        theirs = load_list(other_slots, other->cfg);
    }
    for (const LoadJob& j : l) {
        cairo_surface_t* found = nullptr;
        if (j.both_cfgs)
            for (const LoadJob& o : theirs)
                if (share_key(o, other->cfg) == share_key(j, cfg) && o.both_cfgs) found = *o.dst;
        if (found) *j.dst = cairo_surface_reference(found);
        else todo.push_back(&j);
    }
    auto job = [&](int i) { *todo[size_t(i)]->dst = load_png(dir + "/" + todo[size_t(i)]->name); };
    if (threaded) run_jobs(int(todo.size()), job);
    else for (int i = 0; i < int(todo.size()); ++i) job(i);
    bool ok = true;
    for (const LoadJob& j : l) ok &= sized(*j.dst, j.w, j.h);
    a.ok = ok;
    return a;
}

// --- the images, scaled to the window ------------------------------------------

// A pasted region at window scale: its frames, where they sit, and the soft mask they go through.
struct Region {
    std::vector<cairo_surface_t*> frame;
    cairo_surface_t* mask = nullptr;       // A8, feathered edges (rectangles) -- null for knobs
    int x = 0, y = 0;
};

// `f` = window pixels per asset pixel (0.5 at scale 1), for the cfg the assets hold.
struct Scaled {
    double f = -1.0;
    int    cfg = -1;
    cairo_surface_t* plate[2] = {};
    Region knob[3][2], toggle[2][2], hq, fsw[2];
    // the readout, composed once per model and scale (a text path stroked and filled is the most
    // expensive layer of a frame, and it only changes with the model)
    mutable std::string readout_key;
    mutable cairo_surface_t* readout = nullptr;
};

inline void free_region(Region& r)
{
    for (auto* p : r.frame) cairo_surface_destroy(p);
    if (r.mask) cairo_surface_destroy(r.mask);
    r = Region{};
}

inline void free_scaled(Scaled& s)
{
    for (auto*& p : s.plate) if (p) cairo_surface_destroy(p);
    for (auto& k : s.knob) for (auto& r : k) free_region(r);
    for (auto& t : s.toggle) for (auto& r : t) free_region(r);
    free_region(s.hq);
    for (auto& r : s.fsw) free_region(r);
    if (s.readout) cairo_surface_destroy(s.readout);
    s = Scaled{};
}

// One region of a source image, resampled onto its own surface whose top-left sits at the integer
// window pixel (px0, py0). Each window pixel samples the source at the same place the whole-plate
// resample does, so a pasted frame lines up with the plate under it.
inline cairo_surface_t* resample(cairo_surface_t* src, double sx, double sy, double sw, double sh,
                                 double ax, double ay, double f, int px0, int py0, int pw, int ph)
{
    cairo_surface_t* sub = cairo_surface_create_for_rectangle(src, sx, sy, sw, sh);
    cairo_surface_t* dst = cairo_image_surface_create(CAIRO_FORMAT_RGB24, pw, ph);
    cairo_t* cr = cairo_create(dst);
    cairo_translate(cr, -px0, -py0);
    cairo_scale(cr, f, f);
    cairo_set_source_surface(cr, sub, ax, ay);
    cairo_pattern_set_filter(cairo_get_source(cr), f == 1.0 ? CAIRO_FILTER_NEAREST : CAIRO_FILTER_GOOD);
    cairo_pattern_set_extend(cairo_get_source(cr), CAIRO_EXTEND_PAD);
    cairo_paint(cr);
    cairo_destroy(cr);
    cairo_surface_destroy(sub);
    cairo_surface_flush(dst);
    return dst;
}

// A feathered rectangle in window pixels: opaque inside, ramping to 0 over FEATHER asset px at the
// edges -- the outer ring is where a region and its plate show the same thing, so it may blend.
inline cairo_surface_t* feather_mask(int pw, int ph, double fpx)
{
    cairo_surface_t* m = cairo_image_surface_create(CAIRO_FORMAT_A8, pw, ph);
    cairo_surface_flush(m);
    unsigned char* d = cairo_image_surface_get_data(m);
    const int stride = cairo_image_surface_get_stride(m);
    for (int y = 0; y < ph; ++y)
        for (int x = 0; x < pw; ++x) {
            const double e = std::min({ x + 0.5, pw - x - 0.5, y + 0.5, ph - y - 0.5 });
            const double t = std::clamp(e / std::max(fpx, 1.0), 0.0, 1.0);
            d[y * stride + x] = (unsigned char)std::lround(255.0 * t * t * (3.0 - 2.0 * t));
        }
    cairo_surface_mark_dirty(m);
    return m;
}

inline Region scale_region(cairo_surface_t* src, int x0, int y0, int x1, int y1, int frames, double f,
                           bool mask)
{
    Region r;
    const int w = x1 - x0, h = y1 - y0;
    r.x = int(std::floor(x0 * f)); r.y = int(std::floor(y0 * f));
    const int pw = int(std::ceil(x1 * f)) - r.x, ph = int(std::ceil(y1 * f)) - r.y;
    for (int i = 0; i < frames; ++i)
        r.frame.push_back(resample(src, 0, double(i * h), w, h, x0, y0, f, r.x, r.y, pw, ph));
    if (mask) r.mask = feather_mask(pw, ph, FEATHER * f);
    return r;
}

// The images at window scale `f`. One job per image (a whole plate, or one region with its frames),
// longest first: each reads its own source image and writes its own surfaces, so no cairo object is
// shared between two threads. `threaded` false scales on the calling thread alone.
inline void ensure_scaled(Scaled& s, const Assets& a, double f, bool threaded = true)
{
    if (!a.ok || (s.f == f && s.cfg == a.cfg)) return;
    free_scaled(s);
    s.f = f; s.cfg = a.cfg;
    const int W = int(std::ceil(ASSET_W * f)), H = int(std::ceil(ASSET_H * f));
    struct Job { cairo_surface_t* src; cairo_surface_t** plate; Region* region; IBox box; int frames; bool mask; };
    std::vector<Job> l;
    for (int st = 0; st < 2; ++st)
        for (int k = 0; k < n_knobs(a.cfg); ++k) {
            const KnobGeo& kg = knob_geo(a.cfg, k);
            l.push_back({ a.knob[k][st], nullptr, &s.knob[k][st], { kg.x0, kg.y0, kg.x1, kg.y1 }, kg.frames, false });
        }
    for (int st = 0; st < 2; ++st) l.push_back({ a.plate[st], &s.plate[st], nullptr, {}, 0, false });
    for (int st = 0; st < 2; ++st) l.push_back({ a.fsw[st], nullptr, &s.fsw[st], FSW, FSW_FRAMES, true });
    for (int st = 0; st < 2; ++st)
        for (int t = 0; t < 2; ++t) l.push_back({ a.toggle[t][st], nullptr, &s.toggle[t][st], toggle_box(t), 2, true });
    l.push_back({ a.hq_on, nullptr, &s.hq, HQ_LIT, 1, true });
    auto job = [&](int i) {
        const Job& j = l[size_t(i)];
        if (j.plate) *j.plate = resample(j.src, 0, 0, ASSET_W, ASSET_H, 0, 0, f, 0, 0, W, H);
        else *j.region = scale_region(j.src, j.box.x0, j.box.y0, j.box.x1, j.box.y1, j.frames, f, j.mask);
    };
    if (threaded) run_jobs(int(l.size()), job);
    else for (int i = 0; i < int(l.size()); ++i) job(i);
}

// --- what is drawn ---------------------------------------------------------------

struct View {
    double sustain = 0.8, tone = 0.5, volume = 0.6;   // knob positions, 0..1
    int    controls = 0;      // VARIANT_CONTROLS of the model shown
    double lit = 1.0;         // 1 active, 0 bypass; in between while crossfading
    double fsw = 0.0;         // 0 rest .. 1 pressed
    bool   toggle_up = false; // the toggle's position (its port > 0.5)
    bool   hq = false;        // the `hq` port
    bool   display_lit = true;   // the readout shows (active, hovered or the list open)
    const char* variant = "";
    cairo_font_face_t* face = nullptr;   // the dot-matrix face; null falls back to monospace
};

inline int knob_frame(const KnobGeo& kg, double v01)
{
    const double deg = KNOB_MIN_DEG + (KNOB_MAX_DEG - KNOB_MIN_DEG) * std::clamp(v01, 0.0, 1.0);
    const double t = (deg - kg.deg_first) / (kg.deg_last - kg.deg_first);
    return std::clamp(int(std::lround(t * (kg.frames - 1))), 0, kg.frames - 1);
}

inline void set_hex(cairo_t* cr, unsigned hex, double alpha)
{
    cairo_set_source_rgba(cr, ((hex >> 16) & 0xff) / 255.0, ((hex >> 8) & 0xff) / 255.0,
                          (hex & 0xff) / 255.0, alpha);
}

// The readout: a backlit dot-matrix line in the colour the scene renders it, with a faint glow.
// The dots are stroked as well as filled: at scale 1 they are thinner than a pixel and a fill alone
// reads dark (the Valve Howler's finding). Font size calibrated on the scene's own render of the
// plate: ink 22 asset px tall on a 48 mm glass.
inline constexpr double TEXT_MM = 3.82;
inline constexpr double GLASS_MM = 48.0;
inline constexpr double DOT_STROKE_PX = 0.9;

inline void draw_readout(cairo_t* cr, const View& v, double f)
{
    const Box& d = DISPLAY;
    const double w = (d.x1 - d.x0) * f;
    const double mm = w / GLASS_MM;
    const double cx = (d.x0 + d.x1) / 2.0 * f, cy = (d.y0 + d.y1) / 2.0 * f;
    std::string label = v.variant;
    for (char& c : label) c = char(std::toupper(static_cast<unsigned char>(c)));

    cairo_save(cr);
    if (v.face) cairo_set_font_face(cr, v.face);
    else cairo_select_font_face(cr, "monospace", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
    cairo_set_font_size(cr, TEXT_MM * mm);
    cairo_text_extents_t e;
    cairo_text_extents(cr, label.c_str(), &e);
    const double max_w = w - 12.0 * mm;
    if (e.width > max_w) {
        cairo_set_font_size(cr, TEXT_MM * mm * max_w / e.width);
        cairo_text_extents(cr, label.c_str(), &e);
    }
    const double tx = cx - 2.0 * mm - e.width / 2.0 - e.x_bearing;
    const double ty = cy - e.height / 2.0 - e.y_bearing;
    cairo_move_to(cr, tx, ty);
    cairo_text_path(cr, label.c_str());
    const double kx = cx + (w / 2.0) - 3.2 * mm;     // caret at the right end of the window
    cairo_move_to(cr, kx - 1.3 * mm, cy - 0.7 * mm);
    cairo_line_to(cr, kx + 1.3 * mm, cy - 0.7 * mm);
    cairo_line_to(cr, kx, cy + 0.8 * mm);
    cairo_close_path(cr);
    cairo_set_line_join(cr, CAIRO_LINE_JOIN_ROUND);
    set_hex(cr, READOUT_HEX, 0.07);
    cairo_set_line_width(cr, 0.6 * mm);
    cairo_stroke_preserve(cr);
    set_hex(cr, READOUT_HEX, 1.0);
    cairo_set_line_width(cr, DOT_STROKE_PX * f * 2.0);    // f * 2 = the UI scale
    cairo_stroke_preserve(cr);
    cairo_fill(cr);
    cairo_restore(cr);
}

// A region's frame through its feathered mask, at `alpha`.
inline void paste(cairo_t* cr, const Region& r, int frame, double alpha)
{
    if (alpha <= 0.0 || r.frame.empty()) return;
    cairo_surface_t* fr = r.frame[size_t(frame)];
    cairo_save(cr);
    // clip to the region first: a group is as large as the clip, and an unclipped one is the whole
    // window -- allocated and composited again for every pasted part, every frame
    cairo_rectangle(cr, r.x, r.y, cairo_image_surface_get_width(fr), cairo_image_surface_get_height(fr));
    cairo_clip(cr);
    if (alpha >= 1.0) {
        cairo_set_source_surface(cr, fr, r.x, r.y);
        if (r.mask) cairo_mask_surface(cr, r.mask, r.x, r.y);
        else cairo_paint(cr);
    } else {
        cairo_push_group(cr);
        cairo_set_source_surface(cr, fr, r.x, r.y);
        if (r.mask) cairo_mask_surface(cr, r.mask, r.x, r.y);
        else cairo_paint(cr);
        cairo_pop_group_to_source(cr);
        cairo_paint_with_alpha(cr, alpha);
    }
    cairo_restore(cr);
}

// A knob's frame through a soft circular mask over its base circle, at `alpha`.
inline void paste_knob(cairo_t* cr, const Region& r, const KnobGeo& kg, int frame, double f, double alpha)
{
    if (alpha <= 0.0 || r.frame.empty()) return;
    cairo_surface_t* fr = r.frame[size_t(frame)];
    cairo_save(cr);
    cairo_rectangle(cr, r.x, r.y, cairo_image_surface_get_width(fr), cairo_image_surface_get_height(fr));
    cairo_clip(cr);                      // the mask is a pattern over the whole window otherwise
    cairo_set_source_surface(cr, fr, r.x, r.y);
    const double r_in = (kg.r + 1.0) * f, r_out = (kg.r + 3.0) * f;
    cairo_pattern_t* m = cairo_pattern_create_radial(kg.cx * f, kg.cy * f, 0, kg.cx * f, kg.cy * f, r_out);
    const double a = std::min(alpha, 1.0);
    cairo_pattern_add_color_stop_rgba(m, 0.0, 0, 0, 0, a);
    cairo_pattern_add_color_stop_rgba(m, r_in / r_out, 0, 0, 0, a);
    cairo_pattern_add_color_stop_rgba(m, 1.0, 0, 0, 0, 0);
    cairo_mask(cr, m);
    cairo_pattern_destroy(m);
    cairo_restore(cr);
}

// The readout through its cache: the glass's own pixels from the current state's plate with the text
// drawn on them, composed the first time a model is shown in that state at this scale, then pasted.
// Drawn on the plate's pixels (not on transparency) so the result is the same, bit for bit, as drawing
// it directly: antialiased dot edges composited twice through 8-bit alpha differ by up to 27 levels.
// While crossfading (`both`) the background changes every frame: drawn directly.
inline void paint_readout(cairo_t* cr, const Scaled& s, const View& v, bool both, int state)
{
    if (both) { draw_readout(cr, v, s.f); return; }
    const std::string key = std::string(v.variant) + (v.face ? "#d" : "#m") + char('0' + state);
    const Box& d = DISPLAY;
    const int x0 = int(std::floor(d.x0 * s.f)), y0 = int(std::floor(d.y0 * s.f));
    if (s.readout_key != key || !s.readout) {
        if (s.readout) cairo_surface_destroy(s.readout);
        const int w = int(std::ceil(d.x1 * s.f)) - x0, h = int(std::ceil(d.y1 * s.f)) - y0;
        s.readout = cairo_image_surface_create(CAIRO_FORMAT_RGB24, w, h);
        cairo_t* c = cairo_create(s.readout);
        cairo_translate(c, -x0, -y0);
        cairo_set_source_surface(c, s.plate[state], 0, 0);   // nothing else covers the glass
        cairo_paint(c);
        draw_readout(c, v, s.f);
        cairo_destroy(c);
        cairo_surface_flush(s.readout);
        s.readout_key = key;
    }
    cairo_set_source_surface(cr, s.readout, x0, y0);
    cairo_paint(cr);
}

// The whole panel, at `f` window pixels per asset pixel, onto `cr` in window pixels.
inline void draw_panel(cairo_t* cr, const Scaled& s, const View& v)
{
    if (s.f <= 0.0) {                       // assets missing: a flat placeholder, never a half panel
        set_hex(cr, BACKDROP_HEX, 1.0);
        cairo_paint(cr);
        return;
    }
    const double f = s.f, lit = std::clamp(v.lit, 0.0, 1.0);
    // Each layer exists in both states. Outside the crossfade only the current state's image is
    // painted: painting the other underneath at full cover doubles the cost of every frame for nothing
    // (measured: two whole plates per frame where one does).
    const bool both = lit > 0.0 && lit < 1.0;
    const int  only = lit >= 1.0 ? ST_ACTIVE : ST_BYPASS;
    cairo_save(cr);
    cairo_set_source_surface(cr, s.plate[both ? ST_BYPASS : only], 0, 0);
    cairo_paint(cr);
    if (both) {
        cairo_set_source_surface(cr, s.plate[ST_ACTIVE], 0, 0);
        cairo_paint_with_alpha(cr, lit);
    }
    const double values3[3] = { v.sustain, v.tone, v.volume };
    for (int k = 0; k < n_knobs(s.cfg); ++k) {
        const KnobGeo& kg = knob_geo(s.cfg, k);
        const int fr = knob_frame(kg, s.cfg == CFG_ONE ? v.volume : values3[k]);
        paste_knob(cr, s.knob[k][both ? ST_BYPASS : only], kg, fr, f, 1.0);
        if (both) paste_knob(cr, s.knob[k][ST_ACTIVE], kg, fr, f, lit);
    }
    if (v.controls != 0) {
        const int t = set_of(v.controls), fr = v.toggle_up ? 1 : 0;
        paste(cr, s.toggle[t][both ? ST_BYPASS : only], fr, 1.0);
        if (both) paste(cr, s.toggle[t][ST_ACTIVE], fr, lit);
    }
    if (v.hq) paste(cr, s.hq, 0, lit);           // lit HQ only with the pedal on
    const int fi = std::clamp(int(std::lround(v.fsw * (FSW_FRAMES - 1))), 0, FSW_FRAMES - 1);
    if (fi > 0) {                                // at rest both plates already show it
        paste(cr, s.fsw[both ? ST_BYPASS : only], fi, 1.0);
        if (both) paste(cr, s.fsw[ST_ACTIVE], fi, lit);
    }
    cairo_restore(cr);
    if (v.display_lit) paint_readout(cr, s, v, both, only);
}

// Hit tests in ASSET pixels (window pixels / f). Returns the knob index in the cfg's order, or -1.
inline int knob_at(int cfg, double ax, double ay)
{
    for (int k = 0; k < n_knobs(cfg); ++k) {
        const KnobGeo& kg = knob_geo(cfg, k);
        const double dx = ax - kg.cx, dy = ay - kg.cy;
        if (dx * dx + dy * dy <= kg.r * kg.r) return k;
    }
    return -1;
}

inline bool in_box(const Box& b, double ax, double ay)
{
    return ax >= b.x0 && ax <= b.x1 && ay >= b.y0 && ay <= b.y1;
}

}  // namespace photo
}  // namespace ui
}  // namespace tocake
