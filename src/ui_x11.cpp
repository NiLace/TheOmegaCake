// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// ui_x11.cpp — LV2 X11 UI for The Omega Cake: the photoreal pedal.
//
// This is the window-and-events layer. All painting lives in `ui_photo.h` (Blender renders composited
// with Cairo); here lives:
//   - a child X11 window inside the host's parent, with the parent's visual and depth,
//   - the bundle's fonts; the panel images, decoded once per process into a cache every instance
//     shares, and scaled from it once per window size and layout,
//   - the X events pumped from the LV2 idle callback: knob drags and the wheel, the footswitch (it
//     sinks while held and toggles `enabled` on release over it), the toggle on the models that have
//     one (`tone_bypass` or `tone_switch`), the HQ button (only while the pedal is on), and the model
//     list under the display,
//   - the host's port changes reflected back (port_event), and `ui:touch` around every gesture so
//     touch-mode automation and MIDI-learn behave.
//
// Which controls a model has comes from VARIANT_CONTROLS (generated from the engine kind the DSP
// dispatches on). The geometry is
// the Blender scene's projection (`ui_photo_layout.h`), never measured off an image.
//
// No Pugl, no SVG renderer: only Xlib + cairo-xlib + FreeType. X11-only by design (runs under XWayland).

#include "ui_photo.h"
#include "ui_fonts.h"
#include "variant_labels_generated.h"

#include <lv2/core/lv2.h>
#include <lv2/ui/ui.h>
#include <lv2/atom/atom.h>
#include <lv2/options/options.h>
#include <lv2/urid/urid.h>

#include <cairo/cairo.h>
#include <cairo/cairo-xlib.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <mutex>
#include <string>

#define TOCAKE_UI_URI "https://nylarea.com/plugins/omegacake#ui"

namespace {

// Port indices — must match src/lv2.cpp / the .ttl.
enum PortIndex {
    PORT_SUSTAIN     = 2,
    PORT_TONE        = 3,
    PORT_VOLUME      = 4,
    PORT_HQ          = 5,
    PORT_VARIANT     = 7,
    PORT_TONE_BYPASS = 8,
    PORT_ENABLED     = 9,
    PORT_TONE_SWITCH = 10,
};

// A host may send a non-finite control value (a broken automation lane, a restored session). One NaN
// reaching the drawing puts the cairo context in a sticky error state and the panel stops painting;
// one NaN in the panel state is written back to the plugin on the next drag. So port_event ignores
// non-finite values, and clamp01 maps NaN to 0 rather than passing it.

inline double clamp01(double v) { return v > 0.0 ? (v < 1.0 ? v : 1.0) : 0.0; }
inline int clampi(int v, int lo, int hi) { return v < lo ? lo : (v > hi ? hi : v); }

// HQ is clickable only with the pedal on.

namespace photo = tocake::ui::photo;

// The knobs' values live by PORT ("slot": 0 sustain, 1 tone, 2 volume); a layout lists them in its own
// order (photo::KNOBS_THREE: sustain, tone, volume; KNOBS_ONE: volume).
constexpr uint32_t SLOT_PORT[3] = { PORT_SUSTAIN, PORT_TONE, PORT_VOLUME };
inline int slot_of(int cfg, int k) { return cfg == photo::CFG_ONE ? 2 : k; }

double seconds_now()
{
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}

// Drag travel: 200 logical px for 0->1; the wheel moves 2 % a notch.
constexpr double DRAG_PX    = 200.0;
constexpr double WHEEL_STEP = 0.02;

// --- the model list ------------------------------------------------------------
struct Popup {
    Window           win  = 0;
    cairo_surface_t* surf = nullptr;
    cairo_t*         cr   = nullptr;
    cairo_surface_t* buf   = nullptr;   // off-screen double buffer
    cairo_t*         bufcr = nullptr;
    int  w = 0, h = 0, row_h = 0;
    int  visible = 0;                   // rows shown at once (the models do not fit a screen)
    int  top = 0;                       // first visible row
    int  hover = -1;                    // row under the mouse (absolute), or -1
    int  frame = 0;                     // the bevelled frame around the rows
    bool open = false;
    Window focus_at_open = 0;           // who had the keyboard focus when it opened
};

struct UI {
    Display* dpy    = nullptr;
    Window   win    = 0;
    Window   parent = 0;
    Window   root   = 0;
    Visual*  visual = nullptr;
    Colormap cmap   = 0;                // ours, when the parent's visual is not the default
    int      depth  = 0;                // the window's (the parent's); the popup uses it too (B21)
    cairo_surface_t* surf = nullptr;
    cairo_t*         cr   = nullptr;
    cairo_surface_t* buf   = nullptr;   // off-screen double buffer: composed here, blitted once
    cairo_t*         bufcr = nullptr;

    int    w = 0, h = 0;
    double scale = 1.0;                 // the host's `ui:scaleFactor`: the window's size at open
    double fit = 0.5;                   // window pixels per asset pixel: the pedal fitted to the window
    int    ox = 0, oy = 0;              // where the fitted pedal starts in the window (centred)

    LV2UI_Write_Function write = nullptr;
    LV2UI_Controller     ctl   = nullptr;
    const LV2UI_Touch*   touch = nullptr;

    tocake::ui::Fonts fonts;
    std::string       photo_dir;        // the bundle's photo/
    photo::Scaled     scaled;           // the images at this window's scale
    bool              cache_user = false;   // counted among the shared image cache's users
    photo::Assets     own[2];           // [cfg] decoded for this instance alone: another bundle's images
    int    failed_cfg = -1;             // a layout whose images did not load, at
    double failed_fit = -1.0;           // this scale: not retried until either changes

    // What the panel shows: the ports' values.
    double knob[3] = { 0.8, 0.5, 0.6 };   // by slot: sustain, tone, volume (the plugin's defaults)
    bool   knob_dirty[3] = {};             // moved, not yet written to the port
    bool   on = true;                      // `enabled`
    bool   hq = false;
    bool   tone_bypass = false, tone_switch = false;
    int    cur_variant = 0;
    bool   display_hover = false;
    bool   pressed = false;                // the footswitch held down

    // The two animations, driven by the clock (not by how often the host idles).
    double lit = 1.0, lit_from = 1.0, lit_t0 = -1.0;     // the look crossfade on the footswitch
    double fsw = 0.0, fsw_from = 0.0, fsw_t0 = -1.0;     // plunger travel

    int    drag = -1;                      // knob index in the current layout, or -1
    double drag_y0 = 0.0, drag_v0 = 0.0;

    Popup popup;
    bool  needs_redraw = true;
    mutable std::string readout;           // what the display shows: year and name (view_of composes it)
};

int display_pos_of(int model);
int menu_row_of(int model);
int controls_of(const UI* ui) { return tocake::VARIANT_CONTROLS[ui->cur_variant]; }
int cfg_of(const UI* ui) { return photo::cfg_of(controls_of(ui)); }

// A drag moves the knob on every motion event, but the port is written once per idle turn with the
// latest value: the host gets the same final value without a write per pixel of travel.
void set_knob(UI* ui, int slot, double v01)
{
    ui->knob[slot] = clamp01(v01);
    ui->knob_dirty[slot] = true;
    ui->needs_redraw = true;
}

void flush_knobs(UI* ui)
{
    for (int s = 0; s < 3; ++s) {
        if (!ui->knob_dirty[s]) continue;
        ui->knob_dirty[s] = false;
        if (ui->write) {
            float pv = float(ui->knob[s]);
            ui->write(ui->ctl, SLOT_PORT[s], sizeof(float), 0, &pv);
        }
    }
}

void notify_touch(UI* ui, uint32_t port, bool held)
{
    if (ui->touch && ui->touch->touch) ui->touch->touch(ui->touch->handle, port, held);
}

void write_port(UI* ui, uint32_t port, float v)
{
    notify_touch(ui, port, true);
    if (ui->write) ui->write(ui->ctl, port, sizeof(float), 0, &v);
    notify_touch(ui, port, false);
}

// The same factor as a UI scale (1 = the logical window), for sizes given in logical px.
double ui_scale(const UI* ui) { return ui->fit * 2.0; }

void update_fit(UI* ui)
{
    const double f = std::fmin(ui->w / photo::ASSET_W, ui->h / photo::ASSET_H);
    ui->fit = f > 0.0 ? f : ui->scale * 0.5;
    ui->ox = int((ui->w - std::ceil(photo::ASSET_W * ui->fit)) / 2.0);
    ui->oy = int((ui->h - std::ceil(photo::ASSET_H * ui->fit)) / 2.0);
    if (ui->ox < 0) ui->ox = 0;
    if (ui->oy < 0) ui->oy = 0;
}

// Window coordinates -> ASSET coordinates, the grid the hit zones are given in.
void to_asset(const UI* ui, double px, double py, double* x, double* y)
{
    *x = (px - ui->ox) / ui->fit;
    *y = (py - ui->oy) / ui->fit;
}

int hit_knob(const UI* ui, double px, double py)
{
    double x, y; to_asset(ui, px, py, &x, &y);
    return photo::knob_at(cfg_of(ui), x, y);
}

bool hit(const UI* ui, const photo::Box& b, double px, double py)
{
    double x, y; to_asset(ui, px, py, &x, &y);
    return photo::in_box(b, x, y);
}

// The toggle with its two legends: the region its state frames are pasted over.
bool hit_toggle_region(const UI* ui, double px, double py)
{
    const photo::IBox& r = photo::toggle_box(photo::set_of(controls_of(ui)));
    return hit(ui, photo::Box{ double(r.x0), double(r.y0), double(r.x1), double(r.y1) }, px, py);
}

// Start an animation from wherever it is now, so a reversal mid-way does not jump.
void start_fade(UI* ui) { ui->lit_from = ui->lit; ui->lit_t0 = seconds_now(); ui->needs_redraw = true; }
void start_fsw(UI* ui)  { ui->fsw_from = ui->fsw; ui->fsw_t0 = seconds_now(); ui->needs_redraw = true; }

// An animation that ENDS in this call still needs its last frame painted: the paint is asked for by
// animate()'s return, and a finished animation used to return "not moving", so the plunger and the
// crossfade stayed on the frame before their last until another event repainted. The call that finishes
// an animation returns true once more.

constexpr bool kPaintFinalFrame = true;

// Advance both animations to `now`; true while either still moves, and on the call that finishes one.
bool animate(UI* ui, double now)
{
    bool moving = false;
    const double lit_to = ui->on ? 1.0 : 0.0;
    if (ui->lit_t0 >= 0.0) {
        const double t = (now - ui->lit_t0) / photo::FADE_S;
        if (t >= 1.0) { ui->lit = lit_to; ui->lit_t0 = -1.0; moving |= kPaintFinalFrame; }
        else { ui->lit = ui->lit_from + (lit_to - ui->lit_from) * std::fmax(0.0, t); moving = true; }
    }
    const double fsw_to = ui->pressed ? 1.0 : 0.0;
    if (ui->fsw_t0 >= 0.0) {
        const double t = (now - ui->fsw_t0) / photo::FSW_S;
        if (t >= 1.0) { ui->fsw = fsw_to; ui->fsw_t0 = -1.0; moving |= kPaintFinalFrame; }
        else { ui->fsw = ui->fsw_from + (fsw_to - ui->fsw_from) * std::fmax(0.0, t); moving = true; }
    }
    return moving;
}

void set_variant(UI* ui, int idx)
{
    idx = clampi(idx, 0, tocake::VARIANT_LABEL_COUNT - 1);
    ui->cur_variant = idx;
    write_port(ui, PORT_VARIANT, float(idx));
    ui->needs_redraw = true;
}

// The order the model list shows (and the wheel walks): display position -> model index. Until the list
// is regrouped it is the model order; the port value is always the MODEL index, whatever the display.
int display_pos_of(int model)
{
    for (int i = 0; i < tocake::VARIANT_LABEL_COUNT; ++i) if (tocake::VARIANT_DISPLAY_ORDER[i] == model) return i;
    return 0;
}

// The menu row (headers included) that shows a model.
int menu_row_of(int model)
{
    for (int r = 0; r < tocake::VARIANT_MENU_COUNT; ++r) if (tocake::VARIANT_MENU[r].model == model) return r;
    return 0;
}

// The wheel over the display steps to the next model DOWN the list (Button5) or up (Button4), and stops
// at either end.
void step_variant(UI* ui, int dir)
{
    const int pos = display_pos_of(ui->cur_variant) + dir;
    if (pos < 0 || pos >= tocake::VARIANT_LABEL_COUNT) return;
    set_variant(ui, tocake::VARIANT_DISPLAY_ORDER[pos]);
}

void set_enabled(UI* ui, bool on)
{
    ui->on = on;
    write_port(ui, PORT_ENABLED, on ? 1.0f : 0.0f);
    start_fade(ui);
}

// --- the model list ------------------------------------------------------------

// The list continues the display: dark glass inside a bevelled frame, the dot-matrix face in the
// readout's lit purple. The models do not fit: a window of rows that scrolls with the wheel.
void popup_render(UI* ui)
{
    Popup& p = ui->popup;
    if (!p.bufcr) return;
    cairo_t* cr = p.bufcr;
    const double fw = p.frame;
    photo::set_hex(cr, 0x2b2c2e, 1.0);
    cairo_paint(cr);
    photo::set_hex(cr, 0x55575a, 1.0);
    cairo_rectangle(cr, 0, 0, p.w, 1.0);
    cairo_fill(cr);
    photo::set_hex(cr, 0x151516, 1.0);
    cairo_rectangle(cr, 0, p.h - 1.0, p.w, 1.0);
    cairo_fill(cr);
    photo::set_hex(cr, 0x050505, 1.0);
    cairo_rectangle(cr, fw - 1.0, fw - 1.0, p.w - 2.0 * fw + 2.0, p.h - 2.0 * fw + 2.0);
    cairo_fill(cr);
    photo::set_hex(cr, 0x0b0b0c, 1.0);
    cairo_rectangle(cr, fw, fw, p.w - 2.0 * fw, p.h - 2.0 * fw);
    cairo_fill(cr);

    if (ui->fonts.doto) cairo_set_font_face(cr, ui->fonts.doto);
    else cairo_select_font_face(cr, "monospace", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
    cairo_set_font_size(cr, p.row_h * 0.58);
    const double dot_stroke = photo::DOT_STROKE_PX * ui_scale(ui);
    const double indent = p.row_h * 0.75;
    const int last = clampi(p.top + p.visible, 0, tocake::VARIANT_MENU_COUNT);
    for (int r = p.top; r < last; ++r) {           // r: a row of the menu (a header or a model)
        const tocake::VariantMenuRow& row = tocake::VARIANT_MENU[r];
        const double y = fw + (r - p.top) * p.row_h;
        const double x0 = fw + p.row_h * 0.45 + row.level * indent;
        std::string text;
        if (row.model < 0) {
            // a header: dimmer, with a rule from its text to the right edge
            text = row.header;
            for (char& c : text) c = char(std::toupper(static_cast<unsigned char>(c)));
            cairo_text_extents_t e;
            cairo_text_extents(cr, text.c_str(), &e);
            const double ty = y + p.row_h / 2.0 - e.height / 2.0 - e.y_bearing;
            photo::set_hex(cr, photo::READOUT_HEX, row.level == 0 ? 0.75 : 0.55);
            cairo_move_to(cr, x0 - e.x_bearing, ty);
            cairo_text_path(cr, text.c_str());
            cairo_set_line_width(cr, dot_stroke);
            cairo_stroke_preserve(cr);
            cairo_fill(cr);
            photo::set_hex(cr, photo::READOUT_HEX, row.level == 0 ? 0.45 : 0.25);
            cairo_rectangle(cr, x0 + e.width + p.row_h * 0.4, y + p.row_h / 2.0, p.w - fw - 8.0 - (x0 + e.width + p.row_h * 0.4), 1.0);
            cairo_fill(cr);
            continue;
        }
        const bool cur = row.model == ui->cur_variant;
        if (cur || r == p.hover) {
            photo::set_hex(cr, photo::READOUT_HEX, cur ? 0.22 : 0.10);
            cairo_rectangle(cr, fw, y, p.w - 2.0 * fw, p.row_h);
            cairo_fill(cr);
        }
        text = std::to_string(tocake::VARIANT_YEARS[row.model]) + "  " + tocake::VARIANT_LABELS[row.model];
        for (char& c : text) c = char(std::toupper(static_cast<unsigned char>(c)));
        cairo_text_extents_t e;
        cairo_text_extents(cr, text.c_str(), &e);
        photo::set_hex(cr, photo::READOUT_HEX, cur || r == p.hover ? 1.0 : 0.85);
        cairo_move_to(cr, x0 - e.x_bearing, y + p.row_h / 2.0 - e.height / 2.0 - e.y_bearing);
        cairo_text_path(cr, text.c_str());
        cairo_set_line_width(cr, dot_stroke);
        cairo_stroke_preserve(cr);
        cairo_fill(cr);
    }
    if (tocake::VARIANT_MENU_COUNT > p.visible) {      // where the window of rows sits in the list
        const double track = p.h - 2.0 * fw;
        photo::set_hex(cr, photo::READOUT_HEX, 0.35);
        cairo_rectangle(cr, p.w - fw - 3.0, fw + track * p.top / tocake::VARIANT_MENU_COUNT, 3.0,
                        track * p.visible / tocake::VARIANT_MENU_COUNT);
        cairo_fill(cr);
    }
    cairo_surface_flush(p.buf);
    cairo_set_source_surface(p.cr, p.buf, 0, 0);
    cairo_paint(p.cr);
    cairo_surface_flush(p.surf);
}

void popup_close(UI* ui)
{
    Popup& p = ui->popup;
    if (!p.open) return;
    XUngrabPointer(ui->dpy, CurrentTime);
    XUngrabKeyboard(ui->dpy, CurrentTime);
    if (p.bufcr) { cairo_destroy(p.bufcr); p.bufcr = nullptr; }
    if (p.buf)   { cairo_surface_destroy(p.buf); p.buf = nullptr; }
    if (p.cr)    { cairo_destroy(p.cr); p.cr = nullptr; }
    if (p.surf)  { cairo_surface_destroy(p.surf); p.surf = nullptr; }
    if (p.win)   { XDestroyWindow(ui->dpy, p.win); p.win = 0; }
    p.open = false;
    p.hover = -1;
    ui->needs_redraw = true;            // the display's light follows the list
    XFlush(ui->dpy);
}

void popup_open(UI* ui)
{
    Popup& p = ui->popup;
    if (p.open) { popup_close(ui); return; }

    const int scr = DefaultScreen(ui->dpy);
    const photo::Box& gl = photo::DISPLAY;          // the black glass
    const photo::Box& bz = photo::DISPLAY_BEZEL;    // its frame
    const double f = ui->fit;
    p.row_h = int(26.0 * ui_scale(ui));
    if (p.row_h < 16) p.row_h = 16;
    p.frame = int(std::lround(4.0 * ui_scale(ui)));
    p.visible = tocake::VARIANT_MENU_COUNT < 18 ? tocake::VARIANT_MENU_COUNT : 18;
    // wider than the glass: three levels of indentation and "1976  JUMBO COLOUR CAKE" must fit; centred on it
    p.w = std::max(int(std::lround((gl.x1 - gl.x0) * f)), int(std::lround(330.0 * ui_scale(ui))));
    p.h = p.row_h * p.visible + 2 * p.frame;
    p.top = clampi(menu_row_of(ui->cur_variant) - p.visible / 2, 0, tocake::VARIANT_MENU_COUNT - p.visible);
    p.hover = -1;

    // The display sits at the bottom of the pedal: the list opens UPWARDS from its frame, and downwards
    // only if the screen has no room above.
    int rx = 0, ry = 0; Window child = 0;
    const int glass_cx = ui->ox + int(std::lround((gl.x0 + gl.x1) / 2.0 * f));
    XTranslateCoordinates(ui->dpy, ui->win, ui->root, glass_cx - p.w / 2,
                          ui->oy + int(bz.y0 * f) - p.h, &rx, &ry, &child);
    if (rx < 0) rx = 0;
    if (ry < 0) ry = ry + p.h + int((bz.y1 - bz.y0) * f);
    const int sh = DisplayHeight(ui->dpy, scr);
    if (ry + p.h > sh) ry = sh - p.h;
    if (ry < 0) ry = 0;

    // B21: the popup is a child of the ROOT but draws with the panel's visual, so it needs that
    // visual's DEPTH and colormap too: under a 32-bit ARGB parent the root's depth does not match
    // (BadMatch, and the list never opened).
    XSetWindowAttributes attr{};
    attr.override_redirect = True;
    attr.background_pixel = BlackPixel(ui->dpy, scr);
    attr.border_pixel = 0;
    attr.event_mask = ExposureMask | ButtonPressMask | ButtonReleaseMask | PointerMotionMask |
                      LeaveWindowMask | KeyPressMask;
    unsigned long mask = CWOverrideRedirect | CWBackPixel | CWBorderPixel | CWEventMask;
    if (ui->cmap) { attr.colormap = ui->cmap; mask |= CWColormap; }
    p.win = XCreateWindow(ui->dpy, ui->root, rx, ry, p.w, p.h, 0,
                          ui->depth ? ui->depth : CopyFromParent, InputOutput, ui->visual, mask, &attr);
    p.surf = cairo_xlib_surface_create(ui->dpy, p.win, ui->visual, p.w, p.h);
    p.cr   = cairo_create(p.surf);
    p.buf  = cairo_surface_create_similar(p.surf, CAIRO_CONTENT_COLOR, p.w, p.h);
    p.bufcr = cairo_create(p.buf);

    XMapRaised(ui->dpy, p.win);
    // owner_events False: every pointer event goes to the popup in its own coordinates, so a click
    // outside arrives out of range and closes. Without the grab that routing is false: the panel's
    // own clicks would be read in the popup's coordinates and select a model. So no grab, no popup.
    p.open = true;
    if (XGrabPointer(ui->dpy, p.win, False, ButtonPressMask | ButtonReleaseMask | PointerMotionMask,
                     GrabModeAsync, GrabModeAsync, None, None, CurrentTime) != GrabSuccess) {
        popup_close(ui);
        return;
    }
    int revert = 0;
    XGetInputFocus(ui->dpy, &p.focus_at_open, &revert);
    // the keyboard too, so Escape closes the list; a failed keyboard grab leaves it usable by mouse
    XGrabKeyboard(ui->dpy, p.win, False, GrabModeAsync, GrabModeAsync, CurrentTime);
    XFlush(ui->dpy);
    popup_render(ui);
    ui->needs_redraw = true;
}

// Returns true if the event was consumed by the popup.
bool popup_handle(UI* ui, XEvent* ev)
{
    Popup& p = ui->popup;
    if (!p.open) return false;
    switch (ev->type) {
        case Expose:
            popup_render(ui);
            return true;
        case MotionNotify: {
            const int x = ev->xmotion.x, y = ev->xmotion.y;
            int s = -1;
            if (x >= p.frame && x < p.w - p.frame && y >= p.frame && y < p.h - p.frame)
                s = p.top + (y - p.frame) / p.row_h;
            if (s >= 0 && (s >= tocake::VARIANT_MENU_COUNT || tocake::VARIANT_MENU[s].model < 0)) s = -1;   // headers do not light
            if (s != p.hover) { p.hover = s; popup_render(ui); }
            return true;
        }
        case ButtonPress: {
            const unsigned b = ev->xbutton.button;
            const int x = ev->xbutton.x, y = ev->xbutton.y;
            if (b == Button4 || b == Button5) {
                const int maxtop = tocake::VARIANT_MENU_COUNT - p.visible;
                p.top = clampi(p.top + (b == Button4 ? -1 : 1), 0, maxtop < 0 ? 0 : maxtop);
                popup_render(ui);
                return true;
            }
            if (b == Button1) {
                // the frame is not a row: (y - frame) / row_h would truncate a click on the top edge
                const bool inside = x >= p.frame && x < p.w - p.frame && y >= p.frame && y < p.h - p.frame;
                if (inside) {
                    const int r = p.top + (y - p.frame) / p.row_h;
                    if (r >= 0 && r < tocake::VARIANT_MENU_COUNT && tocake::VARIANT_MENU[r].model >= 0)
                        set_variant(ui, tocake::VARIANT_MENU[r].model);
                    else return true;        // a header: the list stays open
                }
                popup_close(ui);              // inside or out, one click closes
            }
            return true;
        }
        case KeyPress:
            if (XLookupKeysym(&ev->xkey, 0) == XK_Escape) popup_close(ui);
            return true;
        default:
            return true;   // while it holds the pointer grab it eats everything
    }
}

// --- the decoded images, shared by every instance in the process -------------------

// Decoding the bundle's PNGs is most of the cost of opening the editor, and every instance decodes
// the same files. So they are decoded once, per layout, into a process-wide cache under a mutex
// (hosts may open UIs from different threads; an instance opening during another's decode waits for
// it instead of decoding again). The cache counts its users and is NOT freed when the last one
// closes: a host that closes and reopens the editor reopens it without decoding. It is freed when
// this library is unloaded (dlclose, or the process ending), by the cache's destructor; no instance
// is alive then, because a host destroys its UIs before unloading their library. One cache holds one
// bundle's images: an instance from another bundle, while the cache is in use, decodes its own copy.

struct ImageCache {
    std::mutex    mx;
    std::string   dir;
    photo::Assets set[2];               // [cfg]
    int           users = 0;
    ~ImageCache()
    {
        for (auto& a : set) photo::free_assets(a);
    }
};
ImageCache g_images;

// The images of `cfg` scaled into ui->scaled at ui->fit, decoded only if this process has not
// decoded them yet. False when they do not load (the panel draws its placeholder).
bool scale_images(UI* ui, int cfg)
{
    {   // the process-wide cache, under its lock
        std::lock_guard<std::mutex> lock(g_images.mx);
        if (g_images.users == 0 && g_images.dir != ui->photo_dir) {
            for (auto& a : g_images.set) photo::free_assets(a);
            g_images.dir = ui->photo_dir;
        }
        if (g_images.dir == ui->photo_dir) {
            photo::Assets& a = g_images.set[cfg];
            if (!a.ok) {
                photo::free_assets(a);
                a = photo::load_assets(ui->photo_dir, cfg, &g_images.set[1 - cfg]);
                if (!a.ok) { photo::free_assets(a); return false; }
            }
            if (!ui->cache_user) { ++g_images.users; ui->cache_user = true; }
            // scaled under the lock: the decoded surfaces are shared, and a cairo object is not used
            // from two threads at once
            photo::ensure_scaled(ui->scaled, a, ui->fit);
            return true;
        }
    }
    photo::Assets& a = ui->own[cfg];
    if (!a.ok) {
        photo::free_assets(a);
        a = photo::load_assets(ui->photo_dir, cfg, &ui->own[1 - cfg]);
        if (!a.ok) { photo::free_assets(a); return false; }
    }
    photo::ensure_scaled(ui->scaled, a, ui->fit);
    return true;
}

void release_images(UI* ui)
{
    if (ui->cache_user) {
        std::lock_guard<std::mutex> lock(g_images.mx);
        --g_images.users;
        ui->cache_user = false;
    }
    for (auto& a : ui->own) photo::free_assets(a);
}

// --- painting --------------------------------------------------------------------

void make_buffer(UI* ui)
{
    if (ui->bufcr) { cairo_destroy(ui->bufcr); ui->bufcr = nullptr; }
    if (ui->buf)   { cairo_surface_destroy(ui->buf); ui->buf = nullptr; }
    ui->buf = cairo_surface_create_similar(ui->surf, CAIRO_CONTENT_COLOR, ui->w, ui->h);
    ui->bufcr = cairo_create(ui->buf);
}

photo::View view_of(const UI* ui)
{
    photo::View v;
    v.sustain = ui->knob[0]; v.tone = ui->knob[1]; v.volume = ui->knob[2];
    v.controls = controls_of(ui);
    v.lit = ui->lit;
    v.fsw = ui->fsw;
    v.toggle_up = v.controls == 2 ? ui->tone_switch : ui->tone_bypass;
    v.hq = ui->hq;
    // In bypass the display is dark, and lights while hovered or while the list is open, so a model
    // is never changed blind.
    v.display_lit = ui->on || ui->display_hover || ui->popup.open;
    // the display shows the year and the name
    ui->readout = std::to_string(tocake::VARIANT_YEARS[ui->cur_variant]) + " " + tocake::VARIANT_LABELS[ui->cur_variant];
    v.variant = ui->readout.c_str();
    v.face = ui->fonts.doto;
    return v;
}

void render(UI* ui)
{
    if (!ui->cr || !ui->bufcr) return;
    const photo::View v = view_of(ui);
    const int cfg = photo::cfg_of(v.controls);
    // Painting needs only the scaled images. A new layout or a new window size scales them again
    // from the decoded ones, which the process keeps (see ImageCache).
    const bool failed_here = ui->failed_cfg == cfg && ui->failed_fit == ui->fit;
    if ((ui->scaled.cfg != cfg || ui->scaled.f != ui->fit) && !failed_here) {
        photo::free_scaled(ui->scaled);
        // missing or damaged images: the placeholder, and no new attempt per paint at this cfg and size
        if (!scale_images(ui, cfg)) { ui->failed_cfg = cfg; ui->failed_fit = ui->fit; }
    }

    cairo_t* cr = ui->bufcr;
    cairo_identity_matrix(cr);
    cairo_set_antialias(cr, CAIRO_ANTIALIAS_GOOD);
    if (ui->ox > 0 || ui->oy > 0) {        // the margins around the fitted pedal
        photo::set_hex(cr, photo::BACKDROP_HEX, 1.0);
        cairo_paint(cr);
    }
    cairo_translate(cr, ui->ox, ui->oy);   // integer: the frames stay on the plate's pixel grid
    photo::draw_panel(cr, ui->scaled, v);
    cairo_surface_flush(ui->buf);
    // a single blit to the window: it is never seen half-drawn
    cairo_identity_matrix(ui->cr);
    cairo_set_source_surface(ui->cr, ui->buf, 0, 0);
    cairo_paint(ui->cr);
    cairo_surface_flush(ui->surf);
    XFlush(ui->dpy);
}

void end_drag(UI* ui)
{
    if (ui->drag < 0) return;
    flush_knobs(ui);
    notify_touch(ui, SLOT_PORT[slot_of(cfg_of(ui), ui->drag)], false);
    ui->drag = -1;
}

void handle_event(UI* ui, XEvent* ev)
{
    switch (ev->type) {
        case Expose:
            if (ev->xexpose.count == 0) ui->needs_redraw = true;
            break;

        case ConfigureNotify: {
            const int nw = ev->xconfigure.width, nh = ev->xconfigure.height;
            if (nw != ui->w || nh != ui->h) {
                ui->w = nw; ui->h = nh;
                update_fit(ui);
                cairo_xlib_surface_set_size(ui->surf, nw, nh);
                make_buffer(ui);
                ui->needs_redraw = true;
            }
            break;
        }

        case ButtonPress: {
            const unsigned b = ev->xbutton.button;
            const double px = ev->xbutton.x, py = ev->xbutton.y;
            const int cfg = cfg_of(ui);
            if (b == Button1) {
                // No zone overlaps another (the scene's projection), so the order of the tests does
                // not decide anything.
                const int i = hit_knob(ui, px, py);
                if (i >= 0) {
                    ui->drag = i;
                    ui->drag_y0 = py;
                    ui->drag_v0 = ui->knob[slot_of(cfg, i)];
                    notify_touch(ui, SLOT_PORT[slot_of(cfg, i)], true);
                } else if (hit(ui, photo::FSW_HIT, px, py)) {
                    ui->pressed = true;              // it sinks while held down
                    start_fsw(ui);
                } else if (controls_of(ui) != 0 && hit_toggle_region(ui, px, py)) {
                    // The lever flips; a LEGEND selects its own position: the upper
                    // one (BYPASS / BRIGHT) up, the lower one (TONE / DARK) down.
                    double ax, ay; to_asset(ui, px, py, &ax, &ay);
                    const bool has_tone_switch = controls_of(ui) == 2;
                    bool& state = has_tone_switch ? ui->tone_switch : ui->tone_bypass;
                    const bool want = ay < photo::TOGGLE_HIT.y0 ? true : ay > photo::TOGGLE_HIT.y1 ? false : !state;
                    if (want != state) {
                        state = want;
                        write_port(ui, has_tone_switch ? PORT_TONE_SWITCH : PORT_TONE_BYPASS, want ? 1.0f : 0.0f);
                        ui->needs_redraw = true;
                    }
                } else if (ui->on && hit(ui, photo::HQ_HIT, px, py)) {
                    // HQ is selectable only with the pedal on
                    ui->hq = !ui->hq;
                    write_port(ui, PORT_HQ, ui->hq ? 1.0f : 0.0f);
                    ui->needs_redraw = true;
                } else if (hit(ui, photo::DISPLAY_BEZEL, px, py)) {
                    popup_open(ui);
                }
            } else if ((b == Button4 || b == Button5) && hit(ui, photo::DISPLAY_BEZEL, px, py)) {
                step_variant(ui, b == Button5 ? +1 : -1);
            } else if (b == Button4 || b == Button5) {
                const int i = hit_knob(ui, px, py);
                if (i >= 0) {
                    const int s = slot_of(cfg, i);
                    // the wheel is a gesture too: bracketed with touch, or the host sees an orphan write
                    notify_touch(ui, SLOT_PORT[s], true);
                    set_knob(ui, s, ui->knob[s] + (b == Button4 ? WHEEL_STEP : -WHEEL_STEP));
                    flush_knobs(ui);
                    notify_touch(ui, SLOT_PORT[s], false);
                }
            }
            break;
        }

        case ButtonRelease:
            if (ev->xbutton.button == Button1) {
                end_drag(ui);
                if (ui->pressed) {
                    ui->pressed = false;
                    start_fsw(ui);
                    // the footswitch toggles on release, and only released on top: dragging away
                    // and releasing cancels, like any button
                    if (hit(ui, photo::FSW_HIT, ev->xbutton.x, ev->xbutton.y)) set_enabled(ui, !ui->on);
                }
            }
            break;

        case MotionNotify:
            if (ui->drag >= 0 && !(ev->xmotion.state & Button1Mask)) {
                // the button is up and the release never arrived (another client grabbed the pointer
                // mid-drag): end the drag here, or the knob stays stuck to the mouse
                end_drag(ui);
            } else if (ui->drag >= 0) {
                const double dy = ui->drag_y0 - ev->xmotion.y;
                set_knob(ui, slot_of(cfg_of(ui), ui->drag), ui->drag_v0 + dy / (DRAG_PX * ui_scale(ui)));
            } else {
                const bool over = hit(ui, photo::DISPLAY_BEZEL, ev->xmotion.x, ev->xmotion.y);
                if (over != ui->display_hover) { ui->display_hover = over; ui->needs_redraw = true; }
            }
            break;

        case LeaveNotify:
            if (ui->display_hover) { ui->display_hover = false; ui->needs_redraw = true; }
            break;

        case UnmapNotify:          // the pedal's window hidden: its list goes with it
            popup_close(ui);
            break;

        default:
            break;
    }
}

// --- LV2UI entry points -------------------------------------------------------

// X errors. Xlib's default handler exit()s the whole process -- the HOST -- on any protocol error,
// and a parent window with a visual we do not match (a 32-bit ARGB parent) or one the server does not
// know yet raises one. While any instance of this UI is alive, a handler records errors on OUR
// connections and hands every other error to whatever handler was installed before (the host's).
// The handler is process-wide by Xlib's design, so it is installed with the first instance and
// restored with the last.

int (*g_prev_x_handler)(Display*, XErrorEvent*) = nullptr;
int  g_ui_instances = 0;
Display* g_ui_displays[16] = {};
bool g_x_error = false;

int ui_x_error(Display* d, XErrorEvent* e)
{
    for (Display* ours : g_ui_displays)
        if (ours == d) { g_x_error = true; return 0; }
    return g_prev_x_handler ? g_prev_x_handler(d, e) : 0;
}

void track_display(Display* d, bool add)
{
    for (Display*& slot : g_ui_displays)
        if (add ? slot == nullptr : slot == d) { slot = add ? d : nullptr; return; }
}

void cleanup(LV2UI_Handle handle);

LV2UI_Handle instantiate(const LV2UI_Descriptor*, const char*, const char* bundle_path,
                         LV2UI_Write_Function write_function, LV2UI_Controller controller,
                         LV2UI_Widget* widget, const LV2_Feature* const* features)
{
    Window parent = 0;
    const LV2UI_Resize* resize_host = nullptr;
    const LV2UI_Touch*  touch = nullptr;
    const LV2_Options_Option* options = nullptr;
    const LV2_URID_Map* map = nullptr;
    for (int i = 0; features && features[i]; ++i) {
        if (!std::strcmp(features[i]->URI, LV2_UI__parent))
            parent = static_cast<Window>(reinterpret_cast<uintptr_t>(features[i]->data));
        else if (!std::strcmp(features[i]->URI, LV2_UI__resize))
            resize_host = static_cast<const LV2UI_Resize*>(features[i]->data);
        else if (!std::strcmp(features[i]->URI, LV2_UI__touch))
            touch = static_cast<const LV2UI_Touch*>(features[i]->data);
        else if (!std::strcmp(features[i]->URI, LV2_OPTIONS__options))
            options = static_cast<const LV2_Options_Option*>(features[i]->data);
        else if (!std::strcmp(features[i]->URI, LV2_URID__map))
            map = static_cast<const LV2_URID_Map*>(features[i]->data);
    }

    auto* ui = new UI();
    ui->write = write_function;
    ui->ctl   = controller;
    ui->touch = touch;

    // `ui:scaleFactor`: the host states what scale it wants the window at on HiDPI screens; 1 if not.
    ui->scale = 1.0;
    if (options && map) {
        const LV2_URID urid_scale = map->map(map->handle, LV2_UI__scaleFactor);
        const LV2_URID urid_float = map->map(map->handle, LV2_ATOM__Float);
        for (const LV2_Options_Option* o = options; o->key; ++o)
            if (o->key == urid_scale && o->type == urid_float && o->value) {
                const float f = *static_cast<const float*>(o->value);
                if (f > 0.1f && f < 8.0f) ui->scale = double(f);
            }
    }
    ui->w = int(std::ceil(photo::LOGICAL_W * ui->scale));
    ui->h = int(std::ceil(photo::LOGICAL_H * ui->scale));
    update_fit(ui);

    ui->dpy = XOpenDisplay(nullptr);
    if (!ui->dpy) { delete ui; return nullptr; }

    const int screen = DefaultScreen(ui->dpy);
    ui->visual = DefaultVisual(ui->dpy, screen);
    ui->root   = RootWindow(ui->dpy, screen);
    ui->parent = parent ? parent : ui->root;

    if (g_ui_instances++ == 0) { g_x_error = false; g_prev_x_handler = XSetErrorHandler(ui_x_error); }
    track_display(ui->dpy, true);
    g_x_error = false;
    // The window takes the PARENT's visual and depth, and so does the cairo surface: a child of
    // another depth cannot be drawn into its parent.
    int depth = DefaultDepth(ui->dpy, screen);
    XWindowAttributes pa;
    if (XGetWindowAttributes(ui->dpy, ui->parent, &pa) && !g_x_error) {
        ui->visual = pa.visual;
        depth = pa.depth;
    }
    XSetWindowAttributes wa{};
    unsigned long mask = CWBackPixel | CWBorderPixel;
    wa.background_pixel = 0;
    wa.border_pixel = 0;
    if (ui->visual != DefaultVisual(ui->dpy, screen)) {
        ui->cmap = XCreateColormap(ui->dpy, ui->root, ui->visual, AllocNone);
        wa.colormap = ui->cmap;
        mask |= CWColormap;
    }
    ui->depth = depth;
    ui->win = XCreateWindow(ui->dpy, ui->parent, 0, 0, (unsigned)ui->w, (unsigned)ui->h, 0, depth,
                            InputOutput, ui->visual, mask, &wa);
    XSync(ui->dpy, False);
    if (!ui->win || g_x_error) {   // e.g. a parent the server does not know
        ui->win = 0;
        cleanup(ui);
        return nullptr;
    }

    // Size hints for a host that resizes the plugin window: keep the pedal's proportions, never below
    // half its size. The drawing fits whatever size it gets.
    if (XSizeHints* sh = XAllocSizeHints()) {
        sh->flags = PMinSize | PAspect | PBaseSize;
        sh->min_width  = int(photo::LOGICAL_W / 2);
        sh->min_height = int(photo::LOGICAL_H / 2);
        sh->base_width = ui->w; sh->base_height = ui->h;
        sh->min_aspect.x = sh->max_aspect.x = int(photo::ASSET_W);
        sh->min_aspect.y = sh->max_aspect.y = int(photo::ASSET_H);
        XSetWMNormalHints(ui->dpy, ui->win, sh);
        XFree(sh);
    }

    XSelectInput(ui->dpy, ui->win, ExposureMask | StructureNotifyMask | ButtonPressMask |
                 ButtonReleaseMask | PointerMotionMask | ButtonMotionMask | LeaveWindowMask);

    ui->surf = cairo_xlib_surface_create(ui->dpy, ui->win, ui->visual, ui->w, ui->h);
    ui->cr   = cairo_create(ui->surf);
    make_buffer(ui);

    // The fonts and the panel images travel inside the bundle: loaded from there so the panel looks
    // the same in any host, whatever it has installed.
    // LV2 hands the bundle path WITH its trailing slash; every join below adds one (audit CA-C14).
    std::string bundle = bundle_path ? bundle_path : "";
    while (bundle.size() > 1 && bundle.back() == '/') bundle.pop_back();
    ui->fonts = tocake::ui::load_fonts(bundle);
    ui->photo_dir = (bundle.empty() ? std::string(".") : bundle) + "/photo";   // loaded at the first paint

    XMapWindow(ui->dpy, ui->win);
    XFlush(ui->dpy);
    if (resize_host && resize_host->ui_resize) resize_host->ui_resize(resize_host->handle, ui->w, ui->h);

    *widget = reinterpret_cast<LV2UI_Widget>(static_cast<uintptr_t>(ui->win));
    return ui;
}

void cleanup(LV2UI_Handle handle)
{
    auto* ui = static_cast<UI*>(handle);
    popup_close(ui);
    photo::free_scaled(ui->scaled);
    release_images(ui);
    tocake::ui::free_fonts(ui->fonts);
    if (ui->bufcr) cairo_destroy(ui->bufcr);
    if (ui->buf)   cairo_surface_destroy(ui->buf);
    if (ui->cr)    cairo_destroy(ui->cr);
    if (ui->surf)  cairo_surface_destroy(ui->surf);
    if (ui->win)   XDestroyWindow(ui->dpy, ui->win);
    if (ui->cmap)  XFreeColormap(ui->dpy, ui->cmap);
    if (ui->dpy) {
        XSync(ui->dpy, False);
        XCloseDisplay(ui->dpy);
        track_display(ui->dpy, false);
        if (--g_ui_instances == 0) XSetErrorHandler(g_prev_x_handler);
    }
    delete ui;
}

void port_event(LV2UI_Handle handle, uint32_t port, uint32_t buffer_size, uint32_t format, const void* buffer)
{
    if (format != 0 || buffer_size < sizeof(float)) return;
    auto* ui = static_cast<UI*>(handle);
    const float v = *static_cast<const float*>(buffer);
    if (!std::isfinite(v)) return;
    switch (port) {
        // a knob with a write still pending keeps the user's newer value
        case PORT_SUSTAIN: if (!ui->knob_dirty[0]) ui->knob[0] = clamp01(v); break;
        case PORT_TONE:    if (!ui->knob_dirty[1]) ui->knob[1] = clamp01(v); break;
        case PORT_VOLUME:  if (!ui->knob_dirty[2]) ui->knob[2] = clamp01(v); break;
        case PORT_HQ:          ui->hq = v >= 0.5f; break;
        case PORT_TONE_BYPASS: ui->tone_bypass = v >= 0.5f; break;
        case PORT_TONE_SWITCH: ui->tone_switch = v >= 0.5f; break;
        case PORT_ENABLED:
            if ((v >= 0.5f) != ui->on) { ui->on = v >= 0.5f; start_fade(ui); }
            break;
        case PORT_VARIANT: {
            // Clamp BEFORE converting: a float beyond int's range is UB.
            const float hi = static_cast<float>(tocake::VARIANT_LABEL_COUNT - 1);
            const int idx = static_cast<int>(std::lround(v < 0.0f ? 0.0f : (v > hi ? hi : v)));
            ui->cur_variant = idx;
            break;
        }
        default: return;
    }
    ui->needs_redraw = true;
}

// LV2 idle interface: pump X events, advance the animations, repaint when due.
int ui_idle(LV2UI_Handle handle)
{
    auto* ui = static_cast<UI*>(handle);
    if (!ui->dpy) return 1;
    while (XPending(ui->dpy)) {
        XEvent ev;
        XNextEvent(ui->dpy, &ev);
        if (ui->popup.open && (ev.xany.window == ui->popup.win || ev.type == ButtonPress ||
                               ev.type == ButtonRelease || ev.type == MotionNotify)) {
            popup_handle(ui, &ev);          // the popup holds the pointer grab
        } else {
            handle_event(ui, &ev);
        }
    }
    // The list closes when the keyboard focus moves away from where it was when it opened: an
    // override-redirect list has no window manager to hide it. Polled: under XWayland the grabs do
    // not stop a Wayland window taking the focus.
    if (ui->popup.open) {
        Window fw = 0; int revert = 0;
        XGetInputFocus(ui->dpy, &fw, &revert);
        if (fw != ui->popup.focus_at_open) popup_close(ui);
    }
    flush_knobs(ui);
    if (animate(ui, seconds_now())) ui->needs_redraw = true;
    if (ui->needs_redraw) {
        ui->needs_redraw = false;
        render(ui);
    }
    return 0;
}

const LV2UI_Idle_Interface idle_iface = { ui_idle };

const void* extension_data(const char* uri)
{
    if (!std::strcmp(uri, LV2_UI__idleInterface)) return &idle_iface;
    return nullptr;
}

const LV2UI_Descriptor ui_descriptor = {
    TOCAKE_UI_URI,
    instantiate,
    cleanup,
    port_event,
    extension_data,
};

} // namespace

extern "C" LV2_SYMBOL_EXPORT
const LV2UI_Descriptor* lv2ui_descriptor(uint32_t index)
{
    return (index == 0) ? &ui_descriptor : nullptr;
}
