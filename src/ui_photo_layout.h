// GENERATED FILE — DO NOT EDIT BY HAND.
//
// Where the photoreal GUI's parts land in its renders, in ASSET pixels (the renders' own
// grid, top-left origin), projected from the Blender scene. Regenerate with
// gui/scene/build_assets.sh then gui/scene/pack_assets.py.
#pragma once

namespace tocake {
namespace ui {
namespace photo {

inline constexpr double ASSET_W = 900;
inline constexpr double ASSET_H = 1350;

// A knob: its base circle (the hit zone and the paste mask) and its strip's box.
struct KnobGeo {
    const char* name;
    double cx, cy, r;            // base circle
    int    x0, y0, x1, y1;       // strip frame box
    int    frames;
    double deg_first, deg_last;  // pointer angle of the first and last frame
};
struct Box { double x0, y0, x1, y1; };
struct IBox { int x0, y0, x1, y1; };

inline constexpr KnobGeo KNOBS_THREE[3] = {
    { "sustain", 599.56, 335.49, 67.59, 525, 261, 674, 410, 64, -135.00, 135.00 },
    { "tone", 450.00, 506.09, 67.59, 376, 432, 524, 580, 64, -135.00, 135.00 },
    { "volume", 300.44, 335.49, 67.59, 226, 261, 375, 410, 64, -135.00, 135.00 },
};
inline constexpr KnobGeo KNOBS_ONE[1] = {
    { "volume", 450.00, 419.06, 85.61, 358, 328, 542, 511, 64, -135.00, 135.00 },
};

// pasted regions: the toggle with its legends (per legend set), the lit HQ button, the 3PDT
inline constexpr IBox TOGGLE_BYPASS = { 395, 798, 504, 967 };
inline constexpr IBox TOGGLE_TONE_SWITCH = { 395, 798, 504, 971 };
inline constexpr IBox HQ_LIT = { 551, 1036, 664, 1123 };
inline constexpr IBox FSW = { 366, 977, 534, 1160 };
inline constexpr int FSW_FRAMES = 5;

// hit zones and the display
inline constexpr Box FSW_HIT = { 403.34, 1026.18, 496.66, 1127.64 };
inline constexpr Box TOGGLE_HIT = { 432.78, 830.40, 467.22, 907.01 };
inline constexpr Box HQ_HIT = { 559.95, 1044.29, 655.32, 1114.78 };
inline constexpr Box DISPLAY = { 233.83, 1176.36, 666.17, 1253.25 };
inline constexpr Box DISPLAY_BEZEL = { 223.05, 1165.17, 676.95, 1263.99 };

}  // namespace photo
}  // namespace ui
}  // namespace tocake
