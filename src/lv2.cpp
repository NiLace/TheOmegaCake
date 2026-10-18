// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) The NL Sounds contributors
//
// LV2 wrapper for The Omega Cake: 45 models of the four-stage fuzz behind one selector.
//
// It owns one engine per model slot (built in instantiate(), so switching does
// not allocate), maps the ports to each engine's knob API (tapers, the tone
// switches), runs the active engine through the oversampler at the factor
// src/os_policy.h gives the slot (the HQ switch raises every slot to 4x), and
// handles the knob glide and the footswitch with a latency-matched dry path. Signal convention:
// the audio port float is volts at the pedal's input jack (the models' native
// unit); the output is volts at its output jack, times the slot's makeup.

#include "omega_cake.h"
#include "cake_variant.h"
#include "variant_registry_generated.h"
#include "os_policy.h"              // oversampling factor per slot and the HQ switch
#include "op_cake.h"              // op-amp ("IC") version -- a separate DSP path
#include "pyramid_master.h"       // Pyramid Master PM-3 -- bespoke topology, separate DSP
#include "coloured_goat.h"        // V2 Coloured Goat -- clean first stage + C14 tone coupling
#include "distortion_cake.h"      // V2 Distortion & Cake -- standard circuit + C14 clip2->tone
#include "simple_tone_cake.h"     // Vm Precision / German Cake -- simplified tone + emitter bypass
#include "dreamy_cake.h"           // Vm Dreamy Cake -- grounded-emitter clone
#include "jumbo_goat.h"           // V2 Jumbo Goat -- 6-transistor dual-path + germanium gate
#include "jumbo_open_goat.h"      // V2 Jumbo Open Goat -- Jumbo Goat without the gate
#include "goat_s_od.h"            // V2 Goat&S OD -- a Japanese clone (active tone)
#include "dirty_goat.h"           // V2 Dirty Goat -- an early clone (tone-bottom cap)
#include "raw_cake.h"        // V3 Raw Cake -- V3 tone-bypass board (MPSA18 + switch)
#include "reversed_cake.h"        // V6 Reversed Cake -- V6 tone-bypass board (2N5088 + switch)
#include "jumbo_colour_cake.h"   // V3 Jumbo Colour Cake -- a British three-transistor fuzz
#include "lost_pumpkin.h"        // V4-5 Lost Pumpkin -- the two-dual-op-amp version (its own engine)
#include "output_makeup_generated.h"  // tocake::output_makeup_for (MEASURED level table)
#include "tone_level_generated.h"     // tocake::tone_level_db_for (MEASURED level by tone)
#include "denormal_guard.h"
#include "isa_dispatch.h"
#include "oversampler.h"

#include <lv2/core/lv2.h>

#include <algorithm>
#include <new>
#include <cmath>
#include <cstdint>
#include <cstring>

#define TOCAKE_URI "https://nylarea.com/plugins/omegacake"

namespace {

enum PortIndex {
    PORT_IN        = 0,
    PORT_OUT       = 1,
    PORT_SUSTAIN   = 2,
    PORT_TONE      = 3,
    PORT_VOLUME    = 4,
    PORT_HQ         = 5,  // HQ switch: 4x everywhere (os_policy.h)
    PORT_LATENCY    = 6,  // output: reported latency (samples), for host PDC
    PORT_VARIANT    = 7,  // model selector. Index = spreadsheet tab position; the
                          // label/order/kind come from variant_registry_generated.h
                          // (see VARIANT_SLOTS).
    PORT_TONE_BYPASS = 8, // Tone-Bypass switch (0/1): takes the tone stack out on
                          // the models that had one (the op-amp Pumpkin and
                          // Guru, Raw Cake, Reversed). Ignored by the others.
    PORT_ENABLED    = 9,  // footswitch / bypass: 1 = engaged (process),
                          // 0 = bypassed (dry passthrough). Default 1.
    PORT_TONE_SWITCH = 10,// the Little Pumpkins' two-way tone switch (0 = the
                          // darker pot end, 1 = the brighter). Ignored by the
                          // others, which have a tone knob.
};

// The selector roster (names, order, kinds) is generated from the variant
// spreadsheet into tocake::VARIANT_SLOTS -- index = tab position, label = tab name. The
// VariantKind picks the engine class (TOCAKE_ENGINES below: the Standard DK
// engine for most slots, a circuit of its own for the others, OpCake for the
// op-amp models); build_dsp() switches on it.

// Snap the model-selector port to a valid variant index. Order and count are
// whatever the spreadsheet's tab layout produced (tocake::VARIANT_COUNT); the default
// (.ttl) is the index of "V1 The Pyramid" (tab 0).
// A non-finite value keeps `current`; the value is clamped BEFORE converting,
// because a float beyond int's range converts to an undefined value.
inline int snap_variant(float v, int current) {
    if (!std::isfinite(v)) return current;
    const float hi = static_cast<float>(tocake::VARIANT_COUNT - 1);
    return static_cast<int>(std::lround(v < 0.0f ? 0.0f : (v > hi ? hi : v)));
}

// (Variant selection is table-driven from tocake::VARIANT_SLOTS -- see build_dsp.)

// The input stage (trim, level meter, pickup emulation) is not the pedal's: it
// lives in the host chain. The plugin takes the input as volts at its jack.
// An input sample beyond INPUT_MAX is no signal a guitar or a line can carry;
// it is treated like a non-finite one, as silence.
constexpr float    INPUT_MAX = 1.0e4f;
inline float sanitise_input(float x) noexcept { return (std::isfinite(x) && std::fabs(x) < INPUT_MAX) ? x : 0.0f; }

// One engine per VariantKind; exactly one pointer is non-null. The Plugin
// derives from this, so p->dsp / p->is_opamp name the ACTIVE engine; the bank
// below holds one set per slot, built at instantiate().
struct EngineView {
    tocake::OmegaCake* dsp       = nullptr;   // VariantKind::Standard slots (the transistor pedals)
    tocake::OpCake*  op_dsp    = nullptr;   // the op-amp slots (Pumpkin, Little Pumpkin, Guru)
    tocake::PyramidMaster* pm3_dsp = nullptr; // VariantKind::PM3 slot ("V1 Master Pyramid")
    tocake::ColouredGoat*  cg_dsp  = nullptr; // VariantKind::ColouredGoat slot ("V2 Coloured Goat")
    tocake::DistortionCake* dc_dsp = nullptr; // VariantKind::DistortionCake slot ("V2 Distortion & Cake")
    tocake::SimpleToneCake* st_dsp = nullptr; // VariantKind::SimpleToneCake slots ("Vm Precision/German Cake")
    tocake::DreamyCake*    dr_dsp = nullptr; // VariantKind::DreamyCake slot ("Vm Dreamy Cake")
    tocake::JumboGoat*     jg_dsp  = nullptr; // VariantKind::JumboGoat slot ("V2 Jumbo Goat")
    tocake::JumboOpenGoat* jog_dsp = nullptr; // VariantKind::JumboOpenGoat slot ("V2 Jumbo Open Goat")
    tocake::GoatSOD*       gsod_dsp = nullptr; // VariantKind::GoatSOD slot ("V2 Goat&S OD")
    tocake::DirtyGoat*     dg_dsp   = nullptr; // VariantKind::DirtyGoat slot ("V2 Dirty Goat")
    tocake::RawCake*   rmc_dsp  = nullptr; // VariantKind::RawCake slot ("V3 Raw Cake")
    tocake::ReversedCake*  rc_dsp   = nullptr; // VariantKind::ReversedCake slot ("V6 Reversed Cake")
    tocake::JumboColourCake* jcc_dsp = nullptr; // VariantKind::JumboColourCake slot ("V3 Jumbo Colour Cake")
    tocake::LostPumpkin*   lp_dsp   = nullptr; // VariantKind::LostPumpkin slot ("V4-5 Lost Pumpkin")
    bool             is_opamp  = false;     // op-amp DSP active (Pumpkin or Little)
    bool             is_little = false;     // Little Pumpkin (fixed fuzz/tone) active
    bool             is_pm3    = false;     // Pyramid Master DSP active
    bool             is_cg     = false;     // Coloured Goat DSP active
    bool             is_dc     = false;     // Distortion & Cake DSP active
    bool             is_st     = false;     // Simplified-tone circuit (Precision/German) active
    bool             is_dreamy = false;     // Dreamy Cake (grounded-emitter circuit) active
    bool             is_jg     = false;     // Jumbo Goat DSP active
    bool             is_jog    = false;     // Jumbo Open Goat DSP active
    bool             is_gsod   = false;     // Goat&S OD DSP active
    bool             is_dg     = false;     // Dirty Goat DSP active
    bool             is_rmc    = false;     // Raw Cake DSP active
    bool             is_rc     = false;     // Reversed Cake DSP active
    bool             is_jcc    = false;     // Jumbo Colour Cake DSP active
    bool             is_lp     = false;     // Lost Pumpkin DSP active
};

// Every slot's engine is constructed in instantiate(), twice: `factory` is
// never run and keeps the freshly constructed state; `live` is the one run()
// drives. Selecting a slot copies factory -> live and re-prepares: the effect
// of delete/new/prepare/reset without the heap. run() must not allocate
// (lv2:hardRTCapable), and a model change happens inside run().
struct SlotEngines { EngineView live, factory; };

struct Plugin : EngineView {
    SlotEngines bank[tocake::VARIANT_COUNT];
    const float*    in        = nullptr;
    float*          out       = nullptr;
    const float*    sustain   = nullptr;
    const float*    tone      = nullptr;
    const float*    volume    = nullptr;
    const float*    hq_port    = nullptr;
    float*          latency    = nullptr;
    const float*    variant_port = nullptr;
    const float*    tone_bypass  = nullptr;   // Tone-Bypass switch (PORT_TONE_BYPASS)
    const float*    tone_switch  = nullptr;   // Little Pumpkin tone switch (PORT_TONE_SWITCH)
    const float*    enabled      = nullptr;   // footswitch: 1 = on, 0 = bypass

    tocake::Oversampler os;
    bool            cur_hq     = false;  // HQ switch as last applied (os_policy.h)
    int             cur_variant = tocake::DEFAULT_VARIANT;   // active model index (the TTL's default)
    float           out_makeup = 1.0f;   // per-variant post-model output gain
    // A post-model level that follows the glided TONE knob (tone_level_generated.h; only the models
    // that have a table: the Lost Pumpkin). Smoothed per sample by a one-pole of
    // kToneGainTau toward the value at the knob's current glide step, which moves on the stream grid:
    // no steps when the knob turns, and the output does not depend on the host's block size. A model
    // change starts it fresh.
    static constexpr double kToneGainTau = 0.005;   // s
    const float*    tone_level = nullptr;
    float           tone_gain  = 1.0f;
    float           tone_gain_a = 0.0f;             // per-sample coefficient (set_rate)
    bool            tone_gain_snap = true;

    // Knob glide (see knob_tick). knob_snap: the next tick jumps to the
    // ports (after instantiate, activate and a model change: nothing to glide from).
    static constexpr uint32_t kKnobStep = 16;          // samples between applications
#ifndef TOCAKE_KNOB_TAU_MS
#define TOCAKE_KNOB_TAU_MS 20
#endif
    static constexpr double   kKnobTau  = TOCAKE_KNOB_TAU_MS * 1e-3;   // s, one-pole glide
    static constexpr float    kKnobSnap = 1.0e-4f;     // land exactly within this
    float           knob_now[3] = { 0.8f, 0.5f, 0.6f };  // sustain, tone, volume: the TTL defaults
    uint32_t        knob_countdown = 0;
    bool            knob_snap = true;
    double          knob_a = 0.0;     // glide coefficient per tick at sample_rate (set_rate)

    double          sample_rate = 48000.0;

    // The host rate is fixed for the instance's life: everything derived from
    // it alone is computed here, once.
    void set_rate(double rate) {
        sample_rate = rate;
        knob_a = 1.0 - std::exp(-static_cast<double>(kKnobStep) / (kKnobTau * rate));
        tone_gain_a = static_cast<float>(1.0 - std::exp(-1.0 / (kToneGainTau * rate)));
    }

    // ── THE FOOTSWITCH ────────────────────────────────────────────────
    // 1. The declared latency does NOT change with the footswitch. The host
    //    compensates the latency it was told and does not re-ask on `enabled`,
    //    so a bypass passing `in` straight through would come out that many
    //    samples EARLY: a comb filter in any parallel chain, and a jump on the
    //    switch. The dry path is delayed by exactly the published latency.
    // 2. No hard switching: wet <-> dry fades over kFadeSeconds, and while it
    //    fades the engine keeps running (the wet being faded must exist). Fully
    //    bypassed, the engine stops -- that is the footswitch's CPU saving --
    //    and its state stays frozen, like a real pedal with its signal cut.
    // The first span after activate() does not fade: it jumps to the port.
    static constexpr uint32_t kSpan    = 1024;     // run() is cut into spans of at most this
    static constexpr uint32_t kDryRing = 2048;     // power of two, >= kSpan + any latency
    static constexpr uint32_t kDryMask = kDryRing - 1;
    // In SECONDS: a fade tuned in samples is a different fade at every rate.
    static constexpr double   kFadeSeconds = 512.0 / 48000.0;   // ~10.7 ms
#ifndef TOCAKE_PREROLL_MS
#define TOCAKE_PREROLL_MS 10.0
#endif
    static constexpr double   kPrerollSeconds = TOCAKE_PREROLL_MS / 1000.0;   // engine warm-up before the fade-in
    // The ring must hold a whole span plus the delay, or the dry of a span
    // would read samples the same span already overwrote -- with no error.
    static_assert(tocake::Oversampler::TPP - 1 + kSpan <= kDryRing,
                  "the dry ring does not cover a span plus the oversampler's latency");
    float           dry[kDryRing] = {};
    uint32_t        dry_w      = 0;
    double          mix        = 1.0;    // 1 = wet, 0 = dry
    bool            fade_start = true;   // first span after activate(): jump, do not fade
    // The engine is FROZEN while bypassed (the footswitch's CPU saving), so on re-engage it resumes from a state
    // half a second old and its own transient would land in the fade-in. The PRE-ROLL runs the engine for this many samples with the output still dry, then
    // the fade starts. Counted per SAMPLE, so it is block-size independent. Armed while bypassed, spent on re-engage.
    uint32_t        preroll    = 0;
};

// Per-variant OUTPUT makeup gain (linear), applied post-model in the plugin
// layer -- NOT a circuit change. With Volume at full every model has the same
// median loudness over its sustain x tone grid on a real guitar DI, and none peaks above
// -2 dBTP at any sustain/tone/switch setting, so switching models does not jump
// levels (the real pedals differ widely; e.g. V1 Vintage Pyramid's authentic
// R6=1k output stage is far quieter than the rest). The circuits and the
// SPICE-validated models are untouched; this only scales the post-clipping
// output, so the clipping character of each variant is unchanged.
// Table is MEASURED (tocake::output_makeup_for); regenerate with
// `scripts/gen_output_makeup.py --write`.

// The engine kinds, one row each: (type, pointer, active flag). Every place that
// has to touch "whichever engine is live" goes through this list.
#define TOCAKE_ENGINES(X)                                   \
    X(tocake::OpCake,          op_dsp,   is_opamp)          \
    X(tocake::PyramidMaster,   pm3_dsp,  is_pm3)            \
    X(tocake::ColouredGoat,    cg_dsp,   is_cg)             \
    X(tocake::DistortionCake,  dc_dsp,   is_dc)             \
    X(tocake::SimpleToneCake,  st_dsp,   is_st)             \
    X(tocake::DreamyCake,      dr_dsp,   is_dreamy)         \
    X(tocake::JumboGoat,       jg_dsp,   is_jg)             \
    X(tocake::JumboOpenGoat,   jog_dsp,  is_jog)            \
    X(tocake::GoatSOD,         gsod_dsp, is_gsod)           \
    X(tocake::DirtyGoat,       dg_dsp,   is_dg)             \
    X(tocake::RawCake,     rmc_dsp,  is_rmc)            \
    X(tocake::ReversedCake,    rc_dsp,   is_rc)             \
    X(tocake::JumboColourCake, jcc_dsp,  is_jcc)          \
    X(tocake::LostPumpkin,     lp_dsp,   is_lp)

// Construct slot `idx`'s engine into `e`. ALLOCATES: instantiate() only.
void construct_engine(EngineView& e, int idx)
{
    const tocake::VariantSlot& slot = tocake::VARIANT_SLOTS[idx];
    const bool little = (slot.kind == tocake::VariantKind::OpAmpLittle);
    const bool guru   = (slot.kind == tocake::VariantKind::OpAmpGuru);
    const bool opamp  = (slot.kind == tocake::VariantKind::OpAmp) || little || guru;
    e = EngineView{};
    if (opamp) {
        e.op_dsp = new tocake::OpCake{little ? tocake::TheLittlePumpkin
                                      : guru ? tocake::TheGuruCake : tocake::TheOpCake};
        e.is_opamp = true;
        e.is_little = little;
    } else if (slot.kind == tocake::VariantKind::PM3) {
        e.pm3_dsp = new tocake::PyramidMaster{*slot.var};      e.is_pm3 = true;
    } else if (slot.kind == tocake::VariantKind::ColouredGoat) {
        e.cg_dsp = new tocake::ColouredGoat{*slot.var};        e.is_cg = true;
    } else if (slot.kind == tocake::VariantKind::DistortionCake) {
        e.dc_dsp = new tocake::DistortionCake{*slot.var};      e.is_dc = true;
    } else if (slot.kind == tocake::VariantKind::SimpleToneCake) {
        e.st_dsp = new tocake::SimpleToneCake{*slot.var};      e.is_st = true;
    } else if (slot.kind == tocake::VariantKind::DreamyCake) {
        e.dr_dsp = new tocake::DreamyCake{*slot.var};          e.is_dreamy = true;
    } else if (slot.kind == tocake::VariantKind::JumboGoat) {
        e.jg_dsp = new tocake::JumboGoat{*slot.var};           e.is_jg = true;
    } else if (slot.kind == tocake::VariantKind::JumboOpenGoat) {
        e.jog_dsp = new tocake::JumboOpenGoat{*slot.var};      e.is_jog = true;
    } else if (slot.kind == tocake::VariantKind::GoatSOD) {
        e.gsod_dsp = new tocake::GoatSOD{*slot.var};           e.is_gsod = true;
    } else if (slot.kind == tocake::VariantKind::DirtyGoat) {
        e.dg_dsp = new tocake::DirtyGoat{*slot.var};           e.is_dg = true;
    } else if (slot.kind == tocake::VariantKind::RawCake) {
        e.rmc_dsp = new tocake::RawCake{*slot.var};        e.is_rmc = true;
    } else if (slot.kind == tocake::VariantKind::ReversedCake) {
        e.rc_dsp = new tocake::ReversedCake{*slot.var};        e.is_rc = true;
    } else if (slot.kind == tocake::VariantKind::JumboColourCake) {
        e.jcc_dsp = new tocake::JumboColourCake{*slot.var};    e.is_jcc = true;
    } else if (slot.kind == tocake::VariantKind::LostPumpkin) {
        e.lp_dsp = new tocake::LostPumpkin{tocake::TheLostPumpkin}; e.is_lp = true;
    } else {
        e.dsp = new tocake::OmegaCake{*slot.var};
    }
}

void destroy_engine(EngineView& e)
{
#define TOCAKE_DELETE(T, ptr, flag) delete e.ptr; e.ptr = nullptr;
    TOCAKE_ENGINES(TOCAKE_DELETE)
#undef TOCAKE_DELETE
    delete e.dsp; e.dsp = nullptr;
}

// Return `live` to the state `factory` was constructed in: destroy it and
// copy-construct it again IN PLACE. The engines hold only fixed-size members
// (and a reference to their variant, which is why this is not an assignment),
// so a placement copy does not touch the heap.
template <class T>
T* restore_in_place(T* live, const T* factory)
{
    live->~T();
    return ::new (static_cast<void*>(live)) T(*factory);
}

void restore_factory(EngineView& live, const EngineView& factory)
{
#define TOCAKE_RESTORE(T, ptr, flag) \
    if (live.ptr) live.ptr = restore_in_place(live.ptr, factory.ptr);
    TOCAKE_ENGINES(TOCAKE_RESTORE)
#undef TOCAKE_RESTORE
    if (live.dsp) live.dsp = restore_in_place(live.dsp, factory.dsp);
}

void apply_knobs(Plugin* p, float sus, float tone, float vol);

inline float knob_port(const float* port, float held);

// Make slot `idx` the active engine, prepared at the current (effective)
// oversampled rate. The same sequence as constructing it afresh -- factory
// state, prepare, reset (a click, fine on a deliberate switch) -- but on
// engines built in instantiate(), so it does not allocate.
void build_dsp(Plugin* p, int idx)
{
    const tocake::VariantSlot& slot = tocake::VARIANT_SLOTS[idx];
    p->out_makeup = tocake::output_makeup_for(slot.ident);
    p->knob_snap = true;          // a fresh engine: no glide from the old knobs
    p->tone_level = tocake::tone_level_db_for(slot.ident);
    p->tone_gain_snap = true;
    p->knob_countdown = 0;
    SlotEngines& s = p->bank[idx];
    restore_factory(s.live, s.factory);
    static_cast<EngineView&>(*p) = s.live;
    // The knobs BEFORE reset(): the rest point depends on them (the transistor
    // circuits' tone network carries DC from Q2's collector to R5), so a slot
    // entered with its knobs elsewhere than the factory 0.5 would otherwise
    // start with a bias transient that an instance already there does not have.
    for (int i = 0; i < 3; ++i) p->knob_now[i] = knob_port(i == 0 ? p->sustain : i == 1 ? p->tone : p->volume, p->knob_now[i]);
    apply_knobs(p, p->knob_now[0], p->knob_now[1], p->knob_now[2]);
    const int eff = tocake::slot_factor(slot.kind, p->cur_hq);
    p->os.set_factor(eff);
    const double sr = p->sample_rate * eff;
    if (p->dsp) { p->dsp->prepare(sr); p->dsp->reset(); }
#define TOCAKE_PREPARE(T, ptr, flag) if (p->ptr) { p->ptr->prepare(sr); p->ptr->reset(); }
    TOCAKE_ENGINES(TOCAKE_PREPARE)
#undef TOCAKE_PREPARE
}

LV2_Handle instantiate(const LV2_Descriptor*, double rate,
                       const char*, const LV2_Feature* const*)
{
    auto* p = new Plugin();
    p->set_rate(rate);
    p->cur_hq      = false;                   // HQ off: each slot's default factor
    p->os.prepare(Plugin::kSpan);             // run() hands the engine at most one span
    p->os.set_factor(2);                      // build_dsp sets the slot's own
    for (int i = 0; i < tocake::VARIANT_COUNT; ++i) {
        construct_engine(p->bank[i].live, i);
        construct_engine(p->bank[i].factory, i);
    }
    build_dsp(p, p->cur_variant);             // model runs at the oversampled rate
    return p;
}

void connect_port(LV2_Handle handle, uint32_t port, void* data)
{
    auto* p = static_cast<Plugin*>(handle);
    switch (port) {
        case PORT_IN:      p->in      = static_cast<const float*>(data); break;
        case PORT_OUT:     p->out     = static_cast<float*>(data);       break;
        case PORT_SUSTAIN:   p->sustain   = static_cast<const float*>(data); break;
        case PORT_TONE:      p->tone      = static_cast<const float*>(data); break;
        case PORT_VOLUME:    p->volume    = static_cast<const float*>(data); break;
        case PORT_HQ:         p->hq_port   = static_cast<const float*>(data); break;
        case PORT_LATENCY:    p->latency   = static_cast<float*>(data);       break;
        case PORT_VARIANT:    p->variant_port = static_cast<const float*>(data); break;
        case PORT_TONE_BYPASS: p->tone_bypass = static_cast<const float*>(data); break;
        case PORT_ENABLED:     p->enabled     = static_cast<const float*>(data); break;
        case PORT_TONE_SWITCH: p->tone_switch = static_cast<const float*>(data); break;
        default: break;
    }
}

void activate(LV2_Handle handle)
{
    auto* p = static_cast<Plugin*>(handle);
    // The ports' knobs BEFORE reset(), as build_dsp() does: the rest point
    // depends on them, and a session restored with the tone away from its
    // factory position would otherwise start with a bias thump.
    for (int i = 0; i < 3; ++i) p->knob_now[i] = knob_port(i == 0 ? p->sustain : i == 1 ? p->tone : p->volume, p->knob_now[i]);
    apply_knobs(p, p->knob_now[0], p->knob_now[1], p->knob_now[2]);
    // The active engine, whichever class it is (TOCAKE_ENGINES; else Standard).
#define TOCAKE_RESET(T, ptr, flag) if (p->flag) p->ptr->reset(); else
    TOCAKE_ENGINES(TOCAKE_RESET) p->dsp->reset();
#undef TOCAKE_RESET
    p->os.reset();
    std::fill(std::begin(p->dry), std::end(p->dry), 0.0f);
    p->dry_w      = 0;
    p->fade_start = true;
    p->preroll    = 0;
    p->knob_snap  = true;
    p->knob_countdown = 0;
}

// Model and oversampling-factor ports. Applied at the top of EVERY span, engaged
// or bypassed: the bypassed dry path is delayed by the oversampler's latency
// and the latency port reports it, so a factor change made while bypassed has
// to land then, not when the pedal is next engaged (the host would re-align
// the track on every bypass).
void apply_model_and_factor(Plugin* p)
{
    // Model (variant) change: switch to the newly-selected model's
    // pre-built engine, back in its factory state and re-prepared at the
    // current oversampled rate. No allocation; it re-biases (a click), which is
    // fine on a deliberate model switch.
    const int want_variant = p->variant_port ? snap_variant(*p->variant_port, p->cur_variant)
                                             : p->cur_variant;
    if (want_variant != p->cur_variant) {
        p->cur_variant = want_variant;
        build_dsp(p, want_variant);
    }

    // HQ switched: re-tune the resampler and re-prepare the model at the new
    // (oversampled) rate. Rare, user-driven; the re-prepare resets cap states
    // (a click) and re-biases -- acceptable on a quality switch.
    const bool want_hq = p->hq_port ? (*p->hq_port > 0.5f) : p->cur_hq;
    if (want_hq != p->cur_hq) {
        p->cur_hq = want_hq;
        const int eff = tocake::slot_factor(tocake::VARIANT_SLOTS[p->cur_variant].kind, want_hq);
        p->os.set_factor(eff);
        const double sr = p->sample_rate * eff;
#define TOCAKE_REPREPARE(T, ptr, flag) if (p->flag) p->ptr->prepare(sr); else
        TOCAKE_ENGINES(TOCAKE_REPREPARE) p->dsp->prepare(sr);
#undef TOCAKE_REPREPARE
    }

}

// The Little Pumpkins' tone switch.
inline bool read_tone_switch(const Plugin* p)
{
    return p->tone_switch && *p->tone_switch > 0.5f;
}

void apply_knobs(Plugin* p, float sus, float tone, float vol)
{
    // Knob positions (0..1) to each engine's physical wiper fractions: the
    // engines take fractions, the tapers of each model are applied here (or by
    // the OmegaCake wrapper itself for the Standard slots).
    if (p->is_opamp) {
        // Op-amp circuit tapers (the real pedal's pots): LINEAR Sustain (the
        // schematic's "B10K SUSTAIN" -- B = linear), log Volume, linear Tone.
        // Taper maps knob->wiper fraction above the SPICE-validated solve, so
        // fidelity is unaffected.
        using tocake::apply_taper; using tocake::PotTaper;
        if (p->is_little) {
            // Little Pumpkin: Fuzz FIXED at max, Tone = 2-position switch at the
            // tone-pot extremes (its own port), Volume the only live knob. No
            // Sustain/Tone pots on the real pedal.
            p->op_dsp->set_sustain(1.0);                       // fuzz pinned at max
            const bool tsw = read_tone_switch(p);
            p->op_dsp->set_tone(tsw ? 1.0 : 0.0);             // 2-way tone switch (extremes)
            p->op_dsp->set_bypass(false);                      // tone stack always in
            if (p->volume) p->op_dsp->set_volume(apply_taper(vol, PotTaper::Log));
        } else {
            if (p->sustain) p->op_dsp->set_sustain(apply_taper(sus, PotTaper::Linear));
            if (p->tone)    p->op_dsp->set_tone   (static_cast<double>(tone));
            if (p->volume)  p->op_dsp->set_volume (apply_taper(vol, PotTaper::Log));
            if (p->tone_bypass) p->op_dsp->set_bypass(*p->tone_bypass > 0.5f);
        }
    } else if (p->is_pm3) {
        // PM-3 (PyramidMaster) takes physical wiper fractions, so apply the
        // variant's tapers here (the shared OmegaCake wrapper does this itself).
        using tocake::apply_taper;
        const auto& v = *tocake::VARIANT_SLOTS[p->cur_variant].var;
        if (p->sustain) p->pm3_dsp->set_sustain(apply_taper(sus, v.taper_sustain));
        if (p->tone)    p->pm3_dsp->set_tone   (apply_taper(tone,    v.taper_tone));
        if (p->volume)  p->pm3_dsp->set_volume (apply_taper(vol,  v.taper_volume));
    } else if (p->is_cg) {
        // Coloured Goat (ColouredGoat) also takes physical wiper fractions, so
        // apply the variant's tapers here (as for PM-3).
        using tocake::apply_taper;
        const auto& v = *tocake::VARIANT_SLOTS[p->cur_variant].var;
        if (p->sustain) p->cg_dsp->set_sustain(apply_taper(sus, v.taper_sustain));
        if (p->tone)    p->cg_dsp->set_tone   (apply_taper(tone,    v.taper_tone));
        if (p->volume)  p->cg_dsp->set_volume (apply_taper(vol,  v.taper_volume));
    } else if (p->is_dc) {
        // Distortion & Cake (bespoke: standard circuit + C14) -- physical wipers.
        using tocake::apply_taper;
        const auto& v = *tocake::VARIANT_SLOTS[p->cur_variant].var;
        if (p->sustain) p->dc_dsp->set_sustain(apply_taper(sus, v.taper_sustain));
        if (p->tone)    p->dc_dsp->set_tone   (apply_taper(tone,    v.taper_tone));
        if (p->volume)  p->dc_dsp->set_volume (apply_taper(vol,  v.taper_volume));
    } else if (p->is_st) {
        // Simplified-tone circuit (Precision/German) -- physical wipers.
        using tocake::apply_taper;
        const auto& v = *tocake::VARIANT_SLOTS[p->cur_variant].var;
        if (p->sustain) p->st_dsp->set_sustain(apply_taper(sus, v.taper_sustain));
        if (p->tone)    p->st_dsp->set_tone   (apply_taper(tone,    v.taper_tone));
        if (p->volume)  p->st_dsp->set_volume (apply_taper(vol,  v.taper_volume));
    } else if (p->is_dreamy) {
        // Dreamy Cake -- physical wipers. The tone pot is wired in REVERSE on
        // the real pedal: the engine and the SPICE deck take the straight wiper
        // position, so the knob is inverted here to turn darker clockwise, as
        // the circuit does.
        using tocake::apply_taper;
        const auto& v = *tocake::VARIANT_SLOTS[p->cur_variant].var;
        if (p->sustain) p->dr_dsp->set_sustain(apply_taper(sus, v.taper_sustain));
        if (p->tone)    p->dr_dsp->set_tone   (1.0 - apply_taper(tone, v.taper_tone));
        if (p->volume)  p->dr_dsp->set_volume (apply_taper(vol,  v.taper_volume));
    } else if (p->is_jg) {
        // Jumbo Goat (JumboGoat) also takes physical wiper fractions, so apply
        // the variant's tapers here (as for PM-3 / Coloured Goat).
        using tocake::apply_taper;
        const auto& v = *tocake::VARIANT_SLOTS[p->cur_variant].var;
        if (p->sustain) p->jg_dsp->set_sustain(apply_taper(sus, v.taper_sustain));
        if (p->tone)    p->jg_dsp->set_tone   (apply_taper(tone,    v.taper_tone));
        if (p->volume)  p->jg_dsp->set_volume (apply_taper(vol,  v.taper_volume));
    } else if (p->is_jog) {
        // Jumbo Open Goat -- physical wiper fractions, variant tapers (as JG).
        using tocake::apply_taper;
        const auto& v = *tocake::VARIANT_SLOTS[p->cur_variant].var;
        if (p->sustain) p->jog_dsp->set_sustain(apply_taper(sus, v.taper_sustain));
        if (p->tone)    p->jog_dsp->set_tone   (apply_taper(tone,    v.taper_tone));
        if (p->volume)  p->jog_dsp->set_volume (apply_taper(vol,  v.taper_volume));
    } else if (p->is_gsod) {
        // Goat&S OD (GoatSOD) -- physical wiper fractions, variant tapers.
        using tocake::apply_taper;
        const auto& v = *tocake::VARIANT_SLOTS[p->cur_variant].var;
        if (p->sustain) p->gsod_dsp->set_sustain(apply_taper(sus, v.taper_sustain));
        if (p->tone)    p->gsod_dsp->set_tone   (apply_taper(tone,    v.taper_tone));
        if (p->volume)  p->gsod_dsp->set_volume (apply_taper(vol,  v.taper_volume));
    } else if (p->is_dg) {
        // Dirty Goat (DirtyGoat) -- physical wiper fractions, variant tapers.
        using tocake::apply_taper;
        const auto& v = *tocake::VARIANT_SLOTS[p->cur_variant].var;
        if (p->sustain) p->dg_dsp->set_sustain(apply_taper(sus, v.taper_sustain));
        if (p->tone)    p->dg_dsp->set_tone   (apply_taper(tone,    v.taper_tone));
        if (p->volume)  p->dg_dsp->set_volume (apply_taper(vol,  v.taper_volume));
    } else if (p->is_rmc) {
        // Raw Cake (RawCake) -- physical wiper fractions + the tone-bypass
        // switch on the shared Tone-Bypass port (true = flat R28 bypass).
        using tocake::apply_taper;
        const auto& v = *tocake::VARIANT_SLOTS[p->cur_variant].var;
        if (p->sustain) p->rmc_dsp->set_sustain(apply_taper(sus, v.taper_sustain));
        if (p->tone)    p->rmc_dsp->set_tone   (apply_taper(tone,    v.taper_tone));
        if (p->volume)  p->rmc_dsp->set_volume (apply_taper(vol,  v.taper_volume));
        if (p->tone_bypass) p->rmc_dsp->set_bypass(*p->tone_bypass > 0.5f);
    } else if (p->is_rc) {
        using tocake::apply_taper;
        const auto& v = *tocake::VARIANT_SLOTS[p->cur_variant].var;
        if (p->sustain) p->rc_dsp->set_sustain(apply_taper(sus, v.taper_sustain));
        if (p->tone)    p->rc_dsp->set_tone   (apply_taper(tone,    v.taper_tone));
        if (p->volume)  p->rc_dsp->set_volume (apply_taper(vol,  v.taper_volume));
        if (p->tone_bypass) p->rc_dsp->set_bypass(*p->tone_bypass > 0.5f);
    } else if (p->is_jcc) {
        using tocake::apply_taper;
        const auto& v = *tocake::VARIANT_SLOTS[p->cur_variant].var;
        if (p->sustain) p->jcc_dsp->set_sustain(apply_taper(sus, v.taper_sustain));
        // The tone knob turns the OTHER way round on this circuit's capture: its pot lug 1 (the engine's
        // and the deck's wiper fraction 1) sits on the BASS arm, where every other netlist of the pedal
        // has the treble arm. This circuit's Bass-Treble works like the other models' (clockwise =
        // treble), so the knob drives the wiper from the far end. A SPICE null cannot see this: it
        // compares engine and deck at the same wiper.
        if (p->tone)    p->jcc_dsp->set_tone   (apply_taper(1.0f - tone, v.taper_tone));
        if (p->volume)  p->jcc_dsp->set_volume (apply_taper(vol,  v.taper_volume));
    } else if (p->is_lp) {
        // Lost Pumpkin: all three pots are B (linear) on both drawings of the board; the engine takes wiper
        // fractions with 1 = clockwise (most sustain, treble, most volume), as its deck does. No tone-bypass switch.
        if (p->sustain) p->lp_dsp->set_sustain(static_cast<double>(sus));
        if (p->tone)    p->lp_dsp->set_tone   (static_cast<double>(tone));
        if (p->volume)  p->lp_dsp->set_volume (static_cast<double>(vol));
    } else if (tocake::VARIANT_SLOTS[p->cur_variant].fixed_fuzz_tone_switch) {
        // Transistor Little Pumpkin: no sustain pot (pinned; the wiper sits on
        // its node whatever the value), tone a two-way switch at the pot ends on
        // its own port, volume the only knob.
        p->dsp->set_sustain(1.0);
        const bool tsw = read_tone_switch(p);
        p->dsp->set_tone(tsw ? 1.0 : 0.0);
        if (p->volume)  p->dsp->set_volume (static_cast<double>(vol));
    } else {
        if (p->sustain) p->dsp->set_sustain(static_cast<double>(sus));
        if (p->tone)    p->dsp->set_tone   (static_cast<double>(tone));
        if (p->volume)  p->dsp->set_volume (static_cast<double>(vol));
    }

}

// Knob glide: a stepped knob re-tunes the circuit at once, and at a low input
// level nothing hides the step, which clicks. Each knob glides toward its port value as a one-pole of
// kKnobTau, applied every kKnobStep samples (each application rebuilds the
// engine's matrices, only while moving), and lands EXACTLY on the port value
// once within kKnobSnap -- so a settled knob is the same engine as before.

// A knob port read for the glide. A non-finite value (a broken automation lane)
// keeps the knob where it is: taken as a target, it would make the glide state
// NaN for good -- the knob stops responding and the engine keeps rebuilding every
// tick.
inline float knob_port(const float* port, float held)
{
    if (!port || !std::isfinite(*port)) return held;
    return std::clamp(*port, 0.0f, 1.0f);
}

bool knobs_moving(const Plugin* p)
{
    return p->knob_now[0] != knob_port(p->sustain, p->knob_now[0])
        || p->knob_now[1] != knob_port(p->tone,    p->knob_now[1])
        || p->knob_now[2] != knob_port(p->volume,  p->knob_now[2]);
}

void knob_tick(Plugin* p)
{
    const float target[3] = { knob_port(p->sustain, p->knob_now[0]),
                              knob_port(p->tone,    p->knob_now[1]),
                              knob_port(p->volume,  p->knob_now[2]) };
    const double a = p->knob_a;
    for (int i = 0; i < 3; ++i) {
        if (p->knob_snap) { p->knob_now[i] = target[i]; continue; }
        const float d = target[i] - p->knob_now[i];
        p->knob_now[i] = (std::fabs(d) <= Plugin::kKnobSnap)
            ? target[i] : static_cast<float>(p->knob_now[i] + a * d);
    }
    p->knob_snap = false;
    // Still gliding after this tick: the engine may rebuild without refinement;
    // the tick that lands rebuilds exactly (dk_engine.h, set_coarse).
    const bool gliding = p->knob_now[0] != target[0] || p->knob_now[1] != target[1]
                      || p->knob_now[2] != target[2];
    if (p->dsp) p->dsp->set_coarse(gliding);
#define TOCAKE_COARSE(T, ptr, flag) if (p->ptr) tocake::dk_set_coarse(*p->ptr, gliding, 0);
    TOCAKE_ENGINES(TOCAKE_COARSE)
#undef TOCAKE_COARSE
    apply_knobs(p, p->knob_now[0], p->knob_now[1], p->knob_now[2]);
}

// The level table at a tone position: linear in dB between its kToneLevelSteps points.
float tone_level_gain(const float* table_db, float tone)
{
    const float pos = std::clamp(tone, 0.0f, 1.0f) * static_cast<float>(tocake::kToneLevelSteps - 1);
    const int i = std::min(static_cast<int>(pos), tocake::kToneLevelSteps - 2);
    const float db = table_db[i] + (pos - static_cast<float>(i)) * (table_db[i + 1] - table_db[i]);
    return std::pow(10.0f, db / 20.0f);
}

void run_engine(Plugin* p, float* buf, uint32_t n)
{
    // Run the model at the oversampled rate. In-place is safe: within each
    // chunk the upsampler consumes all input before the downsampler writes
    // output. At factor 1 this is an exact per-sample pass-through.
#define TOCAKE_RUN(T, ptr, flag) if (p->flag) { T* d = p->ptr; tocake::run_block(p->os, buf, n, d); } else
    TOCAKE_ENGINES(TOCAKE_RUN) {
#undef TOCAKE_RUN
        tocake::OmegaCake* d = p->dsp;
        tocake::run_block(p->os, buf, n, d);
    }

}

void run_engaged(Plugin* p, const float* in_span, float* out_span, uint32_t n_samples)
{

    // The input into `out` (base rate); the model runs on it afterwards at the
    // oversampled rate. A non-finite or absurd INPUT sample is silence here,
    // before the engine sees it: one NaN from a plugin upstream must not stop the
    // pedal.
    const float* in  = in_span;
    float*       out = out_span;
    for (uint32_t i = 0; i < n_samples; ++i) out[i] = sanitise_input(in[i]);

    // Run the model at the oversampled rate, the knobs gliding on a grid of
    // kKnobStep samples counted along the STREAM (not from the block start), so
    // the output does not depend on the host's block size. A settled knob
    // costs nothing: the span runs in one piece.
    // The switch ports are read here once per span; a tick due at the span start
    // (knob_tick ends in apply_knobs) reads them anyway, so it is not repeated.
    if (p->knob_countdown != 0)
        apply_knobs(p, p->knob_now[0], p->knob_now[1], p->knob_now[2]);
    for (uint32_t done = 0; done < n_samples; ) {
        if (p->knob_countdown == 0) {
            knob_tick(p);
            p->knob_countdown = Plugin::kKnobStep;
        }
        const uint32_t rem = n_samples - done;
        const uint32_t len = knobs_moving(p) ? std::min(p->knob_countdown, rem) : rem;
        run_engine(p, out + done, len);
        if (p->tone_level) {
            const float tgt = p->out_makeup * tone_level_gain(p->tone_level, p->knob_now[1]);
            if (p->tone_gain_snap) { p->tone_gain = tgt; p->tone_gain_snap = false; }
            float g = p->tone_gain;
            for (uint32_t i = 0; i < len; ++i) { g += p->tone_gain_a * (tgt - g); out[done + i] *= g; }
            p->tone_gain = g;
        }
        done += len;
        const uint32_t cd = p->knob_countdown;   // keep the grid's phase along the stream
        p->knob_countdown = (len < cd) ? cd - len
            : (Plugin::kKnobStep - (len - cd) % Plugin::kKnobStep) % Plugin::kKnobStep;
    }

    // Per-variant output makeup (post-model loudness match, see above). A model with a tone-level
    // table had it applied in the loop above, inside its smoothed gain.
    if (!p->tone_level && p->out_makeup != 1.0f)
        for (uint32_t i = 0; i < n_samples; ++i) out[i] *= p->out_makeup;

    if (p->latency) *p->latency = static_cast<float>(p->os.latency_samples());
}

// One span with the footswitch applied (see Plugin::kSpan). An unconnected
// `enabled` port reads as engaged, so the offline tools that never connect it
// run the engaged path unchanged.
void run_span(Plugin* p, const float* in, float* out, uint32_t n)
{
    // The input goes into the ring FIRST: `out` may be the same buffer. Sanitised like the wet path's: the fade
    // computes g * wet + (1 - g) * dry, and (1 - g) * NaN is NaN even at g = 1.
    const uint32_t w0 = p->dry_w;
    for (uint32_t i = 0; i < n; ++i) p->dry[(w0 + i) & Plugin::kDryMask] = sanitise_input(in[i]);
    p->dry_w = (w0 + n) & Plugin::kDryMask;

    apply_model_and_factor(p);
    const bool   engaged = !p->enabled || *p->enabled >= 0.5f;
    const double target  = engaged ? 1.0 : 0.0;
    if (p->fade_start) { p->mix = target; p->fade_start = false; }

    if (p->mix == 1.0 && target == 1.0) {              // engaged: the usual path, untouched
        run_engaged(p, in, out, n);
        return;
    }
    const uint32_t lat = static_cast<uint32_t>(p->os.latency_samples());
    auto dry_at = [p, w0, lat](uint32_t i) noexcept {
        return p->dry[(w0 + i + Plugin::kDryRing - lat) & Plugin::kDryMask];
    };
    if (p->mix == 0.0 && target == 0.0) {              // bypassed: delayed dry, engine frozen
        p->preroll = static_cast<uint32_t>(std::lround(Plugin::kPrerollSeconds * p->sample_rate));
        for (uint32_t i = 0; i < n; ++i) out[i] = dry_at(i);
        if (p->latency) *p->latency = static_cast<float>(lat);
        return;
    }
    // Fading: the engine runs, and the dry is blended in at the same latency.
    run_engaged(p, in, out, n);
    const double step = 1.0 / (Plugin::kFadeSeconds * p->sample_rate);
    for (uint32_t i = 0; i < n; ++i) {
        if (p->preroll > 0 && target > p->mix) { --p->preroll; out[i] = dry_at(i); continue; }   // pre-roll
        p->mix = (target > p->mix) ? std::min(1.0, p->mix + step) : std::max(0.0, p->mix - step);
        // The GAIN follows a smoothstep of the linear position. A linear gain has a corner at each end of the
        // fade, and the corner multiplies the signal's own slope: a click on a low note. The ends are exact
        // (s(0) = 0, s(1) = 1), so engaged and bypassed stay bit-exact.
        const double g = p->mix * p->mix * (3.0 - 2.0 * p->mix);
        out[i] = static_cast<float>(g * out[i] + (1.0 - g) * dry_at(i));
    }
}

void run(LV2_Handle handle, uint32_t n_samples)
{
    // RAII, not a bare set: MXCSR is per-THREAD, so leaving flush-to-zero on
    // would change the arithmetic of every plugin after us in the chain and
    // of the host's own summing -- an undeclared side effect on state we do
    // not own. The guard restores it when the block ends.
    tocake::DenormalGuard denormal_guard;

    auto* p = static_cast<Plugin*>(handle);
    // Written on every call, a zero-length one included (the host may read it right after activate()); a factor
    // change inside a span rewrites it there.
    if (p->latency) *p->latency = static_cast<float>(p->os.latency_samples());
    if (!p->in || !p->out) return;
    for (uint32_t off = 0; off < n_samples; off += Plugin::kSpan) {
        const uint32_t n = std::min(Plugin::kSpan, n_samples - off);
        run_span(p, p->in + off, p->out + off, n);
    }
}

void deactivate(LV2_Handle) {}

void cleanup(LV2_Handle handle)
{
    auto* p = static_cast<Plugin*>(handle);
    for (auto& slot : p->bank) {           // the active view only aliases bank[cur]
        destroy_engine(slot.live);
        destroy_engine(slot.factory);
    }
    delete p;
}

const void* extension_data(const char*) { return nullptr; }

const LV2_Descriptor descriptor = {
    TOCAKE_URI,
    instantiate,
    connect_port,
    activate,
    run,
    deactivate,
    cleanup,
    extension_data,
};

} // namespace

extern "C" LV2_SYMBOL_EXPORT
const LV2_Descriptor* lv2_descriptor(uint32_t index)
{
    return (index == 0) ? &descriptor : nullptr;
}
