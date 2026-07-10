module;

#include <vector>
#include <memory>
#include <cmath>
#include <cstddef>
#include <algorithm>
#include <functional>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>

export module Geo.Strategy;

export import Geo.Chain;

// =====================================================================
// The 2D SLICE STRATEGY -- toolpathing for one planar slice. Prepared
// chains in, final toolpath out, and nothing else in the universe.
//
// The contract:
//
//   * The CALLER prepares the input. Doctrine orientation (CCW bounds the
//     region to clear, CW bounds keep-islands), open-air marks on free
//     edges -- all established before the strategy sees anything. In the
//     Sketch app the user's drawn directions are that preparation; in the
//     CAM app, orientByNesting + markOpenAir condition the mesh section.
//     The strategy neither knows nor cares which.
//
//   * The STRATEGY returns everything: the generation series (profiles),
//     per-step stats and sanity verdicts (steps), and the final toolpath --
//     whose chain order and travel directions ARE the tool's motion, with
//     direction (order) and milling sense (handedness) already resolved.
//
//   * The caller displays / executes the result BLINDLY. Any further
//     processing it performs is strictly presentational.
// =====================================================================
export namespace Geo {

    enum class StrategyKind { Profile, Hatch };

    struct SliceParams {

        StrategyKind kind = StrategyKind::Profile;

        // HATCH only: the direction the raster strokes run, degrees CCW from +x.
        // Scan lines advance perpendicular to it; 0 = strokes along +x (fewest
        // reversals when the long axis of the region is +x).
        float hatchAngle = 0.0f;

        float toolRadius = 1.0f;      // the boundary clearance: generation 1 offsets by exactly this
        float stepover = 1.0f;        // generation advance as a FRACTION of the tool radius
                                      // (gen 2+ offset by toolRadius * stepover; <= 2.0 for full
                                      // coverage, typically 0.5 .. 1.0)
        int maxGenerations = 32;      // recursion cap (outsets grow forever)

        // Finishing pass: a single thin inset taken right after the boundary
        // clearance ring, so the ladder runs 1*R -> finishWidth*R -> stepover*R...
        // -- leaving a very thin pass that a consumer can finish at its own feed.
        bool  finishPass = true;      // insert the thin finishing generation
        float finishWidth = 0.1f;     // its inset, as a FRACTION of the tool radius

        // Emit the OUTER profile (generation 1, the wall-touching ring)? When off,
        // it is generated but NOT emitted/linked/led -- the wall is deferred to a
        // later finishing pass that runs with this on. (Used by the CAM operation's
        // global finishing scheme; per-slice runs leave it on.)
        bool  emitOuter = true;

        // Post-processing axes (independent):
        bool reverse = false;         // chain ORDER: execute the path backwards
        bool climb = true;            // chain HANDEDNESS: climb keeps the method's
                                      // travel (material on the left); conventional
                                      // reverses every chain

        // Lead-in / lead-out: engagement moves woven around RETRACT steps so the
        // tool eases off the cut before retracting and eases back on after arriving,
        // never landing cold on the cut path. Each lead follows a small INSET of the
        // target itself (the same offset method the generations run on), preserving
        // the target's own motion as it approaches; the consumer also reads a lead as
        // a DESCENT, so the inset loop doubles as a helix ramp.
        bool  lead = false;           // weave lead-in/out chains around retract steps
        float leadInset = 0.25f;      // the "safe offset" inset, as a FRACTION of the
                                      // tool radius (the lead rides this far inside the cut)

        // How LONG the lead runs along the inset. The lead is the horizontal run of
        // the ramp the consumer will descend: to drop `cuttingDepth` at `plungeSlope`
        // degrees off horizontal, the run is cuttingDepth / tan(slope). That arc
        // length of inset is what we emit; the CAM app extrudes it into the 3D slope.
        float cuttingDepth = 4.0f;    // depth the lead ramp descends (model units)
        float plungeSlope = 23.0f;    // lead-IN ramp angle off horizontal (deg); shallower = longer
        float retractSlope = 75.0f;   // lead-OUT ramp angle (deg); steep -- it climbs, not plunges
    };

    // One generation's stats and its sanity verdict.
    struct SliceStep {
        size_t profile = 0;           // index into SliceResult::profiles
        size_t chains = 0;
        float area = 0.0f;            // signed area of the whole system at this step
        float cleared = 0.0f;         // |area| removed since the previous step
        float length = 0.0f;          // total cutting length
        bool included = true;         // negative (unbounded-claim) systems are skipped
    };

    struct SliceResult {

        // The generation series: profiles[0] is the prepared seed (after the
        // strategy's own preparation, e.g. the open-air push); profiles[n+1] is
        // profiles[n] offset by one step.
        std::vector<Profile> profiles;

        // One step per generated profile (1..n).
        std::vector<SliceStep> steps;

        // The final toolpath: flattened, ordered, ACTUAL tool motion.
        std::vector<Chain> toolpath;

        // Non-fatal advisories raised while slicing (e.g. a finishing pass that had
        // to be narrowed). The caller may surface or aggregate these; an empty list
        // means a clean run. Strings are kept human-readable on purpose.
        std::vector<std::string> warnings;

        SliceResult() = default;
        SliceResult(SliceResult&&) = default;
        SliceResult& operator=(SliceResult&&) = default;

        void clear() { profiles.clear(); steps.clear(); toolpath.clear(); warnings.clear(); }
    };

    // The strategy MANAGER: configure it (params), feed it (ingest), run it.
    // It stores the ingested seed and the full result, and run() also RETURNS
    // the result -- by reference into the manager's own storage, so the caller
    // may read it in place or std::move the pieces out.
    struct SliceStrategy {

        // Configuration -- set before run().
        SliceParams params;

        // Ingested input: the prepared seed (the caller's half of the contract).
        Profile seed;

        // Stored output: everything run() produced.
        SliceResult result;

        SliceStrategy() = default;
        explicit SliceStrategy(SliceParams p) : params(p) {}

        // Feed the prepared chains. The seed is consumed by run() (the strategy
        // applies its own preparation to it, e.g. the open-air push).
        void ingest(Profile prepared) {
            seed = std::move(prepared);
        }

        // Run the whole pipeline: prepare -> generations -> steps -> toolpath.
        // The result is stored on the manager AND returned.
        SliceResult& run() {

            result.clear();

            // CLOSED CHAINS ONLY: offsetting has no concept of an open chain. Open
            // chains in the seed are dropped at the door (the only open chains in
            // a RESULT are the link segments the strategy itself weaves).
            std::erase_if(seed.chains, [](const Chain& c) { return !c.closed; });

            // The strategy's own preparation: open-air edges pre-offset into the
            // free region. For PROFILE the push funds BOTH sides of generation 1's
            // inset (a wall-backed sliver erodes by 2R) plus margin off the exact-
            // tangency degeneracy -- hence past 2R. For HATCH it is what lets the
            // raster OVERRUN the open sides: the inset that follows brings walls in
            // by a radius while the pushed-out air edges stay proud, so strokes
            // begin and end in open air (no leads, no entries into material).
            if (params.kind == StrategyKind::Profile || params.kind == StrategyKind::Hatch) {
                for (Chain& c : seed.chains) {
                    c = c.withOpenAirOffset(params.toolRadius * 2.1f);
                }
            }

            result.profiles.push_back(std::move(seed));
            seed = Profile();

            if (params.kind == StrategyKind::Hatch) {
                hatchGenerations();
                buildSteps();
                buildHatchToolpath();
            }
            else {
                profileGenerations();
                buildSteps();
                buildToolpath();
            }

            return result;
        }

        // Stages
        //--------------------------------------------------

        // The PROFILE strategy's generations: recursively offset until extinction
        // (or the cap). Generation 1 offsets by EXACTLY one tool radius -- the
        // boundary clearance pass. When a finishing pass is asked for, generation 2
        // is a single THIN inset (finishWidth * R), tagged Finish; every generation
        // after that advances by toolRadius * stepover, the ring advance.
        void profileGenerations() {

            if (std::fabs(params.toolRadius) <= 1e-6f) { return; }

            const size_t cap = static_cast<size_t>(std::max(1, params.maxGenerations));
            const float step = params.toolRadius * std::max(0.05f, params.stepover);
            const float finish = params.toolRadius * std::max(0.01f, params.finishWidth);

            while (result.profiles.size() <= cap) {

                const size_t n = result.profiles.size();   // profiles so far (seed = 1)

                // The thin finishing ring (generation 2). A finishing inset that is
                // too wide for a NARROW region MITOSES -- the ring splits into more
                // chains than its parent -- which destroys the single continuous skin
                // and forces the finishing pass to plough into unroughed stock. So
                // narrow the inset until it no longer splits, and warn that we did.
                if (n == 2 && params.finishPass) {

                    const Profile& parent = result.profiles.back();   // generation 1
                    float f = finish;
                    Profile fin;
                    bool narrowed = false;

                    for (int attempt = 0; attempt < 12; attempt++) {
                        fin = parent.offsetBy(f);
                        if (fin.empty()) { break; }
                        if (fin.chains.size() <= parent.chains.size()) { break; }  // no split
                        f *= 0.7f; narrowed = true;                                 // try narrower
                    }

                    if (fin.empty()) { break; }
                    if (narrowed) {
                        result.warnings.push_back(
                            "Sliced profile was too narrow to create a finishing pass with the "
                            "given parameters. The finishing pass stepover was reduced to account "
                            "for this.");
                    }
                    for (Chain& c : fin.chains) { c.link = LinkKind::Finish; }
                    result.profiles.push_back(std::move(fin));
                    continue;
                }

                const float amount = (n == 1) ? params.toolRadius : step;
                Profile next = result.profiles.back().offsetBy(amount);
                if (next.empty()) { break; }
                result.profiles.push_back(std::move(next));
            }
        }

        // The HATCH strategy is a bounded raster fill -- no generation series, just
        // the tool-centre reachable region, filled with parallel strokes. This step
        // resolves that region: the prepared seed inset by one tool radius (walls in
        // by R, air edges left proud by the earlier open-air push, so strokes overrun
        // the open sides). buildHatchToolpath scan-fills it.
        void hatchGenerations() {

            if (std::fabs(params.toolRadius) <= 1e-6f) { return; }

            Profile region = result.profiles.front().offsetBy(params.toolRadius);
            if (!region.empty()) { result.profiles.push_back(std::move(region)); }
        }

        // The HATCH toolpath: sweep a scan line across the reachable region at the
        // requested angle, spaced by the stepover; each scan line's transversal
        // crossings pair into interior spans (even-odd), and each span is one stroke.
        // Alternate lines flip so the tool serpentines (turns at the end) instead of
        // returning -- the zig-zag that trades climb purity for half the travel, the
        // correct roughing trade for air-bounded regions. Step-overs that stay over
        // cleared floor are Cut links; a hop that crosses the boundary (a hole, a
        // concavity) is a Retract. Straight strokes, so no leads: the region's own
        // open-air overrun is the free lead-in. Post axis `reverse` flips the order;
        // milling sense is the zig-zag itself (climb is a finishing-pass concern).
        void buildHatchToolpath() {

            result.toolpath.clear();
            if (result.profiles.size() < 2) {
                result.warnings.push_back(
                    "The region collapsed under the tool radius -- it is narrower than the "
                    "cutter, so nothing was hatched.");
                return;
            }

            const Profile& region = result.profiles.back();   // the tool-centre reachable region

            const float a = params.hatchAngle * 3.14159265358979f / 180.0f;
            const Pos d(std::cos(a), std::sin(a));    // strokes run along d
            const Pos p(-std::sin(a), std::cos(a));   // scan lines advance along p

            // The clip edges, and the region's extent in the (d, p) hatch frame
            std::vector<const Stoicheion*> clip;
            float dMin = 1e30f, dMax = -1e30f, pMin = 1e30f, pMax = -1e30f;
            for (const Chain& c : region.chains) {
                if (c.edges.empty()) { continue; }
                for (const auto& e : c.edges) { clip.push_back(e.get()); }
                std::vector<Pos> pts; c.tessellate(pts);
                for (const Pos& q : pts) {
                    const float qd = q.dot(d), qp = q.dot(p);
                    dMin = std::min(dMin, qd); dMax = std::max(dMax, qd);
                    pMin = std::min(pMin, qp); pMax = std::max(pMax, qp);
                }
            }
            if (clip.empty() || dMax <= dMin || pMax <= pMin) { return; }

            const float spacing = std::max(params.toolRadius * std::max(0.05f, params.stepover), 1e-3f);
            const float pad = spacing + 1.0f;

            // Each scan line: intersect the region, pair interior spans into strokes
            struct Stroke { Pos a, b; };
            std::vector<std::vector<Stroke>> lines;

            const int n = std::max(1, static_cast<int>(std::ceil((pMax - pMin) / spacing)));
            for (int k = 0; k < n; k++) {

                const float s = pMin + spacing * (static_cast<float>(k) + 0.5f);
                if (s > pMax) { break; }

                const Segment2 scan(d * (dMin - pad) + p * s, d * (dMax + pad) + p * s);

                std::vector<float> ts;
                for (const Stoicheion* e : clip) {
                    std::vector<Pos> hits;
                    Chain::edgeCross(scan, *e, hits);
                    for (const Pos& h : hits) { ts.push_back(h.dot(d)); }
                }
                if (ts.size() < 2) { continue; }

                std::sort(ts.begin(), ts.end());
                std::vector<float> u;   // dedupe coincident crossings (a shared vertex hit twice)
                for (float t : ts) { if (u.empty() || t - u.back() > 1e-4f) { u.push_back(t); } }

                std::vector<Stroke> strokes;
                for (std::size_t i = 0; i + 1 < u.size(); i += 2) {
                    strokes.push_back({ d * u[i] + p * s, d * u[i + 1] + p * s });
                }
                if (!strokes.empty()) { lines.push_back(std::move(strokes)); }
            }

            if (lines.empty()) {
                result.warnings.push_back(
                    "No hatch strokes fit the region at this stepover -- it is narrower than a "
                    "single pass.");
                return;
            }
            if (params.reverse) { std::reverse(lines.begin(), lines.end()); }

            // Does the straight hop a -> b cross the region boundary (=> a retract)?
            auto crosses = [&](Pos from, Pos to) -> bool {
                if ((to - from).pythag() <= 1e-4f) { return false; }
                const Segment2 hop(from, to);
                for (const Stoicheion* e : clip) {
                    std::vector<Pos> hits; Chain::edgeCross(hop, *e, hits);
                    if (!hits.empty()) { return true; }
                }
                return false;
            };

            // A stay-down (Cut) link is only safe over the just-cleared swath: the hop
            // must be SHORT (within a stepover-plus-radius of the last cut) AND cross no
            // wall. This is the profile strategy's reach+clearance rule, mirrored -- a
            // longer in-region hop would slot full-width at depth, so it retracts.
            const float stayBudget = spacing + params.toolRadius;

            bool haveCursor = false;
            Pos cursor;
            for (std::size_t li = 0; li < lines.size(); li++) {

                std::vector<Stroke>& strokes = lines[li];
                if (li % 2 == 1) {   // serpentine: reverse this line's strokes end-for-end
                    std::reverse(strokes.begin(), strokes.end());
                    for (Stroke& st : strokes) { std::swap(st.a, st.b); }
                }

                for (Stroke& st : strokes) {

                    if (haveCursor && (st.a - cursor).pythag() > 1e-4f) {
                        const bool stayDown =
                            (st.a - cursor).pythag() <= stayBudget && !crosses(cursor, st.a);
                        Chain link; link.closed = false; link.id = newId();
                        link.link = stayDown ? LinkKind::Cut : LinkKind::Retract;
                        link.edges.push_back(std::make_unique<Segment2>(cursor, st.a));
                        result.toolpath.push_back(std::move(link));
                    }

                    Chain cut; cut.closed = false; cut.id = newId(); cut.link = LinkKind::None;
                    cut.edges.push_back(std::make_unique<Segment2>(st.a, st.b));
                    result.toolpath.push_back(std::move(cut));

                    cursor = st.b; haveCursor = true;
                }
            }

            // the inset region boundary is itself a real toolpath line: the tool
            // centre riding it lands the cutting edge on the finished walls. the
            // centred raster stops half a stepover short of it, so trace the
            // boundary LAST as a half-feed entrance ring -- the finisher that
            // brings every walled side to size. same link rule as a raster hop:
            // a short in-region step stays down, otherwise it retracts on.
            for (const Chain& c : region.chains) {
                if (c.edges.empty()) { continue; }

                const Pos start = Chain::eStart(*c.edges.front());

                if (haveCursor && (start - cursor).pythag() > 1e-4f) {
                    const bool stayDown =
                        (start - cursor).pythag() <= stayBudget && !crosses(cursor, start);
                    Chain link; link.closed = false; link.id = newId();
                    link.link = stayDown ? LinkKind::Cut : LinkKind::Retract;
                    link.edges.push_back(std::make_unique<Segment2>(cursor, start));
                    result.toolpath.push_back(std::move(link));
                }

                Chain perim = c.clone();
                perim.link = LinkKind::Entrance;
                result.toolpath.push_back(std::move(perim));

                cursor = start; haveCursor = true;   // a closed loop returns to its start
            }
        }

        // Per-generation stats and the FINAL SANITY RULE: a step whose total
        // signed area is negative claims an unbounded region (an outside grown
        // past an inside without touching -- or a lone outer trace); it is
        // recorded, measured, and skipped by the toolpath.
        void buildSteps() {

            result.steps.clear();
            if (result.profiles.empty()) { return; }

            float prevArea = 0.0f;
            for (const Chain& c : result.profiles.front().chains) { prevArea += c.signedArea(); }

            for (size_t g = 1; g < result.profiles.size(); g++) {
                SliceStep s;
                s.profile = g;
                s.chains = result.profiles[g].chains.size();
                for (const Chain& c : result.profiles[g].chains) {
                    s.area += c.signedArea();
                    for (const auto& e : c.edges) { s.length += Chain::edgeLength(*e); }
                }
                s.cleared = std::fabs(prevArea) - std::fabs(s.area);
                s.included = s.area >= 0.0f;
                prevArea = s.area;
                result.steps.push_back(s);
            }
        }

        // Flatten the ancestral forest DEPTH-FIRST: once a chain splits, one
        // child lineage is followed to extinction before backing up -- just far
        // enough to the nearest unaccounted sibling. Skipped generations are
        // traversed but not emitted; an orphan sweep appends anything whose
        // lineage broke.
        //
        // The OUTER profile (generation 1, the boundary-clearance ring whose tool
        // touches the part wall) is EMITTED only when emitOuter is set. With it
        // off, generation 1 is still traversed (its children must be reached) but
        // never emitted -- so the wall is left untouched and nothing links or leads
        // to it. A caller (e.g. a roughing slice) defers the wall to a later
        // finishing pass that runs with emitOuter on.
        std::vector<const Chain*> flattenDepthFirst() const {

            std::vector<const Chain*> out;
            if (result.profiles.size() < 2) { return out; }

            auto included = [&](size_t g) {
                return (g - 1 < result.steps.size()) ? result.steps[g - 1].included : true;
            };

            // Emittable: a valid generation that is not the deferred outer ring.
            auto emit = [&](size_t g) {
                if (!params.emitOuter && g == 1) { return false; }
                return included(g);
            };

            std::unordered_set<const Chain*> visited;

            std::function<void(const Chain&, size_t)> visit = [&](const Chain& c, size_t g) {
                visited.insert(&c);
                if (emit(g)) { out.push_back(&c); }
                if (g + 1 >= result.profiles.size()) { return; }
                for (const Chain& child : result.profiles[g + 1].chains) {
                    if (child.parent == c.id) { visit(child, g + 1); }
                }
            };

            for (const Chain& c : result.profiles[1].chains) { visit(c, 1); }

            for (size_t g = 1; g < result.profiles.size(); g++) {
                if (!emit(g)) { continue; }
                for (const Chain& c : result.profiles[g].chains) {
                    if (!visited.count(&c)) { out.push_back(&c); }
                }
            }

            return out;
        }

        // The final toolpath: pure post-processing on clones. Direction = chain
        // ORDER only; milling = chain HANDEDNESS only. Then the LINKING pass
        // weaves tagged connector chains between consecutive cutting chains.
        //
        // Whether the tool may stay down for a hop is MEASURED, never assumed
        // (ancestry proved nothing: a parent's other children can lie anywhere).
        // A cut link requires BOTH:
        //
        //   * REACH -- a circle of one tool radius at the tool's position must
        //     intersect the next chain: the hop distance to the nearest entry
        //     is at most ONE TOOL RADIUS. That is the whole budget.
        //
        //   * CLEARANCE -- the straight (x, y) hop must not cross ANY of the
        //     source geometry (the prepared seed). Crossing the source means
        //     cutting through a wall; that hop must retract instead.
        //
        // Cut link: the next loop is re-seated to begin at its point NEAREST
        // the tool, keeping the flow of motion. Retract link: with the start
        // free, the next loop is re-seated to the middle of its longest edge
        // (the most stable entry).
        //
        // Links carry no Z -- they are planar segments with a TAG; what retracting
        // means in 3D is the consumer's business. From here on, chain travel IS
        // tool motion.
        void buildToolpath() {

            result.toolpath.clear();

            std::vector<Chain> run;
            for (const Chain* c : flattenDepthFirst()) { run.push_back(c->clone()); }

            if (params.reverse) {
                std::reverse(run.begin(), run.end());
            }

            if (!params.climb) {
                for (Chain& c : run) { c = c.reversed(); }
            }

            // ENTRANCE tagging -- execution-order truth about engagement. A cutting
            // chain run while NEITHER its parent NOR all of its children are already
            // cut has no cleared swath beside it: the tool meets full-width material
            // (a suspected full slot). Outside-in that is each region's wall ring;
            // inside-out, the innermost leaf; a deferred outer (emitOuter off) leaves
            // its child correctly suspect. Only EMITTED chains count as relief -- an
            // unemitted child cleared nothing. A leaf with an uncut parent is always
            // suspect (zero children never vacuously relieve).
            {
                std::unordered_map<Id, size_t> children;      // emitted children per id
                for (const Chain& c : run) { children[c.parent]++; }

                std::unordered_set<Id> cut;
                std::unordered_map<Id, size_t> cutChildren;
                for (Chain& c : run) {
                    const auto kids = children.find(c.id);
                    const bool parentCut   = cut.count(c.parent) > 0;
                    const bool childrenCut = kids != children.end() && cutChildren[c.id] == kids->second;
                    if (!parentCut && !childrenCut) { c.link = LinkKind::Entrance; }
                    cut.insert(c.id);
                    cutChildren[c.parent]++;
                }
            }

            // Does the straight hop a -> b cross any source chain geometry?
            auto crossesSource = [&](Pos a, Pos b) -> bool {

                if ((b - a).pythag() <= 1e-4f) { return false; }
                if (result.profiles.empty()) { return false; }

                Segment2 hop(a, b);

                for (const Chain& c : result.profiles.front().chains) {
                    for (const auto& e : c.edges) {
                        std::vector<Pos> hits;
                        Chain::edgeCross(hop, *e, hits);
                        if (!hits.empty()) { return true; }
                    }
                }

                return false;
            };

            // Offset a single chain by a SIGNED amount through the very same offset
            // method the generations run on (so it preserves topology: a circle stays a
            // circle). Positive insets, negative outsets. Returns the offset chain
            // nearest `near`, with the foot point and its edge. Empty when the chain is
            // too small to take the offset (it collapses).
            auto offsetNear = [&](const Chain& src, Pos near, float amount,
                                  Pos& oPt, size_t& oK) -> std::optional<Chain> {
                Profile p; p.chains.push_back(src.clone());
                Profile off = p.offsetBy(amount);
                const Chain* best = nullptr; float bestD = 1e30f; Pos bestPt; size_t bestK = 0;
                for (const Chain& c : off.chains) {
                    if (!c.closed || c.edges.empty()) { continue; }
                    Pos pt; size_t k = c.nearestPoint(near, pt);
                    float d = (pt - near).pythag();
                    if (d < bestD) { bestD = d; best = &c; bestPt = pt; bestK = k; }
                }
                if (!best) { return {}; }
                oPt = bestPt; oK = bestK;
                return best->clone();
            };

            // A sub-portion of a CLOSED loop by ARC LENGTH, measured along travel:
            // head=true returns the first `len` from the loop's seam (its start);
            // head=false returns the last `len`, ending at the seam. `len` past the
            // perimeter spirals AROUND: the lead wraps the inset as many full times
            // as the run demands (the small-island helix), plus a partial wrap for the
            // remainder. This is how much of the inset topology the lead actually
            // covers.
            auto portion = [&](const Chain& seated, float len, bool head)
                                 -> std::vector<std::unique_ptr<Stoicheion>> {
                std::vector<std::unique_ptr<Stoicheion>> out;
                if (seated.edges.empty() || len <= 1e-4f) { return out; }

                float total = 0.0f;
                for (const auto& e : seated.edges) { total += Chain::edgeLength(*e); }
                if (total <= 1e-4f) { return out; }

                // Split the run into whole laps + a leftover partial lap (capped, so
                // pathological depth/slope can't spawn unbounded geometry).
                const int maxLaps = 256;
                int laps = std::min(maxLaps, static_cast<int>(len / total));
                const float rem = len - static_cast<float>(laps) * total;

                auto wholeLap = [&]() {
                    for (const auto& e : seated.edges) { out.push_back(e->clone()); }
                };

                // One sub-perimeter run: the first `plen` (head) or the last `plen`
                // (tail) of a single seam-seated lap.
                auto partial = [&](float plen, bool h) {
                    if (plen <= 1e-4f) { return; }
                    const float skip = h ? 0.0f : (total - plen);
                    float acc = 0.0f;
                    for (const auto& e : seated.edges) {
                        const float el = Chain::edgeLength(*e);
                        const float endAcc = acc + el;
                        if (h) {
                            if (endAcc <= plen + 1e-4f) { out.push_back(e->clone()); }
                            else {
                                std::unique_ptr<Stoicheion> a, b;
                                Chain::splitEdge(*e, Chain::edgePointAt(*e, (plen - acc) / el), a, b);
                                out.push_back(std::move(a));
                                break;
                            }
                        }
                        else {
                            if (endAcc <= skip + 1e-4f) { /* before the kept region */ }
                            else if (acc >= skip - 1e-4f) { out.push_back(e->clone()); }
                            else {
                                std::unique_ptr<Stoicheion> a, b;
                                Chain::splitEdge(*e, Chain::edgePointAt(*e, (skip - acc) / el), a, b);
                                out.push_back(std::move(b));
                            }
                        }
                        acc = endAcc;
                    }
                };

                // Head leads OUT of the seam: whole laps then the partial. Tail leads
                // INTO the seam: the partial first, then the whole laps -- every lap is
                // seam-to-seam, so the joins are continuous either way.
                if (head) {
                    for (int i = 0; i < laps; i++) { wholeLap(); }
                    partial(rem, true);
                }
                else {
                    partial(rem, false);
                    for (int i = 0; i < laps; i++) { wholeLap(); }
                }
                return out;
            };

            // A MORPHING lead: take `run` of the seam-seated base loop and blend each
            // sampled point toward its mate (nearest foot) on `other`. Both come from
            // the SAME offset family -- the inset and an offset OF the inset -- so they
            // share topology and the in-between never distorts (morphing against the
            // original profile does, since its inset can split into children that no
            // longer correspond to it). The blend ramps from fully on `other` at the
            // seam-far end to fully on the base at the seam. `seamAtEnd` puts the base
            // end last (lead-in, blend 1->0) or first (lead-out, blend 0->1). A null
            // `other` leaves the lead tracing the base in place -- the valid in-place
            // helix. Tessellated to segments: a varying blend of a curve is not a curve.
            auto morph = [&](const Chain& seated, const Chain* other,
                             float run, bool seamAtEnd) -> Chain {
                Chain lead; lead.closed = false; lead.id = newId();   // kind set by the caller
                auto base = portion(seated, run, /*head=*/!seamAtEnd);
                float total = 0.0f;
                for (const auto& e : base) { total += Chain::edgeLength(*e); }
                if (base.empty() || total <= 1e-4f) { return lead; }

                const float step = std::max(0.05f, params.toolRadius * 0.15f);
                std::vector<Pos> pts;
                float accStart = 0.0f;
                for (size_t ei = 0; ei < base.size(); ei++) {
                    const Stoicheion& e = *base[ei];
                    const float el = Chain::edgeLength(e);
                    const int k = std::max(1, static_cast<int>(std::ceil(el / step)));
                    for (int i = (ei == 0 ? 0 : 1); i <= k; i++) {
                        const float f = static_cast<float>(i) / k;
                        const Pos P = Chain::edgePointAt(e, f);                // on the base loop
                        const float g = accStart + el * f;                    // from base start
                        const float blend = (seamAtEnd ? (total - g) : g) / total;   // 1 other .. 0 base
                        Pos M = P;
                        if (other && blend > 1e-4f) {
                            Pos foot; other->nearestPoint(P, foot);           // its mate on `other`
                            M = P + (foot - P) * blend;                       // lerp base -> other
                        }
                        if (pts.empty() || (M - pts.back()).pythag() > 1e-5f) { pts.push_back(M); }
                    }
                    accStart += el;
                }
                for (size_t i = 1; i < pts.size(); i++) {
                    lead.edges.push_back(std::make_unique<Segment2>(pts[i - 1], pts[i]));
                }
                return lead;
            };

            // The two morph curves: the inset `I` (the safe child) and `O`, an OUTWARD
            // offset of `I` that returns ~onto the cut. Both are offsets of the single
            // chain `I`, so they share topology -- the in-between never distorts (which
            // morphing against the original profile does, since its inset can split).
            // The lead's landing curve is `O` (~the cut), so it approaches smoothly --
            // not the inset with a straight stub. Null when the feature can't inset.
            // Conventional milling reversed every run chain before linking, which flips the
            // side a positive offset lands on — sign the inset so it still lands on the
            // cleared side. Keyed on params.climb, NOT per-chain winding: a keep-island ring
            // is legitimately CW in climb mode and its +offset is already the safe side.
            // (Band-aid until leads are rebuilt on doctrine-oriented chains.)
            const float safe = (params.climb ? 1.0f : -1.0f)
                             * params.toolRadius * std::max(0.02f, params.leadInset);
            auto leadCurves = [&](const Chain& target, Pos near,
                                  std::optional<Chain>& inner, std::optional<Chain>& outer,
                                  Pos& outerPt, size_t& outerK) -> bool {
                Pos iPt; size_t iK;
                inner = offsetNear(target, near, +safe, iPt, iK);
                if (!inner) { return false; }
                outer = offsetNear(*inner, near, -safe, outerPt, outerK);   // back ~onto the cut
                if (!outer) {                          // can't outset: land on the inset itself
                    outer = inner->clone();
                    outerPt = iPt; outerK = iK;
                }
                return true;
            };

            // Lead-IN: ride `O` (~the cut) and morph from the deep inset `I` out onto
            // it over `run`, then a tiny merge fixes the residual O->entry offset so the
            // cut still begins exactly at the entry. `startOut` is where the rapid lands.
            auto buildLeadIn = [&](const Chain& target, Pos entry, float run, Pos& startOut) -> Chain {
                std::optional<Chain> inner, outer; Pos oPt; size_t oK;
                if (!leadCurves(target, entry, inner, outer, oPt, oK)) {
                    Chain lead = morph(target, nullptr, run, /*seamAtEnd=*/true);   // in-place helix
                    lead.link = LinkKind::LeadIn;
                    startOut = lead.edges.empty() ? entry : Chain::eStart(*lead.edges.front());
                    return lead;
                }
                Chain Oseated = outer->startedAt(oK, oPt);             // land here, nearest the cut
                Chain lead = morph(Oseated, &*inner, run, /*seamAtEnd=*/true);
                lead.link = LinkKind::LeadIn;

                Pos seamPt = lead.edges.empty() ? oPt : Chain::eEnd(*lead.edges.back());
                if ((entry - seamPt).pythag() > 1e-4f) {              // close the residual O->cut gap
                    lead.edges.push_back(std::make_unique<Segment2>(seamPt, entry));
                }
                startOut = lead.edges.empty() ? entry : Chain::eStart(*lead.edges.front());
                return lead;
            };

            // Lead-OUT: the mirror -- a tiny merge off the cut onto `O`, then morph from
            // `O` inward to the deep inset `I`, leaving the tool deep so the retract
            // lifts AWAY from the wall. `endOut` is where the retract begins.
            auto buildLeadOut = [&](const Chain& source, Pos exitPt, float run, Pos& endOut) -> Chain {
                std::optional<Chain> inner, outer; Pos oPt; size_t oK;
                if (!leadCurves(source, exitPt, inner, outer, oPt, oK)) {
                    Chain lead = morph(source, nullptr, run, /*seamAtEnd=*/false);
                    lead.link = LinkKind::LeadOut;
                    endOut = lead.edges.empty() ? exitPt : Chain::eEnd(*lead.edges.back());
                    return lead;
                }
                Chain Oseated = outer->startedAt(oK, oPt);
                Chain body = morph(Oseated, &*inner, run, /*seamAtEnd=*/false);

                Chain lead; lead.closed = false; lead.id = newId(); lead.link = LinkKind::LeadOut;
                if ((oPt - exitPt).pythag() > 1e-4f) {                // close the residual cut->O gap
                    lead.edges.push_back(std::make_unique<Segment2>(exitPt, oPt));
                }
                for (auto& e : body.edges) { lead.edges.push_back(std::move(e)); }
                endOut = lead.edges.empty() ? exitPt : Chain::eEnd(*lead.edges.back());
                return lead;
            };

            // How far each lead runs: the horizontal run of a ramp that covers
            // `cuttingDepth` at the given angle off horizontal (run = depth / tan).
            // The lead-OUT climbs rather than plunges, so it runs its own, steeper
            // angle -- a much shorter ramp.
            auto runFor = [&](float slopeDeg) -> float {
                if (params.cuttingDepth <= 0.0f) { return 0.0f; }
                const float r = std::clamp(slopeDeg, 1.0f, 89.0f) * 3.14159265358979f / 180.0f;
                return params.cuttingDepth / std::tan(r);
            };
            const float leadRun = runFor(params.plungeSlope);
            const float leadOutRun = runFor(params.retractSlope);

            // The outer boundary rings (generation 1, the wall-touching profile --
            // the part's perimeter AND every hole). Their START/END seam leaves a
            // faint witness bump on the cut wall, so each is OVERSHOT: after closing
            // the loop the tool carries on along the chain by one tool radius,
            // re-cutting the seam clean. Identified by id (clones / re-seats preserve
            // it). Inner rings never witness a wall, so they are left alone.
            std::unordered_set<Id> outerIds;
            if (result.profiles.size() > 1) {
                for (const Chain& c : result.profiles[1].chains) { outerIds.insert(c.id); }
            }

            Pos cursor;
            bool haveCursor = false;
            std::optional<Chain> prevCut;   // the cut just emitted, for its lead-out

            for (Chain& next : run) {

                if (next.edges.empty()) { continue; }

                // Decide the link kind by measurement, and re-seat the entry.
                bool cutLink = false;

                if (next.closed) {

                    if (haveCursor) {

                        Pos entry;
                        size_t k = next.nearestPoint(cursor, entry);

                        const bool reachable =
                            (entry - cursor).pythag() <= params.toolRadius + 1e-4f;

                        if (reachable && !crossesSource(cursor, entry)) {
                            cutLink = true;
                            next = next.startedAt(k, entry);
                        }
                    }

                    if (!cutLink) {
                        size_t k = next.longestEdgeIndex();
                        next = next.startedAt(k, Chain::edgeMidpoint(*next.edges[k]));
                    }
                }

                // The cut begins at `entry`. Leads are woven ONLY at retract steps (a
                // cut-link already rides down into the next ring -- it needs none).
                // Such a step is ALWAYS paired: a lead-OUT eases off the previous cut
                // before the retract lifts, and a lead-IN eases onto this one after it
                // lands -- so the tool moves AWAY from a chain before retracting and
                // TOWARD a chain after arriving, never up or down the cut wall.
                Pos entry = Chain::eStart(*next.edges.front());
                Pos linkTarget = entry;
                const bool doLead = params.lead && next.closed && !cutLink;

                // 1. Lead-OUT of the previous cut, easing the tool off it before the
                //    retract. This advances the cursor to the lift point.
                if (doLead && haveCursor && prevCut) {
                    Pos outEnd;
                    Chain lo = buildLeadOut(*prevCut, cursor, leadOutRun, outEnd);
                    if (!lo.edges.empty()) {
                        result.toolpath.push_back(std::move(lo));
                        cursor = outEnd;
                    }
                }

                // 2. Lead-IN to the target, built first so the link can deliver the
                //    tool to where the lead-in begins (its inset, not the cut). An
                //    ENTRANCE ring's ramp descends into full-width material, so its
                //    run is floored at one full lap of the ring -- the descent is
                //    spread around the whole loop, never plunged down a short arc.
                std::optional<Chain> leadIn;
                if (doLead) {
                    float runLen = leadRun;
                    if (next.link == LinkKind::Entrance) {
                        float perim = 0.0f;
                        for (const auto& e : next.edges) { perim += Chain::edgeLength(*e); }
                        runLen = std::max(runLen, perim);
                    }
                    Pos inStart;
                    Chain li = buildLeadIn(next, entry, runLen, inStart);
                    if (!li.edges.empty()) { linkTarget = inStart; leadIn = std::move(li); }
                }

                // 3. The link (rapid/retract, or stay-down cut) to the link target.
                if (haveCursor && (linkTarget - cursor).pythag() > 1e-4f) {
                    Chain link;
                    link.closed = false;
                    link.id = newId();
                    link.link = cutLink ? LinkKind::Cut : LinkKind::Retract;
                    link.edges.push_back(std::make_unique<Segment2>(cursor, linkTarget));
                    result.toolpath.push_back(std::move(link));
                }

                // 4. The lead-in, then the cut itself.
                if (leadIn) {
                    result.toolpath.push_back(std::move(*leadIn));
                }

                // The CLOSED loop, kept for the next chain's lead-out (which insets
                // it -- offsetting needs a closed chain), BEFORE any overshoot opens
                // it.
                prevCut = next.clone();

                // Witness-mark overshoot: an outer boundary ring carries on past its
                // seam by one tool radius along its own start, re-cutting the join so
                // no bump is left where the loop opened and closed.
                if (next.closed && params.toolRadius > 1e-4f && outerIds.count(next.id)) {
                    for (auto& e : portion(next, params.toolRadius, /*head=*/true)) {
                        next.edges.push_back(std::move(e));
                    }
                    next.closed = false;
                }

                cursor = Chain::eEnd(*next.edges.back());
                haveCursor = true;

                result.toolpath.push_back(std::move(next));
            }

            // The DEGENERATE SWEEP -- the last rung below extinction. A seed chain
            // that was GIVEN open air and still left no first-generation child has
            // measured its own verdict: the region is narrower than the tool even
            // after the push funded both sides of the inset. But that same verdict
            // is a guarantee -- a region thinner than the tool is cleared entirely
            // by ONE pass riding a single radius off its walls, overhanging into
            // the air. So the fallback insets the chain's non-open-air runs
            // BLINDLY (each wall edge one radius left of travel, gaps bridged
            // straight) and emits each run as an OPEN cut chain, capped by
            // straight tangent leads: both ends hang in free air, so the ramps
            // descend and climb over nothing. No-touch stays sovereign: a sweep
            // that would ride within one radius of any OTHER wall (a thin slot
            // walled on both sides) is refused with a warning, not emitted. Only
            // chains that had air qualify -- a childless chain walled all round is
            // genuinely uncuttable by this tool. Like generation 1, sweeps are the
            // wall pass, so they obey emitOuter.
            if (params.kind == StrategyKind::Profile && params.emitOuter &&
                params.toolRadius > 1e-4f && !result.profiles.empty()) {

                const Profile& seed0 = result.profiles.front();

                // Children count as coverage only if generation 1 actually CUT: a
                // gen-1 the sanity rule disqualified (negative, unbounded claim)
                // cleared nothing, so its chains must not veto the sweep.
                auto hasChildren = [&](Id id) {
                    if (result.profiles.size() < 2) { return false; }
                    if (!result.steps.empty() && !result.steps.front().included) { return false; }
                    for (const Chain& g : result.profiles[1].chains) {
                        if (g.parent == id) { return true; }
                    }
                    return false;
                };

                for (const Chain& c : seed0.chains) {

                    if (!c.closed || c.edges.empty()) { continue; }
                    // Regions only, never islands -- by SIGNED AREA: a sliver seed's
                    // tangency cusps make turningSign unreliable, and the sweep exists
                    // precisely for such slivers.
                    if (c.signedArea() <= 0.0f) { continue; }

                    bool anyAir = false, anyWall = false;
                    for (const auto& e : c.edges) {
                        if (e->openAir) { anyAir = true; } else { anyWall = true; }
                    }
                    if (!anyAir || !anyWall) { continue; }
                    if (hasChildren(c.id)) { continue; }           // it generated: the ladder handled it

                    // Wall runs: maximal consecutive non-open-air stretches, scanned
                    // from an open-air edge so no run splits across the seam.
                    const size_t n = c.edges.size();
                    size_t first = 0;
                    while (first < n && !c.edges[first]->openAir) { first++; }

                    std::vector<std::vector<const Stoicheion*>> walls;
                    bool inRun = false;
                    for (size_t k = 0; k < n; k++) {
                        const Stoicheion* e = c.edges[(first + k) % n].get();
                        if (e->openAir) { inRun = false; continue; }
                        if (!inRun) { walls.emplace_back(); inRun = true; }
                        walls.back().push_back(e);
                    }

                    for (const std::vector<const Stoicheion*>& wall : walls) {

                        // The blind inset: the sweep is the wall run carried one tool
                        // radius into the region, gaps bridged straight. A sweep is by
                        // definition a first engagement of uncut material: Entrance.
                        Chain sweep;
                        sweep.closed = false;
                        sweep.id = newId();
                        sweep.parent = c.id;
                        sweep.link = LinkKind::Entrance;
                        for (const Stoicheion* e : wall) {
                            std::unique_ptr<Stoicheion> o = Chain::offsetLeft(*e, params.toolRadius);
                            if (Chain::edgeLength(*o) < 1e-4f) { continue; }
                            if (!sweep.edges.empty()) {
                                Pos a = Chain::eEnd(*sweep.edges.back());
                                Pos b = Chain::eStart(*o);
                                if ((b - a).pythag() > 1e-4f) {
                                    sweep.edges.push_back(std::make_unique<Segment2>(a, b));
                                }
                            }
                            sweep.edges.push_back(std::move(o));
                        }
                        if (sweep.edges.empty()) { continue; }

                        // The no-touch guard, measured: the sweep must clear every
                        // wall it is not itself riding by a full radius.
                        std::vector<const Stoicheion*> others;
                        for (const Chain& oc : seed0.chains) {
                            for (const auto& oe : oc.edges) {
                                if (oe->openAir) { continue; }
                                bool own = false;
                                for (const Stoicheion* w : wall) { if (w == oe.get()) { own = true; break; } }
                                if (!own) { others.push_back(oe.get()); }
                            }
                        }
                        bool gouges = false;
                        for (const auto& e : sweep.edges) {
                            for (int i = 0; i <= 8 && !gouges; i++) {
                                Pos p = Chain::edgePointAt(*e, static_cast<float>(i) / 8.0f);
                                for (const Stoicheion* w : others) {
                                    if (w->distanceTo(p) < params.toolRadius - 1e-3f) { gouges = true; break; }
                                }
                            }
                            if (gouges) { break; }
                        }
                        if (gouges) {
                            result.warnings.push_back(
                                "A region narrower than the tool was left uncut: its single-sweep "
                                "pass would ride into an opposing wall.");
                            continue;
                        }

                        if (!params.climb) { sweep = sweep.reversed(); }

                        Pos s = Chain::eStart(*sweep.edges.front());
                        Pos t = Chain::eEnd(*sweep.edges.back());
                        Pos sd = Chain::eStartDir(*sweep.edges.front());
                        Pos td = Chain::eEndDir(*sweep.edges.back());

                        // Straight tangent caps, both ends in open air.
                        Pos linkTarget = s;
                        std::optional<Chain> leadIn, leadOut;
                        if (params.lead && leadRun > 1e-4f) {
                            Chain li; li.closed = false; li.id = newId(); li.link = LinkKind::LeadIn;
                            li.edges.push_back(std::make_unique<Segment2>(s - sd * leadRun, s));
                            linkTarget = s - sd * leadRun;
                            leadIn = std::move(li);
                        }
                        if (params.lead && leadOutRun > 1e-4f) {
                            Chain lo; lo.closed = false; lo.id = newId(); lo.link = LinkKind::LeadOut;
                            lo.edges.push_back(std::make_unique<Segment2>(t, t + td * leadOutRun));
                            leadOut = std::move(lo);
                        }

                        // Ease off the previous closed cut before retracting to the
                        // sweep, the same pairing the main loop keeps.
                        if (params.lead && haveCursor && prevCut) {
                            Pos outEnd;
                            Chain lo = buildLeadOut(*prevCut, cursor, leadOutRun, outEnd);
                            if (!lo.edges.empty()) {
                                result.toolpath.push_back(std::move(lo));
                                cursor = outEnd;
                            }
                            prevCut.reset();
                        }

                        if (haveCursor && (linkTarget - cursor).pythag() > 1e-4f) {
                            Chain link;
                            link.closed = false;
                            link.id = newId();
                            link.link = LinkKind::Retract;
                            link.edges.push_back(std::make_unique<Segment2>(cursor, linkTarget));
                            result.toolpath.push_back(std::move(link));
                        }

                        if (leadIn) { result.toolpath.push_back(std::move(*leadIn)); }
                        result.toolpath.push_back(std::move(sweep));
                        cursor = t;
                        if (leadOut) {
                            cursor = t + td * leadOutRun;
                            result.toolpath.push_back(std::move(*leadOut));
                        }
                        haveCursor = true;
                    }
                }
            }
        }
    };

    // Convenience: one-shot form for callers with nothing to hold onto.
    inline SliceResult runSliceStrategy(Profile seed, const SliceParams& params) {
        SliceStrategy strategy(params);
        strategy.ingest(std::move(seed));
        strategy.run();
        return std::move(strategy.result);
    }
}
