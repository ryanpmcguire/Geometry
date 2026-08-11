module;

#include <vector>
#include <memory>
#include <string>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <algorithm>
#include <chrono>

export module Geo.Frontier;

export import Geo.Chain;

// Geo::Frontier -- the frontier-extraction method, reconstructed. A parallel, spec-driven reimplementation of the offset / crossing / merge machinery in Geo::Chain, written against the method documents in external/Geometry/FrontierExtraction/Method (01..06). The original was written in discovery mode and hard-codes one question -- "the internal frontier of an inward offset" -- as constants (minimum crossing level, CCW-only inversion trial, auto-detected merge sign). Here the question is a parameter: * frontier selection and the inversion guard both key on the SENSE sigma = sign(amount * h) -- the offset sign RELATIVE to the chain's handedness -- not on either alone (Method/03); * the merge is an explicit winding-rule difference: material is { winding >= 1 } with doctrine orientation asserted upstream, and minuend/subtrahend roles decide which clusters may yield material at all (Method/04, 06); * mitosis (one chain splitting against itself) and merging (resolved chains meeting) are kept distinct: offsetValid vs merge. Stoicheia are IMPORTED, not reimplemented -- and so are the edge-level numerical kernels (offsetLeft, edgeCross, splitEdgeAtPoints, endpoint / tangent probes), reused via Geo::Chain's statics under the alias Kern. That is deliberate: those are pure stoicheion-level functions carrying the hard-won robustness lore (double-precision crossings, tangency fusion, absolute-distance interiority, transversal-kiss rejection), and sharing them means a differential test against the old kernel compares the METHOD, not two transcriptions of the same arithmetic. SELECTION RULE (corrected 2026-07-11; see offsetValid): the frontier extreme keys on sign(amount) -- the side the curve was pushed -- NOT on sigma. A mirror argument forbids sigma-keyed selection (reflection preserves sigma but negates the level field); production anchors pin both halves (legacy engine: min at amount > 0; the crescent seed: max at amount < 0). Sigma gates only the inversion guard. KNOWN LIMITATION: a growing chain's legitimate secondary births (a pinched-off void between two arms) are dropped with the non-extreme levels. Affects efficiency, not safety.
export namespace Geo::Frontier {

    // The proven edge-level kernel host (Geo::Chain's statics).
    using Kern = Geo::Chain;

    // Which side of the difference a chain argues for. Minuend chains bound the region we are keeping (fill); subtrahend chains bound geometry we are removing from it (avoid). Merge drops any cluster with no minuend: subtrahend geometry alone can never yield material (Method/04).
    enum class Role : std::uint8_t { Minuend, Subtrahend };

    // Chain -- an oriented run of stoicheia; travel is meaning, material is on the LEFT of travel (Method/02, 05).
    struct Chain {

        std::vector<std::unique_ptr<Stoicheion>> edges;   // oriented: travel is start -> end
        bool closed = false;

        Id id = 0;
        Id parent = 0;
        Role role = Role::Minuend;

        // An offset arc whose radius lands below this is a collapsed feature (a fillet at the offset distance, within numerical slop): culled outright
        static constexpr float RadiusCull = 1e-2f;

        // Below this signed area (mm^2) a closed loop bounds no cuttable region: it is stitching noise (cusp slivers off a merged lobe). Judged nowhere NEAR real geometry (smallest legitimate frontier pockets are step^2 scale, ~1 mm^2) and vital everywhere: an un-culled zero-area CW sliver offset-DILATES into an immortal spurious ring (see offsetValid).
        static constexpr float SignificanceAreaFloor = 0.05f;

        Chain() = default;
        Chain(Chain&&) = default;
        Chain& operator=(Chain&&) = default;

        Chain clone() const {
            Chain c;
            c.closed = closed; c.id = id; c.parent = parent; c.role = role;
            for (const auto& e : edges) { c.edges.push_back(e->clone()); }
            return c;
        }

        // The opposite journey: edge order flipped, every edge reversed. For CONDITIONING inputs whose winding is arbitrary -- never for reorienting the method's own outputs (Method/05).
        Chain reversed() const {
            Chain r;
            r.closed = closed; r.id = id; r.parent = parent; r.role = role;
            for (size_t i = edges.size(); i-- > 0; ) { r.edges.push_back(Kern::reversedEdge(*edges[i])); }
            return r;
        }

        void tessellate(std::vector<Pos>& out) const {
            bool first = true;
            for (const auto& e : edges) { e->tessellate(out, !first); first = false; }
        }
        void tessellateDirected(std::vector<Pos>& out) const {
            bool first = true;
            for (const auto& e : edges) { e->tessellateDirected(out, !first); first = false; }
        }

        // Measures

        // Exact signed area (Green's theorem, closed-form per edge). This is the judge for ORIENTING inputs: exact, immune to the tangency-cusp turning flips that bit the old kernel three times.
        float signedArea() const {
            float twice = 0.0f;
            for (const auto& e : edges) { twice += e->signedAreaTerm(); }
            return 0.5f * twice;
        }

        // Handedness as a REGION fact: +1 CCW (bounds material inside), -1 CW (an island), 0 degenerate. Judged by area, not turning -- the sense sigma is about what the chain's region does under the offset, an area question.
        int handedness() const {
            float a = signedArea();
            return (a > 1e-4f) ? 1 : (a < -1e-4f ? -1 : 0);
        }

        // The accumulated tangent rotation: 2*pi * rotation index for a closed path. This is the judge for OFFSET OUTPUT children (the inversion guard): burr lobes can outweigh a thin loop's area and flip its sign, but they can never flip total turning, because offsetting preserves the tangent structure of its source.
        float totalTurning() const {
            const size_t n = edges.size();
            if (n == 0) { return 0.0f; }
            double total = 0.0;
            for (const auto& e : edges) {
                Pos c; float r, a0, sweep; int chir;
                if (Kern::circularOf(*e, c, r, a0, sweep, chir)) {
                    total += static_cast<double>(chir) * sweep;
                }
                else if (Kern::isOffEll(*e)) {
                    const OffsetEllipse2* o = Kern::asOff(*e);
                    const int N = 64;
                    dPos prev(o->tangentAt(0.0f));
                    for (int i = 1; i <= N; i++) {
                        dPos t(o->tangentAt(static_cast<float>(i) / N));
                        total += std::atan2(prev.cross(t), prev.dot(t));
                        prev = t;
                    }
                }
                else if (Kern::isSpline(*e)) {
                    const Spline2* o = Kern::asSpline(*e);
                    const int N = 64;
                    dPos prev(o->tangentAt(0.0f));
                    for (int i = 1; i <= N; i++) {
                        dPos t(o->tangentAt(static_cast<float>(i) / N));
                        total += std::atan2(prev.cross(t), prev.dot(t));
                        prev = t;
                    }
                }
            }
            const size_t corners = closed ? n : (n - 1);
            for (size_t k = 0; k < corners; k++) {
                dPos a(Kern::eEndDir(*edges[k]));
                dPos b(Kern::eStartDir(*edges[(k + 1) % n]));
                total += std::atan2(a.cross(b), a.dot(b));   // signed exterior angle
            }
            return static_cast<float>(total);
        }
        int turningSign() const {
            float t = totalTurning();
            return (t > 1e-3f) ? 1 : (t < -1e-3f) ? -1 : 0;
        }

        // Winding number around a probe, by ACCUMULATED SUBTENDED ANGLE -- deliberately not ray casting: the angle sum is continuous in every vertex position and can only destabilise when the probe sits ON the curve (Method/04 scar tissue).
        int windingAround(Pos probe) const {
            const dPos P(probe);
            double total = 0.0;

            auto chord = [&](dPos a, dPos b) {
                dPos u = a - P, v = b - P;
                return std::atan2(u.cross(v), u.dot(v));   // wrapped to (-pi, pi]
            };

            for (const auto& e : edges) {
                Pos c; float r, a0, sweep; int chir;
                if (Kern::circularOf(*e, c, r, a0, sweep, chir)) {
                    dPos C(c);
                    double R = r;
                    dPos p0 = C + dPos::fromAngle(a0) * R;
                    dPos p1 = C + dPos::fromAngle(static_cast<double>(a0) + sweep) * R;
                    double wrapped = chord(p0, p1);          // CCW-span bearing delta, mod 2pi
                    double d;
                    if ((P - C).len() < R) {
                        // Inside the circle the bearing advances monotonically CCW, so the true advance lies in (0, 2pi): unique branch.
                        d = (wrapped > 0.0) ? wrapped : wrapped + dTAU;
                    }
                    else {
                        d = wrapped;                         // outside: |advance| < pi
                    }
                    total += (chir > 0) ? d : -d;            // travel direction signs it
                }
                else if (Kern::isOffEll(*e)) {
                    const OffsetEllipse2* o = Kern::asOff(*e);
                    const int N = 256;
                    dPos prev(o->pointAt(0.0f));
                    for (int i = 1; i <= N; i++) {
                        dPos cur(o->pointAt(static_cast<float>(i) / N));
                        total += chord(prev, cur);
                        prev = cur;
                    }
                }
                else if (Kern::isSpline(*e)) {
                    const Spline2* o = Kern::asSpline(*e);
                    const int N = 256;
                    dPos prev(o->pointAt(0.0f));
                    for (int i = 1; i <= N; i++) {
                        dPos cur(o->pointAt(static_cast<float>(i) / N));
                        total += chord(prev, cur);
                        prev = cur;
                    }
                }
                else {
                    total += chord(dPos(Kern::eStart(*e)), dPos(Kern::eEnd(*e)));
                }
            }
            return static_cast<int>(std::lround(total / dTAU));
        }

        // The point of greatest x anywhere on the loop: provably unenclosed within its own cluster, and the containment probe (on this loop, far from another unless they genuinely touch).
        static Pos maxXPoint(const Chain& L) {
            Pos best(-1e30f, 0.0f);
            for (const auto& e : L.edges) {
                Pos p = Kern::edgeMaxXPoint(*e);
                if (p.x > best.x) { best = p; }
            }
            return best;
        }

        // A point strictly INSIDE the loop, never on it: nudge inward (-x) of the max-x point, verified against the loop's own winding. A boundary point's winding against a touching sibling is a coin flip; an interior point is unambiguous under tangency.
        static Pos interiorPoint(const Chain& L) {
            const Pos mx = maxXPoint(L);
            for (float nudge : { 1e-3f, 1e-2f, 1e-1f }) {
                const Pos q(mx.x - nudge, mx.y);
                if (L.windingAround(q) != 0) { return q; }
            }
            return mx;   // degenerate sliver: the boundary point is all there is
        }

        // Hard-weld: enforce EXACT junction coincidence by snapping both neighbours to their mean -- never a micro bridge segment (a sub-eps edge is invisible to crossing detection and poison to the walk, and a marked-closed loop with a real gap breaks coincidence matching). assemble() connects within eps, so every interior joint carries up-to-eps slop until this runs. An edge that cannot absorb the snap (degenerate under mutation) falls back to one exact bridge.
        static void hardWeld(Chain& c, float eps) {
            const size_t n = c.edges.size();
            if (n == 0) { return; }
            const size_t joints = c.closed ? n : (n - 1);
            std::vector<size_t> bridgeAfter;
            for (size_t k = 0; k < joints; k++) {
                Stoicheion& A = *c.edges[k];
                Stoicheion& B = *c.edges[(k + 1) % n];
                const Pos a = Kern::eEnd(A);
                const Pos b = Kern::eStart(B);
                const float gap = (b - a).pythag();
                if (gap <= 1e-9f) { continue; }      // already exact
                if (gap > eps)    { continue; }      // a real gap: not ours to hide
                const Pos mid = (a + b) * 0.5f;
                const bool aOk = Kern::setEdgeEnd(A, mid);
                const bool bOk = Kern::setEdgeStart(B, mid);
                if (!aOk || !bOk) { bridgeAfter.push_back(k); }
            }
            for (size_t bi = bridgeAfter.size(); bi > 0; bi--) {
                const size_t k = bridgeAfter[bi - 1];
                const Pos a = Kern::eEnd(*c.edges[k]);
                const Pos b = Kern::eStart(*c.edges[(k + 1) % n]);
                if ((b - a).pythag() > 1e-9f) {
                    c.edges.insert(c.edges.begin() + static_cast<std::ptrdiff_t>(k) + 1,
                                   std::make_unique<Segment2>(a, b));
                }
            }
        }

        // The nesting primitive, by the algebra: loops that INTERSECT are neither neighbors nor nested (the crossing machinery owns that case); loops that are crossing-free are nested-or-neighbors, and that is a GLOBAL fact decidable at ONE sample -- so the effort goes not into many samples (a vote averages around the problem) but into choosing the SAFEST one: the point of the inner loop with MAXIMUM clearance from the outer's curve -- the deepest point of the SDF, where the side reading cannot be ambiguous.
        static bool containedIn(const Chain& inner, const Chain& outer) {
            Pos best{}; float bestD = -1.0f;
            for (const auto& e : inner.edges) {
                const Pos m = Kern::edgeMidpoint(*e);
                float d = 1e30f;
                for (const auto& o : outer.edges) { d = std::min(d, o->distanceTo(m)); }
                if (d > bestD) { bestD = d; best = m; }
            }
            if (bestD < 0.0f) { return false; }
            return outer.windingAround(best) != 0;
        }

        // The parity cull (the nesting rule). Within ONE RESOLVED ring system, alternation is law: the frontier of { winding >= 1 } can never contain a ring nested inside a SAME-handed ring -- the inner ring would bound material on both sides or neither. Any such survivor is a doubled bound or a fabricated artifact loop (the non-tangent-joint micro loop), and every downstream winding it touches is corrupted. CULL the inner ring; never reorient -- flipping fabricates a hole that does not exist. Orientation is established once at consumption and blindly applied; this is the enforcement pass, not a repair.
        static void cullSameHandedNested(std::vector<Chain>& rings) {
            const size_t n = rings.size();
            if (n < 2) { return; }
            std::vector<int>  turn(n, 0);
            std::vector<Pos>  probe(n);
            for (size_t i = 0; i < n; i++) {
                if (!rings[i].closed || rings[i].edges.empty()) { continue; }
                int h = rings[i].handedness();
                if (h == 0) { h = rings[i].turningSign(); }   // degenerate area: turning judges
                turn[i]  = h;
                probe[i] = interiorPoint(rings[i]);
            }
            std::vector<char> cull(n, 0);
            for (size_t i = 0; i < n; i++) {
                if (turn[i] == 0) { continue; }
                for (size_t j = 0; j < n; j++) {
                    if (j == i || turn[j] != turn[i]) { continue; }
                    // i strictly inside j, same handedness -> i is the violator (containment by vote -- no single-probe coin flips)
                    if (containedIn(rings[i], rings[j])
                        && !containedIn(rings[j], rings[i])) {
                        cull[i] = 1;
                        break;
                    }
                }
            }
            size_t w = 0;
            for (size_t i = 0; i < n; i++) {
                if (!cull[i]) { rings[w++] = std::move(rings[i]); }
            }
            rings.resize(w);
        }

        // Beat one: the blind offset (Method/02) Every edge slides to the LEFT of its own travel by `amount`; every corner is bridged by the join known by construction. No winding test, no preferred side -- inset vs outset EMERGES from each chain's own chirality. The result may self-intersect; those crossings are beat two's input.

        // EXPERIMENT SWITCH: true = the legacy always-arc corner joins; false = arcs only where the push side OPENS A GAP (convex relative to the offset), straight lines on the overlap side. The overlap-side join is transient scar tissue the crossing walk excises anyway, and a line cannot wrap, cannot flip its side, and crosses cleanly exactly where the arc construction's side pick lives on the epsilon knife edge.
        static constexpr bool JoinAlwaysArcs = false;

        // The corner bridge between two offset edges. `arcs` selects the mode; tA/tB are the ORIGINAL edges' tangents at the original vertex (the convexity judges), offDirA the offset edge's end tangent (the arc's side pick). Null = tangent-continuous / degenerate: no bridge needed.
        static std::unique_ptr<Stoicheion> join(bool arcs, Pos corner, Pos endA, Pos startB,
                                                Pos tA, Pos tB, Pos offDirA, float amount) {
            if (!arcs) {
                // Gap opens only when the push side is the OUTSIDE of the turn: cross(tA, tB) * amount < 0 under the left-push convention. Overlap-side AND near-tangent (cross ~ 0) corners take the line.
                if (static_cast<double>(tA.cross(tB)) * static_cast<double>(amount) >= 0.0) {
                    return std::make_unique<Segment2>(endA, startB);
                }
            }

            // All angle math in DOUBLE: a nearly-complete around-the-bend arc subtracts two nearly-equal angles.
            dPos C(corner);
            dPos radial = dPos(endA) - C;
            double r = radial.len();
            if (r < 1e-12) { return nullptr; }
            int s = (dPos(offDirA).dot(radial.perp()) >= 0.0) ? 1 : -1;
            double angA = radial.angle(), angB = (dPos(startB) - C).angle();
            double span = (s > 0) ? dWrap(angB - angA) : dWrap(angA - angB);
            // The side dot above is EPSILON-SCALE at a near-tangent corner: a coin-flipped side wraps the span to 2pi - true and the bridge becomes a near-full-circle arc (the opposite-handed mini loop; the same sign-flip-wraps-full-turn disease as the firmware G2/G3 bug). A corner's turn can never exceed pi, so span > pi IS the proof the side flipped: take the complement. Threshold-free, exact.
            if (span > 3.14159265358979324 + 1e-3) {
                s = -s;
                span = 6.28318530717958648 - span;
            }
            double midAng = angA + static_cast<double>(s) * span * 0.5;
            return std::make_unique<Arc2>(corner, endA, startB,
                                          (C + dPos::fromAngle(midAng) * r).f());
        }

        Chain offsetRaw(float amount) const {
            Chain result;
            if (edges.empty()) { return result; }
            result.closed = closed;

            // A lone prime closed curve is its own loop: no corners, just the blind concentric left-offset.
            if (edges.size() == 1 && Kern::isCircle(*edges.front())) {
                std::unique_ptr<Stoicheion> o = Kern::offsetLeft(*edges.front(), amount);
                if (Kern::pruneSmallFeatures && Kern::edgeLength(*o) < Kern::PruneLength) { return result; }
                if (Kern::asCirc(*o)->radius() < RadiusCull) { return result; }   // driven to (near) zero: extinct
                result.edges.push_back(std::move(o));
                return result;
            }
            if (edges.size() == 1 && Kern::isOffEll(*edges.front()) && Kern::asOff(*edges.front())->fullLoop()) {
                std::unique_ptr<Stoicheion> o = Kern::offsetLeft(*edges.front(), amount);
                if (Kern::pruneSmallFeatures && Kern::edgeLength(*o) < Kern::PruneLength) { return result; }
                result.edges.push_back(std::move(o));
                return result;
            }
            if (edges.size() == 1 && Kern::isSpline(*edges.front()) && Kern::asSpline(*edges.front())->fullLoop()) {
                std::unique_ptr<Stoicheion> o = Kern::offsetLeft(*edges.front(), amount);
                if (Kern::pruneSmallFeatures && Kern::edgeLength(*o) < Kern::PruneLength) { return result; }
                result.edges.push_back(std::move(o));
                return result;
            }

            // Guard: a chain holding anything we don't offset passes through.
            for (const auto& e : edges) {
                SKind t = e->type();
                if (t != SKind::Segment && t != SKind::Arc && t != SKind::OffsetEllipse && t != SKind::Spline) {
                    for (const auto& src : edges) { result.edges.push_back(src->clone()); }
                    return result;
                }
            }

            // The push: offset each edge alone; an edge that collapses under the offset generates nothing (its neighbours join across the hole).
            std::vector<std::unique_ptr<Stoicheion>> off;
            std::vector<size_t> src;                   // off[i] came from edges[src[i]]
            off.reserve(edges.size());
            src.reserve(edges.size());
            for (size_t i = 0; i < edges.size(); i++) {
                std::unique_ptr<Stoicheion> o = Kern::offsetLeft(*edges[i], amount);
                if (Kern::pruneSmallFeatures && Kern::edgeLength(*o) < Kern::PruneLength) { continue; }

                // An arc driven to (near) zero radius is a collapsed feature, not a boundary. The length prune alone misses it: through-zero inversion keeps the through-point's bearing, so the junk re-emerges with the COMPLEMENT sweep -- long enough to survive, chirality flipped. At fillet-radius == tool-radius (with section slop) this is the EXPECTED case, not the corner: cull by RADIUS, outright, and let the corner construction join the neighbours directly.
                if (Kern::isArc(*o) && Kern::asArc(*o)->radius() < RadiusCull) { continue; }

                off.push_back(std::move(o));
                src.push_back(i);
            }

            const size_t n = off.size();
            if (n == 0) { return result; }
            const size_t corners = closed ? n : (n > 0 ? n - 1 : 0);

            // The uniform corner construction. Near-tangent seams (the STEP signature) are FUSED, not bridged -- a sliver join loop there is numerical poison (Method/02 scar tissue).
            std::vector<std::unique_ptr<Stoicheion>> joinArc(corners);
            for (size_t k = 0; k < corners; k++) {
                Pos corner = Kern::eEnd(*edges[src[k]]);              // the original corner
                Pos endA = Kern::eEnd(*off[k]);
                Pos startB = Kern::eStart(*off[(k + 1) % n]);
                if ((endA - startB).pythag() < 1e-4f) { continue; }   // tangent-continuous: already joined

                // A hole in the source (an edge pruned or radius-culled between these two) voids the corner-arc construction -- it presumes true adjacency, and across a hole the arc's far end misses startB by up to the culled radius, shredding the chain at kernel eps. Bridge the hole EXACTLY: a straight segment, endpoints verbatim.
                if (((src[k] + 1) % edges.size()) != src[(k + 1) % n]) {
                    joinArc[k] = std::make_unique<Segment2>(endA, startB);
                    continue;
                }

                if (Kern::fuseCloseTangencies) {
                    Pos ta = Kern::eEndDir(*off[k]);
                    Pos tb = Kern::eStartDir(*off[(k + 1) % n]);
                    if (std::fabs(ta.cross(tb)) < Kern::FuseSin && ta.dot(tb) > 0.0f) {
                        Pos mid = (endA + startB) * 0.5f;
                        bool aOk = Kern::setEdgeEnd(*off[k], mid);
                        bool bOk = Kern::setEdgeStart(*off[(k + 1) % n], mid);
                        if (!aOk || !bOk) {
                            joinArc[k] = std::make_unique<Segment2>(
                                Kern::eEnd(*off[k]), Kern::eStart(*off[(k + 1) % n]));
                        }
                        continue;
                    }
                }

                joinArc[k] = join(JoinAlwaysArcs, corner, endA, startB,
                                  Kern::eEndDir(*edges[src[k]]),
                                  Kern::eStartDir(*edges[src[(k + 1) % n]]),
                                  Kern::eEndDir(*off[k]), amount);
            }

            for (size_t i = 0; i < n; i++) {
                result.edges.push_back(std::move(off[i]));
                if (i < corners && joinArc[i]) { result.edges.push_back(std::move(joinArc[i])); }
            }
            return result;
        }

        // Beat two: the crossing-number field (Method/03) Split the raw offset at every self-crossing and walk it in travel order, accumulating a level: at each crossing, the OTHER strand passing right-to-left counts +1, left-to-right -1. Only DIFFERENCES are meaningful (the walk's start is arbitrary), so levels are normalised max -> 0. Which extreme is the true frontier is beat three's decision, NOT the field's.
        struct Field {
            std::vector<std::unique_ptr<Stoicheion>> edges;   // travel order, direction preserved
            std::vector<int> level;                           // accumulated level per sub-edge
        };

        Field field(float eps = 1e-3f) const {
            const size_t n = edges.size();

            // 1. every self-crossing, remembering which two edges meet there. (edgeCross is transversal-only: a tangential kiss is rejected at detection -- its +-1 would be float noise in the walk.)
            struct Rec { Pos p; size_t i, j; bool firedI = false, firedJ = false; };
            std::vector<Rec> recs;
            std::vector<std::vector<Pos>> pts(n);
            for (size_t i = 0; i < n; i++) {
                for (size_t j = i + 1; j < n; j++) {
                    std::vector<Pos> c; Kern::edgeCross(*edges[i], *edges[j], c);
                    for (const Pos& p : c) { recs.push_back({ p, i, j }); pts[i].push_back(p); pts[j].push_back(p); }
                }
            }

            // 2. split every edge at its cuts, keeping chain order; each sub-edge remembers its original edge.
            Field run;
            std::vector<size_t> from;
            for (size_t i = 0; i < n; i++) {
                Kern::splitEdgeAtPoints(*edges[i], pts[i], run.edges);
                while (from.size() < run.edges.size()) { from.push_back(i); }
            }

            // 3. walk strictly in order. A sub-edge wears the level in effect as it is entered; the level updates at its end vertex if that vertex is a crossing (the OTHER strand's travel direction there).
            run.level.resize(run.edges.size(), 0);
            int count = 0;
            for (size_t k = 0; k < run.edges.size(); k++) {
                run.level[k] = count;
                Pos v = Kern::eEnd(*run.edges[k]);
                for (Rec& r : recs) {
                    if ((r.p - v).pythag() >= eps) { continue; }
                    size_t other = (from[k] == r.i) ? r.j
                                 : (from[k] == r.j) ? r.i : static_cast<size_t>(-1);
                    if (other == static_cast<size_t>(-1)) { continue; }
                    // ONE increment per crossing per strand: two crossings within eps of each other put a sub-edge end inside BOTH records' capture radius, and double-firing corrupts the level of every sub-edge after it (the close-intersection case).
                    bool& fired = (from[k] == r.i) ? r.firedI : r.firedJ;
                    if (fired) { continue; }
                    fired = true;
                    // Both tangents from the ORIGINAL curves at the true crossing point: a sub-edge between two close cuts is micro, and its eEndDir is noise -- exactly where the +-1 must not coin-flip.
                    Pos tOur = Kern::travelDirAt(*edges[from[k]], r.p);
                    Pos tOther = Kern::travelDirAt(*edges[other], r.p);
                    count += (tOur.cross(tOther) > 0.0f) ? 1 : -1;
                }
            }

            // Normalise relative to the max: start-invariant with no external measurement. Levels read 0 (max) downward.
            if (!run.level.empty()) {
                int maxN = run.level.front();
                for (int v : run.level) { maxN = std::max(maxN, v); }
                for (int& v : run.level) { v -= maxN; }
            }
            return run;
        }

        // A field fragment: a well-formed sub-chain plus the level its pieces share. (Holds a Chain by value: defined out of line, below.)
        struct Fragment;
        std::vector<Fragment> fragments(float eps = 1e-3f) const;

        // Beat three: sigma-keyed extraction (Method/03)
        std::vector<Chain> offsetValid(float amount, float eps = 1e-3f) const;

        // Merging (Method/04, 06)
        static std::vector<Chain> merge(std::vector<Chain> allChains, float eps = 1e-3f);

        // Building

        // Grow one chain out of the pool, connecting end-to-start until it closes or runs out. `allowReverse` permits flipping an edge to make a connection -- legitimate when assembling unordered input strokes, FORBIDDEN when reassembling oriented frontier pieces (a silent reversal there fabricates wrong-handed loops; the piece set must fail visibly instead).
        static Chain buildOne(std::vector<std::unique_ptr<Stoicheion>>& pool, float eps,
                              bool allowReverse) {
            Chain c;
            if (pool.empty()) { return c; }
            c.edges.push_back(std::move(pool.front()));
            pool.erase(pool.begin());
            Pos loopStart = Kern::eStart(*c.edges.front());

            while (!pool.empty()) {
                Pos cur = Kern::eEnd(*c.edges.back());
                int best = -1; bool rev = false; float bestD = eps;
                for (size_t i = 0; i < pool.size(); i++) {
                    float dS = (cur - Kern::eStart(*pool[i])).pythag();
                    if (dS < bestD) { best = static_cast<int>(i); rev = false; bestD = dS; }
                    if (allowReverse) {
                        float dE = (cur - Kern::eEnd(*pool[i])).pythag();
                        if (dE < bestD) { best = static_cast<int>(i); rev = true; bestD = dE; }
                    }
                }
                if (best < 0) { break; }
                std::unique_ptr<Stoicheion> next = rev ? Kern::reversedEdge(*pool[best]) : std::move(pool[best]);
                pool.erase(pool.begin() + best);
                c.edges.push_back(std::move(next));
                if ((Kern::eEnd(*c.edges.back()) - loopStart).pythag() <= eps) { c.closed = true; break; }
            }
            hardWeld(c, eps);
            return c;
        }

        // Assemble chains from raw stoicheia (reversal allowed: the input is unordered and unoriented). Prime closed curves are their own loops; an ellipse enters EXACTLY as zero-distance offset-ellipse spans cut at quadrant boundaries (curvature extrema, where offset cusps form).
        static std::vector<Chain> assemble(const std::vector<std::unique_ptr<Stoicheion>>& src,
                                           float eps = 1e-3f) {
            std::vector<Chain> out;
            std::vector<std::unique_ptr<Stoicheion>> pool;

            auto ellipseToEdges = [&pool](Pos c, Pos u, Pos v, float t0, float span) {
                int spans = std::max(2, static_cast<int>(std::ceil(std::fabs(span) / (TAU * 0.25f))));
                for (int i = 0; i < spans; i++) {
                    float ta = t0 + span * (static_cast<float>(i) / spans);
                    float tb = t0 + span * (static_cast<float>(i + 1) / spans);
                    pool.push_back(std::make_unique<OffsetEllipse2>(c, u, v, 0.0f, ta, tb));
                }
            };

            for (const auto& e : src) {
                if (!e || e->locked || e->construction || e->isPoint()) { continue; }
                switch (e->type()) {
                    case SKind::Circle: {
                        Chain c; c.closed = true; c.edges.push_back(e->clone()); out.push_back(std::move(c));
                        break;
                    }
                    case SKind::Segment:
                    case SKind::Arc: { pool.push_back(e->clone()); break; }
                    case SKind::Spline: {
                        if (static_cast<const Spline2*>(e.get())->fullLoop()) {
                            Chain c; c.closed = true; c.edges.push_back(e->clone()); out.push_back(std::move(c));
                        }
                        else { pool.push_back(e->clone()); }
                        break;
                    }
                    case SKind::Ellipse: {
                        const Ellipse2* el = static_cast<const Ellipse2*>(e.get());
                        ellipseToEdges(el->c, el->u, el->v, 0.0f, TAU);
                        break;
                    }
                    case SKind::EllipseArc: {
                        const EllipseArc2* ea = static_cast<const EllipseArc2*>(e.get());
                        ellipseToEdges(ea->c, ea->u, ea->v, ea->a0, ea->span());
                        break;
                    }
                    default: break;
                }
            }

            while (!pool.empty()) {
                Chain c = buildOne(pool, eps, /*allowReverse=*/true);
                if (c.edges.empty()) { break; }
                out.push_back(std::move(c));
            }

            for (Chain& c : out) { c.id = newId(); c.parent = 0; }
            return out;
        }
    };

    // A field fragment (out of line: holds a Chain by value).
    struct Chain::Fragment { Chain chain; int level = 0; };

    inline std::vector<Chain::Fragment> Chain::fragments(float eps) const {
        Field run = field(eps);
        std::vector<Fragment> out;

        for (size_t k = 0; k < run.edges.size(); k++) {
            bool fresh = out.empty()
                || out.back().level != run.level[k]
                || (Kern::eEnd(*out.back().chain.edges.back()) - Kern::eStart(*run.edges[k])).pythag() > eps;
            if (fresh) { out.push_back(Fragment{ Chain{}, run.level[k] }); }
            out.back().chain.edges.push_back(std::move(run.edges[k]));
        }

        // Wrap merge: on a closed walk the first and last fragments can be two halves of the same stretch.
        if (closed && out.size() > 1) {
            Fragment& head = out.front();
            Fragment& tail = out.back();
            if (head.level == tail.level &&
                (Kern::eEnd(*tail.chain.edges.back()) - Kern::eStart(*head.chain.edges.front())).pythag() <= eps) {
                for (auto& e : head.chain.edges) { tail.chain.edges.push_back(std::move(e)); }
                head.chain.edges = std::move(tail.chain.edges);
                out.pop_back();
            }
        }

        for (Fragment& nc : out) {
            if (nc.chain.edges.empty()) { continue; }
            nc.chain.closed =
                (Kern::eEnd(*nc.chain.edges.back()) - Kern::eStart(*nc.chain.edges.front())).pythag() <= eps;
        }
        return out;
    }

    // The method end to end for ONE chain: blind offset, crossing field, side-keyed extreme, oriented stitch, sigma-gated inversion guard. TWO different keys, and they are not the same variable: * The EXTREME is keyed on sign(amount) -- the SIDE the curve was pushed, matching the field's own travel-left orientation. Leftward (amount > 0): the frontier is the MINIMUM level. Rightward (amount < 0): the MAXIMUM. It cannot key on sigma: mirror the plane and a proven min-case becomes its opposite-handed twin with every crossing increment negated -- same sigma, opposite extreme label. Anchors: the legacy engine (all chains at amount > 0, min, proven in production, CW-growing islands included) and the crescent seed (CCW growing at amount < 0, max, verified on a real part). * The GUARD is keyed on sigma = sign(amount * h) -- a REGION question: only a shrinking chain (sigma > 0) can collapse through zero, so only there are handedness-flipped children the fold rather than a boundary. Growth (sigma < 0) cannot invert; the guard is off. [Growing secondary births -- a void pinched off between growing arms -- are still dropped with the non-extreme levels: a known, safe-direction limitation.] h is judged by exact signed AREA (the sense is a region question); the guard compares TURNING (offsetting preserves tangent structure, so turning is the burr-immune judge of an offset child's handedness).
    inline std::vector<Chain> Chain::offsetValid(float amount, float eps) const {
        std::vector<Chain> out;

        // Closed chains only: offsetting is a statement about a bounded region's boundary; an open chain bounds nothing.
        if (!closed || edges.empty()) { return out; }

        // Geometric-significance floor. A near-zero-area loop is boundary NOISE (a cusp sliver shed by a merge), and offsetting one is how noise turns IMMORTAL: area ~0 defeats handedness, turningSign takes over, and a CW-turning sliver reads as a growth-side chain -- it DILATES every generation into a spurious ring that outlives the real geometry and cascades to the generation cap. Below the floor there is no region to bound: refuse to offset it at all.
        if (std::fabs(signedArea()) < SignificanceAreaFloor) { return out; }

        int h = handedness();
        if (h == 0) { h = turningSign(); }        // degenerate area: fall back to turning
        if (h == 0) { return out; }               // no orientation at all: nothing to say
        const float sigma = amount * static_cast<float>(h);

        Chain raw = offsetRaw(amount);
        if (raw.edges.empty()) { return out; }

        std::vector<Fragment> frags = raw.fragments(eps);
        if (frags.empty()) { return out; }

        int lo = frags.front().level, hi = frags.front().level;
        for (const Fragment& f : frags) { lo = std::min(lo, f.level); hi = std::max(hi, f.level); }
        const int target = (amount > 0.0f) ? lo : hi;   // the SIDE, not sigma (see above)

        std::vector<std::unique_ptr<Stoicheion>> pool;
        for (Fragment& f : frags) {
            if (f.level != target) { continue; }
            for (auto& e : f.chain.edges) { pool.push_back(std::move(e)); }
        }
        while (!pool.empty()) {
            Chain m = buildOne(pool, eps, /*allowReverse=*/false);   // oriented pieces: never reverse
            if (m.edges.empty()) { break; }
            m.id = newId();
            m.parent = id;
            m.role = role;
            out.push_back(std::move(m));
        }

        // The inversion guard, UNCONDITIONAL. The method is create-then-cull: opposite-handed loops are expected intermediates (through-zero re-emergence, pinched voids, close-crossing level slips) and the cull is the contract. The old sigma gate guarded only the SHRINKING side -- "growth cannot invert" -- but a GROWING opposite-handed chain (an annulus seed's hole family grows every generation) ran with the guard off, and any inverted survivor that reached the kept extreme level sailed through: the un-culled opposite-handed loop class. At the extreme level an opposite-handed child is never a legitimate frontier (legitimate pinched voids live on non-extreme levels by the walk's own accounting), so the cull needs no gate.
        {
            const int srcH = h;
            std::erase_if(out, [srcH](const Chain& m) {
                return m.closed && m.turningSign() == -srcH;
            });
        }
        (void)sigma;

        // The parity cull: one chain's valid offset children may split into disjoint siblings, never into same-handed nesting -- a same-handed child inside another is the non-tangent-joint artifact loop (a level the crossing walk failed to separate). Culled at birth, before it can poison a merge or a cascade winding.
        cullSameHandedNested(out);
        return out;
    }

    // The merge: resolved loops meet; the material region is selected by an EXPLICIT winding rule and its frontier extracted (Method/04, 06). Doctrine orientation is asserted upstream (material on the left of every edge), so material is simply { winding >= 1 } -- no auto-detected system sign, which is exactly what an enclosing subtrahend perimeter used to hijack. Roles close the remaining hole: a cluster containing no MINUEND loop can bound no material and yields nothing, so disjoint avoid geometry can never survive as a spurious frontier. Envelopment without intersection needs NO special case: clustering is by crossing OR CONTAINMENT, so a grown avoid loop enveloped by the grown fill joins its cluster, and the winding keeps the avoid's loop as the frontier while discarding the fill's (the bore: the seed IS the avoid's standoff ring -- "we moved away from the avoid region").
    inline std::vector<Chain> Chain::merge(std::vector<Chain> allChains, float eps) {

        std::vector<Chain> out;
        std::vector<Chain> loops;
        for (Chain& ch : allChains) {
            if (ch.closed && !ch.edges.empty()) { loops.push_back(std::move(ch)); }
            else if (!ch.edges.empty())         { out.push_back(std::move(ch)); }   // opens pass through
        }

        const size_t n = loops.size();
        if (n == 0) { return out; }

        // Cluster by INTERACTION (cross, or one contains the other). Disjoint profiles are independent worlds.
        auto interacts = [&](size_t a, size_t b) {
            std::vector<Pos> x;
            for (const auto& ea : loops[a].edges) {
                for (const auto& eb : loops[b].edges) {
                    Kern::edgeCross(*ea, *eb, x);
                    if (!x.empty()) { return true; }
                }
            }
            // containment by VOTE (see containedIn): no single-probe coin flips
            if (containedIn(loops[a], loops[b])) { return true; }
            if (containedIn(loops[b], loops[a])) { return true; }
            return false;
        };
        std::vector<int> cluster(n, -1);
        int clusters = 0;
        for (size_t a = 0; a < n; a++) {
            if (cluster[a] >= 0) { continue; }
            cluster[a] = clusters;
            std::vector<size_t> queue{ a };
            while (!queue.empty()) {
                size_t cur = queue.back(); queue.pop_back();
                for (size_t b = 0; b < n; b++) {
                    if (cluster[b] >= 0 || !interacts(cur, b)) { continue; }
                    cluster[b] = clusters;
                    queue.push_back(b);
                }
            }
            clusters++;
        }

        for (int c = 0; c < clusters; c++) {
            std::vector<size_t> mem;
            for (size_t a = 0; a < n; a++) { if (cluster[a] == c) { mem.push_back(a); } }

            // The minuend rule: subtrahend geometry alone bounds nothing.
            bool hasMinuend = false;
            for (size_t a : mem) { if (loops[a].role == Role::Minuend) { hasMinuend = true; break; } }
            if (!hasMinuend) { continue; }

            // A lone loop has nothing to overlap: through untouched.
            if (mem.size() == 1) { out.push_back(std::move(loops[mem[0]])); continue; }

            // LEVEL-WALK CLASSIFICATION -- the field algebra, cross-chain. Side probes classified each piece independently at noise-scale offsets and could not distinguish, among several same-handed result chains, the true frontier from non-extremal junk (the annulus that never merges with its inset ring). Crossing numbers can: the left-side winding of a piece changes ONLY at crossings, by a sign read off two exact tangents. 1. record every crossing, BOTH strands (transversal + vertex-landing) 2. split each loop's edges at its crossings, travel order kept 3. per loop: read the other loops' winding ONCE, at the safest piece (max clearance), then ACCUMULATE along the walk 4. closed-walk invariant: one full lap returns to its starting level; violation = a swallowed crossing -> legacy side probes for this cluster (fail toward the proven path) 5. keep pieces whose LEFT level is exactly 1 (material-left)

            // -- 1. crossing records
            struct XRec {
                Pos p; size_t la, lb;
                const Stoicheion* ea; const Stoicheion* eb;
                bool usedA = false, usedB = false;
            };
            std::vector<XRec> recs;
            for (size_t ia = 0; ia < mem.size(); ia++) {
                for (size_t ib = ia + 1; ib < mem.size(); ib++) {
                    const size_t a = mem[ia], b = mem[ib];
                    for (const auto& eA : loops[a].edges) {
                        for (const auto& eB : loops[b].edges) {
                            std::vector<Pos> x;
                            Kern::edgeCross(*eA, *eB, x);
                            for (const Pos& p : x) { recs.push_back({ p, a, b, eA.get(), eB.get() }); }
                            // vertex-landing: a genuine crossing ON a vertex is invisible to edgeCross (strict interiority)
                            const Pos vB = Kern::eStart(*eB);
                            if (eA->distanceTo(vB) < 2e-3f) { recs.push_back({ vB, a, b, eA.get(), eB.get() }); }
                            const Pos vA = Kern::eStart(*eA);
                            if (eB->distanceTo(vA) < 2e-3f) { recs.push_back({ vA, a, b, eA.get(), eB.get() }); }
                        }
                    }
                }
            }
            // duplicates (an edgeCross hit re-found by the vertex sweep) would double a transition; genuine distinct crossings sit further apart
            std::erase_if(recs, [&recs](const XRec& r) {
                for (const XRec& s : recs) {
                    if (&s == &r) { return false; }              // keep the first of a pair
                    if (s.la == r.la && s.lb == r.lb
                        && (s.p - r.p).pythag() < 5e-4f) { return true; }
                }
                return false;
            });

            // -- 2. per-loop pieces, travel-ordered
            struct TaggedPiece { std::unique_ptr<Stoicheion> e; size_t loop; };
            std::vector<TaggedPiece> pieces;
            std::vector<std::pair<size_t, size_t>> range(n, { 0, 0 });   // per-loop [first, count)
            for (size_t a : mem) {
                const size_t first = pieces.size();
                for (const auto& edge : loops[a].edges) {
                    std::vector<Pos> cuts;
                    for (const XRec& r : recs) {
                        if (r.la == a && r.ea == edge.get()) { cuts.push_back(r.p); }
                        if (r.lb == a && r.eb == edge.get()) { cuts.push_back(r.p); }
                    }
                    std::vector<std::unique_ptr<Stoicheion>> sub;
                    Kern::splitEdgeAtPoints(*edge, cuts, sub);
                    for (auto& s : sub) { pieces.push_back({ std::move(s), a }); }
                }
                range[a] = { first, pieces.size() - first };
            }

            // -- 3+4. per-piece side classification. Crossing, winding, turning and intersection are DIFFERENT things, and only crossing -- which side of each other loop this span lies on -- decides the frontier. The cluster loops are simple (post-orient, post-cull), so "which side of B" is well-defined everywhere off B's boundary, and a piece lies strictly BETWEEN crossings: its B-side is CONSTANT over its span and readable directly at its own interior points. No tangent signs, no accumulation state -- the crossings' only job was the cuts that bound the pieces. Three samples must agree; disagreement = a crossing was missed inside this piece (a straddler) -> fall back.
            std::vector<int> level(pieces.size(), 0);
            bool levelOk = true;
            for (size_t a : mem) {
                const auto [first, count] = range[a];
                if (count == 0) { continue; }

                int h = loops[a].handedness();
                if (h == 0) { h = loops[a].turningSign(); }
                if (h == 0) { levelOk = false; break; }
                const int ownLeft = (h > 0) ? 1 : 0;   // own contribution left of travel

                for (size_t k = first; k < first + count && levelOk; k++) {
                    int w[3]; std::size_t wi = 0;
                    for (float t : { 0.25f, 0.5f, 0.75f }) {
                        const Pos q = Kern::edgePointAt(*pieces[k].e, t);
                        int sum = 0;
                        for (size_t b : mem) { if (b != a) { sum += loops[b].windingAround(q); } }
                        w[wi++] = sum;
                    }
                    if (w[0] != w[1] || w[1] != w[2]) { levelOk = false; break; }   // straddler: missed cut
                    level[k] = ownLeft + w[1];
                }
                if (!levelOk) { break; }
            }

            // -- 5. keep the frontier; on a parity violation fall back to the proven side-probe classification for this cluster
            std::vector<std::unique_ptr<Stoicheion>> kept;
            if (levelOk) {
                for (auto& tp : pieces) {
                    if (level[&tp - pieces.data()] == 1) {
                        tp.e->id = loops[tp.loop].id;   // provenance
                        kept.push_back(std::move(tp.e));
                    }
                }
            }
            else {
                for (auto& tp : pieces) {
                    const Pos mid = Kern::edgeMidpoint(*tp.e);
                    const Pos t = Kern::travelDirAt(*tp.e, mid);
                    float clearance = 1e30f;
                    for (size_t b : mem) {
                        for (const auto& le : loops[b].edges) {
                            clearance = std::min(clearance, le->distanceTo(mid));
                        }
                    }
                    const float reach = std::max(1e-4f, std::min(1e-2f * Kern::edgeLength(*tp.e), 0.4f * clearance));
                    int wl = 0, wr = 0;
                    for (size_t b : mem) {
                        wl += loops[b].windingAround(mid + Kern::perpCCW(t) * reach);
                        wr += loops[b].windingAround(mid - Kern::perpCCW(t) * reach);
                    }
                    // doctrine keep: material strictly on the LEFT of travel
                    if (wl >= 1 && wr < 1) {
                        tp.e->id = loops[tp.loop].id;
                        kept.push_back(std::move(tp.e));
                    }
                }
            }

            // Stitch the surviving pieces back into loops. Oriented pieces: reversal is forbidden -- if ends don't meet forward, the piece set is wrong and must FAIL VISIBLY (an open remnant), never be silently flipped into a wrong-handed loop.
            while (!kept.empty()) {
                Chain stitched = buildOne(kept, eps, /*allowReverse=*/false);
                if (stitched.edges.empty()) { break; }

                size_t dominant = mem.front(); float bestLen = -1.0f;
                for (size_t a : mem) {
                    float len = 0.0f;
                    for (const auto& e : stitched.edges) {
                        if (e->id == loops[a].id) { len += Kern::edgeLength(*e); }
                    }
                    if (len > bestLen) { bestLen = len; dominant = a; }
                }
                stitched.id = newId();
                stitched.parent = loops[dominant].parent;
                stitched.role = Role::Minuend;   // a stitched frontier bounds material by construction

                out.push_back(std::move(stitched));
            }
        }

        // Geometric-significance floor on the OUTPUT too: stitched cusp slivers (near-zero area, 2-3 edges) are noise the parity cull cannot judge (no handedness) and the next offsetValid would immortalize. See the input floor above.
        std::erase_if(out, [](const Chain& m) {
            return m.closed && std::fabs(m.signedArea()) < SignificanceAreaFloor;
        });

        // Geometric-significance floor: stitched cusp slivers are noise the parity cull cannot judge (no handedness) and the next offsetValid would immortalize (see SignificanceAreaFloor).
        std::erase_if(out, [](const Chain& m) {
            return m.closed && std::fabs(m.signedArea()) < SignificanceAreaFloor;
        });

        // The parity cull: enforce nesting alternation on the RESOLVED output. A frontier can never legitimately carry a ring nested inside a same-handed ring; any such survivor is a doubled bound or an artifact and corrupts every downstream winding. Cull, never reorient.
        cullSameHandedNested(out);
        return out;
    }

    // Profile -- one GENERATION: the chains bounding the clearable material at one offset depth. The complete method, one generation forward: stage A per chain (offsetValid), stage B across chains (merge). An exhausted profile returns empty -- extinction, not error.
    struct Profile {

        std::vector<Chain> chains;

        // The sanity verdict: a generation whose total signed area goes negative claims an unbounded region (transient garbage -- e.g. a surviving inverted burr). It is recorded and SKIPPED by consumers, never emitted -- and it never halts the ladder: recursion continues and the garbage extinctifies on its own (the legacy engine's concession, kept).
        bool included = true;

        Profile() = default;
        Profile(Profile&&) = default;
        Profile& operator=(Profile&&) = default;

        Profile clone() const {
            Profile p;
            p.included = included;
            for (const Chain& c : chains) { p.chains.push_back(c.clone()); }
            return p;
        }

        bool empty() const { return chains.empty(); }

        float signedArea() const {
            float a = 0.0f;
            for (const Chain& c : chains) { if (c.closed) { a += c.signedArea(); } }
            return a;
        }

        Profile offsetBy(float amount, float eps = 1e-3f) const {
            Profile next;
            std::vector<Chain> valid;
            for (const Chain& c : chains) {
                if (!c.closed) { continue; }
                for (Chain& m : c.offsetValid(amount, eps)) { valid.push_back(std::move(m)); }
            }
            next.chains = Chain::merge(std::move(valid), eps);

            // The cascade contract: NEGATIVE LOBES ALWAYS GET FILTERED. A negative (CW) chain is legitimate ONLY as a hole -- contained in a positive sibling (an annulus hole-family ring). A lone CW chain with no positive container bounds no material; un-filtered it cascades as a growing inverted ring and EMITS while positive geometry keeps its generation's sum positive. Per-chain judgment, not the generation sum.
            std::erase_if(next.chains, [&next](const Chain& c) {
                if (!c.closed) { return false; }
                int h = c.handedness();
                if (h == 0) { h = c.turningSign(); }
                if (h >= 0) { return false; }
                for (const Chain& b : next.chains) {
                    if (&b == &c || !b.closed) { continue; }
                    int hb = b.handedness();
                    if (hb == 0) { hb = b.turningSign(); }
                    if (hb > 0 && Chain::containedIn(c, b)) { return false; }   // a real hole
                }
                return true;                                                    // lone inverted: garbage
            });
            return next;
        }
    };

    // Raster -- parallel strokes filling a region (normally the SEED: the legal tool-centre region, so clipped strokes are safe by construction). Each scan line's transversal crossings pair into interior spans by even-odd; strokes come back OPEN, ordered and DIRECTED: `serpentine` reverses alternate lines end-for-end (hatch -- turns at the ends); without it every stroke runs the same way (scallop -- constant engagement sense, the consumer retracts between strokes). Links, retracts and Z are the consumer's grammar, never woven here.
    inline std::vector<Chain> raster(const Profile& region, float angleDeg, float spacing,
                                     bool serpentine, bool reverse = false) {

        std::vector<Chain> out;

        const float ang = angleDeg * (TAU / 360.0f);
        const Pos d(std::cos(ang), std::sin(ang));    // strokes run along d
        const Pos p(-std::sin(ang), std::cos(ang));   // scan lines advance along p

        // The clip edges and the region's extent in the (d, p) frame
        std::vector<const Stoicheion*> clip;
        float dMin = 1e30f, dMax = -1e30f, pMin = 1e30f, pMax = -1e30f;
        for (const Chain& c : region.chains) {
            if (c.edges.empty()) { continue; }
            for (const auto& e : c.edges) { clip.push_back(e.get()); }
            std::vector<Pos> pts;
            c.tessellate(pts);
            for (const Pos& q : pts) {
                const float qd = q.dot(d), qp = q.dot(p);
                dMin = std::min(dMin, qd); dMax = std::max(dMax, qd);
                pMin = std::min(pMin, qp); pMax = std::max(pMax, qp);
            }
        }
        if (clip.empty() || dMax <= dMin || pMax <= pMin) { return out; }

        spacing = std::max(spacing, 1e-3f);
        const float pad = spacing + 1.0f;

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
                Kern::edgeCross(scan, *e, hits);
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

        if (reverse) { std::reverse(lines.begin(), lines.end()); }   // scan-advance order

        for (std::size_t li = 0; li < lines.size(); li++) {
            std::vector<Stroke>& strokes = lines[li];
            if (serpentine && (li % 2 == 1)) {   // hatch: alternate lines end-for-end
                std::reverse(strokes.begin(), strokes.end());
                for (Stroke& st : strokes) { std::swap(st.a, st.b); }
            }
            for (const Stroke& st : strokes) {
                Chain c;
                c.closed = false;
                c.id = newId();
                c.edges.push_back(std::make_unique<Segment2>(st.a, st.b));
                out.push_back(std::move(c));
            }
        }
        return out;
    }

    // Parallel -- the seed construction and the inlay (Method/06). Input: two sets of raw stoicheia. FILL encloses area where cutting is encouraged (the tool must maximally touch it); AVOID encloses area where no part of the tool may enter. Both arrive with NO orientation promise. The strategy: 1. assemble each set into endpoint-coincident chains; 2. orient by NESTING PARITY with opposite bases -- fill outers CCW (islands CW), avoid outers CW (islands CCW): one shared nesting tree, handedness following depth parity, the avoid being the fill's island (Method/05, 06); 3. GROW BOTH by one tool radius -- fill by offsetLeft(-R), avoid by offsetLeft(+R); opposite signs, both regions expanding; 4. merge: the seed is the frontier of { winding >= 1 } = (fill (+) R) \ (avoid (+) R) -- the maximum extent of the tool CENTRE that reaches the fill and never violates the avoid; 5. inlay: recursively inset the seed by the stepover until extinction. profiles[0] IS the outermost cuttable ring (the one-radius standoff was built into the seed, so there is no "offset by one radius first" -- that rule dissolved into the construction).
    struct Parallel {

        struct Params {
            float toolRadius = 1.0f;
            float stepover = 1.0f;        // generation advance as a FRACTION of the tool radius
            int   maxGenerations = 64;    // inlay cap
            float eps = 1e-3f;            // KERNEL tolerance: crossings, stitching, fragment
                                          // coincidence -- pure geometry math, never input slop
            float weldEps = 1e-2f;        // almost-closed loops within this gap are welded shut
                                          // (BRep-section endpoints carry approximation slop)

            // ASSEMBLY tolerances, PER SYSTEM, derived from each sectioned solid's own tolerance (max vertex tolerance; a spec solid ~1e-7, a healed / measured solid ~1e-2). Two endpoints within a solid's own tolerance are the SAME point per its model, so assembling there closes measured-fragmented loops without ever fusing features the solid keeps genuinely apart. The cut and the avoid are DIFFERENT solids (a spec delta against a measured part is routine), so each system assembles at its own figure. 0 = fall back to eps. Assembly tolerance is input conditioning; the kernel eps above never widens.
            float fillEps = 0.0f;
            float avoidEps = 0.0f;

            // The extended (lead) region's two expansions, as FRACTIONS of the tool radius. The survey of the open air lead motion may actually use: fill grown by a full diameter (centres within feasible reach of the cut), avoid grown PAST the ring standoff -- rings graze avoid + R, but approach moves keep an extra half radius off kept walls (over-avoidance is the safe direction; the seed's exact-R bound is untouched). leadReach = 0 disables the construction.
            float leadReach = 2.0f;   // extended fill expansion
            float leadClear = 1.5f;   // extended avoid expansion (never below 1.0)

            // The TRUE-open-air survey's expansions (fractions of R): fill by one full tool diameter, avoid by two thirds of one. An independent, parallel construction (nothing in the seed/inlay pipeline reads it); its product is the openAir member. 0 disables.
            float airReach = 2.0f;
            float airClear = 4.0f / 3.0f;

            // Generic-position symmetry breaker. The fill is grown by R + fillBias while the avoid is grown by EXACTLY R -- so grown boundaries that would otherwise lie coincident for long stretches (a feature wall a whole multiple of R from a kept wall -- tool-radius-dependent, so it strikes "randomly" per tool) are separated in guaranteed transversal position, beyond every kernel tolerance and the section slop. The asymmetry is the principle, not a fudge: the avoid inflation is the SAFETY BOUND (dimensional truth, exactly R -- never biased), while the fill expansion is only a REACH CLAIM, and over-claiming reach is harmless -- the claim is still clipped by the exact avoid term, and the extra overrun hangs in open air. Greed on the minuend is free; the subtrahend must stay exact.
            float fillBias = 0.05f;       // model units (absolute: kernel tolerances are absolute)

            // Seed reach knob, as a MULTIPLE of the tool radius. The base fill grow is seedOffset*R + fillBias (default 1 = exactly R + fillBias, the plain tool-radius dilation). Larger claims more reach -- the tool-centre region extends further past the true material, buying open-air room. Reach only: the grow is still clipped by the exact avoid term in the merge (avoid stays +R, the hard keep-out), so over-claiming here only hangs extra centres in air, never drives the tool into kept material. Only the BASE fill reads it; the lead/air reach variants keep their own factors.
            float seedOffset = 1.0f;

            // Clip-seed FILL reach, a MULTIPLE of R (a hair over 1). This is a generic-position PROXY of the true tool-centre wall -- exact-boundary clipping (fill (+) exactly R) is numerically fragile (Weiler-Atherton on coincident edges), so the clip region is displaced a hair OUTWARD (0.1R) from the true wall. When seedOffset runs FAT (corner relief / early exit into free-edge air), the cascade rings are GENERATED from the fat seed but must OBEY this proxy of the true tool-centre region -- (fill (+) clipSeedOffset*R) \ (avoid (+) clipAvoidOffset*R). The strategy clips each fat CASCADE ring to it, dropping the far-air arclength. This proxy is the FINISH clip target; the strategy insets it by wallStock to form the ROUGH clip target (the inset-seed proxy). The 0.1R over 1R displaces the clip boundary past kernel eps from any fat-seed-generated edge (generic position, same rationale as fillBias). wallStock is NEVER on this seed -- it re-enters only on the ROUGH clip proxy, strategy-side. Only built when the fat seed is meaningfully larger (seedOffset > clipSeedOffset + eps); otherwise the two regions coincide and the clip is a no-op -- clipSeed stays empty.
            float clipSeedOffset = 1.1f;

            // Clip-seed AVOID grow, a MULTIPLE of R (a hair UNDER 1). The FAT seed grows avoid by exactly +R (the real hard keep-out), so a fat ring's legitimate wall-finish stretch (tool tangent to the wall) lies exactly on avoid (+) R. If the clip seed also grew avoid by +R, that stretch would sit ON the clip's avoid boundary, where winding classification is ill-defined -- the wall pass could be dropped and the wall left uncut. Growing the clip's avoid a hair LESS (0.9R) pushes its boundary 0.1R inside the fat wall stretch, so that stretch is UNAMBIGUOUSLY kept (generic position, mirroring the fill side's displacement). Safe: the fat seed already enforced the true +R limit, so a looser clip avoid can only KEEP already- legal wall pieces, never admit a path into avoid (+) R.
            float clipAvoidOffset = 0.9f;

            // DEPRECATED / UNUSED by the engine. Wall stock is NOT a seed term: wallStock is a BINARY FINAL-EMISSION choice made STRATEGY-side. The engine builds two seed regions -- the (fat) offset seed (profiles[0], the free-form cascade origin) and the TRUE seed (trueSeed, the real tool-centre wall at cut R / avoid R) -- and exposes trueSeed so the strategy emits its outermost ring as either the true seed (finish) or the true seed inset by the wall skin (rough). No cascade geometry ever carries the skin. Retained only so existing initializers compile; the engine seed no longer reads it.
            float wallStock = 0.0f;
        };

        Params params;

        Parallel() = default;
        explicit Parallel(Params p) : params(p) {}

        // Ingested input (consumed by run).
        std::vector<std::unique_ptr<Stoicheion>> fill;
        std::vector<std::unique_ptr<Stoicheion>> avoid;

        // Output: profiles[0] is the seed; each later profile is one stepover deeper. Warnings are plain-language deviations, empty on a clean run.
        std::vector<Profile> profiles;
        std::vector<std::string> warnings;
        double timingPreparationMs = 0.0;
        double timingInlayMs = 0.0;

        // The oriented, welded FILL section (material boundary at this slice), kept after run() for consumers: winding >= 1 against these chains is material, below is open air -- the test that decides where a ring may be entered by a straight plunge (the tool centre column hangs in air).
        std::vector<Chain> fillBoundary;

        // The oriented, welded AVOID section, likewise kept -- primarily so a consumer can SEE what actually bounded the seed (an empty list on a level with kept material is the tell that the avoid never arrived).
        std::vector<Chain> avoidBoundary;

        // The grown systems, preserved pre-merge (the merge consumes its inputs): the per-chain offsetValid results of fill and avoid.
        std::vector<Chain> grownFill;
        std::vector<Chain> grownAvoid;

        // The EXTENDED region: (fill (+) leadReach*R) minus (avoid (+) leadClear*R) -- the survey of the open air lead motion may actually use, built by the seed's own construction on the same preserved originals. Membership is a wall-safety certificate WITH MARGIN: a lead arc lying wholly inside it stays within feasible reach of the cut and keeps an extra half radius off every kept wall (rings graze avoid + R; leads don't).
        Profile extended;

        // The CLIP seed: the true tool-centre region -- (fill (+) clipSeedOffset*R) \ (avoid (+) clipAvoidOffset*R) -- that a FAT-seed run must obey. Rings are generated from the reach-inflated seed (corner relief, early free- edge exit) then clipped back to this region so the far-air arclength is dropped. The avoid term is a hair UNDER +R (clipAvoidOffset < 1) so a fat wall-hugging stretch (at the fat seed's real +R limit) falls inside this clip and is kept, not dropped on an ill-defined boundary. The true +R keep-out is still enforced by the FAT seed, so the looser clip avoid only KEEPS legal travel, never licenses the tool past a wall. Empty unless the fat seed is meaningfully larger than the clip (see clipSeedOffset).
        Profile clipSeed;

        // The TRUE seed: the real tool-centre wall region -- (fill (+) R) \ (avoid (+) R), i.e. seedOffset = 1 and NO wall stock. Built ALWAYS (ungated), by the SAME seed machinery as profiles[0] but at the true radius, from the same preserved originals. This is the AUTHORITATIVE dimensional wall: the strategy emits it (finish) or its uniform inset by the wall skin (rough) as the single outermost ring, and bounds link containment against it. profiles[0] (the possibly-fat offset seed) stays the free-form cascade origin; when seedOffset == 1 the two coincide, so the strategy drops profiles[0] and emits trueSeed once (no duplicate).
        Profile trueSeed;

        // TRUE OPEN AIR -- the entry planner's menu of certified sites. The seed treatment run once more at (fill (+) airReach*R) \ (avoid (+) airClear*R), and the resulting boundary CLASSIFIED against the original sections: pieces lying outside BOTH regions are kept, as OPEN chains (segments of the air ring). On such a segment the whole tool disk hangs a full radius clear of the cut material and a third of a diameter off every kept wall -- a plunge there touches nothing. (The extended region cannot certify this: fill (+) 2R CONTAINS the fill, so membership there proves reach and wall safety, never disk clearance.)
        std::vector<Chain> openAir;

        // Pre-oriented avoid path: when the caller has already assembled and oriented the keep-out as closed loops, feed them whole and SKIP the engine's assemble/weld/dropOpen/orientByParity for the avoid side.
        std::vector<Chain> avoidLoops;
        bool avoidPreOriented = false;

        void ingest(std::vector<std::unique_ptr<Stoicheion>> fillSet,
                    std::vector<std::unique_ptr<Stoicheion>> avoidSet) {
            fill = std::move(fillSet);
            avoid = std::move(avoidSet);
        }

        // Opt-in: supply the avoid as pre-oriented CLOSED loops. INVARIANT -- every loop MUST be closed with role == Role::Subtrahend; the weld / dropOpen safety net and the joint orientByParity are all skipped for these, so an open or wrong-role loop is not corrected here.
        void ingestAvoidLoops(std::vector<Chain> loops) {
            avoidLoops = std::move(loops);
            avoidPreOriented = true;
        }

        // Orient a system by nesting parity: depth = how many OTHER loops of the same system contain this one (probed at the loop's own extreme point -- on the loop, far from the others unless they truly touch). Desired handedness alternates from `base` with depth. Judged and forced by exact SIGNED AREA (input orientation is a seed question; turning is the coin-flip on tangency-cusped slivers). Degenerate zero-area loops are dropped with a warning.
        void orientByParity(std::vector<Chain>& loops, int base, Role role) {

            // DECIDE FIRST, MOVE SECOND. The old single pass std::moved each processed loop into `kept`, so every b < a was a hollow husk during a's containment test -- depth silently read 0 whenever a container happened to sort before its contained loop. Assembly order is arbitrary and mirror-dependent, which made this the one-sided island bug: same feature, opposite faces, opposite verdicts.
            const size_t n = loops.size();
            std::vector<int> cur(n, 0);
            std::vector<int> desired(n, 0);

            for (size_t a = 0; a < n; a++) {
                float area = loops[a].signedArea();
                cur[a] = (area > 1e-4f) ? 1 : (area < -1e-4f ? -1 : 0);
                if (cur[a] == 0) { continue; }

                int depth = 0;
                for (size_t b = 0; b < n; b++) {
                    if (b == a) { continue; }
                    if (Chain::containedIn(loops[a], loops[b])) { depth++; }
                }
                desired[a] = ((depth % 2) == 0) ? base : -base;

                {   // orientation confession: names WHY a loop kept/flipped its sense
                    char buf[160];
                    std::snprintf(buf, sizeof(buf),
                        "[orientByParity] loop #%zu area %+.2f depth %d desired %+d%s",
                        a, static_cast<double>(area), depth, desired[a],
                        (cur[a] == desired[a]) ? "" : " FLIP");
                    warnings.push_back(buf);
                }
            }

            std::vector<Chain> kept;
            for (size_t a = 0; a < n; a++) {
                if (cur[a] == 0) {
                    warnings.push_back("A degenerate (zero-area) input loop was dropped.");
                    continue;
                }
                Chain c = (cur[a] == desired[a]) ? std::move(loops[a]) : loops[a].reversed();
                c.role = role;
                kept.push_back(std::move(c));
            }
            loops = std::move(kept);
        }

        // Run the whole construction: assemble -> orient -> grow both -> merge (the seed) -> inlay to extinction. Returns false (with a warning) when no seed exists -- the avoid swallows the fill, or the fill never assembled into a closed region.
        bool run() {

            using Clock = std::chrono::steady_clock;
            const auto runStarted = Clock::now();
            auto elapsedMs = [](Clock::time_point a, Clock::time_point b) {
                return std::chrono::duration<double, std::milli>(b - a).count();
            };
            timingPreparationMs = 0.0;
            timingInlayMs = 0.0;

            profiles.clear();
            warnings.clear();
            fillBoundary.clear();
            avoidBoundary.clear();
            grownFill.clear();
            grownAvoid.clear();
            extended = Profile();
            clipSeed = Profile();
            openAir.clear();
            // NB: do NOT clear avoidPreOriented/avoidLoops here -- ingestAvoidLoops runs before run(). They are consumed (moved-from) below; a fresh Parallel per level means no stale carry-over.

            const float R = params.toolRadius;
            if (R <= 1e-6f) { warnings.push_back("Tool radius is zero: nothing to construct."); return false; }

            // 1. Assemble. Closed chains only past this door: offsetting is a statement about bounded regions (welding almost-closed input is the caller's conditioning job).
            const float fEps = std::max(params.eps, params.fillEps);
            const float aEps = std::max(params.eps, params.avoidEps);
            std::vector<Chain> fills  = Chain::assemble(fill, fEps);
            std::vector<Chain> avoids = avoidPreOriented ? std::move(avoidLoops)
                                                         : Chain::assemble(avoid, aEps);

            // HARD-weld almost-closed loops shut: section endpoints carry slop past the assembly eps, and a loop that is closed in fact but open by a hair must not be dropped at the door. Closure is enforced by SNAPPING both ends to their mean -- never a micro bridge segment (a sub-eps edge is invisible to crossing detection). Then every INTERIOR joint is snapped exact too: assemble connects within eps, so a "closed" chain still carries up-to-eps slop at every vertex until Chain::hardWeld runs.
            auto weld = [&](std::vector<Chain>& v, float sysEps) {
                const float gapEps = std::max(params.weldEps, sysEps);   // a loop's own residual
                for (Chain& c : v) {                                     // gap scales with its solid
                    if (c.edges.empty()) { continue; }
                    if (!c.closed) {
                        const Pos a = Kern::eEnd(*c.edges.back());
                        const Pos b = Kern::eStart(*c.edges.front());
                        const float gap = (b - a).pythag();
                        if (gap > gapEps) { continue; }
                        if (gap > 1e-9f) {
                            const Pos mid = (a + b) * 0.5f;
                            if (!Kern::setEdgeEnd(*c.edges.back(), mid)
                                || !Kern::setEdgeStart(*c.edges.front(), mid)) {
                                c.edges.push_back(std::make_unique<Segment2>(
                                    Kern::eEnd(*c.edges.back()), Kern::eStart(*c.edges.front())));
                            }
                        }
                        c.closed = true;
                    }
                    Chain::hardWeld(c, gapEps);
                }
            };
            weld(fills, fEps);
            if (!avoidPreOriented) { weld(avoids, aEps); }   // pre-oriented loops are already closed

            // Dropped opens are the assembly's failure report: say HOW each one failed (edge count + closure gap), so a consumer can tell endpoint slop (gap just past the weld) from garbage sections (huge gaps -- the tell of a tilted or drifted slicing frame).
            auto dropOpen = [&](std::vector<Chain>& v, const char* which) {
                int reported = 0;
                for (const Chain& c : v) {
                    if (c.closed || c.edges.empty()) { continue; }
                    if (reported++ >= 4) { break; }
                    const float gap = (Kern::eStart(*c.edges.front()) - Kern::eEnd(*c.edges.back())).pythag();
                    float length = 0.0f;
                    for (const auto& e : c.edges) { length += Kern::edgeLength(*e); }
                    char buf[160];
                    std::snprintf(buf, sizeof(buf),
                        "Open %s chain dropped: %zu edge(s), length %.3f, closure gap %.4f",
                        which, c.edges.size(), length, gap);
                    warnings.push_back(buf);
                }
                std::erase_if(v, [](const Chain& c) { return !c.closed; });
            };
            dropOpen(fills, "fill");
            if (!avoidPreOriented) { dropOpen(avoids, "avoid"); }   // pre-oriented loops guaranteed closed

            if (fills.empty()) {
                warnings.push_back("No closed fill boundary: nothing to cut.");
                return false;
            }

            // 2. Orient: one nesting tree, opposite bases (Method/05, 06).
            orientByParity(fills,  +1, Role::Minuend);      // fill: outers CCW, islands CW
            // avoid: pre-oriented loops keep the caller's own-parity handedness; re-deriving JOINT parity here would flip deep nested walls (gouge).
            if (!avoidPreOriented) { orientByParity(avoids, -1, Role::Subtrahend); }   // avoid: outers CW, islands CCW

            // 3. Grow BOTH by one tool radius. Fill by -R (left of CCW travel is the interior, so -R expands); avoid by +R (left of CW travel is the exterior, so +R expands). Islands of each system shrink under the same value -- dilation of the region, by one blind sign per system (Method/06). The fill takes the extra fillBias (see Params): reach is claimed greedily, the safety bound stays exact, and grown boundaries never sit coincident. A SHRINKING chain (fill outer, avoid island -- positive signed area in its system's parity) may legitimately collapse under the tool radius. A DILATING chain (fill hole, avoid outer -- negative signed area) grows the region: its offset can NEVER legitimately vanish, so an empty offsetValid there is a kernel failure -- and swallowing it silently drops a hard bound from the merge (an island vanishes from the seed, the cascade runs outer-only, and only the independently-built clip stands between the rings and the kept material). Loud + refusing. Exception: a chain under SignificanceAreaFloor bounds no region and offsetValid refuses it BY DESIGN -- its empty offset is the promised outcome, not a lost bound. Refusing the level on one is a false alarm on merge noise (a cusp sliver is not a hard bound).
            std::vector<Chain> grown;
            bool lostBound = false;
            auto growInto = [&](Chain& c, float amount, std::vector<Chain>& keep,
                                const char* system) {
                std::vector<Chain> off = c.offsetValid(amount, params.eps);
                if (off.empty() && !c.closed) {
                    // an open chain bounds nothing and its signedArea is meaningless -- it must not read as a dilating bound. On the pre-oriented path this is a caller contract breach: loud.
                    char buf[160];
                    std::snprintf(buf, sizeof(buf),
                        "%s OPEN chain reached growth (%zu edge(s)) -- "
                        "contract breach upstream; dropped", system, c.edges.size());
                    warnings.push_back(buf);
                    return;
                }
                if (off.empty() && std::fabs(c.signedArea()) < Chain::SignificanceAreaFloor) {
                    return;   // sub-floor noise: unoffsettable by definition, drop
                }
                if (off.empty() && c.signedArea() < 0.0f) {
                    char buf[160];
                    std::snprintf(buf, sizeof(buf),
                        "%s bound failed to offset (dilating side, %zu edge(s), "
                        "area %+.4f) -- seed would lose a hard bound; level refused",
                        system, c.edges.size(), static_cast<double>(c.signedArea()));
                    warnings.push_back(buf);
                    lostBound = true;
                    return;
                }
                for (Chain& m : off) {
                    keep.push_back(m.clone());          // preserved: the merge consumes its input
                    grown.push_back(std::move(m));
                }
            };
            for (Chain& c : fills) {
                growInto(c, -(params.seedOffset * R + params.fillBias), grownFill, "fill");
            }
            for (Chain& c : avoids) {
                growInto(c, +R, grownAvoid, "avoid");
            }
            if (lostBound) { return false; }

            // The sections have served the growth; keep them as the consumer's oracles (fillBoundary: material vs air; avoidBoundary: what bounded).
            fillBoundary = std::move(fills);
            avoidBoundary = std::move(avoids);

            // 4. The seed: one merge. { winding >= 1 } with roles is exactly (fill (+) R) \ (avoid (+) R). Envelopment-without-intersection (the bore: the grown avoid ring sitting whole inside the grown fill) is handled by the merge's containment clustering -- the avoid's ring survives as the seed, the fill's copy is discarded.
            Profile seed;
            seed.chains = Chain::merge(std::move(grown), params.eps);

            // The seed bounds the clearable region; the construction roles are spent. Everything from here on is minuend.
            for (Chain& c : seed.chains) { c.role = Role::Minuend; }

            if (seed.empty() || seed.signedArea() <= 0.0f) {
                warnings.push_back("No seed: the avoid region leaves no reachable extent of the fill "
                                   "(or the fill collapsed under the tool radius).");
                return false;
            }

            // 4b. Merge validation, path-independent: the interior of every DILATING bound (fill hole, avoid outer -- the hard keep-outs) must read OUTSIDE the merged seed, whichever merge path (crossing walk / containment clustering) handled it. A large island whose grown ring INTERSECTS the fat-inset outer takes the crossing walk; if the merge loses it there, the cascade runs outer-only over kept material with just the clip standing in the way. Probe: nudge inward of the ring's max-X vertex (interior of any simple loop lies -X of its rightmost point), verified against the ring's own winding before it may condemn the seed.
            {
                auto retainsInterior = [&](const Chain& c) -> bool {
                    if (!c.closed || c.edges.empty() || c.signedArea() >= 0.0f) { return false; }
                    Pos mx = Kern::eStart(*c.edges.front());
                    for (const auto& e : c.edges) {
                        const Pos a = Kern::eStart(*e);
                        const Pos b = Kern::eEnd(*e);
                        if (a.x > mx.x) { mx = a; }
                        if (b.x > mx.x) { mx = b; }
                    }
                    for (float nudge : { 0.02f * R, 0.1f * R, 0.5f * R }) {
                        const Pos probe{ mx.x - nudge, mx.y };
                        if (c.windingAround(probe) == 0) { continue; }   // not interior; try deeper
                        int w = 0;
                        for (const Chain& s : seed.chains) { w += s.windingAround(probe); }
                        return w >= 1;                                    // kept-material interior inside the seed
                    }
                    return false;
                };
                for (const Chain& c : grownAvoid) {
                    if (retainsInterior(c)) {
                        warnings.push_back("Seed retains an avoid bound's interior -- the merge lost "
                                           "an intersecting island ring; level refused.");
                        return false;
                    }
                }
                for (const Chain& c : grownFill) {
                    if (retainsInterior(c)) {
                        warnings.push_back("Seed retains a fill hole's interior -- the merge lost an "
                                           "intersecting hole ring; level refused.");
                        return false;
                    }
                }
            }

            // The extended region: the identical construction from the same preserved originals, at the lead expansions (see Params). The fill side carries the bias (generic position, greed is free); the avoid side is a clean 1.5R -- above R is the safe direction.
            if (params.leadReach > 1e-4f) {
                const float extF = params.leadReach * R + params.fillBias;
                const float extA = std::max(1.0f, params.leadClear) * R;
                std::vector<Chain> g2;
                for (const Chain& c : fillBoundary) {
                    for (Chain& m : c.offsetValid(-extF, params.eps)) { g2.push_back(std::move(m)); }
                }
                for (const Chain& c : avoidBoundary) {
                    for (Chain& m : c.offsetValid(+extA, params.eps)) { g2.push_back(std::move(m)); }
                }
                extended.chains = Chain::merge(std::move(g2), params.eps);
                for (Chain& c : extended.chains) { c.role = Role::Minuend; }
            }

            // The clip seed: the seed's OWN construction once more, at the true tool radius (clipSeedOffset, a hair over 1R) rather than the fat seedOffset. The fill side carries the same fillBias (generic position); the avoid side grows by clipAvoidOffset*R -- a hair UNDER +R, so the fat seed's wall-hugging stretches (at avoid (+) R, the real limit the fat seed enforced) sit unambiguously INSIDE this clip region and are kept, instead of coincident with its boundary (where winding is ill-defined). The true +R keep-out is still enforced by the FAT seed; this looser clip avoid can only keep already-legal wall pieces, never admit avoid (+) R. Built only when the fat seed is meaningfully larger; when they coincide the clip would be a no-op, so it is skipped and clipSeed left empty.
            if (params.clipSeedOffset > 1e-4f
                && params.seedOffset > params.clipSeedOffset + 0.05f) {
                const float clipF = params.clipSeedOffset * R + params.fillBias;
                const float clipA = params.clipAvoidOffset * R;
                std::vector<Chain> gC;
                for (const Chain& c : fillBoundary) {
                    for (Chain& m : c.offsetValid(-clipF, params.eps)) { gC.push_back(std::move(m)); }
                }
                for (const Chain& c : avoidBoundary) {
                    for (Chain& m : c.offsetValid(+clipA, params.eps)) { gC.push_back(std::move(m)); }
                }
                clipSeed.chains = Chain::merge(std::move(gC), params.eps);
                for (Chain& c : clipSeed.chains) { c.role = Role::Minuend; }
            }

            // The TRUE seed: the seed's OWN construction once more, at the true tool radius (seedOffset = 1) and NO wall stock -- (fill (+) R) \ (avoid (+) R). Built UNCONDITIONALLY (not gated on a fat seed): it is the authoritative dimensional wall the strategy emits as its outermost ring (finish = this; rough = this inset by the wall skin) and bounds link containment against. The fill side keeps the same fillBias for generic position; the avoid side is a clean +R (the real hard keep-out). Distinct from clipSeed, whose avoid is a hair UNDER +R for winding-safe clipping -- this is the exact wall.
            if (std::fabs(params.seedOffset - 1.0f) <= 1e-6f) {
                // The default seed construction is already exactly (fill +(R+bias)) \ (avoid +R). Reuse that result instead of performing the same offsets and winding merge a second time.
                for (const Chain& c : seed.chains) { trueSeed.chains.push_back(c.clone()); }
            }
            else {
                const float trueF = R + params.fillBias;
                std::vector<Chain> gT;
                for (const Chain& c : fillBoundary) {
                    for (Chain& m : c.offsetValid(-trueF, params.eps)) { gT.push_back(std::move(m)); }
                }
                for (const Chain& c : avoidBoundary) {
                    for (Chain& m : c.offsetValid(+R, params.eps)) { gT.push_back(std::move(m)); }
                }
                trueSeed.chains = Chain::merge(std::move(gT), params.eps);
                for (Chain& c : trueSeed.chains) { c.role = Role::Minuend; }
            }

            // The true-open-air survey: the identical construction once more, at the air expansions, from the same preserved originals. The merged ring is part over the cut region (its avoid-derived stretches can run inside the fill) and part in genuine air; each edge is split at its crossings with the ORIGINAL section boundaries and the pieces classified by winding -- outside the fill region (w < 1, CCW base) AND outside the avoid region (w > -1, CW base) is true air. The survivors stitch into open segments, oriented, never reversed.
            if (params.airReach > 1e-4f) {

                const float airF = params.airReach * R + params.fillBias;
                const float airA = params.airClear * R;
                std::vector<Chain> g3;
                for (const Chain& c : fillBoundary) {
                    for (Chain& m : c.offsetValid(-airF, params.eps)) { g3.push_back(std::move(m)); }
                }
                for (const Chain& c : avoidBoundary) {
                    for (Chain& m : c.offsetValid(+airA, params.eps)) { g3.push_back(std::move(m)); }
                }
                std::vector<Chain> ring = Chain::merge(std::move(g3), params.eps);

                std::vector<std::unique_ptr<Stoicheion>> kept;
                for (Chain& rc : ring) {
                    for (auto& e : rc.edges) {

                        std::vector<Pos> cuts;
                        for (const Chain& c : fillBoundary) {
                            for (const auto& oe : c.edges) { Kern::edgeCross(*e, *oe, cuts); }
                        }
                        for (const Chain& c : avoidBoundary) {
                            for (const auto& oe : c.edges) { Kern::edgeCross(*e, *oe, cuts); }
                        }

                        std::vector<std::unique_ptr<Stoicheion>> pieces;
                        Kern::splitEdgeAtPoints(*e, cuts, pieces);

                        for (auto& piece : pieces) {
                            const Pos mid = Kern::edgeMidpoint(*piece);
                            int wF = 0, wA = 0;
                            for (const Chain& c : fillBoundary)  { wF += c.windingAround(mid); }
                            for (const Chain& c : avoidBoundary) { wA += c.windingAround(mid); }
                            if (wF < 1 && wA > -1) { kept.push_back(std::move(piece)); }
                        }
                    }
                }

                while (!kept.empty()) {
                    Chain seg = Chain::buildOne(kept, params.eps, /*allowReverse=*/false);
                    if (seg.edges.empty()) { break; }
                    seg.id = newId();
                    openAir.push_back(std::move(seg));
                }
            }

            profiles.push_back(std::move(seed));

            timingPreparationMs = elapsedMs(runStarted, Clock::now());
            const auto inlayStarted = Clock::now();

            // 5. The inlay: recursively inset until extinction or the cap. profiles[0] is already the outermost CUTTABLE ring -- the seed stands one radius off every wall by construction -- so every advance is a plain stepover. A generation whose total signed area goes non-positive claims an unbounded region (outsides without insides): extinction, per the sanity rule.
            const float step = R * std::max(0.05f, params.stepover);
            while (static_cast<int>(profiles.size()) <= params.maxGenerations) {
                Profile next = profiles.back().offsetBy(step, params.eps);
                if (next.empty()) { break; }
                next.included = next.signedArea() > 0.0f;   // skip downstream, never halt
                profiles.push_back(std::move(next));
            }
            timingInlayMs = elapsedMs(inlayStarted, Clock::now());
            return true;
        }
    };

    // The geometry-side result consumed by a motion strategy.  `parent` indexes
    // the FINAL, already-permuted ring vector. In ring travel the ordering is
    // linkStart -> linkCenter -> linkEnd. Emission begins at linkEnd, traverses
    // the complementary span to linkStart, then links start -> next end.
    // Derived views never alter the canonical ring.
    struct RingPlan {
        std::vector<int> parent;
        std::vector<float> linkCenterDistance;
        std::vector<float> linkStartDistance;
        std::vector<float> linkEndDistance;
        std::vector<Pos> linkCenter;
        std::vector<Pos> linkStart; // outbound departure; BEFORE center in ring travel
        std::vector<Pos> linkEnd;   // inbound landing / this ring's lap entry
        std::vector<Geo::Chain> cycle;   // full cyclic view seated at linkEnd
        std::vector<Geo::Chain> lap;     // directed traversal linkEnd -> linkStart
        std::vector<Geo::Chain> overlap;
        std::vector<char> linkToNext;    // exact projected stay-down link exists
    };

    inline bool inside(const std::vector<Geo::Chain>& region, Pos p) {
        int w = 0;
        for (const Geo::Chain& c : region) { w += c.windingAround(p); }
        return w >= 1;
    }

    // Intersect one cascade loop with the legal seed while retaining a LOOP.
    // Pieces of the cascade inside the seed and pieces of the seed inside the
    // cascade are consumed into one pool; Geo::Chain::buildOne now hard-welds
    // every accepted joint before the loop escapes this layer.
    inline std::vector<Geo::Chain> clipRing(const Geo::Chain& ring,
                                            const std::vector<Geo::Chain>& clip,
                                            float eps, bool& pureBoundary) {
        pureBoundary = false;
        if (ring.edges.empty()) { return {}; }
        if (clip.empty()) {
            std::vector<Geo::Chain> out;
            out.push_back(ring.clone());
            return out;
        }

        // Boolean intersection is region topology, not milling travel. Normalize
        // this one bounded ring to material-left (CCW) before mixing it with the
        // already region-oriented clip system. The source travel sense is restored
        // only after complete result loops exist.
        Geo::Chain workRing = ring.turningSign() < 0 ? ring.reversed() : ring.clone();

        std::vector<std::unique_ptr<Stoicheion>> pool;
        bool ringOutside = false, keptRingPiece = false;
        for (const auto& e : workRing.edges) {
            std::vector<Pos> cuts;
            for (const Geo::Chain& c : clip) {
                for (const auto& b : c.edges) { Geo::Chain::edgeCross(*e, *b, cuts); }
            }
            std::vector<std::unique_ptr<Stoicheion>> pieces;
            Geo::Chain::splitEdgeAtPoints(*e, cuts, pieces);
            for (auto& p : pieces) {
                if (inside(clip, Geo::Chain::edgeMidpoint(*p))) {
                    keptRingPiece = true;
                    pool.push_back(std::move(p));
                }
                else { ringOutside = true; }
            }
        }
        if (!ringOutside) {
            std::vector<Geo::Chain> out;
            out.push_back(ring.clone());
            return out;
        }

        for (const Geo::Chain& boundary : clip) {
            for (const auto& e : boundary.edges) {
                std::vector<Pos> cuts;
                for (const auto& b : workRing.edges) { Geo::Chain::edgeCross(*e, *b, cuts); }
                std::vector<std::unique_ptr<Stoicheion>> pieces;
                Geo::Chain::splitEdgeAtPoints(*e, cuts, pieces);
                for (auto& p : pieces) {
                    if (workRing.windingAround(Geo::Chain::edgeMidpoint(*p)) != 0) {
                        p->openAir = true;
                        pool.push_back(std::move(p));
                    }
                }
            }
        }

        std::vector<Geo::Chain> loops;
        while (!pool.empty()) {
            // All pieces now carry compatible region orientation. Reversing an
            // individual edge here would manufacture a split family.
            Geo::Chain loop = Geo::Chain::buildOne(pool, eps, false);
            if (loop.edges.empty()) { break; }
            if (!loop.closed) {
                // The cascade and clip share the seed's boundary conditions; an
                // open result is precisely the surviving run between two open-air
                // clip points. Close it explicitly across that true-seed air edge.
                // This occurs before link arclengths are established.
                const Pos end = Geo::Chain::eEnd(*loop.edges.back());
                const Pos start = Geo::Chain::eStart(*loop.edges.front());
                if ((start - end).pythag() > eps) {
                    auto closure = std::make_unique<Segment2>(end, start);
                    closure->openAir = true;
                    loop.edges.push_back(std::move(closure));
                }
                loop.closed = true;
                Geo::Chain::hardWeld(loop, eps);
            }
            loop.id = ring.id;
            loop.parent = ring.parent;
            loop.link = ring.link;
            loops.push_back(std::move(loop));
        }

        // Reassembly hygiene before nesting is read.  Coincident traces are the
        // seed-boundary copy of the same frontier; same-handed inner loops are
        // stitch artifacts, not legitimate alternating region boundaries.
        auto interiorPoint = [](const Geo::Chain& c) {
            Pos p = Geo::Chain::edgeMidpoint(*c.edges.front());
            const Pos t = Geo::Chain::eStartDir(*c.edges.front());
            const Pos left(-t.y, t.x);
            for (float d : { 1e-3f, 1e-2f, 1e-1f }) {
                const Pos a = p + left * d, b = p - left * d;
                if (c.windingAround(a) != 0) { return a; }
                if (c.windingAround(b) != 0) { return b; }
            }
            return p;
        };
        auto ridesOn = [](const Geo::Chain& a, const Geo::Chain& b) {
            constexpr float tolerance = 3e-3f;
            for (const auto& e : a.edges) {
                for (float u : { 0.0f, 0.5f, 1.0f }) {
                    const Pos p = Geo::Chain::edgePointAt(*e, u);
                    float distance = 1e30f;
                    for (const auto& be : b.edges) distance = std::min(distance, be->distanceTo(p));
                    if (distance > tolerance) { return false; }
                }
            }
            return true;
        };
        std::vector<char> drop(loops.size(), 0);
        for (size_t i = 0; i < loops.size(); i++) {
            if (loops[i].turningSign() == 0) { drop[i] = 1; continue; }
            for (size_t j = 0; j < loops.size() && !drop[i]; j++) {
                if (i != j && !drop[j] && ridesOn(loops[i], loops[j])) { drop[i] = 1; }
            }
        }
        for (size_t i = 0; i < loops.size(); i++) {
            if (drop[i]) { continue; }
            const Pos p = interiorPoint(loops[i]);
            for (size_t j = 0; j < loops.size(); j++) {
                if (i == j || drop[j] || loops[i].turningSign() != loops[j].turningSign()) { continue; }
                if (loops[j].windingAround(p) != 0
                    && loops[i].windingAround(interiorPoint(loops[j])) == 0) {
                    drop[i] = 1;
                    break;
                }
            }
        }
        size_t write = 0;
        for (size_t i = 0; i < loops.size(); i++) if (!drop[i]) loops[write++] = std::move(loops[i]);
        loops.resize(write);

        // Direction is a cutting contract.  Reassembly may reverse seed pieces,
        // so restore the source sense by nesting parity after all split children
        // are known.  This is also what makes a clipped annulus emit its hole in
        // the opposite direction.
        const int sourceTurn = ring.turningSign();
        for (size_t i = 0; i < loops.size(); i++) {
            Pos probe = Geo::Chain::edgeMidpoint(*loops[i].edges.front());
            const Pos t = Geo::Chain::eStartDir(*loops[i].edges.front());
            const Pos left(-t.y, t.x);
            for (float d : { 1e-3f, 1e-2f, 1e-1f }) {
                const Pos a = probe + left * d, b = probe - left * d;
                if (loops[i].windingAround(a) != 0) { probe = a; break; }
                if (loops[i].windingAround(b) != 0) { probe = b; break; }
            }
            size_t depth = 0;
            for (size_t j = 0; j < loops.size(); j++) {
                if (i != j && loops[j].windingAround(probe) != 0) { depth++; }
            }
            const int wanted = (depth & 1) ? -sourceTurn : sourceTurn;
            if (sourceTurn != 0 && loops[i].turningSign() != wanted) {
                loops[i] = loops[i].reversed();
            }
        }
        pureBoundary = ringOutside && !keptRingPiece && !loops.empty();
        return loops;
    }

    struct RingStation {
        float s = 0.0f;
        Pos p{};
        float cost = 0.0f;
    };

    inline float ringLength(const Geo::Chain& c) {
        float n = 0.0f;
        for (const auto& e : c.edges) { n += Geo::Chain::edgeLength(*e); }
        return n;
    }

    inline std::vector<RingStation> sampleRing(const Geo::Chain& c, float cadence) {
        std::vector<RingStation> out;
        float acc = 0.0f;
        for (size_t k = 0; k < c.edges.size(); k++) {
            const float len = Geo::Chain::edgeLength(*c.edges[k]);
            const int count = std::max(1, static_cast<int>(std::ceil(len / std::max(cadence, 1e-2f))));
            for (int i = 0; i < count; i++) {
                const float u = (static_cast<float>(i) + 0.5f) / static_cast<float>(count);
                out.push_back({ acc + u * len, Geo::Chain::edgePointAt(*c.edges[k], u), 0.0f });
            }
            acc += len;
        }
        return out;
    }

    inline float cyclicDistance(float a, float b, float length) {
        const float d = std::fabs(a - b);
        return std::min(d, std::max(length - d, 0.0f));
    }

    inline float linkSeparation(float ringLength, float toolRadius) {
        // Two radii is the preferred corridor width. A short ring is valid
        // geometry, so saturate at the only possible upper bound instead of
        // treating an unsatisfiable preference as a topology failure.
        return std::min(2.0f * toolRadius, 0.5f * ringLength);
    }

    inline Geo::Chain ringPrefix(const Geo::Chain& c, float length) {
        Geo::Chain out;
        out.id = newId();
        float acc = 0.0f;
        for (const auto& e : c.edges) {
            const float n = Geo::Chain::edgeLength(*e);
            if (acc + n <= length + 1e-4f) {
                out.edges.push_back(e->clone());
                acc += n;
                if (acc >= length - 1e-4f) { break; }
            }
            else if (n > 1e-6f) {
                std::unique_ptr<Stoicheion> a, b;
                Geo::Chain::splitEdge(*e, Geo::Chain::edgePointAt(*e, (length - acc) / n), a, b);
                if (a) { out.edges.push_back(std::move(a)); }
                break;
            }
        }
        Geo::Chain::hardWeld(out, 1e-3f);
        return out;
    }

    struct NormalHit {
        bool found = false;
        bool exact = false;
        Pos point{};
        size_t edge = 0;
        float distance = 0.0f;
    };

    // Cast from a ring link point along its orientation-aware geometric normal.
    // CCW loops have their bounded interior on the left, so their outward normal
    // is right; CW loops reverse that fact. `outward` selects whether the target
    // is the enclosing parent or an enclosed child. The opposite normal is only a
    // fallback for degenerate/reoriented inputs whose turning sign is unavailable.
    inline NormalHit normalProjection(const Geo::Chain& source, Pos origin, Pos tangent,
                                      const Geo::Chain& target, bool outward,
                                      float eps = 1e-3f) {
        NormalHit miss;
        if (source.edges.empty() || target.edges.empty()) { return miss; }
        const float tangentLength = tangent.pythag();
        if (tangentLength <= 1e-8f) {
            Pos closest;
            const size_t edge = target.nearestPoint(origin, closest);
            return { true, false, closest, edge, (closest - origin).pythag() };
        }
        tangent = tangent / tangentLength;

        int turn = source.turningSign();
        if (turn == 0) { turn = 1; }
        // right normal for CCW, left normal for CW
        Pos normal(tangent.y * static_cast<float>(turn),
                   -tangent.x * static_cast<float>(turn));
        if (!outward) { normal = normal * -1.0f; }

        float reach = 1.0f;
        for (const auto& e : target.edges) {
            reach = std::max(reach, (Geo::Chain::eStart(*e) - origin).pythag()
                + Geo::Chain::edgeLength(*e));
            reach = std::max(reach, (Geo::Chain::eEnd(*e) - origin).pythag()
                + Geo::Chain::edgeLength(*e));
        }
        reach *= 2.0f;

        auto cast = [&](Pos rayOrigin, Pos direction) {
            NormalHit best;
            best.distance = 1e30f;
            Segment2 ray(rayOrigin, rayOrigin + direction * reach);
            for (size_t k = 0; k < target.edges.size(); k++) {
                std::vector<Pos> hits;
                Geo::Chain::edgeCross(ray, *target.edges[k], hits);
                for (Pos p : hits) {
                    const Pos delta = p - rayOrigin;
                    const float along = delta.x * direction.x + delta.y * direction.y;
                    const float off = std::fabs(delta.x * direction.y - delta.y * direction.x);
                    if (along <= eps || off > 4.0f * eps || along >= best.distance) { continue; }
                    best = { true, true, p, k, (p - origin).pythag() };
                }
            }
            return best;
        };

        auto rotate = [](Pos v, float angle) {
            const float c = std::cos(angle), s = std::sin(angle);
            return Pos(v.x * c - v.y * s, v.x * s + v.y * c);
        };

        // Exact normal first. Vertex landings and near-tangencies can disappear
        // under crossing tolerances, so retry parallel rays from nearby positions,
        // then minimally perturbed directions, before trying the opposite normal.
        for (float sign : { 1.0f, -1.0f }) {
            const Pos base = normal * sign;
            NormalHit hit = cast(origin, base);
            if (hit.found) { return hit; }
            for (float shift : { 4.0f * eps, -4.0f * eps, 16.0f * eps, -16.0f * eps }) {
                hit = cast(origin + tangent * shift, base);
                if (hit.found) { return hit; }
            }
            for (float angle : { 1e-4f, -1e-4f, 1e-3f, -1e-3f, 1e-2f, -1e-2f }) {
                hit = cast(origin, rotate(base, angle));
                if (hit.found) { return hit; }
            }
        }

        // A closed lineage target always has a usable closest point even if its
        // analytic crossing predicates are numerically inconclusive. Preserve the
        // stay-down transition and let this be the path-wise compromise.
        Pos closest;
        const size_t edge = target.nearestPoint(origin, closest);
        return { true, false, closest, edge, (closest - origin).pythag() };
    }

    // Geometric correspondence between clipped generations. Parallel offsets
    // remain nearest along a local normal even when clipping destroys containment.
    struct NormalAdjacency {
        bool rayHit = false;
        float distance = 1e30f;
    };
    inline NormalAdjacency normalAdjacency(const Geo::Chain& source, const Geo::Chain& target,
                                           float eps = 1e-3f) {
        NormalAdjacency best;
        if (source.edges.empty() || target.edges.empty()) { return best; }
        for (const auto& edge : source.edges) {
            for (float u : { 0.0f, 0.5f, 1.0f }) {
                const Pos p = Geo::Chain::edgePointAt(*edge, u);
                Pos tangent;
                if (u <= 0.0f) { tangent = Geo::Chain::eStartDir(*edge); }
                else if (u >= 1.0f) { tangent = Geo::Chain::eEndDir(*edge); }
                else {
                    tangent = (Geo::Chain::edgePointAt(*edge, 0.51f)
                             - Geo::Chain::edgePointAt(*edge, 0.49f)).normalized();
                }
                const NormalHit hit = normalProjection(source, p, tangent, target, true, eps);
                if (!hit.found) { continue; }
                if ((hit.exact && !best.rayHit)
                    || (hit.exact == best.rayHit && hit.distance < best.distance)) {
                    best.rayHit = hit.exact;
                    best.distance = hit.distance;
                }
            }
        }
        return best;
    }

    // Complete the post-cascade geometry.  Clipping can split one ancestral ring
    // into several loops and invalidate id genealogy, so hierarchy is reconstructed
    // exclusively by winding containment over the FINAL loops.  The function then
    // chooses both link points and returns a leaf-first post-order traversal.
    inline RingPlan planRings(std::vector<Geo::Chain>& rings,
                              std::vector<char>& wall,
                              std::vector<char>& flip,
                              const std::vector<Geo::Chain>& clip,
                              const std::vector<Geo::Chain>& fill,
                              const std::vector<Geo::Chain>& legal,
                              const std::vector<Geo::Chain>& avoid,
                              float toolRadius, float eps = 1e-3f) {
        // Clip only cascade rings.  A pure clip-boundary duplicate is emitted once;
        // an authoritative wall pre-claims that slot.
        if (!clip.empty()) {
            std::vector<Geo::Chain> next;
            std::vector<char> nextWall, nextFlip;
            bool boundaryEmitted = std::any_of(wall.begin(), wall.end(), [](char v) { return v != 0; });
            for (size_t i = 0; i < rings.size(); i++) {
                const char w = i < wall.size() ? wall[i] : 0;
                const char f = i < flip.size() ? flip[i] : 0;
                if (w || !rings[i].closed || rings[i].edges.empty()) {
                    next.push_back(std::move(rings[i])); nextWall.push_back(w); nextFlip.push_back(f);
                    continue;
                }
                bool pure = false;
                std::vector<Geo::Chain> split = clipRing(rings[i], clip, eps, pure);
                if (pure && boundaryEmitted) { continue; }
                if (pure) { boundaryEmitted = true; }
                for (Geo::Chain& c : split) {
                    // openAir is a disk predicate, not merely a centre predicate.
                    for (auto& e : c.edges) {
                        if (!e->openAir) { continue; }
                        for (Pos p : { Geo::Chain::eStart(*e), Geo::Chain::edgeMidpoint(*e), Geo::Chain::eEnd(*e) }) {
                            bool clear = !inside(fill, p);
                            for (const Geo::Chain& b : fill) {
                                for (const auto& be : b.edges) {
                                    if (be->distanceTo(p) < 2.0f * toolRadius) { clear = false; }
                                }
                            }
                            if (!clear) { e->openAir = false; break; }
                        }
                    }
                    next.push_back(std::move(c)); nextWall.push_back(w); nextFlip.push_back(f);
                }
            }
            rings = std::move(next); wall = std::move(nextWall); flip = std::move(nextFlip);
        }

        for (Geo::Chain& r : rings) {
            if (r.edges.empty()) { continue; }
            if (!r.closed) {
                const Pos a = Geo::Chain::eEnd(*r.edges.back());
                const Pos b = Geo::Chain::eStart(*r.edges.front());
                if ((b - a).pythag() > 1e-9f) { r.edges.push_back(std::make_unique<Segment2>(a, b)); }
                r.closed = true;
            }
            Geo::Chain::hardWeld(r, eps);
        }

        const size_t n = rings.size();
        RingPlan plan;
        if (n == 0) { return plan; }

        // Reconstruct generation branches BEFORE link seating. One ancestral ring
        // can yield several clipped pieces carrying the same id; choose the piece
        // belonging to each child by containment first, normal adjacency second.
        std::vector<int> generationParent(n, -1);
        for (size_t child = 0; child < n; child++) {
            if (rings[child].parent == 0 || rings[child].edges.empty()) { continue; }
            int best = -1;
            int bestGeometry = -1;
            float bestDistance = 1e30f;
            const Pos probe = Geo::Chain::edgeMidpoint(*rings[child].edges.front());
            for (size_t candidate = 0; candidate < n; candidate++) {
                if (candidate == child || rings[candidate].id != rings[child].parent
                    || rings[candidate].edges.empty()) { continue; }
                const bool contains = rings[candidate].windingAround(probe) != 0;
                const NormalAdjacency ray = normalAdjacency(rings[child], rings[candidate], eps);
                // IDs constrain only the adjacent generation. Geometry chooses
                // the branch: exact normal intersection is primary, containment
                // corroborates it, and distance resolves equal evidence.
                const int geometry = (ray.rayHit ? 2 : 0) + (contains ? 1 : 0);
                if (best < 0 || geometry > bestGeometry
                    || (geometry == bestGeometry && ray.distance < bestDistance)) {
                    best = static_cast<int>(candidate);
                    bestGeometry = geometry;
                    bestDistance = ray.distance;
                }
            }
            generationParent[child] = best;
        }

        // If ancestry metadata was unavailable, strict containment still gives a
        // topological parent. Choose the smallest enclosing loop (the immediate
        // container), not whichever loop happened to be produced first.
        for (size_t child = 0; child < n; child++) {
            if (generationParent[child] >= 0 || rings[child].edges.empty()) { continue; }
            const Pos probe = Geo::Chain::edgeMidpoint(*rings[child].edges.front());
            int best = -1;
            float bestArea = 1e30f;
            for (size_t candidate = 0; candidate < n; candidate++) {
                if (candidate == child || rings[candidate].edges.empty()
                    || rings[candidate].windingAround(probe) == 0) { continue; }
                const float area = std::fabs(rings[candidate].signedArea());
                if (area < bestArea) { best = static_cast<int>(candidate); bestArea = area; }
            }
            generationParent[child] = best;
        }

        // Preserve requested outside-in/inside-out sense, but make each branch
        // contiguous. Existing relative parent/child positions tell us whether the
        // caller requested parent-first or child-first traversal.
        int directionVote = 0;
        for (size_t i = 0; i < n; i++) if (generationParent[i] >= 0) {
            directionVote += generationParent[i] < static_cast<int>(i) ? 1 : -1;
        }
        const bool parentFirst = directionVote >= 0;
        std::vector<std::vector<size_t>> generationChildren(n);
        std::vector<size_t> roots;
        for (size_t i = 0; i < n; i++) {
            if (generationParent[i] >= 0) {
                generationChildren[static_cast<size_t>(generationParent[i])].push_back(i);
            }
            else { roots.push_back(i); }
        }
        std::vector<size_t> order;
        auto visit = [&](auto&& self, size_t node) -> void {
            if (parentFirst) { order.push_back(node); }
            for (size_t child : generationChildren[node]) { self(self, child); }
            if (!parentFirst) { order.push_back(node); }
        };
        for (size_t root : roots) { visit(visit, root); }
        if (order.size() == n) {
            std::vector<size_t> newIndex(n);
            for (size_t i = 0; i < n; i++) { newIndex[order[i]] = i; }
            std::vector<Geo::Chain> sortedRings;
            std::vector<char> sortedWall, sortedFlip;
            sortedRings.reserve(n); sortedWall.reserve(n); sortedFlip.reserve(n);
            for (size_t old : order) {
                sortedRings.push_back(std::move(rings[old]));
                sortedWall.push_back(old < wall.size() ? wall[old] : 0);
                sortedFlip.push_back(old < flip.size() ? flip[old] : 0);
            }
            std::vector<int> sortedParent(n, -1);
            for (size_t old = 0; old < n; old++) if (generationParent[old] >= 0) {
                sortedParent[newIndex[old]] = static_cast<int>(newIndex[static_cast<size_t>(generationParent[old])]);
            }
            rings = std::move(sortedRings); wall = std::move(sortedWall); flip = std::move(sortedFlip);
            generationParent = std::move(sortedParent);
        }

        plan.parent = generationParent;
        plan.linkCenterDistance.resize(n); plan.linkStartDistance.resize(n); plan.linkEndDistance.resize(n);
        plan.linkCenter.resize(n);
        plan.linkStart.resize(n); plan.linkEnd.resize(n);
        plan.cycle.resize(n); plan.lap.resize(n); plan.overlap.resize(n);

        std::vector<size_t> depth(n, 0);
        for (size_t i = 0; i < n; i++) {
            if (rings[i].edges.empty()) { continue; }
            const Pos p = Geo::Chain::edgeMidpoint(*rings[i].edges.front());
            for (size_t j = 0; j < n; j++) {
                if (i != j && !rings[j].edges.empty() && rings[j].windingAround(p) != 0) { depth[i]++; }
            }
        }
        for (size_t i = 0; i < n; i++) {
            if (rings[i].edges.empty() || plan.parent[i] >= 0) { continue; }
            const Pos p = Geo::Chain::edgeMidpoint(*rings[i].edges.front());
            int enclosing = -1;
            for (size_t j = 0; j < n; j++) {
                if (i == j || rings[j].edges.empty() || rings[j].windingAround(p) == 0) { continue; }
                if (enclosing < 0 || depth[j] > depth[static_cast<size_t>(enclosing)]) {
                    enclosing = static_cast<int>(j);
                }
            }
            plan.parent[i] = enclosing;
        }

        auto clearance = [&](Pos p) {
            float d = 1e30f;
            for (const auto* set : { &fill, &avoid }) {
                for (const Geo::Chain& c : *set) for (const auto& e : c.edges) d = std::min(d, e->distanceTo(p));
            }
            return d;
        };
        auto clearedAir = [&](Pos p) { return !inside(fill, p) && inside(legal, p); };

        std::vector<std::vector<RingStation>> stations(n);
        std::vector<float> lengths(n, 0.0f), seatArc(n, 0.0f);
        std::vector<size_t> byDepth;
        for (size_t i = 0; i < n; i++) {
            if (rings[i].edges.empty()) { continue; }
            lengths[i] = ringLength(rings[i]);
            stations[i] = sampleRing(rings[i], 0.5f * toolRadius);
            if (!stations[i].empty()) { byDepth.push_back(i); }
        }
        std::sort(byDepth.begin(), byDepth.end(), [&](size_t a, size_t b) { return depth[a] < depth[b]; });

        for (size_t i : byDepth) {
            auto& s = stations[i];
            const int parent = plan.parent[i];
            if (parent >= 0 && !stations[static_cast<size_t>(parent)].empty()) {
                const size_t p = static_cast<size_t>(parent);
                for (RingStation& q : s) {
                    q.cost = 1e30f;
                    for (const RingStation& a : stations[p]) {
                        q.cost = std::min(q.cost, a.cost
                            + 0.25f * cyclicDistance(a.s, seatArc[p], lengths[p])
                            + (q.p - a.p).pythag());
                    }
                }
            }
            else {
                // A root seeks the CENTRE of a genuinely long cleared-air run;
                // disconnected runs are never flattened into one sample list.
                // With no qualifying run, use the smoothest full link window.
                const size_t count = s.size();
                const float window = 4.0f * toolRadius;
                std::vector<char> isAir(count, 0);
                bool allAir = true, anyAir = false;
                for (size_t k = 0; k < count; k++) {
                    isAir[k] = clearedAir(s[k].p) ? 1 : 0;
                    allAir = allAir && isAir[k]; anyAir = anyAir || isAir[k];
                }
                std::vector<float> centers;
                if (allAir) {
                    size_t best = 0;
                    for (size_t k = 1; k < count; k++)
                        if (clearance(s[k].p) > clearance(s[best].p)) best = k;
                    centers.push_back(s[best].s);
                }
                else if (anyAir) {
                    size_t material = 0;
                    while (material < count && isAir[material]) material++;
                    size_t runStart = count;
                    for (size_t step = 1; step <= count; step++) {
                        const size_t k = (material + step) % count;
                        if (isAir[k] && runStart == count) { runStart = k; }
                        if ((!isAir[k] || step == count) && runStart != count) {
                            const size_t last = (material + step - 1) % count;
                            float span = s[last].s - s[runStart].s;
                            if (span < 0.0f) { span += lengths[i]; }
                            // Samples sit at interval midpoints; include one cadence
                            // so the test measures the represented run, not only the
                            // distance between its first and last midpoint.
                            span += std::min(0.5f * toolRadius, lengths[i]);
                            if (span >= window) {
                                float center = s[runStart].s + 0.5f * span;
                                if (center >= lengths[i]) { center -= lengths[i]; }
                                centers.push_back(center);
                            }
                            runStart = count;
                        }
                    }
                }
                if (centers.empty()) {
                    std::vector<float> curvature(count, 0.0f), rough(count, 0.0f);
                    for (size_t k = 0; k < count; k++) {
                        const Pos a = s[k].p - s[(k + count - 1) % count].p;
                        const Pos b = s[(k + 1) % count].p - s[k].p;
                        const float la = a.pythag(), lb = b.pythag();
                        if (la > 1e-6f && lb > 1e-6f) {
                            curvature[k] = std::atan2(std::fabs(a.x * b.y - a.y * b.x),
                                a.x * b.x + a.y * b.y) / (0.5f * (la + lb));
                        }
                    }
                    for (size_t k = 0; k < count; k++) {
                        const size_t before = (k + count - 1) % count;
                        const float ds = std::max(cyclicDistance(s[k].s, s[before].s, lengths[i]), 1e-6f);
                        rough[k] = curvature[k] * toolRadius
                            + std::fabs(curvature[k] - curvature[before]) / ds * toolRadius * toolRadius;
                    }
                    size_t best = 0; float bestSum = 1e30f;
                    for (size_t k = 0; k < count; k++) {
                        float sum = 0.0f;
                        for (size_t j = 0; j < count; j++)
                            if (cyclicDistance(s[k].s, s[j].s, lengths[i]) <= 0.5f * window) sum += rough[j];
                        if (sum < bestSum) { best = k; bestSum = sum; }
                    }
                    centers.push_back(s[best].s);
                }
                for (RingStation& q : s) {
                    q.cost = 1e30f;
                    for (float center : centers)
                        q.cost = std::min(q.cost, cyclicDistance(q.s, center, lengths[i]));
                }
            }

            size_t best = 0;
            for (size_t k = 1; k < s.size(); k++) {
                if (s[k].cost + 1e-4f < s[best].cost
                    || (std::fabs(s[k].cost - s[best].cost) <= 1e-4f
                        && clearance(s[k].p) > clearance(s[best].p))) { best = k; }
            }
            seatArc[i] = s[best].s;
            const float separation = linkSeparation(lengths[i], toolRadius);
            const float halfGap = 0.5f * separation;
            auto wrap = [&](float distance) {
                distance = std::fmod(distance, lengths[i]);
                return distance < 0.0f ? distance + lengths[i] : distance;
            };
            plan.linkCenterDistance[i] = wrap(s[best].s);
            plan.linkStartDistance[i] = wrap(s[best].s - halfGap);
            plan.linkEndDistance[i] = wrap(s[best].s + halfGap);
            plan.linkCenter[i] = rings[i].pointAtDistance(plan.linkCenterDistance[i]);
            plan.linkStart[i] = rings[i].pointAtDistance(plan.linkStartDistance[i]);
            plan.linkEnd[i] = rings[i].pointAtDistance(plan.linkEndDistance[i]);
            plan.cycle[i] = rings[i].startedAtDistance(plan.linkEndDistance[i]);
            plan.lap[i] = ringPrefix(plan.cycle[i], lengths[i] - separation);
            plan.overlap[i] = {};
        }

        // Preserve the established march order. The caller has already applied
        // inside-out/outside-in policy (including mixed air-rooted lineages);
        // clipping preserves that relative order. Geometry rebuilds the tree only
        // to decide which consecutive transitions are lineage links.
        std::vector<std::vector<size_t>> children(n);
        for (size_t i = 0; i < n; i++) {
            if (plan.parent[i] >= 0) { children[static_cast<size_t>(plan.parent[i])].push_back(i); }
        }

        // Entrance is execution-order topology, not a fixed leaf property. The
        // first member reached in either direction is full-width unless its parent
        // or all of its children have already been emitted.
        std::vector<char> emitted(n, 0);
        for (size_t i = 0; i < n; i++) {
            const bool parentCut = plan.parent[i] >= 0
                && emitted[static_cast<size_t>(plan.parent[i])] != 0;
            bool childrenCut = !children[i].empty();
            for (size_t child : children[i]) { childrenCut = childrenCut && emitted[child] != 0; }
            rings[i].link = (!parentCut && !childrenCut) ? LinkKind::Entrance : LinkKind::None;
            plan.cycle[i].link = rings[i].link;
            plan.lap[i].link = rings[i].link;
            emitted[i] = 1;
        }


        plan.linkToNext.assign(rings.size(), 0);

        // Resolve stay-down links in actual marching order. Endpoint placement is
        // already complete: the hierarchy-wide field selected every linkCenter,
        // and START/CENTER/END were derived symmetrically on that ring. Do NOT
        // overwrite the next ring from the previous ring's ray. That recurrence
        // changes B's center from A's shifted start, then C's from B's newly
        // shifted start, mechanically rotating the corridor once per generation.
        //
        //     current.linkStart  ------------>  next.linkEnd
        //
        // Both operands are fixed Cartesian points. A normal projection may be a
        // future scoring term for choosing the common corridor, but is not allowed
        // to re-seat either ring after the global seating pass. Branch jumps remain
        // independently seated because they retract.
        for (size_t next = 1; next < rings.size(); next++) {
            const size_t current = next - 1;
            const bool toParent = plan.parent[current] == static_cast<int>(next);
            const bool toChild = plan.parent[next] == static_cast<int>(current);
            if ((!toParent && !toChild) || rings[current].edges.empty() || rings[next].edges.empty()) { continue; }
            plan.linkToNext[current] = 1;
        }

        // Bind metadata to the exact operands Strategy will emit.
        for (size_t i = 0; i < plan.lap.size(); i++) {
            if (plan.lap[i].edges.empty()) { continue; }
            plan.linkEnd[i] = Geo::Chain::eStart(*plan.lap[i].edges.front());
            plan.linkStart[i] = Geo::Chain::eEnd(*plan.lap[i].edges.back());
        }
        return plan;
    }

    // Convenience: one-shot form.
    inline std::vector<Profile> runParallel(std::vector<std::unique_ptr<Stoicheion>> fill,
                                            std::vector<std::unique_ptr<Stoicheion>> avoid,
                                            const Parallel::Params& params) {
        Parallel p(params);
        p.ingest(std::move(fill), std::move(avoid));
        p.run();
        return std::move(p.profiles);
    }
}
