module;

#include <vector>
#include <memory>
#include <string>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <algorithm>

export module Geo.Frontier;

export import Geo.Chain;

// =====================================================================
// Geo::Frontier -- the frontier-extraction method, reconstructed.
//
// A parallel, spec-driven reimplementation of the offset / crossing /
// merge machinery in Geo::Chain, written against the method documents in
// external/Geometry/FrontierExtraction/Method (01..06). The original was
// written in discovery mode and hard-codes one question -- "the internal
// frontier of an inward offset" -- as constants (minimum crossing level,
// CCW-only inversion trial, auto-detected merge sign). Here the question
// is a parameter:
//
//   * frontier selection and the inversion guard both key on the SENSE
//     sigma = sign(amount * h) -- the offset sign RELATIVE to the chain's
//     handedness -- not on either alone (Method/03);
//   * the merge is an explicit winding-rule difference: material is
//     { winding >= 1 } with doctrine orientation asserted upstream, and
//     minuend/subtrahend roles decide which clusters may yield material
//     at all (Method/04, 06);
//   * mitosis (one chain splitting against itself) and merging (resolved
//     chains meeting) are kept distinct: offsetValid vs merge.
//
// Stoicheia are IMPORTED, not reimplemented -- and so are the edge-level
// numerical kernels (offsetLeft, edgeCross, splitEdgeAtPoints, endpoint /
// tangent probes), reused via Geo::Chain's statics under the alias Kern.
// That is deliberate: those are pure stoicheion-level functions carrying
// the hard-won robustness lore (double-precision crossings, tangency
// fusion, absolute-distance interiority, transversal-kiss rejection), and
// sharing them means a differential test against the old kernel compares
// the METHOD, not two transcriptions of the same arithmetic.
//
// SELECTION RULE (corrected 2026-07-11; see offsetValid): the frontier
// extreme keys on sign(amount) -- the side the curve was pushed -- NOT on
// sigma. A mirror argument forbids sigma-keyed selection (reflection
// preserves sigma but negates the level field); production anchors pin
// both halves (legacy engine: min at amount > 0; the crescent seed: max
// at amount < 0). Sigma gates only the inversion guard.
//
// KNOWN LIMITATION: a growing chain's legitimate secondary births (a
// pinched-off void between two arms) are dropped with the non-extreme
// levels. Affects efficiency, not safety.
// =====================================================================
export namespace Geo::Frontier {

    // The proven edge-level kernel host (Geo::Chain's statics).
    using Kern = Geo::Chain;

    // Which side of the difference a chain argues for. Minuend chains bound
    // the region we are keeping (fill); subtrahend chains bound geometry we
    // are removing from it (avoid). Merge drops any cluster with no minuend:
    // subtrahend geometry alone can never yield material (Method/04).
    enum class Role : std::uint8_t { Minuend, Subtrahend };

    // =================================================================
    // Chain -- an oriented run of stoicheia; travel is meaning, material
    // is on the LEFT of travel (Method/02, 05).
    // =================================================================
    struct Chain {

        std::vector<std::unique_ptr<Stoicheion>> edges;   // oriented: travel is start -> end
        bool closed = false;

        Id id = 0;
        Id parent = 0;
        Role role = Role::Minuend;

        // An offset arc whose radius lands below this is a collapsed feature (a
        // fillet at the offset distance, within numerical slop): culled outright
        static constexpr float RadiusCull = 1e-2f;

        Chain() = default;
        Chain(Chain&&) = default;
        Chain& operator=(Chain&&) = default;

        Chain clone() const {
            Chain c;
            c.closed = closed; c.id = id; c.parent = parent; c.role = role;
            for (const auto& e : edges) { c.edges.push_back(e->clone()); }
            return c;
        }

        // The opposite journey: edge order flipped, every edge reversed. For
        // CONDITIONING inputs whose winding is arbitrary -- never for
        // reorienting the method's own outputs (Method/05).
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
        //--------------------------------------------------

        // Exact signed area (Green's theorem, closed-form per edge). This is
        // the judge for ORIENTING inputs: exact, immune to the tangency-cusp
        // turning flips that bit the old kernel three times.
        float signedArea() const {
            float twice = 0.0f;
            for (const auto& e : edges) { twice += e->signedAreaTerm(); }
            return 0.5f * twice;
        }

        // Handedness as a REGION fact: +1 CCW (bounds material inside),
        // -1 CW (an island), 0 degenerate. Judged by area, not turning --
        // the sense sigma is about what the chain's region does under the
        // offset, an area question.
        int handedness() const {
            float a = signedArea();
            return (a > 1e-4f) ? 1 : (a < -1e-4f ? -1 : 0);
        }

        // The accumulated tangent rotation: 2*pi * rotation index for a
        // closed path. This is the judge for OFFSET OUTPUT children (the
        // inversion guard): burr lobes can outweigh a thin loop's area and
        // flip its sign, but they can never flip total turning, because
        // offsetting preserves the tangent structure of its source.
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

        // Winding number around a probe, by ACCUMULATED SUBTENDED ANGLE --
        // deliberately not ray casting: the angle sum is continuous in every
        // vertex position and can only destabilise when the probe sits ON the
        // curve (Method/04 scar tissue).
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
                        // Inside the circle the bearing advances monotonically
                        // CCW, so the true advance lies in (0, 2pi): unique branch.
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

        // The point of greatest x anywhere on the loop: provably unenclosed
        // within its own cluster, and the containment probe (on this loop,
        // far from another unless they genuinely touch).
        static Pos maxXPoint(const Chain& L) {
            Pos best(-1e30f, 0.0f);
            for (const auto& e : L.edges) {
                Pos p = Kern::edgeMaxXPoint(*e);
                if (p.x > best.x) { best = p; }
            }
            return best;
        }

        // Beat one: the blind offset (Method/02)
        //--------------------------------------------------
        // Every edge slides to the LEFT of its own travel by `amount`; every
        // corner is bridged by the join arc known by construction (centred on
        // the original corner, radius |amount|). No winding test, no preferred
        // side -- inset vs outset EMERGES from each chain's own chirality. The
        // result may self-intersect; those crossings are beat two's input.
        Chain offsetRaw(float amount) const {
            Chain result;
            if (edges.empty()) { return result; }
            result.closed = closed;

            // A lone prime closed curve is its own loop: no corners, just the
            // blind concentric left-offset.
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

            // The push: offset each edge alone; an edge that collapses under
            // the offset generates nothing (its neighbours join across the hole).
            std::vector<std::unique_ptr<Stoicheion>> off;
            std::vector<size_t> src;                   // off[i] came from edges[src[i]]
            off.reserve(edges.size());
            src.reserve(edges.size());
            for (size_t i = 0; i < edges.size(); i++) {
                std::unique_ptr<Stoicheion> o = Kern::offsetLeft(*edges[i], amount);
                if (Kern::pruneSmallFeatures && Kern::edgeLength(*o) < Kern::PruneLength) { continue; }

                // An arc driven to (near) zero radius is a collapsed feature, not a
                // boundary. The length prune alone misses it: through-zero inversion
                // keeps the through-point's bearing, so the junk re-emerges with the
                // COMPLEMENT sweep -- long enough to survive, chirality flipped. At
                // fillet-radius == tool-radius (with section slop) this is the
                // EXPECTED case, not the corner: cull by RADIUS, outright, and let
                // the corner construction join the neighbours directly.
                if (Kern::isArc(*o) && Kern::asArc(*o)->radius() < RadiusCull) { continue; }

                off.push_back(std::move(o));
                src.push_back(i);
            }

            const size_t n = off.size();
            if (n == 0) { return result; }
            const size_t corners = closed ? n : (n > 0 ? n - 1 : 0);

            // The uniform corner construction. Near-tangent seams (the STEP
            // signature) are FUSED, not bridged -- a sliver join loop there is
            // numerical poison (Method/02 scar tissue).
            std::vector<std::unique_ptr<Stoicheion>> joinArc(corners);
            for (size_t k = 0; k < corners; k++) {
                Pos corner = Kern::eEnd(*edges[src[k]]);              // the original corner
                Pos endA = Kern::eEnd(*off[k]);
                Pos startB = Kern::eStart(*off[(k + 1) % n]);
                if ((endA - startB).pythag() < 1e-4f) { continue; }   // tangent-continuous: already joined

                // A hole in the source (an edge pruned or radius-culled between
                // these two) voids the corner-arc construction -- it presumes true
                // adjacency, and across a hole the arc's far end misses startB by
                // up to the culled radius, shredding the chain at kernel eps.
                // Bridge the hole EXACTLY: a straight segment, endpoints verbatim.
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

                // All angle math in DOUBLE: a nearly-complete around-the-bend
                // arc subtracts two nearly-equal angles.
                dPos C(corner);
                dPos radial = dPos(endA) - C;
                double r = radial.len();
                if (r < 1e-12) { continue; }
                int s = (dPos(Kern::eEndDir(*off[k])).dot(radial.perp()) >= 0.0) ? 1 : -1;
                double angA = radial.angle(), angB = (dPos(startB) - C).angle();
                double span = (s > 0) ? dWrap(angB - angA) : dWrap(angA - angB);
                double midAng = angA + static_cast<double>(s) * span * 0.5;
                joinArc[k] = std::make_unique<Arc2>(corner, endA, startB,
                                                    (C + dPos::fromAngle(midAng) * r).f());
            }

            for (size_t i = 0; i < n; i++) {
                result.edges.push_back(std::move(off[i]));
                if (i < corners && joinArc[i]) { result.edges.push_back(std::move(joinArc[i])); }
            }
            return result;
        }

        // Beat two: the crossing-number field (Method/03)
        //--------------------------------------------------
        // Split the raw offset at every self-crossing and walk it in travel
        // order, accumulating a level: at each crossing, the OTHER strand
        // passing right-to-left counts +1, left-to-right -1. Only DIFFERENCES
        // are meaningful (the walk's start is arbitrary), so levels are
        // normalised max -> 0. Which extreme is the true frontier is beat
        // three's decision, NOT the field's.
        struct Field {
            std::vector<std::unique_ptr<Stoicheion>> edges;   // travel order, direction preserved
            std::vector<int> level;                           // accumulated level per sub-edge
        };

        Field field(float eps = 1e-3f) const {
            const size_t n = edges.size();

            // 1. every self-crossing, remembering which two edges meet there.
            //    (edgeCross is transversal-only: a tangential kiss is rejected
            //    at detection -- its +-1 would be float noise in the walk.)
            struct Rec { Pos p; size_t i, j; };
            std::vector<Rec> recs;
            std::vector<std::vector<Pos>> pts(n);
            for (size_t i = 0; i < n; i++) {
                for (size_t j = i + 1; j < n; j++) {
                    std::vector<Pos> c; Kern::edgeCross(*edges[i], *edges[j], c);
                    for (const Pos& p : c) { recs.push_back({ p, i, j }); pts[i].push_back(p); pts[j].push_back(p); }
                }
            }

            // 2. split every edge at its cuts, keeping chain order; each
            //    sub-edge remembers its original edge.
            Field run;
            std::vector<size_t> from;
            for (size_t i = 0; i < n; i++) {
                Kern::splitEdgeAtPoints(*edges[i], pts[i], run.edges);
                while (from.size() < run.edges.size()) { from.push_back(i); }
            }

            // 3. walk strictly in order. A sub-edge wears the level in effect
            //    as it is entered; the level updates at its end vertex if that
            //    vertex is a crossing (the OTHER strand's travel direction there).
            run.level.resize(run.edges.size(), 0);
            int count = 0;
            for (size_t k = 0; k < run.edges.size(); k++) {
                run.level[k] = count;
                Pos v = Kern::eEnd(*run.edges[k]);
                for (const Rec& r : recs) {
                    if ((r.p - v).pythag() >= eps) { continue; }
                    size_t other = (from[k] == r.i) ? r.j
                                 : (from[k] == r.j) ? r.i : static_cast<size_t>(-1);
                    if (other == static_cast<size_t>(-1)) { continue; }
                    Pos tOur = Kern::eEndDir(*run.edges[k]);
                    Pos tOther = Kern::travelDirAt(*edges[other], r.p);
                    count += (tOur.cross(tOther) > 0.0f) ? 1 : -1;
                }
            }

            // Normalise relative to the max: start-invariant with no external
            // measurement. Levels read 0 (max) downward.
            if (!run.level.empty()) {
                int maxN = run.level.front();
                for (int v : run.level) { maxN = std::max(maxN, v); }
                for (int& v : run.level) { v -= maxN; }
            }
            return run;
        }

        // A field fragment: a well-formed sub-chain plus the level its pieces
        // share. (Holds a Chain by value: defined out of line, below.)
        struct Fragment;
        std::vector<Fragment> fragments(float eps = 1e-3f) const;

        // Beat three: sigma-keyed extraction (Method/03)
        //--------------------------------------------------
        std::vector<Chain> offsetValid(float amount, float eps = 1e-3f) const;

        // Merging (Method/04, 06)
        //--------------------------------------------------
        static std::vector<Chain> merge(std::vector<Chain> allChains, float eps = 1e-3f);

        // Building
        //--------------------------------------------------

        // Grow one chain out of the pool, connecting end-to-start until it
        // closes or runs out. `allowReverse` permits flipping an edge to make
        // a connection -- legitimate when assembling unordered input strokes,
        // FORBIDDEN when reassembling oriented frontier pieces (a silent
        // reversal there fabricates wrong-handed loops; the piece set must
        // fail visibly instead).
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
            return c;
        }

        // Assemble chains from raw stoicheia (reversal allowed: the input is
        // unordered and unoriented). Prime closed curves are their own loops;
        // an ellipse enters EXACTLY as zero-distance offset-ellipse spans cut
        // at quadrant boundaries (curvature extrema, where offset cusps form).
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

        // Wrap merge: on a closed walk the first and last fragments can be two
        // halves of the same stretch.
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

    // The method end to end for ONE chain: blind offset, crossing field,
    // side-keyed extreme, oriented stitch, sigma-gated inversion guard.
    //
    // TWO different keys, and they are not the same variable:
    //
    //   * The EXTREME is keyed on sign(amount) -- the SIDE the curve was
    //     pushed, matching the field's own travel-left orientation. Leftward
    //     (amount > 0): the frontier is the MINIMUM level. Rightward
    //     (amount < 0): the MAXIMUM. It cannot key on sigma: mirror the plane
    //     and a proven min-case becomes its opposite-handed twin with every
    //     crossing increment negated -- same sigma, opposite extreme label.
    //     Anchors: the legacy engine (all chains at amount > 0, min, proven
    //     in production, CW-growing islands included) and the crescent seed
    //     (CCW growing at amount < 0, max, verified on a real part).
    //
    //   * The GUARD is keyed on sigma = sign(amount * h) -- a REGION question:
    //     only a shrinking chain (sigma > 0) can collapse through zero, so only
    //     there are handedness-flipped children the fold rather than a boundary.
    //     Growth (sigma < 0) cannot invert; the guard is off. [Growing secondary
    //     births -- a void pinched off between growing arms -- are still dropped
    //     with the non-extreme levels: a known, safe-direction limitation.]
    //
    // h is judged by exact signed AREA (the sense is a region question);
    // the guard compares TURNING (offsetting preserves tangent structure, so
    // turning is the burr-immune judge of an offset child's handedness).
    inline std::vector<Chain> Chain::offsetValid(float amount, float eps) const {
        std::vector<Chain> out;

        // Closed chains only: offsetting is a statement about a bounded
        // region's boundary; an open chain bounds nothing.
        if (!closed || edges.empty()) { return out; }

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

        // The inversion guard, sigma-gated.
        if (sigma > 0.0f) {
            const int srcH = h;
            std::erase_if(out, [srcH](const Chain& m) {
                return m.closed && m.turningSign() == -srcH;
            });
        }
        return out;
    }

    // The merge: resolved loops meet; the material region is selected by an
    // EXPLICIT winding rule and its frontier extracted (Method/04, 06).
    //
    // Doctrine orientation is asserted upstream (material on the left of
    // every edge), so material is simply { winding >= 1 } -- no auto-detected
    // system sign, which is exactly what an enclosing subtrahend perimeter
    // used to hijack. Roles close the remaining hole: a cluster containing no
    // MINUEND loop can bound no material and yields nothing, so disjoint
    // avoid geometry can never survive as a spurious frontier.
    //
    // Envelopment without intersection needs NO special case: clustering is
    // by crossing OR CONTAINMENT, so a grown avoid loop enveloped by the
    // grown fill joins its cluster, and the winding keeps the avoid's loop as
    // the frontier while discarding the fill's (the bore: the seed IS the
    // avoid's standoff ring -- "we moved away from the avoid region").
    inline std::vector<Chain> Chain::merge(std::vector<Chain> allChains, float eps) {

        std::vector<Chain> out;
        std::vector<Chain> loops;
        for (Chain& ch : allChains) {
            if (ch.closed && !ch.edges.empty()) { loops.push_back(std::move(ch)); }
            else if (!ch.edges.empty())         { out.push_back(std::move(ch)); }   // opens pass through
        }

        const size_t n = loops.size();
        if (n == 0) { return out; }

        // Cluster by INTERACTION (cross, or one contains the other). Disjoint
        // profiles are independent worlds.
        auto interacts = [&](size_t a, size_t b) {
            std::vector<Pos> x;
            for (const auto& ea : loops[a].edges) {
                for (const auto& eb : loops[b].edges) {
                    Kern::edgeCross(*ea, *eb, x);
                    if (!x.empty()) { return true; }
                }
            }
            if (loops[b].windingAround(maxXPoint(loops[a])) != 0) { return true; }
            if (loops[a].windingAround(maxXPoint(loops[b])) != 0) { return true; }
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

            std::vector<std::unique_ptr<Stoicheion>> kept;
            for (size_t a : mem) {
                for (const auto& edge : loops[a].edges) {

                    // Cut this edge wherever any other cluster loop crosses it.
                    std::vector<Pos> cuts;
                    for (size_t b : mem) {
                        if (b == a) { continue; }
                        for (const auto& eb : loops[b].edges) { Kern::edgeCross(*edge, *eb, cuts); }
                    }
                    std::vector<std::unique_ptr<Stoicheion>> pieces;
                    Kern::splitEdgeAtPoints(*edge, cuts, pieces);

                    // Frontier test: material on exactly one side. The probe
                    // reach is ABSOLUTE distance, bounded by real clearance to
                    // the other loops -- never a parameter tolerance.
                    for (auto& piece : pieces) {
                        Pos mid = Kern::edgeMidpoint(*piece);
                        Pos t = Kern::travelDirAt(*piece, mid);

                        float clearance = 1e30f;
                        for (size_t b : mem) {
                            for (const auto& le : loops[b].edges) {
                                if (le.get() == edge.get()) { continue; }
                                clearance = std::min(clearance, le->distanceTo(mid));
                            }
                        }
                        float reach = std::max(1e-4f, std::min(1e-2f * Kern::edgeLength(*piece), 0.4f * clearance));

                        auto totalWinding = [&](Pos probe) {
                            int w = 0;
                            for (size_t b : mem) { w += loops[b].windingAround(probe); }
                            return w;
                        };
                        bool leftIn  = totalWinding(mid + Kern::perpCCW(t) * reach) >= 1;
                        bool rightIn = totalWinding(mid - Kern::perpCCW(t) * reach) >= 1;
                        if (leftIn != rightIn) {
                            piece->id = loops[a].id;   // provenance
                            kept.push_back(std::move(piece));
                        }
                    }
                }
            }

            // Stitch the surviving pieces back into loops. Oriented pieces:
            // reversal is forbidden -- if ends don't meet forward, the piece
            // set is wrong and must FAIL VISIBLY (an open remnant), never be
            // silently flipped into a wrong-handed loop.
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
        return out;
    }

    // =================================================================
    // Profile -- one GENERATION: the chains bounding the clearable
    // material at one offset depth. The complete method, one generation
    // forward: stage A per chain (offsetValid), stage B across chains
    // (merge). An exhausted profile returns empty -- extinction, not error.
    // =================================================================
    struct Profile {

        std::vector<Chain> chains;

        // The sanity verdict: a generation whose total signed area goes negative
        // claims an unbounded region (transient garbage -- e.g. a surviving
        // inverted burr). It is recorded and SKIPPED by consumers, never emitted
        // -- and it never halts the ladder: recursion continues and the garbage
        // extinctifies on its own (the legacy engine's concession, kept).
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
            return next;
        }
    };

    // =================================================================
    // Raster -- parallel strokes filling a region (normally the SEED: the
    // legal tool-centre region, so clipped strokes are safe by construction).
    // Each scan line's transversal crossings pair into interior spans by
    // even-odd; strokes come back OPEN, ordered and DIRECTED: `serpentine`
    // reverses alternate lines end-for-end (hatch -- turns at the ends);
    // without it every stroke runs the same way (scallop -- constant
    // engagement sense, the consumer retracts between strokes). Links,
    // retracts and Z are the consumer's grammar, never woven here.
    // =================================================================
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

    // =================================================================
    // Parallel -- the seed construction and the inlay (Method/06).
    //
    // Input: two sets of raw stoicheia. FILL encloses area where cutting
    // is encouraged (the tool must maximally touch it); AVOID encloses
    // area where no part of the tool may enter. Both arrive with NO
    // orientation promise.
    //
    // The strategy:
    //   1. assemble each set into endpoint-coincident chains;
    //   2. orient by NESTING PARITY with opposite bases -- fill outers
    //      CCW (islands CW), avoid outers CW (islands CCW): one shared
    //      nesting tree, handedness following depth parity, the avoid
    //      being the fill's island (Method/05, 06);
    //   3. GROW BOTH by one tool radius -- fill by offsetLeft(-R), avoid
    //      by offsetLeft(+R); opposite signs, both regions expanding;
    //   4. merge: the seed is the frontier of { winding >= 1 } =
    //      (fill (+) R) \ (avoid (+) R) -- the maximum extent of the tool
    //      CENTRE that reaches the fill and never violates the avoid;
    //   5. inlay: recursively inset the seed by the stepover until
    //      extinction. profiles[0] IS the outermost cuttable ring (the
    //      one-radius standoff was built into the seed, so there is no
    //      "offset by one radius first" -- that rule dissolved into the
    //      construction).
    // =================================================================
    struct Parallel {

        struct Params {
            float toolRadius = 1.0f;
            float stepover = 1.0f;        // generation advance as a FRACTION of the tool radius
            int   maxGenerations = 64;    // inlay cap
            float eps = 1e-3f;            // KERNEL tolerance: crossings, stitching, fragment
                                          // coincidence -- pure geometry math, never input slop
            float weldEps = 1e-2f;        // almost-closed loops within this gap are welded shut
                                          // (BRep-section endpoints carry approximation slop)

            // ASSEMBLY tolerances, PER SYSTEM, derived from each sectioned solid's
            // own tolerance (max vertex tolerance; a spec solid ~1e-7, a healed /
            // measured solid ~1e-2). Two endpoints within a solid's own tolerance
            // are the SAME point per its model, so assembling there closes
            // measured-fragmented loops without ever fusing features the solid
            // keeps genuinely apart. The cut and the avoid are DIFFERENT solids
            // (a spec delta against a measured part is routine), so each system
            // assembles at its own figure. 0 = fall back to eps. Assembly
            // tolerance is input conditioning; the kernel eps above never widens.
            float fillEps = 0.0f;
            float avoidEps = 0.0f;

            // The extended (lead) region's two expansions, as FRACTIONS of the
            // tool radius. The survey of the open air lead motion may actually
            // use: fill grown by a full diameter (centres within feasible reach
            // of the cut), avoid grown PAST the ring standoff -- rings graze
            // avoid + R, but approach moves keep an extra half radius off kept
            // walls (over-avoidance is the safe direction; the seed's exact-R
            // bound is untouched). leadReach = 0 disables the construction.
            float leadReach = 2.0f;   // extended fill expansion
            float leadClear = 1.5f;   // extended avoid expansion (never below 1.0)

            // The TRUE-open-air survey's expansions (fractions of R): fill by one
            // full tool diameter, avoid by two thirds of one. An independent,
            // parallel construction (nothing in the seed/inlay pipeline reads it);
            // its product is the openAir member. 0 disables.
            float airReach = 2.0f;
            float airClear = 4.0f / 3.0f;

            // Generic-position symmetry breaker. The fill is grown by R + fillBias
            // while the avoid is grown by EXACTLY R -- so grown boundaries that
            // would otherwise lie coincident for long stretches (a feature wall a
            // whole multiple of R from a kept wall -- tool-radius-dependent, so it
            // strikes "randomly" per tool) are separated in guaranteed transversal
            // position, beyond every kernel tolerance and the section slop.
            //
            // The asymmetry is the principle, not a fudge: the avoid inflation is
            // the SAFETY BOUND (dimensional truth, exactly R -- never biased), while
            // the fill expansion is only a REACH CLAIM, and over-claiming reach is
            // harmless -- the claim is still clipped by the exact avoid term, and
            // the extra overrun hangs in open air. Greed on the minuend is free;
            // the subtrahend must stay exact.
            float fillBias = 0.05f;       // model units (absolute: kernel tolerances are absolute)

            // Seed reach knob, as a MULTIPLE of the tool radius. The base fill
            // grow is seedOffset*R + fillBias (default 1 = exactly R + fillBias,
            // the plain tool-radius dilation). Larger claims more reach -- the
            // tool-centre region extends further past the true material, buying
            // open-air room. Reach only: the grow is still clipped by the exact
            // avoid term in the merge (avoid stays +R, the hard keep-out), so
            // over-claiming here only hangs extra centres in air, never drives
            // the tool into kept material. Only the BASE fill reads it; the
            // lead/air reach variants keep their own factors.
            float seedOffset = 1.0f;

            // Clip-seed FILL reach, a MULTIPLE of R (a hair over 1). This is a
            // generic-position PROXY of the true tool-centre wall -- exact-boundary
            // clipping (fill (+) exactly R) is numerically fragile (Weiler-Atherton
            // on coincident edges), so the clip region is displaced a hair OUTWARD
            // (0.1R) from the true wall. When seedOffset runs FAT (corner relief /
            // early exit into free-edge air), the cascade rings are GENERATED from
            // the fat seed but must OBEY this proxy of the true tool-centre region --
            // (fill (+) clipSeedOffset*R) \ (avoid (+) clipAvoidOffset*R). The
            // strategy clips each fat CASCADE ring to it, dropping the far-air
            // arclength. This proxy is the FINISH clip target; the strategy insets it
            // by wallStock to form the ROUGH clip target (the inset-seed proxy). The
            // 0.1R over 1R displaces the clip boundary past kernel eps from any
            // fat-seed-generated edge (generic position, same rationale as fillBias).
            // wallStock is NEVER on this seed -- it re-enters only on the ROUGH clip
            // proxy, strategy-side. Only built when the fat seed is meaningfully
            // larger (seedOffset > clipSeedOffset + eps); otherwise the two regions
            // coincide and the clip is a no-op -- clipSeed stays empty.
            float clipSeedOffset = 1.1f;

            // Clip-seed AVOID grow, a MULTIPLE of R (a hair UNDER 1). The FAT
            // seed grows avoid by exactly +R (the real hard keep-out), so a fat
            // ring's legitimate wall-finish stretch (tool tangent to the wall)
            // lies exactly on avoid (+) R. If the clip seed also grew avoid by
            // +R, that stretch would sit ON the clip's avoid boundary, where
            // winding classification is ill-defined -- the wall pass could be
            // dropped and the wall left uncut. Growing the clip's avoid a hair
            // LESS (0.9R) pushes its boundary 0.1R inside the fat wall stretch,
            // so that stretch is UNAMBIGUOUSLY kept (generic position, mirroring
            // the fill side's displacement). Safe: the fat seed already enforced
            // the true +R limit, so a looser clip avoid can only KEEP already-
            // legal wall pieces, never admit a path into avoid (+) R.
            float clipAvoidOffset = 0.9f;

            // DEPRECATED / UNUSED by the engine. Wall stock is NOT a seed term:
            // wallStock is a BINARY FINAL-EMISSION choice made STRATEGY-side. The
            // engine builds two seed regions -- the (fat) offset seed (profiles[0],
            // the free-form cascade origin) and the TRUE seed (trueSeed, the real
            // tool-centre wall at cut R / avoid R) -- and exposes trueSeed so the
            // strategy emits its outermost ring as either the true seed (finish) or
            // the true seed inset by the wall skin (rough). No cascade geometry ever
            // carries the skin. Retained only so existing initializers compile; the
            // engine seed no longer reads it.
            float wallStock = 0.0f;
        };

        Params params;

        Parallel() = default;
        explicit Parallel(Params p) : params(p) {}

        // Ingested input (consumed by run).
        std::vector<std::unique_ptr<Stoicheion>> fill;
        std::vector<std::unique_ptr<Stoicheion>> avoid;

        // Output: profiles[0] is the seed; each later profile is one stepover
        // deeper. Warnings are plain-language deviations, empty on a clean run.
        std::vector<Profile> profiles;
        std::vector<std::string> warnings;

        // The oriented, welded FILL section (material boundary at this slice),
        // kept after run() for consumers: winding >= 1 against these chains is
        // material, below is open air -- the test that decides where a ring may
        // be entered by a straight plunge (the tool centre column hangs in air).
        std::vector<Chain> fillBoundary;

        // The oriented, welded AVOID section, likewise kept -- primarily so a
        // consumer can SEE what actually bounded the seed (an empty list on a
        // level with kept material is the tell that the avoid never arrived).
        std::vector<Chain> avoidBoundary;

        // The grown systems, preserved pre-merge (the merge consumes its
        // inputs): the per-chain offsetValid results of fill and avoid.
        std::vector<Chain> grownFill;
        std::vector<Chain> grownAvoid;

        // The EXTENDED region: (fill (+) leadReach*R) minus (avoid (+) leadClear*R)
        // -- the survey of the open air lead motion may actually use, built by
        // the seed's own construction on the same preserved originals. Membership
        // is a wall-safety certificate WITH MARGIN: a lead arc lying wholly
        // inside it stays within feasible reach of the cut and keeps an extra
        // half radius off every kept wall (rings graze avoid + R; leads don't).
        Profile extended;

        // The CLIP seed: the true tool-centre region -- (fill (+) clipSeedOffset*R)
        // \ (avoid (+) clipAvoidOffset*R) -- that a FAT-seed run must obey. Rings
        // are generated from the reach-inflated seed (corner relief, early free-
        // edge exit) then clipped back to this region so the far-air arclength is
        // dropped. The avoid term is a hair UNDER +R (clipAvoidOffset < 1) so a
        // fat wall-hugging stretch (at the fat seed's real +R limit) falls inside
        // this clip and is kept, not dropped on an ill-defined boundary. The true
        // +R keep-out is still enforced by the FAT seed, so the looser clip avoid
        // only KEEPS legal travel, never licenses the tool past a wall. Empty unless
        // the fat seed is meaningfully larger than the clip (see clipSeedOffset).
        Profile clipSeed;

        // The TRUE seed: the real tool-centre wall region -- (fill (+) R) \
        // (avoid (+) R), i.e. seedOffset = 1 and NO wall stock. Built ALWAYS
        // (ungated), by the SAME seed machinery as profiles[0] but at the true
        // radius, from the same preserved originals. This is the AUTHORITATIVE
        // dimensional wall: the strategy emits it (finish) or its uniform inset
        // by the wall skin (rough) as the single outermost ring, and bounds link
        // containment against it. profiles[0] (the possibly-fat offset seed) stays
        // the free-form cascade origin; when seedOffset == 1 the two coincide, so
        // the strategy drops profiles[0] and emits trueSeed once (no duplicate).
        Profile trueSeed;

        // TRUE OPEN AIR -- the entry planner's menu of certified sites. The seed
        // treatment run once more at (fill (+) airReach*R) \ (avoid (+) airClear*R),
        // and the resulting boundary CLASSIFIED against the original sections:
        // pieces lying outside BOTH regions are kept, as OPEN chains (segments of
        // the air ring). On such a segment the whole tool disk hangs a full radius
        // clear of the cut material and a third of a diameter off every kept wall
        // -- a plunge there touches nothing. (The extended region cannot certify
        // this: fill (+) 2R CONTAINS the fill, so membership there proves reach
        // and wall safety, never disk clearance.)
        std::vector<Chain> openAir;

        // Pre-oriented avoid path: when the caller has already assembled and
        // oriented the keep-out as closed loops, feed them whole and SKIP the
        // engine's assemble/weld/dropOpen/orientByParity for the avoid side.
        std::vector<Chain> avoidLoops;
        bool avoidPreOriented = false;

        void ingest(std::vector<std::unique_ptr<Stoicheion>> fillSet,
                    std::vector<std::unique_ptr<Stoicheion>> avoidSet) {
            fill = std::move(fillSet);
            avoid = std::move(avoidSet);
        }

        // Opt-in: supply the avoid as pre-oriented CLOSED loops. INVARIANT --
        // every loop MUST be closed with role == Role::Subtrahend; the weld /
        // dropOpen safety net and the joint orientByParity are all skipped for
        // these, so an open or wrong-role loop is not corrected here.
        void ingestAvoidLoops(std::vector<Chain> loops) {
            avoidLoops = std::move(loops);
            avoidPreOriented = true;
        }

        // Orient a system by nesting parity: depth = how many OTHER loops of
        // the same system contain this one (probed at the loop's own extreme
        // point -- on the loop, far from the others unless they truly touch).
        // Desired handedness alternates from `base` with depth. Judged and
        // forced by exact SIGNED AREA (input orientation is a seed question;
        // turning is the coin-flip on tangency-cusped slivers). Degenerate
        // zero-area loops are dropped with a warning.
        void orientByParity(std::vector<Chain>& loops, int base, Role role) {

            std::vector<Chain> kept;
            for (size_t a = 0; a < loops.size(); a++) {

                float area = loops[a].signedArea();
                int cur = (area > 1e-4f) ? 1 : (area < -1e-4f ? -1 : 0);
                if (cur == 0) {
                    warnings.push_back("A degenerate (zero-area) input loop was dropped.");
                    continue;
                }

                int depth = 0;
                Pos probe = Chain::maxXPoint(loops[a]);
                for (size_t b = 0; b < loops.size(); b++) {
                    if (b == a) { continue; }
                    if (loops[b].windingAround(probe) != 0) { depth++; }
                }

                int desired = ((depth % 2) == 0) ? base : -base;
                Chain c = (cur == desired) ? std::move(loops[a]) : loops[a].reversed();
                c.role = role;
                kept.push_back(std::move(c));
            }
            loops = std::move(kept);
        }

        // Run the whole construction: assemble -> orient -> grow both ->
        // merge (the seed) -> inlay to extinction. Returns false (with a
        // warning) when no seed exists -- the avoid swallows the fill, or
        // the fill never assembled into a closed region.
        bool run() {

            profiles.clear();
            warnings.clear();
            fillBoundary.clear();
            avoidBoundary.clear();
            grownFill.clear();
            grownAvoid.clear();
            extended = Profile();
            clipSeed = Profile();
            openAir.clear();
            // NB: do NOT clear avoidPreOriented/avoidLoops here -- ingestAvoidLoops
            // runs before run(). They are consumed (moved-from) below; a fresh
            // Parallel per level means no stale carry-over.

            const float R = params.toolRadius;
            if (R <= 1e-6f) { warnings.push_back("Tool radius is zero: nothing to construct."); return false; }

            // 1. Assemble. Closed chains only past this door: offsetting is a
            // statement about bounded regions (welding almost-closed input is
            // the caller's conditioning job).
            const float fEps = std::max(params.eps, params.fillEps);
            const float aEps = std::max(params.eps, params.avoidEps);
            std::vector<Chain> fills  = Chain::assemble(fill, fEps);
            std::vector<Chain> avoids = avoidPreOriented ? std::move(avoidLoops)
                                                         : Chain::assemble(avoid, aEps);

            // Weld almost-closed loops shut: section endpoints carry slop past
            // the assembly eps, and a loop that is closed in fact but open by a
            // hair must not be dropped at the door. The gap is bridged with a
            // real segment (a marked-closed loop with a geometric hole would
            // break the crossing walk's coincidence matching downstream).
            auto weld = [&](std::vector<Chain>& v, float sysEps) {
                const float gapEps = std::max(params.weldEps, sysEps);   // a loop's own residual
                for (Chain& c : v) {                                     // gap scales with its solid
                    if (c.closed || c.edges.empty()) { continue; }
                    Pos a = Kern::eEnd(*c.edges.back());
                    Pos b = Kern::eStart(*c.edges.front());
                    if ((b - a).pythag() > gapEps) { continue; }
                    if ((b - a).pythag() > 1e-6f) { c.edges.push_back(std::make_unique<Segment2>(a, b)); }
                    c.closed = true;
                }
            };
            weld(fills, fEps);
            if (!avoidPreOriented) { weld(avoids, aEps); }   // pre-oriented loops are already closed

            // Dropped opens are the assembly's failure report: say HOW each one
            // failed (edge count + closure gap), so a consumer can tell endpoint
            // slop (gap just past the weld) from garbage sections (huge gaps --
            // the tell of a tilted or drifted slicing frame).
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
            // avoid: pre-oriented loops keep the caller's own-parity handedness;
            // re-deriving JOINT parity here would flip deep nested walls (gouge).
            if (!avoidPreOriented) { orientByParity(avoids, -1, Role::Subtrahend); }   // avoid: outers CW, islands CCW

            // 3. Grow BOTH by one tool radius. Fill by -R (left of CCW travel
            // is the interior, so -R expands); avoid by +R (left of CW travel
            // is the exterior, so +R expands). Islands of each system shrink
            // under the same value -- dilation of the region, by one blind
            // sign per system (Method/06). The fill takes the extra fillBias
            // (see Params): reach is claimed greedily, the safety bound stays
            // exact, and grown boundaries never sit coincident.
            std::vector<Chain> grown;
            for (Chain& c : fills) {
                for (Chain& m : c.offsetValid(-(params.seedOffset * R + params.fillBias), params.eps)) {
                    grownFill.push_back(m.clone());     // preserved: the merge consumes its input
                    grown.push_back(std::move(m));
                }
            }
            for (Chain& c : avoids) {
                for (Chain& m : c.offsetValid(+R, params.eps)) {
                    grownAvoid.push_back(m.clone());
                    grown.push_back(std::move(m));
                }
            }

            // The sections have served the growth; keep them as the consumer's
            // oracles (fillBoundary: material vs air; avoidBoundary: what bounded).
            fillBoundary = std::move(fills);
            avoidBoundary = std::move(avoids);

            // 4. The seed: one merge. { winding >= 1 } with roles is exactly
            // (fill (+) R) \ (avoid (+) R). Envelopment-without-intersection
            // (the bore: the grown avoid ring sitting whole inside the grown
            // fill) is handled by the merge's containment clustering -- the
            // avoid's ring survives as the seed, the fill's copy is discarded.
            Profile seed;
            seed.chains = Chain::merge(std::move(grown), params.eps);

            // The seed bounds the clearable region; the construction roles are
            // spent. Everything from here on is minuend.
            for (Chain& c : seed.chains) { c.role = Role::Minuend; }

            if (seed.empty() || seed.signedArea() <= 0.0f) {
                warnings.push_back("No seed: the avoid region leaves no reachable extent of the fill "
                                   "(or the fill collapsed under the tool radius).");
                return false;
            }

            // The extended region: the identical construction from the same
            // preserved originals, at the lead expansions (see Params). The
            // fill side carries the bias (generic position, greed is free);
            // the avoid side is a clean 1.5R -- above R is the safe direction.
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

            // The clip seed: the seed's OWN construction once more, at the true
            // tool radius (clipSeedOffset, a hair over 1R) rather than the fat
            // seedOffset. The fill side carries the same fillBias (generic
            // position); the avoid side grows by clipAvoidOffset*R -- a hair
            // UNDER +R, so the fat seed's wall-hugging stretches (at avoid (+) R,
            // the real limit the fat seed enforced) sit unambiguously INSIDE
            // this clip region and are kept, instead of coincident with its
            // boundary (where winding is ill-defined). The true +R keep-out is
            // still enforced by the FAT seed; this looser clip avoid can only
            // keep already-legal wall pieces, never admit avoid (+) R. Built
            // only when the fat seed is meaningfully larger; when they coincide
            // the clip would be a no-op, so it is skipped and clipSeed left empty.
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

            // The TRUE seed: the seed's OWN construction once more, at the true
            // tool radius (seedOffset = 1) and NO wall stock -- (fill (+) R) \
            // (avoid (+) R). Built UNCONDITIONALLY (not gated on a fat seed): it
            // is the authoritative dimensional wall the strategy emits as its
            // outermost ring (finish = this; rough = this inset by the wall skin)
            // and bounds link containment against. The fill side keeps the same
            // fillBias for generic position; the avoid side is a clean +R (the
            // real hard keep-out). Distinct from clipSeed, whose avoid is a hair
            // UNDER +R for winding-safe clipping -- this is the exact wall.
            {
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

            // The true-open-air survey: the identical construction once more, at
            // the air expansions, from the same preserved originals. The merged
            // ring is part over the cut region (its avoid-derived stretches can
            // run inside the fill) and part in genuine air; each edge is split at
            // its crossings with the ORIGINAL section boundaries and the pieces
            // classified by winding -- outside the fill region (w < 1, CCW base)
            // AND outside the avoid region (w > -1, CW base) is true air. The
            // survivors stitch into open segments, oriented, never reversed.
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

            // 5. The inlay: recursively inset until extinction or the cap.
            // profiles[0] is already the outermost CUTTABLE ring -- the seed
            // stands one radius off every wall by construction -- so every
            // advance is a plain stepover. A generation whose total signed
            // area goes non-positive claims an unbounded region (outsides
            // without insides): extinction, per the sanity rule.
            const float step = R * std::max(0.05f, params.stepover);
            while (static_cast<int>(profiles.size()) <= params.maxGenerations) {
                Profile next = profiles.back().offsetBy(step, params.eps);
                if (next.empty()) { break; }
                next.included = next.signedArea() > 0.0f;   // skip downstream, never halt
                profiles.push_back(std::move(next));
            }
            return true;
        }
    };

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
