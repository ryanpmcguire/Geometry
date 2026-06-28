module;

#include <vector>
#include <memory>
#include <string>
#include <utility>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <cmath>
#include <algorithm>

#include <nlohmann/json.hpp>

export module Geo.Geometry;

import Rev.Core.Pos;

// The shared 2D geometry kernel -- stoicheia, relations, the sketch geometry
// container -- consumed by both the Sketch app and the CAM app. Application
// modules re-export these names into their own namespaces (see the
// Sketch.App.Geometry shim on the Sketch side).
export namespace Geo {

    using Json = nlohmann::json;
    using Rev::Core::Pos;

    inline constexpr float PI  = 3.14159265358979323846f;
    inline constexpr float TAU = 2.0f * PI;

    // Stable identity
    //--------------------------------------------------
    // Every entity and every relation carries a unique 64-bit id, persisted as a
    // 16-digit hex string. Ids are the durable handle relations reference by; the
    // pointer graph between them is a derived cache rebuilt from ids on load.

    using Id = std::uint64_t;

    // A fresh random id (never 0, since 0 means "unassigned").
    inline Id newId() {
        static std::mt19937_64 rng(std::random_device{}());
        Id v = 0;
        while (v == 0) { v = rng(); }
        return v;
    }

    inline std::string idToHex(Id id) {
        char buf[17];
        std::snprintf(buf, sizeof(buf), "%016llx", static_cast<unsigned long long>(id));
        return std::string(buf);
    }

    inline Id idFromHex(const std::string& s) {
        return s.empty() ? 0 : static_cast<Id>(std::strtoull(s.c_str(), nullptr, 16));
    }

    // A JSON-serialisable 2D position: Rev's Pos plus a clean round-trip to
    // { "x", "y" }. The stored point type for every geometric entity, so the
    // geometry gets Pos's vector algebra (b - a, normalize, dot, cross, ...).
    struct Pos2 : public Pos {

        using Pos::Pos;                                   // inherit (x,y) constructors
        Pos2() = default;
        Pos2(const Pos& p) : Pos(p) {}

        Json toJson() const { return Json{ { "x", x }, { "y", y } }; }

        static Pos2 fromJson(const Json& j) {
            return Pos2(j.value("x", 0.0f), j.value("y", 0.0f));
        }
    };

    // Ellipse math (shared by entities and tools)
    //--------------------------------------------------

    // World point at parameter t: P(t) = C + cos(t)*U + sin(t)*V.
    inline Pos ellipsePointAt(const Pos& c, const Pos& u, const Pos& v, float t) {
        return c + u * std::cos(t) + v * std::sin(t);
    }

    // The *true* offset (parallel-curve) point at parameter t: P(t) displaced by `d`
    // along its real outward unit normal. The normal is perpendicular to the tangent
    // dP/dt = -sin t U + cos t V (correct even for conjugate, non-orthogonal axes),
    // oriented outward (away from centre). `d` is signed: +grows, -insets. This is
    // analytically equidistant -- the exact offset curve, which is NOT an ellipse.
    inline Pos ellipseOffsetPoint(const Pos& c, const Pos& u, const Pos& v, float t, float d) {
        Pos p = ellipsePointAt(c, u, v, t);
        Pos tan = u * (-std::sin(t)) + v * std::cos(t);
        Pos n(-tan.y, tan.x);
        float nl = n.pythag();
        if (nl > 1e-9f) { n = n / nl; }
        if (n.dot(p - c) < 0.0f) { n = n * -1.0f; }   // outward (away from centre)
        return p + n * d;
    }

    // Parameter t of a world point. With orthogonal conjugate axes, cos t / sin t
    // are the projections onto U and V normalised by their squared lengths.
    inline float ellipseParamOf(const Pos& c, const Pos& u, const Pos& v, const Pos& w) {
        Pos d = w - c;
        float lu = u.dot(u), lv = v.dot(v);
        float cs = (lu > 1e-12f) ? d.dot(u) / lu : 0.0f;
        float sn = (lv > 1e-12f) ? d.dot(v) / lv : 0.0f;
        return std::atan2(sn, cs);
    }

    // Fit the unique ellipse from a centre, a major-axis point A (fixing U = A - C)
    // and a perimeter point B (fixing the conjugate semi-axis V, perpendicular to U).
    struct EllipseFit { Pos u, v; };

    inline EllipseFit fitEllipse(const Pos& c, const Pos& a, const Pos& b) {
        EllipseFit f;
        f.u = a - c;
        float al = f.u.pythag();
        if (al < 1e-6f) { return f; }
        Pos m = f.u / al, perp = m.normal();
        Pos d = b - c;
        float p = d.dot(m), q = d.dot(perp);
        float k = 1.0f - (p / al) * (p / al);
        float bb = (k > 1e-6f) ? std::fabs(q) / std::sqrt(k) : std::fabs(q);
        f.v = perp * bb;
        return f;
    }

    // Circumcentre of three points (the centre of the unique circle through them).
    // Falls back to the centroid when they are (near) collinear.
    inline Pos circumcentre(Pos p0, Pos p1, Pos p2) {
        float d = 2.0f * (p0.x * (p1.y - p2.y) + p1.x * (p2.y - p0.y) + p2.x * (p0.y - p1.y));
        if (std::fabs(d) < 1e-9f) { return (p0 + p1 + p2) / 3.0f; }
        float a2 = p0.dot(p0), b2 = p1.dot(p1), c2 = p2.dot(p2);
        float ux = (a2 * (p1.y - p2.y) + b2 * (p2.y - p0.y) + c2 * (p0.y - p1.y)) / d;
        float uy = (a2 * (p2.x - p1.x) + b2 * (p0.x - p2.x) + c2 * (p1.x - p0.x)) / d;
        return Pos(ux, uy);
    }

    // Wrap an angle delta into [0, TAU).
    inline float wrapTau(float x) { while (x < 0.0f) { x += TAU; } while (x >= TAU) { x -= TAU; } return x; }

    // Sampling helpers (write a polyline approximation of a curve into `out`).
    //--------------------------------------------------

    inline int spanSteps(float span) {
        return std::max(2, static_cast<int>(std::ceil(std::fabs(span) / (TAU / 128.0f))));
    }

    inline void sampleArc(std::vector<Pos>& out, Pos c, float r, float a0, float span, int steps, bool skipFirst = false) {
        for (int i = skipFirst ? 1 : 0; i <= steps; i++) {
            float a = a0 + span * (static_cast<float>(i) / steps);
            out.push_back(c + Pos::fromAngle(a) * r);
        }
    }

    inline void sampleEllipse(std::vector<Pos>& out, Pos c, Pos u, Pos v, float t0, float span, int steps, bool skipFirst = false) {
        for (int i = skipFirst ? 1 : 0; i <= steps; i++) {
            float t = t0 + span * (static_cast<float>(i) / steps);
            out.push_back(ellipsePointAt(c, u, v, t));
        }
    }

    inline Pos closestOnSegment(Pos p, Pos a, Pos b) {
        Pos d = b - a;
        float len2 = d.dot(d);
        float t = (len2 > 1e-9f) ? std::clamp((p - a).dot(d) / len2, 0.0f, 1.0f) : 0.0f;
        return a + d * t;
    }

    // Exact intersection of straight segments / circular arcs. These operate on the
    // *true* analytic forms only -- nothing is ever quantised. Each appends 0..2
    // exact crossing points to `out`. (A full circle is the arc with sweep == TAU.)

    inline bool arcContains(Pos x, Pos c, float a0, float sweep) {
        if (sweep >= TAU - 1e-4f) { return true; }
        return wrapTau((x - c).angle() - a0) <= sweep + 1e-4f;
    }

    inline void isectSegSeg(Pos p1, Pos p2, Pos p3, Pos p4, std::vector<Pos>& out) {
        Pos r = p2 - p1, s = p4 - p3;
        float rxs = r.cross(s);
        if (std::fabs(rxs) < 1e-12f) { return; }                    // parallel / degenerate
        Pos qp = p3 - p1;
        float t = qp.cross(s) / rxs, u = qp.cross(r) / rxs;
        if (t < -1e-4f || t > 1.0f + 1e-4f || u < -1e-4f || u > 1.0f + 1e-4f) { return; }
        out.push_back(p1 + r * t);
    }

    inline void isectSegArc(Pos p1, Pos p2, Pos c, float rad, float a0, float sweep, std::vector<Pos>& out) {
        Pos d = p2 - p1, f = p1 - c;
        float A = d.dot(d);
        if (A < 1e-12f) { return; }
        float B = 2.0f * f.dot(d);
        float C = f.dot(f) - rad * rad;
        float disc = B * B - 4.0f * A * C;
        if (disc < 0.0f) { return; }
        disc = std::sqrt(disc);
        float ts[2] = { (-B - disc) / (2.0f * A), (-B + disc) / (2.0f * A) };
        for (float t : ts) {
            if (t < -1e-4f || t > 1.0f + 1e-4f) { continue; }
            Pos x = p1 + d * t;                                     // exact point on the line *and* circle
            if (arcContains(x, c, a0, sweep)) { out.push_back(x); }
        }
    }

    inline void isectArcArc(Pos c1, float r1, float a1, float s1,
                            Pos c2, float r2, float a2, float s2, std::vector<Pos>& out) {
        Pos d = c2 - c1;
        float dd = d.pythag();
        if (dd < 1e-6f) { return; }                                 // concentric
        if (dd > r1 + r2 + 1e-4f) { return; }                       // too far apart
        if (dd < std::fabs(r1 - r2) - 1e-4f) { return; }            // one inside the other
        float a = (r1 * r1 - r2 * r2 + dd * dd) / (2.0f * dd);
        float h2 = r1 * r1 - a * a;
        float h = (h2 > 0.0f) ? std::sqrt(h2) : 0.0f;
        Pos mid = c1 + d / dd * a;
        Pos perp = Pos(-d.y, d.x) / dd;
        for (float sgn : { 1.0f, -1.0f }) {
            Pos x = mid + perp * (h * sgn);
            if (arcContains(x, c1, a1, s1) && arcContains(x, c2, a2, s2)) { out.push_back(x); }
        }
    }

    // A single possible place the cursor could snap to, with a score. The view
    // gathers these from every nearby entity (plus intersections) each frame and
    // simply picks the highest-scoring one within reach -- and can draw them.
    //
    // Priority, highest first: a defined *point* (endpoint, centre, midpoint,
    // quadrant, vertex) > an *intersection* of two curves > the perpendicular
    // *foot* on a single curve. The tier (axis 2 > real 1 > construction 0) is a
    // minor tiebreak within a category.
    struct SnapCandidate {
        enum Kind { Point, Intersection, OnCurve };
        Pos pos;
        Kind kind = Point;
        int score = 0;
    };

    inline int snapScore(SnapCandidate::Kind kind, int tier) {
        int base = (kind == SnapCandidate::Point) ? 300
                 : (kind == SnapCandidate::Intersection) ? 200 : 100;
        return base + tier;
    }

    // The concrete stoicheion types, as an enum: the cheap, string-free type
    // probe (kind() remains only as the serialisation tag).
    enum class SKind : std::uint8_t {
        Point, Segment, Circle, Arc, Ellipse, EllipseArc, OffsetEllipse, Spline
    };

    // The display / category group of a stoicheion in generated (offset) layers:
    // an enum stamped at creation, never strings. Crossing carries its level in
    // the stoicheion's groupLevel.
    enum class DisplayGroup : std::uint8_t {
        None,
        Profile0,        // the prepared seed chain (yellow)
        ValidPos,        // valid chain, positive signed area (red)
        ValidNeg,        // valid chain, negative signed area (green)
        ToolpathSkip,    // valid by the method, skipped by the toolpath sanity rule
        Crossing,        // debug: crossing-number level colouring (level in groupLevel)
        Minimum,         // debug: the extracted minimum-level path (magenta)
        Intersections,   // debug: crossing dots (yellow)
        Start            // debug: the walk's origin (pink)
    };

    // ===============================================================
    // Stoicheion (στοιχεῖον) -- the fundamental geometric element.
    //
    // Euclid's "Elements" is the Stoicheia; a single primitive is a Stoicheion.
    // Each concrete type stores its own data and knows how to: serialise itself,
    // draw itself as a polyline (tessellate), measure distance to a probe point,
    // expose its anchor (snap) points and its mutable control points, and project
    // a point onto its body. Selection / rendering / snapping / moving in the view
    // are then a single loop of virtual calls -- no per-type branching.
    // ===============================================================
    struct Stoicheion {

        Id id = 0;                   // stable unique identity (assigned on commit)
        bool construction = false;   // reference geometry, not real output
        bool locked = false;         // datum (origin / axes): can't move or delete

        // Open air: beyond this edge lies free space, not part material -- the
        // tool may run off it. A profile strategy pre-pushes open-air edges
        // outward before recursing, so corners against the open region are fully
        // covered instead of left with rounded internal corners.
        bool openAir = false;

        // Group tag: which display category within its layer this stoicheion
        // belongs to. Stamped at creation by whatever operation produced the
        // entity -- never derived downstream; the view maps tag -> style blindly,
        // the same way construction -> grey. Crossing levels ride in groupLevel.
        DisplayGroup group = DisplayGroup::None;
        int groupLevel = 0;

        virtual ~Stoicheion() = default;

        // Identity / serialisation
        virtual SKind type() const = 0;         // the enum probe (cheap, hot paths)
        virtual const char* kind() const = 0;   // the serialisation tag only
        virtual std::unique_ptr<Stoicheion> clone() const = 0;
        virtual Json data() const = 0;          // type-specific payload
        virtual void load(const Json& j) = 0;   // read type-specific payload

        Json toJson() const {
            Json j = data();
            j["id"] = idToHex(id);
            j["kind"] = kind();
            j["construction"] = construction;
            j["locked"] = locked;
            if (openAir) { j["openAir"] = true; }
            if (group != DisplayGroup::None) {
                j["group"] = static_cast<int>(group);
                j["groupLevel"] = groupLevel;
            }
            return j;
        }

        // Rendering: append this entity's polyline approximation (world space).
        // Point-like entities append a single position, drawn as a dot. `skipFirst`
        // omits the entity's first point — used when chaining edges end-to-end so a
        // shared corner isn't emitted twice (no zero-length segments at the joints).
        virtual void tessellate(std::vector<Pos>& out, bool skipFirst = false) const = 0;
        virtual bool isPoint() const { return false; }

        // This edge's exact contribution to the closed loop's signed-area integral,
        // ∮ (x dy - y dx), evaluated analytically along the travel direction -- NOT by
        // sampling. Summed around a closed chain it is twice the signed area, and its
        // sign is the winding (+ CCW, - CW). Each primitive integrates in closed form
        // (a straight chord term, an arc's r^2*dtheta, an ellipse's (uxv)*dt); a point
        // encloses nothing. This is pure geometry, exact to floating point.
        virtual float signedAreaTerm() const { return 0.0f; }

        // Hit-testing: distance from p to the drawn curve.
        virtual float distanceTo(Pos p) const = 0;

        // Point-priority anchors (endpoints, centre) -- used for *selection*.
        virtual void anchors(std::vector<Pos>& out) const {}

        // Feature snap points: every salient position to magnetise placement onto
        // (endpoints, centres, midpoints, circle quadrants, ellipse axis-ends, ...).
        // Defaults to the selection anchors; entities override to add their extras.
        virtual void snapPoints(std::vector<Pos>& out) const { anchors(out); }

        // Exact-intersection forms, used to cross one entity with another. An
        // entity that is *genuinely* straight yields its segment(s); one that is
        // *genuinely* circular yields its circle/arc. Anything else (ellipses)
        // yields neither -- it produces no intersections rather than being
        // approximated. Nothing here is ever quantised.
        virtual void asLineSegments(std::vector<std::pair<Pos, Pos>>& out) const {}
        virtual bool asCircle(Pos& c, float& r, float& a0, float& sweep) const { return false; }

        // Yield this entity's own snap candidates for a cursor at `mouse`: every
        // feature point (scored Point) and the perpendicular foot (scored OnCurve).
        // Intersections are between *pairs* of entities, so the view adds those.
        virtual void snapCandidates(Pos mouse, std::vector<SnapCandidate>& out) const {
            int tier = locked ? 2 : (construction ? 0 : 1);
            std::vector<Pos> sps;
            snapPoints(sps);
            for (const Pos& p : sps) { out.push_back({ p, SnapCandidate::Point, snapScore(SnapCandidate::Point, tier) }); }
            Pos foot;
            if (footOnCurve(mouse, foot)) {
                out.push_back({ foot, SnapCandidate::OnCurve, snapScore(SnapCandidate::OnCurve, tier) });
            }
        }

        // Nearest point on the body for snapping onto the curve; false if none.
        virtual bool footOnCurve(Pos p, Pos& out) const { return false; }

        // Mutable control points -- every absolute point defining the entity (used
        // both for moving and for coincidence detection between entities).
        virtual void controlPoints(std::vector<Pos2*>& out) = 0;

        // Offsetting. `offset` moves this curve sideways by `amount`, *in place*:
        // positive grows / outsets, negative insets (toward the curve's interior --
        // a smaller circle, a parallel line). `offsetBy` is the value-returning twin:
        // it leaves this one untouched and hands back a fresh offset copy. The pair
        // mirrors Pos::normalize / Pos::normalized. A directionless element (a point)
        // does not move. No chaining, no constraints -- each stoicheion offsets alone.
        virtual void offset(float amount) {}
        virtual std::unique_ptr<Stoicheion> offsetBy(float amount) const {
            std::unique_ptr<Stoicheion> copy = clone();
            copy->offset(amount);
            return copy;
        }

        // Offset into a (possibly multi-piece) result. A line / circle / arc offsets
        // exactly to one shape, so the default just emits offsetBy(). Curves whose
        // true offset is not the same kind of curve (ellipses) override this to emit
        // the analytically-equidistant offset as a polyline.
        virtual void offsetInto(float amount, std::vector<std::unique_ptr<Stoicheion>>& out) const {
            out.push_back(offsetBy(amount));
        }

        // Chirality marker housekeeping. A circle / arc carries its orientation as an
        // extra point (b / d). That point is NOT a solver variable and is NOT one of
        // the entity's "defining" (solvable) points -- it is just consistently stored
        // and re-seated to its canonical spot after any change, so it always reads at
        // the proper location and rides along when the entity moves.
        //   chiralitySlot : the anchor/control-point slot of that marker, or -1.
        //   normalize     : re-seat the marker canonically (preserving its side).
        virtual int chiralitySlot() const { return -1; }
        virtual void normalize() {}
    };

    // Exact mutual intersection of two entities. Uses only their true analytic
    // forms (straight segments / circular arcs); if either is neither (an ellipse),
    // no points are produced -- never an approximation.
    inline void intersect(const Stoicheion& A, const Stoicheion& B, std::vector<Pos>& out) {

        std::vector<std::pair<Pos, Pos>> as, bs;
        A.asLineSegments(as);
        B.asLineSegments(bs);

        Pos ac, bc; float ar, br, aa0, asw, ba0, bsw;
        bool aCirc = A.asCircle(ac, ar, aa0, asw);
        bool bCirc = B.asCircle(bc, br, ba0, bsw);

        for (const auto& sa : as) { for (const auto& sb : bs) { isectSegSeg(sa.first, sa.second, sb.first, sb.second, out); } }
        if (bCirc) { for (const auto& sa : as) { isectSegArc(sa.first, sa.second, bc, br, ba0, bsw, out); } }
        if (aCirc) { for (const auto& sb : bs) { isectSegArc(sb.first, sb.second, ac, ar, aa0, asw, out); } }
        if (aCirc && bCirc) { isectArcArc(ac, ar, aa0, asw, bc, br, ba0, bsw, out); }
    }

    // Point
    //--------------------------------------------------
    struct Point2 : public Stoicheion {

        Pos2 p;

        Point2() = default;
        Point2(const Pos& pos) : p(pos) {}

        SKind type() const override { return SKind::Point; }
        const char* kind() const override { return "point"; }
        std::unique_ptr<Stoicheion> clone() const override { return std::make_unique<Point2>(*this); }
        Json data() const override { return Json{ { "p", p.toJson() } }; }
        void load(const Json& j) override { p = Pos2::fromJson(j.value("p", Json::object())); }

        void tessellate(std::vector<Pos>& out, bool skipFirst = false) const override { if (!skipFirst) { out.push_back(p); } }
        bool isPoint() const override { return true; }
        float distanceTo(Pos q) const override { return (q - p).pythag(); }
        void anchors(std::vector<Pos>& out) const override { out.push_back(p); }
        void controlPoints(std::vector<Pos2*>& out) override { out.push_back(&p); }
    };

    // Segment
    //--------------------------------------------------
    struct Segment2 : public Stoicheion {

        Pos2 a, b;

        Segment2() = default;
        Segment2(const Pos& a, const Pos& b) : a(a), b(b) {}

        SKind type() const override { return SKind::Segment; }
        const char* kind() const override { return "segment"; }
        std::unique_ptr<Stoicheion> clone() const override { return std::make_unique<Segment2>(*this); }
        Json data() const override { return Json{ { "a", a.toJson() }, { "b", b.toJson() } }; }
        void load(const Json& j) override {
            a = Pos2::fromJson(j.value("a", Json::object()));
            b = Pos2::fromJson(j.value("b", Json::object()));
        }

        void tessellate(std::vector<Pos>& out, bool skipFirst = false) const override { if (!skipFirst) { out.push_back(a); } out.push_back(b); }
        // Exact: integral of (x dy - y dx) along the straight chord a -> b.
        float signedAreaTerm() const override { return a.x * b.y - a.y * b.x; }
        float distanceTo(Pos p) const override { return (p - closestOnSegment(p, a, b)).pythag(); }
        void anchors(std::vector<Pos>& out) const override { out.push_back(a); out.push_back(b); }
        void snapPoints(std::vector<Pos>& out) const override {
            out.push_back(a); out.push_back(b); out.push_back((a + b) * 0.5f);   // + midpoint
        }
        void asLineSegments(std::vector<std::pair<Pos, Pos>>& out) const override { out.push_back({ a, b }); }
        bool footOnCurve(Pos p, Pos& out) const override { out = closestOnSegment(p, a, b); return true; }
        void controlPoints(std::vector<Pos2*>& out) override { out.push_back(&a); out.push_back(&b); }

        // Slide the line sideways by `amount` along its left-hand normal (the +90deg
        // rotation of a->b); a negative amount slides the other way.
        void offset(float amount) override {
            Pos d = b - a; float len = d.pythag();
            if (len < 1e-9f) { return; }
            Pos n(-d.y / len, d.x / len);
            a -= n * amount;
            b -= n * amount;
        }
    };

    // Circle -- centre c, radius point a (|a - c| is the radius), and chirality
    // marker b: a point on the circle a quarter-turn (pi/4) from a. The *side* b sits
    // on is the orientation -- b a quarter-turn CCW of a is a CCW circle, CW is a CW
    // circle. Exactly like an arc (c, a, b), orientation lives in the points: no
    // signed radius, no bool. So a CW circle offsets opposite to a CCW one, and an
    // inset that drives the radius through zero turns the circle inside-out -- the
    // radius re-emerges positive on the far side and the chirality flips.
    struct Circle2 : public Stoicheion {

        static constexpr float Right = TAU * 0.25f;   // pi/2 (a quarter turn)

        Pos2 c, a, b;

        Circle2() = default;
        // Centre + radius: a CCW circle (a on +x, b a quarter-turn CCW of it).
        Circle2(const Pos& c, float r)
            : c(c), a(c + Pos(r, 0.0f)), b(c + Pos::fromAngle(Right) * r) {}
        Circle2(const Pos& c, const Pos& a, const Pos& b) : c(c), a(a), b(b) {}

        // The unique circle through three points; its chirality follows their winding
        // (p0 -> p1 -> p2 counterclockwise gives a CCW circle).
        static Circle2 ThroughThreePoints(Pos p0, Pos p1, Pos p2) {
            Pos o = circumcentre(p0, p1, p2);
            Pos ra = p0 - o;
            bool ccw = (p1 - p0).cross(p2 - p0) >= 0.0f;
            Pos marker = ccw ? (o + Pos(-ra.y, ra.x)) : (o + Pos(ra.y, -ra.x));   // b at +/-90 deg
            return Circle2(o, p0, marker);
        }
        // The circle having `a`-`b` as a diameter (CCW).
        static Circle2 Diameter(Pos a, Pos b) {
            return Circle2((a + b) * 0.5f, (a - b).pythag() * 0.5f);
        }

        float radius() const { return (a - c).pythag(); }
        int chirality() const { return ((a - c).cross(b - c) >= 0.0f) ? 1 : -1; }   // +1 CCW, -1 CW

        SKind type() const override { return SKind::Circle; }
        const char* kind() const override { return "circle"; }
        std::unique_ptr<Stoicheion> clone() const override { return std::make_unique<Circle2>(*this); }
        Json data() const override { return Json{ { "c", c.toJson() }, { "a", a.toJson() }, { "b", b.toJson() } }; }
        void load(const Json& j) override {
            c = Pos2::fromJson(j.value("c", Json::object()));
            a = Pos2::fromJson(j.value("a", Json::object()));
            b = Pos2::fromJson(j.value("b", Json::object()));
        }

        void tessellate(std::vector<Pos>& out, bool skipFirst = false) const override {
            sampleArc(out, c, radius(), (a - c).angle(), TAU, spanSteps(TAU), skipFirst);
        }
        // Exact: a full turn contributes (chirality) * r^2 * TAU (the centre terms
        // vanish since start == end). Twice the area is +-pi*r^2 * 2, i.e. r^2*TAU.
        float signedAreaTerm() const override {
            float r = radius();
            return static_cast<float>(chirality()) * r * r * TAU;
        }
        float distanceTo(Pos p) const override { return std::fabs((p - c).pythag() - radius()); }
        void anchors(std::vector<Pos>& out) const override { out.push_back(c); out.push_back(a); out.push_back(b); }
        void snapPoints(std::vector<Pos>& out) const override {
            float r = radius();
            out.push_back(c); out.push_back(a); out.push_back(b);       // centre + handles
            out.push_back(c + Pos(r, 0.0f)); out.push_back(c + Pos(-r, 0.0f));   // + quadrants
            out.push_back(c + Pos(0.0f, r)); out.push_back(c + Pos(0.0f, -r));
        }
        bool asCircle(Pos& cc, float& rr, float& a0, float& sweep) const override {
            cc = c; rr = radius(); a0 = (a - c).angle(); sweep = TAU; return true;
        }
        bool footOnCurve(Pos p, Pos& out) const override {
            float r = radius();
            Pos d = p - c; float dd = d.pythag();
            out = (dd > 1e-6f) ? (c + d / dd * r) : (c + Pos(r, 0.0f));
            return true;
        }
        void controlPoints(std::vector<Pos2*>& out) override { out.push_back(&c); out.push_back(&a); out.push_back(&b); }

        // b is the chirality marker, not a solved point, and is kept canonical: a
        // quarter-turn (pi/2) from a, on the side it currently sits. Re-seated after
        // any move/solve so it always reads "at the proper location".
        int chiralitySlot() const override { return 2; }
        void normalize() override {
            Pos ra = a - c;
            b = (chirality() > 0) ? (c + Pos(-ra.y, ra.x))    // +90 deg (CCW)
                                  : (c + Pos(ra.y, -ra.x));   // -90 deg (CW)
        }

        // Concentric offset, chirality-aware. A CCW circle shrinks for a negative
        // amount; a CW circle grows -- opposite, just as a CW chain offsets opposite
        // to a CCW one. Driving the radius through zero turns the circle inside-out:
        // it re-emerges on the far side with the chirality flipped (b swaps sides).
        void offset(float amount) override {
            float r = radius();
            if (r < 1e-9f) { return; }
            int s = chirality();
            float nr = r + amount * static_cast<float>(s);
            float ang = (a - c).angle();
            if (nr < 0.0f) { nr = -nr; ang += TAU * 0.5f; s = -s; }     // through the singularity
            a = c + Pos::fromAngle(ang) * nr;
            b = c + Pos::fromAngle(ang + static_cast<float>(s) * Right) * nr;
        }
    };

    // Arc -- centre c, endpoints a and b, and a *through-point* d that fixes which
    // way round we travel: the arc runs "from A to B through D". Whether D sits on the
    // CCW or CW side of the chord (about the centre) is the chirality -- so the same
    // pair of endpoints can be the CCW arc or the CW arc, and it offsets accordingly,
    // just like a circle's marker point. Only D's *angle* about the centre matters;
    // its distance is free.
    struct Arc2 : public Stoicheion {

        Pos2 c, a, b, d;

        Arc2() = default;
        // Two-endpoint form: a CCW arc, D placed at the CCW midpoint.
        Arc2(const Pos& c, const Pos& a, const Pos& b) : c(c), a(a), b(b) {
            float ang = (a - c).angle() + wrapTau((b - c).angle() - (a - c).angle()) * 0.5f;
            d = c + Pos::fromAngle(ang) * (a - c).pythag();
        }
        Arc2(const Pos& c, const Pos& a, const Pos& b, const Pos& d) : c(c), a(a), b(b), d(d) {}

        // Construct the arc through `a` and `b` that is tangent to `tanA` at a and
        // `tanB` at b. The centre is where the perpendiculars to the two tangents
        // meet (well-defined unless the tangents are parallel) -- and for a fillet at
        // an offset corner that meeting point is exactly the original (pre-inset)
        // corner, with a and b equidistant from it. The travel direction (D / the
        // sweep) is chosen so the arc leaves a heading along tanA.
        static Arc2 Tangent(Pos a, Pos b, Pos tanA, Pos tanB) {
            Pos na(-tanA.y, tanA.x), nb(-tanB.y, tanB.x);   // perpendiculars to the tangents
            float rxs = na.cross(nb);
            Pos o = (std::fabs(rxs) > 1e-9f) ? (a + na * ((b - a).cross(nb) / rxs))
                                             : ((a + b) * 0.5f);   // parallel: degenerate fallback
            float radius = (a - o).pythag();
            Pos radial = a - o;
            Pos ccwTan(-radial.y, radial.x);                       // CCW tangent at a
            int s = (tanA.dot(ccwTan) >= 0.0f) ? 1 : -1;           // does tanA head CCW?
            float angA = radial.angle();
            float span = (s > 0) ? wrapTau((b - o).angle() - angA) : wrapTau(angA - (b - o).angle());
            float midAng = angA + static_cast<float>(s) * span * 0.5f;
            return Arc2(o, a, b, o + Pos::fromAngle(midAng) * radius);
        }

        // The arc that passes through three points, from p0 to p2 by way of p1 -- p1
        // is the through-point, so it fixes both which arc and its chirality.
        static Arc2 ThreePoint(Pos p0, Pos p1, Pos p2) {
            return Arc2(circumcentre(p0, p1, p2), p0, p2, p1);
        }

        float radius() const { return (a - c).pythag(); }

        // +1 if we travel CCW from A to B (D is on the CCW side), -1 if CW.
        int chirality() const {
            float angA = (a - c).angle();
            float ccwB = wrapTau((b - c).angle() - angA);
            float ccwD = wrapTau((d - c).angle() - angA);
            return (ccwD <= ccwB) ? 1 : -1;
        }

        // The covered arc as a CCW start angle + positive sweep (the half through D).
        // Drawing / hit-testing / intersection care about the point set, not the
        // travel direction, so they use this.
        void range(float& a0, float& sweep) const {
            float angA = (a - c).angle();
            float ccwB = wrapTau((b - c).angle() - angA);
            float ccwD = wrapTau((d - c).angle() - angA);
            if (ccwD <= ccwB) { a0 = angA; sweep = ccwB; }                    // CCW: A -> B through D
            else              { a0 = (b - c).angle(); sweep = TAU - ccwB; }   // CW: the complementary half
        }

        SKind type() const override { return SKind::Arc; }
        const char* kind() const override { return "arc"; }
        std::unique_ptr<Stoicheion> clone() const override { return std::make_unique<Arc2>(*this); }
        Json data() const override {
            return Json{ { "c", c.toJson() }, { "a", a.toJson() }, { "b", b.toJson() }, { "d", d.toJson() } };
        }
        void load(const Json& j) override {
            c = Pos2::fromJson(j.value("c", Json::object()));
            a = Pos2::fromJson(j.value("a", Json::object()));
            b = Pos2::fromJson(j.value("b", Json::object()));
            d = Pos2::fromJson(j.value("d", Json::object()));
        }

        void tessellate(std::vector<Pos>& out, bool skipFirst = false) const override {
            float r = radius(); if (r <= 0.0f) { return; }
            float a0, sweep; range(a0, sweep);
            sampleArc(out, c, r, a0, sweep, spanSteps(sweep), skipFirst);
        }
        // Exact: integral of (x dy - y dx) over the arc A -> B through D. Closed form
        // is c.x*(b.y-a.y) - c.y*(b.x-a.x) + r^2 * dtheta, where dtheta is the SIGNED
        // travel sweep (+ when chirality is CCW, - when CW).
        float signedAreaTerm() const override {
            float r = radius();
            float angA = (a - c).angle();
            float dtheta = (chirality() > 0) ? wrapTau((b - c).angle() - angA)
                                             : -wrapTau(angA - (b - c).angle());
            return c.x * (b.y - a.y) - c.y * (b.x - a.x) + r * r * dtheta;
        }
        float distanceTo(Pos p) const override {
            float r = radius();
            if (r <= 0.0f) { return (p - a).pythag(); }
            float a0, sweep; range(a0, sweep);
            float aP = wrapTau((p - c).angle() - a0);
            if (aP <= sweep) { return std::fabs((p - c).pythag() - r); }
            return std::min((p - a).pythag(), (p - b).pythag());
        }
        // Anchors mirror controlPoints (the solver and drag rely on the order):
        // centre, the two endpoints, and the through-point d.
        void anchors(std::vector<Pos>& out) const override {
            out.push_back(c); out.push_back(a); out.push_back(b); out.push_back(d);
        }
        void snapPoints(std::vector<Pos>& out) const override {
            out.push_back(c); out.push_back(a); out.push_back(b);
            float r = radius();
            if (r > 0.0f) { float a0, sweep; range(a0, sweep); out.push_back(c + Pos::fromAngle(a0 + sweep * 0.5f) * r); }
        }
        bool asCircle(Pos& cc, float& rr, float& a0, float& sweep) const override {
            float r = radius();
            if (r <= 0.0f) { return false; }
            cc = c; rr = r; range(a0, sweep); return true;
        }
        bool footOnCurve(Pos p, Pos& out) const override {
            float r = radius(); if (r <= 0.0f) { return false; }
            float a0, sweep; range(a0, sweep);
            float aP = wrapTau((p - c).angle() - a0);
            if (aP > sweep) { return false; }
            Pos dir = p - c; float dd = dir.pythag();
            out = (dd < 1e-6f) ? Pos(a) : (c + dir / dd * r);
            return true;
        }
        void controlPoints(std::vector<Pos2*>& out) override {
            out.push_back(&c); out.push_back(&a); out.push_back(&b); out.push_back(&d);
        }

        // Concentric offset, chirality-aware. The chirality is the travel direction
        // A -> B through D: CCW (+1) behaves like a CCW boundary, so a negative amount
        // shrinks it; CW (-1) reverses the sense, so the same amount grows it. Driving
        // the radius through zero flips the endpoints to their antipodes AND flips
        // the chirality: the through-point keeps its ORIGINAL bearing (|nr|), landing
        // on the complementary side of the chord -- the inverted arc travels the
        // other way round, exactly as a segment spline of the arc would.
        void offset(float amount) override {
            float rad = radius();
            if (rad < 1e-9f) { return; }
            float nr = rad + amount * static_cast<float>(chirality());
            a = c + (a - c) / rad * nr;
            b = c + (b - c) / rad * nr;
            float ld = (d - c).pythag();
            if (ld > 1e-9f) { d = c + (d - c) / ld * std::fabs(nr); }   // |nr|: chirality flips through zero
        }

        // d is the chirality marker, not a solved point, and is kept canonical:
        // exactly the midpoint of the half it currently indicates, on the circle.
        // Re-seated after any move/solve so it stays "at the proper location".
        int chiralitySlot() const override { return 3; }
        void normalize() override {
            float r = radius();
            if (r < 1e-9f) { return; }
            float a0, sweep; range(a0, sweep);
            d = c + Pos::fromAngle(a0 + sweep * 0.5f) * r;
        }
    };

    // Ellipse (centre + conjugate semi-axes U, V; only the centre is an absolute,
    // movable point -- translating the entity must not disturb the axis vectors)
    //--------------------------------------------------
    struct Ellipse2 : public Stoicheion {

        Pos2 c, u, v;

        Ellipse2() = default;
        Ellipse2(const Pos& c, const Pos& u, const Pos& v) : c(c), u(u), v(v) {}

        SKind type() const override { return SKind::Ellipse; }
        const char* kind() const override { return "ellipse"; }
        std::unique_ptr<Stoicheion> clone() const override { return std::make_unique<Ellipse2>(*this); }
        Json data() const override { return Json{ { "c", c.toJson() }, { "u", u.toJson() }, { "v", v.toJson() } }; }
        void load(const Json& j) override {
            c = Pos2::fromJson(j.value("c", Json::object()));
            u = Pos2::fromJson(j.value("u", Json::object()));
            v = Pos2::fromJson(j.value("v", Json::object()));
        }

        void tessellate(std::vector<Pos>& out, bool skipFirst = false) const override { sampleEllipse(out, c, u, v, 0.0f, TAU, spanSteps(TAU), skipFirst); }
        // Exact: a full ellipse contributes (u x v) * TAU (twice its area, signed by
        // handedness); the centre terms vanish since start == end.
        float signedAreaTerm() const override { return u.cross(v) * TAU; }
        float distanceTo(Pos p) const override {
            Pos foot; footNearest(p, 0.0f, TAU, foot); return (p - foot).pythag();
        }
        void anchors(std::vector<Pos>& out) const override { out.push_back(c); }
        void snapPoints(std::vector<Pos>& out) const override {
            out.push_back(c);                            // centre + the four axis ends
            out.push_back(c + u); out.push_back(c - u);
            out.push_back(c + v); out.push_back(c - v);
        }
        bool footOnCurve(Pos p, Pos& out) const override { footNearest(p, 0.0f, TAU, out); return true; }
        void controlPoints(std::vector<Pos2*>& out) override { out.push_back(&c); }

        // Rough in-place offset (kept for previews/quick edits): grow each conjugate
        // semi-axis, scaled by handedness sign(u x v). This only matches the true
        // offset at the four vertices -- offsetInto() is the analytically correct one.
        void offset(float amount) override {
            float d = amount * (u.cross(v) >= 0.0f ? 1.0f : -1.0f);
            float lu = u.pythag(), lv = v.pythag();
            if (lu > 1e-9f) { u += u / lu * d; }
            if (lv > 1e-9f) { v += v / lv * d; }
        }

        // The analytically-equidistant offset: every point pushed along its real
        // normal by `amount` (handedness-scaled: CCW insets, CW grows), emitted as a
        // closed polyline -- the exact parallel curve, which is not an ellipse.
        void offsetInto(float amount, std::vector<std::unique_ptr<Stoicheion>>& out) const override {
            float d = amount * (u.cross(v) >= 0.0f ? 1.0f : -1.0f);
            int steps = spanSteps(TAU);
            Pos prev = ellipseOffsetPoint(c, u, v, 0.0f, d);
            for (int i = 1; i <= steps; i++) {
                float t = TAU * (static_cast<float>(i) / static_cast<float>(steps));
                Pos cur = ellipseOffsetPoint(c, u, v, t, d);
                out.push_back(std::make_unique<Segment2>(prev, cur));
                prev = cur;
            }
        }

      protected:
        // Coarse nearest-point search over a parameter span (no closed form).
        void footNearest(Pos p, float t0, float span, Pos& out) const {
            const int N = 96; float best = 1e30f;
            for (int i = 0; i <= N; i++) {
                Pos e = ellipsePointAt(c, u, v, t0 + span * (static_cast<float>(i) / N));
                float d = (p - e).pythag();
                if (d < best) { best = d; out = e; }
            }
        }
    };

    // Elliptical arc (swept by increasing parameter a0 -> a1, i.e. CCW)
    //--------------------------------------------------
    struct EllipseArc2 : public Stoicheion {

        Pos2 c, u, v;
        float a0 = 0.0f, a1 = 0.0f;

        EllipseArc2() = default;
        EllipseArc2(const Pos& c, const Pos& u, const Pos& v, float a0, float a1)
            : c(c), u(u), v(v), a0(a0), a1(a1) {}

        // Construct the elliptical quarter-arc through a and b, tangent to tanA at a
        // and tanB at b. The centre O is placed so that (a-O) is parallel to tanB and
        // (b-O) to tanA -- i.e. a and b are the ends of *conjugate* semi-diameters, so
        // U = a-O and V = b-O are conjugate axes and the arc from parameter 0 (a) to
        // pi/2 (b) is tangent to both edges. Unlike a circular arc this is unique only
        // up to the conjugate-diameter choice; this is the natural canonical one.
        static EllipseArc2 Tangent(Pos a, Pos b, Pos tanA, Pos tanB) {
            float rxs = tanB.cross(tanA);
            Pos o = (std::fabs(rxs) > 1e-9f) ? (a + tanB * ((b - a).cross(tanA) / rxs))
                                             : ((a + b) * 0.5f);   // parallel tangents: fallback
            return EllipseArc2(o, a - o, b - o, 0.0f, TAU * 0.25f);
        }

        float span() const { return wrapTau(a1 - a0); }
        Pos startPoint() const { return ellipsePointAt(c, u, v, a0); }
        Pos endPoint() const { return ellipsePointAt(c, u, v, a0 + span()); }

        SKind type() const override { return SKind::EllipseArc; }
        const char* kind() const override { return "ellipseArc"; }
        std::unique_ptr<Stoicheion> clone() const override { return std::make_unique<EllipseArc2>(*this); }
        Json data() const override {
            return Json{ { "c", c.toJson() }, { "u", u.toJson() }, { "v", v.toJson() }, { "a0", a0 }, { "a1", a1 } };
        }
        void load(const Json& j) override {
            c = Pos2::fromJson(j.value("c", Json::object()));
            u = Pos2::fromJson(j.value("u", Json::object()));
            v = Pos2::fromJson(j.value("v", Json::object()));
            a0 = j.value("a0", 0.0f); a1 = j.value("a1", 0.0f);
        }

        void tessellate(std::vector<Pos>& out, bool skipFirst = false) const override {
            float sp = span(); sampleEllipse(out, c, u, v, a0, sp, spanSteps(sp), skipFirst);
        }
        // Exact: integral of (x dy - y dx) over the elliptical arc, swept by parameter
        // a0 -> a0+span(). Closed form is c.x*(p1.y-p0.y) - c.y*(p1.x-p0.x) plus the
        // central (u x v)*dt term; the handedness of u,v carries the chirality sign.
        float signedAreaTerm() const override {
            Pos p0 = startPoint(), p1 = endPoint();
            return c.x * (p1.y - p0.y) - c.y * (p1.x - p0.x) + u.cross(v) * span();
        }
        float distanceTo(Pos p) const override {
            Pos foot; footNearest(p, foot); return (p - foot).pythag();
        }
        void anchors(std::vector<Pos>& out) const override {
            out.push_back(c); out.push_back(startPoint()); out.push_back(endPoint());
        }
        bool footOnCurve(Pos p, Pos& out) const override { footNearest(p, out); return true; }
        void controlPoints(std::vector<Pos2*>& out) override { out.push_back(&c); }

        // Rough in-place offset (kept for previews): grow/shrink the semi-axes,
        // handedness-scaled. offsetInto() is the analytically correct one.
        void offset(float amount) override {
            float d = amount * (u.cross(v) >= 0.0f ? 1.0f : -1.0f);
            float lu = u.pythag(), lv = v.pythag();
            if (lu > 1e-9f) { u += u / lu * d; }
            if (lv > 1e-9f) { v += v / lv * d; }
        }

        // The analytically-equidistant offset over the arc's parameter span, emitted
        // as a polyline -- the exact parallel curve (handedness-scaled: CCW insets).
        void offsetInto(float amount, std::vector<std::unique_ptr<Stoicheion>>& out) const override {
            float d = amount * (u.cross(v) >= 0.0f ? 1.0f : -1.0f);
            float sp = span();
            int steps = spanSteps(sp);
            Pos prev = ellipseOffsetPoint(c, u, v, a0, d);
            for (int i = 1; i <= steps; i++) {
                float t = a0 + sp * (static_cast<float>(i) / static_cast<float>(steps));
                Pos cur = ellipseOffsetPoint(c, u, v, t, d);
                out.push_back(std::make_unique<Segment2>(prev, cur));
                prev = cur;
            }
        }

      protected:
        void footNearest(Pos p, Pos& out) const {
            const int N = 96; float sp = span(), best = 1e30f;
            for (int i = 0; i <= N; i++) {
                Pos e = ellipsePointAt(c, u, v, a0 + sp * (static_cast<float>(i) / N));
                float d = (p - e).pythag();
                if (d < best) { best = d; out = e; }
            }
        }
    };

    // Offset ellipse -- the EXACT parallel curve of an ellipse: every point of the
    // base ellipse displaced distance `d` along its true outward normal. This is
    // the curve an ellipse offsets to (it is NOT an ellipse, nor any Bezier --
    // offsets of non-Pythagorean-hodograph curves are non-rational), represented
    // exactly by its defining data: base (c, u, v), distance d, and the parameter
    // span. Travel runs a0 -> a1 in PARAMETER, either direction -- so reversal is
    // swapping a0/a1 and splitting is interval arithmetic, both exact. The family
    // is closed under offsetting: offsetting by delta just changes d.
    struct OffsetEllipse2 : public Stoicheion {

        Pos2 c, u, v;                      // the base ellipse (conjugate semi-axes)
        float d = 0.0f;                    // outward offset distance (signed)
        float a0 = 0.0f, a1 = 0.0f;        // travel: parameter a0 -> a1

        OffsetEllipse2() = default;
        OffsetEllipse2(const Pos& c, const Pos& u, const Pos& v, float d, float a0, float a1)
            : c(c), u(u), v(v), d(d), a0(a0), a1(a1) {}

        // Normalised travel coordinate s in [0,1] -> parameter / point / tangent.
        float paramAt(float s) const { return a0 + (a1 - a0) * s; }
        Pos pointAt(float s) const { return ellipseOffsetPoint(c, u, v, paramAt(s), d); }

        // Unit travel tangent by guarded central difference. (The analytic form
        // involves the unit-normal derivative; the difference is accurate to
        // ~1e-7 of scale and degrades gracefully at offset cusps, where the true
        // tangent is genuinely discontinuous anyway.)
        Pos tangentAt(float s) const {
            const float h = 1e-3f;
            Pos a = pointAt(std::max(0.0f, s - h));
            Pos b = pointAt(std::min(1.0f, s + h));
            Pos t = b - a; float l = t.pythag();
            return (l > 1e-12f) ? t / l : Pos(1.0f, 0.0f);
        }

        bool fullLoop() const { return std::fabs(std::fabs(a1 - a0) - TAU) < 1e-5f; }

        SKind type() const override { return SKind::OffsetEllipse; }
        const char* kind() const override { return "offsetEllipse"; }
        std::unique_ptr<Stoicheion> clone() const override { return std::make_unique<OffsetEllipse2>(*this); }
        Json data() const override {
            return Json{ { "c", c.toJson() }, { "u", u.toJson() }, { "v", v.toJson() },
                         { "d", d }, { "a0", a0 }, { "a1", a1 } };
        }
        void load(const Json& j) override {
            c = Pos2::fromJson(j.value("c", Json::object()));
            u = Pos2::fromJson(j.value("u", Json::object()));
            v = Pos2::fromJson(j.value("v", Json::object()));
            d = j.value("d", 0.0f); a0 = j.value("a0", 0.0f); a1 = j.value("a1", 0.0f);
        }

        void tessellate(std::vector<Pos>& out, bool skipFirst = false) const override {
            int steps = spanSteps(std::fabs(a1 - a0));
            for (int i = skipFirst ? 1 : 0; i <= steps; i++) { out.push_back(pointAt(static_cast<float>(i) / steps)); }
        }

        // Quadrature of ∮(x dy - y dx): the chord shoelace over a fine sampling of
        // the exact curve (no elementary closed form exists for the parallel
        // curve's area integral -- it is elliptic).
        float signedAreaTerm() const override {
            const int N = 256;
            float sum = 0.0f;
            Pos prev = pointAt(0.0f);
            for (int i = 1; i <= N; i++) {
                Pos p = pointAt(static_cast<float>(i) / N);
                sum += prev.x * p.y - prev.y * p.x;
                prev = p;
            }
            return sum;
        }

        float distanceTo(Pos p) const override {
            Pos foot; footNearest(p, foot); return (p - foot).pythag();
        }
        void anchors(std::vector<Pos>& out) const override {
            out.push_back(pointAt(0.0f)); out.push_back(pointAt(1.0f));
        }
        bool footOnCurve(Pos p, Pos& out) const override { footNearest(p, out); return true; }
        void controlPoints(std::vector<Pos2*>& out) override {}   // derived geometry: nothing solvable

        // Stoicheion::offset semantics (+grows): outward is simply more d.
        void offset(float amount) override { d += amount; }

      protected:
        void footNearest(Pos p, Pos& out) const {
            const int N = 96; float best = 1e30f;
            for (int i = 0; i <= N; i++) {
                Pos e = pointAt(static_cast<float>(i) / N);
                float dist = (p - e).pythag();
                if (dist < best) { best = dist; out = e; }
            }
        }
    };

    // Spline -- a (possibly offset) NURBS / B-spline curve, the catch-all for every
    // freeform edge OpenCASCADE hands back from a planar BRep section (torus and
    // freeform-surface sections, and conics when the section is approximated). Like
    // OffsetEllipse2 it is a parametric curve handled entirely NUMERICALLY, and it
    // carries the same offset trick: the true offset of a spline is not a spline, so
    // it is represented as base curve + signed distance `d`, evaluated as
    // base(u) + d * (left normal). The family is closed under offsetting -- offsetting
    // by delta just changes d -- exactly as for OffsetEllipse2.
    //
    // The base curve is its OWN control points / knots / weights (no OpenCASCADE
    // handle): Geo.Geometry stays OCC-free, and the De Boor evaluator below is the
    // whole dependency. Travel runs a0 -> a1 in the base KNOT parameter, either
    // direction -- so reversal is swapping a0/a1 and splitting is interval arithmetic
    // on [a0,a1], both leaving the control net untouched. A non-rational spline omits
    // weights (treated as all 1); a rational one (OCC's exact conics) carries them.
    struct Spline2 : public Stoicheion {

        std::vector<Pos2>  ctrl;           // control points
        std::vector<float> knots;          // flat knot vector (size = ctrl.size() + degree + 1)
        std::vector<float> weights;        // rational weights (empty => non-rational); size = ctrl.size()
        int   degree = 3;
        float d  = 0.0f;                   // outward offset distance (signed)
        float a0 = 0.0f, a1 = 0.0f;        // travel: base parameter a0 -> a1

        Spline2() = default;
        Spline2(std::vector<Pos2> ctrl, std::vector<float> knots, std::vector<float> weights, int degree)
            : ctrl(std::move(ctrl)), knots(std::move(knots)), weights(std::move(weights)), degree(degree) {
            a0 = domainMin(); a1 = domainMax();
        }
        Spline2(std::vector<Pos2> ctrl, std::vector<float> knots, std::vector<float> weights,
                int degree, float d, float a0, float a1)
            : ctrl(std::move(ctrl)), knots(std::move(knots)), weights(std::move(weights)),
              degree(degree), d(d), a0(a0), a1(a1) {}

        // The clamped-knot domain [knots[degree], knots[numCtrl]] the curve is defined over.
        float domainMin() const { return (static_cast<int>(knots.size()) > degree) ? knots[degree] : 0.0f; }
        float domainMax() const {
            const int n = static_cast<int>(ctrl.size());
            return (static_cast<int>(knots.size()) > n) ? knots[n] : 1.0f;
        }

        // The knot span k with knots[k] <= u < knots[k+1], clamped to [degree, numCtrl-1]
        // (Piegl & Tiller, "The NURBS Book").
        int findSpan(float u) const {
            const int n = static_cast<int>(ctrl.size()) - 1;
            if (u >= knots[n + 1]) { return n; }
            if (u <= knots[degree]) { return degree; }
            int low = degree, high = n + 1, mid = (low + high) / 2;
            while (u < knots[mid] || u >= knots[mid + 1]) {
                if (u < knots[mid]) { high = mid; } else { low = mid; }
                mid = (low + high) / 2;
            }
            return mid;
        }

        // De Boor evaluation of the BASE (un-offset) curve at knot parameter u. Works
        // in homogeneous coordinates so rational (NURBS) curves evaluate exactly.
        Pos basePointAt(float u) const {
            const int n = static_cast<int>(ctrl.size());
            if (n == 0) { return Pos(); }
            if (n == 1) { return ctrl[0]; }
            const int p = degree;
            if (p < 1 || p > 24 || static_cast<int>(knots.size()) < n + p + 1) { return ctrl[0]; }
            u = std::clamp(u, knots[p], knots[n]);
            const int k = findSpan(u);
            const bool rational = (weights.size() == ctrl.size());
            float dx[26], dy[26], dw[26];
            for (int j = 0; j <= p; j++) {
                const int idx = k - p + j;
                const float w = rational ? weights[idx] : 1.0f;
                dx[j] = ctrl[idx].x * w; dy[j] = ctrl[idx].y * w; dw[j] = w;
            }
            for (int r = 1; r <= p; r++) {
                for (int j = p; j >= r; j--) {
                    const int idx = k - p + j;
                    const float left = knots[idx], right = knots[idx + p - r + 1];
                    const float den = right - left;
                    const float a = (den > 1e-12f) ? (u - left) / den : 0.0f;
                    dx[j] = (1.0f - a) * dx[j - 1] + a * dx[j];
                    dy[j] = (1.0f - a) * dy[j - 1] + a * dy[j];
                    dw[j] = (1.0f - a) * dw[j - 1] + a * dw[j];
                }
            }
            const float w = dw[p];
            return (std::fabs(w) > 1e-12f) ? Pos(dx[p] / w, dy[p] / w) : Pos(dx[p], dy[p]);
        }

        // Unit base tangent by guarded central difference (degrees of the rational
        // derivative aren't worth carrying; the difference is ample at section scale).
        Pos baseTangentAt(float u) const {
            const float lo = domainMin(), hi = domainMax();
            const float h = std::max(1e-4f, (hi - lo) * 1e-3f);
            Pos t = basePointAt(std::min(hi, u + h)) - basePointAt(std::max(lo, u - h));
            float l = t.pythag();
            return (l > 1e-12f) ? t / l : Pos(1.0f, 0.0f);
        }

        // The realized (offset) point at base parameter u: base displaced by d along
        // the left normal of the base tangent (no displacement when d == 0).
        Pos offsetPointAt(float u) const {
            Pos base = basePointAt(u);
            if (std::fabs(d) < 1e-12f) { return base; }
            Pos t = baseTangentAt(u);
            return base + Pos(-t.y, t.x) * d;
        }

        // Normalised travel coordinate s in [0,1] -> base parameter / point / tangent.
        float paramAt(float s) const { return a0 + (a1 - a0) * s; }
        Pos   pointAt(float s) const { return offsetPointAt(paramAt(s)); }
        Pos tangentAt(float s) const {
            const float h = 1e-3f;
            Pos t = pointAt(std::min(1.0f, s + h)) - pointAt(std::max(0.0f, s - h));
            float l = t.pythag();
            return (l > 1e-12f) ? t / l : Pos(1.0f, 0.0f);
        }

        // A closed curve whose travel covers the whole domain: no endpoints (like a
        // full circle / full-span offset ellipse).
        bool fullLoop() const {
            if (ctrl.size() < 2) { return false; }
            const float lo = domainMin(), hi = domainMax();
            if (std::fabs(std::fabs(a1 - a0) - (hi - lo)) > 1e-4f * std::max(1.0f, hi - lo)) { return false; }
            return (basePointAt(lo) - basePointAt(hi)).pythag() < 1e-4f;
        }

        SKind type() const override { return SKind::Spline; }
        const char* kind() const override { return "spline"; }
        std::unique_ptr<Stoicheion> clone() const override { return std::make_unique<Spline2>(*this); }
        Json data() const override {
            Json cps = Json::array(); for (const Pos2& p : ctrl)    { cps.push_back(p.toJson()); }
            Json kns = Json::array(); for (float k : knots)         { kns.push_back(k); }
            Json wts = Json::array(); for (float w : weights)       { wts.push_back(w); }
            return Json{ { "ctrl", cps }, { "knots", kns }, { "weights", wts },
                         { "degree", degree }, { "d", d }, { "a0", a0 }, { "a1", a1 } };
        }
        void load(const Json& j) override {
            ctrl.clear(); knots.clear(); weights.clear();
            if (auto it = j.find("ctrl");    it != j.end() && it->is_array()) { for (const Json& p : *it) { ctrl.push_back(Pos2::fromJson(p)); } }
            if (auto it = j.find("knots");   it != j.end() && it->is_array()) { for (const Json& k : *it) { knots.push_back(k.get<float>()); } }
            if (auto it = j.find("weights"); it != j.end() && it->is_array()) { for (const Json& w : *it) { weights.push_back(w.get<float>()); } }
            degree = j.value("degree", 3);
            d = j.value("d", 0.0f); a0 = j.value("a0", 0.0f); a1 = j.value("a1", 0.0f);
        }

        void tessellate(std::vector<Pos>& out, bool skipFirst = false) const override {
            int steps = sampleSteps();
            for (int i = skipFirst ? 1 : 0; i <= steps; i++) { out.push_back(pointAt(static_cast<float>(i) / steps)); }
        }
        // Quadrature of ∮(x dy - y dx) over the travel span -- the chord shoelace over a
        // fine sampling of the exact curve (no elementary closed form for a spline's
        // area integral).
        float signedAreaTerm() const override {
            const int N = 256;
            float sum = 0.0f;
            Pos prev = pointAt(0.0f);
            for (int i = 1; i <= N; i++) {
                Pos p = pointAt(static_cast<float>(i) / N);
                sum += prev.x * p.y - prev.y * p.x;
                prev = p;
            }
            return sum;
        }
        float distanceTo(Pos p) const override { Pos foot; footNearest(p, foot); return (p - foot).pythag(); }
        void anchors(std::vector<Pos>& out) const override {
            out.push_back(pointAt(0.0f)); out.push_back(pointAt(1.0f));
        }
        bool footOnCurve(Pos p, Pos& out) const override { footNearest(p, out); return true; }
        void controlPoints(std::vector<Pos2*>& out) override {}   // derived geometry: nothing solvable

        // Stoicheion::offset semantics (+grows): outward is simply more d. The family
        // is closed under offsetting, so offsetBy()/offsetInto() default forms stay
        // exact (one Spline2 with a larger d) -- no tessellation needed.
        void offset(float amount) override { d += amount; }

      protected:
        // Tessellation density from a coarse arc-length estimate (chord target ~0.25
        // model units, matching the section sampler); bounded so trims and giants behave.
        int sampleSteps() const {
            const int C = 16; float len = 0.0f; Pos prev = pointAt(0.0f);
            for (int i = 1; i <= C; i++) { Pos p = pointAt(static_cast<float>(i) / C); len += (p - prev).pythag(); prev = p; }
            return std::clamp(static_cast<int>(std::ceil(len / 0.25f)), 24, 256);
        }
        void footNearest(Pos p, Pos& out) const {
            const int N = 96; float best = 1e30f;
            for (int i = 0; i <= N; i++) {
                Pos e = pointAt(static_cast<float>(i) / N);
                float dist = (p - e).pythag();
                if (dist < best) { best = dist; out = e; }
            }
        }
    };

    // Factory: reconstruct an entity from its JSON envelope.
    //--------------------------------------------------
    inline std::unique_ptr<Stoicheion> stoicheionFromJson(const Json& j) {

        std::string kind = j.value("kind", std::string());
        std::unique_ptr<Stoicheion> e;

        if      (kind == "point")      { e = std::make_unique<Point2>(); }
        else if (kind == "segment")    { e = std::make_unique<Segment2>(); }
        else if (kind == "circle")     { e = std::make_unique<Circle2>(); }
        else if (kind == "arc")        { e = std::make_unique<Arc2>(); }
        else if (kind == "ellipse")    { e = std::make_unique<Ellipse2>(); }
        else if (kind == "ellipseArc") { e = std::make_unique<EllipseArc2>(); }
        else if (kind == "offsetEllipse") { e = std::make_unique<OffsetEllipse2>(); }
        else if (kind == "spline")     { e = std::make_unique<Spline2>(); }
        else { return nullptr; }

        e->load(j);
        e->id = idFromHex(j.value("id", std::string()));
        if (e->id == 0) { e->id = newId(); }
        e->construction = j.value("construction", false);
        e->locked = j.value("locked", false);
        e->openAir = j.value("openAir", false);
        e->group = static_cast<DisplayGroup>(j.value("group", 0));
        e->groupLevel = j.value("groupLevel", 0);
        return e;
    }

    // The full contents of a sketch: the polymorphic list of geometry entities.
    // Constraint relations between entities are intentionally out of scope here --
    // sketch constraint solving is handled externally; this is purely the entity
    // document (the committed geometry, or a tool's transient preview).
    struct SketchGeometry {

        std::vector<std::unique_ptr<Stoicheion>> entities;

        SketchGeometry() = default;
        SketchGeometry(SketchGeometry&&) = default;
        SketchGeometry& operator=(SketchGeometry&&) = default;
        SketchGeometry(const SketchGeometry&) = delete;
        SketchGeometry& operator=(const SketchGeometry&) = delete;

        void clear() { entities.clear(); }
        bool empty() const { return entities.empty(); }

        // Resolve an entity by its stable id (linear for now; a cached id->pointer
        // map comes later).
        Stoicheion* byId(Id id) const {
            for (const auto& e : entities) { if (e && e->id == id) { return e.get(); } }
            return nullptr;
        }

        // Re-seat every entity's chirality marker (circle b / arc d) to its canonical
        // spot. The markers ride along with moves but are never independent points, so
        // this is run after any bulk edit that may have shifted control points.
        void normalize() {
            for (const auto& e : entities) { if (e) { e->normalize(); } }
        }

        // Append an entity, assigning it a stable id if it doesn't have one yet.
        // Returns the new index (so a tool can keep extending it).
        size_t add(std::unique_ptr<Stoicheion> e) {
            if (e->id == 0) { e->id = newId(); }
            entities.push_back(std::move(e));
            return entities.size() - 1;
        }

        Json toJson() const {
            Json out;
            out["entities"] = Json::array();
            for (const auto& e : entities) { if (e) { out["entities"].push_back(e->toJson()); } }
            return out;
        }

        static SketchGeometry fromJson(const Json& j) {
            SketchGeometry g;
            if (auto it = j.find("entities"); it != j.end() && it->is_array()) {
                for (const Json& e : *it) {
                    if (auto ent = stoicheionFromJson(e)) { g.entities.push_back(std::move(ent)); }
                }
            }
            return g;
        }
    };
}
