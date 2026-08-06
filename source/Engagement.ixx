module;

#include <vector>
#include <cmath>
#include <algorithm>
#include <cstddef>

export module Geo.Engagement;

import Rev.Core.Pos;

export import Geo.Chain;

// =====================================================================
// Engagement -- TEST analytic engagement estimation.
//
// The engagement field answers, as a function of arc length along a toolpath
// chain: how much of the tool's swept circle meets material that is still
// standing? The Minkowski insight makes this a pure distance query, no
// boolean geometry required: a chain's cleared region is the Minkowski sum of
// the chain and the tool circle, so a point is CLEARED iff it lies within
// toolR of any already-travelled tool center. At arclength s the RIM of the
// tool circle is sampled; a rim point engages iff it is in material (not in
// the air region) and not cleared by (a) any previously cut chain or (b) this
// chain's own prefix [0, s). The prefix test is exclusive and strict
// (< R - eps), which needs no positional epsilon: earlier centers sit
// strictly closer than R to the trailing rim and strictly farther from the
// leading rim, so the trailing half reads cleared and the front reads intact
// by construction.
//
// This is deliberately a bounded ESTIMATOR: chains are tessellated once into
// polylines with bounding boxes, sampling is capped, and every query is an
// early-out point-to-segment sweep. Consumers treat the result as a field to
// derive feed influences from -- and can disable the whole mechanism.
// =====================================================================
export namespace Geo::Engagement {

    using Rev::Core::Pos;

    // A chain prepared for distance queries: directed polyline + cumulative
    // arclength + inflated bounding box.
    struct Swept {

        std::vector<Pos>   pts;
        std::vector<float> cum;      // arclength at each point
        float length = 0.0f;
        float x0 = 0.0f, y0 = 0.0f, x1 = 0.0f, y1 = 0.0f;

        static Swept of(const Chain& c) {
            Swept s;
            c.tessellateDirected(s.pts);
            s.cum.resize(s.pts.size(), 0.0f);
            if (!s.pts.empty()) {
                s.x0 = s.x1 = s.pts[0].x;
                s.y0 = s.y1 = s.pts[0].y;
            }
            for (std::size_t i = 0; i < s.pts.size(); i++) {
                if (i > 0) { s.cum[i] = s.cum[i - 1] + (s.pts[i] - s.pts[i - 1]).pythag(); }
                s.x0 = std::min(s.x0, s.pts[i].x); s.x1 = std::max(s.x1, s.pts[i].x);
                s.y0 = std::min(s.y0, s.pts[i].y); s.y1 = std::max(s.y1, s.pts[i].y);
            }
            s.length = s.cum.empty() ? 0.0f : s.cum.back();
            return s;
        }

        // Point at arclength (clamped)
        Pos at(float s) const {
            if (pts.empty()) { return {}; }
            if (s <= 0.0f) { return pts.front(); }
            if (s >= length) { return pts.back(); }
            std::size_t i = std::size_t(std::upper_bound(cum.begin(), cum.end(), s) - cum.begin());
            if (i < 1) { i = 1; }
            if (i >= pts.size()) { i = pts.size() - 1; }
            const float seg = cum[i] - cum[i - 1];
            const float t = seg > 1e-9f ? (s - cum[i - 1]) / seg : 0.0f;
            return pts[i - 1] + (pts[i] - pts[i - 1]) * t;
        }

        // Squared distance from q to the polyline, restricted to arclength < limit
        // (limit < 0 = the whole chain). `stop2` allows an early out: once a
        // distance at or below it is found, that value returns immediately.
        float dist2(Pos q, float limit = -1.0f, float stop2 = 0.0f) const {

            if (pts.empty()) { return 1e30f; }
            // bbox reject against the caller's interest radius (sqrt(best) work
            // saved by the conservative box test)
            float best = 1e30f;
            for (std::size_t i = 0; i + 1 < pts.size(); i++) {
                if (limit >= 0.0f && cum[i] >= limit) { break; }
                Pos a = pts[i];
                Pos b = pts[i + 1];
                float segLen = cum[i + 1] - cum[i];
                if (limit >= 0.0f && cum[i + 1] > limit && segLen > 1e-9f) {
                    b = a + (b - a) * ((limit - cum[i]) / segLen);   // clip the last segment at the limit
                }
                const Pos ab = b - a;
                const float d2ab = ab.dot(ab);
                float t = d2ab > 1e-12f ? (q - a).dot(ab) / d2ab : 0.0f;
                t = std::clamp(t, 0.0f, 1.0f);
                const Pos foot = a + ab * t;
                const float d2 = (q - foot).dot(q - foot);
                if (d2 < best) {
                    best = d2;
                    if (best <= stop2) { return best; }
                }
            }
            return best;
        }

        bool nearBox(Pos q, float r) const {
            return q.x >= x0 - r && q.x <= x1 + r && q.y >= y0 - r && q.y <= y1 + r;
        }
    };

    struct Sample { float s = 0.0f; float engagement = 0.0f; };

    struct Field {
        std::vector<Sample> samples;

        // clamped lerp over the samples
        float at(float s) const {
            if (samples.empty()) { return 0.0f; }
            if (s <= samples.front().s) { return samples.front().engagement; }
            if (s >= samples.back().s) { return samples.back().engagement; }
            std::size_t i = 1;
            while (i < samples.size() && samples[i].s < s) { i++; }
            const Sample& a = samples[i - 1];
            const Sample& b = samples[i];
            const float d = b.s - a.s;
            const float t = d > 1e-9f ? (s - a.s) / d : 0.0f;
            return a.engagement + (b.engagement - a.engagement) * t;
        }
    };

    // Estimate the engagement field of `path` cut by a tool of radius toolR.
    //   cleared  : previously cut chains (each fully swept by the tool)
    //   material : ORIENTED FILL boundary chains -- winding >= 1 means the point
    //              stands in material (this matches the slicer's fill-section
    //              convention; the region OUTSIDE is open air). Empty = everything
    //              is material (fail-safe: engagement reads high, never low).
    //   step     : sample spacing along the path (<= 0 -> toolR / 4)
    //   rim      : rim samples per circle
    // Engagement is the fraction of the FULL rim standing in uncleared material:
    // ~0.5 = full slot (the entire leading arc), ~0 = cutting air.
    inline Field estimate(const Chain& path,
                          const std::vector<const Chain*>& cleared,
                          const std::vector<const Chain*>& material,
                          float toolR, float step = 0.0f, int rim = 16) {

        Field field;
        if (toolR <= 1e-5f) { return field; }
        if (step <= 1e-6f) { step = 0.25f * toolR; }

        Swept self = Swept::of(path);
        if (self.length <= 1e-6f || self.pts.size() < 2) { return field; }

        std::vector<Swept> others;
        others.reserve(cleared.size());
        for (const Chain* c : cleared) {
            if (!c || c->edges.empty()) { continue; }
            Swept s = Swept::of(*c);
            if (!s.pts.empty()) { others.push_back(std::move(s)); }
        }

        // bounded work: never more than this many circle probes per chain
        constexpr int MaxSamples = 2048;
        int n = int(std::ceil(self.length / step));
        if (n > MaxSamples) { n = MaxSamples; step = self.length / float(n); }
        if (n < 1) { n = 1; }
        if (rim < 6) { rim = 6; }

        const float eps = 1e-4f;
        const float rIn2 = (toolR - eps) * (toolR - eps);
        constexpr float TAU = 6.28318530718f;

        field.samples.reserve(std::size_t(n) + 1);
        for (int i = 0; i <= n; i++) {

            const float s = std::min(self.length, float(i) * step);
            const Pos p = self.at(s);

            int engaged = 0;
            for (int k = 0; k < rim; k++) {
                const float a = TAU * (float(k) + 0.5f) / float(rim);
                const Pos q = p + Pos(std::cos(a), std::sin(a)) * toolR;

                // material test: winding >= 1 against the oriented fill = material
                // stands here; outside the fill is open air (nothing to cut)
                if (!material.empty()) {
                    int w = 0;
                    for (const Chain* mc : material) { if (mc) { w += mc->windingAround(q); } }
                    if (w < 1) { continue; }
                }

                // cleared by this path's own prefix [0, s)? strict < R:
                // trailing rim reads cleared, leading rim reads intact
                if (self.dist2(q, s, rIn2) < rIn2) { continue; }

                // cleared by an earlier chain's full sweep?
                bool clearedByOther = false;
                for (const Swept& o : others) {
                    if (!o.nearBox(q, toolR)) { continue; }
                    if (o.dist2(q, -1.0f, rIn2) < rIn2) { clearedByOther = true; break; }
                }
                if (clearedByOther) { continue; }

                engaged++;
            }

            field.samples.push_back({ s, float(engaged) / float(rim) });
        }

        return field;
    }
}
