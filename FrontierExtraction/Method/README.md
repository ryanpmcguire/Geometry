# Frontier Extraction — the Method

This folder is the **written contract** for a deliberate, from-scratch reconstruction of the
Geometry library's offset / crossing / merge machinery. The original implementation
(`Chain.ixx`, `offsetValid`, `mitose`) was written in *discovery* mode — we were finding out
whether the crossing-number + winding-frontier idea worked at all, and it does, wonderfully.
Because there was only ever one question being asked of it (*"what is the internal frontier of
an inward offset?"*), the answer to that one question got welded into the core as a constant
rather than exposed as a choice.

These documents state the method **generally** — the way it should be re-expressed once the
tacit understanding is made explicit. They are the target the parallel implementation writes
to, and the checklist the differential harness verifies against.

## The documents, in reading order

1. **[01-Stoicheia-And-Handedness.md](01-Stoicheia-And-Handedness.md)** — why geometry is
   represented as *points*, never bools or per-type sign variables; how each primitive carries
   its handedness *in its point construction*, and why that makes it respond to an offset by its
   own nature.
2. **[02-Offsetting-Semantics.md](02-Offsetting-Semantics.md)** — what offsetting *is*, and why
   a counter-clockwise and a clockwise chain respond oppositely to the *same* signed offset
   value — the emergent, chirality-driven behaviour the whole method rests on.
3. **[03-The-Frontier.md](03-The-Frontier.md)** — the concept of a *frontier*; how a single
   blindly-offset chain self-intersects, and how the true frontier is extracted as an **extremal
   level-set** of a crossing-number field — with the extreme (minimum *or* maximum) chosen by the
   offset's *sense relative to handedness*.
4. **[04-Mitosis-And-Merging.md](04-Mitosis-And-Merging.md)** — the same frontier idea across
   *many* chains: **mitosis** (one region splitting into several) and **merging** (regions
   combining), the winding-rule that selects the material region, and the handedness laws that
   tie it together — including why two insides can never create a new outside, but two outsides
   can create a new inside.
5. **[05-Opposite-Handed-Islands.md](05-Opposite-Handed-Islands.md)** — why an internal island
   *must* carry the opposite handedness to the profile containing it: the "material on the left"
   invariant forces it, it is what makes offset crossings resolve into consistently-handed chains,
   and it is why the method applies *one* offset value to *differently-handed* chains rather than
   different values to same-handed ones — the choice from which the whole algebra of chain identity
   emerges. Includes the mirror: neither absolute handedness is privileged, only the disagreement.
6. **[06-The-Seed-Construction.md](06-The-Seed-Construction.md)** — the synthesis. Given
   randomly-oriented **cut** and **avoid** systems, build the seed `(cut ⊕ R) ∖ (avoid ⊕ R)` — the
   true maximum extent of the tool centre — by orienting the avoid as the cut's island, growing both
   by the tool radius, and merging. Contains the hole discovery (a hole within material is an island,
   so its handedness flips, extending the parity rule across complementary systems) and the warning
   that the whole construction lives in the `σ < 0` regime and so depends on 03's unproven growing-
   direction extraction.

## The discipline that governs the reconstruction

Three rules, because a reconstruction of working, safety-critical geometry fails in specific
ways:

- **The existing code is the oracle.** For the entire proven envelope — inward offsets, the real
  slice corpus the strategies already run — the reconstruction must reproduce the current output
  chain-for-chain, within tolerance, *before anything relies on it*. Elegance earns nothing here.
- **The new regions have no oracle.** Outward offsets, the maximum-frontier selection, and the
  difference-merge are precisely the questions the old code never answered, so they cannot be
  checked against it. They are pinned with **constructed analytic fixtures** — cases whose exact
  answer is written by hand — never by "it looks right."
- **The scar tissue is content, not clutter.** The old code's double-precision crossing kernels,
  tangency fusion, small-feature pruning, absolute-distance interiority, and transversal-kiss
  rejection are each a real failure mode already paid for. "Beautify" must never mean "assume the
  edge cases away." They are catalogued here and carried across faithfully.

Where a document states a rule the current code *proves* by running, it says so. Where it states
a rule the current code merely *never contradicts* because the case is under-exercised, it says
that too, and marks it a fixture obligation. The distinction is the whole point.
