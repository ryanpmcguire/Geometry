# The Frontier — extracting one chain's true offset

## What a frontier is

A **frontier** is the boundary between material and void. When we offset a chain and ask "what is
the new boundary of the region?", we are asking for the frontier of the offset. The difficulty is
that the blind offset ([02-Offsetting-Semantics.md](02-Offsetting-Semantics.md)) does not hand us
a clean boundary — it hands us a **self-intersecting tangle**: the true frontier, plus excursion
loops at converging corners, plus wherever the offset drove distant parts of the chain across
each other. Extraction is the act of recovering the true frontier from that tangle.

The key realisation is that the tangle is *not* ambiguous. The self-intersections partition the
plane into faces, and there is an exact, integer-valued field over those faces that says how
"deep" each one is. The true frontier is the boundary of an **extremal level-set** of that field.
Everything reduces to two questions: *how is the field built*, and *which extreme do we take*.

## The crossing-number field

Build the field by a single **travel-ordered walk** of the raw offset. Split the offset at every
self-crossing into sub-edges, keep them in travel order, and walk them, carrying a running
integer **level**. At each crossing the walk passes, look at how the *other* strand crosses this
one: a strand passing from our right to our left counts `+1`, left to right counts `−1`. Each
sub-edge wears the level in force as it is traversed. Fragments that are consecutive and share a
level belong to the same piece of frontier; the level *changes* exactly at a crossing.

Two properties of this field matter for what follows:

- **Its orientation is set by handedness.** The `+1`/`−1` at each crossing is computed from
  travel direction (`ourTangent × otherTangent`). Reverse the chain's handedness and you reverse
  every travel tangent, which flips the sign of the whole field. So the field is not
  handedness-neutral; "high" and "low" are defined *relative to travel*.
- **Only differences are meaningful.** The walk's absolute starting value is arbitrary (it depends
  where the closed walk began), so the field carries a free additive constant. Only the *relative*
  levels — which faces are deeper than which — are real. The implementation removes the constant by
  normalising one extreme to zero; that normalisation is a *choice of which extreme to anchor on*,
  and — this is the trap — the current code anchors on the maximum because it only ever needed the
  minimum. See below.

## Which extreme: the side decides

*(Corrected 2026-07-11. An earlier revision of this section claimed the extreme was selected by
the sense `σ = sign(amount · h)`. That is impossible, and the error surfaced in production as
avoid boundaries failing to bound. The record of why is kept below, because the argument is the
lesson.)*

The true frontier is the level-set at an **extreme** of the field — and which extreme is set by
**`sign(amount)` alone**: the *side* the curve was displaced, in the same travel-left orientation
the field itself is built from.

- **Leftward offset (`amount > 0`)** — the frontier is the **minimum** level.
- **Rightward offset (`amount < 0`)** — the frontier is the **maximum** level.

**Why it cannot key on `σ`: the mirror argument.** Reflect the plane. A CCW chain inset by `+a`
(the production-proven case, frontier at the minimum) becomes a CW chain offset by `−a` — same
geometry, same true frontier fragments. But every increment in the walk is a 2D cross product,
and cross products **negate under reflection**: the whole level field negates, so the fragments
that sat at the minimum now sit at the maximum. `σ` is mirror-invariant (`(−a)·(−h) = a·h`), yet
the correct extreme label flipped. No rule reading `σ` alone can be right. `sign(amount)` flips
under the mirror exactly as the label does.

**The evidence anchors, one per side.** Leftward-min is proven by the entire legacy engine (every
chain, CCW outers and CW islands alike, offset at `amount > 0` and extracted at the minimum, in
production for months). Rightward-max is proven by the crescent-sliver seed on a real part: a CCW
fill grown (`amount < 0`) through its own tangency cusps — a genuinely self-intersecting raw
offset — extracted the correct envelope at the maximum.

The earlier intuition ("shrinking buries the real boundary at the deep end; growing pushes it to
the rim") remains true *geometrically* — but "deep end" is a fact about the plane, while
minimum-vs-maximum is a fact about the field's *labels*, and the labels flip with travel
direction. `σ` captures the geometry; `sign(amount)` captures the labels; selection needs the
labels. The sense `σ` still governs everything that is genuinely a *region* question — shrink
versus grow, and the inversion guard below.

## What the evidence actually covers — and does not

The legacy extraction takes the **minimum** level, unconditionally, and every chain it ever ran —
CCW outers shrinking *and* CW islands growing — was offset at `amount > 0`. Under the corrected
rule that is not a coincidence it survived; it is the rule: **leftward → minimum**, regardless of
handedness, regardless of shrink-or-grow. The legacy engine is a months-long production proof of
the left half of the rule.

The right half (`amount < 0` → maximum) is anchored by the crescent-sliver seed: a CCW fill grown
through its own tangency cusps — a genuinely self-crossing raw offset — extracting the correct
outer envelope at the maximum, verified on a real part.

The remaining under-evidenced corner: a **CW chain at `amount > 0` whose growth self-crosses**
(a keep-island with a concavity tighter than the offset). The rule says minimum; the legacy
engine agrees but ran that case mostly on flat fields (small stepovers rarely self-cross), so the
evidence is thin rather than absent.

- **Fixture obligation.** Construct the thin-corner case with a known answer (a concave CW island
  grown until the concavity closes; the mirror "U" outset as its CCW twin). Confirm the extremes
  match `sign(amount)` on both, and that the max-anchored normalisation is label-consistent. A
  *measured* result, written as a known-answer test — the mirror argument constrains the answer,
  but the fixture is what earns it.

## The inversion law for a single chain

Extraction produces candidate loops; one law decides which survive. A chain being **inset** can be
driven through zero — a feature narrower than twice the offset collapses and, by the positional
representation, re-emerges *inside-out* with its handedness flipped
([01-Stoicheia-And-Handedness.md](01-Stoicheia-And-Handedness.md)). That flipped candidate is an
**inversion**: it is not a real boundary, it is the geometry folded through itself.

The law: **a candidate's handedness may not flip away from its source's in the shrinking
direction.** Concretely, a CCW source (an inside) that produces a CW child (an outside) has
inverted — *discard it*. This is how recursion terminates honestly: a ring too small to inset
again simply produces nothing. Extinction, not error.

Crucially, this guard is governed by the **sense**, not by handedness alone. Inversion-through-
zero is possible only while a chain is **shrinking** (`σ > 0`) — that is the only direction in
which a dimension can be driven negative. A *growing* chain (`σ < 0`) cannot collapse, so the
guard must **not** fire on it; a growing configuration is *allowed* to change topology (it can
legitimately enclose new regions — see [04-Mitosis-And-Merging.md](04-Mitosis-And-Merging.md)).
The current code writes this guard as "source is CCW," which is again a frozen proxy: at the one
tested direction, `σ > 0` *is* "source is CCW." Generalised, the guard is "`σ > 0` and the child's
handedness flipped."

## Carried across to the reconstruction

- Extraction is: build the crossing-number field by a travel-ordered walk, take the frontier at an
  **extreme** level-set, apply the inversion guard.
- The extreme is chosen by **`sign(amount)`** — leftward → minimum, rightward → maximum. It cannot
  key on `σ` (the mirror argument above) and it cannot key on handedness alone.
- The inversion guard fires only when `σ > 0` (shrinking); it discards handedness-flipped children.
  Selection and the guard read **different variables** — conflating them was the original
  hard-coding, and mis-fusing them into `σ` was the reconstruction's one real error.
- Anchors: leftward-min (legacy engine, production), rightward-max (the crescent seed, real part).
  The least-evidenced corner is a **CW chain at `amount > 0` whose growth self-crosses** (min per
  the rule; legacy ran it mostly vacuously) — still worth a dedicated known-answer fixture.
- Preserve the numerical scar tissue verbatim: crossings computed in **double precision**;
  tangent-seam **fusion**; **small-feature pruning**; interiority measured in **absolute distance**
  (never a parameter tolerance, which scales with edge length); and **transversal-kiss rejection**
  (a tangential touch is not a crossing — its `±1` is pure float noise and must not enter the walk).
