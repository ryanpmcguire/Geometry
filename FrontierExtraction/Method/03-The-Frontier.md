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

## Which extreme: the sense decides

The true frontier is the level-set at an **extreme** of the field — but *which* extreme, minimum
or maximum, is not fixed. It is set by the offset's **sense** relative to the chain's handedness,
`σ = sign(amount · h)` — the same scalar that decides shrink-versus-grow in
[02-Offsetting-Semantics.md](02-Offsetting-Semantics.md):

- When the chain's own region **shrinks** (`σ > 0` — a CCW inset, or a CW island being inset), the
  true frontier is driven to the **deep** end of the field. The excursion loops sit *outside* it,
  at shallower levels. The frontier is the **minimum** level.
- When the chain's own region **grows** (`σ < 0` — a CCW outset, or a CW island being grown), the
  true frontier is the **outer envelope**. The excursions now sit *inside* it. The frontier is the
  **opposite** extreme.

The intuition is direct: shrinking buries the real boundary at the bottom of the tangle and throws
the debris outward; growing pushes the real boundary to the outer rim and traps the debris inside.
The field's orientation flips with handedness *and* the "which extreme is real" flips with offset
direction — and because both effects are already fused into `σ`, **the single scalar `σ` selects
the extreme.** This is why the method can do minimum *or* maximum winding analysis and get the
right frontier either way: it is not two algorithms, it is one algorithm reading one sign.

> This corrects an earlier, looser statement that the selection was handedness-*invariant* and
> keyed on `sign(amount)` alone. It is not. It keys on `sign(amount · h)`. The two happen to
> coincide for CCW chains, which is why the confusion was easy — and why the CW case is where it
> would silently go wrong.

## What the current code actually proves — and does not

The existing extraction takes the **minimum** level, unconditionally, and it is correct — for the
cases it runs. Those cases are all `σ > 0`: the recursion only ever insets (positive amount), so
CCW outers shrink (`σ > 0`, minimum, correct) and CW islands, at that same positive amount, *grow*
(`σ < 0`) — but a *growing convex-ish island does not self-intersect*, so its field is flat, its
minimum equals its maximum, and taking the minimum happens to be right by vacuity.

The gap is exact and worth stating plainly:

> **A growing chain that genuinely self-intersects — a CW island with a concavity, grown until the
> concavity closes on itself; equivalently a CCW chain outset until a neck pinches — is the case
> `σ < 0` with a non-flat field. The current uniform-minimum rule extracts the wrong extreme
> there. It has never been caught because that case does not arise in the tested corpus.**

This is not a hypothetical for the reconstruction: the whole reason the reconstruction exists —
the outward-offset seed construction — deliberately drives chains in the `σ < 0` direction. So:

- **Fixture obligation.** Construct a self-intersecting *growing* chain with a known answer (a "U"
  or dumbbell outset until the arms cross; a concave CW island grown until it closes). Confirm the
  extreme that recovers the true envelope, and confirm the normalisation anchor flips with `σ`
  rather than being hard-anchored to the maximum. This is a *measured* result, written as a
  known-answer test — never derived on paper and trusted.

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
- The extreme (min vs max) and the normalisation anchor are chosen by the **sense** `σ =
  sign(amount · h)`, not by `sign(amount)` and not by handedness alone.
- The inversion guard fires only when `σ > 0` (shrinking); it discards handedness-flipped children.
- `σ > 0` reproduces today's behaviour exactly. `σ < 0` with a self-intersecting field is
  **unproven** in the current code and is a known-answer fixture, not a derivation.
- Preserve the numerical scar tissue verbatim: crossings computed in **double precision**;
  tangent-seam **fusion**; **small-feature pruning**; interiority measured in **absolute distance**
  (never a parameter tolerance, which scales with edge length); and **transversal-kiss rejection**
  (a tangential touch is not a crossing — its `±1` is pure float noise and must not enter the walk).
