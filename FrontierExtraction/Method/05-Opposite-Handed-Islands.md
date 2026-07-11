# Opposite-Handed Islands — why the convention is forced

The doctrine says an outer profile is counter-clockwise and an internal island is clockwise
([02-Offsetting-Semantics.md](02-Offsetting-Semantics.md)). This is not a labelling convenience —
it is **forced**, and getting it wrong does not merely mislabel the geometry, it destroys the
merge. This document is the argument for *why* islands must carry the opposite handedness to the
profile that contains them, and why that single choice is what makes chain identity an *algebra*
rather than a pile of special cases.

## The one invariant: material on the left, everywhere

The whole method rests on a single global invariant: **material is on the left of the travel of
every boundary edge.** Extraction, offsetting, and merging are all just consequences of it. So the
representational question is entirely: *what handedness must each loop carry for that invariant to
hold?*

- An **outer profile** bounds material on its *inside*. Walk it counter-clockwise and your left
  hand points inward — toward the material. **Material on the left ⟹ outer is CCW.**
- An **internal island** is keep-out: it bounds material on its *outside* (the island itself is a
  hole in the clearing). Walk it clockwise and your left hand points *outward* — toward the
  material. **Material on the left ⟹ island is CW.**

So opposite handedness is not chosen; it is *dictated* by demanding one consistent side for
material. The outer and its islands wind oppositely precisely so that a single sentence —
"material on the left" — is true of every edge in the system at once. In the language of planar
topology this is exactly `∂` of the material region: the boundary of a 2-chain with holes is its
outer loop plus its holes with *reversed* orientation. The convention is a boundary operator, not
a preference.

## Why this makes merging produce consistent chains

Now offset the system by one blind value ([02](02-Offsetting-Semantics.md)). The CCW outer insets;
the CW island grows; both advance the internal frontier; and — the point — **both still keep
material on their left**, because the blind offset moves each edge sideways without touching its
handedness. The invariant survives the offset.

When the shrinking outer and the growing island finally meet and cross — the clearing pinching
closed around the island — look at the crossing. The outer arrives travelling one way with
material on its left; the island departs travelling the *opposite* rotational sense, also with
material on its left. Their tangents are in **opposition**. That opposition is exactly what a
frontier *is*: a place where two strands run antiparallel because material is on one side and void
on the other. The merge can therefore **splice** them — enter along the outer, hand off to the
island at the crossing, continue — and because material stays on the left across the hand-off, the
spliced loop is **consistently handed end to end.** One new well-formed chain comes out, and its
handedness is unambiguous because every stoicheion in it agrees on which side the material is.

This is the merge working: oppositely-handed inputs cross in opposition, opposition marks a real
frontier, and the frontier reassembles into consistently-handed offspring.

## The counterfactual: same-handed chains, differently-signed offsets

Now suppose we refused the convention and drew *islands counter-clockwise too* — same handedness
as the outer. The geometry can still be made to move correctly: inset the outer with a positive
offset, and *outset* the island with a negative one. Same handedness everywhere; different offset
values per chain. Superficially it clears the same material.

It breaks, and it breaks at exactly the invariant. A CCW island cannot have material on its left —
walk it counter-clockwise and the left hand points *inward*, into the island, which is keep-out,
not material. So a CCW island necessarily carries material on its **right**. The system no longer
has one consistent side for material: the outer's edges have it on the left, the island's on the
right.

Trace the consequence to the crossing. Two same-handed (both CCW) chains, traced each along its own
handedness, meet at their intersection in **agreement** — both circulating the same rotational
sense — not in opposition. There is no antiparallel hand-off, because there is no consistent
frontier to hand off *along*: the two strands do not bound a common region on a common side. Splice
them anyway and you produce a candidate chain whose stoicheia are of **inconsistent handedness** —
some edges carrying material on the left, some on the right — which is not a chain in any usable
sense. The merge has nothing coherent to extract. This is the precise failure the convention exists
to prevent: *there is no longer any consistently-handed frontier.*

The deeper reason the same-handed scheme also fails in spirit: it replaces one blind, global
operation with **per-chain decisions** — this chain gets `+amount`, that one gets `−amount`. That
is exactly the external judgement the blind method was built to abolish
([01-Stoicheia-And-Handedness.md](01-Stoicheia-And-Handedness.md)). The moment the offset sign
depends on which chain you are looking at, correctness depends on classifying chains correctly
*before* offsetting — the very thing we wanted the geometry to decide for itself.

## The principle: same offset to differently-handed chains

Here is the choice, stated as sharply as it deserves:

> **Apply the *same* offset value to *differently-handed* chains — never *different* offset values
> to *same-handed* chains.**

The two are geometrically capable of moving the same boundaries. They are *not* equivalent, and the
difference is total:

- **Same value, differently handed** (the convention). One blind operation, no per-chain decision,
  the material-on-left invariant preserved, crossings in opposition, frontiers consistent, merge
  products consistently handed. A **path-theoretic algebra of chain identity** emerges: chains
  compose and split under one boundary operator, handedness behaves like a conserved sign, and the
  identity of an offspring chain is determined — not chosen — by the handedness algebra of its
  parents.
- **Different values, same handed** (the counterfactual). Per-chain classification restored, the
  invariant broken, crossings in agreement, no consistent frontier, offspring of mixed handedness.
  No algebra — just cases.

The convention is what lets the method be *blind*, and blindness is what lets it be *correct*: one
sign, applied to a system whose handednesses were chosen so that the single blind offset means the
right thing to every chain, and whose crossings therefore always resolve into consistently-handed
descendants. Opposite-handed islands are the price of admission to that algebra, and the algebra is
the whole method.

## The mirror: neither handedness is privileged

Nothing above depended on *which* handedness we called "outer." It depended only on the outer and
the island being **opposite**, and on the offset sense being chosen to advance the internal
frontier. The absolute assignment is a free global choice, and the proof is a clean symmetry.

Reversing a chain's travel does two things at once: it flips the chain's handedness (CCW ↔ CW) and
it swaps its left and right. So for any edge, `offsetLeft(reverse(e), amount)` moves toward what was
`e`'s *right* — which is identical geometry to `offsetLeft(e, −amount)`. Therefore:

> **Reversing every chain's travel and negating the offset value is a symmetry of the whole
> operation — it produces the identical geometry, edge for edge.**

Apply that symmetry to the standard convention and read off its mirror:

- **Standard.** Outer CCW, island CW, offset by **+amount**, material on the **left** of every edge.
- **Mirror.** Outer CW, island CCW, offset by **−amount**, material on the **right** of every edge.

Everything holds in the mirror, term for term. A CW outer is *inset* by the negative offset — its
interior is on its right, so moving right shrinks it — and a CCW island is *grown* by that same
negative offset; both advance the internal frontier, exactly as before. Material sits consistently
on the right of every edge, so the invariant is intact, merely reflected. Oppositely-handed chains
still cross in opposition, so the merge still splices consistently-handed offspring. The clearing
that results is identical geometry; only the names of the orientations changed.

The reason it works so cleanly is that the quantity the method actually reads is invariant under the
mirror. The **sense** `σ = sign(amount · h)` — which chooses the frontier extreme and gates the
inversion law ([03-The-Frontier.md](03-The-Frontier.md)) — is unchanged when both `amount` and `h`
flip sign, since `sign((−amount) · (−h)) = sign(amount · h)`. Min-versus-max extraction, shrink-
versus-grow behaviour, the inversion guard: all identical in both conventions, because none of them
was ever a function of absolute handedness or absolute offset sign — only of `σ`.

So the *forced* part of the convention is **opposition** — a region and its holes must disagree —
and the *free* part is **which one you call CCW**. A system picks one and holds it (orienting by
exact signed area is how it asserts the choice); the two conventions are then a global reflection
apart, under which the method is exactly invariant. There is no privileged handedness, only a
required disagreement between a region and its holes.

## Carried across to the reconstruction

- The outer-CCW / island-CW convention is **forced by "material on the left," not chosen** — it is
  `∂` of the material region. The reconstruction must treat it as an invariant, and it may *assert*
  it (orient by exact signed area) but must never *depend on the caller* having a preferred style.
- Only the **opposition** of region and hole is forced; the **absolute** handedness is free. The
  standard (outer CCW, `+` offset, material-left) and its mirror (outer CW, `−` offset,
  material-right) are a global reflection apart, and `σ = sign(amount · h)` is invariant across it —
  so the reconstruction fixes one convention for definiteness, never because the other is wrong.
- Advance the frontier with **one blind offset value across the whole system**; never branch the
  offset sign per chain. If the implementation finds itself choosing a sign based on a chain's role,
  it has left the algebra and the merge will not close.
- Consistent handedness of merge offspring is the **correctness signal**: a splice that yields a
  chain of mixed handedness is not a warning to smooth over — it means the inputs violated the
  invariant, and the reconstruction should fail visibly there rather than emit a malformed chain.
