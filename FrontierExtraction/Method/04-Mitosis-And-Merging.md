# Mitosis and Merging — the frontier across many chains

[03-The-Frontier.md](03-The-Frontier.md) extracted the true frontier of **one** chain's offset.
But a system is a *set* of chains — an outer wall, its islands, several disjoint features — and
their offsets interact. Two things happen when they do, and they are **different operations** that
the original code, having only needed one combined behaviour, fused into a single function. The
reconstruction separates them.

- **Mitosis** — one region *splits* into several. As an offset advances, a single ring can pinch
  in two where the region between two obstacles narrows past the point one loop can bridge. The
  name is literal: one becomes many.
- **Merging** — several chains, each already resolved against *itself*, are combined to find the
  system's current set of frontiers. This is the cross-chain operation: chains that were separate
  now meet, and we must decide what the shared material region is and where its boundary runs.

Both are the *same frontier idea* one level up: instead of a crossing-number field over one
chain's self-intersections, we work with a **winding-number field** over the whole arrangement of
loops, and the frontier is again the boundary of a selected region of that field.

## The winding field and the material region

Every closed loop contributes to a **winding number** `w(p)` at each point `p`: `+1` for each CCW
loop that encloses `p`, `−1` for each CW loop. `w` is an exact geometric invariant (computed by
accumulated subtended angle, not ray-casting, so it is stable — it can only wobble when `p` sits
*on* a curve, never because distant geometry shifted).

The **material region** is a level-set of `w`. In the proven, inward-only code it is fixed as
`{ w · s ≥ 1 }`, where `s` is a system sign the code *auto-detects* from the outermost loop's
handedness. That auto-detection is exactly the frozen policy this reconstruction unfreezes: it is
the merge's analogue of the single-chain "always take the minimum," and it fails the moment the
caller wants a region the heuristic does not expect (an enclosing keep-out perimeter drives the
detector to the wrong sign, and the material region inverts). See the difference-merge below.

The generalised merge takes the region **as an argument**, not a guess:

- an explicit **region predicate** — which sign of `w`, and at which threshold, counts as material
  (`w ≥ 1` for a union or a difference; `w ≥ 2` for an intersection);
- a **role per input loop** — *minuend* (a positive contributor, the region we are keeping) versus
  *subtrahend* (a negative contributor, geometry we are removing). A cluster that contains **no
  minuend** contributes no material and is dropped, so disjoint subtrahend geometry cannot survive
  as a spurious frontier.

With those two exposed, the same machine computes union, difference, and intersection, and the
inward-only material region is simply the default predicate. The *frontier extraction* underneath
is unchanged: cut every loop where another crosses it, and **keep an edge piece iff the material
region lies on exactly one of its two sides** — material on both sides is interior (discarded),
void on both sides is exterior (discarded), one of each is a true frontier. The survivors stitch
back into loops at the crossing points; because they meet there exactly, this is *reassembly, not
inference*.

## Clustering: only interacting loops judge each other

Loops that neither cross nor contain one another are **independent worlds** and must never judge
each other's handedness. So the merge first clusters loops by interaction (they cross, or one
contains the other), and resolves each cluster alone. This is not an optimisation; it is
correctness — a disjoint feature elsewhere in the section has no standing to flip the sign of the
region we care about here. (It is also the mechanism by which the *role* rule above drops
subtrahend-only clusters: a cluster with no minuend simply yields nothing.)

## The handedness laws: why two insides cannot make an outside, but two outsides can make an inside

Here is the deep asymmetry the whole method obeys, stated first as the rule the single-chain guard
already encodes and then as the topological fact underneath it.

Call a **CCW loop an "inside"** (it bounds material) and a **CW loop an "outside"** (it bounds a
void / keep-island). Offsetting changes topology, and the law is:

> **An outside may gain a new inside, but an inside may never gain a new outside.**

Concretely, on a single source chain: a CW source (an outside) may legitimately give birth to a
CCW child (a new inside) — an *outset pinching off a new cavity* is a real event. But a CCW source
(an inside) producing a CW child (a new outside) is **inversion** — the ring collapsed through
itself — and that child is discarded ([03-The-Frontier.md](03-The-Frontier.md)).

Why the asymmetry is real, and not a convention:

- **Two insides cannot create a new outside.** Insetting *removes* material. When two material
  regions (two insides) advance toward each other, they can only *merge* — the void between them is
  consumed, not manufactured. Nothing about shrinking material can wrap a *new* void out of a place
  that was solid on both sides. The only way an "inside" appears to sprout an "outside" is the
  collapse-through-zero fold, which is an artifact, not a boundary — hence the discard.
- **Two outsides can create a new inside.** Growing regions *enclose*. When two voids (two
  outsides) expand toward each other, the strip of material caught between them can be pinched off
  and survive as a **new isolated inside** — a genuine island born of the growth. This is a real
  topology change, and it is *permitted*: a growing configuration (`σ < 0`, which never collapses)
  is exactly the regime that is allowed to birth new pieces.

So the law is the topological shadow of the sense `σ`. Shrinking (`σ > 0`) can only *destroy*
voids and *split* material (mitosis), and any apparent new void is an inversion to reject. Growing
(`σ < 0`) can *enclose* voids and *pinch off* new material, and those new insides are legitimate.
The single-chain inversion guard and this merge-level law are one principle seen at two scales:
**handedness may advance toward "more inside" freely, but a flip toward "more outside" is only ever
the fold.**

## How it all ties together

The three levers of the method are not independent knobs; they are locked together by the sense:

- **Offset** (`amount`) and **handedness** (`h`) combine into the **sense** `σ = sign(amount · h)`.
- The **sense** chooses the frontier extreme in single-chain extraction (min for `σ > 0`, the outer
  envelope for `σ < 0`) **and** gates the inversion guard (fires only for `σ > 0`).
- The **crossing number** (one chain) and the **winding number** (many chains) are the same
  integer-field idea at two scales, and in both the frontier is the boundary of a region *selected*
  from that field — a selection that, generalised, is an argument rather than a constant.

The old code answered one question — the internal frontier of an inward offset, its material region
auto-detected — and answered it well. The reconstruction keeps that answer as the default and makes
the *question* sayable, so mitosis, union, difference, and intersection are one machine reading
one honest set of parameters.

## Carried across to the reconstruction

- **Mitosis** (single region splitting) and **merging** (combining resolved chains) are distinct
  operations, not one fused function.
- The merge is a **winding-rule** region extractor: expose the **region predicate** (sign +
  threshold) and a **minuend/subtrahend role** per loop; the inward-only `{ w · s ≥ 1 }` with
  auto-detected `s` is the *default*, not the mechanism.
- **Cluster by interaction** before judging; disjoint loops never cross-judge; a minuend-free
  cluster yields nothing (this, not a distance filter, is what makes an enclosing keep-out
  perimeter drop out cleanly).
- The **handedness law** — outside-may-gain-an-inside, inside-may-never-gain-an-outside — is the
  topological form of the sense, and must hold identically in the single-chain guard and the
  cross-chain merge.
- Preserve the merge's scar tissue: **winding by subtended angle**, not ray-casting; interaction
  and containment probed at each loop's **extreme point** (on the loop, far from the others, so the
  answer changes only at a real geometric event); and **oriented reassembly only** — a stitched
  frontier that will not close in travel order must fail visibly, never be silently reversed into a
  wrong-handed loop.
