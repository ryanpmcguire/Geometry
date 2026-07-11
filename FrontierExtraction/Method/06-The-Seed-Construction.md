# The Seed Construction — cut against avoid

This is the synthesis. Everything in documents 01–05 exists to make *this* construction exact:
given two systems of randomly-oriented input chains, produce the **seed** — the true maximum
extent of the tool *centre* that accomplishes as much of the cut as possible while never letting
the tool enter the avoid region. The seed is the boundary the clearing recursion then insets from;
getting it right is getting the whole toolpath's extent right.

## The inputs and the target

Two systems arrive, each a set of closed chains, each **randomly oriented** (the caller has not
promised any handedness):

- **Cut** — chains enclosing area where cutting is *encouraged*: the tool must maximally touch it.
- **Avoid** — chains enclosing area where cutting is *forbidden*: **no part of the circular tool
  of radius `R` may enter it.** (Safety assumption: the avoid system must enclose *all* keep-out.
  Anything enclosed by neither system is treated as free air — cuttable but empty.)

Reason about the tool *centre*. The tool is a disk of radius `R`.

- To keep the tool out of avoid, the centre must stay at least `R` away from it: the centre is
  excluded from **avoid ⊕ R** (the avoid region grown by `R`).
- To touch the cut at all, the centre must be within `R` of some cut point: the useful centres are
  exactly **cut ⊕ R** (the cut region grown by `R`). Its boundary is where the tool is *tangent* to
  the cut edge — the furthest the centre can sit and still reach the cut.

The maximum safe, useful extent of the tool centre is therefore the exact region

> **seed = (cut ⊕ R) ∖ (avoid ⊕ R).**

Both operands are the input **grown by one tool radius**; the seed is their difference. Everything
below is how to compute that difference with nothing but the blind offset and the consistent-
handedness merge — no Minkowski sum, no special-cased boolean.

## The construction

**Step 1 — orient by nesting parity, with opposite bases.** The inputs are unoriented, so assert
orientation from containment (exact signed area, never trusted from the caller):

- The **cut** system is a region: outer loops **CCW**, their holes **CW**, holes-in-holes CCW —
  handedness alternating with nesting depth.
- The **avoid** system is oriented with the **opposite base**: outer loops **CW**, their holes
  **CCW**, and so on by parity.

The opposite base is not a second rule — it is the *same* rule ([05](05-Opposite-Handed-Islands.md))
seen one level up. From the perspective of the clearable material, the whole avoid region is an
**island**: a keep-out sitting inside the space the tool works. "Flip the handedness of islands"
therefore flips the avoid's entire base, and then its own holes flip again by parity. The cut and
avoid are not two conventions meeting; they are **one nesting tree**, cut at the root, avoid a
child, handedness following parity all the way down.

**Step 2 — grow both.** Dilate each system by `R`, using the blind offset. Because the two systems
carry opposite base handedness, growing them takes **opposite `offsetLeft` signs** — this is
[02](02-Offsetting-Semantics.md)'s "CCW and CW respond oppositely to the same value," read
forward:

| system | base handedness | to **grow** (⊕ R) | sense `σ` |
|---|---|---|---|
| cut   | CCW | `offsetLeft(-R)` | `−1` (growing) |
| avoid | CW  | `offsetLeft(+R)` | `−1` (growing) |

Within *each* system the single sign is applied across its own alternating-handed members (a CCW
outer and its CW holes both grow under `offsetLeft(-R)`), exactly as [05](05-Opposite-Handed-Islands.md)
requires. Between the two systems the signs are opposite, because here we **grow both regions** —
and growing opposite-handed regions is opposite work. This does not contradict "same sign to
differently-handed chains": that rule advances *one* region's frontier; the seed instead *expands*
two complementary regions, which is a different combination of the same primitive.

**Step 3 — merge under the material-on-left invariant.** Feed every grown chain — cut and avoid
together — into one merge ([04](04-Mitosis-And-Merging.md)). Because the orientation was chosen so
that clearable material is consistently on the left of every edge, the merge needs no auto-detected
system sign: the material region is simply `{ winding ≥ 1 }`, and that region *is* the difference.
The surviving frontiers, spliced at the crossings of grown-cut against grown-avoid, are the seed.

## Why the merge is the difference (and why handedness replaces the boolean)

With cut CCW (`+1` inside) and avoid CW (`−1` inside the avoid island), the winding at a point is
`+1` for being inside the grown cut and `−1` for being inside the grown avoid:

- inside cut, outside avoid → `+1` → **material** (the reachable, safe seed).
- inside cut, inside avoid → `+1 − 1 = 0` → not material (the `R` standoff collar where avoid
  encroaches — correctly excluded).
- outside cut → `0` or `−1` → not material.

So `{ winding ≥ 1 }` is exactly `(cut ⊕ R) ∖ (avoid ⊕ R)`. The difference was never computed as a
boolean; it *fell out* of choosing the avoid's handedness to be the island of the cut. **Handedness
is the region selector** — which is why this construction succeeds where an auto-detected merge
fails on an enclosing keep-out: an enclosing avoid perimeter (CW, grown) has non-material on both
its sides, so the frontier test discards it with no spurious loop and no sign to guess.

## The crucial scenario: the hole

Sketching cases, the construction holds across a wide variety — until a **hole**: a void inside
material (the everyday internal pocket). Here the avoid system is the kept material, which has the
pocket as a **hole in it**. The temptation is to orient the whole avoid uniformly CW. That is the
one case where it breaks, and the winding shows exactly why. Evaluate deep inside the pocket
(cut = pocket, avoid = perimeter + pocket-hole):

- pocket-hole oriented **CW** (uniform, wrong): `+1(cut) − 1(perimeter) − 1(hole) = −1` → **not
  material.** The pocket — the very thing we mean to cut — is wrongly excluded.
- pocket-hole **flipped to CCW** (the island rule): `+1(cut) − 1(perimeter) + 1(hole) = +1` →
  **material.** Correct.

The resolution is topological, not a patch: **a hole within a material is itself an island**, so by
"flip the handedness of islands" it must be CCW — the *opposite* of the avoid's own CW base. The
rule that flips a region's holes is the *same* rule that flipped the avoid relative to the cut, and
it keeps flipping by nesting depth. This is the important discovery the construction turns on:

> **"Flip the handedness of islands" is a single universal rule of nesting parity, and it holds not
> only within one region but across *complementary* systems — the avoid is the cut's island, and
> the avoid's hole is the avoid's island. Handedness is determined by depth in one shared nesting
> tree, and complementary facts (cut vs. avoid) are just opposite parities in that one tree.**

Once handedness is assigned by global nesting parity, the seed construction is correct in every
scenario tested on paper, holes included, with no case analysis — one dilation sign per system, one
merge.

## Worked examples: the bore and the collar

**A plain bore.** Drilling a round hole: the cut is a single CCW disk — the bore material to remove,
no islands, "a counter-clockwise inner region only." The avoid is the kept part: a CW outer
perimeter with the bore as a **CCW island** (its CW base flipped once by the island rule). Notice
the bore wall circle appears in *both* systems — as the cut's outer boundary and as the avoid's
inner island — and it carries the **same handedness (CCW) in both**. That is not a coincidence, it
is a correctness invariant:

> **Complementary systems must agree on the handedness of a wall they share.** Parity is assigned by
> depth in one nesting tree, so the shared wall gets one answer no matter which system you reach it
> through. Two systems that disagreed on a shared wall could not merge — so an *observed* agreement
> is a cheap sanity check on the orientation pass.

Growth pulls the one circle into two. The cut grows it outward to `bore + R`; the avoid grows its
island inward (kept material encroaching) to `bore − R`. Winding deep in the bore:
`+1(cut) − 1(perimeter) + 1(bore-island) = +1` → material. In the standoff collar,
`bore − R < r < bore + R`: `+1 − 1 + 0 = 0` → not material. The seed is the disk of radius
`bore − R` — the tool centre standing `R` off the wall. The cut's `bore + R` copy is scaffolding,
discarded by the frontier test; the surviving frontier is the avoid's `bore − R` island.

**Why the duplication is the point.** That a shared wall appears twice and grows to two copies is
not redundancy — it is the mechanism that resolves wall-sides and air-sides in one pass, with no
per-edge tagging. Where a cut edge is backed by kept material, the avoid supplies its `−R` copy and
*that* survives: an `R` standoff. Where a cut edge faces open air, there is no avoid there, so only
the cut's `+R` copy exists and it survives: the overrun. The open-air push faked this distinction
with a boolean per edge; here it falls out of two systems drawn from the same walls.

**A collar inside a bore.** Bore out a ring but leave a central pillar — a cylindrical collar of
kept material: the cut is now an annulus, a CCW outer at the bore wall and a **CW inner island**
around the collar ("it would itself also have a clockwise inner island"). The collar is kept, so it
also appears in the avoid as a small region of its own — CW-flavoured, its outer boundary CW. So the
collar wall is **CW whether reached as the cut's island or the avoid's outer**: agreement again, one
level deeper. No special case is needed — the avoid's collar copy grows outward to `collar + R` and
survives as the standoff around the pillar, exactly as the perimeter did around the bore, while the
cut's island copy is scaffolding. Nesting depth alone decided every handedness, and the seed stands
`R` off *both* walls.

## The load-bearing dependency (read before trusting this)

There is one honest catch, and it ties this document back to [03](03-The-Frontier.md). **Both
dilations are growing offsets — `σ = −1` for cut and for avoid alike.** The entire seed construction
lives in the `σ < 0` regime. And a growing chain that self-intersects (dilate any non-convex cut or
avoid loop until a concavity closes) is precisely the single-chain case that document 03 flags as
**unproven**: the current `minLevel` extraction was only ever validated for `σ > 0` (shrinking), and
almost certainly extracts the wrong extreme for `σ < 0`.

So the seed is only as sound as the `σ < 0` self-mitosis it rests on. This construction **cannot be
wired onto the existing extraction as-is** — it would exercise, on real parts, exactly the growing-
self-intersection case the old code never answered. The `σ < 0` extreme (and its normalisation
anchor) must be pinned by known-answer fixture *first*; the seed is the caller that makes that
fixture non-optional.

## Carried across to the reconstruction

- Orient both systems by **nesting parity from exact signed area** — cut base CCW, avoid base CW —
  never trusting input handedness. The opposite base is the island rule applied to the whole avoid.
- The seed is `(cut ⊕ R) ∖ (avoid ⊕ R)`: **grow both by the tool radius** (cut `offsetLeft(-R)`,
  avoid `offsetLeft(+R)` — opposite signs, both growing), then one merge at `{ winding ≥ 1 }`.
- Handedness is the region selector; the difference is not a boolean but a consequence of the avoid
  being the cut's island. An enclosing avoid perimeter drops out by the frontier test, with no
  auto-detected sign and no locality filter.
- **A hole within material is an island → its handedness flips.** Assign every loop's handedness by
  its depth in the single cut/avoid nesting tree. This is the one rule that makes the hole case, and
  therefore ordinary pockets and bores, come out right.
- **Complementary systems agree on shared walls.** A wall the cut and avoid share must receive the
  same handedness from both (parity is one tree). A shared wall the two systems orient *differently*
  is a bug in the orientation pass, not geometry to reconcile — the reconstruction can assert this as
  a cheap invariant check.
- **A shared wall grows to two copies, and that is the mechanism, not waste.** The avoid's copy wins
  at kept walls (standoff `−R`); the cut's copy survives where no avoid backs it (air overrun `+R`).
  Wall-standoff and air-overrun come from one construction with no per-edge boolean.
- **This whole construction is `σ < 0`.** It depends entirely on the growing-direction self-mitosis
  of [03](03-The-Frontier.md), which is the reconstruction's first and highest-priority fixture. Do
  not let a real part ride on it until that is measured, not assumed.
- Degeneracies to fixture: cut and avoid sharing a wall (their dilations overlap in a `2R` band —
  lean on transversal-kiss rejection so the coincident boundary does not manufacture false
  crossings); deep nesting parity; and an avoid that does not reach this slice at all (no avoid
  loops → the seed is the lone dilated cut, the fully-open extent).
