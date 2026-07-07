# Geo — Chains, Profiles & the Slice Strategy

This is the conceptual guide to the Geometry library's toolpathing stack: what a
**stoicheion** is, what a **chain** is, what a **profile** is, and how the **slice
strategy** turns a prepared boundary into actual tool motion. The source files
(`Geometry.ixx`, `Chain.ixx`, `Strategy.ixx`) are the reference; this file is the map —
and, more importantly, the *why*. Each layer has a small number of load-bearing ideas
that look fiddly until you see the problem they solve. We'll do the problems.

The whole stack in one picture, bottom up:

```
Stoicheion   (segment / arc / circle / ellipse / spline ...)
  "one exact curve, carrying its own orientation"
    └─ Chain        an ordered, oriented run of stoicheia — travel IS meaning
         └─ Profile     one GENERATION: the chains bounding material at one offset depth
              └─ SliceStrategy   profiles in series → ordered, linked, led TOOLPATH
```

Each layer only knows the one below it. A stoicheion knows nothing of chains; a chain
knows nothing of generations; a profile knows nothing of tools. The strategy at the top
is the only thing that knows what a *tool radius* is.

---

## Stoicheion — the atom

Euclid's *Elements* is the Στοιχεῖα (Stoicheia); a single primitive is a **stoicheion**.
It is the fundamental geometric element: one exact curve that knows how to serialise
itself, tessellate itself into a polyline, measure distance to a probe point, intersect
with another stoicheion, and — the operation everything above is built on — **offset**
itself sideways.

The concrete roster (`Geometry.ixx`):

- **`Point2`** — a position. Encloses nothing, offsets nowhere; exists for sketching.
- **`Segment2`** — a straight run `a → b`.
- **`Circle2`** — centre `c`, radius point `a`, and a *chirality marker* `b` (see below).
- **`Arc2`** — centre `c`, endpoints `a`, `b`, and a *through-point* `d`: the arc runs
  "from A to B **through D**", which fixes which of the two possible arcs it is *and*
  which way it travels.
- **`Ellipse2` / `EllipseArc2`** — centre + conjugate semi-axes `u`, `v` (+ a parameter
  span for the arc).
- **`OffsetEllipse2`** — the *exact parallel curve* of an ellipse: base `(c, u, v)` plus
  a signed offset distance `d` and a parameter span.
- **`Spline2`** — a NURBS / B-spline (the catch-all for freeform edges OpenCASCADE hands
  back from a BRep section), likewise carried as base curve + signed offset `d`.

Three ideas make this roster a *system* rather than a grab-bag:

### 1. Orientation lives in points, not flags

A circle carries no signed radius and no `bool ccw`. It carries a marker point `b`
sitting a quarter-turn from its radius point `a` — and **which side `b` sits on *is* the
orientation**. An arc's through-point `d` does the same job. This sounds like a quirk
until you offset *through zero*: drive an inset circle's radius through its centre and it
re-emerges on the far side, inside-out, with the marker naturally landed on the other
side — the chirality flips *by geometry*, not by special-cased sign juggling. Every
downstream rule (which way an offset goes, what a loop's handedness is) reads orientation
off the points and therefore survives every degenerate case the same way the geometry
does.

### 2. Exact where a closed form exists; certified numerics where it provably doesn't

Segments, arcs, and circles intersect analytically — nothing is ever quantised. Signed
area contributions (`signedAreaTerm`) are closed-form Green's-theorem integrals per edge.
But the true offset of an ellipse is *not an ellipse* (nor any rational curve), and a
spline's offset is not a spline — so those cases are handled **numerically on the exact
curve**: `OffsetEllipse2` and `Spline2` evaluate `base(t) + d · normal(t)` directly, and
root-finding runs a certified scan-and-bisect (`scanRoots`) that finds every transversal
crossing to machine precision. The doctrine: *approximate the arithmetic, never the
model*. There is no tolerance at which the geometry itself was polygonised.

### 3. The offset families are closed

Offsetting a line gives a line; a circle, a concentric circle. For the curves whose
offsets leave their own kind, the representation absorbs the problem: an
`OffsetEllipse2` offset by `δ` is the same object with `d + δ`; same for `Spline2`.
**Offsetting never degrades the representation** — a circle stays a circle through fifty
generations, a spline stays one spline. This is why the strategy can recurse offsets to
extinction without accumulating tessellation error.

One more flag worth knowing: **`openAir`**. An edge marked open-air says "beyond me lies
free space, not part material — the tool may run off this edge." It is set by the
*caller* (it is a statement about the world, not the geometry) and consumed by the
strategy's preparation step (below).

---

## Chain — travel is meaning

A **`Chain`** is an ordered run of stoicheia connected end to start. It is the unit of
offsetting and, ultimately, the unit of tool motion. The one mental shift to make:

> A chain is not a *shape*. It is a *journey*. Its travel direction and handedness are
> semantic, and everything downstream reads meaning out of them.

The conventions that carry that meaning:

- **Material lies on the LEFT of travel.** A counter-clockwise (CCW) loop therefore
  bounds a region to be cleared (its inside is material-to-remove... its left is its
  interior); a clockwise (CW) loop bounds a *keep-island* growing in the middle of a
  clearing. This CCW-clears / CW-keeps rule is the **doctrine orientation**, and it is
  the *caller's* job to establish it before anything here runs (in the Sketch app the
  user's drawn directions are that preparation; in the CAM app, nesting-parity
  orientation + open-air marking condition the mesh section).
- **`closed`** — only closed chains bound regions, so only closed chains may be offset
  or asked region questions (winding, area). Open chains exist legitimately only as
  *results*: the link and lead segments a strategy weaves.
- **Handedness is measured two ways.** `signedArea()` is the exact Green's-theorem area
  (sign = winding). `totalTurning()` is the accumulated tangent rotation — robust where
  area is not: a burr lobe can outweigh a thin valid loop's *area* and flip its sign, but
  it can never flip the total *turning*, because offsetting preserves the tangent
  structure of its source. Where the two disagree, turning is the truth.
- **Start-point freedom.** A closed loop is the same cycle wherever you enter it, so
  `startedAt(k, p)` re-seats where the loop *begins* without touching travel or
  geometry. Strategies use this constantly (to enter a ring at the point nearest the
  tool). By contrast `reversed()` genuinely changes meaning — it exists for
  *conditioning inputs* whose winding is arbitrary, never for reorienting the method's
  own outputs.

### Identity and ancestry

Every chain carries a random 64-bit `id` and a `parent` id — the chain it descended from
under offsetting (0 = a root). Children reference parents, never the other way, so the
**ancestral forest is implicit in the chains themselves**: profile 1's chains are roots'
children, profile 2's chains point back into profile 1, and a later pass can walk the
whole family tree with nothing but the ids. When an offset ring splits into several
(narrow region — see mitosis), each child gets a fresh id and the same parent: siblings.

### LinkKind — what a chain *is* to the tool

Beyond its geometry, every chain wears a semantic tag:

| `LinkKind` | Meaning |
|---|---|
| `None`    | an ordinary cutting chain |
| `Cut`     | a stay-down hop to the next ring — the tool keeps cutting |
| `Retract` | the tool must fully retract before this hop |
| `LeadIn`  | a ramp easing *onto* a cut after a retract (the consumer descends it) |
| `LeadOut` | the mirror: easing *off* a cut before the retract lifts |
| `Finish`  | a cutting chain flagged as the thin finishing pass (own feed/speed) |
| `Entrance`| a cutting chain making the *first engagement* of its region — a suspected full slot, to be taken slowly |

The tags are semantic, not cosmetic: the CAM app decides retract heights and ramp
descents from them; the Sketch app merely picks colours.

---

## Offsetting — the crossing method

How does one chain become its offset? Not by the classic mitre/clip corner casework.
The method is deliberately **blind**, in two acts:

**Act one — the raw offset (`offsetRaw`).** Every edge slides to the *left of its own
travel* by the amount. No winding test, no preferred side: since material is on the
left, the chain's own chirality alone decides what "left" means — a CCW loop insets, a
CW loop outsets, *emergently*. Every corner is then bridged by the join arc known **by
construction**: centred on the original corner, radius = the offset amount, continuing
the first edge's travel. Converging corners are forced to take the long way round (a
little loop that exits the shape and returns); diverging corners take the short way. No
decision was made anywhere — and the price is that the raw result **self-intersects**.

**Act two — extraction by crossing number.** Those self-crossings are not a bug; they
are the raw material. The raw chain is split at every self-crossing and walked in order,
accumulating a count: each time the walk passes a crossing, the *other* strand's travel
direction there bumps the count up or down. Fragments sharing the minimum accumulated
number are the genuine offset boundary; everything above it is the excursion loops the
old mitre/clip code would have hand-trimmed. The minimum-level fragments meet
end-to-start at the crossing vertices exactly, so stitching them back into loops is
*reassembly, not inference*.

**The handedness law (the trial).** One rule polices the output: *an outside may gain a
new inside, but an inside may never gain a new outside.* A CW source can legitimately
birth a CCW child (an outset pinching off a new cavity), but a CCW source producing a CW
child is **inversion** — the ring collapsed through itself — and that candidate is
discarded. This is how recursion terminates honestly: a ring too small to inset again
simply produces nothing. Extinction, not error.

---

## Mitosis — when loops meet

`offsetValid` runs each chain in its own universe. But a *generation* is a set of chains,
and their offsets can collide: an outer boundary shrinking inward meets an island growing
outward. `Chain::mitose` resolves the encounter, and its name is literal — this is where
one ring **splits into several** when the region between obstacles narrows.

The mechanics, briefly: loops are clustered by interaction (they cross, or one contains
the other — disjoint clusters are independent worlds and never judge each other). Within
a cluster, every loop is cut at its crossings with every other, and a piece survives iff
it is a true **frontier** of the material: total winding says "inside" on exactly one of
its two sides. Both sides material → interior of the union, discarded. Both sides void →
swallowed or escaped, discarded. The survivors stitch back into well-formed loops at the
crossing points. Nothing is assumed: the material region is *measured* by winding, with
the cluster's own orientation read off its loops.

Ancestry threads through: each stitched offspring takes a fresh id, and its `parent` is
the parent of its **dominant contributor** by length — so lineage always lands on the
previous generation, and the family forest stays sound even through splits and merges.

---

## Profile — one generation

A **`Profile`** is one generation of the system: *the set of chains that together bound
the material at one offset depth*. Profile 0 is the prepared input; profile n+1 is
profile n offset by one step. The entire two-act method hides behind one call:

```cpp
Profile next = current.offsetBy(amount);
//   stage A — each chain offset in its own universe (offsetValid:
//             blind offset, crossing walk, minimum level, handedness law)
//   stage B — the validated chains meet (mitose: cluster, frontier test, stitch)
```

So *recursive offsetting is just a loop*, and an exhausted profile (everything inverted
away) returns empty — which is precisely the termination condition the strategy relies
on. The profile series **is** the generational record: keep every profile and you hold
the full ancestral forest, generation by generation, parent ids pointing backward.

---

## SliceStrategy — from generations to motion

`Geo::SliceStrategy` (in `Strategy.ixx`) is the manager that turns a prepared profile
into a final toolpath, for **one planar slice**. Its contract is strict, and worth
internalising because it is what keeps this layer pure:

- **The caller prepares the input.** Doctrine orientation established, open-air edges
  marked — *before* the strategy sees anything. The strategy neither knows nor cares
  whether a human drew the directions or a mesh section was conditioned into them.
- **The strategy returns everything.** The generation series (`result.profiles`),
  per-step stats and verdicts (`result.steps`), and the final toolpath
  (`result.toolpath`) — whose chain order and travel directions **are** the tool's
  motion, with direction (order) and milling sense (handedness) already resolved.
- **The caller consumes the result blindly.** Any further processing is strictly
  presentational. If a consumer finds itself "fixing" a toolpath, the fix belongs here.

Usage is three moves (or one, via `runSliceStrategy`):

```cpp
Geo::SliceStrategy strategy({ .toolRadius = r, .stepover = 0.6f });
strategy.ingest(std::move(preparedSeed));   // the caller's half of the contract
Geo::SliceResult& out = strategy.run();     // prepare → generations → steps → toolpath
```

### The generation ladder

`run()` drops open chains at the door (offsetting has no concept of them), then applies
the strategy's own preparation: every open-air edge is **offset** out into the free
region — each alone, as the true parallel curve of itself (with material on the left of
travel, the air is the right, so this is exactly the *opposite* of the inset the
generations run). The gaps the push opens are rejoined with **arcs tangent to what they
join**. Between two *pushed* neighbours: the same arc-known-by-construction the
recursive solver builds at every offset corner — centred on the original corner, radius
= the push, tangent to both pushed edges (a corner point offsets to an arc; a chamfer
is not the parallel boundary). At a *wall/air transition*: the tangent **fillet** from
the perimeter's unmoved endpoint onto the pushed edge — at a tangency-reversal tip this
is the half-circle of radius push/2, whose generation-1 inset grazes the tip vertex at
exactly one tool radius. The fillet passes *through* the vertex leaving along the
wall's own tangent, staying inside the wedge the tagging certified as air (a cap
*centred* on the vertex would bulge across the wall — that construction is forbidden).
Straight bridges remain only where an arc would be wrong: converging corners (the
long-way-round loop self-crosses, and a seed must stay a simple loop) and lateral steps
at mid-edge wall/air splits (a tangent join there is an S-biarc; the fillet
construction self-reports the mismatch and stands down). Every join lies in vacated
ground, so it is marked open-air itself. The push amount must fund **both**
sides of generation 1's inset: a wall-backed sliver erodes by `2R` (R from the wall side
*and* R from the air side — the offset method is uniform and grants the air no
discount), so the push runs past `2R` with margin, not merely past the exact-tangency
`R`. (Deciding to do this at all is *toolpathy* logic, which is why the call lives in
the strategy, not the chain method.)

Then the ladder, for the Profile kind:

1. **Generation 1** offsets by **exactly one tool radius** — the boundary-clearance
   ring, whose tool grazes the part wall.
2. **Generation 2**, when `finishPass` is on, is a single **thin** inset
   (`finishWidth × R`), tagged `Finish` — so a consumer can run the wall's final skim at
   its own feed. A finishing inset too wide for a narrow region would *mitose* — split
   into more chains than its parent — destroying the single continuous skin; the
   strategy narrows it until it no longer splits, and records a warning.
3. **Every later generation** advances by `stepover × R`, recursing until extinction or
   the `maxGenerations` cap.

### Steps and the sanity rule

`buildSteps` measures each generation (chain count, signed area, cleared area, cutting
length) and applies the final verdict: a step whose **total signed area is negative**
claims an unbounded region — an outside grown past an inside without touching, or a lone
outer trace. It is recorded and measured, but `included = false`: the toolpath skips it.

### Ordering the rings

`flattenDepthFirst` linearises the ancestral forest **depth-first**: once a ring splits,
one child lineage is followed to extinction before backing up — just far enough to the
nearest unaccounted sibling. Skipped generations are traversed but not emitted, an
orphan sweep appends anything whose lineage broke, and the outer ring obeys `emitOuter`
(a roughing pass can defer the wall entirely to a later finishing run). Two
post-processing axes then apply, independent by design:

- **`reverse`** flips chain *order* only — outside-in becomes inside-out.
- **`climb`** flips chain *handedness* only — climb keeps the method's travel (material
  on the left); conventional reverses every chain.

This ordering stage is deliberately a seam: it decides *sequence*, and nothing else —
smarter orderings slot in here without touching linking or leads.

### Entrance tagging — the full-slot suspects

Once the final execution order is fixed, each cutting chain is judged for *engagement*:
a chain run while **neither its parent nor all of its children are already cut** has no
cleared swath beside it — the tool meets full-width material, a suspected full slot —
and is retagged `Entrance`. The rule is execution-order truth, not ancestry alone:
outside-in it catches each region's wall ring; inside-out (`reverse`) the innermost
leaf; a deferred outer (`emitOuter` off) leaves its child correctly suspect. Only
*emitted* chains count as relief — an unemitted child cleared nothing — and a leaf with
an uncut parent is always suspect (zero children never vacuously relieve). Degenerate
sweeps are `Entrance` by definition. A consumer gives entrance chains their own slow
feed, exactly as it does `Finish`.

### Linking — measured, never assumed

Between consecutive cutting chains the strategy weaves connector chains, and the central
discipline is that **whether the tool may stay down is measured, never inferred from
ancestry** (ancestry proves nothing: a parent's other children can lie anywhere). A
stay-down `Cut` link requires *both*:

- **Reach** — the hop to the next chain's nearest entry is at most **one tool radius**.
  A circle of one tool radius at the tool's position must touch the next chain; that is
  the whole budget.
- **Clearance** — the straight hop must not cross **any of the source geometry**.
  Crossing the seed means cutting through a wall; that hop must retract instead.

On a cut link, the next loop is re-seated (`startedAt`) to begin at its point nearest
the tool, keeping the flow of motion. On a retract link, with the start free, it is
re-seated to the midpoint of its longest edge — the most stable entry. Links carry **no
Z**: they are planar segments with a *tag*; what retracting means in 3D is the
consumer's business.

### Leads — easing on and off the cut

With `lead` on, engagement moves are woven **only around retract steps** (a cut link
already rides down into the next ring). They come paired: a **lead-out** eases the tool
off the previous cut before the retract lifts, and a **lead-in** eases onto the new one
after the rapid lands — the tool always moves *away* from a wall before retracting and
*toward* a cut after arriving, never landing cold on the cut path.

The construction reuses the one operation this library trusts: offsetting. A lead rides
a small **inset of its own target** (`leadInset × R`), morphing between the inset and an
offset of that inset — both members of one offset family, so they share topology and the
in-between never distorts. The lead's arc length is the horizontal run of the ramp the
consumer will descend (`cuttingDepth / tan(plungeSlope)`), so the same chain doubles as
a **helix ramp**: on an island too small for the run, the lead simply wraps the inset in
whole laps. Lead-out runs its own, much steeper angle (`retractSlope`) — it climbs, it
doesn't plunge.

An **`Entrance`** target's lead-in gets one extra rule: its ramp descends into
full-width material, so the run is floored at **one full lap** of the ring — the
descent is spread around the whole loop, never plunged down a short arc.

### The witness-mark overshoot

One last courtesy to the part: an outer boundary ring's start/end seam would leave a
faint witness bump on the finished wall. So every generation-1 ring **overshoots**:
after closing its loop, the tool carries on along the chain by one tool radius,
re-cutting the seam clean. Inner rings never witness a wall and are left alone.

### The degenerate sweep — the rung below extinction

A seed chain that was *given* open air and still left no generation-1 child has measured
its own verdict: the region is narrower than the tool even after the push funded both
sides of the inset. But that verdict is simultaneously a guarantee — a region thinner
than the tool is cleared **entirely by one pass** riding a single radius off its walls,
overhanging into the air. So instead of silently emitting nothing, the strategy insets
the chain's non-open-air runs *blindly* (each wall edge one radius left of travel, gaps
bridged straight) and emits each run as an **open** cut chain, capped by straight
tangent leads — both ends hang in free air, so the ramps descend and climb over
nothing.

Two rules keep it honest. Only chains that *had* air qualify: a childless chain walled
all the way round is genuinely uncuttable by this tool, and stays uncut. And no-touch
stays sovereign — the sweep is *measured* against every wall it isn't itself riding,
and one that would pass within a tool radius of an opposing wall (a thin slot walled on
both sides) is refused with a warning rather than emitted. Like generation 1, sweeps
are a wall pass, so they respect `emitOuter`.

---

## How not to use this layer

- **Don't feed unoriented input.** Doctrine orientation (CCW clears, CW keeps) is the
  caller's contract. The strategy will faithfully toolpath whatever the directions
  *say*; if they're arbitrary, the output is meaningless, not merely wrong-handed.
- **Don't offset open chains, or expect them to cut.** Only closed chains bound
  regions. Open chains in a result are links and leads — the strategy's own weaving.
- **Don't reorient the method's outputs with `reversed()`.** It exists to *condition
  inputs*. Use `climb` for milling sense and `reverse` for execution order.
- **Don't "fix" a toolpath downstream.** The result is consumed blindly; if the motion
  is wrong, the strategy (or the preparation) is where truth lives.
- **Don't infer stay-down hops from ancestry.** A parent's other children can be
  anywhere. Reach and clearance are *measured* per hop; that measurement is the rule.
- **Don't read links as 3D moves.** Links and leads are planar chains with tags. Z
  belongs to the consumer (retract heights, ramp descents).
- **Don't quantise the model.** If you extend a stoicheion, keep the doctrine: exact
  closed forms where they exist, certified numerics on the exact curve where they don't
  — and keep orientation in the points, so degeneracies flip chirality by geometry.
- **Don't treat extinction as an error.** A ring too small to offset produces nothing;
  a profile with nothing left returns empty. That is the recursion *ending correctly*.
