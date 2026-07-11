# Stoicheia and Handedness

## The atom

A **stoicheion** (στοιχεῖον — Euclid's *Elements* is the *Stoicheia*; one element is a
stoicheion) is a single exact curve. The roster is small and closed: `Point`, `Segment`,
`Circle`, `Arc`, `Ellipse`, `EllipseArc`, `OffsetEllipse`, `Spline`. Each knows how to
tessellate itself, measure distance to a probe, intersect another stoicheion, and — the
operation this whole library is built on — **offset itself sideways**.

Everything below is about one decision: **orientation and shape live in the points, not in
flags.** That decision is what makes the offset method degeneracy-proof, and it is the reason a
stoicheion "responds to an offset by its own nature" rather than by a switch statement.

## Points, not bools and per-type variables

The tempting representation is the obvious one: a circle is a centre, a radius (a scalar), and a
`bool clockwise`; an arc is a centre, a radius, a start angle, a sweep, and a direction flag.
Handedness becomes a boolean; the "amount of curve" becomes a signed scalar. This is how most
CAD kernels do it, and it is a trap, because **every flag is a special case waiting to be
handled at every degeneracy.**

We represent the same information *positionally* instead:

- A **circle** is a centre `c`, a **radius point** `a` (so the radius is `|a − c|`, never a
  separate scalar), and a **chirality marker** `b`, a point a quarter-turn from `a`. *Which side
  `b` sits on is the orientation.* There is no `bool clockwise`.
- An **arc** is a centre `c`, two endpoints `a`, `b`, and a **through-point** `d`: the arc runs
  "from `a` to `b` **through** `d`." `d` fixes *which* of the two possible arcs it is **and**
  which way it travels. There is no direction flag and no signed sweep.

The marker point is not a solver variable and is not one of the entity's defining degrees of
freedom — it is a *carried* point, re-seated to its canonical position after any change so it
always reads correctly and rides along when the entity moves. Its whole job is to *hold the
handedness geometrically*.

## Why positional handedness is not a stylistic choice

Consider offsetting a circle inward — insetting it — past its own centre. With a signed-scalar
radius you now have a negative radius, which is nonsense, so you write special-case code:
detect the sign flip, negate the radius, flip the `clockwise` bool, rotate the start angle by π.
Every primitive needs its own version of this, and every one is a place to get it subtly wrong.

With the positional representation there is **nothing to special-case.** Drive the inset through
the centre and the radius point `a` crosses to the far side of `c` on its own; the marker `b`,
kept canonical, lands on the *other* side; and the circle re-emerges inside-out with its
chirality flipped — *by geometry*, because that is what those points do when you move them. The
inversion is not detected and handled; it simply happens, correctly, the way the physical
situation happens.

This is the invariant the whole method leans on:

> **Every downstream rule — which way an offset goes, what a loop's handedness is, whether a
> candidate inverted — reads orientation off the points. So every degeneracy flips chirality the
> same way the geometry does, with no rule to author and no case to miss.**

A boolean cannot do this. A boolean has to be *told* to flip, by code that has to *notice* the
degeneracy first. Points flip because they moved.

## Stoicheometric handedness: why a primitive responds "by its nature"

Because handedness is carried in the points, an individual stoicheion has an **implied
handedness** the moment it is constructed — the side its marker sits on, the side its
through-point sits on. And offsetting is defined uniformly for *all* stoicheia as a single
operation: **move the defining points to the left of travel by the offset amount** (see
[02-Offsetting-Semantics.md](02-Offsetting-Semantics.md)).

The consequence is that the *same* offset operation produces *opposite* geometric motion
depending on the primitive's implied handedness, and it does so with no branch:

- A counter-clockwise circle, offset "left" by a positive amount, **shrinks** — because for
  CCW travel the left side is the interior.
- A clockwise circle, offset "left" by the *same* positive amount, **grows** — because for CW
  travel the left side is the exterior.

Neither circle was asked which way to move. Each moved the way its own point construction
dictates. The offset primitive is blind to handedness; handedness is a property the points
already carry, so the *response* to a blind operation is where handedness expresses itself.

That is what "responds to offset by its own nature" means, and it is the seed of everything in
the next three documents: a chain is a run of such stoicheia, a chain's handedness is the sum of
its edges' implied handedness, and the entire frontier method is the study of what happens when
you apply *one* blind offset to a whole chain and let each edge, and the chain as a whole,
respond according to the handedness it already holds.

## Carried across to the reconstruction

- Orientation and magnitude are **points**, never a `bool` or a signed scalar. No primitive
  carries a `clockwise` flag or a signed radius.
- The chirality marker / through-point is a **carried, canonicalised** point — re-seated after
  every mutation, never solved for, never a degree of freedom.
- Offsetting is **one uniform operation** on the points; per-handedness behaviour is *emergent*,
  never branched. If the reconstruction finds itself writing `if (clockwise)` inside an offset,
  it has left the method.
