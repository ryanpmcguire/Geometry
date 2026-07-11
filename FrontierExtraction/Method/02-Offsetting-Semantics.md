# Offsetting Semantics

## What offsetting is

To **offset** a curve is to produce its *parallel curve* at a signed distance: every point
displaced sideways by the same amount along the local normal. For a straight segment the offset
is a parallel segment; for a circle, a concentric circle; for an arc, a concentric arc. For a
whole **chain** — an ordered run of stoicheia joined end to start — offsetting slides every edge
sideways and then reconciles what happens at the corners.

The method's offset is deliberately **blind**: every edge slides to the **left of its own
travel** by the signed amount, and *nothing else is consulted* — no winding test, no "which side
is inside," no preferred direction. This blindness is not laziness; it is the source of the
method's correctness, because it lets handedness — which the edges already carry in their points
([01-Stoicheia-And-Handedness.md](01-Stoicheia-And-Handedness.md)) — decide the meaning of the
operation instead of a fragile external judgement.

## The convention, stated once

Fix the vocabulary the rest of the documents use:

- **`offsetLeft(edge, amount)`** displaces the edge to the **left of travel** by `amount`. A
  positive `amount` moves left; a negative `amount` moves right. This is the *only* offset
  primitive; everything else is built from it.
- **Handedness** `h` of a closed chain is `+1` for counter-clockwise (CCW) travel, `−1` for
  clockwise (CW). It is read off the geometry (the exact signed area, Green's theorem — never a
  flag), so it is available for any chain at any time.
- **Material is on the left of travel.** This is the doctrine every consumer agrees to before
  the method runs. It means a CCW loop bounds a region of material (its interior is on its left);
  a CW loop bounds a *keep-island* — a hole of material sitting inside a clearing (its interior
  is on its right).

## Why CCW and CW respond oppositely to the same value

This is the crux of the whole system, and it follows from the convention with no further
machinery. `offsetLeft` moves an edge toward the chain's **left**. Whether the chain's *interior*
lies to its left depends entirely on handedness:

- For a **CCW** chain (`h = +1`), the interior is on the left. So a **positive** amount moves the
  boundary *toward* the interior — it **insets** (shrinks the region).
- For a **CW** chain (`h = −1`), the interior is on the right. So the very same **positive**
  amount moves the boundary *away* from the interior — it **outsets** (grows the region).

One operation, one sign, opposite geometric effect. Nobody decided this per chain; it fell out of
"move left" meeting two different handednesses. This is why the offset is called
**chirality-emergent**: inset-versus-outset is *never inferred* — it *emerges* from the chain's
own orientation.

There is a single scalar that captures it. The chain's **own region shrinks** exactly when the
offset moves toward its interior, i.e.

> **own region shrinks ⟺ `amount · h > 0`.**

Check it: CCW (`h=+1`) with positive amount → `amount·h > 0` → shrinks ✓. CW (`h=−1`) with
positive amount → `amount·h < 0` → grows ✓. This product — call it the offset's **sense** —
is the variable that governs frontier extraction and the inversion law in the next two documents.
It is *not* the raw sign of the amount; it is the sign of the amount **relative to the chain's
handedness**, and conflating the two is the single most consequential mistake the reconstruction
can make.

## The unified meaning: advancing the internal frontier

There is a reason the method applies *one* sign to a whole mixed system of CCW outers and CW
islands. At a positive amount, the CCW outer *insets* and the CW island *grows* — and **both of
those advance the internal frontier of the clearable material.** Insetting the outer wall and
growing the keep-island both shrink the material that remains to be cleared. So a single positive
offset marches the whole system's internal frontier inward, regardless of the mix of
handednesses. That is why recursive insetting is "just a loop": one signed amount, applied
blindly, means the same *systemic* thing to every chain even as it means opposite *geometric*
things to CCW and CW edges.

The negative direction is the mirror: a negative amount grows the CCW outer and shrinks the CW
island — both **retreat** the internal frontier, i.e. grow the clearable region outward. The
existing library never marches this way; it is a first-class direction in the reconstruction, and
the first place the old assumptions must be unlearned.

## Corners: the join known by construction

Sliding every edge left opens gaps and overlaps at the corners. The method does **not** solve
corners by mitre/clip casework. At every corner it inserts the **join arc known by
construction**: centred on the *original* (pre-offset) corner, radius `|amount|`, continuing the
first edge's travel into the second. A diverging corner takes the short way round (a rounded
join); a converging corner is forced the long way round (a loop that exits and returns). No
decision is made anywhere — the join is *built*, not chosen.

The price of blindness is that the raw offset **self-intersects**: converging corners produce
little excursion loops, and distant parts of the chain that the offset pushed into each other
cross. Those self-intersections are not a defect to avoid; they are the **raw material** from
which the true frontier is extracted. That extraction is the subject of
[03-The-Frontier.md](03-The-Frontier.md).

## The offset families are closed

One more property the representation buys, and the reconstruction must preserve: offsetting never
degrades the *kind* of a curve. A line offsets to a line, a circle to a concentric circle, and
the parametric curves whose true offset is not their own kind (the offset-ellipse, the spline)
absorb the offset into a stored signed distance `d` rather than tessellating. So a circle stays a
single circle through fifty generations and a spline stays one spline — recursion to extinction
accumulates **no tessellation error**. This is why the frontier method can be exact: the geometry
is never quantised, only the arithmetic is approximated where a closed form provably does not
exist.

## Carried across to the reconstruction

- One offset primitive, `offsetLeft`, moving toward the chain's left; sign moves left/right.
- The governing scalar for meaning is the **sense** `sign(amount · h)`, not `sign(amount)`.
- Corners are the **construction join** (original corner, radius `|amount|`), never mitre/clip.
- The raw offset is *expected* to self-intersect; that is input to extraction, not an error.
- Both offset directions are first-class. The negative (frontier-retreating / region-growing)
  direction is the one the old code never exercised and the reconstruction must build honestly.
