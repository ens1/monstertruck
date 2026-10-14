# Native B-spline cut parameter precision

The native curve `Cut` implementation previously applied the shared geometric
`TOLERANCE` (a model length) to dimensionless parameters. A requested cut close
to a knot or domain end could silently move to that knot. A fast or large curve
could therefore move its boundary by much more than the caller's world-space
residual budget. Knot insertion also treated small positive parameter spans as
zero through `inv_or_zero`, so simply deleting the final snap was insufficient.

`BsplineCurve::try_cut` now retains the requested finite `f64` value literally.
Its left and right active domains meet at exactly that stored parameter. Exact
knot multiplicity, exact domain-end comparisons, direct dimensionless insertion
ratios, and component interpolation replace the geometric near/zero decisions
in the split path. The coefficient interpolation applies to generic Euclidean
and homogeneous controls without projecting or resampling the curve. Large
opposite finite controls do not require an overflowing difference to be formed.
Finite knot intervals that overflow during subtraction use a common scaled
ratio; tiny positive intervals use direct division rather than an overflowing
reciprocal.

`try_cut` validates source shape, finite knots and controls, sorted knots, the
active range `[knot[degree], knot[control_count]]`, and complete literal repeated
knot blocks. Overfull interior knot blocks and unsupported collapsed source
domains are refused. Inactive endpoint padding from native derivative nets or
endpoint knot insertion remains supported; it does not change the active domain. All insertion and output checks happen on a clone before the source is
replaced. A failure returns `BsplineCutError` and leaves the source unchanged.
Cuts outside the active range now fail instead of synthesizing exterior origin
pieces. The existing infallible `Cut` trait retains its signature and explicitly
panics with that error for unsupported inputs; it does not publish a snapped or
unchanged substitute as successful output.

At an exact active-domain end, the empty side is a degree-zero boundary sentinel
with repeated equal knots and the actual evaluated homogeneous boundary control.
Native curve evaluation explicitly returns that stored control for a zero-width
range. Native surface evaluation similarly retains the true stored boundary
row/column when one axis collapses. These sentinel sides do not claim a nonzero
CAD edge or surface area.

`BsplineSurface::try_cut_u` and `try_cut_v` split all native control columns or rows
through the same curve contract. A surface mutates only after every control
curve succeeds. Unclamped endpoint cuts evaluate the actual source boundary
rather than copying a first or last row that may not lie on that boundary.
The infallible surface split methods explicitly panic on a rejected split.
Surface knot insertion uses the same direct ratios and component interpolation.

The B-spline basis value recurrence now also forms direct ratios, so successful
subnormal positive-span position cuts do not immediately evaluate as
`infinity * zero`. Explicit derivative-curve and derivative-surface construction
uses only exact-zero spans and divides each coefficient difference directly by
its span; valid narrow spans no longer acquire zero tangents merely because of
a geometric tolerance. Truly unrepresentable derivatives can still be nonfinite.

Clamping inspection and `BsplineCurve::clamp` use literal multiplicity. A nearby
but distinct endpoint knot does not qualify a source for Bernstein proofs.
Legacy `sectional_curve` now uses the correct axis components and exact active
range when deciding whether to trim. Its broader degree/knot synchronization,
normalization, and diagonal construction remain legacy utilities; this change
does not certify that helper or its parameterization.

This correction bounds the previous parameter relocation to **zero**. Ordinary
native coefficient arithmetic still rounds, and this API is not a certificate
of world-space locus error, rational denominator positivity, or full derivative
conditioning. Geometry admission, fresh continuation tubes, certified parameter
lifting, and final boundary/vertex residual checks remain separate requirements.
Generic construction/normalization helpers outside this split dependency path
retain their existing policies.

Validation performed: targeted Rust formatting, scoped `git diff --check`, and
`cargo clippy -p monstertruck-geometry --lib -- -D warnings` through a workspace
using the local fork dependency. The lint command compiled the geometry library.
No code tests, geometry probes, examples, or interactive smoke checks were run.
