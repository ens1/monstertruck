# Native intersection branches and bounded continuation

The classic Boolean backend now retains both face-local UV traces alongside the
native implicit surface intersection. Splitting retains the producing face
identities and the leader parameter frame. Face division uses the producing
carrier's native UV evaluations, rather than projecting a sampled 3D polyline
onto an arbitrary carrier. Exact affine-plane boundaries remain available;
interpolated traces on nonlinear carriers are hints, never authoritative pcurves.

`certify_intersection_continuation_with_halfspaces` admits a **discovered** branch
over every complete leader interval using parameterized Krawczyk contractions.
It consumes the original native carrier controls and knots with outward interval
arithmetic. Supported maps are affine planes and continuous clamped polynomial
or positive-weight rational splines of degree 1–16. The guide must be a continuous
degree-one spline in the paired traces' parameter frame. Angular analytic maps,
singular junctions, unsupported guides and ill-conditioned proofs are reported
explicitly. The convenience `certify_intersection_continuation` supplies no
additional retaining halfspaces.

The result establishes a unique intersection/cross-section root in each supplied
UV tube for every parameter of its span. Successive one-sided root enclosures
have a bounded gap at the requested model-length precision. Up to 64 retaining
halfspaces are checked by native scalar enclosures over whole spans with slack
no greater than the smaller of that precision and native `TOLERANCE`. Carrier
control-point projection preserves useful plane-residual correlations instead
of relying on independent 3D coordinate boxes. Inconclusive retaining planes
subdivide under the same bounded proof policy.

The consumer installs the returned tubes on the same `SurfaceCurve`. Evaluation
rejects a final native projection outside all applicable tubes or outside its
retaining halfspaces. Newton iterations are not constrained to a tube, so the
certificate does not guarantee numerical convergence. A fallible projection
failure stays a failure; a sampled leader is not substituted as boundary geometry.

Tubes are transient admission data, omitted from serde and geometric content
hashes. Traces and retaining halfspaces remain part of the native definition.
Transforms and mutable leader/carrier access invalidate tubes. Faithful cuts
restrict them to the actual retained parameter ranges; differently parameterized
cuts invalidate them. Inversion reverses their parameter intervals. Consumers
requiring admitted geometry must reprove and reinstall restrictions after source
changes or loading, using geometric precision separately from display/search
tolerance. The classic branch constructor invokes this admission for supported
carriers at native geometric precision. Unsupported angular maps continue to use
their existing fallible native checks and are not described as certified.

All work is bounded: 32,768 guide/trace controls, 65,536 surface controls, 8,192
subdivision nodes/tubes, depth 24, and 64 million interval-work units. Exhausting
a limit returns an explicit resource or unresolved error. An unresolved proof is
not evidence of invalid geometry or nonexistence.

## Remaining foundation work

This proves local continuation of discovered branches, not complete branch
discovery, exact global continuity, global injectivity, self-separation, or
complete Boolean classification. Mesh interference still seeds discovery, and
mesh absence does not prove native intersection absence. A complete replacement
requires bounded **per-use trim domains**, coverage/root-exclusion over those
domains, contact/tangency handling, and native inside/outside classification.
Natural carrier domains alone cannot establish absence for an arbitrarily
trimmed face. The existing trimmed-face/shell/solid lane should supply that
domain authority; another parallel sampled trim representation would not solve
the gap.

## Verification for this checkpoint

Only library compilation, formatting and strict linting are performed. No code
tests, assertion probes or runtime smoke checks are introduced or executed.
