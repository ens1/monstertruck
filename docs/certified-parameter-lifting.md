# Certified algebraic parameter-curve lifting

`Curve::try_lift_parameter_curve_with_error` returns a rational 3D curve and a
complete world-distance bound against the original parameter curve on its
original native carrier. The caller owns the permitted residual and must check
the returned bound before admitting or exporting the replacement. The original
parameter interval and direction are preserved. No polyline search guide,
sampled fit or post-hoc agreement of sampled points supplies the certificate.

The supported inputs are clamped positive rational parameter curves on planes,
or clamped positive rational tensor-product spline cells. Planes use their
original three defining points; outward subtraction retains axis arithmetic
error. Polynomial surface nets retain their native UV frame. Analytic angular
surfaces, non-rational UV maps and unresolved surface-knot crossings refuse.

Rational lifting starts from the original source controls and knot values.
Outward knot insertion extracts the native Bezier pieces without treating
rounded library subdivision outputs as exact inputs. Positive rational UV
control bounds must fit a single original surface knot cell. Degree-one UV
pieces can be split across cells when the outward original knot-plane root
identifies an exactly representable native parameter. An uncertain crossing
cannot select a neighboring polynomial patch by a rounded containment vote.

Tensor-product Bernstein composition uses outward arithmetic for coordinate
normalization, basis powers, products and coefficient sums. The emitted
positive denominator is bounded away from zero. The final world-distance
bound includes homogeneous coefficient rounding, numerator/denominator
conditioning and any shared-piece homogeneous scaling and endpoint join.
The result may exceed a caller's permitted residual even when the algebraic
family is supported; this is a conservative refusal rather than a fitted
replacement or an assertion that the original geometry is invalid.

`Curve::try_lift_parameter_curve` is a convenience wrapper requiring the bound
to fit the native `TOLERANCE`. It does not return an unchecked algebraic lift.
`Curve::try_parameter_distance_bound` compares corresponding normalized
parameters. It includes each parameter lift's error, original knot extraction
and outward interval subdivision in the cross-multiplied Bernstein bound.
A large upper bound does not establish different geometric loci: two valid
traversals can use different speeds, and the enclosure may be conservative.
The caller must not label such an upper bound as proof of invalid geometry.
Unresolved normalized knot domains return `None`.

The implementation limits degree to 32, pieces/cells to 256 and native controls
to 16,384. Deterministic polynomial work estimates additionally bound source
knot extraction, composition and synchronized distance calculations. Nonfinite
arithmetic, loss of a positive interval denominator, discontinuous internal
knots, unclamped input or exhausted bounds return `None`.

Implementation verification uses compilation, formatting, static review and
Clippy. No tests, assertion probes, example executions or manual smoke checks
were performed. These checks do not establish interactive modeling acceptance.
