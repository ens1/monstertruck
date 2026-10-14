//! Algebraic lifting of rational parameter curves, without sampling a 3D guide.

use super::*;

const MAX_DEGREE: usize = 32;
const MAX_SPANS: usize = 256;
const MAX_CONTROLS: usize = 16_384;
const MAX_PROOF_WORK: usize = 8_000_000;

fn finite3(p: Vector3) -> bool { [p.x, p.y, p.z].into_iter().all(f64::is_finite) }
fn finite4(p: Vector4) -> bool { [p.x, p.y, p.z, p.w].into_iter().all(f64::is_finite) }

fn homogeneous(curve: &Curve2D) -> Option<BsplineCurve<Vector3>> {
    let result = match curve {
        Curve2D::Line(Line(a, b)) => BsplineCurve::new(
            KnotVector::from(vec![0.0, 0.0, 1.0, 1.0]),
            vec![Vector3::new(a.x, a.y, 1.0), Vector3::new(b.x, b.y, 1.0)],
        ),
        Curve2D::Polyline(points) if points.len() >= 2 && points.len() <= MAX_SPANS => {
            let end = (points.len() - 1) as f64;
            BsplineCurve::new(
                KnotVector::from(
                    iter::once(0.0)
                        .chain((0..points.len()).map(|i| i as f64))
                        .chain(iter::once(end))
                        .collect::<Vec<_>>(),
                ),
                points.iter().map(|p| Vector3::new(p.x, p.y, 1.0)).collect(),
            )
        }
        Curve2D::BsplineCurve(curve) => BsplineCurve::new(
            curve.knot_vector().clone(),
            curve
                .control_points()
                .iter()
                .map(|p| Vector3::new(p.x, p.y, 1.0))
                .collect(),
        ),
        Curve2D::NurbsCurve(curve) => curve.non_rationalized().clone(),
        _ => return None,
    };
    (result.degree() <= MAX_DEGREE
        && result.control_points().len() <= MAX_CONTROLS
        && result.knot_vector().iter().all(|k| k.is_finite())
        && result
            .control_points()
            .iter()
            .all(|p| finite3(*p) && p.z > 0.0))
    .then_some(result)
}

mod bernstein;
mod certified;
mod interval;

impl Curve {
    /// Lifts an authoritative parameter curve and bounds its full world-space error.
    ///
    /// The bound includes original knot extraction, algebraic composition,
    /// denominator conditioning, emitted control rounding and shared-piece joins.
    /// Supports positive rational planar curves and clamped rational tensor-product
    /// cells. Uncertain knot-cell crossings, nonpositive interval denominators,
    /// unsupported native UV frames and exceeded work limits return `None`.
    /// Parameter interval and direction are preserved. A search guide is never
    /// used as replacement geometry. Callers must compare the returned distance
    /// against their explicit modeling residual before using the lifted curve.
    pub fn try_lift_parameter_curve_with_error(&self) -> Option<(Curve, f64)> {
        certified::lift(self)
    }

    /// Lifts a parameter curve only when its certified error fits native precision.
    ///
    /// Use [`Self::try_lift_parameter_curve_with_error`] when the caller owns an
    /// explicit modeling residual. This convenience method requires the complete
    /// bound to be at most the native `TOLERANCE`; it does not return an unchecked
    /// algebraic approximation.
    pub fn try_lift_parameter_curve(&self) -> Option<Curve> {
        self.try_lift_parameter_curve_with_error()
            .filter(|(_, error)| *error <= TOLERANCE)
            .map(|(curve, _)| curve)
    }

    /// Bounds distance at corresponding normalized parameters over the full curves.
    ///
    /// Original knot extraction, interval subdivision and optional parameter-curve
    /// lifting errors are included. An upper bound larger than a caller's residual
    /// does not prove different loci: the curves may have different traversal
    /// speeds, or the interval certificate may be too conservative. Unsupported
    /// representations or unresolved normalized knot domains return `None`.
    pub fn try_parameter_distance_bound(&self, other: &Self) -> Option<f64> {
        certified::distance(self, other)
    }
}
