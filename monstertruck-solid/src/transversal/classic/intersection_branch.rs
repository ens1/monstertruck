//! Temporary intersection geometry with exact producing face-use identities.

use super::{ShapeOpsSurface, finite_point};
use monstertruck_geometry::prelude::*;
use monstertruck_topology::FaceId;

pub(super) type NativeIntersection<S> = SurfaceCurve<
    BsplineCurve<Point3>,
    S,
    S,
    ParameterCurve<BoundaryCurve2D, S>,
    ParameterCurve<BoundaryCurve2D, S>,
>;

#[derive(Clone, Debug)]
pub(super) struct IntersectionBranch<S> {
    pub(super) curve: NativeIntersection<S>,
    source_faces: Option<[FaceId<S>; 2]>,
}

impl<S: ShapeOpsSurface> IntersectionBranch<S> {
    pub(super) fn new(curve: NativeIntersection<S>) -> Self {
        Self {
            curve,
            source_faces: None,
        }
    }
    pub(super) fn bind(&mut self, faces: [FaceId<S>; 2]) { self.source_faces = Some(faces); }
    fn side(&self, face: FaceId<S>) -> Option<usize> {
        self.source_faces?.iter().position(|source| *source == face)
    }
    pub(super) fn parameter_on(&self, face: FaceId<S>, t: f64) -> Option<Point2> {
        let faces = self.source_faces?;
        let (_, uv0, uv1) = self.curve.search_triple(t, 100)?;
        if faces[0] == face {
            Some(uv0)
        } else if faces[1] == face {
            Some(uv1)
        } else {
            None
        }
    }
    pub(super) fn parameter_samples(&self, face: FaceId<S>, tolerance: f64) -> Option<Vec<Point2>> {
        let side = self.side(face)?;
        let surface = if side == 0 {
            self.curve.surface0()
        } else {
            self.curve.surface1()
        };
        let range = self.range_tuple();
        let mut knots = vec![range.0, range.1];
        knots.extend(
            self.curve.parameter_traces()?[side]
                .knot_vector()
                .iter()
                .copied()
                .filter(|t| *t > range.0 && *t < range.1),
        );
        knots.sort_by(f64::total_cmp);
        knots.dedup();
        let mut output = Vec::with_capacity(knots.len());
        let mut work = 0;
        for span in knots.windows(2) {
            let mut pending = vec![(
                span[0],
                span[1],
                self.parameter_on(face, span[0])?,
                self.parameter_on(face, span[1])?,
                0,
            )];
            if output.is_empty() {
                output.push(self.parameter_on(face, span[0])?);
            }
            while let Some((a, b, uv0, uv1, depth)) = pending.pop() {
                work += 1;
                if work > 65_536 {
                    return None;
                }
                let mut flat = true;
                let mut middle = None;
                for fraction in [0.25, 0.5, 0.75] {
                    let t = a * (1.0 - fraction) + b * fraction;
                    let (point, p0, p1) = self.curve.search_triple(t, 100)?;
                    let uv = if side == 0 { p0 } else { p1 };
                    let chord_uv = uv0 + (uv1 - uv0) * fraction;
                    let chord_point = surface.subs(chord_uv.x, chord_uv.y);
                    if !finite_point(chord_point) {
                        return None;
                    }
                    flat &= point.distance(chord_point) <= tolerance;
                    if fraction == 0.5 {
                        middle = Some((t, uv));
                    }
                }
                if flat {
                    output.push(uv1);
                } else {
                    let (t, uv) = middle?;
                    if depth == 24 || t <= a || t >= b {
                        return None;
                    }
                    pending.extend([(t, b, uv, uv1, depth + 1), (a, t, uv0, uv, depth + 1)]);
                }
            }
        }
        Some(output)
    }
    pub(super) fn snap_endpoint(&mut self, front: bool, point: Point3) -> Option<()> {
        let mut curve = self.curve.clone();
        let mut traces = curve.parameter_traces()?.clone();
        let range = curve.range_tuple();
        let t = if front { range.0 } else { range.1 };
        let uv = |surface: &S, hint: Point2| {
            let parameter = surface.search_parameter(point, Some(hint.into()), 100)?;
            let axis = |value: f64, seed: f64, period: Option<f64>| match period
                .filter(|period| period.is_finite() && *period > 0.0)
            {
                Some(period) => value + ((seed - value) / period).round() * period,
                None => value,
            };
            let parameter = Point2::new(
                axis(parameter.0, hint.x, surface.period_u()),
                axis(parameter.1, hint.y, surface.period_v()),
            );
            let on_surface = surface.subs(parameter.x, parameter.y);
            (finite_point(on_surface) && point.near(&on_surface)).then_some(parameter)
        };
        let values = [
            uv(curve.surface0(), traces[0].subs(t))?,
            uv(curve.surface1(), traces[1].subs(t))?,
        ];
        let index = if front {
            0
        } else {
            curve.leader().control_points().len() - 1
        };
        *curve.leader_mut().control_point_mut(index) = point;
        for (trace, value) in traces.iter_mut().zip(values) {
            let index = if front {
                0
            } else {
                trace.control_points().len() - 1
            };
            *trace.control_point_mut(index) = value;
        }
        if let Some(boundary) = curve.boundary0_mut() {
            *boundary = ParameterCurve::new(
                BoundaryCurve2D::BsplineCurve(traces[0].clone()),
                boundary.surface().clone(),
            );
        }
        if let Some(boundary) = curve.boundary1_mut() {
            *boundary = ParameterCurve::new(
                BoundaryCurve2D::BsplineCurve(traces[1].clone()),
                boundary.surface().clone(),
            );
        }
        self.curve = curve.with_parameter_traces(traces)?;
        Some(())
    }
}

pub(super) fn trace<S: ShapeOpsSurface>(
    surface0: S,
    surface1: S,
    points: Vec<Point3>,
    tolerance: f64,
) -> Option<IntersectionBranch<S>> {
    if points.len() < 2 || points.len().checked_mul(2)? > 32_768 {
        return None;
    }
    let denominator = (points.len() - 1) as f64;
    let knot_vector = |parameters: &[f64]| {
        KnotVector::from(
            std::iter::once(parameters[0])
                .chain(parameters.iter().copied())
                .chain(std::iter::once(parameters[parameters.len() - 1]))
                .collect::<Vec<_>>(),
        )
    };
    let parameters = (0..points.len())
        .map(|index| index as f64 / denominator)
        .collect::<Vec<_>>();
    let raw = BsplineCurve::new(knot_vector(&parameters), points);
    let intersection = IntersectionCurve::new(surface0.clone(), surface1.clone(), raw);
    let mut recovered: Vec<(Point3, Point2, Point2)> = Vec::with_capacity(parameters.len() * 2 - 1);
    let mut trace_parameters = Vec::with_capacity(recovered.capacity());
    let unwrap = |uv: Point2, previous: Point2, surface: &S| {
        let axis = |value: f64, hint: f64, period: Option<f64>| match period {
            Some(period) if period.is_finite() && period > 0.0 => {
                value + ((hint - value) / period).round() * period
            }
            _ => value,
        };
        Point2::new(
            axis(uv.x, previous.x, surface.period_u()),
            axis(uv.y, previous.y, surface.period_v()),
        )
    };
    let values = parameters
        .iter()
        .enumerate()
        .flat_map(|(index, parameter)| {
            std::iter::once(*parameter).chain(
                parameters
                    .get(index + 1)
                    .map(|next| (parameter + next) * 0.5),
            )
        });
    for t in values {
        let seed = intersection.leader().subs(t);
        let hint0 = recovered.last().map(|triple| triple.1.into());
        let hint1 = recovered.last().map(|triple| triple.2.into());
        let (point, mut uv0, mut uv1) = intersection
            .search_triple_with_hints(t, hint0, hint1, 100)
            .or_else(|| intersection.search_nearest_point(seed, hint0, hint1, 100))?;
        if let Some(previous) = recovered.last() {
            uv0 = unwrap(uv0, previous.1, &surface0);
            uv1 = unwrap(uv1, previous.2, &surface1);
        }
        let q0 = surface0.subs(uv0.x, uv0.y);
        let q1 = surface1.subs(uv1.x, uv1.y);
        if !finite_point(point)
            || !finite_point(q0)
            || !finite_point(q1)
            || !uv0.x.is_finite()
            || !uv0.y.is_finite()
            || !uv1.x.is_finite()
            || !uv1.y.is_finite()
            || point.distance(q0) > TOLERANCE
            || point.distance(q1) > TOLERANCE
            || q0.distance(q1) > TOLERANCE
            || point.distance(seed) > 4.0 * tolerance
            || recovered
                .last()
                .is_some_and(|previous| previous.0.distance(point) <= f64::EPSILON)
        {
            return None;
        }
        recovered.push((point, uv0, uv1));
        trace_parameters.push(t);
    }
    let knots = knot_vector(&trace_parameters);
    let leader = BsplineCurve::new(
        knots.clone(),
        recovered.iter().map(|triple| triple.0).collect(),
    );
    let traces = [
        BsplineCurve::new(
            knots.clone(),
            recovered.iter().map(|triple| triple.1).collect(),
        ),
        BsplineCurve::new(knots, recovered.iter().map(|triple| triple.2).collect()),
    ];
    // An analytic plane with no exact patch-domain parameterization is affine
    // in its own UV frame. Both affine carriers make these degree-one pcurves
    // actual boundaries. General planar spline carriers can have nonlinear UV
    // frames and deliberately do not enter this constructor.
    let affine_plane = |surface: &S| {
        !surface.supports_exact_patch_domains()
            && matches!(
                surface.try_into_analytic_surface_kind(),
                Some(AnalyticSurfaceKind::Plane(_))
            )
    };
    let boundaries = if affine_plane(&surface0) && affine_plane(&surface1) {
        [
            Some(ParameterCurve::new(
                BoundaryCurve2D::BsplineCurve(traces[0].clone()),
                surface0.clone(),
            )),
            Some(ParameterCurve::new(
                BoundaryCurve2D::BsplineCurve(traces[1].clone()),
                surface1.clone(),
            )),
        ]
    } else {
        [None, None]
    };
    let [boundary0, boundary1] = boundaries;
    let curve = SurfaceCurve::with_boundaries(surface0, surface1, leader, boundary0, boundary1)
        .with_parameter_traces(traces)?;
    Some(IntersectionBranch::new(curve))
}

impl<S: ShapeOpsSurface> ParametricCurve for IntersectionBranch<S> {
    type Point = Point3;
    type Vector = Vector3;
    fn evaluate(&self, t: f64) -> Point3 { self.curve.evaluate(t) }
    fn derivative(&self, t: f64) -> Vector3 { self.curve.derivative(t) }
    fn derivative_2(&self, t: f64) -> Vector3 { self.curve.derivative_2(t) }
    fn derivative_n(&self, n: usize, t: f64) -> Vector3 { self.curve.derivative_n(n, t) }
    fn derivatives(&self, n: usize, t: f64) -> CurveDerivatives<Vector3> {
        self.curve.derivatives(n, t)
    }
    fn parameter_range(&self) -> ParameterRange { self.curve.parameter_range() }
}

impl<S: ShapeOpsSurface> BoundedCurve for IntersectionBranch<S> {}

impl<S: ShapeOpsSurface> Cut for IntersectionBranch<S> {
    fn cut(&mut self, t: f64) -> Self {
        Self {
            curve: self.curve.cut(t),
            source_faces: self.source_faces,
        }
    }
}

impl<S: ShapeOpsSurface> Invertible for IntersectionBranch<S> {
    fn invert(&mut self) { self.curve.invert(); }
}

impl<S: ShapeOpsSurface> ParameterDivision1D for IntersectionBranch<S> {
    type Point = Point3;
    fn try_parameter_division(
        &self,
        range: (f64, f64),
        tol: f64,
    ) -> Option<(Vec<f64>, Vec<Point3>)> {
        self.curve.try_parameter_division(range, tol)
    }
    fn parameter_division(&self, range: (f64, f64), tol: f64) -> (Vec<f64>, Vec<Point3>) {
        self.curve.parameter_division(range, tol)
    }
}

impl<S: ShapeOpsSurface> SearchParameter<CurveParameter> for IntersectionBranch<S> {
    type Point = Point3;
    fn search_parameter<H: Into<SearchParameterHint1D>>(
        &self,
        point: Point3,
        hint: H,
        trials: usize,
    ) -> Option<f64> {
        self.curve.search_parameter(point, hint, trials)
    }
}

impl<S: ShapeOpsSurface> SearchNearestParameter<CurveParameter> for IntersectionBranch<S> {
    type Point = Point3;
    fn search_nearest_parameter<H: Into<SearchParameterHint1D>>(
        &self,
        point: Point3,
        hint: H,
        trials: usize,
    ) -> Option<f64> {
        self.curve.search_nearest_parameter(point, hint, trials)
    }
}

impl<S: ShapeOpsSurface> ParameterBoundary2D<S> for IntersectionBranch<S> {
    fn parameter_boundary_2d(&self, surface: &S, tolerance: f64) -> Option<Vec<Point2>> {
        let parameters = self
            .curve
            .try_parameter_division(self.range_tuple(), tolerance)?
            .0;
        // This trait feeds search tessellation, not native face construction.
        // Native division uses `parameter_on` and producing face identity instead.
        let side = |first: bool| {
            parameters
                .iter()
                .map(|t| {
                    let (point, uv0, uv1) = self.curve.search_triple(*t, 100)?;
                    let uv = if first { uv0 } else { uv1 };
                    let on_surface = surface.subs(uv.x, uv.y);
                    (finite_point(on_surface) && point.near(&on_surface)).then_some(uv)
                })
                .collect::<Option<Vec<_>>>()
        };
        match (side(true), side(false)) {
            (Some(first), Some(second)) if first == second => Some(first),
            (Some(_), Some(_)) => None,
            (first, second) => first.or(second),
        }
    }
}
