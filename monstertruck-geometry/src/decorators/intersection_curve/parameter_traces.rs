//! Bounded, non-geometric continuation hints for an implicit intersection branch.

use super::*;
use serde::{Deserializer, de::Error};

const MAX_TRACE_CONTROLS: usize = 32_768;
type BoxedParameterTraces = Box<[BsplineCurve<Point2>; 2]>;

pub(super) fn validate_tubes(tubes: &[SurfaceParameterTube]) -> bool {
    (1..=8_192).contains(&tubes.len())
        && tubes.iter().all(|tube| {
            tube.parameter_range.0.is_finite()
                && tube.parameter_range.1.is_finite()
                && tube.parameter_range.0 < tube.parameter_range.1
                && tube
                    .uv_bounds
                    .iter()
                    .all(|(low, high)| low.is_finite() && high.is_finite() && low < high)
        })
        && tubes
            .windows(2)
            .all(|pair| pair[0].parameter_range.1 == pair[1].parameter_range.0)
}

pub(super) fn validate(traces: &[BsplineCurve<Point2>; 2]) -> bool {
    let valid = |trace: &BsplineCurve<Point2>| {
        let controls = trace.control_points();
        let knots = trace.knot_vector();
        (2..=MAX_TRACE_CONTROLS).contains(&controls.len())
            && knots.len() == controls.len() + 2
            && controls
                .iter()
                .all(|point| point.x.is_finite() && point.y.is_finite())
            && knots.iter().all(|parameter| parameter.is_finite())
            && knots.windows(2).all(|pair| pair[0] <= pair[1])
            && knots[0] == knots[1]
            && knots[knots.len() - 1] == knots[knots.len() - 2]
            && knots[1] < knots[knots.len() - 2]
            && knots[1..knots.len() - 1]
                .windows(2)
                .all(|pair| pair[0] < pair[1])
    };
    traces.iter().all(valid) && traces[0].range_tuple() == traces[1].range_tuple()
}

pub(in crate::decorators) fn deserialize_parameter_traces<'de, D>(
    deserializer: D,
) -> std::result::Result<Option<BoxedParameterTraces>, D::Error>
where D: Deserializer<'de> {
    let traces = Option::<[BsplineCurve<Point2>; 2]>::deserialize(deserializer)?;
    if traces.as_ref().is_none_or(validate) {
        Ok(traces.map(Box::new))
    } else {
        Err(D::Error::custom(
            "invalid paired intersection continuation traces",
        ))
    }
}

fn align(parameter: Point2, hint: Point2, surface: &impl ParametricSurface3D) -> Option<Point2> {
    let component = |value: f64, seed: f64, period: Option<f64>| match period
        .filter(|period| period.is_finite() && *period > 0.0)
    {
        Some(period) => value + ((seed - value) / period).round() * period,
        None => value,
    };
    let parameter = Point2::new(
        component(parameter.x, hint.x, surface.period_u()),
        component(parameter.y, hint.y, surface.period_v()),
    );
    (parameter.x.is_finite() && parameter.y.is_finite()).then_some(parameter)
}

pub(super) fn align_triple(
    (point, uv0, uv1): (Point3, Point2, Point2),
    hints: Option<[Point2; 2]>,
    surface0: &impl ParametricSurface3D,
    surface1: &impl ParametricSurface3D,
) -> Option<(Point3, Point2, Point2)> {
    let (uv0, uv1) = match hints {
        Some([hint0, hint1]) => (align(uv0, hint0, surface0)?, align(uv1, hint1, surface1)?),
        None => (uv0, uv1),
    };
    let (p0, p1) = (surface0.subs(uv0.x, uv0.y), surface1.subs(uv1.x, uv1.y));
    (finite_point(p0) && finite_point(p1) && point.near(&p0) && point.near(&p1))
        .then_some((point, uv0, uv1))
}
