//! Bounded display tessellation of the evaluated intersection, never its leader.

use super::*;

const MAX_DEPTH: usize = 24;
const MAX_NODES: usize = 16_384;

#[derive(Clone, Copy)]
struct Span {
    range: (f64, f64),
    ends: (Point3, Point3),
    depth: usize,
}

pub(super) fn try_divide(
    curve: &impl ParametricCurve3D,
    range: (f64, f64),
    tolerance: f64,
) -> Option<(Vec<f64>, Vec<Point3>)> {
    if !tolerance.is_finite()
        || tolerance < TOLERANCE
        || !range.0.is_finite()
        || !range.1.is_finite()
        || range.0 > range.1
    {
        return None;
    }
    let ends = (curve.evaluate(range.0), curve.evaluate(range.1));
    if !finite_point(ends.0) || !finite_point(ends.1) {
        return None;
    }
    if range.0 == range.1 {
        return Some((vec![range.0], vec![ends.0]));
    }
    let mut pending = vec![Span {
        range,
        ends,
        depth: 0,
    }];
    let mut parameters = vec![range.0];
    let mut points = vec![ends.0];
    let mut visited = 0;
    while let Some(span) = pending.pop() {
        visited += 1;
        if visited > MAX_NODES {
            return None;
        }
        let parameter = |fraction| span.range.0 * (1.0 - fraction) + span.range.1 * fraction;
        let middle_parameter = parameter(0.5);
        let middle = curve.evaluate(middle_parameter);
        if !finite_point(middle) {
            return None;
        }
        let samples = [
            (0.25, curve.evaluate(parameter(0.25))),
            (0.5, middle),
            (0.75, curve.evaluate(parameter(0.75))),
        ];
        if samples.iter().any(|(_, point)| !finite_point(*point)) {
            return None;
        }
        let flat = samples.iter().all(|(fraction, point)| {
            let chord = span.ends.0 + (span.ends.1 - span.ends.0) * *fraction;
            point.distance2(chord) <= tolerance * tolerance
        });
        if flat {
            parameters.push(span.range.1);
            points.push(span.ends.1);
        } else {
            if span.depth == MAX_DEPTH
                || middle_parameter <= span.range.0
                || middle_parameter >= span.range.1
            {
                return None;
            }
            pending.extend([
                Span {
                    range: (middle_parameter, span.range.1),
                    ends: (middle, span.ends.1),
                    depth: span.depth + 1,
                },
                Span {
                    range: (span.range.0, middle_parameter),
                    ends: (span.ends.0, middle),
                    depth: span.depth + 1,
                },
            ]);
        }
    }
    Some((parameters, points))
}
