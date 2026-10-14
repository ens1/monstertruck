//! Certified algebraic lifting from original native coefficients.
use super::{
    bernstein::{self, Bezier},
    interval::Interval,
    *,
};

fn exact<const N: usize>(values: [f64; N]) -> Option<[Interval; N]> {
    values
        .into_iter()
        .map(Interval::exact)
        .collect::<Option<Vec<_>>>()?
        .try_into()
        .ok()
}
fn vector(point: Vector4) -> [f64; 4] { [point.x, point.y, point.z, point.w] }
fn choose_controls(controls: &[[Interval; 4]]) -> Option<Vec<[f64; 4]>> {
    controls
        .iter()
        .map(|control| {
            control
                .iter()
                .copied()
                .map(Interval::midpoint)
                .collect::<Option<Vec<_>>>()?
                .try_into()
                .ok()
        })
        .collect()
}
fn native(controls: Vec<[f64; 4]>, knots: Vec<f64>) -> Option<Curve> {
    let points = controls
        .into_iter()
        .map(|p| Vector4::new(p[0], p[1], p[2], p[3]))
        .collect();
    BsplineCurve::try_new(KnotVector::from(knots), points)
        .ok()
        .map(|curve| Curve::NurbsCurve(NurbsCurve::new(curve)))
}
fn planar(uv: &BsplineCurve<Vector3>, plane: &Plane) -> Option<(Curve, f64)> {
    let knots = uv.knot_vector();
    if !uv.is_clamped() || uv.range_tuple() != (*knots.first()?, *knots.last()?) {
        return None;
    }
    let [origin, one, another] = plane.defining_points();
    let axes = (0..3)
        .map(|axis| {
            Some((
                Interval::exact(origin[axis])?,
                Interval::exact(one[axis])?.sub(Interval::exact(origin[axis])?)?,
                Interval::exact(another[axis])?.sub(Interval::exact(origin[axis])?)?,
            ))
        })
        .collect::<Option<Vec<_>>>()?;
    let controls = uv
        .control_points()
        .iter()
        .map(|p| {
            let u = Interval::exact(p.x)?;
            let v = Interval::exact(p.y)?;
            let w = Interval::exact(p.z)?;
            let coordinates = axes
                .iter()
                .map(|(origin, x, y)| origin.mul(w)?.add(x.mul(u)?)?.add(y.mul(v)?))
                .collect::<Option<Vec<_>>>()?;
            Some([coordinates[0], coordinates[1], coordinates[2], w])
        })
        .collect::<Option<Vec<_>>>()?;
    let rounded = choose_controls(&controls)?;
    let error = bernstein::rounding_error(&controls, &rounded)?;
    let curve = native(rounded, uv.knot_vector().to_vec())?;
    Some((curve, error))
}

struct Patch {
    u: [f64; 2],
    v: [f64; 2],
    controls: Vec<Vec<[Interval; 4]>>,
}

fn patches(surface: &BsplineSurface<Vector4>) -> Option<Vec<Patch>> {
    let controls = surface.control_points();
    let columns = controls.first()?.len();
    if controls.len().checked_mul(columns)? > MAX_CONTROLS
        || controls.iter().any(|row| row.len() != columns)
        || controls
            .iter()
            .flatten()
            .any(|p| !finite4(*p) || p.w <= 0.0)
    {
        return None;
    }
    let extraction_score = controls
        .len()
        .checked_mul(columns)?
        .checked_mul(4)?
        .checked_mul(
            surface
                .udegree()
                .checked_mul(surface.knot_vector_u().len())?
                .checked_add(
                    surface
                        .vdegree()
                        .checked_mul(surface.knot_vector_v().len())?,
                )?,
        )?;
    if extraction_score > MAX_PROOF_WORK {
        return None;
    }
    let u = (0..columns)
        .map(|column| {
            let points = controls
                .iter()
                .map(|row| exact(vector(row[column])))
                .collect::<Option<Vec<_>>>()?;
            bernstein::extract(surface.knot_vector_u(), points, surface.udegree())
        })
        .collect::<Option<Vec<_>>>()?;
    let spans = u.first()?.len();
    let mut result = Vec::new();
    for i in 0..spans {
        let rows = (0..=surface.udegree())
            .map(|row| {
                let points = u.iter().map(|column| column[i].controls[row]).collect();
                bernstein::extract(surface.knot_vector_v(), points, surface.vdegree())
            })
            .collect::<Option<Vec<_>>>()?;
        for j in 0..rows.first()?.len() {
            result.push(Patch {
                u: [u[0][i].start, u[0][i].end],
                v: [rows[0][j].start, rows[0][j].end],
                controls: rows.iter().map(|row| row[j].controls.clone()).collect(),
            });
            if result.len() > MAX_SPANS {
                return None;
            }
        }
    }
    Some(result)
}

fn cell<'a>(patches: &'a [Patch], curve: &Bezier<3>) -> Option<&'a Patch> {
    if curve.controls.iter().any(|p| p[2].low <= 0.0) {
        return None;
    }
    let bounds = (0..2)
        .map(|axis| {
            let coordinates = curve
                .controls
                .iter()
                .map(|p| p[axis].div(p[2]))
                .collect::<Option<Vec<_>>>()?;
            Some((
                coordinates
                    .iter()
                    .map(|p| p.low)
                    .fold(f64::INFINITY, f64::min),
                coordinates
                    .iter()
                    .map(|p| p.high)
                    .fold(f64::NEG_INFINITY, f64::max),
            ))
        })
        .collect::<Option<Vec<_>>>()?;
    patches.iter().find(|patch| {
        bounds[0].0 >= patch.u[0]
            && bounds[0].1 <= patch.u[1]
            && bounds[1].0 >= patch.v[0]
            && bounds[1].1 <= patch.v[1]
    })
}
fn power(values: &[Interval], exponent: usize) -> Option<Vec<Interval>> {
    (0..exponent).try_fold(vec![Interval::exact(1.0)?], |result, _| {
        bernstein::product(&result, values)
    })
}
fn basis(
    values: &[Interval],
    complement: &[Interval],
    degree: usize,
) -> Option<Vec<Vec<Interval>>> {
    (0..=degree)
        .map(|i| {
            let coefficient = bernstein::choose(degree, i)?;
            bernstein::product(&power(values, i)?, &power(complement, degree - i)?)?
                .into_iter()
                .map(|value| value.mul(coefficient))
                .collect()
        })
        .collect()
}
fn compose(patch: &Patch, curve: &Bezier<3>) -> Option<Bezier<4>> {
    let m = patch.controls.len().checked_sub(1)?;
    let n = patch.controls.first()?.len().checked_sub(1)?;
    let degree = curve
        .controls
        .len()
        .checked_sub(1)?
        .checked_mul(m.checked_add(n)?)?;
    if degree == 0 || degree > MAX_DEGREE {
        return None;
    }
    let coordinates = |axis: usize, range: [f64; 2]| {
        let low = Interval::exact(range[0])?;
        let width = Interval::exact(range[1])?.sub(low)?;
        let values = curve
            .controls
            .iter()
            .map(|p| p[axis].sub(low.mul(p[2])?)?.div(width))
            .collect::<Option<Vec<_>>>()?;
        let complement = values
            .iter()
            .zip(&curve.controls)
            .map(|(value, p)| p[2].sub(*value))
            .collect::<Option<Vec<_>>>()?;
        Some((values, complement))
    };
    let (u, uc) = coordinates(0, patch.u)?;
    let (v, vc) = coordinates(1, patch.v)?;
    let ub = basis(&u, &uc, m)?;
    let vb = basis(&v, &vc, n)?;
    let terms = patch
        .controls
        .iter()
        .enumerate()
        .flat_map(|(i, row)| row.iter().enumerate().map(move |(j, p)| (i, j, p)))
        .map(|(i, j, p)| Some((*p, bernstein::product(&ub[i], &vb[j])?)))
        .collect::<Option<Vec<_>>>()?;
    let controls = (0..=degree)
        .map(|k| {
            let values = (0..4)
                .map(|axis| {
                    terms
                        .iter()
                        .try_fold(Interval::exact(0.0)?, |sum, (point, coefficient)| {
                            sum.add(point[axis].mul(coefficient[k])?)
                        })
                })
                .collect::<Option<Vec<_>>>()?;
            values.try_into().ok()
        })
        .collect::<Option<Vec<_>>>()?;
    Some(Bezier {
        start: curve.start,
        end: curve.end,
        controls,
    })
}

fn split_linear(source: Bezier<3>, patches: &[Patch]) -> Option<Vec<Bezier<3>>> {
    if source.controls.len() != 2 {
        return Some(vec![source]);
    }
    let mut cuts = Vec::new();
    for (axis, ranges) in [
        (0, patches.iter().map(|p| p.u).collect::<Vec<_>>()),
        (1, patches.iter().map(|p| p.v).collect()),
    ] {
        for knot in ranges.into_iter().flatten() {
            let knot = Interval::exact(knot)?;
            let front = source.controls[0][axis].sub(knot.mul(source.controls[0][2])?)?;
            let back = source.controls[1][axis].sub(knot.mul(source.controls[1][2])?)?;
            if (front.high < 0.0 && back.low > 0.0) || (front.low > 0.0 && back.high < 0.0) {
                let fraction = front.div(front.sub(back)?)?;
                let parameter = Interval::exact(source.start)?
                    .add(
                        Interval::exact(source.end)?
                            .sub(Interval::exact(source.start)?)?
                            .mul(fraction)?,
                    )?
                    .point()?;
                if parameter > source.start && parameter < source.end {
                    cuts.push(parameter);
                }
            }
        }
    }
    cuts.sort_by(f64::total_cmp);
    cuts.dedup();
    if cuts.len() >= MAX_SPANS {
        return None;
    }
    let mut tail = source;
    let mut result = Vec::new();
    for cut in cuts {
        let (left, right) = bernstein::split(&tail, cut)?;
        result.push(left);
        tail = right;
    }
    result.push(tail);
    Some(result)
}

fn assemble(pieces: Vec<Bezier<4>>) -> Option<(Curve, f64)> {
    let degree = pieces.first()?.controls.len().checked_sub(1)?;
    let first = pieces.first()?.start;
    let mut controls: Vec<[f64; 4]> = Vec::new();
    let mut knots = vec![first; degree + 1];
    let mut error = 0.0_f64;
    let mut previous_end = first;
    for piece in pieces {
        if piece.controls.len() != degree + 1
            || piece.start != previous_end
            || piece.end <= piece.start
        {
            return None;
        }
        let original = if let Some(previous) = controls.last() {
            let weight = piece.controls[0][3].midpoint()?;
            let scale = previous[3] / weight;
            if !scale.is_finite() || scale <= 0.0 {
                return None;
            }
            let scale = Interval::exact(scale)?;
            piece
                .controls
                .iter()
                .map(|p| {
                    p.iter()
                        .copied()
                        .map(|value| value.mul(scale))
                        .collect::<Option<Vec<_>>>()?
                        .try_into()
                        .ok()
                })
                .collect::<Option<Vec<[Interval; 4]>>>()?
        } else {
            piece.controls.clone()
        };
        let mut rounded = choose_controls(&original)?;
        if let Some(previous) = controls.last() {
            rounded[0] = *previous;
            knots.extend(iter::repeat_n(piece.start, degree));
        }
        error = error.max(bernstein::rounding_error(&original, &rounded)?);
        if controls.is_empty() {
            controls.extend(rounded);
        } else {
            controls.extend(rounded.into_iter().skip(1));
        }
        if controls.len() > MAX_CONTROLS {
            return None;
        }
        previous_end = piece.end;
    }
    knots.extend(iter::repeat_n(previous_end, degree + 1));
    Some((native(controls, knots)?, error))
}

pub(super) fn lift(curve: &Curve) -> Option<(Curve, f64)> {
    let Curve::ParameterCurve(parameter) = curve else {
        return None;
    };
    let uv = homogeneous(parameter.curve())?;
    if let Surface::Plane(plane) = parameter.surface().as_ref() {
        return planar(&uv, plane);
    }
    let surface = match parameter.surface().as_ref() {
        Surface::BsplineSurface(surface) => surface.try_into_homogeneous_bspline_surface()?,
        Surface::NurbsSurface(surface) => surface.non_rationalized().clone(),
        _ => return None,
    };
    let patches = patches(&surface)?;
    let controls = uv
        .control_points()
        .iter()
        .map(|p| exact([p.x, p.y, p.z]))
        .collect::<Option<Vec<_>>>()?;
    let spans = bernstein::extract(uv.knot_vector(), controls, uv.degree())?;
    let pieces = spans
        .into_iter()
        .map(|span| split_linear(span, &patches))
        .collect::<Option<Vec<_>>>()?
        .into_iter()
        .flatten()
        .collect::<Vec<_>>();
    if pieces.len() > MAX_SPANS {
        return None;
    }
    let mut work = 0usize;
    let lifted = pieces
        .iter()
        .map(|piece| {
            let patch = cell(&patches, piece)?;
            let degree = piece
                .controls
                .len()
                .checked_sub(1)?
                .checked_mul(patch.controls.len() + patch.controls[0].len() - 2)?;
            let score = (degree + 1)
                .checked_pow(3)?
                .checked_mul(patch.controls.len())?
                .checked_mul(patch.controls[0].len())?;
            work = work.checked_add(score)?;
            if work > MAX_PROOF_WORK {
                return None;
            }
            compose(patch, piece)
        })
        .collect::<Option<Vec<_>>>()?;
    assemble(lifted)
}

fn rational(curve: &Curve) -> Option<(BsplineCurve<Vector4>, f64)> {
    match curve {
        Curve::Line(_) | Curve::BsplineCurve(_) | Curve::NurbsCurve(_) => {
            Some((curve.try_into_homogeneous_bspline_curve()?, 0.0))
        }
        Curve::ParameterCurve(_) => {
            let (lifted, error) = lift(curve)?;
            Some((lifted.try_into_homogeneous_bspline_curve()?, error))
        }
        _ => None,
    }
}
fn normalized(curve: &BsplineCurve<Vector4>) -> Option<Vec<Bezier<4>>> {
    if curve.degree() > MAX_DEGREE
        || curve
            .control_points()
            .iter()
            .any(|p| !finite4(*p) || p.w <= 0.0)
    {
        return None;
    }
    let controls = curve
        .control_points()
        .iter()
        .map(|p| exact(vector(*p)))
        .collect::<Option<Vec<_>>>()?;
    let spans = bernstein::extract(curve.knot_vector(), controls, curve.degree())?;
    let start = Interval::exact(spans.first()?.start)?;
    let width = Interval::exact(spans.last()?.end)?.sub(start)?;
    spans
        .into_iter()
        .map(|piece| {
            Some(Bezier {
                start: Interval::exact(piece.start)?
                    .sub(start)?
                    .div(width)?
                    .point()?,
                end: Interval::exact(piece.end)?
                    .sub(start)?
                    .div(width)?
                    .point()?,
                controls: piece.controls,
            })
        })
        .collect()
}
fn restrict(source: &Bezier<4>, start: f64, end: f64) -> Option<Bezier<4>> {
    let tail = if start > source.start {
        bernstein::split(source, start)?.1
    } else {
        source.clone()
    };
    if end < tail.end {
        Some(bernstein::split(&tail, end)?.0)
    } else {
        Some(tail)
    }
}
fn span_distance(a: &[[Interval; 4]], b: &[[Interval; 4]]) -> Option<f64> {
    let aw = a.iter().map(|p| p[3]).collect::<Vec<_>>();
    let bw = b.iter().map(|p| p[3]).collect::<Vec<_>>();
    let denominator = Interval::exact(aw.iter().map(|w| w.low).fold(f64::INFINITY, f64::min))?
        .mul(Interval::exact(
            bw.iter().map(|w| w.low).fold(f64::INFINITY, f64::min),
        )?)?;
    if denominator.low <= 0.0 {
        return None;
    }
    let bounds = (0..3)
        .map(|axis| {
            let lhs = bernstein::product(&a.iter().map(|p| p[axis]).collect::<Vec<_>>(), &bw)?;
            let rhs = bernstein::product(&b.iter().map(|p| p[axis]).collect::<Vec<_>>(), &aw)?;
            let bound = lhs
                .into_iter()
                .zip(rhs)
                .map(|(x, y)| x.sub(y).map(Interval::magnitude))
                .collect::<Option<Vec<_>>>()?
                .into_iter()
                .fold(0.0, f64::max);
            Interval::exact(bound)?.div(denominator)
        })
        .collect::<Option<Vec<_>>>()?;
    bernstein::norm(bounds.try_into().ok()?)
}
pub(super) fn distance(a: &Curve, b: &Curve) -> Option<f64> {
    let (a, ae) = rational(a)?;
    let (b, be) = rational(b)?;
    let a = normalized(&a)?;
    let b = normalized(&b)?;
    let mut knots = a
        .iter()
        .chain(&b)
        .flat_map(|p| [p.start, p.end])
        .collect::<Vec<_>>();
    knots.sort_by(f64::total_cmp);
    knots.dedup();
    if knots.len() > MAX_SPANS + 1 {
        return None;
    }
    let mut work = 0usize;
    let bound = knots
        .windows(2)
        .map(|range| {
            let left = a
                .iter()
                .find(|p| p.start <= range[0] && p.end >= range[1])?;
            let right = b
                .iter()
                .find(|p| p.start <= range[0] && p.end >= range[1])?;
            let score = (left.controls.len() + right.controls.len())
                .checked_pow(2)?
                .checked_mul(left.controls.len())?
                .checked_mul(right.controls.len())?;
            work = work.checked_add(score)?;
            if work > MAX_PROOF_WORK {
                return None;
            }
            span_distance(
                &restrict(left, range[0], range[1])?.controls,
                &restrict(right, range[0], range[1])?.controls,
            )
        })
        .collect::<Option<Vec<_>>>()?
        .into_iter()
        .fold(0.0, f64::max);
    Interval::exact(bound)?
        .add(Interval::exact(ae)?)?
        .add(Interval::exact(be)?)
        .map(|value| value.high)
}
