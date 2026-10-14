//! Original-knot extraction and positive rational Bernstein bounds.
use super::{MAX_CONTROLS, MAX_DEGREE, MAX_PROOF_WORK, MAX_SPANS, interval::Interval};

#[derive(Clone)]
pub(super) struct Bezier<const N: usize> {
    pub(super) start: f64,
    pub(super) end: f64,
    pub(super) controls: Vec<[Interval; N]>,
}

fn insert<const N: usize>(
    knots: &mut Vec<f64>,
    controls: &mut Vec<[Interval; N]>,
    degree: usize,
    knot: f64,
) -> Option<()> {
    let k = knots.iter().rposition(|value| *value <= knot)?;
    let multiplicity = knots.iter().filter(|value| **value == knot).count();
    if multiplicity >= degree || k < degree || k < multiplicity {
        return None;
    }
    let next = (0..=controls.len())
        .map(|i| {
            if i <= k - degree {
                return Some(controls[i]);
            }
            if i > k - multiplicity {
                return Some(controls[i - 1]);
            }
            let alpha = Interval::exact(knot)?
                .sub(Interval::exact(knots[i])?)?
                .div(Interval::exact(knots[i + degree])?.sub(Interval::exact(knots[i])?)?)?;
            let complement = Interval::exact(1.0)?.sub(alpha)?;
            let values = (0..N)
                .map(|axis| {
                    controls[i][axis]
                        .mul(alpha)?
                        .add(controls[i - 1][axis].mul(complement)?)
                })
                .collect::<Option<Vec<_>>>()?;
            values.try_into().ok()
        })
        .collect::<Option<Vec<_>>>()?;
    knots.insert(k + 1, knot);
    *controls = next;
    Some(())
}

/// Extracts clamped original B-spline spans with interval knot-insertion error.
pub(super) fn extract<const N: usize>(
    source_knots: &[f64],
    source_controls: Vec<[Interval; N]>,
    degree: usize,
) -> Option<Vec<Bezier<N>>> {
    if degree == 0
        || degree > MAX_DEGREE
        || source_controls.len() > MAX_CONTROLS
        || source_knots.len() != source_controls.len() + degree + 1
        || source_knots.iter().any(|k| !k.is_finite())
        || source_knots.windows(2).any(|k| k[0] > k[1])
    {
        return None;
    }
    let first = *source_knots.first()?;
    let last = *source_knots.last()?;
    if last <= first
        || source_knots.iter().take(degree + 1).any(|k| *k != first)
        || source_knots
            .iter()
            .rev()
            .take(degree + 1)
            .any(|k| *k != last)
    {
        return None;
    }
    let mut unique = source_knots.to_vec();
    unique.dedup();
    if unique.len() > MAX_SPANS + 1 {
        return None;
    }
    if source_controls
        .len()
        .checked_mul(degree)?
        .checked_mul(unique.len())?
        .checked_mul(N)?
        > MAX_PROOF_WORK
    {
        return None;
    }
    let mut knots = source_knots.to_vec();
    let mut controls = source_controls;
    for knot in unique.iter().copied().filter(|k| *k > first && *k < last) {
        let multiplicity = knots.iter().filter(|k| **k == knot).count();
        if multiplicity > degree {
            return None;
        }
        for _ in multiplicity..degree {
            insert(&mut knots, &mut controls, degree, knot)?;
        }
        if controls.len() > MAX_CONTROLS {
            return None;
        }
    }
    unique
        .windows(2)
        .enumerate()
        .map(|(i, span)| {
            let start = i * degree;
            Some(Bezier {
                start: span[0],
                end: span[1],
                controls: controls.get(start..=start + degree)?.to_vec(),
            })
        })
        .collect()
}

pub(super) fn split<const N: usize>(
    source: &Bezier<N>,
    parameter: f64,
) -> Option<(Bezier<N>, Bezier<N>)> {
    if parameter <= source.start || parameter >= source.end {
        return None;
    }
    let alpha = Interval::exact(parameter)?
        .sub(Interval::exact(source.start)?)?
        .div(Interval::exact(source.end)?.sub(Interval::exact(source.start)?)?)?;
    let complement = Interval::exact(1.0)?.sub(alpha)?;
    let mut levels = vec![source.controls.clone()];
    while levels.last()?.len() > 1 {
        let next = levels
            .last()?
            .windows(2)
            .map(|pair| {
                let values = (0..N)
                    .map(|axis| {
                        pair[0][axis]
                            .mul(complement)?
                            .add(pair[1][axis].mul(alpha)?)
                    })
                    .collect::<Option<Vec<_>>>()?;
                values.try_into().ok()
            })
            .collect::<Option<Vec<_>>>()?;
        levels.push(next);
    }
    Some((
        Bezier {
            start: source.start,
            end: parameter,
            controls: levels.iter().map(|row| row[0]).collect(),
        },
        Bezier {
            start: parameter,
            end: source.end,
            controls: levels.iter().rev().map(|row| row[row.len() - 1]).collect(),
        },
    ))
}

pub(super) fn choose(n: usize, k: usize) -> Option<Interval> {
    (0..k.min(n - k)).try_fold(Interval::exact(1.0)?, |value, i| {
        value
            .mul(Interval::exact((n - i) as f64)?)?
            .div(Interval::exact((i + 1) as f64)?)
    })
}

pub(super) fn product(a: &[Interval], b: &[Interval]) -> Option<Vec<Interval>> {
    let m = a.len().checked_sub(1)?;
    let n = b.len().checked_sub(1)?;
    if m + n > MAX_DEGREE * 2 {
        return None;
    }
    (0..=m + n)
        .map(|k| {
            (k.saturating_sub(n)..=k.min(m)).try_fold(Interval::exact(0.0)?, |sum, i| {
                let factor = choose(m, i)?
                    .mul(choose(n, k - i)?)?
                    .div(choose(m + n, k)?)?;
                sum.add(a[i].mul(b[k - i])?.mul(factor)?)
            })
        })
        .collect()
}

pub(super) fn norm(bounds: [Interval; 3]) -> Option<f64> {
    bounds
        .into_iter()
        .try_fold(Interval::exact(0.0)?, |sum, b| {
            let b = Interval::exact(b.magnitude())?;
            sum.add(b.mul(b)?)
        })?
        .sqrt()
        .map(|b| b.high)
}

/// Controls enclose the original polynomial, while rounded controls are the
/// emitted geometry. Includes numerator, positive denominator and output error.
pub(super) fn rounding_error(original: &[[Interval; 4]], rounded: &[[f64; 4]]) -> Option<f64> {
    if original.len() != rounded.len()
        || rounded
            .iter()
            .any(|p| p.iter().any(|x| !x.is_finite()) || p[3] <= 0.0)
    {
        return None;
    }
    let minimum = original
        .iter()
        .map(|p| p[3].low)
        .fold(f64::INFINITY, f64::min);
    if minimum <= 0.0 {
        return None;
    }
    let errors = (0..4)
        .map(|axis| {
            original
                .iter()
                .zip(rounded)
                .map(|(a, b)| {
                    a[axis]
                        .sub(Interval::exact(b[axis])?)
                        .map(Interval::magnitude)
                })
                .collect::<Option<Vec<_>>>()
                .map(|values| values.into_iter().fold(0.0, f64::max))
        })
        .collect::<Option<Vec<_>>>()?;
    let bounds = (0..3)
        .map(|axis| {
            let magnitude = rounded
                .iter()
                .map(|p| {
                    Interval::exact(p[axis])?
                        .div(Interval::exact(p[3])?)
                        .map(Interval::magnitude)
                })
                .collect::<Option<Vec<_>>>()?
                .into_iter()
                .fold(0.0, f64::max);
            Interval::exact(errors[axis])?
                .add(Interval::exact(magnitude)?.mul(Interval::exact(errors[3])?)?)?
                .div(Interval::exact(minimum)?)
        })
        .collect::<Option<Vec<_>>>()?;
    norm(bounds.try_into().ok()?)
}
