//! Whole-span admission of a discovered implicit surface-intersection branch.
//!
//! Parameterized Krawczyk contractions establish a unique native root in each
//! UV tube for every leader parameter in its span. One-sided roots at leader
//! corners are enclosed separately and their gap is bounded in model units.
//! This does not establish discovery completeness or global self-separation.

mod interval;
mod surface;

use interval::{Interval, Matrix, Parameters, Vector, distance_bound, dot};
use monstertruck_geometry::prelude::*;
use std::result::Result;
use surface::Carrier;
use thiserror::Error;

const MAX_WORK: usize = 64_000_000;
const MAX_NODES: usize = 8_192;
const MAX_DEPTH: usize = 24;

/// Admission evidence for one discovered intersection branch.
///
/// This transient result applies to the supplied geometry, not subsequently
/// modified carriers, traces or leaders. It is not a serialized authority token.
#[derive(Clone, Debug)]
pub struct ContinuationCertificate {
    /// Version of the bounded continuation admission algorithm.
    pub version: u32,
    /// Number of accepted complete leader spans.
    pub spans: usize,
    /// Number of joins with bounded one-sided endpoint roots.
    pub joins: usize,
    /// Maximum certified gap across those joins, in model length units.
    pub maximum_join_gap: f64,
    /// Requested geometric precision, independent of search or display tolerance.
    pub precision: f64,
    /// Number of retaining halfspaces admitted over every complete span.
    pub halfspaces: usize,
    /// Permitted UV tubes that the native evaluator must enforce for this result.
    pub parameter_tubes: Vec<SurfaceParameterTube>,
}

/// Reasons that whole-span continuation admission could not be established.
///
/// An unresolved proof does not establish that the branch is absent or invalid.
#[derive(Clone, Copy, Debug, Error)]
pub enum ContinuationError {
    /// The carrier has no supported exact parameter map.
    #[error(
        "continuation proof requires a native affine plane or clamped polynomial/rational spline with positive weights and degree 1–16"
    )]
    UnsupportedCarrier,
    /// The guide is not an exactly convertible continuous degree-one spline.
    #[error(
        "continuation proof requires a continuous degree-one guide in the retained traces' parameter frame"
    )]
    UnsupportedGuide,
    /// Geometry, traces or precision have malformed or non-finite source data.
    #[error("continuation proof received invalid geometry, traces or precision")]
    InvalidData,
    /// A bounded enclosure or contraction was inconclusive.
    #[error("continuation proof is unresolved during {stage}; this is not proof of nonexistence")]
    Unresolved {
        /// Inconclusive proof stage.
        stage: &'static str,
    },
    /// A deterministic work or subdivision limit was reached.
    #[error("continuation proof exhausted its bounded work during {stage}")]
    ResourceLimit {
        /// Exhausted resource stage.
        stage: &'static str,
    },
}

#[derive(Default)]
struct Work(usize);
impl Work {
    fn spend(&mut self, count: usize) -> Result<(), ContinuationError> {
        self.0 = self
            .0
            .checked_add(count)
            .ok_or(ContinuationError::ResourceLimit {
                stage: "interval carrier work",
            })?;
        if self.0 > MAX_WORK {
            Err(ContinuationError::ResourceLimit {
                stage: "interval carrier work",
            })
        } else {
            Ok(())
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct GuideSpan {
    range: Interval,
    front: Vector,
    back: Vector,
    normal: Vector,
}

impl GuideSpan {
    fn point(self, parameter: Interval) -> Vector {
        let factor = parameter
            .sub(Interval::point(self.range.low))
            .div(Interval::point(self.range.high).sub(Interval::point(self.range.low)));
        std::array::from_fn(|index| {
            self.front[index].add(self.back[index].sub(self.front[index]).mul(factor))
        })
    }
}

#[derive(Clone, Copy, Debug)]
struct CertifiedSpan {
    parameter: Interval,
    parameters: Parameters,
    guide: GuideSpan,
}

fn jacobian(a: surface::Bounds, b: surface::Bounds, normal: &Vector) -> Matrix {
    let columns = [a.u, a.v, b.u.map(Interval::neg), b.v.map(Interval::neg)];
    let final_row = [
        dot(&a.u, normal),
        dot(&a.v, normal),
        dot(&b.u, normal),
        dot(&b.v, normal),
    ]
    .map(|value| value.mul(Interval::point(0.5)));
    std::array::from_fn(|row| {
        std::array::from_fn(|column| {
            if row < 3 {
                columns[column][row]
            } else {
                final_row[column]
            }
        })
    })
}

fn krawczyk(
    carriers: &[Carrier; 2],
    guide: GuideSpan,
    parameter: Interval,
    parameters: Parameters,
    work: &mut Work,
) -> Result<(Parameters, Matrix, [f64; 4]), ContinuationError> {
    work.spend(128)?;
    let center = parameters.map(Interval::middle);
    let centers = center.map(Interval::point);
    let a = carriers[0].bounds([parameters[0], parameters[1]], work)?;
    let b = carriers[1].bounds([parameters[2], parameters[3]], work)?;
    let ac = carriers[0].bounds([centers[0], centers[1]], work)?;
    let bc = carriers[1].bounds([centers[2], centers[3]], work)?;
    let j = jacobian(a, b, &guide.normal);
    let jc = jacobian(ac, bc, &guide.normal);
    let approximate = Matrix4::from_cols(
        Vector4::new(
            jc[0][0].middle(),
            jc[1][0].middle(),
            jc[2][0].middle(),
            jc[3][0].middle(),
        ),
        Vector4::new(
            jc[0][1].middle(),
            jc[1][1].middle(),
            jc[2][1].middle(),
            jc[3][1].middle(),
        ),
        Vector4::new(
            jc[0][2].middle(),
            jc[1][2].middle(),
            jc[2][2].middle(),
            jc[3][2].middle(),
        ),
        Vector4::new(
            jc[0][3].middle(),
            jc[1][3].middle(),
            jc[2][3].middle(),
            jc[3][3].middle(),
        ),
    );
    let inverse = approximate.invert().ok_or(ContinuationError::Unresolved {
        stage: "transversal Jacobian preconditioner",
    })?;
    let c: Matrix = std::array::from_fn(|row| {
        std::array::from_fn(|column| Interval::point(inverse[column][row]))
    });
    if !c.iter().flatten().all(|value| value.valid()) {
        return Err(ContinuationError::Unresolved {
            stage: "finite Jacobian preconditioner",
        });
    }
    let position = guide.point(parameter);
    let midpoint = std::array::from_fn(|index| {
        ac.point[index]
            .add(bc.point[index])
            .mul(Interval::point(0.5))
            .sub(position[index])
    });
    let residual: Parameters = std::array::from_fn(|row| {
        if row < 3 {
            ac.point[row].sub(bc.point[row])
        } else {
            dot(&midpoint, &guide.normal)
        }
    });
    let contraction: Matrix = std::array::from_fn(|row| {
        std::array::from_fn(|column| {
            let product = (0..4).fold(Interval::ZERO, |sum, index| {
                sum.add(c[row][index].mul(j[index][column]))
            });
            Interval::point(if row == column { 1.0 } else { 0.0 }).sub(product)
        })
    });
    let image: Parameters = std::array::from_fn(|row| {
        let correction = (0..4).fold(Interval::ZERO, |sum, column| {
            sum.add(c[row][column].mul(residual[column]))
        });
        let remainder = (0..4).fold(Interval::ZERO, |sum, column| {
            sum.add(contraction[row][column].mul(parameters[column].sub(centers[column])))
        });
        centers[row].sub(correction).add(remainder)
    });
    if !image.iter().all(|value| value.valid()) {
        return Err(ContinuationError::Unresolved {
            stage: "finite parameterized Krawczyk image",
        });
    }
    let radii = std::array::from_fn(|index| parameters[index].sub(centers[index]).absolute());
    Ok((image, contraction, radii))
}

fn certify_box(
    carriers: &[Carrier; 2],
    guide: GuideSpan,
    parameter: Interval,
    parameters: Parameters,
    work: &mut Work,
) -> Result<bool, ContinuationError> {
    let (image, contraction, radii) = krawczyk(carriers, guide, parameter, parameters, work)?;
    let contained = parameters
        .iter()
        .zip(image)
        .all(|(container, image)| container.interior_contains(image));
    let contracts = (0..4).all(|row| {
        let norm = (0..4).fold(Interval::ZERO, |sum, column| {
            sum.add(
                Interval::point(contraction[row][column].absolute())
                    .mul(Interval::point(radii[column]))
                    .div(Interval::point(radii[row])),
            )
        });
        norm.valid() && norm.high < 1.0
    });
    Ok(contained && contracts)
}

fn seed_box(
    carriers: &[Carrier; 2],
    guide: GuideSpan,
    first: [f64; 4],
    second: [f64; 4],
    precision: f64,
    work: &mut Work,
) -> Result<Parameters, ContinuationError> {
    let center = std::array::from_fn::<_, 4, _>(|index| first[index] * 0.5 + second[index] * 0.5);
    let a = carriers[0].bounds(
        [Interval::point(center[0]), Interval::point(center[1])],
        work,
    )?;
    let b = carriers[1].bounds(
        [Interval::point(center[2]), Interval::point(center[3])],
        work,
    )?;
    let j = jacobian(a, b, &guide.normal);
    let matrix = Matrix4::from_cols(
        Vector4::new(
            j[0][0].middle(),
            j[1][0].middle(),
            j[2][0].middle(),
            j[3][0].middle(),
        ),
        Vector4::new(
            j[0][1].middle(),
            j[1][1].middle(),
            j[2][1].middle(),
            j[3][1].middle(),
        ),
        Vector4::new(
            j[0][2].middle(),
            j[1][2].middle(),
            j[2][2].middle(),
            j[3][2].middle(),
        ),
        Vector4::new(
            j[0][3].middle(),
            j[1][3].middle(),
            j[2][3].middle(),
            j[3][3].middle(),
        ),
    );
    let inverse = matrix.invert().ok_or(ContinuationError::Unresolved {
        stage: "seed Jacobian conditioning",
    })?;
    let normal_size = guide
        .normal
        .iter()
        .fold(0.0, |sum, value| (sum + value.absolute()).next_up());
    let parameters = std::array::from_fn(|row| {
        let sensitivity = inverse[0][row].abs()
            + inverse[1][row].abs()
            + inverse[2][row].abs()
            + inverse[3][row].abs() * normal_size;
        let margin = (precision * sensitivity * 32.0)
            .max(center[row].abs().max(1.0) * (512.0 * f64::EPSILON));
        let radius = ((first[row] - second[row]).abs() * 0.75).max(margin);
        Interval::new(
            (center[row] - radius).next_down(),
            (center[row] + radius).next_up(),
        )
    });
    if parameters
        .iter()
        .all(|parameter| parameter.valid() && parameter.low < parameter.high)
    {
        Ok(parameters)
    } else {
        Err(ContinuationError::Unresolved {
            stage: "finite precision-scaled seed tube",
        })
    }
}

fn endpoint(
    carriers: &[Carrier; 2],
    span: CertifiedSpan,
    t: f64,
    work: &mut Work,
) -> Result<Vector, ContinuationError> {
    let mut parameters = span.parameters;
    for _ in 0..16 {
        let (image, _, _) = krawczyk(carriers, span.guide, Interval::point(t), parameters, work)?;
        parameters = std::array::from_fn(|index| parameters[index].intersect(image[index]));
        if !parameters.iter().all(|value| value.valid()) {
            return Err(ContinuationError::Unresolved {
                stage: "shared endpoint root enclosure",
            });
        }
    }
    Ok(carriers[0]
        .bounds([parameters[0], parameters[1]], work)?
        .point)
}

#[expect(
    clippy::too_many_arguments,
    reason = "whole-span retaining-plane admission inputs"
)]
fn halfspace_admission(
    carriers: &[Carrier; 2],
    guide: GuideSpan,
    parameter: Interval,
    mut parameters: Parameters,
    origin: Point3,
    halfspaces: &[(Point3, Vector3)],
    precision: f64,
    work: &mut Work,
) -> Result<bool, ContinuationError> {
    if halfspaces.is_empty() {
        return Ok(true);
    }
    // The contraction already proves a root in the original tube for every t.
    // Each image intersection further encloses all of those roots; it is not a
    // sampled evaluation, nor another assertion of contraction on a zero-width box.
    for _ in 0..8 {
        let (image, _, _) = krawczyk(carriers, guide, parameter, parameters, work)?;
        parameters = std::array::from_fn(|index| parameters[index].intersect(image[index]));
        if !parameters.iter().all(|value| value.valid()) {
            return Err(ContinuationError::Unresolved {
                stage: "retained halfspace root enclosure",
            });
        }
    }
    for (base, normal) in halfspaces {
        let normal = [normal.x, normal.y, normal.z].map(Interval::point);
        let norm_square = dot(&normal, &normal);
        if !norm_square.valid() || norm_square.low <= 0.0 {
            return Err(ContinuationError::InvalidData);
        }
        let magnitude = Interval::new(
            norm_square.low.sqrt().next_down(),
            norm_square.high.sqrt().next_up(),
        );
        let offset = std::array::from_fn(|index| {
            Interval::point(origin[index]).sub(Interval::point(base[index]))
        });
        let first =
            carriers[0].scalar_bound([parameters[0], parameters[1]], &normal, &offset, work)?;
        let second =
            carriers[1].scalar_bound([parameters[2], parameters[3]], &normal, &offset, work)?;
        // Both carriers represent the same root, so either enclosure can sharpen
        // the bound. A disjoint pair is unresolved rather than evidence of absence.
        let signed = first.intersect(second);
        let threshold = Interval::point(-precision.min(TOLERANCE)).mul(magnitude);
        if !signed.valid() || !threshold.valid() || signed.low < threshold.high {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Certifies bounded continuation of one discovered implicit intersection branch.
///
/// Supported carriers are native affine planes and clamped polynomial or
/// positive-weight rational spline maps of degree 1–16. Angular analytic maps,
/// discontinuous carriers and guides other than continuous degree-one splines
/// are explicitly unsupported. Native end-span spline extrapolation is included.
///
/// `precision` is a geometric residual in model length units; do not supply the
/// tessellation or intersection-search tolerance. Entire-span contractions prove
/// existence and local uniqueness in their UV tubes. Join gaps are bounded to
/// this precision; exact global continuity, injectivity, self-separation and
/// absence of undiscovered branches are not claimed.
///
/// `halfspaces` retain roots with normalized plane residual at least
/// `-min(precision, TOLERANCE)`. Their scalar residuals are enclosed over every
/// accepted complete span, rather than checked only at endpoints or samples.
/// At most 64 finite retaining planes are supported. A contact or clipping
/// boundary that cannot be enclosed within the bounded subdivision budget is
/// unresolved, not evidence of a missing branch.
///
/// Install the returned parameter tubes on the same native [`SurfaceCurve`]
/// before publishing it. Native projection rejects a result outside these tubes
/// and checks its retaining halfspaces. Newton iterations remain unrestricted;
/// the certificate does not guarantee that numerical projection will converge.
/// These transient tubes must be reproved after source changes or deserialization.
pub fn certify_intersection_continuation_with_halfspaces<C, S0, S1>(
    leader: &C,
    surface0: &S0,
    surface1: &S1,
    traces: &[BsplineCurve<Point2>; 2],
    halfspaces: &[(Point3, Vector3)],
    precision: f64,
) -> Result<ContinuationCertificate, ContinuationError>
where
    C: ParametricCurve3D + TryIntoHomogeneousBsplineCurve + Clone,
    S0: ParametricSurface3D
        + SearchNearestParameter<SurfaceParameter, Point = Point3>
        + SupportsExactPatchDomains
        + TryIntoHomogeneousBsplineSurface
        + TryIntoAnalyticSurfaceKind
        + Clone,
    S1: ParametricSurface3D
        + SearchNearestParameter<SurfaceParameter, Point = Point3>
        + SupportsExactPatchDomains
        + TryIntoHomogeneousBsplineSurface
        + TryIntoAnalyticSurfaceKind
        + Clone,
{
    if !precision.is_finite()
        || precision <= 0.0
        || halfspaces.len() > 64
        || halfspaces.iter().any(|(base, normal)| {
            !base.x.is_finite()
                || !base.y.is_finite()
                || !base.z.is_finite()
                || !normal.x.is_finite()
                || !normal.y.is_finite()
                || !normal.z.is_finite()
                || normal.magnitude2() == 0.0
                || !normal.magnitude2().is_finite()
        })
    {
        return Err(ContinuationError::InvalidData);
    }
    let source = leader
        .try_into_homogeneous_bspline_curve()
        .ok_or(ContinuationError::UnsupportedGuide)?;
    let controls = source.control_points();
    let knots = source.knot_vector();
    if !(2..=32_768).contains(&controls.len())
        || knots.len() != controls.len() + 2
        || controls
            .iter()
            .any(|h| h.w != 1.0 || !h.x.is_finite() || !h.y.is_finite() || !h.z.is_finite())
        || knots.iter().any(|value| !value.is_finite())
        || knots.windows(2).any(|pair| pair[0] > pair[1])
        || knots[0] != knots[1]
        || knots[knots.len() - 1] != knots[knots.len() - 2]
        || knots[1..knots.len() - 1]
            .windows(2)
            .any(|pair| pair[0] >= pair[1])
    {
        return Err(ContinuationError::UnsupportedGuide);
    }
    let curve = SurfaceCurve::with_boundaries(
        surface0.clone(),
        surface1.clone(),
        leader.clone(),
        Option::<()>::None,
        Option::<()>::None,
    )
    .with_parameter_traces(traces.clone())
    .ok_or(ContinuationError::InvalidData)?;
    if traces[0].range_tuple() != source.range_tuple()
        || traces[1].range_tuple() != source.range_tuple()
    {
        return Err(ContinuationError::InvalidData);
    }
    let origin = Point3::new(controls[0].x, controls[0].y, controls[0].z);
    let carriers = [
        Carrier::new(surface0, origin)?,
        Carrier::new(surface1, origin)?,
    ];
    let local = |point: &Vector4| {
        std::array::from_fn(|index| {
            Interval::point(point[index]).sub(Interval::point(origin[index]))
        })
    };
    let mut work = Work::default();
    let mut accepted = Vec::new();
    let mut nodes = 0;
    for (index, control_pair) in controls.windows(2).enumerate() {
        let guide_range = Interval::new(knots[index + 1], knots[index + 2]);
        let middle = guide_range.middle();
        let normal = leader.derivative(middle);
        if !normal.x.is_finite()
            || !normal.y.is_finite()
            || !normal.z.is_finite()
            || normal.magnitude2() == 0.0
        {
            return Err(ContinuationError::Unresolved {
                stage: "finite leader cross-section",
            });
        }
        let guide = GuideSpan {
            range: guide_range,
            front: local(&control_pair[0]),
            back: local(&control_pair[1]),
            normal: [normal.x, normal.y, normal.z].map(Interval::point),
        };
        let mut pending = vec![(guide_range, 0)];
        while let Some((parameter, depth)) = pending.pop() {
            nodes += 1;
            if nodes > MAX_NODES {
                return Err(ContinuationError::ResourceLimit {
                    stage: "continuation tube subdivisions",
                });
            }
            let qa =
                curve
                    .search_triple(parameter.low, 100)
                    .ok_or(ContinuationError::Unresolved {
                        stage: "native continuation seed",
                    })?;
            let qb =
                curve
                    .search_triple(parameter.high, 100)
                    .ok_or(ContinuationError::Unresolved {
                        stage: "native continuation seed",
                    })?;
            let values_a = [qa.1.x, qa.1.y, qa.2.x, qa.2.y];
            let values_b = [qb.1.x, qb.1.y, qb.2.x, qb.2.y];
            let result = seed_box(&carriers, guide, values_a, values_b, precision, &mut work)
                .and_then(|parameters| {
                    certify_box(&carriers, guide, parameter, parameters, &mut work).and_then(
                        |certified| {
                            if certified {
                                halfspace_admission(
                                    &carriers, guide, parameter, parameters, origin, halfspaces,
                                    precision, &mut work,
                                )
                                .map(|retained| (retained, parameters))
                            } else {
                                Ok((false, parameters))
                            }
                        },
                    )
                });
            match result {
                Ok((true, parameters)) => accepted.push(CertifiedSpan {
                    parameter,
                    parameters,
                    guide,
                }),
                Err(error @ ContinuationError::ResourceLimit { .. }) => return Err(error),
                _ => {
                    let middle = parameter.middle();
                    if depth == MAX_DEPTH || middle <= parameter.low || middle >= parameter.high {
                        return Err(ContinuationError::Unresolved {
                            stage: "whole-span contraction and retained halfspace admission",
                        });
                    }
                    pending.extend([
                        (Interval::new(middle, parameter.high), depth + 1),
                        (Interval::new(parameter.low, middle), depth + 1),
                    ]);
                }
            }
        }
    }
    let mut maximum_join_gap = 0.0_f64;
    for pair in accepted.windows(2) {
        if pair[0].parameter.high != pair[1].parameter.low {
            return Err(ContinuationError::InvalidData);
        }
        let left = endpoint(&carriers, pair[0], pair[0].parameter.high, &mut work)?;
        let right = endpoint(&carriers, pair[1], pair[1].parameter.low, &mut work)?;
        let gap = distance_bound(&left, &right);
        if !gap.is_finite() || gap > precision {
            return Err(ContinuationError::Unresolved {
                stage: "one-sided endpoint join precision",
            });
        }
        maximum_join_gap = maximum_join_gap.max(gap);
    }
    let first = *accepted.first().ok_or(ContinuationError::InvalidData)?;
    let last = *accepted.last().ok_or(ContinuationError::InvalidData)?;
    let front = Point3::new(controls[0].x, controls[0].y, controls[0].z);
    let back_control = controls.last().ok_or(ContinuationError::InvalidData)?;
    let back = Point3::new(back_control.x, back_control.y, back_control.z);
    let mut joins = accepted.len().saturating_sub(1);
    if front.distance(back) <= precision {
        let a = endpoint(&carriers, first, first.parameter.low, &mut work)?;
        let b = endpoint(&carriers, last, last.parameter.high, &mut work)?;
        let gap = distance_bound(&a, &b);
        if !gap.is_finite() || gap > precision {
            return Err(ContinuationError::Unresolved {
                stage: "closed-branch endpoint join precision",
            });
        }
        maximum_join_gap = maximum_join_gap.max(gap);
        joins += 1;
    }
    let parameter_tubes = accepted
        .iter()
        .map(|span| SurfaceParameterTube {
            parameter_range: (span.parameter.low, span.parameter.high),
            uv_bounds: span.parameters.map(|value| (value.low, value.high)),
        })
        .collect();
    Ok(ContinuationCertificate {
        version: 1,
        spans: accepted.len(),
        joins,
        maximum_join_gap,
        precision,
        halfspaces: halfspaces.len(),
        parameter_tubes,
    })
}

/// Certifies a discovered branch without additional retaining halfspaces.
///
/// See [`certify_intersection_continuation_with_halfspaces`] for the supported
/// carriers, geometric precision policy, runtime guard contract and limitations.
pub fn certify_intersection_continuation<C, S0, S1>(
    leader: &C,
    surface0: &S0,
    surface1: &S1,
    traces: &[BsplineCurve<Point2>; 2],
    precision: f64,
) -> Result<ContinuationCertificate, ContinuationError>
where
    C: ParametricCurve3D + TryIntoHomogeneousBsplineCurve + Clone,
    S0: ParametricSurface3D
        + SearchNearestParameter<SurfaceParameter, Point = Point3>
        + SupportsExactPatchDomains
        + TryIntoHomogeneousBsplineSurface
        + TryIntoAnalyticSurfaceKind
        + Clone,
    S1: ParametricSurface3D
        + SearchNearestParameter<SurfaceParameter, Point = Point3>
        + SupportsExactPatchDomains
        + TryIntoHomogeneousBsplineSurface
        + TryIntoAnalyticSurfaceKind
        + Clone,
{
    certify_intersection_continuation_with_halfspaces(
        leader,
        surface0,
        surface1,
        traces,
        &[],
        precision,
    )
}
