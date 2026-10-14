//! Native spline values and first derivatives over complete parameter rectangles.
//!
//! Recurrences use outward arithmetic directly on source controls and knots.
//! Rounded knot insertion or sampled surface fits are not proof inputs.

use super::interval::{Interval, Vector, dot};
use super::{ContinuationError, Work};
use monstertruck_geometry::prelude::*;
use std::result::Result;

#[derive(Debug)]
struct Axis {
    knots: Vec<f64>,
    degree: usize,
    controls: usize,
}

impl Axis {
    fn new(knots: &KnotVector, controls: usize) -> Result<Self, ContinuationError> {
        let degree = knots
            .len()
            .checked_sub(controls + 1)
            .ok_or(ContinuationError::InvalidData)?;
        if !(1..=16).contains(&degree)
            || controls < degree + 1
            || knots.iter().any(|value| !value.is_finite())
            || knots.windows(2).any(|pair| pair[0] > pair[1])
            || knots[0] >= knots[knots.len() - 1]
            || knots[..=degree].iter().any(|knot| *knot != knots[0])
            || knots[knots.len() - degree - 1..]
                .iter()
                .any(|knot| *knot != knots[knots.len() - 1])
        {
            return Err(ContinuationError::UnsupportedCarrier);
        }
        let mut last = knots[0];
        let mut count = 0;
        for knot in knots {
            count = if *knot == last { count + 1 } else { 1 };
            last = *knot;
            if *knot > knots[0] && *knot < knots[knots.len() - 1] && count > degree {
                return Err(ContinuationError::UnsupportedCarrier);
            }
        }
        Ok(Self {
            knots: knots.iter().copied().collect(),
            degree,
            controls,
        })
    }

    fn on_cell(
        &self,
        parameter: Interval,
        cell: usize,
        work: &mut Work,
    ) -> Result<(Vec<Interval>, Vec<Interval>), ContinuationError> {
        let mut previous = (0..self.knots.len() - 1)
            .map(|index| {
                if index == cell {
                    Interval::ONE
                } else {
                    Interval::ZERO
                }
            })
            .collect::<Vec<_>>();
        let mut lower = previous.clone();
        for degree in 1..=self.degree {
            work.spend(previous.len())?;
            let next = (0..previous.len() - 1)
                .map(|index| {
                    let first_span = Interval::point(self.knots[index + degree])
                        .sub(Interval::point(self.knots[index]));
                    let second_span = Interval::point(self.knots[index + degree + 1])
                        .sub(Interval::point(self.knots[index + 1]));
                    let first = if self.knots[index + degree] == self.knots[index] {
                        Interval::ZERO
                    } else {
                        parameter
                            .sub(Interval::point(self.knots[index]))
                            .div(first_span)
                            .mul(previous[index])
                    };
                    let second = if self.knots[index + degree + 1] == self.knots[index + 1] {
                        Interval::ZERO
                    } else {
                        Interval::point(self.knots[index + degree + 1])
                            .sub(parameter)
                            .div(second_span)
                            .mul(previous[index + 1])
                    };
                    first.add(second)
                })
                .collect::<Vec<_>>();
            if degree == self.degree {
                lower = previous;
            }
            previous = next;
        }
        let derivatives = (0..self.controls)
            .map(|index| {
                let left = if self.knots[index + self.degree] == self.knots[index] {
                    Interval::ZERO
                } else {
                    Interval::point(self.degree as f64)
                        .div(
                            Interval::point(self.knots[index + self.degree])
                                .sub(Interval::point(self.knots[index])),
                        )
                        .mul(lower[index])
                };
                let right = if self.knots[index + self.degree + 1] == self.knots[index + 1] {
                    Interval::ZERO
                } else {
                    Interval::point(self.degree as f64)
                        .div(
                            Interval::point(self.knots[index + self.degree + 1])
                                .sub(Interval::point(self.knots[index + 1])),
                        )
                        .mul(lower[index + 1])
                };
                left.sub(right)
            })
            .collect::<Vec<_>>();
        Ok((previous, derivatives))
    }

    fn bounds(
        &self,
        parameter: Interval,
        work: &mut Work,
    ) -> Result<(Vec<Interval>, Vec<Interval>), ContinuationError> {
        if !parameter.valid() {
            return Err(ContinuationError::InvalidData);
        }
        let cells = (self.degree..self.controls)
            .filter(|index| self.knots[*index] < self.knots[*index + 1])
            .collect::<Vec<_>>();
        let first = *cells.first().ok_or(ContinuationError::InvalidData)?;
        let last = *cells.last().ok_or(ContinuationError::InvalidData)?;
        let mut result: Option<(Vec<Interval>, Vec<Interval>)> = None;
        for cell in cells {
            let low = if cell == first {
                parameter.low
            } else {
                parameter.low.max(self.knots[cell])
            };
            let high = if cell == last {
                parameter.high
            } else {
                parameter.high.min(self.knots[cell + 1])
            };
            if low > high {
                continue;
            }
            let values = self.on_cell(Interval::new(low, high), cell, work)?;
            result = Some(match result {
                None => values,
                Some((basis, derivatives)) => (
                    basis
                        .into_iter()
                        .zip(values.0)
                        .map(|(a, b)| a.hull(b))
                        .collect(),
                    derivatives
                        .into_iter()
                        .zip(values.1)
                        .map(|(a, b)| a.hull(b))
                        .collect(),
                ),
            });
        }
        result.ok_or(ContinuationError::Unresolved {
            stage: "native knot-domain enclosure",
        })
    }
}

#[derive(Debug)]
pub(super) struct Spline {
    u: Axis,
    v: Axis,
    controls: Vec<Vec<[Interval; 4]>>,
    constant_weight: Option<f64>,
}

#[derive(Debug)]
pub(super) enum Carrier {
    Plane {
        origin: Vector,
        u: Vector,
        v: Vector,
    },
    Spline(Box<Spline>),
}

#[derive(Clone, Copy, Debug)]
pub(super) struct Bounds {
    pub(super) point: Vector,
    pub(super) u: Vector,
    pub(super) v: Vector,
}

impl Carrier {
    pub(super) fn new<S>(surface: &S, origin: Point3) -> Result<Self, ContinuationError>
    where S: ParametricSurface3D
            + SupportsExactPatchDomains
            + TryIntoHomogeneousBsplineSurface
            + TryIntoAnalyticSurfaceKind {
        if !surface.supports_exact_patch_domains() {
            let Some(AnalyticSurfaceKind::Plane(plane)) = surface.try_into_analytic_surface_kind()
            else {
                return Err(ContinuationError::UnsupportedCarrier);
            };
            let [base, one, another] = plane.defining_points();
            return Ok(Self::Plane {
                origin: [base.x, base.y, base.z]
                    .map(Interval::point)
                    .into_iter()
                    .zip([origin.x, origin.y, origin.z])
                    .map(|(a, b)| a.sub(Interval::point(b)))
                    .collect::<Vec<_>>()
                    .try_into()
                    .map_err(|_| ContinuationError::InvalidData)?,
                u: std::array::from_fn(|index| {
                    Interval::point(one[index]).sub(Interval::point(base[index]))
                }),
                v: std::array::from_fn(|index| {
                    Interval::point(another[index]).sub(Interval::point(base[index]))
                }),
            });
        }
        let net = surface
            .try_into_homogeneous_bspline_surface()
            .ok_or(ContinuationError::UnsupportedCarrier)?;
        let rows = net.control_points();
        let columns = rows.first().ok_or(ContinuationError::InvalidData)?.len();
        if columns == 0
            || rows.len() > 65_536
            || columns > 65_536
            || rows
                .len()
                .checked_mul(columns)
                .is_none_or(|size| size > 65_536)
            || rows.iter().any(|row| row.len() != columns)
            || rows.iter().flatten().any(|h| {
                !h.x.is_finite()
                    || !h.y.is_finite()
                    || !h.z.is_finite()
                    || !h.w.is_finite()
                    || h.w <= 0.0
            })
        {
            return Err(ContinuationError::InvalidData);
        }
        let weight = rows[0][0].w;
        let constant_weight = rows
            .iter()
            .flatten()
            .all(|h| h.w == weight)
            .then_some(weight);
        let controls = rows
            .iter()
            .map(|row| {
                row.iter()
                    .map(|h| {
                        [
                            Interval::point(h.x)
                                .sub(Interval::point(origin.x).mul(Interval::point(h.w))),
                            Interval::point(h.y)
                                .sub(Interval::point(origin.y).mul(Interval::point(h.w))),
                            Interval::point(h.z)
                                .sub(Interval::point(origin.z).mul(Interval::point(h.w))),
                            Interval::point(h.w),
                        ]
                    })
                    .collect()
            })
            .collect();
        Ok(Self::Spline(Box::new(Spline {
            u: Axis::new(net.knot_vector_u(), rows.len())?,
            v: Axis::new(net.knot_vector_v(), columns)?,
            controls,
            constant_weight,
        })))
    }

    pub(super) fn bounds(
        &self,
        uv: [Interval; 2],
        work: &mut Work,
    ) -> Result<Bounds, ContinuationError> {
        let bounds = match self {
            Self::Plane { origin, u, v } => Bounds {
                point: std::array::from_fn(|index| {
                    origin[index]
                        .add(u[index].mul(uv[0]))
                        .add(v[index].mul(uv[1]))
                }),
                u: *u,
                v: *v,
            },
            Self::Spline(spline) => {
                let (u, du) = spline.u.bounds(uv[0], work)?;
                let (v, dv) = spline.v.bounds(uv[1], work)?;
                work.spend(
                    spline
                        .controls
                        .len()
                        .saturating_mul(spline.controls[0].len())
                        .saturating_mul(12),
                )?;
                let mut h = [Interval::ZERO; 4];
                let mut hu = [Interval::ZERO; 4];
                let mut hv = [Interval::ZERO; 4];
                for (i, row) in spline.controls.iter().enumerate() {
                    for (j, control) in row.iter().enumerate() {
                        for component in 0..4 {
                            h[component] = h[component].add(control[component].mul(u[i].mul(v[j])));
                            hu[component] =
                                hu[component].add(control[component].mul(du[i].mul(v[j])));
                            hv[component] =
                                hv[component].add(control[component].mul(u[i].mul(dv[j])));
                        }
                    }
                }
                if let Some(weight) = spline.constant_weight {
                    h[3] = Interval::point(weight);
                    hu[3] = Interval::ZERO;
                    hv[3] = Interval::ZERO;
                }
                if !h[3].valid() || h[3].low <= 0.0 {
                    return Err(ContinuationError::Unresolved {
                        stage: "positive rational denominator",
                    });
                }
                let point = std::array::from_fn(|index| h[index].div(h[3]));
                Bounds {
                    point,
                    u: std::array::from_fn(|index| {
                        hu[index].sub(point[index].mul(hu[3])).div(h[3])
                    }),
                    v: std::array::from_fn(|index| {
                        hv[index].sub(point[index].mul(hv[3])).div(h[3])
                    }),
                }
            }
        };
        if bounds
            .point
            .iter()
            .chain(&bounds.u)
            .chain(&bounds.v)
            .all(|value| value.valid())
        {
            Ok(bounds)
        } else {
            Err(ContinuationError::Unresolved {
                stage: "finite native carrier enclosure",
            })
        }
    }

    /// Bounds a plane residual before expanding coordinates independently.
    /// This preserves exact zero coefficients for a native affine carrier
    /// parallel to the retaining plane, and projects homogeneous spline controls
    /// before summing their basis values.
    pub(super) fn scalar_bound(
        &self,
        uv: [Interval; 2],
        normal: &Vector,
        offset: &Vector,
        work: &mut Work,
    ) -> Result<Interval, ContinuationError> {
        let scalar = match self {
            Self::Plane { origin, u, v } => dot(origin, normal)
                .add(dot(u, normal).mul(uv[0]))
                .add(dot(v, normal).mul(uv[1])),
            Self::Spline(spline) => {
                let (u, _) = spline.u.bounds(uv[0], work)?;
                let (v, _) = spline.v.bounds(uv[1], work)?;
                work.spend(
                    spline
                        .controls
                        .len()
                        .saturating_mul(spline.controls[0].len())
                        .saturating_mul(8),
                )?;
                let mut numerator = Interval::ZERO;
                let mut denominator = Interval::ZERO;
                for (i, row) in spline.controls.iter().enumerate() {
                    for (j, control) in row.iter().enumerate() {
                        let basis = u[i].mul(v[j]);
                        let projected = dot(&[control[0], control[1], control[2]], normal);
                        numerator = numerator.add(projected.mul(basis));
                        denominator = denominator.add(control[3].mul(basis));
                    }
                }
                if let Some(weight) = spline.constant_weight {
                    denominator = Interval::point(weight);
                }
                if !denominator.valid() || denominator.low <= 0.0 {
                    return Err(ContinuationError::Unresolved {
                        stage: "retained halfspace rational denominator",
                    });
                }
                numerator.div(denominator)
            }
        }
        .add(dot(offset, normal));
        if scalar.valid() {
            Ok(scalar)
        } else {
            Err(ContinuationError::Unresolved {
                stage: "finite retained halfspace enclosure",
            })
        }
    }
}
