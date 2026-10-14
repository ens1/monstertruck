//! Classic (0.3.2) face division.
//!
//! Ported verbatim from the published 0.3.2 crate's `transversal::divide_face`:
//! projects each loop's edges into the face parameter domain, splits the face
//! into positively oriented pre-faces with their contained holes, and tags each
//! with its `and`/`or`/`unknown` status.

#![allow(clippy::many_single_char_names)]

use super::faces_classification::FacesClassification;
use super::loops_store::*;
use monstertruck_meshing::prelude::*;
use monstertruck_topology::*;
use rustc_hash::FxHashMap as HashMap;
use std::collections::hash_map::Entry;
use std::ops::Deref;

const MAX_PARAMETER_POINTS: usize = 1_000_000;
type DomainResult<T> = std::result::Result<T, &'static str>;
type FaceDivisionResult<C, S> = std::result::Result<
    FacesClassification<Point3, super::AltCurve<C, S>, S>,
    (usize, &'static str),
>;

fn project_parameter<S>(
    surface: &S,
    point: Point3,
    hint: Option<Point2>,
    tol: f64,
) -> Option<Point2>
where
    S: ParametricSurface3D
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    let seed = hint.map(Into::into);
    let (mut u, mut v) = surface
        .search_parameter(point, seed, 100)
        .or_else(|| surface.search_nearest_parameter(point, seed, 100))?;
    if let Some(hint) = hint {
        if let Some(period) = surface
            .period_u()
            .filter(|period| period.is_finite() && *period > 0.0)
        {
            u += ((hint.x - u) / period).round() * period;
        }
        if let Some(period) = surface
            .period_v()
            .filter(|period| period.is_finite() && *period > 0.0)
        {
            v += ((hint.y - v) / period).round() * period;
        }
    }
    if !u.is_finite() || !v.is_finite() {
        return None;
    }
    let projected = surface.subs(u, v);
    let error = projected.distance(point);
    (error.is_finite() && error <= tol).then_some(Point2::new(u, v))
}

fn create_parameter_boundary<C: super::ShapeOpsCurve<S>, S: super::ShapeOpsSurface>(
    face: &Face<Point3, super::AltCurve<C, S>, S>,
    wire: &Wire<Point3, super::AltCurve<C, S>>,
    polys: &mut HashMap<EdgeId<super::AltCurve<C, S>>, PolylineCurve<Point3>>,
    tol: f64,
) -> DomainResult<PolylineCurve<Point2>> {
    let surface = face.surface();
    let mut values: Vec<Point2> = Vec::new();
    for edge in wire.edge_iter() {
        let curve = edge.curve();
        let mut parameters = match &curve {
            crate::alternative::Alternative::SecondType(branch) => branch
                .parameter_samples(face.id(), tol)
                .ok_or("native branch continuation in its producing face")?,
            crate::alternative::Alternative::FirstType(native) => {
                if let Some(parameters) = native.parameter_boundary_2d(&surface, tol) {
                    parameters
                } else {
                    if let Entry::Vacant(entry) = polys.entry(edge.id()) {
                        let points = native
                            .try_parameter_division(native.range_tuple(), tol)
                            .ok_or("finite native edge sampling")?
                            .1;
                        entry.insert(PolylineCurve(points));
                    }
                    let mut hint = None;
                    polys
                        .get(&edge.id())
                        .ok_or("native edge sampling cache")?
                        .iter()
                        .map(|point| {
                            let parameter = project_parameter(&surface, *point, hint, tol)?;
                            hint = Some(parameter);
                            Some(parameter)
                        })
                        .collect::<Option<Vec<_>>>()
                        .ok_or("native boundary projection")?
                }
            }
        };
        if !edge.orientation() {
            parameters.reverse();
        }
        if parameters.len() < 2
            || values
                .len()
                .checked_add(parameters.len())
                .ok_or("parameter point budget")?
                > MAX_PARAMETER_POINTS
        {
            return Err("invalid or excessive parameter-boundary samples");
        }
        let first = *parameters.first().ok_or("empty parameter boundary")?;
        let last = *parameters.last().ok_or("empty parameter boundary")?;
        let on_first = surface.subs(first.x, first.y);
        let on_last = surface.subs(last.x, last.y);
        if !on_first.near(&edge.front().point()) || !on_last.near(&edge.back().point()) {
            return Err("parameter-boundary endpoints do not match their shared vertices");
        }
        if parameters
            .iter()
            .any(|point| !point.x.is_finite() || !point.y.is_finite())
        {
            return Err("non-finite face-local parameters");
        }
        if let Some(previous) = values.last() {
            let offset = |value: f64, seed: f64, period: Option<f64>| {
                period
                    .filter(|period| period.is_finite() && *period > 0.0)
                    .map_or(0.0, |period| ((seed - value) / period).round() * period)
            };
            let shift = Vector2::new(
                offset(first.x, previous.x, surface.period_u()),
                offset(first.y, previous.y, surface.period_v()),
            );
            parameters.iter_mut().for_each(|point| *point += shift);
            let delta = parameters[0] - previous;
            let error = delta.x.abs() * surface.uder(previous.x, previous.y).magnitude()
                + delta.y.abs() * surface.vder(previous.x, previous.y).magnitude();
            if !error.is_finite() || error > 4.0 * TOLERANCE {
                return Err("inconsistent face-local coedge junction");
            }
            values.extend(parameters.into_iter().skip(1));
        } else {
            values.extend(parameters);
        }
    }
    let (first, last) = (
        *values.first().ok_or("empty parameter loop")?,
        *values.last().ok_or("empty parameter loop")?,
    );
    let delta = last - first;
    let du = surface.uder(first.x, first.y).magnitude();
    let dv = surface.vder(first.x, first.y).magnitude();
    let error = delta.x.abs() * du + delta.y.abs() * dv;
    if !error.is_finite() || du == 0.0 || dv == 0.0 || error > 4.0 * TOLERANCE {
        // A nonzero periodic winding needs a seam-aware face chart, not an
        // invented closing chord in the unwrapped parameter plane.
        return Err("nonzero periodic winding or singular parameter chart");
    }
    *values.last_mut().ok_or("empty parameter loop")? = first;
    Ok(PolylineCurve(values))
}

// UV coordinates have arbitrary units, independent scales and origins. A model
// length tolerance must never be compared to their signed area. Translate and
// normalize both axes, then retain only a numerically resolved orientation.
fn positive_area(poly: &PolylineCurve<Point2>) -> Option<bool> {
    let origin = *poly.first()?;
    let scale = poly
        .iter()
        .fold(Vector2::new(0.0_f64, 0.0_f64), |scale, point| {
            let delta = point - origin;
            Vector2::new(scale.x.max(delta.x.abs()), scale.y.max(delta.y.abs()))
        });
    if scale.x <= 0.0 || scale.y <= 0.0 || !scale.x.is_finite() || !scale.y.is_finite() {
        return None;
    }
    let (mut area, mut correction, mut absolute) = (0.0, 0.0, 0.0);
    for (a, b) in poly
        .iter()
        .zip(poly.iter().cycle().skip(1))
        .take(poly.len())
    {
        let a = a - origin;
        let b = b - origin;
        let left = (a.x / scale.x) * (b.y / scale.y);
        let right = (a.y / scale.y) * (b.x / scale.x);
        super::compensated_add(&mut area, &mut correction, left - right);
        absolute += left.abs() + right.abs();
    }
    (area.is_finite() && area.abs() > absolute * (64.0 * f64::EPSILON)).then_some(area > 0.0)
}

fn contains(poly: &PolylineCurve<Point2>, point: Point2) -> bool {
    let scale = poly
        .iter()
        .fold(Vector2::new(0.0_f64, 0.0_f64), |scale, vertex| {
            let delta = vertex - point;
            Vector2::new(scale.x.max(delta.x.abs()), scale.y.max(delta.y.abs()))
        });
    if scale.x <= 0.0 || scale.y <= 0.0 {
        return false;
    }
    let mut winding = 0_i32;
    for (a, b) in poly
        .iter()
        .zip(poly.iter().cycle().skip(1))
        .take(poly.len())
    {
        let a = a - point;
        let b = b - point;
        let a = Vector2::new(a.x / scale.x, a.y / scale.y);
        let b = Vector2::new(b.x / scale.x, b.y / scale.y);
        let det = a.x * b.y - a.y * b.x;
        let bound = (a.x * b.y).abs() + (a.y * b.x).abs();
        if det.abs() <= 64.0 * f64::EPSILON * bound && a.dot(b) <= 0.0 {
            return true;
        }
        if a.y <= 0.0 && b.y > 0.0 && det > 0.0 {
            winding += 1;
        }
        if a.y > 0.0 && b.y <= 0.0 && det < 0.0 {
            winding -= 1;
        }
    }
    winding != 0
}

#[derive(Clone, Debug)]
struct WireChunk<'a, C> {
    poly: PolylineCurve<Point2>,
    wire: &'a BoundaryWire<Point3, C>,
}

type FaceWithShapesOpStatus<C, S> = (Face<Point3, C, S>, ShapesOpStatus);
fn divide_one_face<C, S>(
    face: &Face<Point3, super::AltCurve<C, S>, S>,
    loops: &Loops<Point3, super::AltCurve<C, S>>,
    tol: f64,
) -> DomainResult<Vec<FaceWithShapesOpStatus<super::AltCurve<C, S>, S>>>
where
    C: super::ShapeOpsCurve<S>,
    S: super::ShapeOpsSurface,
{
    let (mut pre_faces, mut negative_wires) = (Vec::new(), Vec::new());
    let mut map = HashMap::default();
    loops.iter().try_for_each(|wire| {
        let poly = create_parameter_boundary(face, wire, &mut map, tol)?;
        match positive_area(&poly).ok_or("unresolved parameter-loop orientation")? {
            true => pre_faces.push(vec![WireChunk { poly, wire }]),
            false => negative_wires.push(WireChunk { poly, wire }),
        }
        Ok(())
    })?;
    negative_wires.into_iter().try_for_each(|chunk| {
        let point = *chunk.poly.first().ok_or("empty parameter hole")?;
        let mut owner: Option<usize> = None;
        for (index, pre_face) in pre_faces.iter().enumerate() {
            let Some(outer) = pre_face.first() else {
                continue;
            };
            if contains(&outer.poly, point)
                && owner.is_none_or(|old| contains(&pre_faces[old][0].poly, outer.poly.front()))
            {
                owner = Some(index);
            }
        }
        // An orphan hole is a failed split, not an ignorable loop.
        let owner = owner.ok_or("orphan or ambiguous parameter hole")?;
        let outer = &pre_faces[owner][0];
        if outer.wire.deref().inverse() == *chunk.wire.deref() {
            // Cancellation requires the same boundary, not merely equal areas.
            pre_faces[owner].clear();
        } else {
            pre_faces[owner].push(chunk);
        }
        Ok(())
    })?;
    let vec: Vec<_> = pre_faces
        .into_iter()
        .filter(|pre_face| !pre_face.is_empty())
        .map(|pre_face| {
            let surface = face.surface();
            let op = pre_face
                .iter()
                .find(|chunk| chunk.wire.status() != ShapesOpStatus::Unknown);
            let status = match op {
                Some(chunk) => chunk.wire.status(),
                None => ShapesOpStatus::Unknown,
            };
            let wires: Vec<Wire<Point3, super::AltCurve<C, S>>> = pre_face
                .into_iter()
                .map(|chunk| chunk.wire.deref().clone())
                .collect();
            let mut new_face = Face::try_new(wires, surface).ok()?;
            if !face.orientation() {
                new_face.invert();
            }
            Some((new_face, status))
        })
        .collect::<Option<Vec<_>>>()
        .ok_or("native split-face topology")?;
    Ok(vec)
}

pub(super) fn divide_faces<C, S>(
    shell: &Shell<Point3, super::AltCurve<C, S>, S>,
    loops_store: &LoopsStore<Point3, super::AltCurve<C, S>>,
    tol: f64,
) -> FaceDivisionResult<C, S>
where
    C: super::ShapeOpsCurve<S>,
    S: super::ShapeOpsSurface,
{
    if shell.len() != loops_store.len() {
        return Err((0, "input face/loop cardinality"));
    }
    let mut face_index = 0;
    let mut res = FacesClassification::<Point3, super::AltCurve<C, S>, S>::default();
    shell
        .iter()
        .zip(loops_store)
        .enumerate()
        .try_for_each(|(index, (face, loops))| {
            face_index = index;
            if loops
                .iter()
                .all(|wire| wire.status() == ShapesOpStatus::Unknown)
            {
                res.push(face.clone(), ShapesOpStatus::Unknown);
            } else {
                let vec = divide_one_face(face, loops, tol)?;
                vec.into_iter()
                    .for_each(|(face, status)| res.push(face, status));
            }
            Ok(())
        })
        .map_err(|stage| (face_index, stage))?;
    Ok(res)
}
