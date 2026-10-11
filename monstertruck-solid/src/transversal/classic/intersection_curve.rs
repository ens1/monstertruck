//! Classic (0.3.2) marching intersection-curve backend.
//!
//! Extracts the mesh-vs-mesh interference segments for a surface pair, chains
//! them into polylines, and re-samples each polyline onto both surfaces to
//! produce a cleaned [`IntersectionCurve`] leader. Ported verbatim from the
//! published 0.3.2 crate's `transversal::intersection_curve`, except that the
//! parameter-space polylines the 0.3.2 wrapper carried are dropped here: the
//! classic loops-store consumes only the 3D leader.

use monstertruck_core::cgmath64::*;
use monstertruck_geometry::prelude::*;
use monstertruck_meshing::prelude::*;

use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

const MAX_INTERFERENCE_SEGMENTS: usize = 262_144;

type Polyline = PolylineCurve<Point3>;

// Marching produces a segment soup. Only a disjoint union of simple open
// chains and closed loops is admissible here. A branch at a tangent contact or
// singular carrier needs a dedicated resolver, never an arbitrary graph walk.
fn checked_polylines(lines: &[(Point3, Point3)]) -> Option<Vec<Polyline>> {
    if lines.is_empty() {
        return Some(Vec::new());
    }
    if lines.len() > MAX_INTERFERENCE_SEGMENTS {
        return None;
    }
    let origin = lines[0].0;
    if lines
        .iter()
        .any(|(a, b)| !super::finite_point(*a) || !super::finite_point(*b))
    {
        return None;
    }
    let scale = lines
        .iter()
        .flat_map(|(a, b)| [a, b])
        .fold(0.0_f64, |scale, point| {
            let delta = *point - origin;
            scale
                .max(delta.x.abs())
                .max(delta.y.abs())
                .max(delta.z.abs())
        });
    if scale * (128.0 * f64::EPSILON) > TOLERANCE {
        return None;
    }
    let join = TOLERANCE;
    let mut buckets: HashMap<[i64; 3], Vec<usize>> = HashMap::default();
    let mut points: Vec<Point3> = Vec::new();
    let mut edges: Vec<(usize, usize)> = Vec::new();
    let mut edge_set = HashSet::default();
    let mut node = |point: Point3| -> Option<usize> {
        let relative = (point - origin) / join;
        let values = [relative.x.floor(), relative.y.floor(), relative.z.floor()];
        if values
            .iter()
            .any(|value| !value.is_finite() || value.abs() >= (i64::MAX - 2) as f64)
        {
            return None;
        }
        let key = values.map(|value| value as i64);
        let mut found = None;
        for x in -1..=1 {
            for y in -1..=1 {
                for z in -1..=1 {
                    if let Some(candidates) = buckets.get(&[key[0] + x, key[1] + y, key[2] + z]) {
                        for &index in candidates {
                            if points[index].distance(point) <= join {
                                if found.is_some_and(|old| old != index) {
                                    return None;
                                }
                                found = Some(index);
                            }
                        }
                    }
                }
            }
        }
        if let Some(index) = found {
            return Some(index);
        }
        let index = points.len();
        points.push(point);
        buckets.entry(key).or_default().push(index);
        Some(index)
    };
    for &(a, b) in lines {
        let a = node(a)?;
        let b = node(b)?;
        if a != b {
            let key = (a.min(b), a.max(b));
            if edge_set.insert(key) {
                edges.push(key);
            }
        }
    }
    let mut adjacency = vec![Vec::new(); points.len()];
    for (index, &(a, b)) in edges.iter().enumerate() {
        adjacency[a].push(index);
        adjacency[b].push(index);
        if adjacency[a].len() > 2 || adjacency[b].len() > 2 {
            return None;
        }
    }
    let mut used = vec![false; edges.len()];
    let mut result = Vec::new();
    // Open ends first so an open chain is never split by beginning in its middle.
    for start in (0..points.len())
        .filter(|&index| adjacency[index].len() == 1)
        .chain((0..points.len()).filter(|&index| adjacency[index].len() == 2))
    {
        if adjacency[start].iter().all(|&edge| used[edge]) {
            continue;
        }
        let mut current = start;
        let mut line = vec![points[start]];
        while let Some(&edge) = adjacency[current].iter().find(|&&edge| !used[edge]) {
            used[edge] = true;
            let (a, b) = edges[edge];
            current = if a == current { b } else { a };
            line.push(points[current]);
            if current == start {
                break;
            }
        }
        if line.len() < 2 || (current == start && line.len() < 4) {
            return None;
        }
        result.push(PolylineCurve(line));
    }
    used.iter().all(|used| *used).then_some(result)
}

/// Re-sample a search branch onto both native surfaces. Every accepted point
/// must agree with both carriers and remain near its original search branch.
fn build_intersection_curve<S>(
    surface0: S,
    surface1: S,
    poly: Polyline,
    tol: f64,
) -> Option<IntersectionCurve<Polyline, S, S>>
where
    S: ParametricSurface3D
        + Clone
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    let ic = IntersectionCurve::new(surface0.clone(), surface1.clone(), poly);
    let raw = ic.leader();
    let len = raw.len();
    if len < 2 {
        return None;
    }
    let validate = |triple: (Point3, Point2, Point2), seed: Point3| {
        let (point, uv0, uv1) = triple;
        if !super::finite_point(point)
            || !uv0.x.is_finite()
            || !uv0.y.is_finite()
            || !uv1.x.is_finite()
            || !uv1.y.is_finite()
        {
            return None;
        }
        let q0 = surface0.subs(uv0.x, uv0.y);
        let q1 = surface1.subs(uv1.x, uv1.y);
        if !super::finite_point(q0) || !super::finite_point(q1) {
            return None;
        }
        let residual = TOLERANCE;
        if q0.distance(q1) > residual
            || point.distance(q0) > residual
            || point.distance(q1) > residual
            || point.distance(seed) > 4.0 * tol
        {
            return None;
        }
        Some(point)
    };
    let recover = |t: f64| {
        let seed = raw.subs(t);
        ic.search_triple(t, 100)
            .and_then(|triple| validate(triple, seed))
            .or_else(|| {
                ic.search_nearest_point(seed, None, None, 100)
                    .and_then(|triple| validate(triple, seed))
            })
            .or_else(|| {
                let uv0 = surface0
                    .search_parameter(seed, None, 100)
                    .or_else(|| surface0.search_nearest_parameter(seed, None, 100))?;
                let uv1 = surface1
                    .search_parameter(seed, None, 100)
                    .or_else(|| surface1.search_nearest_parameter(seed, None, 100))?;
                let q0 = surface0.subs(uv0.0, uv0.1);
                let q1 = surface1.subs(uv1.0, uv1.1);
                validate((q0.midpoint(q1), uv0.into(), uv1.into()), seed)
            })
    };
    let mut polyline: Polyline = PolylineCurve(Vec::with_capacity(len));
    for index in 0..len {
        let point = recover(index as f64)?;
        if polyline
            .last()
            .is_some_and(|previous| previous.distance(point) <= f64::EPSILON)
        {
            return None;
        }
        polyline.push(point);
        if index + 1 < len {
            // Guard against a leader segment passing between different native
            // branches or a projection singularity hidden between its vertices.
            recover(index as f64 + 0.5)?;
        }
    }
    if raw[0].near(&raw[len - 1]) {
        polyline[len - 1] = polyline[0];
    }
    Some(IntersectionCurve::new(surface0, surface1, polyline))
}

const MAX_TRIANGLE_PAIR_WORK: usize = 64_000_000;

struct SearchTriangle {
    vertices: [Point3; 3],
    normal: Vector3,
    bounds: BoundingBox<Point3>,
}

fn search_triangles(mesh: &PolygonMesh) -> Option<Vec<SearchTriangle>> {
    let mut result = Vec::new();
    for polygon in mesh.face_iter() {
        for index in 2..polygon.len() {
            let vertices = [polygon[0], polygon[index - 1], polygon[index]]
                .map(|vertex| mesh.positions()[vertex.pos]);
            let e0 = vertices[1] - vertices[0];
            let e1 = vertices[2] - vertices[0];
            let scale = [
                e0.x.abs(),
                e0.y.abs(),
                e0.z.abs(),
                e1.x.abs(),
                e1.y.abs(),
                e1.z.abs(),
            ]
            .into_iter()
            .fold(0.0_f64, f64::max);
            if scale == 0.0 {
                continue;
            }
            if !scale.is_finite() {
                return None;
            }
            let normal = (e0 / scale).cross(e1 / scale);
            let length = normal.magnitude();
            if length == 0.0 {
                continue;
            }
            if !length.is_finite() {
                return None;
            }
            result.push(SearchTriangle {
                vertices,
                normal: normal / length,
                bounds: vertices.into_iter().collect(),
            });
        }
    }
    Some(result)
}

fn disjoint_bounds(a: BoundingBox<Point3>, b: BoundingBox<Point3>) -> bool {
    (0..3).any(|axis| a.max()[axis] < b.min()[axis] || b.max()[axis] < a.min()[axis])
}

fn section_points(vertices: [Point3; 3], distances: [f64; 3], epsilon: f64) -> Vec<Point3> {
    let mut result = Vec::with_capacity(4);
    for index in 0..3 {
        let next = (index + 1) % 3;
        if distances[index].abs() <= epsilon {
            result.push(vertices[index]);
        }
        // An endpoint admitted as lying on the plane must not also create a
        // tiny duplicate crossing on the same edge.
        if distances[index].abs() > epsilon
            && distances[next].abs() > epsilon
            && distances[index].is_sign_positive() != distances[next].is_sign_positive()
        {
            let weight = distances[index] / (distances[index] - distances[next]);
            result.push(vertices[index] + (vertices[next] - vertices[index]) * weight);
        }
    }
    result
}

fn triangle_interference(a: &SearchTriangle, b: &SearchTriangle) -> Option<(Point3, Point3)> {
    let direction = a.normal.cross(b.normal);
    let length = direction.magnitude();
    if !length.is_finite() || length <= 64.0 * f64::EPSILON {
        return None;
    }
    let direction = direction / length;
    let origin = a.vertices[0];
    let local_extent = a
        .vertices
        .iter()
        .chain(b.vertices.iter())
        .fold(0.0_f64, |scale, point| {
            let relative = *point - origin;
            scale
                .max(relative.x.abs())
                .max(relative.y.abs())
                .max(relative.z.abs())
        });
    let epsilon = (64.0 * f64::EPSILON * local_extent).max(f64::MIN_POSITIVE);
    let da = a.vertices.map(|p| (p - b.vertices[0]).dot(b.normal));
    let db = b.vertices.map(|p| (p - a.vertices[0]).dot(a.normal));
    let one_side = |d: &[f64; 3]| d.iter().all(|d| *d > epsilon) || d.iter().all(|d| *d < -epsilon);
    if one_side(&da) || one_side(&db) {
        return None;
    }
    let interval = |points: Vec<Point3>| {
        let mut points: Vec<_> = points
            .into_iter()
            .map(|point| (direction.dot(point - origin), point))
            .collect();
        points.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
        Some((*points.first()?, *points.last()?))
    };
    let (a0, a1) = interval(section_points(a.vertices, da, epsilon))?;
    let (b0, b1) = interval(section_points(b.vertices, db, epsilon))?;
    let start = if a0.0 >= b0.0 { a0 } else { b0 };
    let end = if a1.0 <= b1.0 { a1 } else { b1 };
    (end.0 - start.0 > epsilon).then_some((start.1, end.1))
}

fn checked_interferences(
    polygon0: &PolygonMesh,
    polygon1: &PolygonMesh,
    work: &mut usize,
) -> Option<Vec<(Point3, Point3)>> {
    let triangles0 = search_triangles(polygon0)?;
    let mut triangles1 = search_triangles(polygon1)?;
    // A conservative interval sweep reduces work on localized face crossings.
    // The visit budget includes rejected intervals, keeping degenerate cases bounded.
    let bounds: BoundingBox<Point3> = polygon0
        .positions()
        .iter()
        .chain(polygon1.positions())
        .copied()
        .collect();
    let diagonal = bounds.max() - bounds.min();
    let axis = if diagonal.x >= diagonal.y && diagonal.x >= diagonal.z {
        0
    } else if diagonal.y >= diagonal.z {
        1
    } else {
        2
    };
    triangles1.sort_unstable_by(|a, b| a.bounds.min()[axis].total_cmp(&b.bounds.min()[axis]));
    let mut result = Vec::new();
    for a in &triangles0 {
        for b in &triangles1 {
            *work = work.checked_add(1)?;
            if *work > MAX_TRIANGLE_PAIR_WORK {
                return None;
            }
            if b.bounds.min()[axis] > a.bounds.max()[axis] {
                break;
            }
            if !disjoint_bounds(a.bounds, b.bounds)
                && let Some(line) = triangle_interference(a, b)
            {
                if !super::finite_point(line.0) || !super::finite_point(line.1) {
                    return None;
                }
                if result.len() == MAX_INTERFERENCE_SEGMENTS {
                    return None;
                }
                result.push(line);
            }
        }
    }
    Some(result)
}

type IntersectionTuple<S> = (Polyline, IntersectionCurve<Polyline, S, S>);

/// Marching intersection curves between two trimmed surfaces, keyed by their
/// tessellated polygon meshes. Returns, per intersection branch, the raw
/// interference polyline and the cleaned polyline intersection curve.
pub(super) fn intersection_curves<S>(
    surface0: S,
    polygon0: &PolygonMesh,
    surface1: S,
    polygon1: &PolygonMesh,
    tol: f64,
    work: &mut usize,
) -> Option<Vec<IntersectionTuple<S>>>
where
    S: ParametricSurface3D
        + Clone
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    let bounds0: BoundingBox<Point3> = polygon0.positions().iter().copied().collect();
    let bounds1: BoundingBox<Point3> = polygon1.positions().iter().copied().collect();
    if disjoint_bounds(bounds0, bounds1) {
        return Some(Vec::new());
    }
    let interferences = checked_interferences(polygon0, polygon1, work)?;
    checked_polylines(&interferences)?
        .into_iter()
        .map(|polyline| {
            let curve = build_intersection_curve(
                surface0.clone(),
                surface1.clone(),
                polyline.clone(),
                tol,
            )?;
            Some((polyline, curve))
        })
        .collect()
}
