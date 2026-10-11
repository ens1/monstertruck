//! Native straight-boundary planar arrangements. Infinite supporting lines
//! partition each complete trimmed face; only analytic inside/outside decisions
//! select cells. Internal subdivision edges are cancelled before faces are
//! rebuilt, and a common topology pool conforms every remaining shared edge.
//! Curved or unproven inputs remain on the surface-intersection backend.

use super::{ClassicResult, ShapeOpsCurve, ShapeOpsError, ShapeOpsSurface, ShellOrientationHints};
use monstertruck_geometry::prelude::*;
use monstertruck_topology::*;
use rustc_hash::FxHashMap as HashMap;

const EPS: f64 = TOLERANCE;
const ANGULAR: f64 = 128.0 * f64::EPSILON;
const MAX_FACES: usize = 512;
const MAX_EDGES: usize = 16_384;
const MAX_CELLS: usize = 65_536;
const MAX_VERTICES: usize = 65_536;
const MAX_WORK: usize = 64_000_000;

fn failed(stage: &'static str) -> ShapeOpsError { ShapeOpsError::PlanarArrangementFailed { stage } }
fn budget(work: &mut usize, count: usize) -> ClassicResult<()> {
    *work = work
        .checked_add(count)
        .ok_or(ShapeOpsError::ResourceLimit {
            stage: "native planar arrangement",
            limit: MAX_WORK,
        })?;
    if *work > MAX_WORK {
        return Err(ShapeOpsError::ResourceLimit {
            stage: "native planar arrangement",
            limit: MAX_WORK,
        });
    }
    Ok(())
}
fn finite2(p: Point2) -> bool { p.x.is_finite() && p.y.is_finite() }
fn cross(a: Vector2, b: Vector2) -> f64 { a.x * b.y - a.y * b.x }
fn key(a: usize, b: usize) -> (usize, usize) { (a.min(b), a.max(b)) }

#[derive(Clone, Copy)]
struct Frame {
    origin: Point3,
    x: Vector3,
    y: Vector3,
    normal: Vector3,
}
impl Frame {
    fn project(self, p: Point3) -> Point2 {
        let d = p - self.origin;
        Point2::new(d.dot(self.x), d.dot(self.y))
    }
    fn lift(self, p: Point2) -> Point3 { self.origin + self.x * p.x + self.y * p.y }
    fn signed(self, p: Point3) -> f64 { (p - self.origin).dot(self.normal) }
}

struct FaceData<C, S> {
    source: Face<Point3, C, S>,
    frame: Frame,
    loops: Vec<Vec<Point2>>,
    operand: usize,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Location {
    Outside,
    Inside,
    Boundary,
}

fn on_segment(p: Point2, a: Point2, b: Point2) -> bool {
    let delta = b - a;
    let length = delta.magnitude();
    if !length.is_finite() || length <= EPS {
        return false;
    }
    let d = p - a;
    let along = d.dot(delta / length);
    cross(d, delta / length).abs() <= EPS && along >= -EPS && along <= length + EPS
}

fn location(loops: &[Vec<Point2>], p: Point2, work: &mut usize) -> ClassicResult<Location> {
    if !finite2(p) {
        return Err(failed("non-finite planar query"));
    }
    let mut inside = false;
    for polygon in loops {
        budget(work, polygon.len())?;
        for i in 0..polygon.len() {
            let (a, b) = (polygon[i], polygon[(i + 1) % polygon.len()]);
            if on_segment(p, a, b) {
                return Ok(Location::Boundary);
            }
            if (a.y > p.y) != (b.y > p.y) {
                // Translate first; no intersection with a world-origin ray is
                // used, and the half-open endpoint rule counts a vertex once.
                let a = a - p;
                let b = b - p;
                let x = a.x + (b.x - a.x) * (-a.y / (b.y - a.y));
                if !x.is_finite() {
                    return Err(failed("unresolved planar ray"));
                }
                if x > 0.0 {
                    inside = !inside;
                }
            }
        }
    }
    Ok(if inside {
        Location::Inside
    } else {
        Location::Outside
    })
}

fn area(polygon: &[Point2], work: &mut usize) -> ClassicResult<f64> {
    budget(work, polygon.len())?;
    let origin = polygon[0];
    let (mut sum, mut correction, mut absolute) = (0.0, 0.0, 0.0);
    for i in 1..polygon.len() - 1 {
        let value = cross(polygon[i] - origin, polygon[i + 1] - origin);
        super::compensated_add(&mut sum, &mut correction, value);
        absolute += value.abs();
    }
    if !sum.is_finite() || sum.abs() <= absolute * ANGULAR || sum.abs() <= EPS * EPS {
        return Err(failed("collapsed or unresolved planar cell area"));
    }
    Ok(sum * 0.5)
}

fn validate_loops(loops: &[Vec<Point2>], work: &mut usize) -> ClassicResult<()> {
    let mut segments = Vec::new();
    for (loop_index, polygon) in loops.iter().enumerate() {
        if (area(polygon, work)? > 0.0) != (loop_index == 0) {
            return Err(failed("inconsistent outer/hole loop orientation"));
        }
        for i in 0..polygon.len() {
            segments.push((loop_index, i, polygon[i], polygon[(i + 1) % polygon.len()]));
        }
    }
    for (i, &(li, ei, a, b)) in segments.iter().enumerate() {
        let delta = b - a;
        let length = delta.magnitude();
        let u = delta / length;
        for &(lj, ej, c, d) in &segments[i + 1..] {
            budget(work, 1)?;
            let other = d - c;
            let span = other.magnitude();
            let v = other / span;
            let adjacent =
                li == lj && ((ei + 1) % loops[li].len() == ej || (ej + 1) % loops[li].len() == ei);
            let determinant = cross(u, v);
            if determinant.abs() <= ANGULAR {
                if cross(c - a, u).abs() <= EPS {
                    let start = (c - a).dot(u);
                    let end = (d - a).dot(u);
                    let overlap = length.min(start.max(end)) - 0.0_f64.max(start.min(end));
                    if overlap > EPS || (overlap >= -EPS && !adjacent) {
                        return Err(failed("overlapping or touching input trim segments"));
                    }
                }
            } else {
                let t = cross(c - a, v) / determinant;
                let s = cross(c - a, u) / determinant;
                if !t.is_finite() || !s.is_finite() {
                    return Err(failed("unresolved input trim intersection"));
                }
                if t >= -EPS && t <= length + EPS && s >= -EPS && s <= span + EPS && !adjacent {
                    return Err(failed("self-intersecting or touching input trim loops"));
                }
            }
        }
    }
    for i in 1..loops.len() {
        if location(std::slice::from_ref(&loops[0]), loops[i][0], work)? != Location::Inside {
            return Err(failed("input hole lies outside its outer loop"));
        }
        for j in 1..loops.len() {
            if i != j
                && location(std::slice::from_ref(&loops[j]), loops[i][0], work)?
                    != Location::Outside
            {
                return Err(failed("nested or overlapping input holes"));
            }
        }
    }
    Ok(())
}

fn admit<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    shell: &Shell<Point3, C, S>,
    operand: usize,
    work: &mut usize,
) -> ClassicResult<Option<Vec<FaceData<C, S>>>> {
    let mut result = Vec::new();
    let mut count = 0;
    for face in shell {
        let surface = face.oriented_surface();
        let Some(AnalyticSurfaceKind::Plane(plane)) = surface.try_into_analytic_surface_kind()
        else {
            return Ok(None);
        };
        let normal = plane.normal();
        if !super::finite_point(plane.origin()) || !normal.magnitude().is_finite() {
            return Err(failed("invalid native plane"));
        }
        let Some(first) = face.vertex_iter().next() else {
            return Ok(None);
        };
        let origin = first.point();
        let Some((u, v)) = surface.search_parameter(origin, None, 100) else {
            return Ok(None);
        };
        let oriented = surface.normal(u, v);
        let length = oriented.magnitude();
        if !length.is_finite() || length <= 0.0 {
            return Ok(None);
        }
        let normal = oriented / length;
        let x = normal
            .cross(if normal.x.abs() < 0.8 {
                Vector3::unit_x()
            } else {
                Vector3::unit_y()
            })
            .normalize();
        let frame = Frame {
            origin,
            x,
            y: normal.cross(x),
            normal,
        };
        // Analytic identity is rechecked at fixed kernel precision over the
        // complete positive control hull, rather than accepting a fitted plane.
        let Some(net) = surface.try_into_homogeneous_bspline_surface() else {
            return Ok(None);
        };
        if net.control_points().iter().map(Vec::len).sum::<usize>() > MAX_EDGES {
            return Ok(None);
        }
        for h in net.control_points().iter().flatten() {
            if !h.w.is_finite() || h.w <= 0.0 {
                return Ok(None);
            }
            let p = Point3::from_vec(h.truncate() / h.w);
            if !super::finite_point(p) || frame.signed(p).abs() > EPS {
                return Ok(None);
            }
        }
        let mut loops = Vec::new();
        for wire in face.boundaries() {
            if wire.len() < 3 || !wire.is_closed() {
                return Ok(None);
            }
            let mut polygon = Vec::new();
            for edge in wire {
                count += 1;
                if count > MAX_EDGES {
                    return Ok(None);
                }
                let Some((a, b)) = edge.oriented_curve().exact_straight_segment() else {
                    return Ok(None);
                };
                if C::from_straight_segment(a, b).is_none() {
                    return Ok(None);
                }
                if !super::finite_point(a)
                    || !super::finite_point(b)
                    || a.distance(edge.front().point()) > EPS
                    || b.distance(edge.back().point()) > EPS
                    || a.distance(b) <= EPS
                    || frame.signed(a).abs() > EPS
                    || frame.signed(b).abs() > EPS
                {
                    return Err(failed("invalid exact straight boundary"));
                }
                let p = frame.project(edge.front().point());
                let scale = (a - origin).magnitude().max(a.to_vec().magnitude());
                if !finite2(p) || scale * ANGULAR > EPS {
                    return Err(failed("coordinates cannot resolve kernel precision"));
                }
                polygon.push(p);
            }
            area(&polygon, work)?;
            loops.push(polygon);
        }
        if loops.is_empty() {
            return Ok(None);
        }
        result.push(FaceData {
            source: face.clone(),
            frame,
            loops,
            operand,
        });
    }
    Ok(Some(result))
}

#[derive(Clone, Copy)]
struct Line2 {
    normal: Vector2,
    offset: f64,
}
impl Line2 {
    fn signed(self, p: Point2) -> f64 { p.to_vec().dot(self.normal) - self.offset }
}
fn add_line(lines: &mut Vec<Line2>, mut line: Line2, work: &mut usize) -> ClassicResult<()> {
    let length = line.normal.magnitude();
    if !length.is_finite() || !line.offset.is_finite() || length <= 0.0 {
        return Err(failed("invalid partition line"));
    }
    line.normal /= length;
    line.offset /= length;
    if line.normal.x < -ANGULAR || (line.normal.x.abs() <= ANGULAR && line.normal.y < 0.0) {
        line.normal = -line.normal;
        line.offset = -line.offset;
    }
    budget(work, lines.len())?;
    if lines.iter().any(|old| {
        (old.normal - line.normal).magnitude() <= ANGULAR && (old.offset - line.offset).abs() <= EPS
    }) {
        return Ok(());
    }
    if lines.len() >= MAX_EDGES {
        return Err(failed("too many partition lines"));
    }
    lines.push(line);
    Ok(())
}

fn split(
    polygon: &[Point2],
    line: Line2,
    work: &mut usize,
) -> ClassicResult<Option<[Vec<Point2>; 2]>> {
    let distances = polygon.iter().map(|p| line.signed(*p)).collect::<Vec<_>>();
    if distances.iter().any(|d| !d.is_finite()) {
        return Err(failed("non-finite partition distance"));
    }
    if !distances.iter().any(|d| *d > EPS) || !distances.iter().any(|d| *d < -EPS) {
        return Ok(None);
    }
    let mut result = [Vec::new(), Vec::new()];
    for i in 0..polygon.len() {
        let j = (i + 1) % polygon.len();
        let (a, b) = (polygon[i], polygon[j]);
        let (d, e) = (distances[i], distances[j]);
        if d >= -EPS {
            result[0].push(a);
        }
        if d <= EPS {
            result[1].push(a);
        }
        if (d > EPS && e < -EPS) || (d < -EPS && e > EPS) {
            let p = a + (b - a) * (d / (d - e));
            if !finite2(p) || line.signed(p).abs() > EPS {
                return Err(failed("unresolved partition intersection"));
            }
            result[0].push(p);
            result[1].push(p);
        }
    }
    for half in &mut result {
        half.dedup_by(|a, b| a.distance(*b) <= EPS);
        if half.len() > 1 && half[0].distance(*half.last().unwrap()) <= EPS {
            half.pop();
        }
        if half.len() < 3 {
            return Err(failed("partition collapses a bounded cell"));
        }
        area(half, work)?;
    }
    Ok(Some(result))
}

fn cells<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    faces: &[FaceData<C, S>],
    index: usize,
    work: &mut usize,
) -> ClassicResult<Vec<Vec<Point2>>> {
    let face = &faces[index];
    let mut lines = Vec::new();
    for other in faces {
        let cross = face.frame.normal.cross(other.frame.normal).magnitude();
        if cross <= ANGULAR {
            if face.frame.signed(other.frame.origin).abs() > EPS {
                continue;
            }
            for polygon in &other.loops {
                for i in 0..polygon.len() {
                    let a = face.frame.project(other.frame.lift(polygon[i]));
                    let b = face
                        .frame
                        .project(other.frame.lift(polygon[(i + 1) % polygon.len()]));
                    let delta = b - a;
                    let normal = Vector2::new(-delta.y, delta.x);
                    add_line(
                        &mut lines,
                        Line2 {
                            normal,
                            offset: normal.dot(a.to_vec()),
                        },
                        work,
                    )?;
                }
            }
        } else {
            add_line(
                &mut lines,
                Line2 {
                    normal: Vector2::new(
                        other.frame.normal.dot(face.frame.x),
                        other.frame.normal.dot(face.frame.y),
                    ),
                    offset: -other.frame.signed(face.frame.origin),
                },
                work,
            )?;
        }
    }
    let points = face.loops.iter().flatten();
    let mut min = Point2::new(f64::INFINITY, f64::INFINITY);
    let mut max = Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
    for p in points {
        min.x = min.x.min(p.x);
        min.y = min.y.min(p.y);
        max.x = max.x.max(p.x);
        max.y = max.y.max(p.y);
    }
    let mut cells = vec![vec![
        min,
        Point2::new(max.x, min.y),
        max,
        Point2::new(min.x, max.y),
    ]];
    area(&cells[0], work)?;
    for line in lines {
        let mut next = Vec::new();
        for cell in cells {
            budget(work, cell.len())?;
            if let Some([a, b]) = split(&cell, line, work)? {
                next.push(a);
                next.push(b);
            } else {
                next.push(cell);
            }
            if next.len() > MAX_CELLS {
                return Err(failed("too many planar cells"));
            }
        }
        cells = next;
    }
    let mut result = Vec::new();
    for cell in cells {
        let center =
            Point2::from_vec(cell.iter().map(|p| p.to_vec()).sum::<Vector2>() / cell.len() as f64);
        match location(&face.loops, center, work)? {
            Location::Inside => result.push(cell),
            Location::Outside => {}
            Location::Boundary => return Err(failed("cell has no resolved interior")),
        }
    }
    Ok(result)
}

fn inverted<C, S>(faces: &[FaceData<C, S>], operand: usize) -> ClassicResult<bool> {
    let origin = faces
        .iter()
        .find(|f| f.operand == operand)
        .ok_or_else(|| failed("empty operand"))?
        .frame
        .origin;
    let (mut volume, mut correction, mut absolute) = (0.0, 0.0, 0.0);
    for face in faces.iter().filter(|f| f.operand == operand) {
        for polygon in &face.loops {
            let a = face.frame.lift(polygon[0]) - origin;
            for i in 1..polygon.len() - 1 {
                let b = face.frame.lift(polygon[i]) - origin;
                let c = face.frame.lift(polygon[i + 1]) - origin;
                let value = a.dot(b.cross(c));
                super::compensated_add(&mut volume, &mut correction, value);
                absolute += value.abs();
            }
        }
    }
    if !volume.is_finite() || volume.abs() <= absolute * ANGULAR {
        return Err(failed("unresolved material orientation"));
    }
    Ok(volume < 0.0)
}

fn inside<C, S>(
    faces: &[FaceData<C, S>],
    operand: usize,
    point: Point3,
    inverted: bool,
    work: &mut usize,
) -> ClassicResult<bool> {
    // A finite set of deterministic ray directions is only a degeneracy
    // avoidance device. Each accepted classification uses every exact face.
    for attempt in 0..48 {
        let t = attempt as f64 + 0.371;
        let ray = Vector3::new(
            1.0,
            (t * 1.618033988749895).sin(),
            (t * 2.414213562373095).cos(),
        )
        .normalize();
        let mut count = 0_i32;
        let mut ambiguous = false;
        for face in faces.iter().filter(|f| f.operand == operand) {
            budget(work, 1)?;
            let denominator = face.frame.normal.dot(ray);
            if denominator.abs() <= ANGULAR {
                ambiguous = true;
                break;
            }
            let t = -face.frame.signed(point) / denominator;
            if !t.is_finite() {
                return Err(failed("non-finite analytic ray intersection"));
            }
            if t <= EPS {
                continue;
            }
            let p = point + ray * t;
            if !super::finite_point(p) {
                return Err(failed("non-finite analytic ray point"));
            }
            match location(&face.loops, face.frame.project(p), work)? {
                Location::Boundary => {
                    ambiguous = true;
                    break;
                }
                Location::Inside => count += if denominator > 0.0 { 1 } else { -1 },
                Location::Outside => {}
            }
        }
        if !ambiguous {
            return match (inverted, count) {
                (false, 0) | (true, -1) => Ok(false),
                (false, 1) | (true, 0) => Ok(true),
                _ => Err(failed("inconsistent analytic material winding")),
            };
        }
    }
    Err(failed("no nondegenerate analytic classification ray"))
}

fn selection<C, S>(
    faces: &[FaceData<C, S>],
    index: usize,
    point: Point3,
    inverted: [bool; 2],
    work: &mut usize,
) -> ClassicResult<[bool; 2]> {
    let source = &faces[index];
    let mut coincidence = None;
    for other in faces.iter().filter(|f| f.operand != source.operand) {
        budget(work, 1)?;
        if other.frame.signed(point).abs() > EPS {
            continue;
        }
        match location(&other.loops, other.frame.project(point), work)? {
            Location::Outside => continue,
            Location::Boundary => return Err(failed("cell touches an unresolved opposing trim")),
            Location::Inside => {}
        }
        if source.frame.normal.cross(other.frame.normal).magnitude() > ANGULAR {
            return Err(failed("cell interior lies on a transverse plane"));
        }
        let same = source.frame.normal.dot(other.frame.normal) > 0.0;
        if coincidence.is_some() {
            return Err(failed("overlapping input coplanar face interiors"));
        }
        coincidence = Some(same);
    }
    if let Some(same) = coincidence {
        // Same-oriented coincident material boundaries belong to both results,
        // once. Opposed boundaries are interior to union and absent from AND.
        return Ok([same && source.operand == 0; 2]);
    }
    let other = 1 - source.operand;
    let inside = inside(faces, other, point, inverted[other], work)?;
    Ok([inside, !inside])
}

struct Cell {
    source: usize,
    selected: [bool; 2],
    points: Vec<Point3>,
}
struct Vertices {
    points: Vec<Vertex<Point3>>,
    members: Vec<Vec<Point3>>,
}
impl Vertices {
    fn insert(
        &mut self,
        point: Point3,
        preferred: Option<Vertex<Point3>>,
        work: &mut usize,
    ) -> ClassicResult<usize> {
        let mut found = None;
        for (i, members) in self.members.iter().enumerate() {
            budget(work, members.len())?;
            if members.iter().any(|p| p.distance(point) <= EPS) {
                if found.is_some() {
                    return Err(failed("ambiguous nontransitive vertex weld"));
                }
                if members.iter().any(|p| p.distance(point) > EPS) {
                    return Err(failed(
                        "nontransitive vertex cluster exceeds kernel precision",
                    ));
                }
                found = Some(i);
            }
        }
        if let Some(index) = found {
            if !self.members[index].contains(&point) {
                self.members[index].push(point);
            }
            return Ok(index);
        }
        if self.points.len() >= MAX_VERTICES {
            return Err(failed("too many arrangement vertices"));
        }
        let index = self.points.len();
        self.points
            .push(preferred.unwrap_or_else(|| Vertex::new(point)));
        self.members.push(vec![point]);
        Ok(index)
    }
}

fn conformed_segment(
    a: usize,
    b: usize,
    vertices: &Vertices,
    work: &mut usize,
) -> ClassicResult<Vec<usize>> {
    let start = vertices.points[a].point();
    let delta = vertices.points[b].point() - start;
    let length = delta.magnitude();
    if !length.is_finite() || length <= EPS {
        return Err(failed("collapsed arrangement edge"));
    }
    let direction = delta / length;
    let mut points = vec![(0.0, a), (length, b)];
    budget(work, vertices.points.len())?;
    for (index, vertex) in vertices.points.iter().enumerate() {
        if index == a || index == b {
            continue;
        }
        let local = vertex.point() - start;
        let distance = local.dot(direction);
        if distance > EPS && distance < length - EPS && local.cross(direction).magnitude() <= EPS {
            points.push((distance, index));
        }
    }
    points.sort_by(|a, b| a.0.total_cmp(&b.0));
    if points.windows(2).any(|p| p[1].0 - p[0].0 <= EPS) {
        return Err(failed("unresolved edge split order"));
    }
    Ok(points.into_iter().map(|(_, i)| i).collect())
}

fn loops_from_edges(
    edges: &HashMap<(usize, usize), (usize, usize)>,
) -> ClassicResult<Vec<Vec<usize>>> {
    let mut next = HashMap::default();
    let mut previous = HashMap::default();
    for &(a, b) in edges.values() {
        if next.insert(a, b).is_some() || previous.insert(b, a).is_some() {
            return Err(failed("nonmanifold planar region boundary"));
        }
    }
    let mut result = Vec::new();
    while let Some(&start) = next.keys().next() {
        let mut current = start;
        let mut wire = Vec::new();
        loop {
            wire.push(current);
            current = next
                .remove(&current)
                .ok_or_else(|| failed("open planar region boundary"))?;
            if current == start {
                break;
            }
        }
        if wire.len() < 3 {
            return Err(failed("collapsed planar boundary loop"));
        }
        result.push(wire);
    }
    Ok(result)
}

fn assemble<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    faces: &[FaceData<C, S>],
    cells: &[Cell],
    work: &mut usize,
) -> ClassicResult<[Shell<Point3, C, S>; 2]> {
    let mut vertices = Vertices {
        points: Vec::new(),
        members: Vec::new(),
    };
    // Input vertices are immutable candidates; retaining them never mutates or
    // welds the original operands in place.
    for face in faces {
        for vertex in face.source.vertex_iter() {
            vertices.insert(vertex.point(), Some(vertex.clone()), work)?;
        }
    }
    let mut cell_vertices = Vec::new();
    for cell in cells {
        let mut points = Vec::new();
        for &point in &cell.points {
            points.push(vertices.insert(point, None, work)?);
        }
        cell_vertices.push(points);
    }
    let mut curves: HashMap<(usize, usize), Edge<Point3, C>> = HashMap::default();
    for face in faces {
        for edge in face.source.edge_iter() {
            let a = vertices.insert(edge.absolute_front().point(), None, work)?;
            let b = vertices.insert(edge.absolute_back().point(), None, work)?;
            if vertices.points[a].id() == edge.absolute_front().id()
                && vertices.points[b].id() == edge.absolute_back().id()
            {
                curves.entry(key(a, b)).or_insert_with(|| {
                    if edge.orientation() {
                        edge
                    } else {
                        edge.inverse()
                    }
                });
            }
        }
    }
    type Boundary = HashMap<(usize, usize), (usize, usize)>;
    let mut groups: Vec<[Boundary; 2]> = (0..faces.len())
        .map(|_| [HashMap::default(), HashMap::default()])
        .collect();
    for (cell, points) in cells.iter().zip(cell_vertices) {
        for i in 0..points.len() {
            let pieces =
                conformed_segment(points[i], points[(i + 1) % points.len()], &vertices, work)?;
            for pair in pieces.windows(2) {
                let (a, b) = (pair[0], pair[1]);
                for (selected, group) in cell.selected.iter().zip(&mut groups[cell.source]) {
                    if !selected {
                        continue;
                    }
                    if let Some((old_a, old_b)) = group.remove(&key(a, b)) {
                        if old_a != b || old_b != a {
                            return Err(failed("overlapping selected planar cells"));
                        }
                    } else {
                        group.insert(key(a, b), (a, b));
                    }
                }
            }
        }
    }
    let mut output = [Shell::new(), Shell::new()];
    for (face, group) in faces.iter().zip(groups) {
        for (operation, group) in group.into_iter().enumerate() {
            budget(work, group.len())?;
            let wires = loops_from_edges(&group)?;
            let mut outer = Vec::new();
            let mut holes = Vec::new();
            for wire in wires {
                let polygon = wire
                    .iter()
                    .map(|i| face.frame.project(vertices.points[*i].point()))
                    .collect::<Vec<_>>();
                let signed_area = area(&polygon, work)?;
                if signed_area > 0.0 {
                    outer.push((wire, polygon, signed_area, Vec::new()));
                } else {
                    holes.push((wire, polygon));
                }
            }
            for (wire, polygon) in holes {
                let mut owner = None;
                for (i, (_, boundary, size, _)) in outer.iter().enumerate() {
                    if location(std::slice::from_ref(boundary), polygon[0], work)?
                        == Location::Inside
                        && owner.is_none_or(|(_, previous)| *size < previous)
                    {
                        owner = Some((i, *size));
                    }
                }
                let (index, _) = owner.ok_or_else(|| failed("orphan planar hole"))?;
                outer[index].3.push(wire);
            }
            for (wire, _, _, holes) in outer {
                let mut boundaries = Vec::new();
                for indices in std::iter::once(wire).chain(holes) {
                    let mut wire = Wire::new();
                    for i in 0..indices.len() {
                        let (a, b) = (indices[i], indices[(i + 1) % indices.len()]);
                        let edge = if let Some(edge) = curves.get(&key(a, b)) {
                            edge.clone()
                        } else {
                            let curve = C::from_straight_segment(
                                vertices.points[a].point(),
                                vertices.points[b].point(),
                            )
                            .ok_or_else(|| failed("curve type cannot construct exact lines"))?;
                            let edge =
                                Edge::try_new(&vertices.points[a], &vertices.points[b], curve)
                                    .map_err(|_| failed("invalid native arrangement edge"))?;
                            curves.insert(key(a, b), edge.clone());
                            edge
                        };
                        wire.push_back(if edge.front().id() == vertices.points[a].id() {
                            edge
                        } else {
                            edge.inverse()
                        });
                    }
                    boundaries.push(wire);
                }
                let source_boundaries = face.source.boundaries();
                let same = boundaries.len() == source_boundaries.len()
                    && boundaries.iter().zip(&source_boundaries).all(|(a, b)| {
                        a.len() == b.len()
                            && (0..a.len()).any(|shift| {
                                a.iter()
                                    .enumerate()
                                    .all(|(i, e)| e == &b[(i + shift) % b.len()])
                            })
                    });
                output[operation].push(if same {
                    face.source.clone()
                } else {
                    let orientation = face.source.orientation();
                    let boundaries = if orientation {
                        boundaries
                    } else {
                        boundaries.into_iter().map(|wire| wire.inverse()).collect()
                    };
                    let built = Face::try_new(boundaries, face.source.surface())
                        .map_err(|_| failed("invalid native planar face"))?;
                    if orientation { built } else { built.inverse() }
                });
            }
        }
    }
    Ok(output)
}

pub(super) fn partition<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    a: &Shell<Point3, C, S>,
    b: &Shell<Point3, C, S>,
    hints: Option<ShellOrientationHints>,
) -> ClassicResult<Option<[Shell<Point3, C, S>; 2]>> {
    if a.len() + b.len() > MAX_FACES {
        return Ok(None);
    }
    let mut work = 0;
    let Some(mut faces) = admit(a, 0, &mut work)? else {
        return Ok(None);
    };
    let Some(other) = admit(b, 1, &mut work)? else {
        return Ok(None);
    };
    faces.extend(other);
    let signs = if let Some(hints) = hints {
        [hints.first_inverted, hints.second_inverted]
    } else {
        [inverted(&faces, 0)?, inverted(&faces, 1)?]
    };
    for face in &faces {
        validate_loops(&face.loops, &mut work)?;
    }
    let mut result = Vec::new();
    for index in 0..faces.len() {
        for polygon in cells(&faces, index, &mut work)? {
            let face = &faces[index];
            let center = Point2::from_vec(
                polygon.iter().map(|p| p.to_vec()).sum::<Vector2>() / polygon.len() as f64,
            );
            let selected = selection(&faces, index, face.frame.lift(center), signs, &mut work)?;
            if selected.iter().any(|b| *b) {
                if result.len() >= MAX_CELLS {
                    return Err(failed("too many selected arrangement cells"));
                }
                result.push(Cell {
                    source: index,
                    selected,
                    points: polygon.into_iter().map(|p| face.frame.lift(p)).collect(),
                });
            }
        }
    }
    assemble(&faces, &result, &mut work).map(Some)
}
