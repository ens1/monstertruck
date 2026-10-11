//! Classic (0.3.2) boolean pipeline -- the self-contained default backend.
//!
//! This subtree grafts the proven upstream boolean assembly from the published
//! 0.3.2 crate. It is compiled only when the marching SSI is the active backend
//! (default features); a build with an external SSI backend never touches it
//! and runs the current `integrate` pipeline unchanged. The public entry points in
//! `transversal::integrate` cfg-dispatch `and`/`or`/`difference`/
//! `symmetric_difference` here.
//!
//! Graft-boundary adaptations vs. the 0.3.2 source (all documented inline):
//! - The Alternative intersection-curve arm is a raw
//!   `IntersectionCurve<PolylineCurve<Point3>, S, S>` (0.3.2 wrapped it with
//!   parameter polylines the pipeline never consumed).
//! - `altshell_to_shell` converts that arm back to `C` via
//!   `SurfaceCurve::with_boundaries(..).into()` (the conversion the current
//!   `ShapeOpsCurve` bound provides), not `IntersectionCurve::new(..).into()`.
//! - Complete material domains, including cavities and disconnected components,
//!   are classified together using oriented solid-angle winding.
//! - Failures identify their pipeline stage. A valid empty result is a solid
//!   with no boundary shells; ambiguity is never converted into empty geometry.

mod divide_face;
mod faces_classification;
mod intersection_curve;
mod loops_store;
mod planar;

use super::integrate::{ShapeOpsCurve, ShapeOpsError, ShapeOpsSurface, ShellOrientationHints};
use crate::alternative::Alternative;
use monstertruck_geometry::prelude::*;
use monstertruck_meshing::prelude::*;
use monstertruck_topology::*;
use std::iter;

type ClassicResult<T> = std::result::Result<T, ShapeOpsError>;

type AltCurve<C, S> = Alternative<C, IntersectionCurve<PolylineCurve<Point3>, S, S>>;
type AltCurveShell<C, S> = Shell<Point3, AltCurve<C, S>, S>;

// These limits bound the quadratic face-pair stage and the linear classification
// work separately. Tessellation remains a search aid; it never becomes output
// boundary geometry.
const MAX_INPUT_FACES: usize = 16_384;
const MAX_FACE_PAIRS: usize = 4_000_000;
const MAX_SEARCH_TRIANGLES: usize = 1_000_000;
const MAX_CLASSIFICATION_WORK: usize = 64_000_000;
const MAX_FACE_SAMPLES: usize = 12;

type SearchShell = Shell<Point3, PolylineCurve<Point3>, Option<PolygonMesh>>;

fn finite_point(point: Point3) -> bool {
    point.x.is_finite() && point.y.is_finite() && point.z.is_finite()
}

fn compensated_add(sum: &mut f64, correction: &mut f64, term: f64) {
    let corrected = term - *correction;
    let next = *sum + corrected;
    *correction = (next - *sum) - corrected;
    *sum = next;
}

struct SearchDomain {
    triangles: Vec<[Vector3; 3]>,
    origin: Point3,
    scale: f64,
    boundary_distance2: f64,
    inverted: bool,
}

impl SearchDomain {
    fn new(
        shell: &SearchShell,
        tol: f64,
        shell_index: usize,
        orientation_hint: Option<bool>,
    ) -> ClassicResult<Self> {
        let mut meshes = Vec::with_capacity(shell.len());
        let mut bounds = BoundingBox::<Point3>::default();
        let mut triangle_count = 0_usize;
        for (face_index, face) in shell.iter().enumerate() {
            let failed = || ShapeOpsError::ClassificationMeshFailed {
                shell_index,
                face_index,
            };
            let mesh = face.surface().ok_or_else(failed)?;
            if mesh.positions().is_empty() || mesh.positions().iter().any(|p| !finite_point(*p)) {
                return Err(failed());
            }
            for &point in mesh.positions() {
                bounds.push(point);
            }
            let face_triangles = mesh
                .face_iter()
                .map(|polygon| polygon.len().saturating_sub(2))
                .sum::<usize>();
            if face_triangles == 0 {
                return Err(failed());
            }
            triangle_count =
                triangle_count
                    .checked_add(face_triangles)
                    .ok_or(ShapeOpsError::ResourceLimit {
                        stage: "classification mesh triangles",
                        limit: MAX_SEARCH_TRIANGLES,
                    })?;
            if triangle_count > MAX_SEARCH_TRIANGLES {
                return Err(ShapeOpsError::ResourceLimit {
                    stage: "classification mesh triangles",
                    limit: MAX_SEARCH_TRIANGLES,
                });
            }
            meshes.push((mesh, face.orientation()));
        }
        // Scaling and translating before determinants avoids cancellation for
        // parts far from the origin and dependence on the caller's length unit.
        let origin = bounds.min();
        let diagonal = bounds.max() - origin;
        let scale = diagonal.x.max(diagonal.y).max(diagonal.z);
        if !scale.is_finite() || scale <= 0.0 {
            return Err(ShapeOpsError::ShellOrientationFailed { shell_index });
        }
        let mut triangles = Vec::with_capacity(triangle_count);
        let (mut volume, mut correction, mut absolute_volume) = (0.0, 0.0, 0.0);
        for (mesh, orientation) in meshes {
            for polygon in mesh.face_iter() {
                for index in 2..polygon.len() {
                    let mut triangle = [polygon[0], polygon[index - 1], polygon[index]]
                        .map(|vertex| (mesh.positions()[vertex.pos] - origin) / scale);
                    if !orientation {
                        triangle.swap(1, 2);
                    }
                    let term = triangle[0].dot(triangle[1].cross(triangle[2]));
                    absolute_volume += term.abs();
                    compensated_add(&mut volume, &mut correction, term);
                    triangles.push(triangle);
                }
            }
        }
        let inverted = match orientation_hint {
            Some(hint) => hint,
            None if volume.is_finite()
                && volume.abs() > absolute_volume * (64.0 * f64::EPSILON) =>
            {
                volume < 0.0
            }
            None => return Err(ShapeOpsError::ShellOrientationFailed { shell_index }),
        };
        // A point within the tessellation error of the other boundary cannot be
        // classified from that mesh. Reject it rather than choosing the side of
        // a chord instead of the side of the actual curved surface.
        let distance = (2.0 * tol / scale).max(128.0 * f64::EPSILON);
        Ok(Self {
            triangles,
            origin,
            scale,
            boundary_distance2: distance * distance,
            inverted,
        })
    }

    fn inside(&self, point: Point3, work: &mut usize) -> ClassicResult<Option<bool>> {
        *work = work
            .checked_add(self.triangles.len())
            .ok_or(ShapeOpsError::ResourceLimit {
                stage: "face classification triangle visits",
                limit: MAX_CLASSIFICATION_WORK,
            })?;
        if *work > MAX_CLASSIFICATION_WORK {
            return Err(ShapeOpsError::ResourceLimit {
                stage: "face classification triangle visits",
                limit: MAX_CLASSIFICATION_WORK,
            });
        }
        let point = (point - self.origin) / self.scale;
        let (mut angle, mut correction) = (0.0, 0.0);
        for triangle in &self.triangles {
            if triangle_distance2(point, *triangle) <= self.boundary_distance2 {
                return Ok(None);
            }
            let [a, b, c] = triangle.map(|vertex| (vertex - point).normalize());
            let numerator = a.dot(b.cross(c));
            let denominator = 1.0 + a.dot(b) + b.dot(c) + c.dot(a);
            let term = 2.0 * numerator.atan2(denominator);
            if !term.is_finite() {
                return Ok(None);
            }
            compensated_add(&mut angle, &mut correction, term);
        }
        // The signed solid-angle sum is independent of ray direction and has
        // no edge/vertex ray ambiguity. In a valid finite material domain it is
        // 0 outside/cavities and +1 inside; a complemented domain is 0 or -1.
        let winding = angle / (4.0 * std::f64::consts::PI);
        let integer = winding.round();
        if (winding - integer).abs() > 1.0e-5 {
            return Ok(None);
        }
        let valid = if self.inverted {
            integer == 0.0 || integer == -1.0
        } else {
            integer == 0.0 || integer == 1.0
        };
        Ok(valid.then_some(if self.inverted {
            integer == 0.0
        } else {
            integer == 1.0
        }))
    }
}

fn segment_distance2(point: Vector3, a: Vector3, b: Vector3) -> f64 {
    let edge = b - a;
    let length2 = edge.magnitude2();
    if length2 == 0.0 {
        (point - a).magnitude2()
    } else {
        let t = ((point - a).dot(edge) / length2).clamp(0.0, 1.0);
        (point - (a + edge * t)).magnitude2()
    }
}

fn triangle_distance2(point: Vector3, [a, b, c]: [Vector3; 3]) -> f64 {
    let ab = b - a;
    let ac = c - a;
    let normal = ab.cross(ac);
    let normal2 = normal.magnitude2();
    if normal2 > 0.0 {
        let ap = point - a;
        // Cross-product barycentrics avoid subtracting nearly equal squared
        // lengths for thin triangles.
        let u = ap.cross(ac).dot(normal) / normal2;
        let v = ab.cross(ap).dot(normal) / normal2;
        if u >= 0.0 && v >= 0.0 && u + v <= 1.0 {
            return ap.dot(normal).powi(2) / normal2;
        }
    }
    segment_distance2(point, a, b)
        .min(segment_distance2(point, b, c))
        .min(segment_distance2(point, c, a))
}

fn classify_face<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    face: &Face<Point3, AltCurve<C, S>, S>,
    other: &SearchDomain,
    tol: f64,
    work: &mut usize,
) -> ClassicResult<Option<bool>> {
    // Use interior points in the trimmed parameter domain, not a boundary
    // vertex (often shared with the other solid) or a chord-space barycenter.
    let single_face: Shell<_, _, _> = iter::once(face.clone()).collect();
    let triangulated = single_face.triangulation(tol);
    let Some(mesh) = triangulated.first().and_then(Face::surface) else {
        return Ok(None);
    };
    let mut candidates = Vec::with_capacity(MAX_FACE_SAMPLES + 1);
    for polygon in mesh.face_iter() {
        for index in 2..polygon.len() {
            let vertices = [polygon[0], polygon[index - 1], polygon[index]];
            let [a, b, c] = vertices.map(|vertex| mesh.positions()[vertex.pos]);
            let area = (b - a).cross(c - a).magnitude2();
            if !area.is_finite() || area <= 0.0 {
                continue;
            }
            let uv = vertices.map(|vertex| {
                vertex
                    .uv
                    .and_then(|index| mesh.uv_coords().as_slice().get(index))
                    .copied()
            });
            let [Some(u), Some(v), Some(w)] = uv else {
                continue;
            };
            let uv = (u + v + w) / 3.0;
            if !uv.x.is_finite() || !uv.y.is_finite() {
                continue;
            }
            candidates.push((area, uv));
            candidates.sort_unstable_by(|a, b| b.0.total_cmp(&a.0));
            candidates.truncate(MAX_FACE_SAMPLES);
        }
    }
    let surface = face.surface();
    let mut classification = None;
    for (_, uv) in candidates {
        let point = surface.subs(uv.x, uv.y);
        if !finite_point(point) {
            return Ok(None);
        }
        if let Some(current) = other.inside(point, work)? {
            if classification.is_some_and(|previous| previous != current) {
                return Ok(None);
            }
            classification = Some(current);
        }
    }
    Ok(classification)
}

/// Convert an Alternative-curve shell back to a target-curve shell, approximating
/// each intersection-curve arm with a polyline B-spline wrapped in a
/// `SurfaceCurve` (uses the `From<SurfaceCurve<BsplineCurve<Point3>, ..>>`
/// conversion the current `ShapeOpsCurve` bound guarantees).
fn altshell_to_shell<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    altshell: &AltCurveShell<C, S>,
) -> Option<Shell<Point3, C, S>> {
    altshell.try_mapped(
        |p| Some(*p),
        |c| match c {
            Alternative::FirstType(c) => Some(c.clone()),
            Alternative::SecondType(ic) => {
                let surface0 = ic.surface0().clone();
                let surface1 = ic.surface1().clone();
                let mut points = ic.leader().as_slice().to_vec();
                points.dedup_by(|lhs, rhs| (*lhs).near(&*rhs));
                if points.len() < 2 {
                    return None;
                }
                let denominator = (points.len() - 1) as f64;
                let knot_vec = KnotVector::from(
                    iter::once(0.0)
                        .chain((0..points.len()).map(|index| index as f64 / denominator))
                        .chain(iter::once(1.0))
                        .collect::<Vec<_>>(),
                );
                let bspline = BsplineCurve::new(knot_vec, points);
                let boundary0: Option<ParameterCurve<BoundaryCurve2D, S>> = None;
                let boundary1: Option<ParameterCurve<BoundaryCurve2D, S>> = None;
                Some(
                    SurfaceCurve::with_boundaries(
                        surface0, surface1, bspline, boundary0, boundary1,
                    )
                    .into(),
                )
            }
        },
        |s| Some(s.clone()),
    )
}

/// Split two complete material boundaries into intersection and union faces.
/// Inner cavity shells and disconnected components must be classified together.
fn process_one_pair_of_shells<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    shell0: &Shell<Point3, C, S>,
    shell1: &Shell<Point3, C, S>,
    tol: f64,
    orientation_hints: Option<ShellOrientationHints>,
) -> ClassicResult<[Shell<Point3, C, S>; 2]> {
    if let Some(result) = planar::partition(shell0, shell1, orientation_hints)? {
        return Ok(result);
    }
    let pairs = shell0.len().checked_mul(shell1.len());
    if pairs.is_none_or(|pairs| pairs > MAX_FACE_PAIRS) {
        return Err(ShapeOpsError::ResourceLimit {
            stage: "surface intersection face pairs",
            limit: MAX_FACE_PAIRS,
        });
    }
    let poly_shell0 = shell0.triangulation(tol);
    let poly_shell1 = shell1.triangulation(tol);
    let domain0 = SearchDomain::new(
        &poly_shell0,
        tol,
        0,
        orientation_hints.map(|hints| hints.first_inverted),
    )?;
    let domain1 = SearchDomain::new(
        &poly_shell1,
        tol,
        1,
        orientation_hints.map(|hints| hints.second_inverted),
    )?;
    let altshell0: AltCurveShell<C, S> =
        shell0.mapped(|x| *x, |c| Alternative::FirstType(c.clone()), Clone::clone);
    let altshell1: AltCurveShell<C, S> =
        shell1.mapped(|x| *x, |c| Alternative::FirstType(c.clone()), Clone::clone);
    let loops_store::LoopsStoreQuadruple {
        geom_loops_store0: loops_store0,
        geom_loops_store1: loops_store1,
        ..
    } = loops_store::create_loops_stores(&altshell0, &poly_shell0, &altshell1, &poly_shell1, tol)
        .ok_or(ShapeOpsError::IntersectionLoopsFailed)?;
    let mut cls0 = divide_face::divide_faces(&altshell0, &loops_store0, tol)
        .ok_or(ShapeOpsError::DivideFacesFailed { shell_index: 0 })?;
    cls0.integrate_by_component()
        .ok_or(ShapeOpsError::UnknownClassificationFailed { shell_index: 0 })?;
    let mut cls1 = divide_face::divide_faces(&altshell1, &loops_store1, tol)
        .ok_or(ShapeOpsError::DivideFacesFailed { shell_index: 1 })?;
    cls1.integrate_by_component()
        .ok_or(ShapeOpsError::UnknownClassificationFailed { shell_index: 1 })?;
    let [mut and0, mut or0, unknown0] = cls0.and_or_unknown();
    let [mut and1, mut or1, unknown1] = cls1.and_or_unknown();
    let mut work = 0;
    for (shell_index, unknown, domain, and_faces, or_faces) in [
        (0, unknown0, &domain1, &mut and0, &mut or0),
        (1, unknown1, &domain0, &mut and1, &mut or1),
    ] {
        for (face_index, face) in unknown.into_iter().enumerate() {
            let inside = classify_face(&face, domain, tol, &mut work)?.ok_or(
                ShapeOpsError::AmbiguousFaceClassification {
                    shell_index,
                    face_index,
                },
            )?;
            if inside {
                and_faces.push(face);
            } else {
                or_faces.push(face);
            }
        }
    }
    and0.append(&mut and1);
    or0.append(&mut or1);
    Ok([
        altshell_to_shell(&and0)
            .ok_or(ShapeOpsError::AltShellConversionFailed { operation: "and" })?,
        altshell_to_shell(&or0)
            .ok_or(ShapeOpsError::AltShellConversionFailed { operation: "or" })?,
    ])
}

fn validate_inputs<C, S>(
    solid0: &Solid<Point3, C, S>,
    solid1: &Solid<Point3, C, S>,
    tol: f64,
) -> ClassicResult<()> {
    if !tol.is_finite() || tol < TOLERANCE {
        return Err(ShapeOpsError::InvalidTolerance);
    }
    for (operand_index, solid) in [solid0, solid1].into_iter().enumerate() {
        let faces = solid
            .boundaries()
            .iter()
            .map(|shell| shell.len())
            .sum::<usize>();
        if faces > MAX_INPUT_FACES {
            return Err(ShapeOpsError::ResourceLimit {
                stage: "input faces",
                limit: MAX_INPUT_FACES,
            });
        }
        for (shell_index, shell) in solid.boundaries().iter().enumerate() {
            shell
                .check_solid_boundary()
                .map_err(|source| ShapeOpsError::InvalidInputShell {
                    operand_index,
                    shell_index,
                    source,
                })?;
        }
    }
    Ok(())
}

// This shortcut requires identical native boundary identities and trims; no
// approximate geometric equality or shared carrier alone establishes equality.
fn identical_boundaries<C, S>(
    solid0: &Solid<Point3, C, S>,
    solid1: &Solid<Point3, C, S>,
) -> Option<bool> {
    if solid0.boundaries().len() != solid1.boundaries().len() {
        return None;
    }
    let mut reversed = None;
    for (shell0, shell1) in solid0.boundaries().iter().zip(solid1.boundaries()) {
        if shell0.len() != shell1.len() {
            return None;
        }
        for (face0, face1) in shell0.iter().zip(shell1.iter()) {
            if face0.id() != face1.id()
                || face0.absolute_boundaries() != face1.absolute_boundaries()
            {
                return None;
            }
            let current = face0.orientation() != face1.orientation();
            if reversed.is_some_and(|previous| previous != current) {
                return None;
            }
            reversed = Some(current);
        }
    }
    reversed
}

fn finalize<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    operation: &'static str,
    shell: Shell<Point3, C, S>,
) -> ClassicResult<Solid<Point3, C, S>> {
    let boundaries = shell.connected_components();
    Solid::try_new(boundaries)
        .map_err(|source| ShapeOpsError::InvalidOutputShell { operation, source })
}

/// AND operation between two solids (classic backend).
pub(crate) fn and<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    solid0: &Solid<Point3, C, S>,
    solid1: &Solid<Point3, C, S>,
    tol: f64,
) -> ClassicResult<Solid<Point3, C, S>> {
    and_with_orientation_hints(solid0, solid1, None, tol)
}

pub(crate) fn and_with_orientation_hints<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    solid0: &Solid<Point3, C, S>,
    solid1: &Solid<Point3, C, S>,
    hints: Option<ShellOrientationHints>,
    tol: f64,
) -> ClassicResult<Solid<Point3, C, S>> {
    let operation = "and";
    validate_inputs(solid0, solid1, tol)?;
    if solid0.is_empty() || solid1.is_empty() {
        return Ok(Solid::empty());
    }
    if let Some(reversed) = identical_boundaries(solid0, solid1) {
        return Ok(if reversed {
            Solid::empty()
        } else {
            solid0.clone()
        });
    }
    // Cavities are part of the same material domain. Iteratively intersecting
    // outer stock with each cavity boundary inverted the intended semantics.
    let shell0: Shell<_, _, _> = solid0.face_iter().cloned().collect();
    let shell1: Shell<_, _, _> = solid1.face_iter().cloned().collect();
    let [and_shell, _] = process_one_pair_of_shells(&shell0, &shell1, tol, hints)?;
    finalize(operation, and_shell)
}

/// OR operation between two solids (classic backend).
pub(crate) fn or<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    solid0: &Solid<Point3, C, S>,
    solid1: &Solid<Point3, C, S>,
    tol: f64,
) -> ClassicResult<Solid<Point3, C, S>> {
    let operation = "or";
    validate_inputs(solid0, solid1, tol)?;
    if solid0.is_empty() {
        return Ok(solid1.clone());
    }
    if solid1.is_empty() {
        return Ok(solid0.clone());
    }
    if let Some(reversed) = identical_boundaries(solid0, solid1) {
        return if reversed {
            Err(ShapeOpsError::UnboundedResult { operation })
        } else {
            Ok(solid0.clone())
        };
    }
    let shell0: Shell<_, _, _> = solid0.face_iter().cloned().collect();
    let shell1: Shell<_, _, _> = solid1.face_iter().cloned().collect();
    let [_, or_shell] = process_one_pair_of_shells(&shell0, &shell1, tol, None)?;
    // The union of two nonempty material domains cannot be empty. With
    // complemented operands, cancellation of every boundary is universal.
    if or_shell.is_empty() {
        return Err(ShapeOpsError::UnboundedResult { operation });
    }
    finalize(operation, or_shell)
}

/// Difference: the region inside `solid0` but outside `solid1` (classic backend).
pub(crate) fn difference<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    solid0: &Solid<Point3, C, S>,
    solid1: &Solid<Point3, C, S>,
    tol: f64,
) -> ClassicResult<Solid<Point3, C, S>> {
    validate_inputs(solid0, solid1, tol)?;
    if solid0.is_empty() || identical_boundaries(solid0, solid1) == Some(false) {
        return Ok(Solid::empty());
    }
    if solid1.is_empty() {
        return Ok(solid0.clone());
    }
    let mut neg = solid1.clone();
    neg.not();
    and(solid0, &neg, tol)
}

/// Symmetric difference (XOR): the region inside exactly one solid (classic backend).
pub(crate) fn symmetric_difference<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    solid0: &Solid<Point3, C, S>,
    solid1: &Solid<Point3, C, S>,
    tol: f64,
) -> ClassicResult<Solid<Point3, C, S>> {
    let d0 = difference(solid0, solid1, tol)?;
    let d1 = difference(solid1, solid0, tol)?;
    or(&d0, &d1, tol)
}
