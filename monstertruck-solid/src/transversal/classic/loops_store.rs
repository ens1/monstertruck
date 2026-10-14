//! Classic (0.3.2) boolean loops store.
//!
//! Ported verbatim from the published 0.3.2 crate's `transversal::loops_store`,
//! adjusted only for the graft: the marching intersection-curve backend now
//! retains an implicit native surface curve with paired UV continuation hints
//! and producing face identities. Its geometric boundary is never the mesh
//! polyline or an interpolated UV trace.

#![allow(clippy::many_single_char_names)]

use super::intersection_curve;
use monstertruck_core::cgmath64::*;
use monstertruck_geometry::prelude::*;
use monstertruck_meshing::prelude::*;
use monstertruck_topology::{Vertex, *};
use rustc_hash::FxHashMap as HashMap;

type PolylineCurve = monstertruck_meshing::prelude::PolylineCurve<Point3>;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum ShapesOpStatus {
    Unknown,
    And,
    Or,
}

impl ShapesOpStatus {
    fn not(self) -> Self {
        match self {
            Self::Unknown => Self::Unknown,
            Self::And => Self::Or,
            Self::Or => Self::And,
        }
    }
}

#[derive(Clone, Debug)]
pub(super) struct BoundaryWire<P, C> {
    wire: Wire<P, C>,
    status: ShapesOpStatus,
}

impl<P, C> BoundaryWire<P, C> {
    #[inline(always)]
    pub(super) fn new(wire: Wire<P, C>, status: ShapesOpStatus) -> Self { Self { wire, status } }
    #[inline(always)]
    pub(super) fn status(&self) -> ShapesOpStatus { self.status }
    #[inline(always)]
    pub(super) fn inverse(&self) -> Self {
        Self {
            wire: self.wire.inverse(),
            status: self.status.not(),
        }
    }
}

impl ShapesOpStatus {
    fn from_is_curve<C, S0, S1, T0: Clone, T1: Clone>(
        curve: &SurfaceCurve<C, S0, S1, T0, T1>,
    ) -> Option<ShapesOpStatus>
    where
        C: ParametricCurve3D + BoundedCurve,
        S0: ParametricSurface3D + SearchNearestParameter<SurfaceParameter, Point = Point3>,
        S1: ParametricSurface3D + SearchNearestParameter<SurfaceParameter, Point = Point3>, {
        let (t0, t1) = curve.range_tuple();
        let mut resolved = None;
        for fraction in [0.25, 0.5, 0.75] {
            let t = t0 + (t1 - t0) * fraction;
            let (_, pt0, pt1) = curve.search_triple(t, 100)?;
            let der = curve.leader().derivative(t);
            let normal0 = curve.surface0().normal(pt0.x, pt0.y);
            let normal1 = curve.surface1().normal(pt1.x, pt1.y);
            let denominator = normal0.magnitude() * der.magnitude() * normal1.magnitude();
            if !denominator.is_finite() || denominator <= 0.0 {
                return None;
            }
            let direction = normal0.cross(der).dot(normal1) / denominator;
            // Tangent/coincident contacts have no transversal inside/outside
            // direction; assigning them And by default invents material.
            if !direction.is_finite() {
                return None;
            }
            if direction.abs() <= 1.0e-10 {
                continue;
            }
            let status = if direction > 0.0 {
                ShapesOpStatus::Or
            } else {
                ShapesOpStatus::And
            };
            if resolved.is_some_and(|old| old != status) {
                return None;
            }
            resolved = Some(status);
        }
        resolved
    }
}

impl<P, C> std::ops::Deref for BoundaryWire<P, C> {
    type Target = Wire<P, C>;
    #[inline(always)]
    fn deref(&self) -> &Self::Target { &self.wire }
}

impl<P, C> std::ops::DerefMut for BoundaryWire<P, C> {
    #[inline(always)]
    fn deref_mut(&mut self) -> &mut Self::Target { &mut self.wire }
}

#[derive(Clone, Debug)]
pub(super) struct Loops<P, C>(Vec<BoundaryWire<P, C>>);
#[derive(Clone, Debug)]
pub(super) struct LoopsStore<P, C>(Vec<Loops<P, C>>);

impl<P, C> std::ops::Deref for Loops<P, C> {
    type Target = Vec<BoundaryWire<P, C>>;
    #[inline(always)]
    fn deref(&self) -> &Self::Target { &self.0 }
}

impl<P, C> std::ops::DerefMut for Loops<P, C> {
    #[inline(always)]
    fn deref_mut(&mut self) -> &mut Self::Target { &mut self.0 }
}

impl<P, C> std::ops::Deref for LoopsStore<P, C> {
    type Target = Vec<Loops<P, C>>;
    #[inline(always)]
    fn deref(&self) -> &Self::Target { &self.0 }
}

impl<P, C> std::ops::DerefMut for LoopsStore<P, C> {
    #[inline(always)]
    fn deref_mut(&mut self) -> &mut Self::Target { &mut self.0 }
}

impl<P, C> FromIterator<BoundaryWire<P, C>> for Loops<P, C> {
    #[inline(always)]
    fn from_iter<I: IntoIterator<Item = BoundaryWire<P, C>>>(iter: I) -> Self {
        Self(Vec::from_iter(iter))
    }
}

impl<'a, P, C, S> From<&'a Face<P, C, S>> for Loops<P, C> {
    #[inline(always)]
    fn from(face: &'a Face<P, C, S>) -> Loops<P, C> {
        face.absolute_boundaries()
            .iter()
            .map(|wire| BoundaryWire::new(wire.clone(), ShapesOpStatus::Unknown))
            .collect()
    }
}

impl<'a, P: 'a, C: 'a, S: 'a> FromIterator<&'a Face<P, C, S>> for LoopsStore<P, C> {
    fn from_iter<I: IntoIterator<Item = &'a Face<P, C, S>>>(iter: I) -> Self {
        Self(iter.into_iter().map(Loops::from).collect())
    }
}

impl<'a, P, C> IntoIterator for &'a LoopsStore<P, C> {
    type Item = <&'a Vec<Loops<P, C>> as IntoIterator>::Item;
    type IntoIter = <&'a Vec<Loops<P, C>> as IntoIterator>::IntoIter;
    fn into_iter(self) -> Self::IntoIter { self.0.iter() }
}

#[derive(Clone, Debug, Copy, PartialEq)]
enum ParameterKind {
    Front,
    Back,
    Inner(f64),
}

impl<P: Copy, C: Clone> Loops<P, C> {
    fn search_parameter(&self, pt: P) -> Option<(usize, usize, ParameterKind)>
    where
        P: Tolerance,
        C: BoundedCurve<Point = P> + SearchParameter<CurveParameter, Point = P>, {
        self.iter()
            .enumerate()
            .flat_map(move |(i, wire)| wire.iter().enumerate().map(move |(j, edge)| (i, j, edge)))
            .find_map(|(i, j, edge)| {
                let curve = edge.curve();
                curve.search_parameter(pt, None, 100).and_then(|t| {
                    let (start, end) = curve.range_tuple();
                    let kind = if pt.near(&curve.front()) {
                        ParameterKind::Front
                    } else if pt.near(&curve.back()) {
                        ParameterKind::Back
                    } else if t.is_finite() && start < t && t < end {
                        ParameterKind::Inner(t)
                    } else {
                        return None;
                    };
                    Some((i, j, kind))
                })
            })
    }

    fn change_vertex(
        &mut self,
        old_vertex: &Vertex<P>,
        new_vertex: &Vertex<P>,
        emap: &mut HashMap<EdgeId<C>, Edge<P, C>>,
    ) {
        self.iter_mut()
            .flat_map(|wire| wire.iter_mut())
            .for_each(|edge| {
                let mut new_edge = if edge.absolute_front() == old_vertex {
                    emap.entry(edge.id()).or_insert_with(|| {
                        Edge::new(new_vertex, edge.absolute_back(), edge.curve())
                    })
                } else if edge.absolute_back() == old_vertex {
                    emap.entry(edge.id()).or_insert_with(|| {
                        Edge::new(edge.absolute_front(), new_vertex, edge.curve())
                    })
                } else {
                    return;
                }
                .clone();
                if !edge.orientation() {
                    new_edge.invert();
                }
                // Remove the edge from the HashMap when it is no longer there because Id reassignment will occur.
                if edge.count() == 1 {
                    emap.remove(&edge.id());
                }
                *edge = new_edge;
            })
    }

    fn swap_edge_into_wire(&mut self, edge_id: EdgeId<C>, new_wire: &Wire<P, C>) {
        self.iter_mut().for_each(|wire| {
            let mut iter = wire.iter().enumerate();
            if let Some((idx, edge)) = iter.find(|(_, edge)| edge.id() == edge_id) {
                if edge.orientation() {
                    wire.swap_edge_into_wire(idx, new_wire.clone());
                } else {
                    wire.swap_edge_into_wire(idx, new_wire.inverse());
                }
            }
        });
    }

    #[inline(always)]
    fn add_independent_loop(&mut self, r#loop: BoundaryWire<P, C>) {
        self.push(r#loop.inverse());
        self.push(r#loop);
    }

    fn add_edge(
        &mut self,
        edge0: Edge<P, C>,
        status: ShapesOpStatus,
    ) -> [Option<(usize, usize)>; 2] {
        let a = self.iter().enumerate().find_map(|(i, wire)| {
            wire.iter().enumerate().find_map(|(j, edge)| {
                if edge.front() == edge0.back() {
                    Some((i, j))
                } else {
                    None
                }
            })
        });
        let b = self.iter().enumerate().find_map(|(i, wire)| {
            wire.iter().enumerate().find_map(|(j, edge)| {
                if edge.front() == edge0.front() {
                    Some((i, j))
                } else {
                    None
                }
            })
        });
        if let Some((wire_index0, edge_index0)) = a {
            self[wire_index0].rotate_left(edge_index0);
            self[wire_index0].push_front(edge0.clone());
            self[wire_index0].push_back(edge0.inverse());
        }
        match (a, b) {
            (Some((wire_index0, edge_index0)), Some((wire_index1, edge_index1))) => {
                if wire_index0 == wire_index1 {
                    let len = self[wire_index0].len() - 2;
                    let edge_index1 = (len + edge_index1 - edge_index0) % len + 1;
                    let new_wire = self[wire_index0].split_off(edge_index1);
                    self[wire_index0].status = status;
                    self.push(BoundaryWire::new(new_wire, status.not()));
                } else {
                    let mut new_wire0 = self[wire_index1].clone();
                    let mut new_wire1 = new_wire0.split_off(edge_index1);
                    new_wire0.append(&mut self[wire_index0]);
                    new_wire0.append(&mut new_wire1);
                    self[wire_index0] = new_wire0;
                    self.swap_remove(wire_index1);
                }
            }
            (None, Some((wire_index1, edge_index1))) => {
                self[wire_index1].rotate_left(edge_index1);
                self[wire_index1].push_front(edge0.inverse());
                self[wire_index1].push_back(edge0);
            }
            (None, None) => self.push(BoundaryWire::new(
                vec![edge0.inverse(), edge0].into(),
                ShapesOpStatus::Unknown,
            )),
            _ => {}
        }
        [a, b]
    }
}

impl<P: Copy + Tolerance, C: Clone> LoopsStore<P, C> {
    #[inline(always)]
    fn change_vertex(
        &mut self,
        old_vertex: &Vertex<P>,
        new_vertex: &Vertex<P>,
        emap: &mut HashMap<EdgeId<C>, Edge<P, C>>,
    ) {
        self.iter_mut()
            .for_each(|loops| loops.change_vertex(old_vertex, new_vertex, emap));
    }

    #[inline(always)]
    fn swap_edge_into_wire(&mut self, edge_id: EdgeId<C>, new_wire: &Wire<P, C>) {
        self.iter_mut()
            .for_each(|loops| loops.swap_edge_into_wire(edge_id, new_wire))
    }

    fn add_polygon_vertex(
        &mut self,
        loops_index: usize,
        v: &Vertex<P>,
        emap: &mut HashMap<EdgeId<C>, Edge<P, C>>,
    ) -> Option<(usize, usize, ParameterKind)>
    where
        C: Cut<Point = P> + SearchParameter<CurveParameter, Point = P>,
    {
        let pt = v.point();
        let (wire_index, edge_index, kind) = self[loops_index].search_parameter(pt)?;
        match kind {
            ParameterKind::Front => {
                let old_vertex = self[loops_index][wire_index][edge_index]
                    .absolute_front()
                    .clone();
                self.change_vertex(&old_vertex, v, emap);
            }
            ParameterKind::Back => {
                let old_vertex = self[loops_index][wire_index][edge_index]
                    .absolute_back()
                    .clone();
                self.change_vertex(&old_vertex, v, emap);
            }
            ParameterKind::Inner(t) => {
                let edge = self[loops_index][wire_index][edge_index].absolute_clone();
                let edge_id = edge.id();
                let (edge0, edge1) = edge.cut_with_parameter(v, t)?;
                let new_wire: Wire<_, _> = vec![edge0, edge1].into();
                self.swap_edge_into_wire(edge_id, &new_wire);
            }
        }
        Some((wire_index, edge_index, kind))
    }
}

impl<C> LoopsStore<Point3, C> {
    fn add_geom_vertex<S>(
        &mut self,
        (loops_index, wire_index, edge_index): (usize, usize, usize),
        v: &Vertex<Point3>,
        kind: ParameterKind,
        another_surface: &S,
        emap: &mut HashMap<EdgeId<C>, Edge<Point3, C>>,
        tol: f64,
    ) -> Option<()>
    where
        C: Cut<Point = Point3, Vector = Vector3>
            + SearchNearestParameter<CurveParameter, Point = Point3>,
        S: ParametricSurface3D + SearchNearestParameter<SurfaceParameter, Point = Point3>,
    {
        let seed = v.point();
        match kind {
            ParameterKind::Front => {
                let old_vertex = self[loops_index][wire_index][edge_index]
                    .absolute_front()
                    .clone();
                v.set_point(old_vertex.point());
                self.change_vertex(&old_vertex, v, emap);
            }
            ParameterKind::Back => {
                let old_vertex = self[loops_index][wire_index][edge_index]
                    .absolute_back()
                    .clone();
                v.set_point(old_vertex.point());
                self.change_vertex(&old_vertex, v, emap);
            }
            ParameterKind::Inner(_) => {
                let curve = self[loops_index][wire_index][edge_index].curve();
                let (pt, t, _) =
                    curve_surface_projection(&curve, None, another_surface, None, v.point(), 100)?;
                v.set_point(pt);
                let edge = self[loops_index][wire_index][edge_index].absolute_clone();
                let edge_id = edge.id();
                let (edge0, edge1) = edge.cut_with_parameter(v, t)?;
                let new_wire: Wire<_, _> = vec![edge0, edge1].into();
                self.swap_edge_into_wire(edge_id, &new_wire);
            }
        }
        let point = v.point();
        if !super::finite_point(point) || point.distance(seed) > 4.0 * tol {
            return None;
        }
        let uv = another_surface.search_nearest_parameter(point, None, 100)?;
        let on_surface = another_surface.subs(uv.0, uv.1);
        let residual = TOLERANCE;
        (super::finite_point(on_surface) && point.distance(on_surface) <= residual).then_some(())
    }
}

fn curve_surface_projection<C, S>(
    curve: &C,
    curve_hint: Option<f64>,
    surface: &S,
    surface_hint: Option<(f64, f64)>,
    point: Point3,
    trials: usize,
) -> Option<(Point3, f64, Point2)>
where
    C: ParametricCurve3D + BoundedCurve + SearchNearestParameter<CurveParameter, Point = Point3>,
    S: ParametricSurface3D + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    let (start, end) = curve.range_tuple();
    if !start.is_finite() || !end.is_finite() || end <= start || !super::finite_point(point) {
        return None;
    }
    let mut t = curve
        .search_nearest_parameter(point, curve_hint, 100)?
        .clamp(start, end);
    let mut hint = surface_hint;
    for _ in 0..trials.min(100) {
        let on_curve = curve.subs(t);
        if !super::finite_point(on_curve) {
            return None;
        }
        let uv = surface.search_nearest_parameter(on_curve, hint, 100)?;
        if !uv.0.is_finite() || !uv.1.is_finite() {
            return None;
        }
        let on_surface = surface.subs(uv.0, uv.1);
        let distance = on_surface.distance(on_curve);
        let tolerance = TOLERANCE;
        if distance.is_finite() && distance <= tolerance {
            return Some((on_curve, t, uv.into()));
        }
        if !distance.is_finite() {
            return None;
        }
        let tangent = curve.derivative(t);
        let normal = surface.normal(uv.0, uv.1);
        let denominator = tangent.dot(normal);
        let scale = tangent.magnitude() * normal.magnitude();
        if !denominator.is_finite()
            || !scale.is_finite()
            || scale <= 0.0
            || denominator.abs() <= scale * (64.0 * f64::EPSILON)
        {
            return None;
        }
        let mut step = ((on_surface - on_curve).dot(normal) / denominator)
            .clamp(-(end - start) * 0.25, (end - start) * 0.25);
        let mut accepted = None;
        for _ in 0..12 {
            let next = (t + step).clamp(start, end);
            if next == t {
                break;
            }
            let next_point = curve.subs(next);
            if super::finite_point(next_point)
                && let Some(next_uv) = surface.search_nearest_parameter(next_point, Some(uv), 100)
            {
                let next_surface = surface.subs(next_uv.0, next_uv.1);
                let next_distance = next_point.distance(next_surface);
                if next_distance.is_finite() && next_distance < distance {
                    accepted = Some((next, next_uv));
                    break;
                }
            }
            step *= 0.5;
        }
        let (next, next_hint) = accepted?;
        t = next;
        hint = Some(next_hint);
    }
    None
}

fn create_independent_loop<P, C, D>(mut poly_curve0: C) -> Wire<P, D>
where
    C: Cut<Point = P>,
    D: From<C>, {
    let (t0, t1) = poly_curve0.range_tuple();
    let t = (t0 + t1) / 2.0;
    let poly_curve1 = poly_curve0.cut(t);
    let v0 = Vertex::new(poly_curve0.front());
    let v1 = Vertex::new(poly_curve1.front());
    let edge0 = Edge::new(&v0, &v1, poly_curve0.into());
    let edge1 = Edge::new(&v1, &v0, poly_curve1.into());
    wire![edge0, edge1]
}

#[allow(dead_code)]
pub(super) struct LoopsStoreQuadruple<C> {
    pub(super) geom_loops_store0: LoopsStore<Point3, C>,
    pub(super) poly_loops_store0: LoopsStore<Point3, PolylineCurve>,
    pub(super) geom_loops_store1: LoopsStore<Point3, C>,
    pub(super) poly_loops_store1: LoopsStore<Point3, PolylineCurve>,
}

pub(super) fn create_loops_stores<C, S>(
    geom_shell0: &Shell<Point3, C, S>,
    poly_shell0: &Shell<Point3, PolylineCurve, Option<PolygonMesh>>,
    geom_shell1: &Shell<Point3, C, S>,
    poly_shell1: &Shell<Point3, PolylineCurve, Option<PolygonMesh>>,
    tol: f64,
) -> super::ClassicResult<LoopsStoreQuadruple<C>>
where
    C: SearchNearestParameter<CurveParameter, Point = Point3>
        + SearchParameter<CurveParameter, Point = Point3>
        + Cut<Point = Point3, Vector = Vector3>
        + From<super::intersection_branch::IntersectionBranch<S>>,
    S: super::ShapeOpsSurface,
{
    let mut geom_loops_store0: LoopsStore<_, _> = geom_shell0.face_iter().collect();
    let mut poly_loops_store0: LoopsStore<_, _> = poly_shell0.face_iter().collect();
    let mut geom_loops_store1: LoopsStore<_, _> = geom_shell1.face_iter().collect();
    let mut poly_loops_store1: LoopsStore<_, _> = poly_shell1.face_iter().collect();
    let store0_len = geom_loops_store0.len();
    let store1_len = geom_loops_store1.len();
    let mut intersection_work = 0;
    let mut pair = (0, 0);
    let mut stage = "classification mesh admission";
    (0..store0_len)
        .flat_map(move |i| (0..store1_len).map(move |j| (i, j)))
        .try_for_each(|(face_index0, face_index1)| {
            pair = (face_index0, face_index1);
            stage = "classification mesh admission";
            let ori0 = geom_shell0[face_index0].orientation();
            let ori1 = geom_shell1[face_index1].orientation();
            let surface0 = geom_shell0[face_index0].surface();
            let surface1 = geom_shell1[face_index1].surface();
            let polygon0 = poly_shell0[face_index0].surface()?;
            let polygon1 = poly_shell1[face_index1].surface()?;
            stage = "bounded intersection search and paired UV continuation";
            intersection_curve::intersection_curves(
                surface0.clone(),
                &polygon0,
                surface1.clone(),
                &polygon1,
                tol,
                &mut intersection_work,
            )?
            .into_iter()
            .try_for_each(|(polyline, mut intersection_curve)| {
                stage = "transversal branch orientation";
                intersection_curve
                    .bind([geom_shell0[face_index0].id(), geom_shell1[face_index1].id()]);
                let status = ShapesOpStatus::from_is_curve(&intersection_curve.curve)?;
                let (status0, status1) = match (ori0, ori1) {
                    (true, true) => (status, status.not()),
                    (true, false) => (status.not(), status.not()),
                    (false, true) => (status, status),
                    (false, false) => (status.not(), status),
                };
                stage = "shared native boundary endpoints and intersection-loop assembly";
                if polyline.front().near(&polyline.back()) {
                    let poly_wire = create_independent_loop(polyline);
                    poly_loops_store0[face_index0]
                        .add_independent_loop(BoundaryWire::new(poly_wire.clone(), status0));
                    poly_loops_store1[face_index1]
                        .add_independent_loop(BoundaryWire::new(poly_wire, status1));
                    stage = "whole-span continuation admission";
                    match intersection_curve.admit_continuation(TOLERANCE) {
                        Ok(()) | Err(super::super::ContinuationError::UnsupportedCarrier) => {}
                        Err(_) => return None,
                    }
                    let geom_wire = create_independent_loop(intersection_curve);
                    geom_loops_store0[face_index0]
                        .add_independent_loop(BoundaryWire::new(geom_wire.clone(), status0));
                    geom_loops_store1[face_index1]
                        .add_independent_loop(BoundaryWire::new(geom_wire, status1));
                } else {
                    let pv0 = Vertex::new(polyline.front());
                    let pv1 = Vertex::new(polyline.back());
                    let gv0 = Vertex::new(intersection_curve.front());
                    let gv1 = Vertex::new(intersection_curve.back());
                    let mut pemap0 = HashMap::default();
                    let mut pemap1 = HashMap::default();
                    let mut gemap0 = HashMap::default();
                    let mut gemap1 = HashMap::default();
                    let idx00 =
                        poly_loops_store0.add_polygon_vertex(face_index0, &pv0, &mut pemap0);
                    if let Some((wire_index, edge_index, kind)) = idx00 {
                        geom_loops_store0.add_geom_vertex(
                            (face_index0, wire_index, edge_index),
                            &gv0,
                            kind,
                            &surface1,
                            &mut gemap0,
                            tol,
                        )?;
                        intersection_curve.snap_endpoint(true, gv0.point())?;
                    }
                    let idx01 =
                        poly_loops_store0.add_polygon_vertex(face_index0, &pv1, &mut pemap1);
                    if let Some((wire_index, edge_index, kind)) = idx01 {
                        geom_loops_store0.add_geom_vertex(
                            (face_index0, wire_index, edge_index),
                            &gv1,
                            kind,
                            &surface1,
                            &mut gemap1,
                            tol,
                        )?;
                        intersection_curve.snap_endpoint(false, gv1.point())?;
                    }
                    let idx10 =
                        poly_loops_store1.add_polygon_vertex(face_index1, &pv0, &mut pemap0);
                    if let Some((wire_index, edge_index, kind)) = idx10 {
                        geom_loops_store1.add_geom_vertex(
                            (face_index1, wire_index, edge_index),
                            &gv0,
                            kind,
                            &surface0,
                            &mut gemap0,
                            tol,
                        )?;
                        intersection_curve.snap_endpoint(true, gv0.point())?;
                    }
                    let idx11 =
                        poly_loops_store1.add_polygon_vertex(face_index1, &pv1, &mut pemap1);
                    if let Some((wire_index, edge_index, kind)) = idx11 {
                        geom_loops_store1.add_geom_vertex(
                            (face_index1, wire_index, edge_index),
                            &gv1,
                            kind,
                            &surface0,
                            &mut gemap1,
                            tol,
                        )?;
                        intersection_curve.snap_endpoint(false, gv1.point())?;
                    }
                    stage = "whole-span continuation admission";
                    match intersection_curve.admit_continuation(TOLERANCE) {
                        Ok(()) | Err(super::super::ContinuationError::UnsupportedCarrier) => {}
                        Err(_) => return None,
                    }
                    let pedge = Edge::new(&pv0, &pv1, polyline);
                    let gedge = Edge::new(&gv0, &gv1, intersection_curve.into());
                    poly_loops_store0[face_index0].add_edge(pedge.clone(), status0);
                    geom_loops_store0[face_index0].add_edge(gedge.clone(), status0);
                    poly_loops_store1[face_index1].add_edge(pedge, status1);
                    geom_loops_store1[face_index1].add_edge(gedge, status1);
                }
                Some(())
            })
        })
        .ok_or({
            if intersection_work > intersection_curve::MAX_TRIANGLE_PAIR_WORK {
                super::ShapeOpsError::ResourceLimit {
                    stage: "intersection triangle pair visits",
                    limit: intersection_curve::MAX_TRIANGLE_PAIR_WORK,
                }
            } else {
                super::ShapeOpsError::SurfaceIntersectionFailed {
                    first_face_index: pair.0,
                    second_face_index: pair.1,
                    stage,
                }
            }
        })?;
    if geom_loops_store0
        .iter()
        .chain(geom_loops_store1.iter())
        .flat_map(|loops| loops.iter())
        .flat_map(|wire| wire.iter())
        .any(|edge| !edge.is_geometric_consistent())
    {
        return Err(super::ShapeOpsError::SurfaceIntersectionFailed {
            first_face_index: pair.0,
            second_face_index: pair.1,
            stage: "native edge/vertex consistency",
        });
    }
    Ok(LoopsStoreQuadruple {
        geom_loops_store0,
        poly_loops_store0,
        geom_loops_store1,
        poly_loops_store1,
    })
}
