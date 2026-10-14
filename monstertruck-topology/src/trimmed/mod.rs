use crate::compress::{
    CompressDirector, CompressedEdgeIndex, CompressedEdgeUse, CompressedTrimmedFace,
    CompressedTrimmedShell, CompressedTrimmedSolid,
};
use crate::{errors::Error, *};
use std::collections::{HashMap, HashSet};

/// Destructured parts of a [`TrimmedFace`].
pub type TrimmedFaceParts<P, C, S, T> = (Face<P, C, S>, Vec<Vec<Option<T>>>);

/// Runtime face with face-local trim storage per edge-use.
#[derive(Clone, Debug)]
pub struct TrimmedFace<P, C, S, T> {
    face: Face<P, C, S>,
    trims: Vec<Vec<Option<T>>>,
}

/// Runtime shell with face-local trim storage per edge-use.
#[derive(Clone, Debug)]
pub struct TrimmedShell<P, C, S, T> {
    face_list: Vec<TrimmedFace<P, C, S, T>>,
}

/// Runtime solid with face-local trim storage per edge-use.
#[derive(Clone, Debug)]
pub struct TrimmedSolid<P, C, S, T> {
    boundaries: Vec<TrimmedShell<P, C, S, T>>,
    id_allocator: StableIdAllocator,
    attributes: SolidAttributes,
}

impl<P, C, S, T> TrimmedFace<P, C, S, T> {
    /// Creates a trimmed face from a face and per-edge-use trims.
    ///
    /// The `trims` layout must match [`Face::absolute_boundaries`].
    pub fn try_new(face: Face<P, C, S>, trims: Vec<Vec<Option<T>>>) -> Result<Self> {
        let same_wire_count = face.absolute_boundaries().len() == trims.len();
        let same_edge_counts = same_wire_count
            && face
                .absolute_boundaries()
                .iter()
                .zip(trims.iter())
                .all(|(wire, trim_wire)| wire.len() == trim_wire.len());
        if same_edge_counts {
            Ok(Self { face, trims })
        } else {
            Err(Error::TrimLayout)
        }
    }

    /// Creates a trimmed face from a face and per-edge-use trims.
    /// # Panics
    /// Panics if the trim layout does not match [`Face::absolute_boundaries`].
    pub fn new(face: Face<P, C, S>, trims: Vec<Vec<Option<T>>>) -> Self {
        Self::try_new(face, trims).expect("TrimmedFace::new: trim layout mismatch")
    }

    /// Returns the wrapped face.
    #[inline(always)]
    pub const fn face(&self) -> &Face<P, C, S> { &self.face }

    /// Returns the face-local trims matching [`Face::absolute_boundaries`].
    #[inline(always)]
    pub const fn trims(&self) -> &Vec<Vec<Option<T>>> { &self.trims }

    /// Consumes `self` and drops trim storage.
    #[inline(always)]
    pub fn erase_trims(self) -> Face<P, C, S> { self.face }

    /// Destructures into the face and trim storage.
    #[inline(always)]
    pub fn into_parts(self) -> TrimmedFaceParts<P, C, S, T> { (self.face, self.trims) }

    /// Reverses the face orientation without changing absolute edge-use trims.
    pub fn invert(&mut self) { self.face.invert(); }
}

impl<P, C, S, T> From<Face<P, C, S>> for TrimmedFace<P, C, S, T> {
    fn from(face: Face<P, C, S>) -> Self {
        let trims = face
            .absolute_boundaries()
            .iter()
            .map(|wire| wire.iter().map(|_| None).collect())
            .collect();
        Self { face, trims }
    }
}

impl<P, C, S, T> From<TrimmedFace<P, C, S, T>> for Face<P, C, S> {
    fn from(face: TrimmedFace<P, C, S, T>) -> Self { face.erase_trims() }
}

impl<P, C, S, T> FromIterator<TrimmedFace<P, C, S, T>> for TrimmedShell<P, C, S, T> {
    fn from_iter<I: IntoIterator<Item = TrimmedFace<P, C, S, T>>>(iter: I) -> Self {
        Self {
            face_list: iter.into_iter().collect(),
        }
    }
}

impl<P, C, S, T> TrimmedShell<P, C, S, T> {
    /// Returns the faces.
    #[inline(always)]
    pub const fn faces(&self) -> &Vec<TrimmedFace<P, C, S, T>> { &self.face_list }

    /// Consumes `self` and drops trim storage.
    pub fn erase_trims(self) -> Shell<P, C, S> {
        self.face_list
            .into_iter()
            .map(TrimmedFace::erase_trims)
            .collect()
    }

    /// A plain [`Shell`] view of the faces, trims dropped, without consuming
    /// `self`. Face clones are handle clones; no geometry is copied.
    fn to_plain_shell(&self) -> Shell<P, C, S> {
        self.face_list
            .iter()
            .map(|trimmed_face| trimmed_face.face.clone())
            .collect()
    }

    /// The precondition this shell must satisfy to be a boundary of a
    /// [`TrimmedSolid`]; see [`Shell::check_solid_boundary`].
    pub fn check_solid_boundary(&self) -> Result<()> {
        self.to_plain_shell().check_solid_boundary()
    }

    /// Reverses every face while retaining trims in the absolute carrier frame.
    pub fn invert(&mut self) { self.face_list.iter_mut().for_each(TrimmedFace::invert); }

    /// Maps shared topology and each face-local trim without erasing trim storage.
    ///
    /// The caller must map a trim's carrier consistently with `surface_mapping`.
    /// A failed geometry or trim mapping rejects the complete result.
    pub fn try_mapped<Q, D, U, V>(
        &self,
        point_mapping: impl FnMut(&P) -> Option<Q>,
        curve_mapping: impl FnMut(&C) -> Option<D>,
        surface_mapping: impl FnMut(&S) -> Option<U>,
        mut trim_mapping: impl FnMut(&T) -> Option<V>,
    ) -> Option<TrimmedShell<Q, D, U, V>> {
        let mapped =
            self.to_plain_shell()
                .try_mapped(point_mapping, curve_mapping, surface_mapping)?;
        mapped
            .into_iter()
            .zip(&self.face_list)
            .map(|(mut face, original)| {
                // Ordinary Shell::try_mapped preserves shared geometry and
                // absolute use order, but newly constructs edge/face IDs.
                // This mapping is one-to-one topology; retain their stable IDs
                // so sparse attributes continue to refer to the mapped elements.
                face.set_stable_id(original.face.stable_id());
                for (wire, source) in face
                    .boundaries
                    .iter_mut()
                    .zip(original.face.absolute_boundaries())
                {
                    for (edge, source) in wire.iter_mut().zip(source) {
                        edge.set_stable_id(source.stable_id());
                    }
                }
                let trims = original
                    .trims
                    .iter()
                    .map(|wire| {
                        wire.iter()
                            .map(|trim| match trim {
                                Some(trim) => trim_mapping(trim).map(Some),
                                None => Some(None),
                            })
                            .collect::<Option<Vec<_>>>()
                    })
                    .collect::<Option<Vec<_>>>()?;
                TrimmedFace::try_new(face, trims).ok()
            })
            .collect()
    }
}

impl<P, C, S, T> From<Shell<P, C, S>> for TrimmedShell<P, C, S, T> {
    fn from(shell: Shell<P, C, S>) -> Self { shell.into_iter().map(TrimmedFace::from).collect() }
}

impl<P, C, S, T> From<TrimmedShell<P, C, S, T>> for Shell<P, C, S> {
    fn from(shell: TrimmedShell<P, C, S, T>) -> Self { shell.erase_trims() }
}

impl<P, C, S, T> TrimmedSolid<P, C, S, T> {
    /// Creates a trimmed solid from its boundary shells, WITHOUT validating
    /// them.
    ///
    /// # This is `Solid::new_unchecked`, not `Solid::new`
    ///
    /// Despite the shared name, this constructor is the trimmed counterpart of
    /// [`Solid::new_unchecked`] and not of [`Solid::new`]: it checks nothing.
    /// Spec 011's ledger entry for class C11, and the corpus row that pinned
    /// it, both recorded this call as "the infallible `Solid::new`" -- it is
    /// not, it never validated at all, and the abort they were chasing came
    /// from a `Solid::debug_new` two hops downstream in
    /// [`Solid::try_mapped`]. That confusion is the reason this doc comment is
    /// this long.
    ///
    /// Use [`Self::try_new`] whenever the boundaries come from data rather
    /// than from a construction that already established the invariant.
    #[inline(always)]
    pub fn new(boundaries: Vec<TrimmedShell<P, C, S, T>>) -> Self {
        Self {
            boundaries,
            id_allocator: StableIdAllocator::new(),
            attributes: SolidAttributes::new(),
        }
    }

    /// Creates a trimmed solid whose boundaries must be non-empty, connected,
    /// and closed manifold -- the SAME precondition as [`Solid::try_new`],
    /// via the single definition in [`Shell::check_solid_boundary`].
    ///
    /// The plain path has enforced this since forever
    /// (`Solid::extract` ends in `Solid::try_new(shells?)?`); the trimmed path
    /// did not, so a `CompressedTrimmedSolid` that was short a face -- the
    /// shape a typed surface refusal upstream produces -- extracted to `Ok`
    /// and only failed later, in a `debug_assertions`-only construction that
    /// panicked in debug and silently returned an invalid solid in release.
    /// Ledger class C11.
    #[inline(always)]
    pub fn try_new(boundaries: Vec<TrimmedShell<P, C, S, T>>) -> Result<Self> {
        for shell in &boundaries {
            shell.check_solid_boundary()?;
        }
        let mut result = Self::new(boundaries);
        result.reconcile_allocator()?;
        Ok(result)
    }

    /// Returns the boundary shells.
    #[inline(always)]
    pub const fn boundaries(&self) -> &Vec<TrimmedShell<P, C, S, T>> { &self.boundaries }

    /// Consumes `self` and returns the boundary shells.
    #[inline(always)]
    pub fn into_boundaries(self) -> Vec<TrimmedShell<P, C, S, T>> { self.boundaries }

    /// Allocator state retained when adding or removing face-local trims.
    pub const fn id_allocator(&self) -> &StableIdAllocator { &self.id_allocator }
    /// Sparse attributes keyed by the stable IDs of this solid's topology.
    pub const fn attributes(&self) -> &SolidAttributes { &self.attributes }
    /// Edit sparse attributes without replacing topology or trim authority.
    pub fn attributes_mut(&mut self) -> &mut SolidAttributes { &mut self.attributes }

    fn reconcile_allocator(&mut self) -> Result<()> {
        let mut highest = 0;
        let mut face_ids = HashMap::new();
        let mut edge_ids = HashMap::new();
        let mut vertex_ids = HashMap::new();
        fn unique<I: Eq>(ids: &mut HashMap<StableId, I>, stable: StableId, id: I) -> Result<()> {
            if stable.is_assigned()
                && ids
                    .insert(stable, id)
                    .is_some_and(|previous| ids.get(&stable) != Some(&previous))
            {
                return Err(Error::StableIdMetadata(
                    "one stable ID identifies multiple distinct topology elements",
                ));
            }
            Ok(())
        }
        for face in self.boundaries.iter().flat_map(|shell| shell.faces()) {
            unique(&mut face_ids, face.face.stable_id(), face.face.id())?;
            highest = highest.max(face.face.stable_id().raw());
            for edge in face.face.edge_iter() {
                unique(&mut edge_ids, edge.stable_id(), edge.id())?;
                for vertex in [edge.absolute_front(), edge.absolute_back()] {
                    unique(&mut vertex_ids, vertex.stable_id(), vertex.id())?;
                }
                highest = highest
                    .max(edge.stable_id().raw())
                    .max(edge.absolute_front().stable_id().raw())
                    .max(edge.absolute_back().stable_id().raw());
            }
        }
        if highest >= u64::MAX - 1
            || self.id_allocator.peek() == 0
            || self.id_allocator.peek() == u64::MAX
        {
            return Err(Error::StableIdMetadata(
                "allocator or topology IDs exhaust the finite allocation range",
            ));
        }
        self.id_allocator.ensure_above(highest);
        Ok(())
    }

    /// Maps one-to-one native topology and all trim carriers while retaining
    /// stable IDs, sparse attributes, and allocator continuity. Invalid mapped
    /// boundaries or any failed mapping reject the complete result.
    pub fn try_mapped<Q, D, U, V>(
        &self,
        mut point_mapping: impl FnMut(&P) -> Option<Q>,
        mut curve_mapping: impl FnMut(&C) -> Option<D>,
        mut surface_mapping: impl FnMut(&S) -> Option<U>,
        mut trim_mapping: impl FnMut(&T) -> Option<V>,
    ) -> Option<TrimmedSolid<Q, D, U, V>> {
        let boundaries = self
            .boundaries
            .iter()
            .map(|shell| {
                shell.try_mapped(
                    &mut point_mapping,
                    &mut curve_mapping,
                    &mut surface_mapping,
                    &mut trim_mapping,
                )
            })
            .collect::<Option<Vec<_>>>()?;
        let mut mapped = TrimmedSolid::try_new(boundaries).ok()?;
        mapped.id_allocator = self.id_allocator.clone();
        mapped.attributes = self.attributes.clone();
        mapped.reconcile_allocator().ok()?;
        Some(mapped)
    }

    /// Consumes `self` and drops trim storage.
    pub fn erase_trims(self) -> Solid<P, C, S> {
        Solid {
            boundaries: self
                .boundaries
                .into_iter()
                .map(TrimmedShell::erase_trims)
                .collect(),
            id_allocator: self.id_allocator,
            attributes: self.attributes,
        }
    }
}

impl<P, C, S, T> From<Solid<P, C, S>> for TrimmedSolid<P, C, S, T> {
    fn from(solid: Solid<P, C, S>) -> Self {
        let Solid {
            boundaries,
            id_allocator,
            attributes,
        } = solid;
        Self {
            boundaries: boundaries.into_iter().map(TrimmedShell::from).collect(),
            id_allocator,
            attributes,
        }
    }
}

impl<P, C, S, T> From<TrimmedSolid<P, C, S, T>> for Solid<P, C, S> {
    fn from(solid: TrimmedSolid<P, C, S, T>) -> Self { solid.erase_trims() }
}

impl<P, C, S, T> FromIterator<TrimmedShell<P, C, S, T>> for TrimmedSolid<P, C, S, T> {
    fn from_iter<I: IntoIterator<Item = TrimmedShell<P, C, S, T>>>(iter: I) -> Self {
        Self::new(iter.into_iter().collect())
    }
}

impl<P: Clone, C: Clone, S: Clone> Shell<P, C, S> {
    /// Creates a runtime trimmed shell by evaluating `trim_curve` for each edge-use.
    pub fn to_trimmed_with_face_trims<T, F>(&self, mut trim_curve: F) -> TrimmedShell<P, C, S, T>
    where F: FnMut(&Edge<P, C>, &S) -> Option<T> {
        self.iter()
            .map(|face| {
                let surface = face.surface();
                let trims = face
                    .absolute_boundaries()
                    .iter()
                    .map(|wire| wire.iter().map(|edge| trim_curve(edge, &surface)).collect())
                    .collect();
                TrimmedFace::new(face.clone(), trims)
            })
            .collect()
    }

    /// Creates a runtime trimmed shell from exact trims carried by the edge curves.
    pub fn to_trimmed_with_exact_face_trims(
        &self,
    ) -> TrimmedShell<P, C, S, <C as ExactParameterBoundary2D<S>>::BoundaryCurve>
    where C: ExactParameterBoundary2D<S> {
        self.to_trimmed_with_face_trims(|edge, surface| {
            edge.curve().exact_parameter_boundary_2d(surface)
        })
    }
}

impl<P: Clone, C: Clone, S: Clone> Solid<P, C, S> {
    /// Creates a runtime trimmed solid by evaluating `trim_curve` for each edge-use.
    pub fn to_trimmed_with_face_trims<T, F>(&self, mut trim_curve: F) -> TrimmedSolid<P, C, S, T>
    where F: FnMut(&Edge<P, C>, &S) -> Option<T> {
        TrimmedSolid {
            boundaries: self
                .boundaries()
                .iter()
                .map(|shell| shell.to_trimmed_with_face_trims(&mut trim_curve))
                .collect(),
            id_allocator: self.id_allocator.clone(),
            attributes: self.attributes.clone(),
        }
    }

    /// Creates a runtime trimmed solid from exact trims carried by the edge curves.
    pub fn to_trimmed_with_exact_face_trims(
        &self,
    ) -> TrimmedSolid<P, C, S, <C as ExactParameterBoundary2D<S>>::BoundaryCurve>
    where C: ExactParameterBoundary2D<S> {
        self.to_trimmed_with_face_trims(|edge, surface| {
            edge.curve().exact_parameter_boundary_2d(surface)
        })
    }
}

impl<P: Clone, C: Clone, S: Clone, T: Clone> From<&TrimmedShell<P, C, S, T>>
    for CompressedTrimmedShell<P, C, S, T>
{
    fn from(shell: &TrimmedShell<P, C, S, T>) -> Self {
        let mut director = CompressDirector::new();
        let faces = shell
            .face_list
            .iter()
            .map(|face| CompressedTrimmedFace {
                boundaries: face
                    .face
                    .absolute_boundaries()
                    .iter()
                    .zip(face.trims.iter())
                    .map(|(wire, trim_wire)| {
                        wire.iter()
                            .zip(trim_wire.iter())
                            .map(|(edge, trim_curve)| {
                                let CompressedEdgeIndex { index, orientation } =
                                    director.get_eid(edge);
                                CompressedEdgeUse {
                                    index,
                                    orientation,
                                    trim_curve: trim_curve.clone(),
                                }
                            })
                            .collect()
                    })
                    .collect(),
                orientation: face.face.orientation(),
                surface: face.face.surface(),
            })
            .collect();
        let (vertices, edges, vsids, esids) = director.stable_ids();
        let fsids = shell
            .face_list
            .iter()
            .map(|face| face.face.stable_id())
            .collect::<Vec<_>>();
        let assigned = vsids
            .iter()
            .chain(&esids)
            .chain(&fsids)
            .any(|id| id.is_assigned());
        Self {
            vertices,
            edges,
            faces,
            vertex_stable_ids: assigned.then_some(vsids),
            edge_stable_ids: assigned.then_some(esids),
            face_stable_ids: assigned.then_some(fsids),
        }
    }
}

impl<P: Clone, C: Clone, S: Clone, T: Clone> From<&TrimmedSolid<P, C, S, T>>
    for CompressedTrimmedSolid<P, C, S, T>
{
    fn from(solid: &TrimmedSolid<P, C, S, T>) -> Self {
        Self {
            boundaries: solid
                .boundaries
                .iter()
                .map(CompressedTrimmedShell::from)
                .collect(),
            id_allocator: Some(solid.id_allocator.clone()),
            attributes: if solid.attributes.is_empty() {
                None
            } else {
                Some(solid.attributes.clone())
            },
        }
    }
}

impl<P, C, S, T> TryFrom<CompressedTrimmedShell<P, C, S, T>> for TrimmedShell<P, C, S, T> {
    type Error = Error;

    fn try_from(shell: CompressedTrimmedShell<P, C, S, T>) -> Result<Self> {
        let CompressedTrimmedShell {
            vertices,
            edges,
            faces,
            vertex_stable_ids,
            edge_stable_ids,
            face_stable_ids,
        } = shell;
        fn metadata(ids: Option<&[StableId]>, length: usize) -> Result<()> {
            if let Some(ids) = ids {
                if ids.len() != length {
                    return Err(Error::StableIdMetadata(
                        "ID array length does not match its topology collection",
                    ));
                }
                let mut assigned = HashSet::new();
                if ids.iter().any(|id| {
                    id.raw() >= u64::MAX - 1 || (id.is_assigned() && !assigned.insert(*id))
                }) {
                    return Err(Error::StableIdMetadata(
                        "ID array contains a repeated or exhausted assigned identity",
                    ));
                }
            }
            Ok(())
        }
        metadata(vertex_stable_ids.as_deref(), vertices.len())?;
        metadata(edge_stable_ids.as_deref(), edges.len())?;
        metadata(face_stable_ids.as_deref(), faces.len())?;
        let mut vertices: Vec<_> = vertices.into_iter().map(Vertex::new).collect();
        if let Some(ids) = vertex_stable_ids {
            for (vertex, id) in vertices.iter_mut().zip(ids) {
                vertex.set_stable_id(id);
            }
        }
        let mut edges = edges
            .into_iter()
            .map(|edge| {
                let vertex = |index| {
                    vertices.get(index).ok_or(Error::InvalidIndex {
                        kind: "vertex",
                        index,
                        length: vertices.len(),
                    })
                };
                let front = vertex(edge.vertices.0)?;
                let back = vertex(edge.vertices.1)?;
                Edge::try_new(front, back, edge.curve)
            })
            .collect::<Result<Vec<_>>>()?;
        if let Some(ids) = edge_stable_ids {
            for (edge, id) in edges.iter_mut().zip(ids) {
                edge.set_stable_id(id);
            }
        }
        faces
            .into_iter()
            .enumerate()
            .map(|(index, face)| {
                let boundaries = face
                    .boundaries
                    .iter()
                    .map(|wire| {
                        wire.iter()
                            .map(|edge_use| {
                                let edge =
                                    edges.get(edge_use.index).ok_or(Error::InvalidIndex {
                                        kind: "edge",
                                        index: edge_use.index,
                                        length: edges.len(),
                                    })?;
                                Ok(if edge_use.orientation {
                                    edge.clone()
                                } else {
                                    edge.inverse()
                                })
                            })
                            .collect::<Result<Wire<P, C>>>()
                    })
                    .collect::<Result<Vec<Wire<P, C>>>>()?;
                let trims = face
                    .boundaries
                    .into_iter()
                    .map(|wire| {
                        wire.into_iter()
                            .map(|edge_use| edge_use.trim_curve)
                            .collect()
                    })
                    .collect();
                let mut runtime_face = Face::try_new(boundaries, face.surface)?;
                if !face.orientation {
                    runtime_face.invert();
                }
                if let Some(ids) = &face_stable_ids {
                    runtime_face.set_stable_id(ids[index]);
                }
                TrimmedFace::try_new(runtime_face, trims)
            })
            .collect()
    }
}

impl<P, C, S, T> TryFrom<CompressedTrimmedSolid<P, C, S, T>> for TrimmedSolid<P, C, S, T> {
    type Error = Error;

    fn try_from(solid: CompressedTrimmedSolid<P, C, S, T>) -> Result<Self> {
        let mut result = Self::try_new(
            solid
                .boundaries
                .into_iter()
                .map(TrimmedShell::try_from)
                .collect::<Result<Vec<_>>>()?,
        )?;
        result.id_allocator = solid.id_allocator.unwrap_or_default();
        result.attributes = solid.attributes.unwrap_or_default();
        result.reconcile_allocator()?;
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
