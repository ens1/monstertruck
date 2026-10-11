//! Classic (0.3.2) per-face and/or/unknown classification.
//!
//! Ported verbatim from the published 0.3.2 crate's
//! `transversal::faces_classification`. Kept private to the classic subtree so
//! it binds to the classic [`ShapesOpStatus`] (the crate-level
//! `faces_classification` binds to the loops-store status of the same
//! name, a distinct type).

use super::loops_store::ShapesOpStatus;
use monstertruck_topology::*;
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

#[derive(Clone, Debug)]
pub(super) struct FacesClassification<P, C, S> {
    shell: Shell<P, C, S>,
    status: HashMap<FaceId<S>, ShapesOpStatus>,
}

impl<P, C, S> Default for FacesClassification<P, C, S> {
    fn default() -> Self {
        Self {
            shell: Default::default(),
            status: HashMap::default(),
        }
    }
}

impl<P, C, S> FacesClassification<P, C, S> {
    pub(super) fn push(&mut self, face: Face<P, C, S>, status: ShapesOpStatus) {
        self.status.insert(face.id(), status);
        self.shell.push(face);
    }

    pub(super) fn and_or_unknown(&self) -> [Shell<P, C, S>; 3] {
        let [mut and, mut or, mut unknown] = <[Shell<P, C, S>; 3]>::default();
        for face in &self.shell {
            // SAFETY: `push()` inserts every face id into `status`.
            match self.status.get(&face.id()).unwrap() {
                ShapesOpStatus::And => and.push(face.clone()),
                ShapesOpStatus::Or => or.push(face.clone()),
                ShapesOpStatus::Unknown => unknown.push(face.clone()),
            }
        }
        [and, or, unknown]
    }

    pub(super) fn integrate_by_component(&mut self) -> Option<()> {
        let [and, or, unknown] = self.and_or_unknown();
        let and_boundary: HashSet<_> = and
            .extract_boundaries()
            .into_iter()
            .flatten()
            .map(|edge| edge.id())
            .collect();
        let or_boundary: HashSet<_> = or
            .extract_boundaries()
            .into_iter()
            .flatten()
            .map(|edge| edge.id())
            .collect();
        for comp in unknown.connected_components() {
            let boundary = comp.extract_boundaries();
            // Closed untouched components have no boundary to inherit from.
            // Open components must agree at every classified neighboring edge;
            // taking the first edge made classification depend on iteration order.
            let mut inherited = None;
            for edge in boundary.iter().flatten() {
                let status = match (
                    and_boundary.contains(&edge.id()),
                    or_boundary.contains(&edge.id()),
                ) {
                    (true, false) => ShapesOpStatus::And,
                    (false, true) => ShapesOpStatus::Or,
                    (false, false) => continue,
                    (true, true) => return None,
                };
                if inherited.is_some_and(|previous| previous != status) {
                    return None;
                }
                inherited = Some(status);
            }
            if let Some(status) = inherited {
                for face in &comp {
                    *self.status.get_mut(&face.id())? = status;
                }
            }
        }
        Some(())
    }
}
