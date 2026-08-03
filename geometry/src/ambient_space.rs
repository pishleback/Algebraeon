use super::*;
use crate::{
    simplex_collection::LabelledSimplexCollection,
    simplicial_disjoint_union::SimplicialDisjointUnion, vector::Vector,
};
use std::sync::{Arc, atomic::AtomicUsize};

/// An affine space over a field.
/// affine_dimension = 0 => the empty space
/// affine_dimension = 1 => one point space
/// affine_dimension = 2 => a line
/// affine_dimension = 3 => a plane
/// ...
#[derive(Debug, Clone)]
pub struct AffineSpace<FS: FieldSignature> {
    field: Arc<FS>,
    // linear dimension = affine dimension - 1
    affine_dimension: usize,
    ident: usize,
}

impl<FS: FieldSignature> PartialEq for AffineSpace<FS> {
    fn eq(&self, other: &Self) -> bool {
        #[cfg(debug_assertions)]
        if self.ident == other.ident {
            assert_eq!(self.affine_dimension, other.affine_dimension);
            assert_eq!(self.field(), other.field());
        }
        self.ident == other.ident
    }
}

impl<FS: FieldSignature> Eq for AffineSpace<FS> {}

impl<FS: FieldSignature + Hash> Hash for AffineSpace<FS> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.ident.hash(state);
    }
}

impl<FS: FieldSignature> AffineSpace<FS> {
    pub fn new_affine(field: Arc<FS>, affine_dimension: usize) -> Self {
        static COUNTER: AtomicUsize = AtomicUsize::new(0);
        Self {
            field,
            affine_dimension,
            ident: COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed),
        }
    }

    pub fn new_empty(field: Arc<FS>) -> Self {
        Self::new_affine(field, 0)
    }

    pub fn new_linear(field: Arc<FS>, linear_dimension: usize) -> Self {
        Self::new_affine(field, linear_dimension + 1)
    }

    pub fn field(&self) -> &Arc<FS> {
        &self.field
    }

    pub fn origin(&self) -> Option<Vector<FS>> {
        Some(Vector::construct(self, |_| self.field().zero()))
    }

    pub fn empty_subset(&self) -> impl LabelledSimplexCollection<FS, ()>
    where
        FS: OrderedRingSignature,
        FS::Elem: Hash,
    {
        SimplicialDisjointUnion::new_unchecked(self, std::collections::HashSet::new())
    }

    pub fn linear_dimension(&self) -> Option<usize> {
        if self.affine_dimension == 0 {
            None
        } else {
            Some(self.affine_dimension - 1)
        }
    }

    pub fn affine_dimension(&self) -> usize {
        self.affine_dimension
    }
}

pub fn common_space<'s, FS: FieldSignature>(
    space1: &'s AffineSpace<FS>,
    space2: &'s AffineSpace<FS>,
) -> Option<&'s AffineSpace<FS>> {
    if space1 == space2 { Some(space1) } else { None }
}
