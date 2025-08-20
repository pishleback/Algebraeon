use super::*;
use crate::{
    ambient_space::AffineSpace,
    simplex::Simplex,
    simplex_collection::LabelledSimplexCollection,
    simplicial_complex::{LabelledSimplicialComplex, SimplicialComplex},
    simplicial_disjoint_union::LabelledSimplicialDisjointUnion,
};
use std::collections::{HashMap, HashSet};

#[derive(Clone)]
pub struct LabelledPartialSimplicialComplex<
    'f,
    FS: OrderedRingSignature + FieldSignature,
    T: Eq + Clone + Send + Sync,
> {
    ambient_space: AffineSpace<'f, FS>,
    simplexes: HashMap<Simplex<'f, FS>, T>,
}

pub type PartialSimplicialComplex<'f, FS> = LabelledPartialSimplicialComplex<'f, FS, ()>;

impl<'f, FS: OrderedRingSignature + FieldSignature> std::fmt::Debug
    for PartialSimplicialComplex<'f, FS>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PartialSimplicialComplex")
            .field("simplexes", &self.simplexes)
            .finish()
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature, T: Eq + Clone + Send + Sync>
    LabelledSimplexCollection<'f, FS, T> for LabelledPartialSimplicialComplex<'f, FS, T>
where
    FS::Set: Hash,
{
    type WithLabel<S: Eq + Clone + Send + Sync> = LabelledPartialSimplicialComplex<'f, FS, S>;
    type SubsetType = LabelledPartialSimplicialComplex<'f, FS, T>;

    fn try_new_labelled(
        ambient_space: AffineSpace<'f, FS>,
        simplexes: std::collections::HashMap<Simplex<'f, FS>, T>,
    ) -> Result<Self, &'static str> {
        Ok(Self {
            ambient_space,
            simplexes,
        })
    }

    fn new_labelled_unchecked(
        ambient_space: AffineSpace<'f, FS>,
        simplexes: std::collections::HashMap<Simplex<'f, FS>, T>,
    ) -> Self {
        Self::try_new_labelled(ambient_space, simplexes).unwrap()
    }

    fn ambient_space(&self) -> AffineSpace<'f, FS> {
        self.ambient_space
    }

    fn labelled_simplexes(&self) -> std::collections::HashMap<&Simplex<'f, FS>, &T> {
        self.simplexes.iter().collect()
    }

    fn into_labelled_simplexes(self) -> std::collections::HashMap<Simplex<'f, FS>, T> {
        self.simplexes
    }

    fn into_partial_simplicial_complex(self) -> LabelledPartialSimplicialComplex<'f, FS, T> {
        self
    }

    fn to_partial_simplicial_complex(&self) -> LabelledPartialSimplicialComplex<'f, FS, T> {
        self.clone()
    }

    fn into_simplicial_disjoint_union(self) -> LabelledSimplicialDisjointUnion<'f, FS, T> {
        LabelledSimplicialDisjointUnion::new_labelled_unchecked(self.ambient_space, self.simplexes)
    }

    fn to_simplicial_disjoint_union(&self) -> LabelledSimplicialDisjointUnion<'f, FS, T> {
        self.clone().into_simplicial_disjoint_union()
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature, T: Eq + Clone + Send + Sync>
    LabelledPartialSimplicialComplex<'f, FS, T>
where
    FS::Set: Hash,
{
    pub fn try_into_simplicial_complex(
        self,
    ) -> Result<LabelledSimplicialComplex<'f, FS, T>, &'static str> {
        LabelledSimplicialComplex::try_new_labelled(self.ambient_space, self.simplexes)
    }

    pub fn into_labelled_simplicial_complex(&self) -> LabelledSimplicialComplex<'f, FS, Option<T>> {
        let mut simplexes = HashSet::new();
        for spx in self.simplexes.keys() {
            for bdry in spx.sub_simplices_not_null() {
                simplexes.insert(bdry);
            }
        }
        LabelledSimplicialComplex::try_new_labelled(
            self.ambient_space(),
            simplexes
                .into_iter()
                .map(|spx| {
                    let label = self.simplexes.get(&spx).cloned();
                    (spx, label)
                })
                .collect(),
        )
        .unwrap()
    }

    pub fn simplify(&self) -> Self {
        self.into_labelled_simplicial_complex()
            .simplify()
            .subset_by_filter(|label| label.is_some())
            .apply_label_function(|label| label.clone().unwrap())
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> PartialSimplicialComplex<'f, FS>
where
    FS::Set: Hash,
{
    pub fn closure(&self) -> SimplicialComplex<'f, FS> {
        self.into_labelled_simplicial_complex().forget_labels()
    }
}
