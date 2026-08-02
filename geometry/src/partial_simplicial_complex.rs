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
    FS: OrderedRingSignature + FieldSignature,
    T: Eq + Clone + Send + Sync,
> {
    ambient_space: AffineSpace<FS>,
    simplexes: HashMap<Simplex<FS>, T>,
}

pub type PartialSimplicialComplex<FS> = LabelledPartialSimplicialComplex<FS, ()>;

impl<FS: OrderedRingSignature + FieldSignature> std::fmt::Debug for PartialSimplicialComplex<FS> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PartialSimplicialComplex")
            .field("simplexes", &self.simplexes)
            .finish()
    }
}

impl<FS: OrderedRingSignature + FieldSignature, T: Eq + Clone + Send + Sync>
    LabelledSimplexCollection<FS, T> for LabelledPartialSimplicialComplex<FS, T>
where
    FS::Elem: Hash,
{
    type WithLabel<S: Eq + Clone + Send + Sync> = LabelledPartialSimplicialComplex<FS, S>;
    type SubsetType = LabelledPartialSimplicialComplex<FS, T>;

    fn try_new_labelled(
        ambient_space: &AffineSpace<FS>,
        simplexes: std::collections::HashMap<Simplex<FS>, T>,
    ) -> Result<Self, &'static str> {
        Ok(Self {
            ambient_space: ambient_space.clone(),
            simplexes,
        })
    }

    fn new_labelled_unchecked(
        ambient_space: &AffineSpace<FS>,
        simplexes: std::collections::HashMap<Simplex<FS>, T>,
    ) -> Self {
        Self::try_new_labelled(&ambient_space, simplexes).unwrap()
    }

    fn ambient_space(&self) -> &AffineSpace<FS> {
        &self.ambient_space
    }

    fn labelled_simplexes(&self) -> std::collections::HashMap<&Simplex<FS>, &T> {
        self.simplexes.iter().collect()
    }

    fn into_labelled_simplexes(self) -> std::collections::HashMap<Simplex<FS>, T> {
        self.simplexes
    }

    fn into_partial_simplicial_complex(self) -> LabelledPartialSimplicialComplex<FS, T> {
        self
    }

    fn to_partial_simplicial_complex(&self) -> LabelledPartialSimplicialComplex<FS, T> {
        self.clone()
    }

    fn into_simplicial_disjoint_union(self) -> LabelledSimplicialDisjointUnion<FS, T> {
        LabelledSimplicialDisjointUnion::new_labelled_unchecked(&self.ambient_space, self.simplexes)
    }

    fn to_simplicial_disjoint_union(&self) -> LabelledSimplicialDisjointUnion<FS, T> {
        self.clone().into_simplicial_disjoint_union()
    }
}

impl<FS: OrderedRingSignature + FieldSignature, T: Eq + Clone + Send + Sync>
    LabelledPartialSimplicialComplex<FS, T>
where
    FS::Elem: Hash,
{
    pub fn try_into_simplicial_complex(
        self,
    ) -> Result<LabelledSimplicialComplex<FS, T>, &'static str> {
        LabelledSimplicialComplex::try_new_labelled(&self.ambient_space, self.simplexes)
    }

    pub fn into_labelled_simplicial_complex(&self) -> LabelledSimplicialComplex<FS, Option<T>> {
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

impl<FS: OrderedRingSignature + FieldSignature> PartialSimplicialComplex<FS>
where
    FS::Elem: Hash,
{
    pub fn closure(&self) -> SimplicialComplex<FS> {
        self.into_labelled_simplicial_complex().forget_labels()
    }
}
