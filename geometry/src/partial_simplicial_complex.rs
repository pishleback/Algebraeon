use super::*;
use crate::{
    ambient_space::AffineSpace, simplex::Simplex, simplex_collection::LabelledSimplexCollection,
    simplicial_complex::LabelledSimplicialComplex,
};
use std::collections::{HashMap, HashSet};

#[derive(Clone)]
pub struct LabelledPartialSimplicialComplex<
    'f,
    FS: OrderedRingSignature + FieldSignature,
    T: Eq + Clone,
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

impl<'f, FS: OrderedRingSignature + FieldSignature, T: Eq + Clone>
    LabelledSimplexCollection<'f, FS, T> for LabelledPartialSimplicialComplex<'f, FS, T>
where
    FS::Set: Hash,
{
    type WithLabel<S: Eq + Clone> = LabelledPartialSimplicialComplex<'f, FS, S>;
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
}

impl<'f, FS: OrderedRingSignature + FieldSignature, T: Eq + Clone>
    LabelledPartialSimplicialComplex<'f, FS, T>
where
    FS::Set: Hash,
{
    pub fn try_as_simplicial_complex(
        self,
    ) -> Result<LabelledSimplicialComplex<'f, FS, T>, &'static str> {
        LabelledSimplicialComplex::try_new_labelled(self.ambient_space, self.simplexes)
    }

    pub fn closure(&self) -> LabelledSimplicialComplex<'f, FS, Option<T>> {
        let mut simplexes = HashSet::new();
        #[allow(clippy::for_kv_map)]
        for (spx, _label) in &self.simplexes {
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
        #[allow(clippy::redundant_closure_for_method_calls)]
        self.closure()
            .simplify()
            .subset_by_filter(|label| label.is_some())
            .apply_label_function(|label| label.clone().unwrap())
    }
}
