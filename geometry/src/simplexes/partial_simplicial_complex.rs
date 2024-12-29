use std::collections::{HashMap, HashSet};

use super::*;

#[derive(Clone)]
pub struct LabelledPartialSimplicialComplex<
    FS: OrderedRingStructure + FieldStructure,
    SP: Borrow<AffineSpace<FS>> + Clone,
    T: Eq + Clone,
> {
    ambient_space: SP,
    simplexes: HashMap<Simplex<FS, SP>, T>,
}

pub type PartialSimplicialComplex<FS, SP> = LabelledPartialSimplicialComplex<FS, SP, ()>;

impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone> std::fmt::Debug
    for PartialSimplicialComplex<FS, SP>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PartialSimplicialComplex")
            .field("simplexes", &self.simplexes)
            .finish()
    }
}

impl<
        FS: OrderedRingStructure + FieldStructure,
        SP: Borrow<AffineSpace<FS>> + Clone,
        T: Eq + Clone,
    > LabelledSimplexCollection<FS, SP, T> for LabelledPartialSimplicialComplex<FS, SP, T>
where
    FS::Set: Hash,
{
    type WithLabel<S: Eq + Clone> = LabelledPartialSimplicialComplex<FS, SP, S>;
    type SubsetType = LabelledPartialSimplicialComplex<FS, SP, T>;

    fn new_labelled(
        ambient_space: SP,
        simplexes: std::collections::HashMap<Simplex<FS, SP>, T>,
    ) -> Result<Self, &'static str> {
        Ok(Self {
            ambient_space,
            simplexes,
        })
    }

    fn new_labelled_unchecked(
        ambient_space: SP,
        simplexes: std::collections::HashMap<Simplex<FS, SP>, T>,
    ) -> Self {
        Self::new_labelled(ambient_space, simplexes).unwrap()
    }

    fn ambient_space(&self) -> SP {
        self.ambient_space.clone()
    }

    fn labelled_simplexes(&self) -> std::collections::HashMap<&Simplex<FS, SP>, &T> {
        self.simplexes.iter().collect()
    }

    fn into_labelled_simplexes(self) -> std::collections::HashMap<Simplex<FS, SP>, T> {
        self.simplexes
    }

    fn into_partial_simplicial_complex(self) -> LabelledPartialSimplicialComplex<FS, SP, T> {
        self
    }
}

impl<
        FS: OrderedRingStructure + FieldStructure,
        SP: Borrow<AffineSpace<FS>> + Clone,
        T: Eq + Clone,
    > LabelledPartialSimplicialComplex<FS, SP, T>
where
    FS::Set: Hash,
{
    pub fn try_as_simplicial_complex(
        self,
    ) -> Result<LabelledSimplicialComplex<FS, SP, T>, &'static str> {
        LabelledSimplicialComplex::new_labelled(self.ambient_space, self.simplexes)
    }

    pub fn closure(&self) -> LabelledSimplicialComplex<FS, SP, Option<T>> {
        let mut simplexes = HashSet::new();
        for (spx, _label) in &self.simplexes {
            for bdry in spx.sub_simplices_not_null() {
                simplexes.insert(bdry);
            }
        }
        LabelledSimplicialComplex::new_labelled(
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
        self.closure()
            .simplify()
            .subset_by_filter(|label| label.is_some())
            .apply_label_function(|label| label.clone().unwrap())
    }
}
