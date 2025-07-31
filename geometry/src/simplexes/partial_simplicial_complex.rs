use super::*;
use std::collections::{HashMap, HashSet};

#[derive(Clone)]
pub struct LabelledPartialSimplicialComplex<
    'f,
    FS: OrderedRingSignature + FieldSignature,
    SP: Borrow<AffineSpace<'f, FS>> + Clone,
    T: Eq + Clone,
> {
    ambient_space: SP,
    simplexes: HashMap<Simplex<'f, FS, SP>, T>,
}

pub type PartialSimplicialComplex<'f, FS, SP> = LabelledPartialSimplicialComplex<'f, FS, SP, ()>;

impl<'f, FS: OrderedRingSignature + FieldSignature, SP: Borrow<AffineSpace<'f, FS>> + Clone>
    std::fmt::Debug for PartialSimplicialComplex<'f, FS, SP>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PartialSimplicialComplex")
            .field("simplexes", &self.simplexes)
            .finish()
    }
}

impl<
    'f,
    FS: OrderedRingSignature + FieldSignature,
    SP: Borrow<AffineSpace<'f, FS>> + Clone,
    T: Eq + Clone,
> LabelledSimplexCollection<'f, FS, SP, T> for LabelledPartialSimplicialComplex<'f, FS, SP, T>
where
    FS::Set: Hash,
{
    type WithLabel<S: Eq + Clone> = LabelledPartialSimplicialComplex<'f, FS, SP, S>;
    type SubsetType = LabelledPartialSimplicialComplex<'f, FS, SP, T>;

    fn new_labelled(
        ambient_space: SP,
        simplexes: std::collections::HashMap<Simplex<'f, FS, SP>, T>,
    ) -> Result<Self, &'static str> {
        Ok(Self {
            ambient_space,
            simplexes,
        })
    }

    fn new_labelled_unchecked(
        ambient_space: SP,
        simplexes: std::collections::HashMap<Simplex<'f, FS, SP>, T>,
    ) -> Self {
        Self::new_labelled(ambient_space, simplexes).unwrap()
    }

    fn ambient_space(&self) -> SP {
        self.ambient_space.clone()
    }

    fn labelled_simplexes(&self) -> std::collections::HashMap<&Simplex<'f, FS, SP>, &T> {
        self.simplexes.iter().collect()
    }

    fn into_labelled_simplexes(self) -> std::collections::HashMap<Simplex<'f, FS, SP>, T> {
        self.simplexes
    }

    fn into_partial_simplicial_complex(self) -> LabelledPartialSimplicialComplex<'f, FS, SP, T> {
        self
    }
}

impl<
    'f,
    FS: OrderedRingSignature + FieldSignature,
    SP: Borrow<AffineSpace<'f, FS>> + Clone,
    T: Eq + Clone,
> LabelledPartialSimplicialComplex<'f, FS, SP, T>
where
    FS::Set: Hash,
{
    pub fn try_as_simplicial_complex(
        self,
    ) -> Result<LabelledSimplicialComplex<'f, FS, SP, T>, &'static str> {
        LabelledSimplicialComplex::new_labelled(self.ambient_space, self.simplexes)
    }

    pub fn closure(&self) -> LabelledSimplicialComplex<'f, FS, SP, Option<T>> {
        let mut simplexes = HashSet::new();
        #[allow(clippy::for_kv_map)]
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
        #[allow(clippy::redundant_closure_for_method_calls)]
        self.closure()
            .simplify()
            .subset_by_filter(|label| label.is_some())
            .apply_label_function(|label| label.clone().unwrap())
    }
}
