use std::collections::{HashMap, HashSet};

use super::*;

pub trait LabelledSimplexCollection<
    FS: OrderedRingSignature + FieldSignature,
    SP: Borrow<AffineSpace<FS>> + Clone,
    T: Eq + Clone,
>: Sized where
    FS::Set: Hash,
{
    type WithLabel<S: Eq + Clone>: LabelledSimplexCollection<FS, SP, S>;
    type SubsetType: LabelledSimplexCollection<FS, SP, T>;

    fn new(
        ambient_space: SP,
        simplexes: HashSet<Simplex<FS, SP>>,
    ) -> Result<Self::WithLabel<()>, &'static str> {
        Self::WithLabel::<()>::new_labelled(
            ambient_space,
            simplexes.into_iter().map(|spx| (spx, ())).collect(),
        )
    }
    fn new_unchecked(
        ambient_space: SP,
        simplexes: HashSet<Simplex<FS, SP>>,
    ) -> Self::WithLabel<()> {
        Self::WithLabel::<()>::new_labelled_unchecked(
            ambient_space,
            simplexes.into_iter().map(|spx| (spx, ())).collect(),
        )
    }

    fn new_labelled(
        ambient_space: SP,
        simplexes: HashMap<Simplex<FS, SP>, T>,
    ) -> Result<Self, &'static str>;
    fn new_labelled_unchecked(ambient_space: SP, simplexes: HashMap<Simplex<FS, SP>, T>) -> Self;

    fn ambient_space(&self) -> SP;

    fn simplexes<'a>(&'a self) -> HashSet<&'a Simplex<FS, SP>>
    where
        T: 'a,
    {
        self.labelled_simplexes().into_keys().collect()
    }
    fn into_simplexes(self) -> HashSet<Simplex<FS, SP>> {
        self.into_labelled_simplexes().into_keys().collect()
    }

    fn labelled_simplexes(&self) -> HashMap<&Simplex<FS, SP>, &T>;
    fn into_labelled_simplexes(self) -> HashMap<Simplex<FS, SP>, T>;

    fn subset_by_label(
        &self,
        label: &T,
    ) -> <Self::SubsetType as LabelledSimplexCollection<FS, SP, T>>::WithLabel<()> {
        self.subset_by_filter(|spx_label| spx_label == label)
            .forget_labels()
    }
    fn into_subset_by_label(
        self,
        label: &T,
    ) -> <Self::SubsetType as LabelledSimplexCollection<FS, SP, T>>::WithLabel<()> {
        self.into_subset_by_filter(|spx_label| spx_label == label)
            .forget_labels()
    }
    fn subset_by_filter(&self, f: impl Fn(&T) -> bool) -> Self::SubsetType {
        Self::SubsetType::new_labelled_unchecked(
            self.ambient_space(),
            self.labelled_simplexes()
                .into_iter()
                .filter(|(_spx, label)| f(label))
                .map(|(spx, label)| (spx.clone(), label.clone()))
                .collect(),
        )
    }
    fn into_subset_by_filter(self, f: impl Fn(&T) -> bool) -> Self::SubsetType {
        Self::SubsetType::new_labelled_unchecked(
            self.ambient_space(),
            self.into_labelled_simplexes()
                .into_iter()
                .filter(|(_spx, label)| f(label))
                .collect(),
        )
    }

    fn into_partial_simplicial_complex(self) -> LabelledPartialSimplicialComplex<FS, SP, T>;

    fn apply_label_function<S: Eq + Clone>(&self, f: impl Fn(&T) -> S) -> Self::WithLabel<S> {
        LabelledSimplexCollection::new_labelled_unchecked(
            self.ambient_space(),
            self.labelled_simplexes()
                .into_iter()
                .map(|(spx, label)| (spx.clone(), f(label)))
                .collect(),
        )
    }
    fn into_apply_label_function<S: Eq + Clone>(self, f: impl Fn(T) -> S) -> Self::WithLabel<S> {
        LabelledSimplexCollection::new_labelled_unchecked(
            self.ambient_space(),
            self.into_labelled_simplexes()
                .into_iter()
                .map(|(spx, label)| (spx, f(label)))
                .collect(),
        )
    }
    fn forget_labels(&self) -> Self::WithLabel<()> {
        self.apply_label_function(|_| ())
    }
    fn into_forget_labels(self) -> Self::WithLabel<()> {
        self.into_apply_label_function(|_| ())
    }

    fn common_label<'a>(
        &'a self,
        simplexes: impl Iterator<Item = &'a Simplex<FS, SP>>,
    ) -> Option<&'a T>
    where
        FS: 'a,
        SP: 'a,
    {
        let mut label = None;
        for spx in simplexes {
            let spx_label = *self.labelled_simplexes().get(&spx).unwrap();
            match label {
                Some(label) => {
                    if label != spx_label {
                        return None;
                    }
                }
                None => {
                    label = Some(spx_label);
                }
            }
        }
        label
    }
}
