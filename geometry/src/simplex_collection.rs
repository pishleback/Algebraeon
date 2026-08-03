use super::*;
use crate::{
    ambient_space::AffineSpace, partial_simplicial_complex::LabelledPartialSimplicialComplex,
    simplex::Simplex, simplicial_disjoint_union::LabelledSimplicialDisjointUnion,
};
use std::collections::{HashMap, HashSet};

/// A collection of disjoint simplices labelled by T
pub trait LabelledSimplexCollection<
    FS: OrderedRingSignature + FieldSignature,
    T: Eq + Clone + Send + Sync,
>: Sized where
    FS::Elem: Hash,
{
    type WithLabel<S: Eq + Clone + Send + Sync>: LabelledSimplexCollection<FS, S>;
    type SubsetType: LabelledSimplexCollection<FS, T>;

    fn try_new(
        ambient_space: &AffineSpace<FS>,
        simplexes: HashSet<Simplex<FS>>,
    ) -> Result<Self::WithLabel<()>, &'static str> {
        Self::WithLabel::<()>::try_new_labelled(
            ambient_space,
            simplexes.into_iter().map(|spx| (spx, ())).collect(),
        )
    }
    fn new_unchecked(
        ambient_space: &AffineSpace<FS>,
        simplexes: HashSet<Simplex<FS>>,
    ) -> Self::WithLabel<()> {
        Self::WithLabel::<()>::new_labelled_unchecked(
            ambient_space,
            simplexes.into_iter().map(|spx| (spx, ())).collect(),
        )
    }

    fn try_new_labelled(
        ambient_space: &AffineSpace<FS>,
        simplexes: HashMap<Simplex<FS>, T>,
    ) -> Result<Self, &'static str>;
    fn new_labelled_unchecked(
        ambient_space: &AffineSpace<FS>,
        simplexes: HashMap<Simplex<FS>, T>,
    ) -> Self;

    fn ambient_space(&self) -> &AffineSpace<FS>;

    fn simplexes<'a>(&'a self) -> HashSet<&'a Simplex<FS>>
    where
        T: 'a,
    {
        self.labelled_simplexes().into_keys().collect()
    }
    fn into_simplexes(self) -> HashSet<Simplex<FS>> {
        self.into_labelled_simplexes().into_keys().collect()
    }

    fn labelled_simplexes(&self) -> HashMap<&Simplex<FS>, &T>;
    fn into_labelled_simplexes(self) -> HashMap<Simplex<FS>, T>;

    fn subset_by_label(
        &self,
        label: &T,
    ) -> <Self::SubsetType as LabelledSimplexCollection<FS, T>>::WithLabel<()> {
        self.subset_by_filter(|spx_label| spx_label == label)
            .forget_labels()
    }
    fn into_subset_by_label(
        self,
        label: &T,
    ) -> <Self::SubsetType as LabelledSimplexCollection<FS, T>>::WithLabel<()> {
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
            &self.ambient_space().clone(),
            self.into_labelled_simplexes()
                .into_iter()
                .filter(|(_spx, label)| f(label))
                .collect(),
        )
    }

    fn into_partial_simplicial_complex(self) -> LabelledPartialSimplicialComplex<FS, T>;
    fn to_partial_simplicial_complex(&self) -> LabelledPartialSimplicialComplex<FS, T>;

    fn into_simplicial_disjoint_union(self) -> LabelledSimplicialDisjointUnion<FS, T>;
    fn to_simplicial_disjoint_union(&self) -> LabelledSimplicialDisjointUnion<FS, T>;

    fn apply_label_function<S: Eq + Clone + Send + Sync>(
        &self,
        f: impl Fn(&T) -> S,
    ) -> Self::WithLabel<S> {
        LabelledSimplexCollection::new_labelled_unchecked(
            self.ambient_space(),
            self.labelled_simplexes()
                .into_iter()
                .map(|(spx, label)| (spx.clone(), f(label)))
                .collect(),
        )
    }
    fn into_apply_label_function<S: Eq + Clone + Send + Sync>(
        self,
        f: impl Fn(T) -> S,
    ) -> Self::WithLabel<S> {
        LabelledSimplexCollection::new_labelled_unchecked(
            &self.ambient_space().clone(),
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

    fn common_label<'a>(&'a self, simplexes: impl Iterator<Item = &'a Simplex<FS>>) -> Option<&'a T>
    where
        FS: 'a,
        AffineSpace<FS>: 'a,
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum InteriorOrBoundary {
    Interior,
    Boundary,
}

pub trait InteriorOrBoundarySimplexCollection<FS: OrderedRingSignature + FieldSignature>:
    LabelledSimplexCollection<FS, InteriorOrBoundary>
where
    FS::Elem: Hash,
{
    fn interior(
        &self,
    ) -> <Self::SubsetType as LabelledSimplexCollection<FS, InteriorOrBoundary>>::WithLabel<()>
    {
        self.subset_by_label(&InteriorOrBoundary::Interior)
    }

    fn boundary(
        &self,
    ) -> <Self::SubsetType as LabelledSimplexCollection<FS, InteriorOrBoundary>>::WithLabel<()>
    {
        self.subset_by_label(&InteriorOrBoundary::Boundary)
    }
}

impl<
    FS: OrderedRingSignature + FieldSignature,
    S: LabelledSimplexCollection<FS, InteriorOrBoundary>,
> InteriorOrBoundarySimplexCollection<FS> for S
where
    FS::Elem: Hash,
{
}
