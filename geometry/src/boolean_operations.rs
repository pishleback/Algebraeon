use super::*;
use crate::{
    ambient_space::common_space,
    convex_hull::ConvexHull,
    partial_simplicial_complex::LabelledPartialSimplicialComplex,
    simplex::Simplex,
    simplex_collection::LabelledSimplexCollection,
    simplicial_complex::{InteriorBoundaryLabel, SimplicialComplex},
    simplicial_disjoint_union::{LabelledSimplicialDisjointUnion, SimplicialDisjointUnion},
};
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum VennLabel {
    Left,
    Middle,
    Right,
}

impl<'f, FS: OrderedRingSignature + FieldSignature> Simplex<'f, FS>
where
    FS::Set: Hash,
{
    pub fn venn(&self, other: &Self) -> LabelledPartialSimplicialComplex<'f, FS, VennLabel> {
        let ambient_space = common_space(self.ambient_space(), other.ambient_space()).unwrap();

        let overlap = ConvexHull::intersect(
            &ConvexHull::from_simplex(self.clone()),
            &ConvexHull::from_simplex(other.clone()),
        );

        let mut self_ext = overlap.clone();
        for pt in self.points() {
            self_ext.extend_by_point(pt.clone());
        }
        let self_parts = self_ext
            .as_simplicial_complex()
            .subset_by_label(&InteriorBoundaryLabel::Interior)
            .into_simplexes();

        let mut other_ext = overlap.clone();
        for pt in other.points() {
            other_ext.extend_by_point(pt.clone());
        }
        let other_parts = other_ext
            .as_simplicial_complex()
            .subset_by_label(&InteriorBoundaryLabel::Interior)
            .into_simplexes();

        let all_parts = self_parts.union(&other_parts);
        LabelledPartialSimplicialComplex::<'f, FS, VennLabel>::new_labelled_unchecked(
            ambient_space.clone(),
            all_parts
                .into_iter()
                .map(|spx| {
                    let label = match (self_parts.contains(spx), other_parts.contains(spx)) {
                        (true, false) => VennLabel::Left,
                        (true, true) => VennLabel::Middle,
                        (false, true) => VennLabel::Right,
                        (false, false) => {
                            unreachable!()
                        }
                    };
                    (spx.clone(), label)
                })
                .collect(),
        )
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature, T: Eq + Clone>
    LabelledSimplicialDisjointUnion<'f, FS, T>
where
    FS::Set: Hash,
{
    pub fn subtract_raw<S: Eq + Clone>(
        &self,
        other: &LabelledSimplicialDisjointUnion<'f, FS, S>,
    ) -> LabelledSimplicialDisjointUnion<'f, FS, T> {
        let ambient_space = common_space(self.ambient_space(), other.ambient_space()).unwrap();

        Self::new_labelled_unchecked(ambient_space.clone(), {
            let mut simplexes = HashMap::new();
            for (self_spx, self_spx_label) in self.labelled_simplexes() {
                let mut self_leftover = HashSet::from([self_spx.clone()]);
                for other_spx in other.simplexes() {
                    self_leftover = self_leftover
                        .into_iter()
                        .flat_map(|self_leftover_spx| {
                            Simplex::venn(&self_leftover_spx, other_spx)
                                .subset_by_label(&VennLabel::Left)
                                .into_simplexes()
                        })
                        .collect();
                }
                for spx in self_leftover {
                    simplexes.insert(spx, self_spx_label.clone());
                }
            }
            simplexes
        })
    }

    pub fn intersection_raw<S: Eq + Clone>(
        &self,
        other: &LabelledSimplicialDisjointUnion<'f, FS, S>,
    ) -> LabelledSimplicialDisjointUnion<'f, FS, (T, S)> {
        let ambient_space = common_space(self.ambient_space(), other.ambient_space()).unwrap();
        LabelledSimplicialDisjointUnion::new_labelled_unchecked(ambient_space.clone(), {
            let mut simplexes = HashMap::new();
            for (self_spx, self_spx_label) in self.labelled_simplexes() {
                for (other_spx, other_spx_label) in other.labelled_simplexes() {
                    for spx in Simplex::venn(self_spx, other_spx)
                        .subset_by_label(&VennLabel::Middle)
                        .into_simplexes()
                    {
                        simplexes.insert(spx, (self_spx_label.clone(), other_spx_label.clone()));
                    }
                }
            }
            simplexes
        })
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> SimplicialDisjointUnion<'f, FS>
where
    FS::Set: Hash,
{
    pub fn union_raw(&self, other: &Self) -> SimplicialDisjointUnion<'f, FS> {
        let ambient_space = common_space(self.ambient_space(), other.ambient_space()).unwrap();
        let mut simplexes = HashSet::new();
        for spx in Self::subtract_raw(other, self).into_simplexes() {
            simplexes.insert(spx);
        }
        for spx in self.simplexes() {
            simplexes.insert(spx.clone());
        }
        Self::new_unchecked(ambient_space.clone(), simplexes)
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> SimplicialComplex<'f, FS>
where
    FS::Set: Hash,
{
    pub fn union_raw(&self, other: &Self) -> Self {
        LabelledSimplicialDisjointUnion::union_raw(&self.into(), &other.into())
            .refine_to_partial_simplicial_complex()
            .try_as_simplicial_complex()
            .unwrap()
    }

    pub fn union(&self, other: &Self) -> Self {
        self.union_raw(other).simplify()
    }

    pub fn intersection_raw(&self, other: &Self) -> Self {
        LabelledSimplicialDisjointUnion::intersection_raw(&self.into(), &other.into())
            .refine_to_partial_simplicial_complex()
            .apply_label_function(|_| ())
            .try_as_simplicial_complex()
            .unwrap()
    }

    pub fn intersection(&self, other: &Self) -> Self {
        self.intersection_raw(other).simplify()
    }
}

/*
 - Venn dju <T1> and dju <T2> to produce dju <(Option<T1>, Option<T2>)>
 - Replace partial simplicial complex (psc) with labelled simplicial complex <bool>
 - Intersect psc, psc -> psc
 - Union psc, psc -> psc
 - Subtract psc, psc -> psc
 - Have a trait for a collection of labelled simplexes
   - Labelled subset
   - Filtered labelled subset
   - Union -> PartialSimplicialComplex
   - Intersection -> PartialSimplicialComplex
   - Difference -> PartialSimplicialComplex
 - Implement it for:
   - SimplexUnion
   - SimplexDisjointUnion
   - SemiSimplicialComplex
   - SimplicialComplex
*/
