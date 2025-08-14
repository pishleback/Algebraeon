use super::*;
use crate::simplex_overlap::{SimplexOverlapResult, simplex_overlap};
use crate::{
    ambient_space::common_space,
    convex_hull::ConvexHull,
    partial_simplicial_complex::{LabelledPartialSimplicialComplex, PartialSimplicialComplex},
    simplex::Simplex,
    simplex_collection::{InteriorOrBoundarySimplexCollection, LabelledSimplexCollection},
    simplicial_complex::{LabelledSimplicialComplex, SimplicialComplex},
    simplicial_disjoint_union::{LabelledSimplicialDisjointUnion, SimplicialDisjointUnion},
};
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum VennLabel {
    Left,
    Middle,
    Right,
}

fn simplex_venn<'f, FS: OrderedRingSignature + FieldSignature>(
    left_simplex: &Simplex<'f, FS>,
    right_simplex: &Simplex<'f, FS>,
) -> LabelledPartialSimplicialComplex<'f, FS, VennLabel>
where
    FS::Set: Hash,
{
    let ambient_space =
        common_space(left_simplex.ambient_space(), right_simplex.ambient_space()).unwrap();

    // optimization
    if simplex_overlap(left_simplex, right_simplex) != SimplexOverlapResult::Overlap {
        return LabelledPartialSimplicialComplex::<'f, FS, VennLabel>::new_labelled_unchecked(
            ambient_space,
            HashMap::from([
                (left_simplex.clone(), VennLabel::Left),
                (right_simplex.clone(), VennLabel::Right),
            ]),
        );
    }

    let overlap = ConvexHull::intersect(
        &ConvexHull::from_simplex(left_simplex.clone()),
        &ConvexHull::from_simplex(right_simplex.clone()),
    );

    // optimization
    if overlap
        .to_simplicial_complex()
        .interior()
        .simplexes()
        .is_empty()
    {
        return LabelledPartialSimplicialComplex::<'f, FS, VennLabel>::new_labelled_unchecked(
            ambient_space,
            HashMap::from([
                (left_simplex.clone(), VennLabel::Left),
                (right_simplex.clone(), VennLabel::Right),
            ]),
        );
    }

    let mut self_ext = overlap.clone();
    for pt in left_simplex.points() {
        self_ext.extend_by_point(pt.clone());
    }
    let self_parts = self_ext.to_simplicial_complex().interior().into_simplexes();

    let mut other_ext = overlap.clone();
    for pt in right_simplex.points() {
        other_ext.extend_by_point(pt.clone());
    }
    let other_parts = other_ext
        .to_simplicial_complex()
        .interior()
        .into_simplexes();

    let all_parts = self_parts.union(&other_parts);
    LabelledPartialSimplicialComplex::<'f, FS, VennLabel>::new_labelled_unchecked(
        ambient_space,
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

impl<'f, FS: OrderedRingSignature + FieldSignature, T: Eq + Clone>
    LabelledSimplicialDisjointUnion<'f, FS, T>
where
    FS::Set: Hash,
{
    pub(crate) fn subtract_raw<S: Eq + Clone>(
        &self,
        other: &LabelledSimplicialDisjointUnion<'f, FS, S>,
    ) -> LabelledSimplicialDisjointUnion<'f, FS, T> {
        let ambient_space = common_space(self.ambient_space(), other.ambient_space()).unwrap();

        Self::new_labelled_unchecked(ambient_space, {
            let mut simplexes = HashMap::new();
            for (self_spx, self_spx_label) in self.labelled_simplexes() {
                let mut self_leftover = HashSet::from([self_spx.clone()]);
                for other_spx in other.simplexes() {
                    self_leftover = self_leftover
                        .into_iter()
                        .flat_map(|self_leftover_spx| {
                            simplex_venn(&self_leftover_spx, other_spx)
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

    pub(crate) fn intersect_raw<S: Eq + Clone>(
        &self,
        other: &LabelledSimplicialDisjointUnion<'f, FS, S>,
    ) -> LabelledSimplicialDisjointUnion<'f, FS, (T, S)> {
        let ambient_space = common_space(self.ambient_space(), other.ambient_space()).unwrap();
        LabelledSimplicialDisjointUnion::new_labelled_unchecked(ambient_space, {
            let mut simplexes = HashMap::new();
            for (self_spx, self_spx_label) in self.labelled_simplexes() {
                for (other_spx, other_spx_label) in other.labelled_simplexes() {
                    for spx in simplex_venn(self_spx, other_spx)
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

    pub(crate) fn union_raw(&self, other: &Self) -> SimplicialDisjointUnion<'f, FS> {
        let ambient_space = common_space(self.ambient_space(), other.ambient_space()).unwrap();
        let mut simplexes = HashSet::new();
        for spx in Self::subtract_raw(other, self).into_simplexes() {
            simplexes.insert(spx);
        }
        for spx in self.simplexes() {
            simplexes.insert(spx.clone());
        }
        Self::new_unchecked(ambient_space, simplexes)
    }
}

pub trait Difference<Other> {
    type Output;
    fn difference(&self, other: &Other) -> Self::Output;
}

pub trait Intersect<Other> {
    type Output;
    fn intersect(&self, other: &Other) -> Self::Output;
}

pub trait Union<Other> {
    type Output;
    fn union(&self, other: &Other) -> Self::Output;
}

impl<'f, FS: OrderedRingSignature + FieldSignature, T: Eq + Clone, S: Eq + Clone>
    Difference<LabelledSimplicialDisjointUnion<'f, FS, S>>
    for LabelledSimplicialDisjointUnion<'f, FS, T>
where
    FS::Set: Hash,
{
    type Output = LabelledPartialSimplicialComplex<'f, FS, T>;

    fn difference(&self, other: &LabelledSimplicialDisjointUnion<'f, FS, S>) -> Self::Output {
        self.subtract_raw(other)
            .refine_into_partial_simplicial_complex()
            .simplify()
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature, T: Eq + Clone, S: Eq + Clone>
    Difference<LabelledPartialSimplicialComplex<'f, FS, S>>
    for LabelledSimplicialDisjointUnion<'f, FS, T>
where
    FS::Set: Hash,
{
    type Output = LabelledPartialSimplicialComplex<'f, FS, T>;

    fn difference(&self, other: &LabelledPartialSimplicialComplex<'f, FS, S>) -> Self::Output {
        self.difference(&other.clone().into_simplicial_disjoint_union())
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature, T: Eq + Clone, S: Eq + Clone>
    Difference<LabelledSimplicialDisjointUnion<'f, FS, S>>
    for LabelledPartialSimplicialComplex<'f, FS, T>
where
    FS::Set: Hash,
{
    type Output = LabelledPartialSimplicialComplex<'f, FS, T>;

    fn difference(&self, other: &LabelledSimplicialDisjointUnion<'f, FS, S>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .difference(other)
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature, T: Eq + Clone, S: Eq + Clone>
    Difference<LabelledPartialSimplicialComplex<'f, FS, S>>
    for LabelledPartialSimplicialComplex<'f, FS, T>
where
    FS::Set: Hash,
{
    type Output = LabelledPartialSimplicialComplex<'f, FS, T>;

    fn difference(&self, other: &LabelledPartialSimplicialComplex<'f, FS, S>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .difference(&other.clone().into_simplicial_disjoint_union())
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature, T: Eq + Clone, S: Eq + Clone>
    Difference<LabelledSimplicialComplex<'f, FS, S>> for LabelledSimplicialDisjointUnion<'f, FS, T>
where
    FS::Set: Hash,
{
    type Output = LabelledPartialSimplicialComplex<'f, FS, T>;

    fn difference(&self, other: &LabelledSimplicialComplex<'f, FS, S>) -> Self::Output {
        self.difference(&other.clone().into_simplicial_disjoint_union())
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature, T: Eq + Clone, S: Eq + Clone>
    Difference<LabelledSimplicialComplex<'f, FS, S>> for LabelledPartialSimplicialComplex<'f, FS, T>
where
    FS::Set: Hash,
{
    type Output = LabelledPartialSimplicialComplex<'f, FS, T>;

    fn difference(&self, other: &LabelledSimplicialComplex<'f, FS, S>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .difference(&other.clone().into_simplicial_disjoint_union())
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature, T: Eq + Clone, S: Eq + Clone>
    Difference<LabelledSimplicialDisjointUnion<'f, FS, S>> for LabelledSimplicialComplex<'f, FS, T>
where
    FS::Set: Hash,
{
    type Output = LabelledPartialSimplicialComplex<'f, FS, T>;

    fn difference(&self, other: &LabelledSimplicialDisjointUnion<'f, FS, S>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .difference(other)
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature, T: Eq + Clone, S: Eq + Clone>
    Difference<LabelledPartialSimplicialComplex<'f, FS, S>> for LabelledSimplicialComplex<'f, FS, T>
where
    FS::Set: Hash,
{
    type Output = LabelledPartialSimplicialComplex<'f, FS, T>;

    fn difference(&self, other: &LabelledPartialSimplicialComplex<'f, FS, S>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .difference(&other.clone().into_simplicial_disjoint_union())
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature, T: Eq + Clone, S: Eq + Clone>
    Difference<LabelledSimplicialComplex<'f, FS, S>> for LabelledSimplicialComplex<'f, FS, T>
where
    FS::Set: Hash,
{
    type Output = LabelledPartialSimplicialComplex<'f, FS, T>;

    fn difference(&self, other: &LabelledSimplicialComplex<'f, FS, S>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .difference(&other.clone().into_simplicial_disjoint_union())
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> Intersect<SimplicialDisjointUnion<'f, FS>>
    for SimplicialDisjointUnion<'f, FS>
where
    FS::Set: Hash,
{
    type Output = PartialSimplicialComplex<'f, FS>;

    fn intersect(&self, other: &SimplicialDisjointUnion<'f, FS>) -> Self::Output {
        self.intersect_raw(other)
            .forget_labels()
            .refine_into_partial_simplicial_complex()
            .simplify()
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> Intersect<PartialSimplicialComplex<'f, FS>>
    for SimplicialDisjointUnion<'f, FS>
where
    FS::Set: Hash,
{
    type Output = PartialSimplicialComplex<'f, FS>;

    fn intersect(&self, other: &PartialSimplicialComplex<'f, FS>) -> Self::Output {
        self.intersect(&other.clone().into_simplicial_disjoint_union())
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> Intersect<SimplicialDisjointUnion<'f, FS>>
    for PartialSimplicialComplex<'f, FS>
where
    FS::Set: Hash,
{
    type Output = PartialSimplicialComplex<'f, FS>;

    fn intersect(&self, other: &SimplicialDisjointUnion<'f, FS>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .intersect(other)
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> Intersect<PartialSimplicialComplex<'f, FS>>
    for PartialSimplicialComplex<'f, FS>
where
    FS::Set: Hash,
{
    type Output = PartialSimplicialComplex<'f, FS>;

    fn intersect(&self, other: &PartialSimplicialComplex<'f, FS>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .intersect(&other.clone().into_simplicial_disjoint_union())
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> Intersect<SimplicialDisjointUnion<'f, FS>>
    for SimplicialComplex<'f, FS>
where
    FS::Set: Hash,
{
    type Output = PartialSimplicialComplex<'f, FS>;

    fn intersect(&self, other: &SimplicialDisjointUnion<'f, FS>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .intersect(other)
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> Intersect<PartialSimplicialComplex<'f, FS>>
    for SimplicialComplex<'f, FS>
where
    FS::Set: Hash,
{
    type Output = PartialSimplicialComplex<'f, FS>;

    fn intersect(&self, other: &PartialSimplicialComplex<'f, FS>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .intersect(&other.clone().into_simplicial_disjoint_union())
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> Intersect<SimplicialComplex<'f, FS>>
    for SimplicialDisjointUnion<'f, FS>
where
    FS::Set: Hash,
{
    type Output = PartialSimplicialComplex<'f, FS>;

    fn intersect(&self, other: &SimplicialComplex<'f, FS>) -> Self::Output {
        self.intersect(&other.clone().into_simplicial_disjoint_union())
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> Intersect<SimplicialComplex<'f, FS>>
    for PartialSimplicialComplex<'f, FS>
where
    FS::Set: Hash,
{
    type Output = PartialSimplicialComplex<'f, FS>;

    fn intersect(&self, other: &SimplicialComplex<'f, FS>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .intersect(&other.clone().into_simplicial_disjoint_union())
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> Intersect<SimplicialComplex<'f, FS>>
    for SimplicialComplex<'f, FS>
where
    FS::Set: Hash,
{
    type Output = SimplicialComplex<'f, FS>;

    fn intersect(&self, other: &SimplicialComplex<'f, FS>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .intersect(&other.clone().into_simplicial_disjoint_union())
            .try_into_simplicial_complex()
            .unwrap()
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> Union<SimplicialDisjointUnion<'f, FS>>
    for SimplicialDisjointUnion<'f, FS>
where
    FS::Set: Hash,
{
    type Output = PartialSimplicialComplex<'f, FS>;

    fn union(&self, other: &SimplicialDisjointUnion<'f, FS>) -> Self::Output {
        self.union_raw(other)
            .refine_into_partial_simplicial_complex()
            .simplify()
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> Union<PartialSimplicialComplex<'f, FS>>
    for SimplicialDisjointUnion<'f, FS>
where
    FS::Set: Hash,
{
    type Output = PartialSimplicialComplex<'f, FS>;

    fn union(&self, other: &PartialSimplicialComplex<'f, FS>) -> Self::Output {
        self.union(&other.clone().into_simplicial_disjoint_union())
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> Union<SimplicialDisjointUnion<'f, FS>>
    for PartialSimplicialComplex<'f, FS>
where
    FS::Set: Hash,
{
    type Output = PartialSimplicialComplex<'f, FS>;

    fn union(&self, other: &SimplicialDisjointUnion<'f, FS>) -> Self::Output {
        self.clone().into_simplicial_disjoint_union().union(other)
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> Union<PartialSimplicialComplex<'f, FS>>
    for PartialSimplicialComplex<'f, FS>
where
    FS::Set: Hash,
{
    type Output = PartialSimplicialComplex<'f, FS>;

    fn union(&self, other: &PartialSimplicialComplex<'f, FS>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .union(&other.clone().into_simplicial_disjoint_union())
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> Union<SimplicialDisjointUnion<'f, FS>>
    for SimplicialComplex<'f, FS>
where
    FS::Set: Hash,
{
    type Output = PartialSimplicialComplex<'f, FS>;

    fn union(&self, other: &SimplicialDisjointUnion<'f, FS>) -> Self::Output {
        self.clone().into_simplicial_disjoint_union().union(other)
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> Union<PartialSimplicialComplex<'f, FS>>
    for SimplicialComplex<'f, FS>
where
    FS::Set: Hash,
{
    type Output = PartialSimplicialComplex<'f, FS>;

    fn union(&self, other: &PartialSimplicialComplex<'f, FS>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .union(&other.clone().into_simplicial_disjoint_union())
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> Union<SimplicialComplex<'f, FS>>
    for SimplicialDisjointUnion<'f, FS>
where
    FS::Set: Hash,
{
    type Output = PartialSimplicialComplex<'f, FS>;

    fn union(&self, other: &SimplicialComplex<'f, FS>) -> Self::Output {
        self.union(&other.clone().into_simplicial_disjoint_union())
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> Union<SimplicialComplex<'f, FS>>
    for PartialSimplicialComplex<'f, FS>
where
    FS::Set: Hash,
{
    type Output = PartialSimplicialComplex<'f, FS>;

    fn union(&self, other: &SimplicialComplex<'f, FS>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .union(&other.clone().into_simplicial_disjoint_union())
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> Union<SimplicialComplex<'f, FS>>
    for SimplicialComplex<'f, FS>
where
    FS::Set: Hash,
{
    type Output = SimplicialComplex<'f, FS>;

    fn union(&self, other: &SimplicialComplex<'f, FS>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .union(&other.clone().into_simplicial_disjoint_union())
            .try_into_simplicial_complex()
            .unwrap()
    }
}

// impl<'f, FS: OrderedRingSignature + FieldSignature> SimplicialComplex<'f, FS>
// where
//     FS::Set: Hash,
// {
//     pub fn union_raw(&self, other: &Self) -> Self {
//         LabelledSimplicialDisjointUnion::union_raw(&self.into(), &other.into())
//             .refine_into_partial_simplicial_complex()
//             .try_into_simplicial_complex()
//             .unwrap()
//     }

//     pub fn union(&self, other: &Self) -> Self {
//         self.union_raw(other).simplify()
//     }

//     pub fn intersect_raw(&self, other: &Self) -> Self {
//         LabelledSimplicialDisjointUnion::intersect_raw(&self.into(), &other.into())
//             .refine_into_partial_simplicial_complex()
//             .into_forget_labels()
//             .try_into_simplicial_complex()
//             .unwrap()
//     }

//     pub fn intersect(&self, other: &Self) -> Self {
//         self.intersect_raw(other).simplify()
//     }
// }

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
