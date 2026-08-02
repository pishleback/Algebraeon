use super::*;
use crate::simplex_overlap::simplex_interior_overlap;
use crate::{
    ambient_space::common_space,
    convex_hull::ConvexHull,
    partial_simplicial_complex::{LabelledPartialSimplicialComplex, PartialSimplicialComplex},
    simplex::Simplex,
    simplex_collection::{InteriorOrBoundarySimplexCollection, LabelledSimplexCollection},
    simplicial_complex::{LabelledSimplicialComplex, SimplicialComplex},
    simplicial_disjoint_union::{LabelledSimplicialDisjointUnion, SimplicialDisjointUnion},
};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum VennLabel {
    Left,
    Middle,
    Right,
}

fn simplex_venn< FS: OrderedRingSignature + FieldSignature>(
    left_simplex: &Simplex< FS>,
    right_simplex: &Simplex< FS>,
) -> LabelledPartialSimplicialComplex< FS, VennLabel>
where
    FS::Elem: Hash,
{
    let ambient_space =
        common_space(left_simplex.ambient_space(), right_simplex.ambient_space()).unwrap();

    // optimization
    if !simplex_interior_overlap(left_simplex, right_simplex) {
        return LabelledPartialSimplicialComplex::< FS, VennLabel>::new_labelled_unchecked(
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
        return LabelledPartialSimplicialComplex::< FS, VennLabel>::new_labelled_unchecked(
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
    LabelledPartialSimplicialComplex::< FS, VennLabel>::new_labelled_unchecked(
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

impl< FS: OrderedRingSignature + FieldSignature, T: Eq + Clone + Send + Sync>
    LabelledSimplicialDisjointUnion< FS, T>
where
    FS::Elem: Hash,
{
    pub(crate) fn subtract_raw<S: Eq + Clone + Send + Sync>(
        &self,
        other: &LabelledSimplicialDisjointUnion< FS, S>,
    ) -> LabelledSimplicialDisjointUnion< FS, T> {
        let ambient_space = common_space(self.ambient_space(), other.ambient_space()).unwrap();

        Self::new_labelled_unchecked(
            ambient_space,
            self.labelled_simplexes()
                .into_iter()
                .collect::<Vec<_>>()
                .into_par_iter()
                .map(|(self_spx, self_spx_label)| {
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
                    self_leftover
                        .into_iter()
                        .map(|spx| (spx, self_spx_label.clone()))
                        .collect::<Vec<_>>()
                })
                .flatten()
                .collect(),
        )
    }

    pub(crate) fn intersect_raw<S: Eq + Clone + Send + Sync>(
        &self,
        other: &LabelledSimplicialDisjointUnion< FS, S>,
    ) -> LabelledSimplicialDisjointUnion< FS, (T, S)> {
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

    pub(crate) fn union_raw(&self, other: &Self) -> SimplicialDisjointUnion< FS> {
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

impl<
    
    FS: OrderedRingSignature + FieldSignature,
    T: Eq + Clone + Send + Sync,
    S: Eq + Clone + Send + Sync,
> Difference<LabelledSimplicialDisjointUnion< FS, S>>
    for LabelledSimplicialDisjointUnion< FS, T>
where
    FS::Elem: Hash,
{
    type Output = LabelledPartialSimplicialComplex< FS, T>;

    fn difference(&self, other: &LabelledSimplicialDisjointUnion< FS, S>) -> Self::Output {
        self.subtract_raw(other)
            .refine_into_partial_simplicial_complex()
            .simplify()
    }
}

impl<
    
    FS: OrderedRingSignature + FieldSignature,
    T: Eq + Clone + Send + Sync,
    S: Eq + Clone + Send + Sync,
> Difference<LabelledPartialSimplicialComplex< FS, S>>
    for LabelledSimplicialDisjointUnion< FS, T>
where
    FS::Elem: Hash,
{
    type Output = LabelledPartialSimplicialComplex< FS, T>;

    fn difference(&self, other: &LabelledPartialSimplicialComplex< FS, S>) -> Self::Output {
        self.difference(&other.clone().into_simplicial_disjoint_union())
    }
}

impl<
    
    FS: OrderedRingSignature + FieldSignature,
    T: Eq + Clone + Send + Sync,
    S: Eq + Clone + Send + Sync,
> Difference<LabelledSimplicialDisjointUnion< FS, S>>
    for LabelledPartialSimplicialComplex< FS, T>
where
    FS::Elem: Hash,
{
    type Output = LabelledPartialSimplicialComplex< FS, T>;

    fn difference(&self, other: &LabelledSimplicialDisjointUnion< FS, S>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .difference(other)
    }
}

impl<
    
    FS: OrderedRingSignature + FieldSignature,
    T: Eq + Clone + Send + Sync,
    S: Eq + Clone + Send + Sync,
> Difference<LabelledPartialSimplicialComplex< FS, S>>
    for LabelledPartialSimplicialComplex< FS, T>
where
    FS::Elem: Hash,
{
    type Output = LabelledPartialSimplicialComplex< FS, T>;

    fn difference(&self, other: &LabelledPartialSimplicialComplex< FS, S>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .difference(&other.clone().into_simplicial_disjoint_union())
    }
}

impl<
    
    FS: OrderedRingSignature + FieldSignature,
    T: Eq + Clone + Send + Sync,
    S: Eq + Clone + Send + Sync,
> Difference<LabelledSimplicialComplex< FS, S>> for LabelledSimplicialDisjointUnion< FS, T>
where
    FS::Elem: Hash,
{
    type Output = LabelledPartialSimplicialComplex< FS, T>;

    fn difference(&self, other: &LabelledSimplicialComplex< FS, S>) -> Self::Output {
        self.difference(&other.clone().into_simplicial_disjoint_union())
    }
}

impl<
    
    FS: OrderedRingSignature + FieldSignature,
    T: Eq + Clone + Send + Sync,
    S: Eq + Clone + Send + Sync,
> Difference<LabelledSimplicialComplex< FS, S>> for LabelledPartialSimplicialComplex< FS, T>
where
    FS::Elem: Hash,
{
    type Output = LabelledPartialSimplicialComplex< FS, T>;

    fn difference(&self, other: &LabelledSimplicialComplex< FS, S>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .difference(&other.clone().into_simplicial_disjoint_union())
    }
}

impl<
    
    FS: OrderedRingSignature + FieldSignature,
    T: Eq + Clone + Send + Sync,
    S: Eq + Clone + Send + Sync,
> Difference<LabelledSimplicialDisjointUnion< FS, S>> for LabelledSimplicialComplex< FS, T>
where
    FS::Elem: Hash,
{
    type Output = LabelledPartialSimplicialComplex< FS, T>;

    fn difference(&self, other: &LabelledSimplicialDisjointUnion< FS, S>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .difference(other)
    }
}

impl<
    
    FS: OrderedRingSignature + FieldSignature,
    T: Eq + Clone + Send + Sync,
    S: Eq + Clone + Send + Sync,
> Difference<LabelledPartialSimplicialComplex< FS, S>> for LabelledSimplicialComplex< FS, T>
where
    FS::Elem: Hash,
{
    type Output = LabelledPartialSimplicialComplex< FS, T>;

    fn difference(&self, other: &LabelledPartialSimplicialComplex< FS, S>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .difference(&other.clone().into_simplicial_disjoint_union())
    }
}

impl<
    
    FS: OrderedRingSignature + FieldSignature,
    T: Eq + Clone + Send + Sync,
    S: Eq + Clone + Send + Sync,
> Difference<LabelledSimplicialComplex< FS, S>> for LabelledSimplicialComplex< FS, T>
where
    FS::Elem: Hash,
{
    type Output = LabelledPartialSimplicialComplex< FS, T>;

    fn difference(&self, other: &LabelledSimplicialComplex< FS, S>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .difference(&other.clone().into_simplicial_disjoint_union())
    }
}

impl< FS: OrderedRingSignature + FieldSignature> Intersect<SimplicialDisjointUnion< FS>>
    for SimplicialDisjointUnion< FS>
where
    FS::Elem: Hash,
{
    type Output = PartialSimplicialComplex< FS>;

    fn intersect(&self, other: &SimplicialDisjointUnion< FS>) -> Self::Output {
        self.intersect_raw(other)
            .forget_labels()
            .refine_into_partial_simplicial_complex()
            .simplify()
    }
}

impl< FS: OrderedRingSignature + FieldSignature> Intersect<PartialSimplicialComplex< FS>>
    for SimplicialDisjointUnion< FS>
where
    FS::Elem: Hash,
{
    type Output = PartialSimplicialComplex< FS>;

    fn intersect(&self, other: &PartialSimplicialComplex< FS>) -> Self::Output {
        self.intersect(&other.clone().into_simplicial_disjoint_union())
    }
}

impl< FS: OrderedRingSignature + FieldSignature> Intersect<SimplicialDisjointUnion< FS>>
    for PartialSimplicialComplex< FS>
where
    FS::Elem: Hash,
{
    type Output = PartialSimplicialComplex< FS>;

    fn intersect(&self, other: &SimplicialDisjointUnion< FS>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .intersect(other)
    }
}

impl< FS: OrderedRingSignature + FieldSignature> Intersect<PartialSimplicialComplex< FS>>
    for PartialSimplicialComplex< FS>
where
    FS::Elem: Hash,
{
    type Output = PartialSimplicialComplex< FS>;

    fn intersect(&self, other: &PartialSimplicialComplex< FS>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .intersect(&other.clone().into_simplicial_disjoint_union())
    }
}

impl< FS: OrderedRingSignature + FieldSignature> Intersect<SimplicialDisjointUnion< FS>>
    for SimplicialComplex< FS>
where
    FS::Elem: Hash,
{
    type Output = PartialSimplicialComplex< FS>;

    fn intersect(&self, other: &SimplicialDisjointUnion< FS>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .intersect(other)
    }
}

impl< FS: OrderedRingSignature + FieldSignature> Intersect<PartialSimplicialComplex< FS>>
    for SimplicialComplex< FS>
where
    FS::Elem: Hash,
{
    type Output = PartialSimplicialComplex< FS>;

    fn intersect(&self, other: &PartialSimplicialComplex< FS>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .intersect(&other.clone().into_simplicial_disjoint_union())
    }
}

impl< FS: OrderedRingSignature + FieldSignature> Intersect<SimplicialComplex< FS>>
    for SimplicialDisjointUnion< FS>
where
    FS::Elem: Hash,
{
    type Output = PartialSimplicialComplex< FS>;

    fn intersect(&self, other: &SimplicialComplex< FS>) -> Self::Output {
        self.intersect(&other.clone().into_simplicial_disjoint_union())
    }
}

impl< FS: OrderedRingSignature + FieldSignature> Intersect<SimplicialComplex< FS>>
    for PartialSimplicialComplex< FS>
where
    FS::Elem: Hash,
{
    type Output = PartialSimplicialComplex< FS>;

    fn intersect(&self, other: &SimplicialComplex< FS>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .intersect(&other.clone().into_simplicial_disjoint_union())
    }
}

impl< FS: OrderedRingSignature + FieldSignature> Intersect<SimplicialComplex< FS>>
    for SimplicialComplex< FS>
where
    FS::Elem: Hash,
{
    type Output = SimplicialComplex< FS>;

    fn intersect(&self, other: &SimplicialComplex< FS>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .intersect(&other.clone().into_simplicial_disjoint_union())
            .try_into_simplicial_complex()
            .unwrap()
    }
}

impl< FS: OrderedRingSignature + FieldSignature> Union<SimplicialDisjointUnion< FS>>
    for SimplicialDisjointUnion< FS>
where
    FS::Elem: Hash,
{
    type Output = PartialSimplicialComplex< FS>;

    fn union(&self, other: &SimplicialDisjointUnion< FS>) -> Self::Output {
        self.union_raw(other)
            .refine_into_partial_simplicial_complex()
            .simplify()
    }
}

impl< FS: OrderedRingSignature + FieldSignature> Union<PartialSimplicialComplex< FS>>
    for SimplicialDisjointUnion< FS>
where
    FS::Elem: Hash,
{
    type Output = PartialSimplicialComplex< FS>;

    fn union(&self, other: &PartialSimplicialComplex< FS>) -> Self::Output {
        self.union(&other.clone().into_simplicial_disjoint_union())
    }
}

impl< FS: OrderedRingSignature + FieldSignature> Union<SimplicialDisjointUnion< FS>>
    for PartialSimplicialComplex< FS>
where
    FS::Elem: Hash,
{
    type Output = PartialSimplicialComplex< FS>;

    fn union(&self, other: &SimplicialDisjointUnion< FS>) -> Self::Output {
        self.clone().into_simplicial_disjoint_union().union(other)
    }
}

impl< FS: OrderedRingSignature + FieldSignature> Union<PartialSimplicialComplex< FS>>
    for PartialSimplicialComplex< FS>
where
    FS::Elem: Hash,
{
    type Output = PartialSimplicialComplex< FS>;

    fn union(&self, other: &PartialSimplicialComplex< FS>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .union(&other.clone().into_simplicial_disjoint_union())
    }
}

impl< FS: OrderedRingSignature + FieldSignature> Union<SimplicialDisjointUnion< FS>>
    for SimplicialComplex< FS>
where
    FS::Elem: Hash,
{
    type Output = PartialSimplicialComplex< FS>;

    fn union(&self, other: &SimplicialDisjointUnion< FS>) -> Self::Output {
        self.clone().into_simplicial_disjoint_union().union(other)
    }
}

impl< FS: OrderedRingSignature + FieldSignature> Union<PartialSimplicialComplex< FS>>
    for SimplicialComplex< FS>
where
    FS::Elem: Hash,
{
    type Output = PartialSimplicialComplex< FS>;

    fn union(&self, other: &PartialSimplicialComplex< FS>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .union(&other.clone().into_simplicial_disjoint_union())
    }
}

impl< FS: OrderedRingSignature + FieldSignature> Union<SimplicialComplex< FS>>
    for SimplicialDisjointUnion< FS>
where
    FS::Elem: Hash,
{
    type Output = PartialSimplicialComplex< FS>;

    fn union(&self, other: &SimplicialComplex< FS>) -> Self::Output {
        self.union(&other.clone().into_simplicial_disjoint_union())
    }
}

impl< FS: OrderedRingSignature + FieldSignature> Union<SimplicialComplex< FS>>
    for PartialSimplicialComplex< FS>
where
    FS::Elem: Hash,
{
    type Output = PartialSimplicialComplex< FS>;

    fn union(&self, other: &SimplicialComplex< FS>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .union(&other.clone().into_simplicial_disjoint_union())
    }
}

impl< FS: OrderedRingSignature + FieldSignature> Union<SimplicialComplex< FS>>
    for SimplicialComplex< FS>
where
    FS::Elem: Hash,
{
    type Output = SimplicialComplex< FS>;

    fn union(&self, other: &SimplicialComplex< FS>) -> Self::Output {
        self.clone()
            .into_simplicial_disjoint_union()
            .union(&other.clone().into_simplicial_disjoint_union())
            .try_into_simplicial_complex()
            .unwrap()
    }
}

// impl< FS: OrderedRingSignature + FieldSignature> SimplicialComplex< FS>
// where
//     FS::Elem: Hash,
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
