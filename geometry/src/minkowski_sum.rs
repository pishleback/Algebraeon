use crate::{
    ambient_space::common_space,
    partial_simplicial_complex::PartialSimplicialComplex,
    simplex_collection::{InteriorOrBoundarySimplexCollection, LabelledSimplexCollection},
    simplicial_disjoint_union::SimplicialDisjointUnion,
};
use algebraeon_rings::structure::{FieldSignature, OrderedRingSignature};
use itertools::Itertools;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::hash::Hash;

pub trait MinkowskiSumRaw<Other> {
    type Output;
    fn minkowski_sum_raw(&self, other: &Other) -> Self::Output;
}

pub trait MinkowskiSum<Other> {
    type Output;
    fn minkowski_sum(&self, other: &Other) -> Self::Output;
}

impl<'f, FS: OrderedRingSignature + FieldSignature>
    MinkowskiSumRaw<PartialSimplicialComplex<'f, FS>> for PartialSimplicialComplex<'f, FS>
where
    FS::Set: Hash,
{
    type Output = SimplicialDisjointUnion<'f, FS>;

    fn minkowski_sum_raw(&self, other: &PartialSimplicialComplex<'f, FS>) -> Self::Output {
        let space = common_space(self.ambient_space(), other.ambient_space()).unwrap();

        self.simplexes()
            .into_iter()
            .cartesian_product(other.simplexes().into_iter().collect::<Vec<_>>())
            .collect::<Vec<_>>()
            .into_par_iter()
            .map(|(self_spx, other_spx)| {
                let mut points = vec![];
                for p in self_spx.points() {
                    for q in other_spx.points() {
                        points.push(p + q);
                    }
                }
                space
                    .convex_hull(points)
                    .to_simplicial_complex()
                    .interior()
                    .into_simplicial_disjoint_union()
            })
            .reduce(
                || space.empty_subset().into_simplicial_disjoint_union(),
                |left, right| {
                    println!(
                        "{:?} + {:?}",
                        left.simplexes().len(),
                        right.simplexes().len()
                    );
                    let shape = left.union_raw(&right);
                    println!(" = {:?}", shape.simplexes().len());
                    shape
                },
            )
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> MinkowskiSum<PartialSimplicialComplex<'f, FS>>
    for PartialSimplicialComplex<'f, FS>
where
    FS::Set: Hash,
{
    type Output = PartialSimplicialComplex<'f, FS>;

    fn minkowski_sum(&self, other: &PartialSimplicialComplex<'f, FS>) -> Self::Output {
        self.minkowski_sum_raw(other)
            .refine_into_partial_simplicial_complex()
            .simplify()
    }
}
