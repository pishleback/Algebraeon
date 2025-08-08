use crate::{
    ambient_space::common_space,
    boolean_operations::Union,
    partial_simplicial_complex::PartialSimplicialComplex,
    simplex_collection::{InteriorOrBoundarySimplexCollection, LabelledSimplexCollection},
};
use algebraeon_rings::structure::{FieldSignature, OrderedRingSignature};
use std::hash::Hash;

pub trait MinkowskiSum<Other> {
    type Output;
    fn minkowski_sum(&self, other: &Other) -> Self::Output;
}

impl<'f, FS: OrderedRingSignature + FieldSignature> MinkowskiSum<PartialSimplicialComplex<'f, FS>>
    for PartialSimplicialComplex<'f, FS>
where
    FS::Set: Hash,
{
    type Output = PartialSimplicialComplex<'f, FS>;

    fn minkowski_sum(&self, other: &PartialSimplicialComplex<'f, FS>) -> Self::Output {
        let space = common_space(self.ambient_space(), other.ambient_space()).unwrap();
        let mut sum = space.empty_subset().into_partial_simplicial_complex();
        for self_spx in self.simplexes() {
            for other_spx in other.simplexes() {
                let mut points = vec![];
                for p in self_spx.points() {
                    for q in other_spx.points() {
                        points.push(p + q);
                    }
                }
                sum = sum.union(&space.convex_hull(points).to_simplicial_complex().interior());
            }
        }
        sum
    }
}
