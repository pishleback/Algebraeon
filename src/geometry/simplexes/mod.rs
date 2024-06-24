use super::*;

impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone> PartialOrd
    for Vector<FS, SP>
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        let space = common_space(
            self.ambient_space().borrow(),
            other.ambient_space().borrow(),
        )?;
        for i in 0..space.linear_dimension().unwrap() {
            match space
                .ordered_field()
                .ring_cmp(self.coordinate(i), other.coordinate(i))
            {
                std::cmp::Ordering::Less => {
                    return Some(std::cmp::Ordering::Less);
                }
                std::cmp::Ordering::Equal => {}
                std::cmp::Ordering::Greater => {
                    return Some(std::cmp::Ordering::Greater);
                }
            }
        }
        Some(std::cmp::Ordering::Equal)
    }
}

impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone> Ord
    for Vector<FS, SP>
{
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match self.partial_cmp(other) {
            Some(ans) => ans,
            None => panic!(),
        }
    }
}

mod simplex;
pub use simplex::*;

mod convex_hull;
pub use convex_hull::*;

mod simplicial_complex;
pub use simplicial_complex::*;

// mod partial_simplicial_complex {
//     use super::*;

//     #[derive(Debug, Clone)]
//     pub struct PartialSimplicialComplex {}

//     impl PartialSimplicialComplex {
//         // pub fn overlap_components(
//         //     a: &PartialSimplicialComplex,
//         //     b: &PartialSimplicialComplex,
//         // ) -> (
//         //     Rc<SimplicialComplex>,                       // a union b
//         //     SubSimplicialComplex<Rc<SimplicialComplex>>, // a intersect b
//         //     SubSimplicialComplex<Rc<SimplicialComplex>>, // a without b
//         //     SubSimplicialComplex<Rc<SimplicialComplex>>, // b without a
//         // ) {
//         //     todo!()
//         // }
//     }
// }
// pub use partial_simplicial_complex::*;

// fn split<
//     FS: OrderedRingStructure + FieldStructure,
//     SP: Borrow<LinearSpace<FS>> + Clone,
// >(
//     a: &Simplex<FS, SP>,
//     b: &OrientedSimplex<FS, SP>,
// ) -> (
//     Vec<Simplex<FS, SP>>,
//     Vec<Simplex<FS, SP>>,
//     Vec<Simplex<FS, SP>>,
// ) {
//     todo!()
// }

// fn intersect<
//     FS: OrderedRingStructure + FieldStructure,
//     SP: Borrow<LinearSpace<FS>> + Clone,
// >(
//     a: &Simplex<FS, SP>,
//     b: &Simplex<FS, SP>,
// ) -> Vec<Simplex<FS, SP>> {
//     todo!()
// }

/*
Plan

TODO: SC venn SC:
    For each simplex in one paired with a simplex in the other subdivide each and subdivide the SCs as induced such that we have a venn diagram

TODO: Refine SC given refinement of some simplex:
    Delete the interior of the simplex
    Get induced refinement of all facets, ridges, ...
    Replace facets, riges, ... with their refinement and induce refinement of adjacent simplexes in the SC
    Reinsert the refined interior of the simplex

TODO: Intersect a pair of simplexes A, B
    Express the subspace spanned by B as the intersection of hyperplanes
    Intersect A other with each hyperplane
    Restrict to the affine subspace spanned by B
    Intersect A with the interior side of each oriented facet of B
    Return to the ambient space
    Compute the convex hull as an SC

TODO: Intersect a simplex with an oriented hyperplane
    Compute the convex hull of any neutral points and any edges with one positive and one negative endpoint

TODO: Get the component of a simplex on the positive side of a hyperplane
    induct on the dimension of the simplex
    0d: the simplex is a point and it is trivial
    nd: take the interior of the convex hull of the positive component of each facet

*/
