use super::*;

// It is helpful for computational reasons to put an ordering on the vectors
// so that the points of a simplex can be ordered
impl<FS: OrderedRingSignature + FieldSignature, SP: Borrow<AffineSpace<FS>> + Clone> PartialOrd
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
impl<FS: OrderedRingSignature + FieldSignature, SP: Borrow<AffineSpace<FS>> + Clone> Ord
    for Vector<FS, SP>
{
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match self.partial_cmp(other) {
            Some(ans) => ans,
            None => panic!(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum InteriorBoundaryLabel {
    Interior,
    Boundary,
}

mod simplex;
pub use simplex::*;

mod convex_hull;
pub use convex_hull::*;

mod simplex_collection;
pub use simplex_collection::*;

mod simplicial_complex;
pub use simplicial_complex::*;

mod partial_simplicial_complex;
pub use partial_simplicial_complex::*;

mod simplicial_disjoint_union;
pub use simplicial_disjoint_union::*;

mod boolean_operations;
// pub use boolean_operations::*;
