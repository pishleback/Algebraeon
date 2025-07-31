use super::*;

// It is helpful for computational reasons to put an ordering on the vectors
// so that the points of a simplex can be ordered
#[allow(clippy::non_canonical_partial_ord_impl)]
impl<'f, FS: OrderedRingSignature + FieldSignature, SP: Borrow<AffineSpace<'f, FS>> + Clone>
    PartialOrd for Vector<'f, FS, SP>
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        let space = common_space(
            self.ambient_space().borrow(),
            other.ambient_space().borrow(),
        )?;
        for i in 0..space.linear_dimension().unwrap() {
            match space
                .field()
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
impl<'f, FS: OrderedRingSignature + FieldSignature, SP: Borrow<AffineSpace<'f, FS>> + Clone> Ord
    for Vector<'f, FS, SP>
{
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        if let Some(ans) = self.partial_cmp(other) {
            ans
        } else {
            panic!();
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
