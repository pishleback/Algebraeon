use crate::{ambient_space::common_space, simplex::Simplex};
use algebraeon_rings::structure::{FieldSignature, OrderedRingSignature};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SimplexOverlapResult {
    Disjoint, //the closures of the simplexes are disjoint
    Touching, //the closures of the simplexes overlap but the interiors are disjoint
    Overlap,  //the interiors of the simplexes overlap
}

/*
Calculate whether two simplexes in the same space overlap.

The idea is to use the Hyperplane Separation Theorem https://en.wikipedia.org/wiki/Hyperplane_separation_theorem:
    Two convex sets are disjoint if and only if there exists a hyperplane separating them

For a pair of simplexes, it suffices to check for separating hyperplanes parallel to the minkowski sum of the simplexes
The linear cosets defined by the hyperplanar faces of the minkowski sum are given by the linear coset spanned by the k-dimensional subsimplexes of A and the (n-k)-dimensional subsimplexes of B for k=1,2,...,n-1 where n+1 = the dimension of the minkowski sum

In the algorithm, we work with the normal vectors to hyperplanes rather than the hyperplanes themselves
*/
pub fn simplex_overlap<'f, FS: OrderedRingSignature + FieldSignature>(
    a: &Simplex<'f, FS>,
    b: &Simplex<'f, FS>,
) -> SimplexOverlapResult {
    let space = common_space(a.ambient_space(), b.ambient_space()).unwrap();

    if a.n() == 0 || b.n() == 0 {
        return SimplexOverlapResult::Disjoint;
    }
    debug_assert!(space.linear_dimension().is_some()); // non-empty space

    let a_root = a.point(0);
    let b_root = b.point(0);

    let a_vecs = (1..a.n()).map(|i| a.point(i) - a_root).collect::<Vec<_>>();
    let b_vecs = (1..b.n()).map(|i| b.point(i) - b_root).collect::<Vec<_>>();

    let linear_vecs_span = space.affine_subspace_from_root_and_linear_span(
        &space.origin().unwrap(),
        a_vecs.iter().chain(b_vecs.iter()).collect(),
    );

    // the dimension of the minkowski sum of a and b
    let s_dim = linear_vecs_span
        .embedded_space()
        .linear_dimension()
        .unwrap();

    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ambient_space::AffineSpace;
    use algebraeon_nzq::Rational;

    #[test]
    fn something_and_null() {
        let space3 = AffineSpace::new_linear(Rational::structure_ref(), 3);

        assert_eq!(
            simplex_overlap(
                &space3.simplex(vec![space3.vector([1, 2, 3])]).unwrap(),
                &space3.simplex(vec![]).unwrap()
            ),
            SimplexOverlapResult::Disjoint
        );

        assert_eq!(
            simplex_overlap(
                &space3
                    .simplex(vec![space3.vector([1, 2, 3]), space3.vector([4, 5, 6])])
                    .unwrap(),
                &space3.simplex(vec![]).unwrap()
            ),
            SimplexOverlapResult::Disjoint
        );
    }

    #[test]
    fn two_points() {
        let space3 = AffineSpace::new_linear(Rational::structure_ref(), 3);

        // different points
        assert_eq!(
            simplex_overlap(
                &space3.simplex(vec![space3.vector([1, 2, 3])]).unwrap(),
                &space3.simplex(vec![space3.vector([1, 3, 2])]).unwrap()
            ),
            SimplexOverlapResult::Disjoint
        );

        // same point
        assert_eq!(
            simplex_overlap(
                &space3.simplex(vec![space3.vector([1, 2, 3])]).unwrap(),
                &space3.simplex(vec![space3.vector([1, 2, 3])]).unwrap()
            ),
            SimplexOverlapResult::Overlap
        );
    }

    #[test]
    fn point_and_line() {
        let space3 = AffineSpace::new_linear(Rational::structure_ref(), 3);

        // point in the middle of the line
        assert_eq!(
            simplex_overlap(
                &space3
                    .simplex(vec![space3.vector([0, 0, 0]), space3.vector([2, 0, 0])])
                    .unwrap(),
                &space3.simplex(vec![space3.vector([1, 0, 0])]).unwrap()
            ),
            SimplexOverlapResult::Overlap
        );

        // point at the end of the line
        assert_eq!(
            simplex_overlap(
                &space3
                    .simplex(vec![space3.vector([0, 0, 0]), space3.vector([2, 0, 0])])
                    .unwrap(),
                &space3.simplex(vec![space3.vector([2, 0, 0])]).unwrap()
            ),
            SimplexOverlapResult::Touching
        );

        // point disjoint from line
        assert_eq!(
            simplex_overlap(
                &space3
                    .simplex(vec![space3.vector([0, 0, 0]), space3.vector([2, 0, 0])])
                    .unwrap(),
                &space3.simplex(vec![space3.vector([0, 1, 0])]).unwrap()
            ),
            SimplexOverlapResult::Disjoint
        );
        assert_eq!(
            simplex_overlap(
                &space3
                    .simplex(vec![space3.vector([0, 0, 0]), space3.vector([2, 0, 0])])
                    .unwrap(),
                &space3.simplex(vec![space3.vector([3, 0, 0])]).unwrap()
            ),
            SimplexOverlapResult::Disjoint
        );
        assert_eq!(
            simplex_overlap(
                &space3
                    .simplex(vec![space3.vector([0, 0, 0]), space3.vector([2, 0, 0])])
                    .unwrap(),
                &space3.simplex(vec![space3.vector([2, 1, 0])]).unwrap()
            ),
            SimplexOverlapResult::Disjoint
        );
        assert_eq!(
            simplex_overlap(
                &space3
                    .simplex(vec![space3.vector([0, 0, 0]), space3.vector([2, 0, 0])])
                    .unwrap(),
                &space3.simplex(vec![space3.vector([3, 1, 0])]).unwrap()
            ),
            SimplexOverlapResult::Disjoint
        );
    }

    #[test]
    fn point_and_triangle() {
        let space3 = AffineSpace::new_linear(Rational::structure_ref(), 3);

        let triangle = space3
            .simplex(vec![
                space3.vector([6, 0, 0]),
                space3.vector([0, 6, 0]),
                space3.vector([0, 0, 6]),
            ])
            .unwrap();

        assert_eq!(
            simplex_overlap(
                &triangle,
                &space3.simplex(vec![space3.vector([0, 0, 0])]).unwrap()
            ),
            SimplexOverlapResult::Disjoint
        );
        assert_eq!(
            simplex_overlap(
                &triangle,
                &space3.simplex(vec![space3.vector([4, 4, 4])]).unwrap()
            ),
            SimplexOverlapResult::Disjoint
        );

        assert_eq!(
            simplex_overlap(
                &triangle,
                &space3.simplex(vec![space3.vector([2, 2, 2])]).unwrap()
            ),
            SimplexOverlapResult::Overlap
        );

        assert_eq!(
            simplex_overlap(
                &triangle,
                &space3.simplex(vec![space3.vector([0, 6, 0])]).unwrap()
            ),
            SimplexOverlapResult::Touching
        );

        assert_eq!(
            simplex_overlap(
                &triangle,
                &space3.simplex(vec![space3.vector([3, 3, 0])]).unwrap()
            ),
            SimplexOverlapResult::Touching
        );
    }

    #[test]
    fn line_and_triangle() {
        let space3 = AffineSpace::new_linear(Rational::structure_ref(), 3);

        let triangle = space3
            .simplex(vec![
                space3.vector([6, 0, 0]),
                space3.vector([0, 6, 0]),
                space3.vector([0, 0, 6]),
            ])
            .unwrap();

        assert_eq!(
            simplex_overlap(
                &triangle,
                &space3
                    .simplex(vec![space3.vector([0, 0, 0]), space3.vector([1, 1, 1])])
                    .unwrap()
            ),
            SimplexOverlapResult::Disjoint
        );
        assert_eq!(
            simplex_overlap(
                &triangle,
                &space3
                    .simplex(vec![space3.vector([0, 0, 0]), space3.vector([2, 2, 2])])
                    .unwrap()
            ),
            SimplexOverlapResult::Touching
        );
        assert_eq!(
            simplex_overlap(
                &triangle,
                &space3
                    .simplex(vec![space3.vector([0, 0, 0]), space3.vector([3, 3, 3])])
                    .unwrap()
            ),
            SimplexOverlapResult::Overlap
        );

        assert_eq!(
            simplex_overlap(
                &triangle,
                &space3
                    .simplex(vec![space3.vector([6, 0, 0]), space3.vector([7, 0, 0])])
                    .unwrap()
            ),
            SimplexOverlapResult::Touching
        );

        assert_eq!(
            simplex_overlap(
                &triangle,
                &space3
                    .simplex(vec![space3.vector([0, -6, 0]), space3.vector([6, 0, 0])])
                    .unwrap()
            ),
            SimplexOverlapResult::Touching
        );

        assert_eq!(
            simplex_overlap(
                &triangle,
                &space3
                    .simplex(vec![space3.vector([6, -12, 0]), space3.vector([12, -6, 0])])
                    .unwrap()
            ),
            SimplexOverlapResult::Touching
        );
    }

    #[test]
    fn triangle_and_triangle() {
        let space3 = AffineSpace::new_linear(Rational::structure_ref(), 3);

        assert_eq!(
            simplex_overlap(
                &space3
                    .simplex(vec![
                        space3.vector([0, -1, 0]),
                        space3.vector([0, 1, 0]),
                        space3.vector([1, 0, 0])
                    ])
                    .unwrap(),
                &space3
                    .simplex(vec![
                        space3.vector([0, 0, -1]),
                        space3.vector([0, 0, 1]),
                        space3.vector([-1, 0, 0])
                    ])
                    .unwrap()
            ),
            SimplexOverlapResult::Touching
        );

        assert_eq!(
            simplex_overlap(
                &space3
                    .simplex(vec![
                        space3.vector([0, -1, 0]),
                        space3.vector([0, 1, 0]),
                        space3.vector([1, 0, 0])
                    ])
                    .unwrap(),
                &space3
                    .simplex(vec![
                        space3.vector([1, 0, -1]),
                        space3.vector([1, 0, 1]),
                        space3.vector([0, 0, 0])
                    ])
                    .unwrap()
            ),
            SimplexOverlapResult::Overlap
        );

        assert_eq!(
            simplex_overlap(
                &space3
                    .simplex(vec![
                        space3.vector([0, -1, 0]),
                        space3.vector([0, 1, 0]),
                        space3.vector([1, 0, 0])
                    ])
                    .unwrap(),
                &space3
                    .simplex(vec![
                        space3.vector([-1, 0, -1]),
                        space3.vector([-1, 0, 1]),
                        space3.vector([-2, 0, 0])
                    ])
                    .unwrap()
            ),
            SimplexOverlapResult::Disjoint
        );
    }
}
