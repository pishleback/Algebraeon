use crate::{
    ambient_space::{self, common_space},
    simplex::Simplex,
    vector::DotProduct,
};
use algebraeon_nzq::Rational;
use algebraeon_rings::{
    matrix::{Matrix, MatrixStructure},
    module::{
        finitely_free_module::FinitelyFreeModuleStructure,
        finitely_free_submodule::FinitelyFreeSubmoduleStructure,
    },
    structure::{FieldSignature, OrderedRingSignature},
};
use itertools::{Combinations, Itertools};

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
    let m = space.linear_dimension().unwrap();

    let field = space.field();

    let a_root = a.point(0);
    let b_root = b.point(0);

    let a_vecs = (1..a.n()).map(|i| a.point(i) - a_root).collect::<Vec<_>>();
    let b_vecs = (1..b.n()).map(|i| b.point(i) - b_root).collect::<Vec<_>>();

    println!("a_vecs = {:?}", a_vecs);
    println!("b_vecs = {:?}", b_vecs);

    let linear_vecs_span = space.affine_subspace_from_root_and_linear_span(
        &space.origin().unwrap(),
        a_vecs.iter().chain(b_vecs.iter()).collect(),
    );

    let linear_vecs_span_submodule = MatrixStructure::new(space.field().clone()).row_span(
        Matrix::construct(a_vecs.len() + b_vecs.len(), m, |r, c| {
            if r < a_vecs.len() {
                a_vecs[r].coordinate(c).clone()
            } else {
                b_vecs[r - a_vecs.len()].coordinate(c).clone()
            }
        }),
    );
    debug_assert_eq!(linear_vecs_span_submodule.module_rank(), m);

    // the dimension of the minkowski sum of A and B
    let n = linear_vecs_span
        .embedded_space()
        .linear_dimension()
        .unwrap();
    println!("n = {:?}", n);

    let mut just_touching = false;

    for i in 0..n {
        let j = n - i - 1;

        // An i-dim subsimplex of A and a j-dim subsimplex of B linearly span a hyperface of the minkowski sum of A and B
        for a_pts in (0..a.n()).combinations(i + 1) {
            for b_pts in (0..b.n()).combinations(j + 1) {
                let a_sub = a.sub_simplex(a_pts.clone());
                let b_sub = b.sub_simplex(b_pts);
                println!("a_sub = {:?}", a_sub);
                println!("b_sub = {:?}", b_sub);

                let a_sub_root = a_sub.point(0);
                let b_sub_root = b_sub.point(0);

                let a_sub_vecs = (0..i)
                    .map(|k| a_sub.point(k + 1) - a_sub_root)
                    .collect::<Vec<_>>();
                let b_sub_vecs = (0..j)
                    .map(|k| b_sub.point(k + 1) - b_sub_root)
                    .collect::<Vec<_>>();

                println!("a_sub_vecs = {:?}", a_sub_vecs);
                println!("b_sub_vecs = {:?}", b_sub_vecs);

                let normal_space = FinitelyFreeSubmoduleStructure::new(
                    FinitelyFreeModuleStructure::<FS, &'f FS>::new(field, m),
                )
                .intersect(
                    &MatrixStructure::new(space.field().clone()).col_kernel(Matrix::construct(
                        i + j,
                        m,
                        |r, c| {
                            if r < i {
                                a_sub_vecs[r].coordinate(c).clone()
                            } else {
                                b_sub_vecs[r - i].coordinate(c).clone()
                            }
                        },
                    )),
                    &linear_vecs_span_submodule,
                );

                println!("{:?}", normal_space.rank());
                debug_assert!(normal_space.rank() >= 1);

                if normal_space.rank() == 1 {
                    let normal =
                        space.vector(normal_space.basis().first().unwrap().iter().cloned());
                    println!("normal = {:?}", normal);

                    for vec in &a_sub_vecs {
                        debug_assert!(field.is_zero(&normal.dot(vec)));
                    }
                    for vec in &b_sub_vecs {
                        debug_assert!(field.is_zero(&normal.dot(vec)));
                    }

                    let mut a_points = a.points().into_iter();
                    let a_proj = normal.dot(a_points.next().unwrap());
                    let mut min_a_proj = a_proj.clone();
                    let mut max_a_proj = a_proj;
                    for a_pt in a_points {
                        let proj = normal.dot(a_pt);
                        if field.ring_cmp(&proj, &min_a_proj) == std::cmp::Ordering::Less {
                            min_a_proj = proj.clone();
                        }
                        if field.ring_cmp(&proj, &max_a_proj) == std::cmp::Ordering::Greater {
                            max_a_proj = proj.clone();
                        }
                    }

                    let mut b_points = b.points().into_iter();
                    let b_proj = normal.dot(b_points.next().unwrap());
                    let mut min_b_proj = b_proj.clone();
                    let mut max_b_proj = b_proj;
                    for b_pt in b_points {
                        let proj = normal.dot(b_pt);
                        if field.ring_cmp(&proj, &min_b_proj) == std::cmp::Ordering::Less {
                            min_b_proj = proj.clone();
                        }
                        if field.ring_cmp(&proj, &max_b_proj) == std::cmp::Ordering::Greater {
                            max_b_proj = proj.clone();
                        }
                    }

                    println!(
                        "({:?}, {:?}) ({:?}, {:?})",
                        min_a_proj, max_a_proj, min_b_proj, max_b_proj
                    );

                    match (
                        field.ring_cmp(&max_a_proj, &min_b_proj),
                        field.ring_cmp(&max_b_proj, &min_a_proj),
                    ) {
                        (std::cmp::Ordering::Less, _) | (_, std::cmp::Ordering::Less) => {
                            return SimplexOverlapResult::Disjoint;
                        }
                        (std::cmp::Ordering::Greater, std::cmp::Ordering::Greater) => {}
                        (std::cmp::Ordering::Equal, _) | (_, std::cmp::Ordering::Equal) => {
                            println!("TOUCH");
                            just_touching = true;
                        }
                    }
                }
            }
        }
    }

    match just_touching {
        false => SimplexOverlapResult::Overlap,
        true => SimplexOverlapResult::Touching,
    }
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
                &space3.simplex(vec![space3.vector([1, 1, 0])]).unwrap()
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
            SimplexOverlapResult::Disjoint
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
