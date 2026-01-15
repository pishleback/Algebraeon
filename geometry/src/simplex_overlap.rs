use crate::{ambient_space::common_space, simplex::Simplex, vector::DotProduct};
use algebraeon_rings::{
    matrix::{Matrix, MatrixStructure},
    module::{
        finitely_free_module::FinitelyFreeModuleStructure,
        finitely_free_submodule::FinitelyFreeSubmoduleStructure,
    },
    structure::{FieldSignature, OrderedRingSignature, ZeroEqSignature},
};
use itertools::Itertools;
use std::collections::HashSet;
use std::hash::Hash;

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
// If RETURN_IF_TOUCHING=true then this function may return SimplexOverlapResult::Touching even if the simplexes are actually disjoint. That gives a performance boost if it only matters that the interiors are disjoint
fn simplex_overlap_impl<
    'f,
    FS: OrderedRingSignature + FieldSignature,
    const RETURN_IF_TOUCHING: bool,
>(
    a: &Simplex<'f, FS>,
    b: &Simplex<'f, FS>,
) -> SimplexOverlapResult
where
    FS::Set: Hash,
{
    let space = common_space(a.ambient_space(), b.ambient_space()).unwrap();

    // If either A or B is the null-simplex then they are disjoint.
    if a.n() == 0 || b.n() == 0 {
        return SimplexOverlapResult::Disjoint;
    }
    debug_assert!(space.linear_dimension().is_some()); // non-empty space
    let space_dim = space.linear_dimension().unwrap();
    let field = space.field();

    // Now since either A nor B is the null simplex we know they are both non-empty, so can pick a root point for each
    let a_root = a.point(0);
    let b_root = b.point(0);

    // For each point of A and B other than the root, take the vector from the root to that point
    let a_vecs = (1..a.n()).map(|i| a.point(i) - a_root).collect::<Vec<_>>();
    let b_vecs = (1..b.n()).map(|i| b.point(i) - b_root).collect::<Vec<_>>();

    // A vector starting in A and ending in B
    let a_to_b_vec = a_root - b_root;

    // This is the linear subspace obtained by rooting A and B at the origin and taking their combined linear span
    // The sum of A and B lives within a coset of this linear span
    let linear_span_foo = space.affine_subspace_from_root_and_linear_span(
        &space.origin().unwrap(),
        a_vecs.iter().chain(b_vecs.iter()).collect(),
    );
    let linear_span_foo_submodule = MatrixStructure::new(space.field().clone()).row_span(
        Matrix::construct(a_vecs.len() + b_vecs.len(), space_dim, |r, c| {
            if r < a_vecs.len() {
                a_vecs[r].coordinate(c).clone()
            } else {
                b_vecs[r - a_vecs.len()].coordinate(c).clone()
            }
        }),
    );
    debug_assert_eq!(linear_span_foo_submodule.module_rank(), space_dim);

    // This is the linear subspace obtained by rooting the affine span of A and B
    // In typical cases, it is equal to linear_span_foo
    //  In that case, we need only check the (finitely many) normals to the hyperplanes (in linear_span_foo) of the sum of A and B
    // In degenerate cases, it has dimension one greater than linear_span_foo
    //  In that case, in addition to the normals to the hyperplanes of the sum of A and B, we need to also check the normal vector to the hyperplane (in linear_span_bar) given by the sum of A and B
    let linear_span_bar = MatrixStructure::new(space.field().clone()).row_span(Matrix::construct(
        a_vecs.len() + b_vecs.len() + 1,
        space_dim,
        |r, c| {
            if r < a_vecs.len() {
                a_vecs[r].coordinate(c).clone()
            } else if r < a_vecs.len() + b_vecs.len() {
                b_vecs[r - a_vecs.len()].coordinate(c).clone()
            } else {
                a_to_b_vec.coordinate(c).clone()
            }
        },
    ));
    debug_assert_eq!(linear_span_bar.module_rank(), space_dim);
    debug_assert!(
        (linear_span_bar.rank() == linear_span_foo_submodule.rank())
            || (linear_span_bar.rank() == linear_span_foo_submodule.rank() + 1)
    );

    // Accumulate all the normals we need to check, excluding any duplicates
    let mut normals = HashSet::new();
    // Record whether A and B just touch when projected along any of the normals
    // If A and B don't overlap then they are touching iff they are touching when projected on some normal in our list
    let mut just_touching = false;

    let n = linear_span_foo.embedded_space().linear_dimension().unwrap();

    for normal in {
        if linear_span_bar.rank() != linear_span_foo_submodule.rank() {
            // Add the normal to the sum of A and B for the degenerate case

            let normal_space =
                FinitelyFreeSubmoduleStructure::new(
                    FinitelyFreeModuleStructure::<FS, &'f FS>::new(field, space_dim),
                )
                .intersect(
                    MatrixStructure::new(space.field().clone()).col_kernel(Matrix::construct(
                        a_vecs.len() + b_vecs.len(),
                        space_dim,
                        |r, c| {
                            if r < a_vecs.len() {
                                a_vecs[r].coordinate(c).clone()
                            } else {
                                b_vecs[r - a_vecs.len()].coordinate(c).clone()
                            }
                        },
                    )),
                    linear_span_bar.clone(),
                );
            debug_assert!(normal_space.rank() == 1);
            let normal = space.vector(normal_space.basis().first().unwrap().iter().cloned());
            for vec in &a_vecs {
                debug_assert!(field.is_zero(&normal.dot(vec)));
            }
            for vec in &b_vecs {
                debug_assert!(field.is_zero(&normal.dot(vec)));
            }
            Some(normal)
        } else {
            None
        }
    }
    .into_iter()
    // Add the normals comming from the faces of the sum of A and B
    .chain((0..n).map(|i| (i, n - i - 1)).flat_map(|(i, j)| {
        (0..a.n())
            .combinations(i + 1)
            .cartesian_product((0..b.n()).combinations(j + 1))
            .filter_map(|(a_pts, b_pts)| {
                let i = a_pts.len() - 1;
                let j = b_pts.len() - 1;

                let a_sub = a.sub_simplex(a_pts.clone());
                let b_sub = b.sub_simplex(b_pts);

                let a_sub_root = a_sub.point(0);
                let b_sub_root = b_sub.point(0);

                let a_sub_vecs = (0..i)
                    .map(|k| a_sub.point(k + 1) - a_sub_root)
                    .collect::<Vec<_>>();
                let b_sub_vecs = (0..j)
                    .map(|k| b_sub.point(k + 1) - b_sub_root)
                    .collect::<Vec<_>>();

                let normal_space = FinitelyFreeSubmoduleStructure::new(
                    FinitelyFreeModuleStructure::<FS, &'f FS>::new(field, space_dim),
                )
                .intersect(
                    MatrixStructure::new(space.field().clone()).col_kernel(Matrix::construct(
                        i + j,
                        space_dim,
                        |r, c| {
                            if r < i {
                                a_sub_vecs[r].coordinate(c).clone()
                            } else {
                                b_sub_vecs[r - i].coordinate(c).clone()
                            }
                        },
                    )),
                    linear_span_foo_submodule.clone(),
                );

                debug_assert!(normal_space.rank() >= 1);

                // If normal_space.rank() >= 1 then this is a degenerate face, there is no unique normal, and we need not check it
                if normal_space.rank() == 1 {
                    let normal =
                        space.vector(normal_space.basis().first().unwrap().iter().cloned());
                    for vec in &a_sub_vecs {
                        debug_assert!(field.is_zero(&normal.dot(vec)));
                    }
                    for vec in &b_sub_vecs {
                        debug_assert!(field.is_zero(&normal.dot(vec)));
                    }
                    Some(normal)
                } else {
                    None
                }
            })
    })) {
        if !normals.contains(&normal) {
            // Find the min and max projections of points of A
            let mut a_points = a.points().iter();
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

            // Find the min and max projections of points of B
            let mut b_points = b.points().iter();
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

            // Check how the projections of A and B overlap
            match (
                field.ring_cmp(&max_a_proj, &min_b_proj),
                field.ring_cmp(&max_b_proj, &min_a_proj),
            ) {
                (std::cmp::Ordering::Less, _) | (_, std::cmp::Ordering::Less) => {
                    return SimplexOverlapResult::Disjoint;
                }
                (std::cmp::Ordering::Greater, std::cmp::Ordering::Greater) => {}
                (std::cmp::Ordering::Equal, _) | (_, std::cmp::Ordering::Equal) => {
                    if RETURN_IF_TOUCHING {
                        return SimplexOverlapResult::Touching;
                    }
                    just_touching = true;
                }
            }

            normals.insert(normal);
        }
    }

    match just_touching {
        false => SimplexOverlapResult::Overlap,
        true => SimplexOverlapResult::Touching,
    }
}

pub fn simplex_interior_overlap<'f, FS: OrderedRingSignature + FieldSignature>(
    a: &Simplex<'f, FS>,
    b: &Simplex<'f, FS>,
) -> bool
where
    FS::Set: Hash,
{
    match simplex_overlap_impl::<FS, true>(a, b) {
        SimplexOverlapResult::Disjoint => false,
        SimplexOverlapResult::Touching => false,
        SimplexOverlapResult::Overlap => true,
    }
}

pub fn simplex_closure_overlap<'f, FS: OrderedRingSignature + FieldSignature>(
    a: &Simplex<'f, FS>,
    b: &Simplex<'f, FS>,
) -> bool
where
    FS::Set: Hash,
{
    match simplex_overlap_impl::<FS, false>(a, b) {
        SimplexOverlapResult::Disjoint => false,
        SimplexOverlapResult::Touching => true,
        SimplexOverlapResult::Overlap => true,
    }
}

pub fn simplex_overlap<'f, FS: OrderedRingSignature + FieldSignature>(
    a: &Simplex<'f, FS>,
    b: &Simplex<'f, FS>,
) -> SimplexOverlapResult
where
    FS::Set: Hash,
{
    simplex_overlap_impl::<FS, false>(a, b)
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
