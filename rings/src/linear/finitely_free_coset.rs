use crate::structure::ModuleSignature;

use super::{
    finitely_free_affine::FinitelyFreeSubmoduleAffineSubset,
    finitely_free_modules::FinitelyFreeModuleStructure,
    finitely_free_submodule::FinitelyFreeSubmodule,
    matrix::{Matrix, ReducedHermiteAlgorithmSignature, UniqueReducedHermiteAlgorithmSignature},
};
use algebraeon_sets::structure::*;

#[derive(Debug, Clone)]
pub struct FinitelyFreeSubmoduleCoset<Ring: ReducedHermiteAlgorithmSignature> {
    // reduced coset offset
    offset: Vec<Ring::Set>,
    submodule: FinitelyFreeSubmodule<Ring>,
}

impl<Ring: ReducedHermiteAlgorithmSignature> FinitelyFreeSubmoduleCoset<Ring> {
    pub fn ring(&self) -> &Ring {
        self.submodule.ring()
    }

    pub fn module(&self) -> &FinitelyFreeModuleStructure<Ring> {
        self.submodule.module()
    }

    pub fn module_rank(&self) -> usize {
        self.submodule.module_rank()
    }

    pub fn coset_rank(&self) -> usize {
        self.submodule.submodule_rank()
    }

    pub fn offset(&self) -> &Vec<Ring::Set> {
        &self.offset
    }

    pub fn submodule(&self) -> &FinitelyFreeSubmodule<Ring> {
        &self.submodule
    }

    pub fn into_offset_and_row_basis_matrix(self) -> (Vec<Ring::Set>, Matrix<Ring::Set>) {
        (self.offset, self.submodule.into_row_basis_matrix())
    }

    pub fn into_offset_and_col_basis_matrix(self) -> (Vec<Ring::Set>, Matrix<Ring::Set>) {
        (self.offset, self.submodule.into_col_basis_matrix())
    }

    pub fn from_offset_and_module(
        offset: &Vec<Ring::Set>,
        submodule: FinitelyFreeSubmodule<Ring>,
    ) -> Self {
        let (_, offset) = submodule.reduce_element(offset);
        Self { offset, submodule }
    }

    pub fn equal_slow(x: &Self, y: &Self) -> bool {
        let module = common_structure::<FinitelyFreeModuleStructure<Ring>>(x.module(), y.module());
        if !FinitelyFreeSubmodule::equal_slow(&x.submodule, &y.submodule) {
            return false;
        }
        x.submodule
            .contains_element(&module.add(&x.offset, &module.neg(&y.offset)))
    }

    pub fn add(x: &Self, y: &Self) -> Self {
        let module = common_structure::<FinitelyFreeModuleStructure<Ring>>(x.module(), y.module());
        Self::from_offset_and_module(
            &module.add(&x.offset, &y.offset),
            FinitelyFreeSubmodule::add(&x.submodule, &y.submodule),
        )
    }

    pub fn intersect(x: &Self, y: &Self) -> FinitelyFreeSubmoduleAffineSubset<Ring> {
        let module = common_structure::<FinitelyFreeModuleStructure<Ring>>(x.module(), y.module());
        let cols = x.module_rank();
        debug_assert_eq!(cols, y.module_rank());
        let x_offset = x.offset();
        let y_offset = y.offset();
        let x_linear_row_basis = x.submodule().row_basis_matrix();
        let y_linear_row_basis = y.submodule().row_basis_matrix();

        /*
        The trick is to try and find one element in the intersection of x and y to take as the offset.
         - If no such element exists then the intersection is empty
         - If such an element does exist then we can take it as an offset for the intersection of the two submodules

        To test if x intersect y is non-empty we can express each as an affine span
        e.g. in the 2 dimensional case x and y can both be expressed as all sums of 3 points where the coefficients add to 1
        If we forget that all coefficients must add to 1 then this is a linear intersection
        Adding an extra first coordinate lets us track also the sums of the coefficients in the linear intersection
        Now we just need to check if the linear intersection contains a point with first coordinate equal to 1
        In row reduced form this is just checking the top left element of the row reduced matrix
         */

        let x_affine_row_basis = Matrix::construct(x_linear_row_basis.rows() + 1, cols, |r, c| {
            let v = if r == 0 {
                &module.ring().zero()
            } else {
                x_linear_row_basis.at(r - 1, c).unwrap()
            };
            module.ring().add(v, &x_offset[c])
        });

        let y_affine_row_basis = Matrix::construct(y_linear_row_basis.rows() + 1, cols, |r, c| {
            let v = if r == 0 {
                &module.ring().zero()
            } else {
                y_linear_row_basis.at(r - 1, c).unwrap()
            };
            module.ring().add(v, &y_offset[c])
        });

        let linearlized_intersection_row_basis = FinitelyFreeSubmodule::intersect(
            &FinitelyFreeSubmodule::matrix_row_span(
                module.ring().clone(),
                Matrix::join_cols(
                    x_affine_row_basis.rows(),
                    vec![
                        &Matrix::construct(x_affine_row_basis.rows(), 1, |_, _| {
                            module.ring().one()
                        }),
                        &x_affine_row_basis,
                    ],
                ),
            ),
            &FinitelyFreeSubmodule::matrix_row_span(
                module.ring().clone(),
                Matrix::join_cols(
                    y_affine_row_basis.rows(),
                    vec![
                        &Matrix::construct(y_affine_row_basis.rows(), 1, |_, _| {
                            module.ring().one()
                        }),
                        &y_affine_row_basis,
                    ],
                ),
            ),
        )
        .into_row_basis_matrix();
        // the coordinates to keep track of coefficient sum will intersect somewhere
        debug_assert!(linearlized_intersection_row_basis.rows() >= 1);
        debug_assert!(linearlized_intersection_row_basis.cols() >= 1);

        #[cfg(debug_assertions)]
        {
            if module
                .ring()
                .is_unit(linearlized_intersection_row_basis.at(0, 0).unwrap())
            {
                assert!(module.ring().equal(
                    linearlized_intersection_row_basis.at(0, 0).unwrap(),
                    &module.ring().one()
                ));
            }
        }

        if module.ring().equal(
            linearlized_intersection_row_basis.at(0, 0).unwrap(),
            &module.ring().one(),
        ) {
            let offset = (0..cols)
                .map(|c| {
                    linearlized_intersection_row_basis
                        .at(0, c + 1)
                        .unwrap()
                        .clone()
                })
                .collect();
            FinitelyFreeSubmoduleAffineSubset::from_coset(Self::from_offset_and_module(
                &offset,
                FinitelyFreeSubmodule::intersect(x.submodule(), y.submodule()),
            ))
        } else {
            FinitelyFreeSubmoduleAffineSubset::new_empty(module)
        }
    }

    pub fn into_affine_subset(self) -> FinitelyFreeSubmoduleAffineSubset<Ring> {
        FinitelyFreeSubmoduleAffineSubset::from_coset(self)
    }
}

impl<Ring: UniqueReducedHermiteAlgorithmSignature> FinitelyFreeSubmoduleCoset<Ring> {
    pub fn equal(x: &Self, y: &Self) -> bool {
        let module = common_structure::<FinitelyFreeModuleStructure<Ring>>(x.module(), y.module());
        module.equal(&x.offset, &y.offset)
            && FinitelyFreeSubmodule::equal(&x.submodule, &y.submodule)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_nzq::Integer;

    #[test]
    fn test_cosets() {
        let coset1 = Matrix::<Integer>::from_rows(vec![vec![15, 0], vec![0, 10]])
            .row_span()
            .coset(&vec![Integer::from(3), Integer::from(3)]);
        let coset2 = Matrix::<Integer>::from_rows(vec![vec![10, 0], vec![0, 15]])
            .row_span()
            .coset(&vec![Integer::from(4), Integer::from(4)]);

        println!("coset1");
        println!("{:?}", coset1.offset());
        coset1.submodule().clone().into_row_basis_matrix().pprint();

        println!();

        println!("coset2");
        println!("{:?}", coset2.offset());
        coset2.submodule().clone().into_row_basis_matrix().pprint();

        println!();

        let coset1_add_coset2 = FinitelyFreeSubmoduleCoset::add(&coset1, &coset2);
        println!("coset1 + coset2");
        println!("{:?}", coset1_add_coset2.offset());
        coset1_add_coset2
            .submodule()
            .clone()
            .into_row_basis_matrix()
            .pprint();
        assert_eq!(
            coset1_add_coset2.offset(),
            &vec![Integer::from(2), Integer::from(2)]
        );
        assert!(FinitelyFreeSubmodule::equal(
            coset1_add_coset2.submodule(),
            &Matrix::<Integer>::from_rows(vec![vec![5, 0], vec![0, 5]]).row_span()
        ));

        let coset1_intersect_coset2 = FinitelyFreeSubmoduleCoset::intersect(&coset1, &coset2);
        assert!(coset1_intersect_coset2.is_empty());
    }

    #[test]
    fn test_cosets_intersect() {
        let coset1 = Matrix::<Integer>::from_rows(vec![vec![6, 4], vec![0, 3]])
            .row_span()
            .coset(&vec![Integer::from(2), Integer::from(3)]);
        let coset2 = Matrix::<Integer>::from_rows(vec![vec![15, 6], vec![0, 7]])
            .row_span()
            .coset(&vec![Integer::from(5), Integer::from(4)]);

        println!("coset1");
        println!("{:?}", coset1.offset());
        coset1.submodule().clone().into_row_basis_matrix().pprint();

        println!();

        println!("coset2");
        println!("{:?}", coset2.offset());
        coset2.submodule().clone().into_row_basis_matrix().pprint();

        println!();

        let coset1_intersect_coset2 =
            FinitelyFreeSubmoduleCoset::intersect(&coset1, &coset2).unwrap_to_coset();
        println!("coset1 & coset2");
        println!("{:?}", coset1_intersect_coset2.offset());
        coset1_intersect_coset2
            .submodule()
            .clone()
            .into_row_basis_matrix()
            .pprint();

        assert_eq!(
            coset1_intersect_coset2.into_offset_and_row_basis_matrix(),
            (
                vec![Integer::from(20), Integer::from(3)],
                Matrix::from_rows(vec![vec![30, 5], vec![0, 21]])
            )
        );
    }
}
