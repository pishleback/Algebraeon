use super::{finitely_free_affine::*, finitely_free_module::*, finitely_free_submodule::*};
use crate::{
    matrix::{Matrix, ReducedHermiteAlgorithmSignature, UniqueReducedHermiteAlgorithmSignature},
    structure::*,
};
use algebraeon_sets::structure::*;
use std::fmt::Debug;

#[derive(Debug, Clone)]
pub struct FinitelyFreeSubmoduleCoset<Set: Clone + Debug> {
    // reduced coset offset
    offset: Vec<Set>,
    submodule: FinitelyFreeSubmodule<Set>,
}

impl<Set: Clone + Debug> FinitelyFreeSubmoduleCoset<Set> {
    pub fn offset(&self) -> &Vec<Set> {
        &self.offset
    }

    pub fn submodule(&self) -> &FinitelyFreeSubmodule<Set> {
        &self.submodule
    }

    pub fn rank(&self) -> usize {
        self.submodule.rank()
    }

    pub fn into_offset_and_row_basis_matrix(self) -> (Vec<Set>, Matrix<Set>) {
        (self.offset, self.submodule.into_row_basis_matrix())
    }

    pub fn into_offset_and_col_basis_matrix(self) -> (Vec<Set>, Matrix<Set>) {
        (self.offset, self.submodule.into_col_basis_matrix())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FinitelyFreeSubmoduleCosetStructure<
    Ring: ReducedHermiteAlgorithmSignature,
    RingB: BorrowedStructure<Ring>,
> {
    module: FinitelyFreeModuleStructure<Ring, RingB>,
}

impl<Ring: ReducedHermiteAlgorithmSignature, RingB: BorrowedStructure<Ring>>
    FinitelyFreeSubmoduleCosetStructure<Ring, RingB>
{
    pub fn new(module: FinitelyFreeModuleStructure<Ring, RingB>) -> Self {
        Self { module }
    }
}

impl<Ring: ReducedHermiteAlgorithmSignature, RingB: BorrowedStructure<Ring>> Signature
    for FinitelyFreeSubmoduleCosetStructure<Ring, RingB>
{
}

impl<Ring: ReducedHermiteAlgorithmSignature, RingB: BorrowedStructure<Ring>> SetSignature
    for FinitelyFreeSubmoduleCosetStructure<Ring, RingB>
{
    type Set = FinitelyFreeSubmoduleCoset<Ring::Set>;

    fn is_element(&self, _x: &Self::Set) -> Result<(), String> {
        //TODO: better checks here
        Ok(())
    }
}

impl<Ring: ReducedHermiteAlgorithmSignature, RingB: BorrowedStructure<Ring>>
    FinitelyFreeSubmoduleCosetStructure<Ring, RingB>
{
    pub fn ring(&self) -> &Ring {
        self.module.ring()
    }

    pub fn module(&self) -> &FinitelyFreeModuleStructure<Ring, RingB> {
        &self.module
    }

    pub fn from_offset_and_submodule(
        &self,
        offset: &Vec<Ring::Set>,
        submodule: FinitelyFreeSubmodule<Ring::Set>,
    ) -> FinitelyFreeSubmoduleCoset<Ring::Set> {
        debug_assert!(self.module().is_element(offset).is_ok());
        debug_assert!(self.module().submodules().is_element(&submodule).is_ok());
        let (_, offset) = self
            .module()
            .submodules()
            .reduce_element(&submodule, offset);
        FinitelyFreeSubmoduleCoset { offset, submodule }
    }

    pub fn from_submodule(
        &self,
        submodule: FinitelyFreeSubmodule<Ring::Set>,
    ) -> FinitelyFreeSubmoduleCoset<Ring::Set> {
        FinitelyFreeSubmoduleCoset {
            offset: self.module().zero(),
            submodule,
        }
    }

    pub fn equal_slow(
        &self,
        x: &FinitelyFreeSubmoduleCoset<Ring::Set>,
        y: &FinitelyFreeSubmoduleCoset<Ring::Set>,
    ) -> bool {
        debug_assert!(self.is_element(x).is_ok());
        debug_assert!(self.is_element(y).is_ok());
        if !self
            .module()
            .submodules()
            .equal_slow(&x.submodule, &y.submodule)
        {
            return false;
        }
        self.module()
            .submodules()
            .contains_element(&x.submodule, &self.module().sub(&x.offset, &y.offset))
    }

    pub fn sum(
        &self,
        x: &FinitelyFreeSubmoduleCoset<Ring::Set>,
        y: &FinitelyFreeSubmoduleCoset<Ring::Set>,
    ) -> FinitelyFreeSubmoduleCoset<Ring::Set> {
        debug_assert!(self.is_element(x).is_ok());
        debug_assert!(self.is_element(y).is_ok());
        self.from_offset_and_submodule(
            &self.module().add(&x.offset, &y.offset),
            self.module().submodules().sum(&x.submodule, &y.submodule),
        )
    }

    pub fn intersect(
        &self,
        x: &FinitelyFreeSubmoduleCoset<Ring::Set>,
        y: &FinitelyFreeSubmoduleCoset<Ring::Set>,
    ) -> FinitelyFreeSubmoduleAffineSubset<Ring::Set> {
        debug_assert!(self.is_element(x).is_ok());
        debug_assert!(self.is_element(y).is_ok());
        let cols = self.module().rank();
        let x_offset = x.offset();
        let y_offset = y.offset();
        let x_linear_row_basis = x.submodule().row_basis_matrix();
        let y_linear_row_basis = y.submodule().row_basis_matrix();

        /*
        The trick is to find one element in the intersection of x and y to take as the offset.
         - If no such element exists then the intersection is empty
         - If such an element exists then we can take it as an offset for the intersection of the two submodules

        To test if x intersect y is non-empty we can express each as an affine span
        e.g. in the 2 dimensional case x and y can both be expressed as all sums of 3 points where the coefficients add to 1
        If we forget that all coefficients must add to 1 then this is a linear intersection
        Adding an extra first coordinate lets us track also the sums of the coefficients in the linear intersection
        Now we just need to check if the linear intersection contains a point with first coordinate equal to 1
        In row reduced form this is just checking the top left element of the row reduced matrix
         */

        let x_affine_row_basis = Matrix::construct(x_linear_row_basis.rows() + 1, cols, |r, c| {
            let v = if r == 0 {
                &self.ring().zero()
            } else {
                x_linear_row_basis.at(r - 1, c).unwrap()
            };
            self.ring().add(v, &x_offset[c])
        });

        let y_affine_row_basis = Matrix::construct(y_linear_row_basis.rows() + 1, cols, |r, c| {
            let v = if r == 0 {
                &self.ring().zero()
            } else {
                y_linear_row_basis.at(r - 1, c).unwrap()
            };
            self.ring().add(v, &y_offset[c])
        });

        let larger_module = self.ring().free_module_structure(self.module().rank() + 1);

        let linearlized_intersection_row_basis = larger_module
            .submodules()
            .intersect(
                &larger_module
                    .submodules()
                    .matrix_row_span(Matrix::join_cols(
                        x_affine_row_basis.rows(),
                        vec![
                            &Matrix::construct(x_affine_row_basis.rows(), 1, |_, _| {
                                self.ring().one()
                            }),
                            &x_affine_row_basis,
                        ],
                    )),
                &larger_module
                    .submodules()
                    .matrix_row_span(Matrix::join_cols(
                        y_affine_row_basis.rows(),
                        vec![
                            &Matrix::construct(y_affine_row_basis.rows(), 1, |_, _| {
                                self.ring().one()
                            }),
                            &y_affine_row_basis,
                        ],
                    )),
            )
            .into_row_basis_matrix();
        // the coordinates to keep track of coefficient sum will intersect somewhere
        debug_assert!(linearlized_intersection_row_basis.rows() >= 1);
        debug_assert!(linearlized_intersection_row_basis.cols() >= 1);

        #[cfg(debug_assertions)]
        {
            if self
                .ring()
                .is_unit(linearlized_intersection_row_basis.at(0, 0).unwrap())
            {
                assert!(self.ring().equal(
                    linearlized_intersection_row_basis.at(0, 0).unwrap(),
                    &self.ring().one()
                ));
            }
        }

        if self.ring().equal(
            linearlized_intersection_row_basis.at(0, 0).unwrap(),
            &self.ring().one(),
        ) {
            let offset = (0..cols)
                .map(|c| {
                    linearlized_intersection_row_basis
                        .at(0, c + 1)
                        .unwrap()
                        .clone()
                })
                .collect();
            FinitelyFreeSubmoduleAffineSubset::NonEmpty(
                self.from_offset_and_submodule(
                    &offset,
                    self.module()
                        .submodules()
                        .intersect(x.submodule(), y.submodule()),
                ),
            )
        } else {
            FinitelyFreeSubmoduleAffineSubset::Empty
        }
    }

    pub fn contains_element(
        &self,
        x: &FinitelyFreeSubmoduleCoset<Ring::Set>,
        p: &Vec<Ring::Set>,
    ) -> bool {
        debug_assert_eq!(p.len(), self.module().rank());
        self.module().submodules().contains_element(
            x.submodule(),
            &self.module().add(p, &self.module().neg(x.offset())),
        )
    }
}

impl<Ring: UniqueReducedHermiteAlgorithmSignature, RingB: BorrowedStructure<Ring>> EqSignature
    for FinitelyFreeSubmoduleCosetStructure<Ring, RingB>
{
    fn equal(&self, x: &Self::Set, y: &Self::Set) -> bool {
        self.module().equal(&x.offset, &y.offset)
            && self.module().submodules().equal(&x.submodule, &y.submodule)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_nzq::Integer;

    #[test]
    fn test_cosets() {
        let module = Integer::structure().into_free_module_structure(2);

        let coset1 = module.submodules().coset(
            &Matrix::<Integer>::from_rows(vec![vec![15, 0], vec![0, 10]]).row_span(),
            &vec![Integer::from(3), Integer::from(3)],
        );
        let coset2 = module.submodules().coset(
            &Matrix::<Integer>::from_rows(vec![vec![10, 0], vec![0, 15]]).row_span(),
            &vec![Integer::from(4), Integer::from(4)],
        );

        println!("coset1");
        println!("{:?}", coset1.offset());
        coset1.submodule().clone().into_row_basis_matrix().pprint();

        println!();

        println!("coset2");
        println!("{:?}", coset2.offset());
        coset2.submodule().clone().into_row_basis_matrix().pprint();

        println!();

        let coset1_add_coset2 = module.cosets().sum(&coset1, &coset2);
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
        assert!(module.submodules().equal(
            coset1_add_coset2.submodule(),
            &Matrix::<Integer>::from_rows(vec![vec![5, 0], vec![0, 5]]).row_span()
        ));

        let coset1_intersect_coset2 = module.cosets().intersect(&coset1, &coset2);
        assert!(coset1_intersect_coset2.is_empty());
    }

    #[test]
    fn test_cosets_intersect() {
        let module = Integer::structure().into_free_module_structure(2);

        let coset1 = module.submodules().coset(
            &Matrix::<Integer>::from_rows(vec![vec![6, 4], vec![0, 3]]).row_span(),
            &vec![Integer::from(2), Integer::from(3)],
        );
        let coset2 = module.submodules().coset(
            &Matrix::<Integer>::from_rows(vec![vec![15, 6], vec![0, 7]]).row_span(),
            &vec![Integer::from(5), Integer::from(4)],
        );

        println!("coset1");
        println!("{:?}", coset1.offset());
        coset1.submodule().clone().into_row_basis_matrix().pprint();

        println!();

        println!("coset2");
        println!("{:?}", coset2.offset());
        coset2.submodule().clone().into_row_basis_matrix().pprint();

        println!();

        let coset1_intersect_coset2 = module
            .cosets()
            .intersect(&coset1, &coset2)
            .unwrap_to_coset();
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
