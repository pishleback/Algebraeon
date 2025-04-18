use super::{finitely_free_modules::FinitelyFreeModuleStructure, matrix::Matrix};
use crate::{linear::matrix::*, structure::*};
use algebraeon_sets::structure::*;

#[derive(Debug, Clone)]
pub struct FinitelyFreeSubmodule<Ring: ReducedHermiteAlgorithmSignature> {
    ring: Ring,
    // a matrix in reduced hermite normal form with all non-zero rows whose rows form a basis for the submodule
    row_basis: Matrix<Ring::Set>,
    // the columns of the pivots of row_basis
    pivots: Vec<usize>,
}

impl<Ring: ReducedHermiteAlgorithmSignature> FinitelyFreeSubmodule<Ring> {
    pub fn ring(&self) -> &Ring {
        &self.ring
    }

    pub fn module_rank(&self) -> usize {
        self.row_basis.cols()
    }

    pub fn submodule_rank(&self) -> usize {
        self.row_basis.rows()
    }

    pub fn module(&self) -> FinitelyFreeModuleStructure<Ring> {
        FinitelyFreeModuleStructure::new(self.ring().clone(), self.module_rank())
    }

    pub fn submodule(&self) -> FinitelyFreeModuleStructure<Ring> {
        FinitelyFreeModuleStructure::new(self.ring().clone(), self.submodule_rank())
    }

    pub fn into_row_basis(self) -> Matrix<Ring::Set> {
        self.row_basis
    }

    pub fn into_col_basis(self) -> Matrix<Ring::Set> {
        self.row_basis.transpose()
    }

    pub fn basis(&self) -> Vec<Vec<Ring::Set>> {
        (0..self.row_basis.rows())
            .map(|r| {
                (0..self.module().rank())
                    .map(|c| self.row_basis.at(r, c).unwrap().clone())
                    .collect()
            })
            .collect()
    }

    pub fn matrix_row_span(ring: Ring, matrix: Matrix<Ring::Set>) -> Self {
        let cols = matrix.cols();
        let (h, _u, _det, pivots) =
            MatrixStructure::new(ring.clone()).row_reduced_hermite_algorithm(matrix);
        let row_basis = h.submatrix((0..pivots.len()).collect(), (0..cols).collect());
        FinitelyFreeSubmodule {
            ring,
            row_basis,
            pivots,
        }
    }

    pub fn matrix_col_span(ring: Ring, matrix: Matrix<Ring::Set>) -> Self {
        Self::matrix_row_span(ring, matrix.transpose())
    }

    pub fn matrix_row_kernel(ring: Ring, matrix: Matrix<Ring::Set>) -> Self {
        let rows = matrix.rows();
        let (_h, u, _u_det, pivs) =
            MatrixStructure::new(ring.clone()).row_hermite_algorithm(matrix);
        debug_assert_eq!(rows, u.rows());
        debug_assert_eq!(rows, u.cols());
        let ker = u.submatrix((pivs.len()..rows).collect(), (0..rows).collect());
        Self::matrix_row_span(ring, ker)
    }

    pub fn matrix_col_kernel(ring: Ring, matrix: Matrix<Ring::Set>) -> Self {
        Self::matrix_row_kernel(ring, matrix.transpose())
    }

    pub fn reduce_element(&self, element: &Vec<Ring::Set>) -> (Vec<Ring::Set>, Vec<Ring::Set>) {
        debug_assert!(self.module().is_element(&element));
        let mut reduced_element = element.clone();
        let mut offset = vec![];
        for (r, &c) in self.pivots.iter().enumerate() {
            let quo = self
                .ring()
                .quo(&reduced_element[c], self.row_basis.at(r, c).unwrap())
                .unwrap();
            for c2 in 0..self.module_rank() {
                reduced_element[c2] = self.ring().add(
                    &reduced_element[c2],
                    &self
                        .ring()
                        .neg(&self.ring().mul(&quo, &self.row_basis.at(r, c2).unwrap())),
                );
            }
            offset.push(quo);
        }
        (offset, reduced_element)
    }
}

/// The set of all submodules of some module represented by a basis in hermite normal form
/// UNIQUE: The reduced hermite normal form representing a submodule is unique
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FinitelyFreeSubmodulesStructure<Ring: RingSignature, const UNIQUE: bool> {
    module: FinitelyFreeModuleStructure<Ring>,
}

impl<Ring: ReducedHermiteAlgorithmSignature> FinitelyFreeSubmodulesStructure<Ring, false> {
    pub fn new_nonunique_reduction(module: FinitelyFreeModuleStructure<Ring>) -> Self {
        Self { module }
    }
}

impl<Ring: UniqueReducedHermiteAlgorithmSignature> FinitelyFreeSubmodulesStructure<Ring, true> {
    pub fn new(module: FinitelyFreeModuleStructure<Ring>) -> Self {
        Self { module }
    }
}

impl<Ring: ReducedHermiteAlgorithmSignature, const UNIQUE: bool> Signature
    for FinitelyFreeSubmodulesStructure<Ring, UNIQUE>
{
}

impl<Ring: ReducedHermiteAlgorithmSignature, const UNIQUE: bool> SetSignature
    for FinitelyFreeSubmodulesStructure<Ring, UNIQUE>
{
    type Set = FinitelyFreeSubmodule<Ring>;

    fn is_element(&self, sm: &Self::Set) -> bool {
        if sm.ring() != self.ring() {
            return false;
        }
        if sm.row_basis.cols() != self.module().rank() {
            return false;
        }
        // todo check sm.row_basis is in reduced hermite normal form with all non-zero rows
        true
    }
}

impl<Ring: ReducedHermiteAlgorithmSignature, const UNIQUE: bool> EqSignature
    for FinitelyFreeSubmodulesStructure<Ring, UNIQUE>
{
    fn equal(&self, x: &Self::Set, y: &Self::Set) -> bool {
        if UNIQUE {
            MatrixStructure::new(self.ring().clone()).equal(&x.row_basis, &y.row_basis)
        } else {
            self.contains(x, y) && self.contains(y, x)
        }
    }
}

impl<Ring: ReducedHermiteAlgorithmSignature, const UNIQUE: bool>
    SubModuleSignature<Ring, FinitelyFreeModuleStructure<Ring>>
    for FinitelyFreeSubmodulesStructure<Ring, UNIQUE>
{
    fn ring(&self) -> &Ring {
        self.module.ring()
    }

    fn module(&self) -> &FinitelyFreeModuleStructure<Ring> {
        &self.module
    }

    fn improper_submodule(&self) -> Self::Set {
        Self::Set::matrix_row_span(
            self.ring().clone(),
            MatrixStructure::new(self.ring().clone()).ident(self.module().rank()),
        )
    }

    fn add(&self, x: &Self::Set, y: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(&x));
        debug_assert!(self.is_element(&y));
        Self::Set::matrix_row_span(
            self.ring().clone(),
            Matrix::join_rows(
                self.module().rank(),
                vec![x.clone().into_row_basis(), y.clone().into_row_basis()],
            ),
        )
    }

    fn intersect(&self, x: &Self::Set, y: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(&x));
        debug_assert!(self.is_element(&y));
        let x_rows = x.clone().into_row_basis();
        let y_rows = y.clone().into_row_basis();
        let matrix = Matrix::join_rows(self.module().rank(), vec![&x_rows, &y_rows]);
        let matrix_ker = Self::Set::matrix_row_kernel(self.ring().clone(), matrix).into_row_basis();
        let matrix_ker_first_part = matrix_ker.submatrix(
            (0..matrix_ker.rows()).collect(),
            (0..x_rows.rows()).collect(),
        );
        Self::Set::matrix_row_span(
            self.ring().clone(),
            MatrixStructure::new(self.ring().clone())
                .mul(&matrix_ker_first_part, &x_rows)
                .unwrap(),
        )
    }

    fn generated(&self, generators: Vec<&Vec<Ring::Set>>) -> Self::Set {
        for generator in &generators {
            debug_assert!(self.module.is_element(generator));
        }
        let row_span = Matrix::construct(generators.len(), self.module().rank(), |r, c| {
            generators[r][c].clone()
        });
        Self::Set::matrix_row_span(self.ring().clone(), row_span)
    }

    fn contains_element(&self, submodule: &Self::Set, element: &Vec<Ring::Set>) -> bool {
        debug_assert!(self.is_element(&submodule));
        debug_assert!(self.module().is_element(&element));
        let (_offset, element_reduced) = submodule.reduce_element(element);
        element_reduced
            .iter()
            .all(|coeff| self.ring().is_zero(coeff))
    }

    fn contains(&self, x: &Self::Set, y: &Self::Set) -> bool {
        debug_assert!(self.is_element(&x));
        debug_assert!(self.is_element(&y));
        for b in y.basis() {
            if !self.contains_element(&x, &b) {
                return false;
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_nzq::Integer;

    #[test]
    fn test_finitely_free_submodule_reduction() {
        let a = Matrix::from_rows(vec![vec![1, 2, 4, 5], vec![1, 2, 3, 4]]);
        a.pprint();
        let a_reduced =
            FinitelyFreeSubmodule::matrix_row_span(Integer::structure(), a).into_row_basis();
        let a_expected = Matrix::from_rows(vec![vec![1, 2, 0, 1], vec![0, 0, 1, 1]]);
        a_reduced.pprint();
        a_expected.pprint();
        assert_eq!(a_reduced, a_expected);
    }

    #[test]
    fn test_finitely_free_submodule_unreduced_equal() {
        let module = FinitelyFreeModuleStructure::new(Integer::structure(), 4);
        let submodule = FinitelyFreeSubmodulesStructure::new(module.clone());

        assert!(submodule.equal(
            &FinitelyFreeSubmodule::matrix_row_span(
                Integer::structure(),
                Matrix::from_rows(vec![vec![1, 2, 3, 4], vec![0, 0, 1, 1]])
            ),
            &FinitelyFreeSubmodule::matrix_row_span(
                Integer::structure(),
                Matrix::from_rows(vec![vec![1, 2, 0, 1], vec![0, 0, 1, 1]])
            )
        ))
    }

    #[test]
    fn test_finitely_free_submodule_intersect() {
        let module = FinitelyFreeModuleStructure::new(Integer::structure(), 4);
        let submodule = FinitelyFreeSubmodulesStructure::new_nonunique_reduction(module.clone());

        let a = FinitelyFreeSubmodule::matrix_row_span(
            Integer::structure(),
            Matrix::from_rows(vec![
                vec![2, 0, 0, 0],
                vec![0, 2, 0, 0],
                vec![0, 0, 2, 0],
                vec![0, 0, 0, 2],
            ]),
        );

        let b = FinitelyFreeSubmodule::matrix_row_span(
            Integer::structure(),
            Matrix::from_rows(vec![vec![3, 3, 3, 3], vec![3, 3, -3, -3]]),
        );

        let c = FinitelyFreeSubmodule::matrix_row_span(
            Integer::structure(),
            Matrix::from_rows(vec![vec![6, 6, 0, 0], vec![0, 0, 6, 6]]),
        );

        let s = submodule.intersect(&a, &b);

        s.clone().into_row_basis().pprint();

        assert!(submodule.equal(&c, &s));
    }

    #[test]
    fn test_finitely_free_submodule_element_reduction() {
        let a = FinitelyFreeSubmodule::matrix_row_span(
            Integer::structure(),
            Matrix::from_rows(vec![vec![3, 2, 0, 3], vec![0, 14, 3, 1], vec![0, 0, 0, 10]]),
        );
        a.clone().into_row_basis().pprint();

        let element = vec![20, 20, 20, 20]
            .into_iter()
            .map(|x| Integer::from(x))
            .collect::<Vec<_>>();
        println!("element = {:?}", element);

        let (offset, reduced_element) = a.reduce_element(&element);

        println!("offset = {:?}", offset);
        println!("reduced_element = {:?}", reduced_element);
        debug_assert_eq!(
            reduced_element,
            vec![2, 8, 20, 2]
                .into_iter()
                .map(|x| Integer::from(x))
                .collect::<Vec<_>>()
        )
    }
}
