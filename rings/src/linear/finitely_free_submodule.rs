use super::{finitely_free_modules::FinitelyFreeModuleStructure, matrix::Matrix};
use crate::{linear::matrix::*, structure::*};
use algebraeon_sets::structure::*;

#[derive(Debug, Clone)]
pub struct FinitelyFreeSubmodule<Ring: RingSignature> {
    // a matrix in hermite normal form with all non-zero rows whose rows form a basis for the submodule
    // may also be required to be reduced or unique depending on the structure used
    row_basis: Matrix<Ring::Set>,
    // the columns of the pivots of row_basis
    pivots: Vec<usize>,
}

impl<Ring: RingSignature> FinitelyFreeSubmodule<Ring> {
    pub fn into_row_basis(self) -> Matrix<Ring::Set> {
        self.row_basis
    }

    pub fn into_col_basis(self) -> Matrix<Ring::Set> {
        self.row_basis.transpose()
    }
}

/// The set of all submodules of some module represented by a basis in hermite normal form
/// REDUCED: Submodules are represented by matricies in _reduced_ hermite normal form
/// UNIQUE: The hermite normal form representing a submodule is _unique_
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FinitelyFreeSubmoduleStructure<
    Ring: RingSignature,
    const REDUCED: bool,
    const UNIQUE: bool,
> {
    module: FinitelyFreeModuleStructure<Ring>,
}

impl<Ring: RingSignature, const REDUCED: bool, const UNIQUE: bool>
    FinitelyFreeSubmoduleStructure<Ring, REDUCED, UNIQUE>
where
    Self: FinitelyFreeSubmoduleHermiteAlgorithm<Ring, UNIQUE>,
{
    pub fn basis(&self, sm: &FinitelyFreeSubmodule<Ring>) -> Vec<Vec<Ring::Set>> {
        debug_assert!(self.is_element(sm));
        (0..sm.row_basis.rows())
            .map(|r| {
                (0..self.module().rank())
                    .map(|c| sm.row_basis.at(r, c).unwrap().clone())
                    .collect()
            })
            .collect()
    }
}

impl<Ring: HermiteAlgorithmSignature> FinitelyFreeSubmoduleStructure<Ring, false, false> {
    pub fn new_unreduced(module: FinitelyFreeModuleStructure<Ring>) -> Self {
        Self { module }
    }
}

impl<Ring: ReducedHermiteAlgorithmSignature> FinitelyFreeSubmoduleStructure<Ring, true, false> {
    pub fn new_reduced(module: FinitelyFreeModuleStructure<Ring>) -> Self {
        Self { module }
    }
}

impl<Ring: UniqueReducedHermiteAlgorithmSignature>
    FinitelyFreeSubmoduleStructure<Ring, true, true>
{
    pub fn new_uniquely_reduced(module: FinitelyFreeModuleStructure<Ring>) -> Self {
        Self { module }
    }
}

mod hermite_algorithm_impl {
    use super::*;
    pub trait FinitelyFreeSubmoduleHermiteAlgorithm<Ring: RingSignature, const UNIQUE: bool>:
        SetSignature<Set = FinitelyFreeSubmodule<Ring>>
    {
        fn row_hermite_algorithm(
            &self,
            matrix: Matrix<Ring::Set>,
        ) -> (Matrix<Ring::Set>, Matrix<Ring::Set>, Ring::Set, Vec<usize>);
    }
}
use hermite_algorithm_impl::*;

impl<Ring: HermiteAlgorithmSignature, const UNIQUE: bool>
    FinitelyFreeSubmoduleHermiteAlgorithm<Ring, UNIQUE>
    for FinitelyFreeSubmoduleStructure<Ring, false, UNIQUE>
{
    fn row_hermite_algorithm(
        &self,
        matrix: Matrix<<Ring>::Set>,
    ) -> (
        Matrix<<Ring>::Set>,
        Matrix<<Ring>::Set>,
        <Ring>::Set,
        Vec<usize>,
    ) {
        MatrixStructure::new(self.ring().clone()).row_hermite_algorithm(matrix)
    }
}

impl<Ring: ReducedHermiteAlgorithmSignature, const UNIQUE: bool>
    FinitelyFreeSubmoduleHermiteAlgorithm<Ring, UNIQUE>
    for FinitelyFreeSubmoduleStructure<Ring, true, UNIQUE>
{
    fn row_hermite_algorithm(
        &self,
        matrix: Matrix<<Ring>::Set>,
    ) -> (
        Matrix<<Ring>::Set>,
        Matrix<<Ring>::Set>,
        <Ring>::Set,
        Vec<usize>,
    ) {
        MatrixStructure::new(self.ring().clone()).row_reduced_hermite_algorithm(matrix)
    }
}

impl<Ring: RingSignature, const REDUCED: bool, const UNIQUE: bool>
    FinitelyFreeSubmoduleStructure<Ring, REDUCED, UNIQUE>
where
    Self: FinitelyFreeSubmoduleHermiteAlgorithm<Ring, UNIQUE>,
{
    pub fn matrix_row_span(&self, matrix: Matrix<<Ring>::Set>) -> FinitelyFreeSubmodule<Ring> {
        debug_assert_eq!(matrix.cols(), self.module().rank());
        let (h, _u, _det, pivots) = self.row_hermite_algorithm(matrix);
        let row_basis = h.submatrix(
            (0..pivots.len()).collect(),
            (0..self.module().rank()).collect(),
        );
        FinitelyFreeSubmodule { row_basis, pivots }
    }

    pub fn matrix_col_span(&self, matrix: Matrix<Ring::Set>) -> FinitelyFreeSubmodule<Ring> {
        self.matrix_row_span(matrix.transpose())
    }

    pub fn matrix_row_kernel(&self, matrix: Matrix<Ring::Set>) -> FinitelyFreeSubmodule<Ring> {
        debug_assert_eq!(matrix.rows(), self.module().rank());
        let rows = matrix.rows();
        let (_h, u, _u_det, pivs) = self.row_hermite_algorithm(matrix);
        debug_assert_eq!(rows, u.rows());
        debug_assert_eq!(rows, u.cols());
        let ker = u.submatrix((pivs.len()..rows).collect(), (0..rows).collect());
        self.matrix_row_span(ker)
    }

    pub fn matrix_col_kernel(&self, matrix: Matrix<Ring::Set>) -> FinitelyFreeSubmodule<Ring> {
        self.matrix_row_kernel(matrix.transpose())
    }
}

impl<Ring: RingSignature, const REDUCED: bool, const UNIQUE: bool> Signature
    for FinitelyFreeSubmoduleStructure<Ring, REDUCED, UNIQUE>
where
    Self: FinitelyFreeSubmoduleHermiteAlgorithm<Ring, UNIQUE>,
{
}

impl<Ring: RingSignature, const REDUCED: bool, const UNIQUE: bool> SetSignature
    for FinitelyFreeSubmoduleStructure<Ring, REDUCED, UNIQUE>
where
    Self: FinitelyFreeSubmoduleHermiteAlgorithm<Ring, UNIQUE>,
{
    type Set = FinitelyFreeSubmodule<Ring>;

    fn is_element(&self, sm: &Self::Set) -> bool {
        if sm.row_basis.cols() != self.module().rank() {
            return false;
        }
        // todo check sm.row_basis is in hermite normal form with all non-zero rows
        if REDUCED {
            // todo check sm.row_basis is in _reduced_ hermite normal form
        }
        true
    }
}

impl<Ring: RingSignature, const REDUCED: bool, const UNIQUE: bool> EqSignature
    for FinitelyFreeSubmoduleStructure<Ring, REDUCED, UNIQUE>
where
    Self: FinitelyFreeSubmoduleHermiteAlgorithm<Ring, UNIQUE>,
{
    fn equal(&self, x: &Self::Set, y: &Self::Set) -> bool {
        if UNIQUE {
            MatrixStructure::new(self.ring().clone()).equal(&x.row_basis, &y.row_basis)
        } else {
            self.contains(x, y) && self.contains(y, x)
        }
    }
}

impl<Ring: RingSignature, const REDUCED: bool, const UNIQUE: bool>
    SubModuleSignature<Ring, FinitelyFreeModuleStructure<Ring>>
    for FinitelyFreeSubmoduleStructure<Ring, REDUCED, UNIQUE>
where
    Self: FinitelyFreeSubmoduleHermiteAlgorithm<Ring, UNIQUE>,
{
    fn ring(&self) -> &Ring {
        self.module.ring()
    }

    fn module(&self) -> &FinitelyFreeModuleStructure<Ring> {
        &self.module
    }

    fn improper_submodule(&self) -> Self::Set {
        self.matrix_row_span(MatrixStructure::new(self.ring().clone()).ident(self.module().rank()))
    }

    fn add(&self, x: &Self::Set, y: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(&x));
        debug_assert!(self.is_element(&y));
        self.matrix_row_span(Matrix::join_rows(
            self.module().rank(),
            vec![x.clone().into_row_basis(), y.clone().into_row_basis()],
        ))
    }

    fn intersect(&self, x: &Self::Set, y: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(&x));
        debug_assert!(self.is_element(&y));
        let x_rows = x.clone().into_row_basis();
        let y_rows = y.clone().into_row_basis();
        let matrix = Matrix::join_rows(self.module().rank(), vec![&x_rows, &y_rows]);
        let matrix_ker = Self {
            module: FinitelyFreeModuleStructure::new(self.ring().clone(), matrix.rows()),
        }
        .matrix_row_kernel(matrix)
        .into_row_basis();
        let matrix_ker_first_part = matrix_ker.submatrix(
            (0..matrix_ker.rows()).collect(),
            (0..x_rows.rows()).collect(),
        );
        self.matrix_row_span(
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
        self.matrix_row_span(row_span)
    }

    fn contains_element(&self, x: &Self::Set, p: &Vec<Ring::Set>) -> bool {
        debug_assert!(self.is_element(x));
        debug_assert!(self.module().is_element(&p));
        todo!()
    }

    fn contains(&self, x: &Self::Set, y: &Self::Set) -> bool {
        debug_assert!(self.is_element(&x));
        debug_assert!(self.is_element(&y));
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_nzq::Integer;

    #[test]
    fn test_finitely_free_submodule_reduction() {
        let module = FinitelyFreeModuleStructure::new(Integer::structure(), 4);
        {
            let a = Matrix::from_rows(vec![vec![1, 2, 4, 5], vec![1, 2, 3, 4]]);
            a.pprint();
            let a_reduced = FinitelyFreeSubmoduleStructure::new_uniquely_reduced(module.clone())
                .matrix_row_span(a)
                .into_row_basis();
            let a_expected = Matrix::from_rows(vec![vec![1, 2, 0, 1], vec![0, 0, 1, 1]]);
            a_reduced.pprint();
            a_expected.pprint();
            assert_eq!(a_reduced, a_expected);
        }
        println!();
        {
            let a = Matrix::from_rows(vec![vec![1, 2, 4, 5], vec![1, 2, 3, 4]]);
            a.pprint();
            let a_reduced = FinitelyFreeSubmoduleStructure::new_unreduced(module.clone())
                .matrix_row_span(a)
                .into_row_basis();
            let a_expected = Matrix::from_rows(vec![vec![1, 2, 3, 4], vec![0, 0, 1, 1]]);
            a_reduced.pprint();
            a_expected.pprint();
            assert_eq!(a_reduced, a_expected);
        }
    }

    #[test]
    fn test_finitely_free_submodule_unreduced_equal() {
        let module = FinitelyFreeModuleStructure::new(Integer::structure(), 4);
        let submodule = FinitelyFreeSubmoduleStructure::new_unreduced(module.clone());

        assert!(submodule.equal(
            &submodule.matrix_row_span(Matrix::from_rows(vec![vec![1, 2, 3, 4], vec![0, 0, 1, 1]])),
            &submodule.matrix_row_span(Matrix::from_rows(vec![vec![1, 2, 0, 1], vec![0, 0, 1, 1]]))
        ))
    }

    #[test]
    fn test_finitely_free_submodule_intersect() {
        let module = FinitelyFreeModuleStructure::new(Integer::structure(), 4);
        let submodule = FinitelyFreeSubmoduleStructure::new_uniquely_reduced(module.clone());

        let a = submodule.matrix_row_span(Matrix::from_rows(vec![
            vec![2, 0, 0, 0],
            vec![0, 2, 0, 0],
            vec![0, 0, 2, 0],
            vec![0, 0, 0, 2],
        ]));

        let b = submodule.matrix_row_span(Matrix::from_rows(vec![
            vec![3, 3, 3, 3],
            vec![3, 3, -3, -3],
        ]));

        let c =
            submodule.matrix_row_span(Matrix::from_rows(vec![vec![6, 6, 0, 0], vec![0, 0, 6, 6]]));

        let s = submodule.intersect(&a, &b);

        s.clone().into_row_basis().pprint();

        assert!(submodule.equal(&c, &s));
    }
}
