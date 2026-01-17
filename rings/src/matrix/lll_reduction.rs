use std::str::FromStr;

use crate::{
    isolated_algebraic::polynomial::Interleave,
    linear::finitely_free_module::RingToFinitelyFreeModuleSignature,
    matrix::{Matrix, MatrixStructure, RealInnerProduct, SymmetricMatrix},
    structure::{FieldSignature, RealSubsetSignature},
};
use algebraeon_nzq::{Integer, Rational};
use algebraeon_sets::structure::{BorrowedStructure, MetaType, SetSignature};

// fn is_positive_definite(mat: &SymmetricMatrix<Integer>) -> bool {
//     let n = mat.n();
//     let mat = Matrix::construct(n, n, |r, c| mat.get(r, c).unwrap().clone());
//     for i in 1..n {
//         if mat
//             .submatrix((0..i).collect(), (0..i).collect())
//             .det()
//             .unwrap()
//             <= Integer::ZERO
//         {
//             return false;
//         }
//     }
//     true
// }

// fn integral_lll_reduction(
//     row_basis: Matrix<Integer>,
//     gram_matrix: &SymmetricMatrix<Integer>,
// ) -> Result<(), String> {
//     // source: Cohen H - a course in computational algebraic number theory
//     debug_assert!(is_positive_definite(gram_matrix));
//     let n = row_basis.rows(); // number of vectors in the basis
//     let m = row_basis.cols(); // dimension of the ambient space
//     assert_eq!(m, gram_matrix.n());
//     debug_assert_eq!(row_basis.rank(), n);

//     // keep track of the transformatinos performed on the lattice basis
//     // the columns express the linear combination of the input basis required to obtain the output basis
//     let mut h = Matrix::<Integer>::ident(n);
//     let mut k = 1;
//     let mut k_max = 0;
//     let mut d = vec![Integer::ZERO; n];
//     d[0] = todo!("row_basis.get_row(0) dot row_basis.get_row(0)");
//     let mut lambda = todo!();

//     loop {
//         if k > k_max {
//             k_max = k;
//             for j in 0..(k + 1) {
//                 let mut u: Integer = todo!("row_basis.get_row(k) dot row_basis.get_row(j)");
//                 for i in 0..j {
//                     u = (d[i] * u - lambda[k][i] * lambda[j][i])
//                         .try_divide(d[i - 1])
//                         .unwrap();
//                 }
//             }
//         }
//     }

//     todo!()
// }

impl<FS: RealSubsetSignature + FieldSignature, FSB: BorrowedStructure<FS>>
    MatrixStructure<FS, FSB>
{
    /// Take the rows of M = `mat` as the basis for a lattice.
    ///
    /// Perform an LLL reduction on the lattice, returning a matrix H such that the rows of H*M are an LLL reduced basis for the lattice.
    pub fn lll_row_reduction_algorithm(
        &self,
        mat: Matrix<FS::Set>,
        inner_product: &impl RealInnerProduct<FS>,
        delta: &Rational,
    ) -> Matrix<FS::Set> {
        assert!(mat.cols() >= mat.rows());
        debug_assert!(self.rank(mat.clone()) == mat.rows());
        // 1/4 < delta <= 1
        assert!(Rational::ONE < Rational::from(4) * delta && delta <= &Rational::ONE);
        let n = mat.rows(); // size of the basis
        let m = mat.cols(); // dimension of the ambient space

        println!("n={n} m={m}");

        todo!()
    }
}

impl<F: MetaType> Matrix<F>
where
    F::Signature: RealSubsetSignature + FieldSignature,
{
    pub fn lll_row_reduction_algorithm(
        self,
        inner_product: &impl RealInnerProduct<F::Signature>,
        delta: &Rational,
    ) -> Matrix<F> {
        Self::structure().lll_row_reduction_algorithm(self, inner_product, delta)
    }
}

#[cfg(test)]
mod tests {
    // use crate::matrix::{
    //     Matrix, SymmetricMatrix,
    //     lll_reduction::{integral_lll_reduction, is_positive_definite},
    // };
    // use algebraeon_nzq::Integer;

    // #[test]
    // fn test_is_positive_definite() {
    //     let mat = Matrix::from_rows(vec![
    //         vec![Integer::from(2), Integer::from(-1), Integer::from(0)],
    //         vec![Integer::from(-1), Integer::from(2), Integer::from(-1)],
    //         vec![Integer::from(0), Integer::from(-1), Integer::from(0)],
    //     ]);
    //     assert!(is_positive_definite(&mat.try_into().unwrap()));
    // }

    use crate::matrix::{Matrix, StandardInnerProduct};
    use algebraeon_nzq::Rational;
    use algebraeon_sets::structure::MetaType;
    use std::str::FromStr;

    #[test]
    fn test() {
        todo!();
    }
}
