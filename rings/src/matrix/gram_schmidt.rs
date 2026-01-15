use super::*;

// #[derive(Debug, Clone)]
// struct FiniteDimensionalInnerProductSpaceStructure<Ring: ComplexSubsetStructure> {
//     ring: Ring,
//     rows: usize,
//     cols: usize,
//     // Let n = self.rows * self.cols.
//     // Then self.values is an n by n matrix containing the inner products
//     // of the elementary basis vectors when expressed as metamatricies.
//     values: Matrix<Ring::Set>,
// }

// impl<Ring: ComplexSubsetStructure> PartialEq for FiniteDimensionalInnerProductSpaceStructure<Ring> {
//     fn eq(&self, other: &Self) -> bool {
//         self.ring == other.ring
//             && self.rows == other.rows
//             && self.cols == other.cols
//             && MatrixStructure::new(self.ring.clone()).equal(&self.values, &other.values)
//     }
// }

// impl<Ring: ComplexSubsetStructure> Eq for FiniteDimensionalInnerProductSpaceStructure<Ring> {}

impl<FS: ComplexConjugateSignature, FSB: BorrowedStructure<FS>> MatrixStructure<FS, FSB> {
    pub fn conjugate(&self, mat: &Matrix<FS::Set>) -> Matrix<FS::Set> {
        mat.apply_map(|x| self.ring().conjugate(x))
    }

    pub fn conjugate_transpose(&self, mat: &Matrix<FS::Set>) -> Matrix<FS::Set> {
        self.conjugate(mat).transpose()
    }
}

impl<FS: ComplexConjugateSignature + RingSignature, FSB: BorrowedStructure<FS>>
    MatrixStructure<FS, FSB>
{
    // dot product of a and conj(b)
    pub fn inner_product(&self, a: &Matrix<FS::Set>, b: &Matrix<FS::Set>) -> FS::Set {
        self.dot(a, &self.conjugate(b))
    }
}

impl<FS: ComplexConjugateSignature + FieldSignature, FSB: BorrowedStructure<FS>>
    MatrixStructure<FS, FSB>
{
    //return mat=LQ where L is lower triangular and Q is row-orthogonal (not orthonormal)
    pub fn gram_schmidt_row_orthogonalization_algorithm(
        &self,
        mut mat: Matrix<FS::Set>,
    ) -> (Matrix<FS::Set>, Matrix<FS::Set>) {
        #[cfg(debug_assertions)]
        let original_mat = mat.clone();

        let mut lt = self.ident(mat.rows());
        for i in 0..mat.rows() {
            for j in 0..i {
                //col_i = col_i - (projection of col_j onto col_i)
                let lambda = self
                    .ring()
                    .try_divide(
                        &self.inner_product(&mat.get_row_submatrix(i), &mat.get_row_submatrix(j)),
                        &self.inner_product(&mat.get_row_submatrix(j), &mat.get_row_submatrix(j)),
                    )
                    .unwrap();
                //col_i -= lambda col_j
                let row_opp = ElementaryOpp::new_row_opp(
                    self.ring().clone(),
                    ElementaryOppType::AddRowMul {
                        i,
                        j,
                        x: self.ring().neg(&lambda),
                    },
                );
                row_opp.apply(&mut lt);
                row_opp.apply(&mut mat);
            }
        }

        for i in 0..mat.rows() {
            for j in (i + 1)..mat.rows() {
                debug_assert!(self.ring().is_zero(
                    &self.inner_product(&mat.get_row_submatrix(i), &mat.get_row_submatrix(j)),
                ));
            }
        }

        #[cfg(debug_assertions)]
        assert!(self.equal(&mat, &self.mul(&lt, &original_mat).unwrap()));

        (lt, mat)
    }

    //return mat=QR where Q is col-orthogonal (not orthonormal) and R is upper triangular and
    pub fn gram_schmidt_col_orthogonalization_algorithm(
        &self,
        mat: Matrix<FS::Set>,
    ) -> (Matrix<FS::Set>, Matrix<FS::Set>) {
        let (l, q) = self.gram_schmidt_row_orthogonalization_algorithm(mat.transpose());
        (q.transpose(), l.transpose())
    }

    pub fn gram_schmidt_row_orthogonalization(&self, mat: Matrix<FS::Set>) -> Matrix<FS::Set> {
        self.gram_schmidt_row_orthogonalization_algorithm(mat).1
    }

    pub fn gram_schmidt_col_orthogonalization(&self, mat: Matrix<FS::Set>) -> Matrix<FS::Set> {
        self.gram_schmidt_col_orthogonalization_algorithm(mat).0
    }
}

impl<
    FS: ComplexConjugateSignature + PositiveRealNthRootSignature + FieldSignature,
    FSB: BorrowedStructure<FS>,
> MatrixStructure<FS, FSB>
{
    //return L*mat=Q where L is lower triangular and Q is orthonormal
    pub fn lq_decomposition_algorithm(
        &self,
        mat: Matrix<FS::Set>,
    ) -> (Matrix<FS::Set>, Matrix<FS::Set>) {
        let (mut lt, mut mat) = self.gram_schmidt_row_orthogonalization_algorithm(mat);

        for r in 0..mat.rows() {
            let row = mat.get_row_submatrix(r);
            let lensq = self.inner_product(&row, &row);
            let row_opp = ElementaryOpp::new_row_opp(
                self.ring().clone(),
                ElementaryOppType::UnitMul {
                    row: r,
                    unit: self
                        .ring()
                        .try_reciprocal(&self.ring().nth_root(&lensq, 2).unwrap())
                        .unwrap(),
                },
            );
            row_opp.apply(&mut lt);
            row_opp.apply(&mut mat);
        }

        debug_assert!(
            self.equal(
                &self.ident(mat.rows()),
                &self
                    .mul(&mat, &self.conjugate(&mat.transpose_ref()))
                    .unwrap()
            )
        );

        (lt, mat)
    }

    //return mat*R=Q where Q is col-orthogonal (not orthonormal) and R is upper triangular
    pub fn qr_decomposition_algorithm(
        &self,
        mat: Matrix<FS::Set>,
    ) -> (Matrix<FS::Set>, Matrix<FS::Set>) {
        let (l, q) = self.lq_decomposition_algorithm(mat.transpose());
        (q.transpose(), l.transpose())
    }

    pub fn gram_schmidt_row_orthonormalization(&self, mat: Matrix<FS::Set>) -> Matrix<FS::Set> {
        self.lq_decomposition_algorithm(mat).1
    }

    pub fn gram_schmidt_col_orthonormalization(&self, mat: Matrix<FS::Set>) -> Matrix<FS::Set> {
        self.qr_decomposition_algorithm(mat).0
    }
}

impl<F: MetaType> Matrix<F>
where
    F::Signature: ComplexConjugateSignature + FieldSignature,
{
    pub fn gram_schmidt_row_orthogonalization_algorithm(self) -> (Matrix<F>, Matrix<F>) {
        Self::structure().gram_schmidt_row_orthogonalization_algorithm(self)
    }

    pub fn gram_schmidt_col_orthogonalization_algorithm(self) -> (Matrix<F>, Matrix<F>) {
        Self::structure().gram_schmidt_col_orthogonalization_algorithm(self)
    }

    pub fn gram_schmidt_row_orthogonalization(self) -> Matrix<F> {
        Self::structure().gram_schmidt_row_orthogonalization(self)
    }

    pub fn gram_schmidt_col_orthogonalization(self) -> Matrix<F> {
        Self::structure().gram_schmidt_col_orthogonalization(self)
    }
}

impl<F: MetaType> Matrix<F>
where
    F::Signature: ComplexConjugateSignature + PositiveRealNthRootSignature + FieldSignature,
{
    pub fn lq_decomposition_algorithm(self) -> (Matrix<F>, Matrix<F>) {
        Self::structure().lq_decomposition_algorithm(self)
    }

    pub fn qr_decomposition_algorithm(self) -> (Matrix<F>, Matrix<F>) {
        Self::structure().qr_decomposition_algorithm(self)
    }

    pub fn gram_schmidt_row_orthonormalization(self) -> Matrix<F> {
        Self::structure().gram_schmidt_row_orthonormalization(self)
    }

    pub fn gram_schmidt_col_orthonormalization(self) -> Matrix<F> {
        Self::structure().gram_schmidt_col_orthonormalization(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::isolated_algebraic::{ComplexAlgebraic, RealAlgebraic};
    use std::str::FromStr;

    #[test]
    fn rational_gram_schmidt() {
        let mat = Matrix::<Rational>::from_rows(vec![
            vec![Rational::from(1), Rational::from(-1), Rational::from(3)],
            vec![Rational::from(1), Rational::from(0), Rational::from(5)],
            vec![Rational::from(1), Rational::from(2), Rational::from(6)],
        ]);

        let mat_expected_gs = Matrix::from_rows(vec![
            vec![
                Rational::from_str("1").unwrap(),
                Rational::from_str("-4/3").unwrap(),
                Rational::from_str("-3/7").unwrap(),
            ],
            vec![
                Rational::from_str("1").unwrap(),
                Rational::from_str("-1/3").unwrap(),
                Rational::from_str("9/14").unwrap(),
            ],
            vec![
                Rational::from_str("1").unwrap(),
                Rational::from_str("5/3").unwrap(),
                Rational::from_str("-3/14").unwrap(),
            ],
        ]);

        mat.pprint();
        mat_expected_gs.pprint();

        let (mat_actual_gs, mat_actual_ut) =
            mat.clone().gram_schmidt_col_orthogonalization_algorithm();

        mat_actual_gs.pprint();
        mat_actual_ut.pprint();
        Matrix::mul(&mat, &mat_actual_ut).unwrap().pprint();

        assert_eq!(mat.gram_schmidt_col_orthogonalization(), mat_expected_gs);
    }

    #[allow(clippy::erasing_op)]
    #[test]
    fn complex_gram_schmidt() {
        let i = &ComplexAlgebraic::i().into_ergonomic();

        let mat = Matrix::<ComplexAlgebraic>::from_rows(vec![
            vec![(1 + 0 * i).into_verbose(), (1 * i).into_verbose()],
            vec![(1 + 0 * i).into_verbose(), (1 + 0 * i).into_verbose()],
        ]);
        mat.pprint();
        mat.gram_schmidt_col_orthogonalization().pprint();

        let mat = Matrix::<ComplexAlgebraic>::from_rows(vec![
            vec![
                (-2 + 2 * i).into_verbose(),
                (7 + 3 * i).into_verbose(),
                (7 + 3 * i).into_verbose(),
            ],
            vec![
                (3 + 3 * i).into_verbose(),
                (-2 + 4 * i).into_verbose(),
                (6 + 2 * i).into_verbose(),
            ],
            vec![
                (2 + 2 * i).into_verbose(),
                (8 + 0 * i).into_verbose(),
                (-9 + 1 * i).into_verbose(),
            ],
        ]);
        mat.pprint();
        mat.clone().gram_schmidt_col_orthogonalization().pprint();
    }

    #[test]
    fn complex_gram_schmidt_normalized() {
        let one = &RealAlgebraic::one().into_ergonomic();

        let mat = Matrix::<RealAlgebraic>::from_rows(vec![
            vec![(1 * one).into_verbose(), (1 * one).into_verbose()],
            vec![(1 * one).into_verbose(), (2 * one).into_verbose()],
        ]);
        mat.pprint();
        mat.gram_schmidt_col_orthonormalization().pprint();

        let i = &ComplexAlgebraic::i().into_ergonomic();
        let mat = Matrix::<ComplexAlgebraic>::from_rows(vec![
            vec![(-2 + 2 * i).into_verbose(), (-9 + 1 * i).into_verbose()],
            vec![(3 + 3 * i).into_verbose(), (-2 + 4 * i).into_verbose()],
        ]);
        mat.pprint();
        mat.clone().gram_schmidt_col_orthonormalization().pprint();
    }
}
