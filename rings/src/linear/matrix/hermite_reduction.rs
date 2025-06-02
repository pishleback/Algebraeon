use crate::linear::{
    finitely_free_affine::FinitelyFreeSubmoduleAffineSubset,
    finitely_free_module::FinitelyFreeModuleStructure,
    finitely_free_submodule::FinitelyFreeSubmodule,
};

use super::*;

/// Rings for which hermite normal forms can be computed
pub trait HermiteAlgorithmSignature: BezoutDomainSignature {}
impl<Ring: BezoutDomainSignature> HermiteAlgorithmSignature for Ring {}

/// Rings for which reduced hermite normal forms can be computed
pub trait ReducedHermiteAlgorithmSignature:
    HermiteAlgorithmSignature + EuclideanDomainSignature + FavoriteAssociateSignature
{
}
impl<Ring: HermiteAlgorithmSignature + EuclideanDomainSignature + FavoriteAssociateSignature>
    ReducedHermiteAlgorithmSignature for Ring
{
}

/// Rings for which reduced hermite normal forms can be computed and are unique
pub trait UniqueReducedHermiteAlgorithmSignature: ReducedHermiteAlgorithmSignature {}
impl UniqueReducedHermiteAlgorithmSignature for IntegerCanonicalStructure {}
impl<Field: FieldSignature> UniqueReducedHermiteAlgorithmSignature for Field {}

impl<Ring: HermiteAlgorithmSignature, RingB: BorrowedStructure<Ring>> MatrixStructure<Ring, RingB> {
    /// Return (H, U, u_det, pivots) such that
    /// - H is in row hermite normal form, meaning
    /// - U is invertible
    /// - UM=H
    /// - u_det is the determinant of u
    /// - pivots[r] is the column of the rth pivot
    pub fn row_hermite_algorithm(
        &self,
        mut m: Matrix<Ring::Set>,
    ) -> (Matrix<Ring::Set>, Matrix<Ring::Set>, Ring::Set, Vec<usize>) {
        //build up U by applying row opps to the identity as we go
        let mut u = self.ident(m.rows());
        let mut u_det = self.ring().one();
        let mut pivs = vec![];

        let (mut pr, mut pc) = (0, 0);
        'pivot_loop: while pr < m.rows() {
            //find the next pivot row
            //the next pivot row is the next row with a non-zero element below the previous pivot
            'next_pivot_loop: loop {
                if pc == m.cols() {
                    break 'pivot_loop;
                }

                for r in pr..m.rows() {
                    if !self.ring().equal(m.at(r, pc).unwrap(), &self.ring().zero()) {
                        break 'next_pivot_loop;
                    }
                }

                pc += 1;
            }
            pivs.push(pc);

            if pr + 1 < m.rows() {
                //reduce everything below the pivot to zero
                for r in pr + 1..m.rows() {
                    let a = m.at(pr, pc).unwrap();
                    let b = m.at(r, pc).unwrap();
                    //if a=0 and b=0 there is nothing to do. The reduction step would fail because d=0 and we divide by d, so just skip it in this case
                    if !self.ring().equal(a, &self.ring().zero())
                        || !self.ring().equal(b, &self.ring().zero())
                    {
                        let (d, x, y) = self.ring().xgcd(a, b);
                        debug_assert!(
                            self.ring().equal(
                                &self
                                    .ring()
                                    .add(&self.ring().mul(&x, a), &self.ring().mul(&y, b)),
                                &d
                            )
                        );
                        // perform the following row opps on self
                        // / x  -b/d \
                        // \ y   a/d /
                        let row_opp = ElementaryOpp::new_row_opp(
                            self.ring().clone(),
                            ElementaryOppType::TwoInv {
                                i: pr,
                                j: r,
                                a: x,
                                b: y,
                                //TODO: compute b/d and a/d at the same time d is computed?
                                c: self.ring().neg(&self.ring().div(b, &d).unwrap()),
                                d: self.ring().div(a, &d).unwrap(),
                            },
                        );
                        //this will implicitly put the pivot into fav assoc form because that is what gcd does
                        row_opp.apply(&mut m);
                        row_opp.apply(&mut u);
                        self.ring().mul_mut(&mut u_det, &row_opp.det());
                    }
                }
            } else {
                //explicitly put the pivot into fav assoc form
                let (unit, _assoc) = self.ring().factor_fav_assoc(m.at(pr, pc).unwrap());
                let row_opp = ElementaryOpp::new_row_opp(
                    self.ring().clone(),
                    ElementaryOppType::UnitMul {
                        row: pr,
                        unit: self.ring().inv(&unit).unwrap(),
                    },
                );
                //this will implicitly put the pivot into fav assoc form because that is what the gcd returns
                row_opp.apply(&mut m);
                row_opp.apply(&mut u);
                self.ring().mul_mut(&mut u_det, &row_opp.det());
            }

            //should have eliminated everything below the pivot
            for r in pr + 1..m.rows() {
                debug_assert!(self.ring().equal(m.at(r, pc).unwrap(), &self.ring().zero()));
            }
            pr += 1;
        }

        if m.rows() <= 4 {
            debug_assert!(self.ring().equal(&self.det_naive(&u).unwrap(), &u_det));
        }

        (m, u, u_det, pivs)
    }

    pub fn col_hermite_algorithm(
        &self,
        a: Matrix<Ring::Set>,
    ) -> (Matrix<Ring::Set>, Matrix<Ring::Set>, Ring::Set, Vec<usize>) {
        let (rh, ru, u_det, pivs) = self.row_hermite_algorithm(a.transpose());
        (rh.transpose(), ru.transpose(), u_det, pivs)
    }

    fn det_hermite(&self, a: Matrix<Ring::Set>) -> Ring::Set {
        let n = a.rows();
        debug_assert_eq!(n, a.cols());
        let (h, _u, u_det, _pivs) = self.row_hermite_algorithm(a);
        //h = u * self, we know det(u), and h is upper triangular
        let mut h_det = self.ring().one();
        for i in 0..n {
            self.ring().mul_mut(&mut h_det, h.at(i, i).unwrap());
        }
        self.ring().div(&h_det, &u_det).unwrap()
    }

    pub fn det(&self, a: Matrix<Ring::Set>) -> Result<Ring::Set, MatOppErr> {
        let n = a.rows();
        if n != a.cols() {
            Err(MatOppErr::NotSquare)
        } else if n <= 3 {
            //for speed
            Ok(self.det_naive(&a).unwrap())
        } else {
            Ok(self.det_hermite(a))
        }
    }

    pub fn rank(&self, a: Matrix<Ring::Set>) -> usize {
        let (_h, _u, _u_det, pivs) = self.row_hermite_algorithm(a);
        pivs.len()
    }
}

impl<Ring: ReducedHermiteAlgorithmSignature, RingB: BorrowedStructure<Ring>>
    MatrixStructure<Ring, RingB>
{
    /// Returns (H, U, u_det, pivots) such that
    /// - H is in row reduced hermite normal form, meaning entries above pivots have euclidean norm strictly less than the pivot
    /// - U is invertible
    /// - UM=H
    /// - pivots[r] is the column of the rth pivot and pivots.len() == rank(A)
    pub fn row_reduced_hermite_algorithm(
        &self,
        m: Matrix<Ring::Set>,
    ) -> (Matrix<Ring::Set>, Matrix<Ring::Set>, Ring::Set, Vec<usize>) {
        let (mut h, mut u, u_det, pivs) = self.row_hermite_algorithm(m);

        for (pr, pc) in pivs.iter().enumerate() {
            for r in 0..pr {
                //reduce h[r, pc] so that it has norm less than h[pr, pc]
                let a = h.at(r, *pc).unwrap();
                let b = h.at(pr, *pc).unwrap();
                //a = b*q + r
                let q = self.ring().quo(a, b).unwrap();
                let row_opp = ElementaryOpp::new_row_opp(
                    self.ring().clone(),
                    ElementaryOppType::AddRowMul {
                        i: r,
                        j: pr,
                        x: self.ring().neg(&q),
                    },
                );
                row_opp.apply(&mut h);
                row_opp.apply(&mut u);
                // adding a multiple of a row does not change the determinant so no need to update u_det
            }
        }

        (h, u, u_det, pivs)
    }

    pub fn row_reduced_hermite_normal_form(&self, m: Matrix<Ring::Set>) -> Matrix<Ring::Set> {
        self.row_reduced_hermite_algorithm(m).0
    }

    pub fn col_reduced_hermite_algorithm(
        &self,
        m: Matrix<Ring::Set>,
    ) -> (Matrix<Ring::Set>, Matrix<Ring::Set>, Ring::Set, Vec<usize>) {
        let (rh, ru, ru_det, pivs) = self.row_reduced_hermite_algorithm(m.transpose());
        (rh.transpose(), ru.transpose(), ru_det, pivs)
    }

    pub fn col_reduced_hermite_normal_form(&self, m: Matrix<Ring::Set>) -> Matrix<Ring::Set> {
        self.col_reduced_hermite_algorithm(m).0
    }

    pub fn inv(&self, a: Matrix<Ring::Set>) -> Result<Matrix<Ring::Set>, MatOppErr> {
        let n = a.rows();
        if n == a.cols() {
            let (h, u, _u_det, _pivs) = self.row_reduced_hermite_algorithm(a);
            //h = u*a
            if self.equal(&h, &self.ident(n)) {
                Ok(u)
            } else {
                Err(MatOppErr::Singular)
            }
        } else {
            Err(MatOppErr::NotSquare)
        }
    }

    pub fn row_span(&self, matrix: Matrix<Ring::Set>) -> FinitelyFreeSubmodule<Ring::Set> {
        FinitelyFreeModuleStructure::<Ring, _>::new(self.ring(), matrix.cols())
            .into_submodules()
            .matrix_row_span(matrix)
    }

    pub fn col_span(&self, matrix: Matrix<Ring::Set>) -> FinitelyFreeSubmodule<Ring::Set> {
        FinitelyFreeModuleStructure::<Ring, _>::new(self.ring(), matrix.rows())
            .into_submodules()
            .matrix_col_span(matrix)
    }

    pub fn row_kernel(&self, matrix: Matrix<Ring::Set>) -> FinitelyFreeSubmodule<Ring::Set> {
        FinitelyFreeModuleStructure::<Ring, _>::new(self.ring(), matrix.rows())
            .into_submodules()
            .matrix_row_kernel(matrix)
    }

    pub fn col_kernel(&self, matrix: Matrix<Ring::Set>) -> FinitelyFreeSubmodule<Ring::Set> {
        FinitelyFreeModuleStructure::<Ring, _>::new(self.ring(), matrix.cols())
            .into_submodules()
            .matrix_col_kernel(matrix)
    }

    pub fn row_affine_span(
        &self,
        matrix: Matrix<Ring::Set>,
    ) -> FinitelyFreeSubmoduleAffineSubset<Ring::Set> {
        let span = (0..matrix.rows())
            .map(|r| matrix.get_row(r))
            .collect::<Vec<_>>();
        FinitelyFreeModuleStructure::<Ring, _>::new(self.ring(), matrix.cols())
            .affine_subsets()
            .from_affine_span(span.iter().collect())
    }

    pub fn col_affine_span(
        &self,
        matrix: Matrix<Ring::Set>,
    ) -> FinitelyFreeSubmoduleAffineSubset<Ring::Set> {
        self.row_affine_span(matrix.transpose())
    }

    pub fn row_solve(
        &self,
        matrix: Matrix<Ring::Set>,
        y: &Vec<Ring::Set>,
    ) -> Option<Vec<Ring::Set>> {
        let submodules = FinitelyFreeModuleStructure::<Ring, _>::new(self.ring(), matrix.cols())
            .into_submodules();
        let (row_span_submodule, basis_in_terms_of_matrix_rows) =
            submodules.matrix_row_span_and_basis(matrix);
        let (offset, y_reduced) = submodules.reduce_element(&row_span_submodule, y);
        if y_reduced.iter().all(|v| self.ring().is_zero(v)) {
            Some(
                self.mul(
                    &Matrix::from_rows(vec![offset]),
                    &basis_in_terms_of_matrix_rows,
                )
                .unwrap()
                .get_row(0),
            )
        } else {
            None
        }
    }

    pub fn col_solve(
        &self,
        matrix: Matrix<Ring::Set>,
        y: &Vec<Ring::Set>,
    ) -> Option<Vec<Ring::Set>> {
        self.row_solve(matrix.transpose(), y)
    }

    pub fn row_solution_set(
        &self,
        matrix: Matrix<Ring::Set>,
        y: &Vec<Ring::Set>,
    ) -> FinitelyFreeSubmoduleAffineSubset<Ring::Set> {
        let module = FinitelyFreeModuleStructure::<Ring, _>::new(self.ring(), matrix.rows());
        match self.row_solve(matrix.clone(), y) {
            Some(offset) => FinitelyFreeSubmoduleAffineSubset::NonEmpty(
                module
                    .cosets()
                    .from_offset_and_submodule(&offset, self.row_kernel(matrix)),
            ),
            None => FinitelyFreeSubmoduleAffineSubset::Empty,
        }
    }

    pub fn col_solution_set(
        &self,
        matrix: Matrix<Ring::Set>,
        y: &Vec<Ring::Set>,
    ) -> FinitelyFreeSubmoduleAffineSubset<Ring::Set> {
        self.row_solution_set(matrix.transpose(), y)
    }
}

impl<R: MetaType> Matrix<R>
where
    R::Signature: BezoutDomainSignature,
{
    pub fn row_hermite_algorithm(&self) -> (Self, Self, R, Vec<usize>) {
        Self::structure().row_hermite_algorithm(self.clone())
    }

    pub fn col_hermite_algorithm(&self) -> (Self, Self, R, Vec<usize>) {
        Self::structure().col_hermite_algorithm(self.clone())
    }

    pub fn det(&self) -> Result<R, MatOppErr> {
        Self::structure().det(self.clone())
    }

    pub fn rank(&self) -> usize {
        Self::structure().rank(self.clone())
    }
}

impl<R: MetaType> Matrix<R>
where
    R::Signature: EuclideanDomainSignature + BezoutDomainSignature + FavoriteAssociateSignature,
{
    pub fn row_reduced_hermite_algorithm(&self) -> (Self, Self, R, Vec<usize>) {
        Self::structure().row_reduced_hermite_algorithm(self.clone())
    }

    pub fn row_reduced_hermite_normal_form(&self) -> Self {
        Self::structure().row_reduced_hermite_normal_form(self.clone())
    }

    pub fn col_reduced_hermite_algorithm(&self) -> (Self, Self, R, Vec<usize>) {
        Self::structure().col_reduced_hermite_algorithm(self.clone())
    }

    pub fn col_reduced_hermite_normal_form(&self) -> Self {
        Self::structure().col_reduced_hermite_normal_form(self.clone())
    }

    pub fn inv(&self) -> Result<Matrix<R>, MatOppErr> {
        Self::structure().inv(self.clone())
    }

    pub fn row_span(self) -> FinitelyFreeSubmodule<R> {
        Self::structure().row_span(self)
    }

    pub fn col_span(self) -> FinitelyFreeSubmodule<R> {
        Self::structure().col_span(self)
    }

    pub fn row_kernel(self) -> FinitelyFreeSubmodule<R> {
        Self::structure().row_kernel(self)
    }

    pub fn col_kernel(self) -> FinitelyFreeSubmodule<R> {
        Self::structure().col_kernel(self)
    }

    pub fn row_affine_span(self) -> FinitelyFreeSubmoduleAffineSubset<R> {
        Self::structure().row_affine_span(self)
    }

    pub fn col_affine_span(self) -> FinitelyFreeSubmoduleAffineSubset<R> {
        Self::structure().col_affine_span(self)
    }

    pub fn row_solve(self, y: &Vec<R>) -> Option<Vec<R>> {
        Self::structure().row_solve(self, y)
    }

    pub fn col_solve(self, y: &Vec<R>) -> Option<Vec<R>> {
        Self::structure().col_solve(self, y)
    }

    pub fn row_solution_set(self, y: &Vec<R>) -> FinitelyFreeSubmoduleAffineSubset<R> {
        Self::structure().row_solution_set(self, y)
    }

    pub fn col_solution_set(self, y: &Vec<R>) -> FinitelyFreeSubmoduleAffineSubset<R> {
        Self::structure().col_solution_set(self, y)
    }
}

#[cfg(test)]
mod tests {
    use crate::linear::finitely_free_module::RingToFinitelyFreeModuleSignature;

    use super::*;

    #[test]
    fn hermite_algorithm() {
        for a in [
            Matrix::from_rows(vec![
                vec![
                    Integer::from(2),
                    Integer::from(3),
                    Integer::from(6),
                    Integer::from(2),
                ],
                vec![
                    Integer::from(5),
                    Integer::from(6),
                    Integer::from(1),
                    Integer::from(6),
                ],
                vec![
                    Integer::from(8),
                    Integer::from(3),
                    Integer::from(1),
                    Integer::from(1),
                ],
            ]),
            Matrix::from_rows(vec![
                vec![
                    Integer::from(2),
                    Integer::from(3),
                    Integer::from(6),
                    Integer::from(2),
                ],
                vec![
                    Integer::from(5),
                    Integer::from(6),
                    Integer::from(1),
                    Integer::from(6),
                ],
                vec![
                    Integer::from(8),
                    Integer::from(3),
                    Integer::from(1),
                    Integer::from(1),
                ],
            ])
            .transpose(),
            Matrix::from_rows(vec![
                vec![
                    Integer::from(2),
                    Integer::from(3),
                    Integer::from(5),
                    Integer::from(2),
                ],
                vec![
                    Integer::from(5),
                    Integer::from(6),
                    Integer::from(11),
                    Integer::from(6),
                ],
                vec![
                    Integer::from(8),
                    Integer::from(3),
                    Integer::from(11),
                    Integer::from(1),
                ],
            ]),
            Matrix::from_rows(vec![
                vec![Integer::from(0), Integer::from(4), Integer::from(4)],
                vec![Integer::from(1), Integer::from(6), Integer::from(12)],
                vec![Integer::from(1), Integer::from(4), Integer::from(16)],
            ]),
        ] {
            println!();
            println!("hermite reduced row algorithm for");
            a.pprint();
            let (h, u, _u_det, pivs) = a.clone().row_reduced_hermite_algorithm();
            println!("H =");
            h.pprint();
            println!("pivs = {:?}", pivs);
            println!("U =");
            u.pprint();
            assert_eq!(h, Matrix::mul(&u, &a).unwrap());

            //trace the boundary of zeros and check that everything under is zero
            let mut rz = 0;
            for cz in 0..h.cols() {
                if let Some(cp) = pivs.get(rz) {
                    if cp == &cz {
                        rz += 1;
                    }
                }
                for r in rz..h.rows() {
                    assert_eq!(h.at(r, cz).unwrap(), &Integer::zero());
                }
            }

            //check pivot columns
            for (pr, pc) in pivs.iter().enumerate() {
                assert!(h.at(pr, *pc).unwrap() != &Integer::zero());
                for r in 0..h.rows() {
                    #[allow(clippy::comparison_chain)]
                    if r > pr {
                        assert_eq!(h.at(r, *pc).unwrap(), &Integer::zero());
                    } else if r == pr {
                        let (_unit, assoc) = Integer::factor_fav_assoc(h.at(r, *pc).unwrap());
                        assert_eq!(&assoc, h.at(r, *pc).unwrap());
                    } else {
                        assert!(
                            Integer::norm(h.at(r, *pc).unwrap())
                                < Integer::norm(h.at(pr, *pc).unwrap())
                        );
                    }
                }
            }

            println!();
            println!("hermite reduced col algorithm for");
            a.pprint();
            let (h, u, _u_det, pivs) = a.clone().col_reduced_hermite_algorithm();
            println!("H =");
            h.pprint();
            println!("pivs = {:?}", pivs);
            println!("U =");
            u.pprint();

            //trace the boundary of zeros and check that everything to the right is zero
            let mut cz = 0;
            for rz in 0..h.rows() {
                if let Some(rp) = pivs.get(cz) {
                    if rp == &rz {
                        cz += 1;
                    }
                }
                for c in cz..h.cols() {
                    assert_eq!(h.at(rz, c).unwrap(), &Integer::zero());
                }
            }

            //check the pivot rows
            assert_eq!(h, Matrix::mul(&a, &u).unwrap());
            for (pc, pr) in pivs.iter().enumerate() {
                assert!(h.at(*pr, pc).unwrap() != &Integer::zero());
                for c in 0..h.cols() {
                    #[allow(clippy::comparison_chain)]
                    if c > pc {
                        assert_eq!(h.at(*pr, c).unwrap(), &Integer::zero());
                    } else if c == pc {
                        let (_unit, assoc) = Integer::factor_fav_assoc(h.at(*pr, c).unwrap());
                        assert_eq!(&assoc, h.at(*pr, c).unwrap());
                    } else {
                        assert!(
                            Integer::norm(h.at(*pr, c).unwrap())
                                < Integer::norm(h.at(*pr, pc).unwrap())
                        );
                    }
                }
            }
        }

        {
            //integer reduced hermite normal form is unique, so we can fully check an example
            let a = Matrix::<Integer>::from_rows(vec![
                vec![
                    Integer::from(2),
                    Integer::from(3),
                    Integer::from(6),
                    Integer::from(2),
                ],
                vec![
                    Integer::from(5),
                    Integer::from(6),
                    Integer::from(1),
                    Integer::from(6),
                ],
                vec![
                    Integer::from(8),
                    Integer::from(3),
                    Integer::from(1),
                    Integer::from(1),
                ],
            ]);

            let expected_h = Matrix::from_rows(vec![
                vec![
                    Integer::from(1),
                    Integer::from(0),
                    Integer::from(50),
                    Integer::from(-11),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(3),
                    Integer::from(28),
                    Integer::from(-2),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(61),
                    Integer::from(-13),
                ],
            ]);

            let expected_u = Matrix::from_rows(vec![
                vec![Integer::from(9), Integer::from(-5), Integer::from(1)],
                vec![Integer::from(5), Integer::from(-2), Integer::from(0)],
                vec![Integer::from(11), Integer::from(-6), Integer::from(1)],
            ]);

            let (h, u, _u_det, pivs) = a.clone().row_reduced_hermite_algorithm();

            assert_eq!(h, expected_h);
            assert_eq!(u, expected_u);
            assert_eq!(pivs, vec![0, 1, 2]);
        }

        {
            //this one used to cause a dividion by zero error when replacing (a, b) with (gcd, 0) when (a, b) = (0, 0)
            let a = Matrix::<Integer>::from_rows(vec![
                vec![Integer::from(1), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(1), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(1)],
            ]);

            let expected_h = Matrix::from_rows(vec![
                vec![Integer::from(1), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(1), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(1)],
                vec![Integer::from(0), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(0)],
            ]);

            let expected_u = Matrix::from_rows(vec![
                vec![
                    Integer::from(1),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(1),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(1),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(1),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(-1),
                    Integer::from(0),
                    Integer::from(0),
                ],
            ]);

            let (h, u, _u_det, pivs) = a.clone().row_reduced_hermite_algorithm();

            assert_eq!(h, expected_h);
            assert_eq!(u, expected_u);
            assert_eq!(pivs, vec![0, 1, 2]);
        }
    }

    #[test]
    fn invert() {
        let a = Matrix::<Rational>::from_rows(vec![
            vec![Rational::from(2), Rational::from(4), Rational::from(4)],
            vec![Rational::from(-6), Rational::from(6), Rational::from(12)],
            vec![Rational::from(10), Rational::from(7), Rational::from(17)],
        ]);
        a.pprint();
        println!("{:?}", a.rank());
        let b = Matrix::inv(&a).unwrap();
        b.pprint();
        debug_assert_eq!(Matrix::mul(&a, &b).unwrap(), Matrix::ident(3));
        debug_assert_eq!(Matrix::mul(&b, &a).unwrap(), Matrix::ident(3));
    }

    #[test]
    fn span_and_kernel_rank() {
        let mat = Matrix::<Integer>::from_rows(vec![
            vec![Integer::from(1), Integer::from(1), Integer::from(1)],
            vec![Integer::from(1), Integer::from(2), Integer::from(1)],
            vec![Integer::from(1), Integer::from(1), Integer::from(1)],
            vec![Integer::from(1), Integer::from(1), Integer::from(1)],
        ]);

        assert_eq!(mat.clone().row_span().rank(), 2);
        assert_eq!(mat.clone().col_span().rank(), 2);
        assert_eq!(mat.clone().row_kernel().rank(), 2);
        assert_eq!(mat.clone().col_kernel().rank(), 1);
    }

    #[test]
    fn affine_span() {
        {
            let module = Integer::structure().into_free_module(2);

            //row affine span
            let lat1 = Matrix::<Integer>::from_rows(vec![
                vec![Integer::from(1), Integer::from(1)],
                vec![Integer::from(3), Integer::from(1)],
                vec![Integer::from(2), Integer::from(3)],
            ])
            .row_affine_span();

            let lat2 = FinitelyFreeSubmoduleAffineSubset::NonEmpty(
                module.cosets().from_offset_and_submodule(
                    &vec![Integer::from(2), Integer::from(3)],
                    Matrix::<Integer>::from_rows(vec![vec![1, 2], vec![-1, 2]]).row_span(),
                ),
            );

            println!("lat1 = {:?}", lat1);
            println!("lat2 = {:?}", lat2);

            assert!(module.affine_subsets().equal(&lat1, &lat2));
            assert!(module.affine_subsets().equal_slow(&lat1, &lat2));
        }

        {
            let module = Integer::structure().into_free_module(2);

            //column affine span
            let lat1 = Matrix::<Integer>::from_rows(vec![
                vec![Integer::from(1), Integer::from(3), Integer::from(2)],
                vec![Integer::from(1), Integer::from(1), Integer::from(3)],
            ])
            .col_affine_span();

            let lat2 = FinitelyFreeSubmoduleAffineSubset::NonEmpty(
                module.cosets().from_offset_and_submodule(
                    &vec![Integer::from(2), Integer::from(3)],
                    Matrix::<Integer>::from_rows(vec![vec![1, -1], vec![2, 2]]).col_span(),
                ),
            );

            println!("lat1 = {:?}", lat1);
            println!("lat2 = {:?}", lat2);

            assert!(module.affine_subsets().equal(&lat1, &lat2));
            assert!(module.affine_subsets().equal_slow(&lat1, &lat2));
        }
    }

    #[test]
    fn span_and_kernel_points() {
        let module = Integer::structure().into_free_module(4);

        let mat = Matrix::<Integer>::from_rows(vec![
            vec![
                Integer::from(1),
                Integer::from(1),
                Integer::from(0),
                Integer::from(0),
            ],
            vec![
                Integer::from(1),
                Integer::from(1),
                Integer::from(0),
                Integer::from(0),
            ],
            vec![
                Integer::from(0),
                Integer::from(0),
                Integer::from(3),
                Integer::from(5),
            ],
        ]);

        println!("matrix");
        mat.pprint();

        let k = mat.col_kernel();

        assert!(module.submodules().contains_element(
            &k,
            &vec![
                Integer::from(-1),
                Integer::from(1),
                Integer::from(5),
                Integer::from(-3)
            ]
        ));
    }

    #[test]
    fn test_row_solve() {
        let module = Integer::structure().into_free_module(3);

        let matrix =
            Matrix::<Integer>::from_rows(vec![vec![1, 0, 0], vec![1, 0, 1], vec![1, 1, 1]]);
        let x =
            matrix
                .clone()
                .row_solve(&vec![Integer::from(4), Integer::from(4), Integer::from(7)]);
        assert_eq!(
            x.unwrap(),
            vec![Integer::from(-3), Integer::from(3), Integer::from(4),]
        );

        let a = Matrix::<Integer>::from_rows(vec![vec![1, 0, 0], vec![1, 0, 1]])
            .row_solution_set(&vec![Integer::from(2), Integer::from(3), Integer::from(2)]);
        for s in module.affine_subsets().affine_basis(&a) {
            println!("{:?}", s);
        }
    }
}
