use super::*;

/// Rings for which hermite normal forms can be computed
pub trait HermiteAlgorithmSignature: BezoutDomainSignature {}
impl<Ring: BezoutDomainSignature> HermiteAlgorithmSignature for Ring {}

/// Rings for which reduced hermite normal forms can be computed
pub trait ReducedHermiteAlgorithmSignature:
    HermiteAlgorithmSignature + EuclideanDivisionSignature + FavoriteAssociateSignature
{
}
impl<Ring: HermiteAlgorithmSignature + EuclideanDivisionSignature + FavoriteAssociateSignature>
    ReducedHermiteAlgorithmSignature for Ring
{
}

/// Rings for which reduced hermite normal forms can be computed and are unique
pub trait UniqueReducedHermiteAlgorithmSignature: ReducedHermiteAlgorithmSignature {}
impl UniqueReducedHermiteAlgorithmSignature for IntegerCanonicalStructure {}
impl<Field: FieldSignature> UniqueReducedHermiteAlgorithmSignature for Field {}

impl<Ring: BezoutDomainSignature> MatrixStructure<Ring> {
    #[deprecated]
    pub fn row_span(&self, a: Matrix<Ring::Set>) -> LinearSubspace<Ring::Set> {
        LinearSubspaceStructure::new(self.ring().clone()).from_span(
            1,
            a.cols(),
            (0..a.rows())
                .map(|r| a.submatrix(vec![r], (0..a.cols()).collect()))
                .collect(),
        )
    }

    #[deprecated]
    pub fn col_span(&self, a: Matrix<Ring::Set>) -> LinearSubspace<Ring::Set> {
        LinearSubspaceStructure::new(self.ring().clone()).from_span(
            a.rows(),
            1,
            (0..a.cols())
                .map(|c| a.submatrix((0..a.rows()).collect(), vec![c]))
                .collect(),
        )
    }

    #[deprecated]
    pub fn row_affine_span(&self, a: Matrix<Ring::Set>) -> AffineSubspace<Ring::Set> {
        let affine_lattice_structure = AffineSubspaceStructure::new(self.ring().clone());
        if a.rows() == 0 {
            affine_lattice_structure.empty(1, a.cols())
        } else {
            let offset = a.get_row(0);

            let b = Matrix::construct(a.rows() - 1, a.cols(), |r, c| {
                self.ring().add(
                    &self.ring().neg(offset.at(0, c).unwrap()),
                    a.at(r + 1, c).unwrap(),
                )
            });

            let linlat = self.row_span(b);

            affine_lattice_structure.from_offset_and_linear_lattice(1, a.cols(), offset, linlat)
        }
    }

    #[deprecated]
    pub fn col_affine_span(&self, a: Matrix<Ring::Set>) -> AffineSubspace<Ring::Set> {
        let affine_lattice_structure = AffineSubspaceStructure::new(self.ring().clone());
        if a.cols() == 0 {
            affine_lattice_structure.empty(a.rows(), 1)
        } else {
            let offset = a.get_col(0);

            let b = Matrix::construct(a.rows(), a.cols() - 1, |r, c| {
                self.ring().add(
                    &self.ring().neg(offset.at(r, 0).unwrap()),
                    a.at(r, c + 1).unwrap(),
                )
            });

            let linlat = self.col_span(b);

            affine_lattice_structure.from_offset_and_linear_lattice(a.rows(), 1, offset, linlat)
        }
    }

    #[deprecated]
    pub fn row_kernel(&self, a: Matrix<Ring::Set>) -> LinearSubspace<Ring::Set> {
        let (_h, u, _u_det, pivs) = self.row_hermite_algorithm(a);
        LinearSubspaceStructure::new(self.ring().clone()).from_basis(
            1,
            u.cols(),
            (pivs.len()..u.rows())
                .into_iter()
                .map(|r| u.submatrix(vec![r], (0..u.cols()).collect()))
                .collect(),
        )
    }

    #[deprecated]
    pub fn col_kernel(&self, a: Matrix<Ring::Set>) -> LinearSubspace<Ring::Set> {
        let (_h, u, _u_det, pivs) = self.col_hermite_algorithm(a);
        LinearSubspaceStructure::new(self.ring().clone()).from_basis(
            u.rows(),
            1,
            (pivs.len()..u.cols())
                .into_iter()
                .map(|c| u.submatrix((0..u.rows()).collect(), vec![c]))
                .collect(),
        )
    }

    #[deprecated]
    pub fn row_solve(
        &self,
        m: &Matrix<Ring::Set>,
        y: impl Borrow<Matrix<Ring::Set>>,
    ) -> Option<Matrix<Ring::Set>> {
        match self.col_solve(&m.transpose_ref(), &y.borrow().transpose_ref()) {
            Some(x) => Some(x.transpose()),
            None => None,
        }
    }

    #[deprecated]
    pub fn col_solve(
        &self,
        m: &Matrix<Ring::Set>,
        y: impl Borrow<Matrix<Ring::Set>>,
    ) -> Option<Matrix<Ring::Set>> {
        assert_eq!(y.borrow().rows(), m.rows());
        assert_eq!(y.borrow().cols(), 1);
        //the kernel of ext_mat is related to the solution
        let ext_mat = Matrix::join_cols(m.rows(), vec![y.borrow(), m]);
        //we are looking for a point in the column kernel where the first coordinate is 1
        let col_ker = self.col_kernel(ext_mat);

        let first_coords: Vec<&Ring::Set> = (0..LinearSubspaceStructure::new(self.ring().clone())
            .rank(&col_ker))
            .map(|basis_num| {
                LinearSubspaceStructure::new(self.ring().clone())
                    .basis_matrix_element(&col_ker, basis_num, 0, 0)
            })
            .collect();

        let (g, taps) = self.ring().xgcd_list(first_coords);

        if self.ring().is_unit(&g) {
            debug_assert!(self.ring().equal(&g, &self.ring().one()));
        }
        if self.ring().equal(&g, &self.ring().one()) {
            //there is a solution
            //it is given by -(sum(taps * col_ker.basis)) with the first coordinate (equal to 1) removed
            let mut ext_ans = self.zero(m.cols() + 1, 1);
            for basis_num in 0..LinearSubspaceStructure::new(self.ring().clone()).rank(&col_ker) {
                self.add_mut(
                    &mut ext_ans,
                    &self.mul_scalar_ref(
                        &LinearSubspaceStructure::new(self.ring().clone())
                            .basis_matrix(&col_ker, basis_num),
                        &taps[basis_num],
                    ),
                )
                .unwrap();
            }
            debug_assert!(
                self.ring()
                    .equal(ext_ans.at(0, 0).unwrap(), &self.ring().one())
            );
            let x = self.neg(ext_ans.submatrix((1..ext_ans.rows()).collect(), vec![0]));
            debug_assert!(self.equal(&self.mul(m, &x).unwrap(), y.borrow()));
            Some(x)
        } else {
            None //there is no solution
        }
    }

    #[deprecated]
    pub fn row_solution_lattice(
        &self,
        m: &Matrix<Ring::Set>,
        y: impl Borrow<Matrix<Ring::Set>>,
    ) -> AffineSubspace<Ring::Set> {
        match self.row_solve(m, y) {
            Some(x) => AffineSubspaceStructure::new(self.ring().clone())
                .from_offset_and_linear_lattice(1, m.rows(), x, self.row_kernel(m.clone())),
            None => AffineSubspaceStructure::new(self.ring().clone()).empty(1, m.rows()),
        }
    }

    #[deprecated]
    pub fn col_solution_lattice(
        &self,
        m: &Matrix<Ring::Set>,
        y: impl Borrow<Matrix<Ring::Set>>,
    ) -> AffineSubspace<Ring::Set> {
        match self.col_solve(m, y) {
            Some(x) => AffineSubspaceStructure::new(self.ring().clone())
                .from_offset_and_linear_lattice(m.cols(), 1, x, self.col_kernel(m.clone())),
            None => AffineSubspaceStructure::new(self.ring().clone()).empty(m.cols(), 1),
        }
    }
}

impl<Ring: HermiteAlgorithmSignature> MatrixStructure<Ring> {
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

impl<Ring: ReducedHermiteAlgorithmSignature> MatrixStructure<Ring> {
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
        if n != a.cols() {
            Err(MatOppErr::NotSquare)
        } else {
            let (h, u, _u_det, _pivs) = self.row_reduced_hermite_algorithm(a);
            //h = u*a
            if self.equal(&h, &self.ident(n)) {
                Ok(u)
            } else {
                Err(MatOppErr::Singular)
            }
        }
    }
}

impl<R: MetaType> Matrix<R>
where
    R::Signature: BezoutDomainSignature,
{
    pub fn row_span(&self) -> LinearSubspace<R> {
        Self::structure().row_span(self.clone())
    }

    pub fn col_span(&self) -> LinearSubspace<R> {
        Self::structure().col_span(self.clone())
    }

    pub fn row_affine_span(&self) -> AffineSubspace<R> {
        Self::structure().row_affine_span(self.clone())
    }

    pub fn col_affine_span(&self) -> AffineSubspace<R> {
        Self::structure().col_affine_span(self.clone())
    }

    pub fn row_kernel(&self) -> LinearSubspace<R> {
        Self::structure().row_kernel(self.clone())
    }

    pub fn col_kernel(&self) -> LinearSubspace<R> {
        Self::structure().col_kernel(self.clone())
    }

    pub fn row_solve(&self, y: impl Borrow<Self>) -> Option<Self> {
        Self::structure().row_solve(self, y)
    }

    pub fn col_solve(&self, y: impl Borrow<Self>) -> Option<Self> {
        Self::structure().col_solve(self, y)
    }

    pub fn row_solution_lattice(&self, y: impl Borrow<Self>) -> AffineSubspace<R> {
        Self::structure().row_solution_lattice(self, y)
    }

    pub fn col_solution_lattice(&self, y: impl Borrow<Self>) -> AffineSubspace<R> {
        Self::structure().col_solution_lattice(self, y)
    }

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
    R::Signature: EuclideanDivisionSignature + BezoutDomainSignature + FavoriteAssociateSignature,
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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hermite_algorithm() {
        for a in vec![
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
                match pivs.get(rz) {
                    Some(cp) => {
                        if cp == &cz {
                            rz += 1;
                        }
                    }
                    None => {}
                }
                for r in rz..h.rows() {
                    assert_eq!(h.at(r, cz).unwrap(), &Integer::zero());
                }
            }

            //check pivot columns
            for (pr, pc) in pivs.iter().enumerate() {
                assert!(h.at(pr, *pc).unwrap() != &Integer::zero());
                for r in 0..h.rows() {
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
                match pivs.get(cz) {
                    Some(rp) => {
                        if rp == &rz {
                            cz += 1;
                        }
                    }
                    None => {}
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
            //row affine span
            let lat1 = Matrix::<Integer>::from_rows(vec![
                vec![Integer::from(1), Integer::from(1)],
                vec![Integer::from(3), Integer::from(1)],
                vec![Integer::from(2), Integer::from(3)],
            ])
            .row_affine_span();

            let lat2 = AffineSubspace::from_offset_and_linear_lattice(
                1,
                2,
                Matrix::from_rows(vec![vec![Integer::from(2), Integer::from(3)]]),
                Matrix::from_rows(vec![
                    vec![Integer::from(1), Integer::from(2)],
                    vec![Integer::from(-1), Integer::from(2)],
                ])
                .row_span(),
            );

            lat1.pprint();
            lat2.pprint();

            assert_eq!(lat1, lat2);
        }

        {
            //column affine span
            let lat1 = Matrix::<Integer>::from_rows(vec![
                vec![Integer::from(1), Integer::from(3), Integer::from(2)],
                vec![Integer::from(1), Integer::from(1), Integer::from(3)],
            ])
            .col_affine_span();

            let lat2 = AffineSubspace::from_offset_and_linear_lattice(
                2,
                1,
                Matrix::from_rows(vec![vec![Integer::from(2)], vec![Integer::from(3)]]),
                Matrix::from_rows(vec![
                    vec![Integer::from(1), Integer::from(-1)],
                    vec![Integer::from(2), Integer::from(2)],
                ])
                .col_span(),
            );

            lat1.pprint();
            lat2.pprint();

            assert_eq!(&lat1, &lat2);
        }
    }

    #[test]
    fn span_and_kernel_points() {
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
        println!("kernel");
        k.pprint();

        assert!(k.contains_point(Matrix::from_rows(vec![
            vec![Integer::from(-1)],
            vec![Integer::from(1)],
            vec![Integer::from(5)],
            vec![Integer::from(-3)]
        ])));
    }
}
