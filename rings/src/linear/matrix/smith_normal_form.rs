use super::*;

impl<RS: BezoutDomainSignature, RSB : BorrowedStructure<RS>> MatrixStructure<RS, RSB> {
    //return (u, s, v, k) such that self = usv and s is in smith normal form (with diagonal entries their favorite associates) and u, v are invertible and k is the number of non-zero elements in the diagonal of s
    pub fn smith_algorithm(
        &self,
        mut m: Matrix<RS::Set>,
    ) -> (Matrix<RS::Set>, Matrix<RS::Set>, Matrix<RS::Set>, usize) {
        let mut u = self.ident(m.rows());
        let mut v = self.ident(m.cols());

        let mut n = 0;
        'inductive_loop: while n < m.rows() && n < m.cols() {
            //search for a non-zero element to make the new starting point for (n, n)
            //having a non-zero element is necessary later in the algorithm
            'search_for_nonzero_element: {
                if self.ring().equal(m.at(n, n).unwrap(), &self.ring().zero()) {
                    //search the first row to start with
                    for c in n + 1..m.cols() {
                        if !self.ring().equal(m.at(n, c).unwrap(), &self.ring().zero()) {
                            //swap column n and column c
                            let col_opp = ElementaryOpp::new_col_opp(
                                self.ring().clone(),
                                ElementaryOppType::Swap(n, c),
                            );
                            col_opp.apply(&mut m);
                            col_opp.apply(&mut v);

                            break 'search_for_nonzero_element;
                        }
                    }
                    //search all the rows below row n
                    for r in n + 1..m.rows() {
                        for c in n..m.cols() {
                            if !self.ring().equal(m.at(r, c).unwrap(), &self.ring().zero()) {
                                //swap column n and column c
                                let col_opp = ElementaryOpp::new_col_opp(
                                    self.ring().clone(),
                                    ElementaryOppType::Swap(n, c),
                                );
                                col_opp.apply(&mut m);
                                col_opp.apply(&mut v);

                                //swap row n and row r
                                let row_opp = ElementaryOpp::new_col_opp(
                                    self.ring().clone(),
                                    ElementaryOppType::Swap(n, r),
                                );
                                row_opp.apply(&mut m);
                                row_opp.apply(&mut u);

                                break 'search_for_nonzero_element;
                            }
                        }
                    }
                    //everything remaining is zero. The algorithm can terminate here
                    break 'inductive_loop;
                }
            }

            //turn (n, n) into its favorite associate
            let (unit, _assoc) = self.ring().factor_fav_assoc(m.at(n, n).unwrap());
            let row_opp = ElementaryOpp::new_row_opp(
                self.ring().clone(),
                ElementaryOppType::UnitMul {
                    row: n,
                    unit: self.ring().inv(&unit).unwrap(),
                },
            );
            row_opp.apply(&mut m);
            row_opp.apply(&mut u);

            let mut first = true;
            let mut all_divisible;
            'zero_first_row_and_column_loop: loop {
                //replace the first row (a0, a1, ..., ak) with (gcd, 0, ..., 0). Might mess up the first column in the process
                all_divisible = true;
                for c in n + 1..m.cols() {
                    let a = m.at(n, n).unwrap();
                    let b = m.at(n, c).unwrap();
                    match self.ring().div(b, a) {
                        Ok(q) => {
                            //b is a multiple of a
                            //replace (a, b) with (a, 0) by subtracting a multiple of a from b
                            let col_opp = ElementaryOpp::new_col_opp(
                                self.ring().clone(),
                                ElementaryOppType::AddRowMul {
                                    i: c,
                                    j: n,
                                    x: self.ring().neg(&q),
                                },
                            );
                            col_opp.apply(&mut m);
                            col_opp.apply(&mut v);
                        }
                        Err(RingDivisionError::NotDivisible) => {
                            all_divisible = false;
                            //b is not a multiple of a
                            //replace (a, b) with (gcd, 0)
                            let (d, x, y) = self.ring().xgcd(a, b);
                            debug_assert!(
                                self.ring().equal(
                                    &self
                                        .ring()
                                        .add(&self.ring().mul(&x, a), &self.ring().mul(&y, b)),
                                    &d
                                )
                            );
                            let col_opp = ElementaryOpp::new_col_opp(
                                self.ring().clone(),
                                ElementaryOppType::TwoInv {
                                    i: n,
                                    j: c,
                                    a: x,
                                    b: y,
                                    c: self.ring().neg(&self.ring().div(b, &d).unwrap()),
                                    d: self.ring().div(a, &d).unwrap(),
                                },
                            );
                            col_opp.apply(&mut m);
                            col_opp.apply(&mut v);
                        }
                        Err(RingDivisionError::DivideByZero) => {
                            //swap a and b
                            //a=0 so this does have the effect of (a, b) -> (gcd(a, b), 0)
                            let col_opp = ElementaryOpp::new_col_opp(
                                self.ring().clone(),
                                ElementaryOppType::Swap(n, c),
                            );
                            col_opp.apply(&mut m);
                            col_opp.apply(&mut v);
                        }
                    }
                }
                if all_divisible && !first {
                    break 'zero_first_row_and_column_loop;
                }
                first = false;

                //replace the first column (a0, a1, ..., ak) with (gcd, 0, ..., 0). Might mess up the first row in the process
                all_divisible = true;
                for r in n + 1..m.rows() {
                    let a = m.at(n, n).unwrap();
                    let b = m.at(r, n).unwrap();
                    match self.ring().div(b, a) {
                        Ok(q) => {
                            //b is a multiple of a
                            //replace (a, b) with (a, 0) by subtracting a multiple of a from b
                            let col_opp = ElementaryOpp::new_row_opp(
                                self.ring().clone(),
                                ElementaryOppType::AddRowMul {
                                    i: r,
                                    j: n,
                                    x: self.ring().neg(&q),
                                },
                            );
                            col_opp.apply(&mut m);
                            col_opp.apply(&mut u);
                        }
                        Err(RingDivisionError::NotDivisible) => {
                            all_divisible = false;
                            //b is not a multiple of a
                            //replace (a, b) with (gcd, 0)
                            let (d, x, y) = self.ring().xgcd(a, b);
                            debug_assert!(
                                self.ring().equal(
                                    &self
                                        .ring()
                                        .add(&self.ring().mul(&x, a), &self.ring().mul(&y, b)),
                                    &d
                                )
                            );
                            let row_opp = ElementaryOpp::new_row_opp(
                                self.ring().clone(),
                                ElementaryOppType::TwoInv {
                                    i: n,
                                    j: r,
                                    a: x,
                                    b: y,
                                    c: self.ring().neg(&self.ring().div(b, &d).unwrap()),
                                    d: self.ring().div(a, &d).unwrap(),
                                },
                            );
                            row_opp.apply(&mut m);
                            row_opp.apply(&mut u);
                        }
                        Err(RingDivisionError::DivideByZero) => {
                            //swap a and b
                            //a=0 so this does have the effect of (a, b) -> (gcd(a, b), 0)
                            let col_opp = ElementaryOpp::new_row_opp(
                                self.ring().clone(),
                                ElementaryOppType::Swap(n, r),
                            );
                            col_opp.apply(&mut m);
                            col_opp.apply(&mut u);
                        }
                    }
                }
                if all_divisible {
                    break 'zero_first_row_and_column_loop;
                }
            }
            //now the first row and the first column are all zero except the top left element at (n, n) which is non-zero
            debug_assert!(!self.ring().equal(m.at(n, n).unwrap(), &self.ring().zero()));
            //some more fiddling is needed now to make sure the top left element divides everything else
            for r in n + 1..m.rows() {
                //row(n) = row(n) + row(r)
                let row_opp = ElementaryOpp::new_row_opp(
                    self.ring().clone(),
                    ElementaryOppType::AddRowMul {
                        i: n,
                        j: r,
                        x: self.ring().one(),
                    },
                );
                row_opp.apply(&mut m);
                row_opp.apply(&mut u);

                //row(n) goes from (g, a1, a2, ..., an) to (gcd, 0, 0, ..., 0) by applying col opps
                for c in n + 1..m.cols() {
                    let a = m.at(n, n).unwrap();
                    let b = m.at(n, c).unwrap();
                    // println!("a={:?} b={:?}", a, b);
                    //if a=0 and b=0 then there is nothing to do and the following step would fail when dividing by g=0
                    if !self.ring().is_zero(a) || !self.ring().is_zero(b) {
                        //b might not be a multiple of a
                        //replace (a, b) with (gcd, 0) to fix this
                        let (g, x, y) = self.ring().xgcd(a, b);
                        debug_assert!(
                            self.ring().equal(
                                &self
                                    .ring()
                                    .add(&self.ring().mul(&x, a), &self.ring().mul(&y, b)),
                                &g
                            )
                        );
                        let col_opp = ElementaryOpp::new_col_opp(
                            self.ring().clone(),
                            ElementaryOppType::TwoInv {
                                i: n,
                                j: c,
                                a: x,
                                b: y,
                                c: self.ring().neg(&self.ring().div(b, &g).unwrap()),
                                d: self.ring().div(a, &g).unwrap(),
                            },
                        );
                        col_opp.apply(&mut m);
                        col_opp.apply(&mut v);
                    }
                }
            }

            //fix the first column
            for fix_r in n + 1..m.rows() {
                let a = m.at(n, n).unwrap();
                let b = m.at(fix_r, n).unwrap();
                let q = self.ring().div(b, a).unwrap();
                let col_opp = ElementaryOpp::new_row_opp(
                    self.ring().clone(),
                    ElementaryOppType::AddRowMul {
                        i: fix_r,
                        j: n,
                        x: self.ring().neg(&q),
                    },
                );
                col_opp.apply(&mut m);
                col_opp.apply(&mut u);
            }

            if self.ring().equal(m.at(n, n).unwrap(), &self.ring().zero()) {
                //the bottom right submatrix is all zero
                break 'inductive_loop;
            }
            n += 1;
        }

        (u, m, v, n)
    }
}

impl<R: MetaType> Matrix<R>
where
    R::Signature: BezoutDomainSignature,
{
    pub fn smith_algorithm(&self) -> (Self, Self, Self, usize) {
        Self::structure().smith_algorithm(self.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_smith_algorithm() {
        {
            let a = Matrix::<Integer>::from_rows(vec![
                vec![Integer::from(2), Integer::from(4), Integer::from(4)],
                vec![Integer::from(-6), Integer::from(6), Integer::from(12)],
                vec![Integer::from(10), Integer::from(4), Integer::from(16)],
            ]);
            let (u, s, v, k) = a.clone().smith_algorithm();
            assert_eq!(s, Matrix::mul(&Matrix::mul(&u, &a).unwrap(), &v).unwrap());
            assert_eq!(k, 3);
            assert_eq!(
                s,
                Matrix::from_rows(vec![
                    vec![Integer::from(2), Integer::from(0), Integer::from(0)],
                    vec![Integer::from(0), Integer::from(2), Integer::from(0)],
                    vec![Integer::from(0), Integer::from(0), Integer::from(156)]
                ])
            );

            let a = Matrix::<Integer>::from_rows(vec![
                vec![
                    Integer::from(-6),
                    Integer::from(111),
                    Integer::from(-36),
                    Integer::from(6),
                ],
                vec![
                    Integer::from(5),
                    Integer::from(-672),
                    Integer::from(210),
                    Integer::from(74),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(-255),
                    Integer::from(81),
                    Integer::from(24),
                ],
                vec![
                    Integer::from(-7),
                    Integer::from(255),
                    Integer::from(-81),
                    Integer::from(-10),
                ],
            ]);
            let (u, s, v, k) = a.clone().smith_algorithm();
            assert_eq!(s, Matrix::mul(&Matrix::mul(&u, &a).unwrap(), &v).unwrap());
            assert_eq!(k, 3);
            assert_eq!(
                s,
                Matrix::from_rows(vec![
                    vec![
                        Integer::from(1),
                        Integer::from(0),
                        Integer::from(0),
                        Integer::from(0)
                    ],
                    vec![
                        Integer::from(0),
                        Integer::from(3),
                        Integer::from(0),
                        Integer::from(0)
                    ],
                    vec![
                        Integer::from(0),
                        Integer::from(0),
                        Integer::from(21),
                        Integer::from(0)
                    ],
                    vec![
                        Integer::from(0),
                        Integer::from(0),
                        Integer::from(0),
                        Integer::from(0)
                    ]
                ])
            );
        }

        {
            //used to cause a divide by zero
            let a = Matrix::<Integer>::from_rows(vec![
                vec![
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(1),
                ],
            ]);
            let (_u, _s, _v, k) = a.clone().smith_algorithm();
            assert_eq!(k, 1);
        }
    }
}
