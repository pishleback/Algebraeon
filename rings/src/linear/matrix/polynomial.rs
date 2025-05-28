use super::*;

impl<FS: FieldSignature, FSB: BorrowedStructure<FS>> MatrixStructure<FS, FSB> {
    pub fn presentation_matrix(
        &self,
        m: Matrix<FS::Set>,
    ) -> Result<Matrix<Polynomial<FS::Set>>, MatOppErr> {
        let n = m.rows();
        if n != m.cols() {
            Err(MatOppErr::NotSquare)
        } else {
            let poly_ring = PolynomialStructure::new(self.ring().clone());
            let poly_mat_struct = MatrixStructure::new(poly_ring.clone());
            Ok(poly_mat_struct
                .add(
                    &m.apply_map(|x| Polynomial::constant(x.clone())),
                    &poly_mat_struct
                        .neg(poly_mat_struct.diag(&(0..n).map(|_i| poly_ring.var()).collect())),
                )
                .unwrap())
        }
    }

    pub fn minimal_polynomial(&self, m: Matrix<FS::Set>) -> Result<Polynomial<FS::Set>, MatOppErr> {
        match self.presentation_matrix(m) {
            Ok(pres_mat) => {
                let poly_ring = PolynomialStructure::new(self.ring().clone());
                let poly_mat_struct = MatrixStructure::new(poly_ring.clone());
                let (_u, s, _v, k) = poly_mat_struct.smith_algorithm(pres_mat);
                debug_assert!(k > 0); //cant be all zero becasue we are taking SNF of a non-zero matrix
                Ok(s.at(k - 1, k - 1).unwrap().clone())
            }
            Err(MatOppErr::NotSquare) => Err(MatOppErr::NotSquare),
            Err(_) => panic!(),
        }
    }

    pub fn characteristic_polynomial(
        &self,
        m: Matrix<FS::Set>,
    ) -> Result<Polynomial<FS::Set>, MatOppErr> {
        match self.presentation_matrix(m) {
            Ok(pres_mat) => {
                let poly_ring = PolynomialStructure::new(self.ring().clone());
                let poly_mat_struct = MatrixStructure::new(poly_ring.clone());
                let (_u, s, _v, k) = poly_mat_struct.smith_algorithm(pres_mat);
                debug_assert!(k > 0); //cant be all zero becasue we are taking SNF of a non-zero matrix
                let mut char_poly = poly_ring.one();
                for i in 0..k {
                    poly_ring.mul_mut(&mut char_poly, s.at(i, i).unwrap())
                }
                Ok(char_poly)
            }
            Err(MatOppErr::NotSquare) => Err(MatOppErr::NotSquare),
            Err(_) => panic!(),
        }
    }
}

impl<F: MetaType> Matrix<F>
where
    F::Signature: FieldSignature,
{
    pub fn presentation_matrix(&self) -> Result<Matrix<Polynomial<F>>, MatOppErr> {
        Self::structure().presentation_matrix(self.clone())
    }

    pub fn minimal_polynomial(&self) -> Result<Polynomial<F>, MatOppErr> {
        Self::structure().minimal_polynomial(self.clone())
    }

    pub fn characteristic_polynomial(&self) -> Result<Polynomial<F>, MatOppErr> {
        Self::structure().characteristic_polynomial(self.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn min_and_char_polys() {
        {
            let a = Matrix::<Integer>::from_rows(vec![
                vec![Integer::from(0), Integer::from(4), Integer::from(4)],
                vec![Integer::from(1), Integer::from(4), Integer::from(16)],
            ]);
            let (_u, s, _v, _k) = a.smith_algorithm();

            assert_eq!(
                s,
                Matrix::from_rows(vec![
                    vec![Integer::from(1), Integer::from(0), Integer::from(0),],
                    vec![Integer::from(0), Integer::from(4), Integer::from(0)],
                ])
            );
        }

        {
            let a = Matrix::<Rational>::from_rows(vec![
                vec![Rational::from(0), Rational::from(0), Rational::from(0)],
                vec![Rational::from(0), Rational::from(0), Rational::from(1)],
                vec![Rational::from(0), Rational::from(0), Rational::from(0)],
            ]);
            let min_p = a.clone().minimal_polynomial().unwrap();
            let char_p = a.clone().characteristic_polynomial().unwrap();
            assert_eq!(
                &min_p,
                &Polynomial::from_coeffs(vec![
                    Rational::from(0),
                    Rational::from(0),
                    Rational::from(1)
                ])
            );
            assert_eq!(
                &char_p,
                &Polynomial::from_coeffs(vec![
                    Rational::from(0),
                    Rational::from(0),
                    Rational::from(0),
                    Rational::from(1)
                ])
            );
        }

        {
            let a = Matrix::<Rational>::from_rows(vec![
                vec![
                    Rational::from(4),
                    Rational::from(0),
                    Rational::from(0),
                    Rational::from(0),
                ],
                vec![
                    Rational::from(0),
                    Rational::from(4),
                    Rational::from(0),
                    Rational::from(0),
                ],
                vec![
                    Rational::from(0),
                    Rational::from(1),
                    Rational::from(4),
                    Rational::from(0),
                ],
                vec![
                    Rational::from(0),
                    Rational::from(0),
                    Rational::from(1),
                    Rational::from(4),
                ],
            ]);
            let min_p = a.clone().minimal_polynomial().unwrap();
            let char_p = a.clone().characteristic_polynomial().unwrap();
            assert_eq!(
                &min_p,
                &Polynomial::from_coeffs(vec![Rational::from(-4), Rational::from(1),])
                    .nat_pow(&Natural::from(3u8))
            );
            assert_eq!(
                &char_p,
                &Polynomial::from_coeffs(vec![Rational::from(-4), Rational::from(1),])
                    .nat_pow(&Natural::from(4u8))
            );
        }
    }
}
