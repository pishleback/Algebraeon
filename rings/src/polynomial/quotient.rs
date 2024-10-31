use super::polynomial::*;
use crate::{
    linear::matrix::*,
    ring_structure::{quotient::*, structure::*},
    structure::*,
};

impl<FS: FieldStructure, const IS_FIELD: bool> QuotientStructure<PolynomialStructure<FS>, IS_FIELD>
where
    PolynomialStructure<FS>: Structure<Set = Polynomial<FS::Set>> + UniqueFactorizationStructure,
{
    pub fn generator(&self) -> Polynomial<FS::Set> {
        self.ring().var()
    }

    pub fn degree(&self) -> usize {
        self.ring().degree(self.modulus()).unwrap()
    }

    //matrix representing column vector multiplication by a on the left
    pub fn col_multiplication_matrix(&self, a: &Polynomial<FS::Set>) -> Matrix<FS::Set> {
        let poly_ring = self.ring();
        let field = poly_ring.coeff_ring();

        let a_reduced = self.reduce(a);
        let deg = self.degree();
        Matrix::from_cols(
            (0..self.degree())
                .map(|i| {
                    let mut coeffs = self
                        .reduce(&poly_ring.mul(a, &poly_ring.var_pow(i)))
                        .into_coeffs();
                    debug_assert!(coeffs.len() <= deg);
                    while coeffs.len() < deg {
                        coeffs.push(field.zero());
                    }
                    coeffs
                })
                .collect(),
        )
    }
    //matrix representing row vector multiplication by a on the right
    pub fn row_multiplication_matrix(&self, a: &Polynomial<FS::Set>) -> Matrix<FS::Set> {
        self.col_multiplication_matrix(a).transpose()
    }

    pub fn to_col_vector(&self, a: &Polynomial<FS::Set>) -> Matrix<FS::Set> {
        let a_reduced = self.reduce(a);
        Matrix::construct(self.degree(), 1, |r, c| {
            self.ring().coeff(&a_reduced, r).clone()
        })
    }
    pub fn to_row_vector(&self, a: &Polynomial<FS::Set>) -> Matrix<FS::Set> {
        self.to_col_vector(a).transpose()
    }

    pub fn from_col_vector(&self, v: Matrix<FS::Set>) -> Polynomial<FS::Set> {
        assert_eq!(v.cols(), 1);
        assert_eq!(v.rows(), self.degree());
        Polynomial::from_coeffs(
            (0..self.degree())
                .map(|i| v.at(i, 0).unwrap().clone())
                .collect(),
        )
    }
    pub fn from_row_vector(&self, v: Matrix<FS::Set>) -> Polynomial<FS::Set> {
        self.from_col_vector(v.transpose())
    }

    pub fn min_poly(&self, a: &Polynomial<FS::Set>) -> Polynomial<FS::Set> {
        MatrixStructure::new(self.ring().coeff_ring())
            .minimal_polynomial(self.col_multiplication_matrix(a))
            .unwrap()
    }

    pub fn norm(&self, a: &Polynomial<FS::Set>) -> FS::Set {
        MatrixStructure::new(self.ring().coeff_ring())
            .det(self.col_multiplication_matrix(a))
            .unwrap()
    }

    pub fn trace(&self, a: &Polynomial<FS::Set>) -> FS::Set {
        MatrixStructure::new(self.ring().coeff_ring())
            .trace(&self.col_multiplication_matrix(a))
            .unwrap()
    }
}
