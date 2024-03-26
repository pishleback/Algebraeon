use malachite_q::Rational;

use crate::rings::{
    linear::matrix::Matrix,
    polynomial::polynomial::*,
    ring_structure::{cannonical::Ring, quotient::*},
    structure::*,
};

type ANFStructure = QuotientStructure<PolynomialStructure<CannonicalStructure<Rational>>, true>;

fn new_anf(f: Polynomial<Rational>) -> ANFStructure {
    ANFStructure::new(PolynomialStructure::new(Rational::structure()).into(), f)
}

impl ANFStructure {
    pub fn degree(&self) -> usize {
        self.modulus().degree().unwrap()
    }

    //matrix representing column vector multiplication by a on the left
    pub fn col_multiplication_matrix(&self, a: &Polynomial<Rational>) -> Matrix<Rational> {
        let a_reduced = self.reduce(a);
        let deg = self.degree();
        Matrix::from_cols(
            (0..self.degree())
                .map(|i| {
                    let mut coeffs = self
                        .reduce(&Polynomial::mul(a, &Polynomial::var_pow(i)))
                        .into_coeffs();
                    debug_assert!(coeffs.len() <= deg);
                    while coeffs.len() < deg {
                        coeffs.push(Rational::zero());
                    }
                    coeffs
                })
                .collect(),
        )
    }
    //matrix representing row vector multiplication by a on the right
    pub fn row_multiplication_matrix(&self, a: &Polynomial<Rational>) -> Matrix<Rational> {
        self.col_multiplication_matrix(a).transpose()
    }

    pub fn to_col_vector(&self, a: &Polynomial<Rational>) -> Matrix<Rational> {
        let a_reduced = self.reduce(a);
        Matrix::construct(self.degree(), 1, |r, c| a_reduced.coeff(r))
    }
    pub fn to_row_vector(&self, a: &Polynomial<Rational>) -> Matrix<Rational> {
        self.to_col_vector(a).transpose()
    }

    pub fn from_col_vector(&self, v: Matrix<Rational>) -> Polynomial<Rational> {
        assert_eq!(v.cols(), 1);
        assert_eq!(v.rows(), self.degree());
        Polynomial::from_coeffs(
            (0..self.degree())
                .map(|i| v.at(i, 0).unwrap().clone())
                .collect(),
        )
    }
    pub fn from_row_vector(&self, v: Matrix<Rational>) -> Polynomial<Rational> {
        self.from_col_vector(v.transpose())
    }

    pub fn min_poly(&self, a: &Polynomial<Rational>) -> Polynomial<Rational> {
        self.col_multiplication_matrix(a).minimal_polynomial().unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_anf_to_and_from_vector() {
        let x = &Polynomial::<Rational>::var().into_ring();
        let anf = new_anf((x.pow(5) - x + 1).into_set());
        let alpha = (x.pow(9) + 5).into_set();

        println!("{}", alpha);
        println!("{}", anf.reduce(&alpha));
        anf.to_col_vector(&alpha).pprint();

        assert_eq!(
            anf.to_col_vector(&alpha),
            Matrix::from_cols(vec![vec![
                Rational::from(4),
                Rational::from(1),
                Rational::from(0),
                Rational::from(0),
                Rational::from(-1)
            ]])
        );

        assert!(anf.equal(
            &anf.from_col_vector(Matrix::from_cols(vec![vec![
                Rational::from(4),
                Rational::from(1),
                Rational::from(0),
                Rational::from(0),
                Rational::from(-1)
            ]])),
            &alpha
        ));
    }

    #[test]
    fn test_anf() {
        let x = &Polynomial::<Rational>::var().into_ring();
        let anf = new_anf((x.pow(5) - x + 1).into_set());
        println!("{:?}", anf);

        let alpha = (x.pow(1)).into_set();
        println!("{}", anf.min_poly(&alpha));
    }
}
