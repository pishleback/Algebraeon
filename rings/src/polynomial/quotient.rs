use super::polynomial::*;
use crate::{linear::matrix::*, rings::quotient::QuotientStructure, structure::*};
use algebraeon_sets::structure::*;

pub type FieldExtensionByPolynomialQuotientStructure<FS> =
    QuotientStructure<PolynomialStructure<FS>, true>;

impl<FS: FieldSignature, const IS_FIELD: bool> QuotientStructure<PolynomialStructure<FS>, IS_FIELD>
where
    PolynomialStructure<FS>: SetSignature<Set = Polynomial<FS::Set>>,
{
    pub fn generator(&self) -> Polynomial<FS::Set> {
        self.ring().var()
    }

    //matrix representing column vector multiplication by a on the left
    pub fn col_multiplication_matrix(&self, a: &Polynomial<FS::Set>) -> Matrix<FS::Set> {
        let poly_ring = self.ring();
        let field = poly_ring.coeff_ring();
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
        Matrix::construct(self.degree(), 1, |r, _c| {
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
        MatrixStructure::new(self.ring().coeff_ring().clone())
            .minimal_polynomial(self.col_multiplication_matrix(a))
            .unwrap()
    }

    pub fn degree(&self) -> usize {
        self.ring().degree(self.modulus()).unwrap()
    }

    pub fn norm(&self, a: &Polynomial<FS::Set>) -> FS::Set {
        MatrixStructure::new(self.ring().coeff_ring().clone())
            .det(self.col_multiplication_matrix(a))
            .unwrap()
    }

    pub fn trace(&self, a: &Polynomial<FS::Set>) -> FS::Set {
        MatrixStructure::new(self.ring().coeff_ring().clone())
            .trace(&self.col_multiplication_matrix(a))
            .unwrap()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FieldExtensionByPolynomialQuotient<Field: FieldSignature> {
    extension_field: FieldExtensionByPolynomialQuotientStructure<Field>,
}

impl<Field: FieldSignature> FieldExtensionByPolynomialQuotient<Field> {
    pub fn new(extension_field: FieldExtensionByPolynomialQuotientStructure<Field>) -> Self {
        Self { extension_field }
    }
}

impl<Field: FieldSignature> Morphism<Field, FieldExtensionByPolynomialQuotientStructure<Field>>
    for FieldExtensionByPolynomialQuotient<Field>
{
    fn domain(&self) -> &Field {
        self.extension_field.ring().coeff_ring()
    }

    fn range(&self) -> &FieldExtensionByPolynomialQuotientStructure<Field> {
        &self.extension_field
    }
}

impl<Field: FieldSignature> Function<Field, FieldExtensionByPolynomialQuotientStructure<Field>>
    for FieldExtensionByPolynomialQuotient<Field>
{
    fn image(&self, x: &Field::Set) -> Polynomial<Field::Set> {
        Polynomial::constant(x.clone())
    }
}

impl<Field: FieldSignature>
    RingHomomorphism<Field, FieldExtensionByPolynomialQuotientStructure<Field>>
    for FieldExtensionByPolynomialQuotient<Field>
{
}

impl<Field: FieldSignature>
    InjectiveFunction<Field, FieldExtensionByPolynomialQuotientStructure<Field>>
    for FieldExtensionByPolynomialQuotient<Field>
{
    fn try_preimage(&self, x: &Polynomial<Field::Set>) -> Option<Field::Set> {
        PolynomialStructure::new(self.domain().clone()).as_constant(&self.range().reduce(x))
    }
}

impl<Field: FieldSignature>
    FiniteDimensionalFieldExtension<Field, FieldExtensionByPolynomialQuotientStructure<Field>>
    for FieldExtensionByPolynomialQuotient<Field>
{
    fn degree(&self) -> usize {
        self.range().degree()
    }

    fn norm(&self, a: &Polynomial<Field::Set>) -> Field::Set {
        self.range().norm(a)
    }

    fn trace(&self, a: &Polynomial<Field::Set>) -> Field::Set {
        self.range().trace(a)
    }

    fn min_poly(&self, a: &Polynomial<Field::Set>) -> Polynomial<<Field>::Set> {
        self.range().min_poly(a)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_nzq::Rational;
    use std::str::FromStr;

    #[test]
    fn finite_dimensional_field_extension_structure() {
        let x = PolynomialStructure::new(Rational::structure())
            .var()
            .into_ergonomic();
        {
            let p = (x.pow(3) + &x - 1).into_verbose();
            let f = p.algebraic_number_field();
            let ext = FieldExtensionByPolynomialQuotient::new(f);
            assert_eq!(ext.degree(), 3);
            assert_eq!(
                ext.image(&Rational::from_str("4").unwrap()),
                (4 * x.pow(0)).into_verbose()
            );
            assert_eq!(ext.try_preimage(&(3 * x.pow(2) + 1).into_verbose()), None);
            assert_eq!(
                ext.try_preimage(&(x.pow(3) + &x + 1).into_verbose()),
                Some(Rational::from_str("2").unwrap())
            );

            assert_eq!(
                ext.norm(&(5 * x.pow(1) + 2).into_verbose()),
                Rational::from_str("183").unwrap()
            );
        }
        {
            // Z[i]
            let p = (x.pow(2) + 1).into_verbose();
            let f = p.algebraic_number_field();
            let ext = FieldExtensionByPolynomialQuotient::new(f);
            assert_eq!(ext.degree(), 2);
            // a^2 + b^2
            assert_eq!(
                ext.norm(&(3 + 4 * &x).into_verbose()),
                Rational::from_str("25").unwrap()
            );
            // 2a
            assert_eq!(
                ext.trace(&(3 + 4 * &x).into_verbose()),
                Rational::from_str("6").unwrap()
            );
            // min_poly(1+i) = x^2 - 2x + 2
            assert_eq!(
                ext.min_poly(&(1 + &x).into_verbose()),
                (x.pow(2) - 2 * &x + 2).into_verbose()
            );
        }
    }
}
