use std::borrow::Cow;

use super::{Polynomial, polynomial_ring::*};
use crate::{matrix::*, structure::*};
use algebraeon_nzq::Natural;
use algebraeon_sets::structure::*;

pub type FieldExtensionByPolynomialQuotientStructure<FS, FSB, FSPB> =
    QuotientStructure<PolynomialStructure<FS, FSB>, FSPB, true>;

impl<
    FS: FieldSignature + CharacteristicSignature,
    FSB: BorrowedStructure<FS>,
    FSPB: BorrowedStructure<PolynomialStructure<FS, FSB>>,
    const IS_FIELD: bool,
> CharacteristicSignature for QuotientStructure<PolynomialStructure<FS, FSB>, FSPB, IS_FIELD>
where
    PolynomialStructure<FS, FSB>: SetSignature<Set = Polynomial<FS::Set>>,
{
    fn characteristic(&self) -> Natural {
        self.ring().characteristic()
    }
}

impl<
    FS: FieldSignature,
    FSB: BorrowedStructure<FS>,
    FSPB: BorrowedStructure<PolynomialStructure<FS, FSB>>,
> FieldExtensionByPolynomialQuotientStructure<FS, FSB, FSPB>
where
    PolynomialStructure<FS, FSB>: SetSignature<Set = Polynomial<FS::Set>>,
{
    pub fn field_inclusion<'a>(
        &'a self,
    ) -> FieldExtensionByPolynomialQuotientInclusion<FS, FSB, FSPB> {
        FieldExtensionByPolynomialQuotientInclusion::new(self.clone())
    }
}

impl<
    FS: FieldSignature,
    FSB: BorrowedStructure<FS>,
    FSPB: BorrowedStructure<PolynomialStructure<FS, FSB>>,
    const IS_FIELD: bool,
> QuotientStructure<PolynomialStructure<FS, FSB>, FSPB, IS_FIELD>
where
    PolynomialStructure<FS, FSB>: SetSignature<Set = Polynomial<FS::Set>>,
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
    pub fn to_vector(&self, a: &Polynomial<FS::Set>) -> Vec<FS::Set> {
        let mut v = self.reduce(a).into_coeffs();
        debug_assert!(v.len() <= self.degree());
        while v.len() < self.degree() {
            v.push(self.ring().coeff_ring().zero());
        }
        v
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
    pub fn from_vector(&self, v: Vec<FS::Set>) -> Polynomial<FS::Set> {
        Polynomial::from_coeffs(v)
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
pub struct FieldExtensionByPolynomialQuotientInclusion<
    Field: FieldSignature,
    FieldB: BorrowedStructure<Field>,
    FieldPolyB: BorrowedStructure<PolynomialStructure<Field, FieldB>>,
> {
    extension_field: FieldExtensionByPolynomialQuotientStructure<Field, FieldB, FieldPolyB>,
}

impl<
    Field: FieldSignature,
    FieldB: BorrowedStructure<Field>,
    FieldPolyB: BorrowedStructure<PolynomialStructure<Field, FieldB>>,
> FieldExtensionByPolynomialQuotientInclusion<Field, FieldB, FieldPolyB>
{
    pub fn new(
        extension_field: FieldExtensionByPolynomialQuotientStructure<Field, FieldB, FieldPolyB>,
    ) -> Self {
        Self { extension_field }
    }
}

impl<
    Field: FieldSignature,
    FieldB: BorrowedStructure<Field>,
    FieldPolyB: BorrowedStructure<PolynomialStructure<Field, FieldB>>,
> Morphism<Field, FieldExtensionByPolynomialQuotientStructure<Field, FieldB, FieldPolyB>>
    for FieldExtensionByPolynomialQuotientInclusion<Field, FieldB, FieldPolyB>
{
    fn domain(&self) -> &Field {
        self.extension_field.ring().coeff_ring()
    }

    fn range(&self) -> &FieldExtensionByPolynomialQuotientStructure<Field, FieldB, FieldPolyB> {
        &self.extension_field
    }
}

impl<
    Field: FieldSignature,
    FieldB: BorrowedStructure<Field>,
    FieldPolyB: BorrowedStructure<PolynomialStructure<Field, FieldB>>,
> Function<Field, FieldExtensionByPolynomialQuotientStructure<Field, FieldB, FieldPolyB>>
    for FieldExtensionByPolynomialQuotientInclusion<Field, FieldB, FieldPolyB>
{
    fn image(&self, x: &Field::Set) -> Polynomial<Field::Set> {
        Polynomial::constant(x.clone())
    }
}

impl<
    Field: FieldSignature,
    FieldB: BorrowedStructure<Field>,
    FieldPolyB: BorrowedStructure<PolynomialStructure<Field, FieldB>>,
> RingHomomorphism<Field, FieldExtensionByPolynomialQuotientStructure<Field, FieldB, FieldPolyB>>
    for FieldExtensionByPolynomialQuotientInclusion<Field, FieldB, FieldPolyB>
{
}

impl<
    Field: FieldSignature,
    FieldB: BorrowedStructure<Field>,
    FieldPolyB: BorrowedStructure<PolynomialStructure<Field, FieldB>>,
> InjectiveFunction<Field, FieldExtensionByPolynomialQuotientStructure<Field, FieldB, FieldPolyB>>
    for FieldExtensionByPolynomialQuotientInclusion<Field, FieldB, FieldPolyB>
{
    fn try_preimage(&self, x: &Polynomial<Field::Set>) -> Option<Field::Set> {
        PolynomialStructure::new(self.domain().clone()).as_constant(&self.range().reduce(x))
    }
}

impl<
    Field: FieldSignature,
    FieldB: BorrowedStructure<Field>,
    FieldPolyB: BorrowedStructure<PolynomialStructure<Field, FieldB>>,
>
    FiniteRankFreeRingExtension<
        Field,
        FieldExtensionByPolynomialQuotientStructure<Field, FieldB, FieldPolyB>,
    > for FieldExtensionByPolynomialQuotientInclusion<Field, FieldB, FieldPolyB>
{
    type Basis = EnumeratedFiniteSetStructure;

    fn basis_set(&self) -> impl std::borrow::Borrow<Self::Basis> {
        EnumeratedFiniteSetStructure::new(self.range().degree())
    }

    fn to_component<'a>(&self, b: &usize, v: &'a Polynomial<Field::Set>) -> Cow<'a, Field::Set> {
        Cow::Owned(
            self.domain()
                .polynomials()
                .coeff(&self.range().reduce(v), *b)
                .clone(),
        )
    }

    fn from_component(&self, b: &usize, r: &Field::Set) -> Polynomial<Field::Set> {
        self.domain().polynomials().constant_var_pow(r.clone(), *b)
    }
}

impl<
    Field: FieldSignature,
    FieldB: BorrowedStructure<Field>,
    FieldPolyB: BorrowedStructure<PolynomialStructure<Field, FieldB>>,
>
    FiniteDimensionalFieldExtension<
        Field,
        FieldExtensionByPolynomialQuotientStructure<Field, FieldB, FieldPolyB>,
    > for FieldExtensionByPolynomialQuotientInclusion<Field, FieldB, FieldPolyB>
{
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
            let f = p.algebraic_number_field().unwrap();
            let ext = FieldExtensionByPolynomialQuotientInclusion::new(f);
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
            let f = p.algebraic_number_field().unwrap();
            let ext = FieldExtensionByPolynomialQuotientInclusion::new(f);
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
