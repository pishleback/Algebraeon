use std::borrow::Cow;

use super::{Polynomial, polynomial_ring::*};
use crate::{matrix::*, structure::*};
use algebraeon_nzq::{Integer, Natural};
use algebraeon_sets::structure::*;

pub type PolynomialQuotientRingStructure<FS, FSB, FSPB, const IS_FIELD: bool> =
    QuotientStructure<PolynomialStructure<FS, FSB>, FSPB, IS_FIELD>;

impl<
    FS: FieldSignature + CharacteristicSignature,
    FSB: BorrowedStructure<FS>,
    FSPB: BorrowedStructure<PolynomialStructure<FS, FSB>>,
    const IS_FIELD: bool,
> CharacteristicSignature for PolynomialQuotientRingStructure<FS, FSB, FSPB, IS_FIELD>
where
    PolynomialStructure<FS, FSB>: SetSignature<Set = Polynomial<FS::Set>>,
{
    fn characteristic(&self) -> Natural {
        self.ring().characteristic()
    }
}

impl<
    FS: FieldSignature + CharZeroRingSignature,
    FSB: BorrowedStructure<FS>,
    FSPB: BorrowedStructure<PolynomialStructure<FS, FSB>>,
    const IS_FIELD: bool,
> CharZeroRingSignature for PolynomialQuotientRingStructure<FS, FSB, FSPB, IS_FIELD>
where
    PolynomialStructure<FS, FSB>: SetSignature<Set = Polynomial<FS::Set>>,
{
    fn try_to_int(&self, x: &Self::Set) -> Option<Integer> {
        let x_reduced = self.reduce(x);
        self.ring().try_to_int(&x_reduced)
    }
}

impl<
    FS: FieldSignature,
    FSB: BorrowedStructure<FS>,
    FSPB: BorrowedStructure<PolynomialStructure<FS, FSB>>,
    const IS_FIELD: bool,
> PolynomialQuotientRingStructure<FS, FSB, FSPB, IS_FIELD>
where
    PolynomialStructure<FS, FSB>: SetSignature<Set = Polynomial<FS::Set>>,
{
    pub fn coefficient_ring_inclusion<'a>(
        &'a self,
    ) -> PolynomialQuotientRingExtension<FS, FSB, FSPB, IS_FIELD> {
        PolynomialQuotientRingExtension::new(self.clone())
    }
}

impl<
    FS: FieldSignature,
    FSB: BorrowedStructure<FS>,
    FSPB: BorrowedStructure<PolynomialStructure<FS, FSB>>,
    const IS_FIELD: bool,
> PolynomialQuotientRingStructure<FS, FSB, FSPB, IS_FIELD>
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
            self.ring().coeff(&a_reduced, r).into_owned()
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
pub struct PolynomialQuotientRingExtension<
    Field: FieldSignature,
    FieldB: BorrowedStructure<Field>,
    FieldPolyB: BorrowedStructure<PolynomialStructure<Field, FieldB>>,
    const IS_FIELD: bool,
> {
    polynomial_quotient_ring: PolynomialQuotientRingStructure<Field, FieldB, FieldPolyB, IS_FIELD>,
}

impl<
    Field: FieldSignature,
    FieldB: BorrowedStructure<Field>,
    FieldPolyB: BorrowedStructure<PolynomialStructure<Field, FieldB>>,
    const IS_FIELD: bool,
> PolynomialQuotientRingExtension<Field, FieldB, FieldPolyB, IS_FIELD>
{
    pub fn new(
        polynomial_quotient_ring: PolynomialQuotientRingStructure<
            Field,
            FieldB,
            FieldPolyB,
            IS_FIELD,
        >,
    ) -> Self {
        Self {
            polynomial_quotient_ring,
        }
    }
}

impl<
    Field: FieldSignature,
    FieldB: BorrowedStructure<Field>,
    FieldPolyB: BorrowedStructure<PolynomialStructure<Field, FieldB>>,
    const IS_FIELD: bool,
> Morphism<Field, PolynomialQuotientRingStructure<Field, FieldB, FieldPolyB, IS_FIELD>>
    for PolynomialQuotientRingExtension<Field, FieldB, FieldPolyB, IS_FIELD>
{
    fn domain(&self) -> &Field {
        self.polynomial_quotient_ring.ring().coeff_ring()
    }

    fn range(&self) -> &PolynomialQuotientRingStructure<Field, FieldB, FieldPolyB, IS_FIELD> {
        &self.polynomial_quotient_ring
    }
}

impl<
    Field: FieldSignature,
    FieldB: BorrowedStructure<Field>,
    FieldPolyB: BorrowedStructure<PolynomialStructure<Field, FieldB>>,
    const IS_FIELD: bool,
> Function<Field, PolynomialQuotientRingStructure<Field, FieldB, FieldPolyB, IS_FIELD>>
    for PolynomialQuotientRingExtension<Field, FieldB, FieldPolyB, IS_FIELD>
{
    fn image(&self, x: &Field::Set) -> Polynomial<Field::Set> {
        Polynomial::constant(x.clone())
    }
}

impl<
    Field: FieldSignature,
    FieldB: BorrowedStructure<Field>,
    FieldPolyB: BorrowedStructure<PolynomialStructure<Field, FieldB>>,
    const IS_FIELD: bool,
> RingHomomorphism<Field, PolynomialQuotientRingStructure<Field, FieldB, FieldPolyB, IS_FIELD>>
    for PolynomialQuotientRingExtension<Field, FieldB, FieldPolyB, IS_FIELD>
{
}

impl<
    Field: FieldSignature,
    FieldB: BorrowedStructure<Field>,
    FieldPolyB: BorrowedStructure<PolynomialStructure<Field, FieldB>>,
    const IS_FIELD: bool,
> InjectiveFunction<Field, PolynomialQuotientRingStructure<Field, FieldB, FieldPolyB, IS_FIELD>>
    for PolynomialQuotientRingExtension<Field, FieldB, FieldPolyB, IS_FIELD>
{
    fn try_preimage(&self, x: &Polynomial<Field::Set>) -> Option<Field::Set> {
        self.domain()
            .polynomial_ring()
            .as_constant(&self.range().reduce(x))
    }
}

impl<
    'h,
    Field: FieldSignature,
    FieldB: BorrowedStructure<Field>,
    FieldPolyB: BorrowedStructure<PolynomialStructure<Field, FieldB>>,
    const IS_FIELD: bool,
> FreeModuleSignature<Field>
    for RingHomomorphismRangeModuleStructure<
        'h,
        Field,
        PolynomialQuotientRingStructure<Field, FieldB, FieldPolyB, IS_FIELD>,
        PolynomialQuotientRingExtension<Field, FieldB, FieldPolyB, IS_FIELD>,
    >
{
    type Basis = EnumeratedFiniteSetStructure;

    fn basis_set(&self) -> impl std::borrow::Borrow<Self::Basis> {
        EnumeratedFiniteSetStructure::new(self.module().degree())
    }

    fn to_component<'a>(&self, b: &usize, v: &'a Polynomial<Field::Set>) -> Cow<'a, Field::Set> {
        Cow::Owned(
            self.ring()
                .polynomial_ring()
                .coeff(&self.module().reduce(v), *b)
                .into_owned(),
        )
    }

    fn from_component(&self, b: &usize, r: &Field::Set) -> Polynomial<Field::Set> {
        self.ring()
            .polynomial_ring()
            .constant_var_pow(r.clone(), *b)
    }
}

impl<
    'h,
    Field: FieldSignature,
    FieldB: BorrowedStructure<Field>,
    FieldPolyB: BorrowedStructure<PolynomialStructure<Field, FieldB>>,
    const IS_FIELD: bool,
> FinitelyFreeModuleSignature<Field>
    for RingHomomorphismRangeModuleStructure<
        'h,
        Field,
        PolynomialQuotientRingStructure<Field, FieldB, FieldPolyB, IS_FIELD>,
        PolynomialQuotientRingExtension<Field, FieldB, FieldPolyB, IS_FIELD>,
    >
{
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_nzq::Rational;
    use std::str::FromStr;

    #[test]
    fn finite_dimensional_field_extension_structure() {
        let x = Rational::structure()
            .into_polynomial_ring()
            .var()
            .into_ergonomic();
        {
            let p = (x.pow(3) + &x - 1).into_verbose();
            let f = p.algebraic_number_field().unwrap();
            let ext = PolynomialQuotientRingExtension::new(f);
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
            let ext = PolynomialQuotientRingExtension::new(f);
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
