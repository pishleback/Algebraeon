use super::{Polynomial, polynomial_structure::*};
use crate::{matrix::*, structure::*};
use algebraeon_sets::sets::EnumeratedFiniteSetStructure;
use algebraeon_structures::*;
use itertools::Itertools;
use std::{
    borrow::{Borrow, Cow},
    rc::Rc,
};

pub type PolynomialQuotientRingStructure<FS, const IS_FIELD: bool> =
    EuclideanRemainderQuotientStructure<PolynomialStructure<FS>, IS_FIELD>;

impl<FS: FieldSignature + CharacteristicSignature, const IS_FIELD: bool> CharacteristicSignature
    for PolynomialQuotientRingStructure<FS, IS_FIELD>
where
    PolynomialStructure<FS>: SetSignature<Elem = Polynomial<FS::Elem>>,
{
    fn characteristic(self: &Rc<Self>) -> Natural {
        self.ring().characteristic()
    }
}

impl<FS: FieldSignature + CharZeroRingSignature, const IS_FIELD: bool> CharZeroRingSignature
    for PolynomialQuotientRingStructure<FS, IS_FIELD>
where
    PolynomialStructure<FS>: SetSignature<Elem = Polynomial<FS::Elem>>,
{
    fn try_to_int(self: &Rc<Self>, x: &Self::Elem) -> Option<Integer> {
        let x_reduced = self.reduce(x);
        self.ring().try_to_int(&x_reduced)
    }
}

impl<FS: FieldSignature, const IS_FIELD: bool> PolynomialQuotientRingStructure<FS, IS_FIELD>
where
    PolynomialStructure<FS>: SetSignature<Elem = Polynomial<FS::Elem>>,
{
    pub fn coefficient_ring_inclusion(
        self: &Rc<Self>,
    ) -> Rc<PolynomialQuotientRingExtension<FS, IS_FIELD>> {
        PolynomialQuotientRingExtension::new(self.clone())
    }
}

impl<FS: FieldSignature, const IS_FIELD: bool> PolynomialQuotientRingStructure<FS, IS_FIELD>
where
    PolynomialStructure<FS>: SetSignature<Elem = Polynomial<FS::Elem>>,
{
    pub fn generator(self: &Rc<Self>) -> Polynomial<FS::Elem> {
        self.ring().var()
    }

    pub fn col_multiplication_matrix(
        self: &Rc<Self>,
        a: &Polynomial<FS::Elem>,
    ) -> Matrix<FS::Elem> {
        self.coefficient_ring_inclusion()
            .col_multiplication_matrix(a)
    }

    pub fn row_multiplication_matrix(
        self: &Rc<Self>,
        a: &Polynomial<FS::Elem>,
    ) -> Matrix<FS::Elem> {
        self.coefficient_ring_inclusion()
            .row_multiplication_matrix(a)
    }

    pub fn to_col(self: &Rc<Self>, a: &Polynomial<FS::Elem>) -> Matrix<FS::Elem> {
        self.coefficient_ring_inclusion()
            .range_module_structure()
            .to_col(a)
    }

    pub fn to_row(self: &Rc<Self>, a: &Polynomial<FS::Elem>) -> Matrix<FS::Elem> {
        self.coefficient_ring_inclusion()
            .range_module_structure()
            .to_row(a)
    }

    pub fn to_vec(self: &Rc<Self>, a: &Polynomial<FS::Elem>) -> Vec<FS::Elem> {
        self.coefficient_ring_inclusion()
            .range_module_structure()
            .to_vec(a)
    }

    pub fn from_col(self: &Rc<Self>, v: Matrix<FS::Elem>) -> Polynomial<FS::Elem> {
        self.coefficient_ring_inclusion()
            .range_module_structure()
            .from_col(v)
    }

    pub fn from_row(self: &Rc<Self>, v: Matrix<FS::Elem>) -> Polynomial<FS::Elem> {
        self.coefficient_ring_inclusion()
            .range_module_structure()
            .from_row(v)
    }

    pub fn from_vec(self: &Rc<Self>, v: Vec<impl Borrow<FS::Elem>>) -> Polynomial<FS::Elem> {
        self.coefficient_ring_inclusion()
            .range_module_structure()
            .from_vec(v)
    }

    pub fn degree(&self) -> usize {
        self.ring().degree(self.modulus().as_ref()).unwrap()
    }
}

impl<FS: FieldSignature> PolynomialQuotientRingStructure<FS, true>
where
    PolynomialStructure<FS>: SetSignature<Elem = Polynomial<FS::Elem>>,
{
    pub fn min_poly(self: &Rc<Self>, a: &Polynomial<FS::Elem>) -> Polynomial<FS::Elem> {
        self.coefficient_ring_inclusion().min_poly(a)
    }

    pub fn norm(self: &Rc<Self>, a: &Polynomial<FS::Elem>) -> FS::Elem {
        self.coefficient_ring_inclusion().norm(a)
    }

    pub fn trace(self: &Rc<Self>, a: &Polynomial<FS::Elem>) -> FS::Elem {
        self.coefficient_ring_inclusion().trace(a)
    }
}

impl<const IS_FIELD: bool, FS: FieldSignature + FiniteSetSignature> CountableSetSignature
    for PolynomialQuotientRingStructure<FS, IS_FIELD>
where
    PolynomialStructure<FS>: SetSignature<Elem = Polynomial<FS::Elem>>,
{
    fn generate_all_elements(self: Rc<Self>) -> impl Iterator<Item = Self::Elem> {
        let n = self.coefficient_ring_inclusion().degree();
        let ring_elements = self.ring().coeff_ring().list_all_elements();
        let hom = self.coefficient_ring_inclusion().clone();
        (0..n)
            .map(move |_| ring_elements.clone())
            .multi_cartesian_product()
            .map(move |v| hom.from_vec(v))
    }
}

impl<const IS_FIELD: bool, FS: FieldSignature + FiniteSetSignature> FiniteSetSignature
    for PolynomialQuotientRingStructure<FS, IS_FIELD>
where
    PolynomialStructure<FS>: SetSignature<Elem = Polynomial<FS::Elem>>,
{
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PolynomialQuotientRingExtension<Field: FieldSignature, const IS_FIELD: bool> {
    polynomial_quotient_ring: Rc<PolynomialQuotientRingStructure<Field, IS_FIELD>>,
}

impl<Field: FieldSignature, const IS_FIELD: bool> PolynomialQuotientRingExtension<Field, IS_FIELD> {
    pub fn new(
        polynomial_quotient_ring: Rc<PolynomialQuotientRingStructure<Field, IS_FIELD>>,
    ) -> Rc<Self> {
        Self {
            polynomial_quotient_ring,
        }
        .into()
    }
}

impl<Field: FieldSignature, const IS_FIELD: bool>
    Morphism<Field, PolynomialQuotientRingStructure<Field, IS_FIELD>>
    for PolynomialQuotientRingExtension<Field, IS_FIELD>
{
    fn domain(self: &Rc<Self>) -> Rc<Field> {
        self.polynomial_quotient_ring.ring().coeff_ring()
    }

    fn range(self: &Rc<Self>) -> Rc<PolynomialQuotientRingStructure<Field, IS_FIELD>> {
        self.polynomial_quotient_ring.clone()
    }
}

impl<Field: FieldSignature, const IS_FIELD: bool>
    FunctionMorphism<Field, PolynomialQuotientRingStructure<Field, IS_FIELD>>
    for PolynomialQuotientRingExtension<Field, IS_FIELD>
{
    fn image(self: &Rc<Self>, x: &Field::Elem) -> Polynomial<Field::Elem> {
        Polynomial::constant(x.clone())
    }
}

impl<Field: FieldSignature, const IS_FIELD: bool>
    RingHomomorphism<Field, PolynomialQuotientRingStructure<Field, IS_FIELD>>
    for PolynomialQuotientRingExtension<Field, IS_FIELD>
{
}

impl<Field: FieldSignature, const IS_FIELD: bool>
    InjectiveFunctionMorphism<Field, PolynomialQuotientRingStructure<Field, IS_FIELD>>
    for PolynomialQuotientRingExtension<Field, IS_FIELD>
{
    fn try_preimage(self: &Rc<Self>, x: &Polynomial<Field::Elem>) -> Option<Field::Elem> {
        self.domain()
            .polynomials()
            .as_constant(&self.range().reduce(x))
    }
}

impl<Field: FieldSignature, const IS_FIELD: bool>
    FreeModuleSignature<EnumeratedFiniteSetStructure, Field>
    for RingHomomorphismRangeModuleStructure<
        Field,
        PolynomialQuotientRingStructure<Field, IS_FIELD>,
        PolynomialQuotientRingExtension<Field, IS_FIELD>,
    >
{
    fn basis_set(self: &Rc<Self>) -> Rc<EnumeratedFiniteSetStructure> {
        EnumeratedFiniteSetStructure::new(self.module().degree())
    }

    fn to_component<'a>(
        self: &Rc<Self>,
        b: &usize,
        v: &'a Polynomial<Field::Elem>,
    ) -> Cow<'a, Field::Elem> {
        Cow::Owned(
            self.ring()
                .polynomials()
                .coeff(&self.module().reduce(v), *b)
                .into_owned(),
        )
    }

    fn from_component(self: &Rc<Self>, b: &usize, r: &Field::Elem) -> Polynomial<Field::Elem> {
        self.ring().polynomials().constant_var_pow(r.clone(), *b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

    #[test]
    fn finite_dimensional_field_extension_structure() {
        let x = Rational::structure().polynomials().var().into_ergonomic();
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
