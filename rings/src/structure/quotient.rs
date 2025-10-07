use crate::structure::*;
use algebraeon_sets::structure::*;
use std::{borrow::Borrow, marker::PhantomData};

#[derive(Debug, Clone)]
pub struct QuotientStructure<
    RS: EuclideanDomainSignature,
    RSB: BorrowedStructure<RS>,
    const IS_FIELD: bool,
> {
    _ring: PhantomData<RS>,
    ring: RSB,
    modulus: RS::Set,
}

impl<RS: EuclideanDomainSignature, RSB: BorrowedStructure<RS>, const IS_FIELD: bool>
    QuotientStructure<RS, RSB, IS_FIELD>
{
    pub fn new_unchecked(ring: RSB, modulus: RS::Set) -> Self {
        assert!(!ring.borrow().is_zero(&modulus));
        Self {
            _ring: PhantomData,
            ring,
            modulus,
        }
    }

    pub fn ring(&self) -> &RS {
        self.ring.borrow()
    }

    pub fn modulus(&self) -> &RS::Set {
        &self.modulus
    }

    pub fn reduce(&self, a: impl Borrow<RS::Set>) -> RS::Set {
        self.ring().rem(a.borrow(), &self.modulus)
    }
}

impl<RS: EuclideanDomainSignature, RSB: BorrowedStructure<RS>> QuotientStructure<RS, RSB, false> {
    fn new_ring(ring: RSB, modulus: RS::Set) -> Self {
        Self::new_unchecked(ring, modulus)
    }
}

impl<RS: EuclideanDomainSignature + FactorableSignature, RSB: BorrowedStructure<RS>>
    QuotientStructure<RS, RSB, true>
{
    fn new_field_unchecked(ring: RSB, modulus: RS::Set) -> Self {
        debug_assert!(ring.borrow().is_irreducible(&modulus));
        Self::new_unchecked(ring, modulus)
    }

    fn new_field(ring: RSB, modulus: RS::Set) -> Result<Self, ()> {
        if ring.borrow().is_irreducible(&modulus) {
            Ok(Self::new_unchecked(ring, modulus))
        } else {
            Err(())
        }
    }
}

pub trait RingToQuotientRingSignature: EuclideanDomainSignature {
    fn quotient_ring<'a>(&'a self, modulus: Self::Set) -> QuotientStructure<Self, &'a Self, false> {
        QuotientStructure::new_ring(self, modulus)
    }

    fn into_quotient_ring(self, modulus: Self::Set) -> QuotientStructure<Self, Self, false> {
        QuotientStructure::new_ring(self, modulus)
    }
}
impl<Ring: EuclideanDomainSignature> RingToQuotientRingSignature for Ring {}

pub trait RingToQuotientFieldSignature: EuclideanDomainSignature + FactorableSignature {
    fn quotient_field<'a>(
        &'a self,
        modulus: Self::Set,
    ) -> Result<QuotientStructure<Self, &'a Self, true>, ()> {
        QuotientStructure::new_field(self, modulus)
    }

    fn into_quotient_field(
        self,
        modulus: Self::Set,
    ) -> Result<QuotientStructure<Self, Self, true>, ()> {
        QuotientStructure::new_field(self, modulus)
    }

    fn quotient_field_unchecked<'a>(
        &'a self,
        modulus: Self::Set,
    ) -> QuotientStructure<Self, &'a Self, true> {
        QuotientStructure::new_field_unchecked(self, modulus)
    }

    fn into_quotient_field_unchecked(
        self,
        modulus: Self::Set,
    ) -> QuotientStructure<Self, Self, true> {
        QuotientStructure::new_field_unchecked(self, modulus)
    }
}
impl<Ring: EuclideanDomainSignature + FactorableSignature> RingToQuotientFieldSignature for Ring {}

impl<RS: EuclideanDomainSignature, RSB: BorrowedStructure<RS>, const IS_FIELD: bool> PartialEq
    for QuotientStructure<RS, RSB, IS_FIELD>
{
    fn eq(&self, other: &Self) -> bool {
        if self.ring == other.ring {
            debug_assert_eq!(
                self.ring().equal(&self.modulus, &other.modulus),
                other.ring().equal(&self.modulus, &other.modulus)
            );
            self.ring().equal(&self.modulus, &other.modulus)
        } else {
            false
        }
    }
}

impl<RS: EuclideanDomainSignature, RSB: BorrowedStructure<RS>, const IS_FIELD: bool> Eq
    for QuotientStructure<RS, RSB, IS_FIELD>
{
}

impl<RS: EuclideanDomainSignature, RSB: BorrowedStructure<RS>, const IS_FIELD: bool> Signature
    for QuotientStructure<RS, RSB, IS_FIELD>
{
}

impl<RS: EuclideanDomainSignature, RSB: BorrowedStructure<RS>, const IS_FIELD: bool> SetSignature
    for QuotientStructure<RS, RSB, IS_FIELD>
{
    type Set = RS::Set;

    fn is_element(&self, _x: &Self::Set) -> Result<(), String> {
        Ok(())
    }
}

impl<
    RS: EuclideanDomainSignature + ToStringSignature,
    RSB: BorrowedStructure<RS>,
    const IS_FIELD: bool,
> ToStringSignature for QuotientStructure<RS, RSB, IS_FIELD>
{
    fn to_string(&self, elem: &Self::Set) -> String {
        self.ring().to_string(elem)
    }
}

impl<RS: EuclideanDomainSignature, RSB: BorrowedStructure<RS>, const IS_FIELD: bool> EqSignature
    for QuotientStructure<RS, RSB, IS_FIELD>
{
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        self.ring().is_zero(
            &self
                .ring()
                .rem(&self.ring().add(a, &self.ring().neg(b)), &self.modulus),
        )
    }
}

impl<RS: EuclideanDomainSignature, RSB: BorrowedStructure<RS>, const IS_FIELD: bool>
    AdditiveMonoidSignature for QuotientStructure<RS, RSB, IS_FIELD>
{
    fn zero(&self) -> Self::Set {
        self.ring().zero()
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.ring().rem(&self.ring().add(a, b), &self.modulus)
    }
}

impl<RS: EuclideanDomainSignature, RSB: BorrowedStructure<RS>, const IS_FIELD: bool>
    AdditiveGroupSignature for QuotientStructure<RS, RSB, IS_FIELD>
{
    fn neg(&self, a: &Self::Set) -> Self::Set {
        self.ring().neg(a)
    }

    fn sub(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.ring().sub(a, b)
    }
}

impl<RS: EuclideanDomainSignature, RSB: BorrowedStructure<RS>, const IS_FIELD: bool>
    SemiRingSignature for QuotientStructure<RS, RSB, IS_FIELD>
{
    fn one(&self) -> Self::Set {
        self.ring().one()
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.ring().rem(&self.ring().mul(a, b), &self.modulus)
    }
}

impl<RS: EuclideanDomainSignature, RSB: BorrowedStructure<RS>, const IS_FIELD: bool>
    RingSignature for QuotientStructure<RS, RSB, IS_FIELD>
{
    fn is_reduced(&self) -> Result<bool, String> {
        if IS_FIELD {
            return Ok(true);
        }
        Err("not implemented yet: is_reduced for quotient rings that are not fields".to_string())
    }
}

impl<
    RS: EuclideanDomainSignature + FavoriteAssociateSignature,
    RSB: BorrowedStructure<RS>,
    const IS_FIELD: bool,
> SemiRingUnitsSignature for QuotientStructure<RS, RSB, IS_FIELD>
{
    fn inv(&self, x: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        if self.is_zero(x) {
            Err(RingDivisionError::DivideByZero)
        } else {
            let (g, a, b) = self.ring().euclidean_xgcd(x.clone(), self.modulus.clone());
            debug_assert!(
                self.ring().equal(
                    &g,
                    &self
                        .ring()
                        .add(&self.ring().mul(&a, x), &self.ring().mul(&b, &self.modulus))
                )
            );
            if self.equal(&g, &self.one()) {
                Ok(self.reduce(a))
            } else {
                Err(RingDivisionError::NotDivisible)
            }
        }
    }
}

impl<RS: EuclideanDomainSignature + FavoriteAssociateSignature, RSB: BorrowedStructure<RS>>
    IntegralDomainSignature for QuotientStructure<RS, RSB, true>
{
    fn div(&self, top: &Self::Set, bot: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        match self.inv(bot) {
            Ok(bot_inv) => Ok(self.mul(top, &bot_inv)),
            Err(err) => Err(err),
        }
    }
}

impl<RS: EuclideanDomainSignature + FavoriteAssociateSignature, RSB: BorrowedStructure<RS>>
    FieldSignature for QuotientStructure<RS, RSB, true>
{
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_nzq::*;

    #[test]
    fn test() {
        let mod5 = QuotientStructure::new_field(Integer::structure(), Integer::from(5)).unwrap();

        assert!(mod5.equal(
            &mod5.div(&Integer::from(3), &Integer::from(2)).unwrap(),
            &Integer::from(9)
        ));
    }
}
