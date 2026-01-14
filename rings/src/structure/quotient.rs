use crate::structure::*;
use algebraeon_sets::structure::*;
use std::{borrow::Borrow, marker::PhantomData};

/// A quotient of a Euclidean domain by a non-zero element.
#[derive(Debug, Clone)]
pub struct EuclideanRemainderQuotientStructure<
    RS: EuclideanDomainSignature,
    RSB: BorrowedStructure<RS>,
    const IS_FIELD: bool,
> {
    _ring: PhantomData<RS>,
    ring: RSB,
    modulus: RS::Set,
}

impl<RS: EuclideanDomainSignature, RSB: BorrowedStructure<RS>, const IS_FIELD: bool>
    EuclideanRemainderQuotientStructure<RS, RSB, IS_FIELD>
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

impl<RS: EuclideanDomainSignature, RSB: BorrowedStructure<RS>>
    EuclideanRemainderQuotientStructure<RS, RSB, false>
{
    fn new_ring(ring: RSB, modulus: RS::Set) -> Self {
        Self::new_unchecked(ring, modulus)
    }
}

impl<RS: EuclideanDomainSignature + FactoringMonoidSignature, RSB: BorrowedStructure<RS>>
    EuclideanRemainderQuotientStructure<RS, RSB, true>
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
    fn quotient_ring(
        &self,
        modulus: Self::Set,
    ) -> EuclideanRemainderQuotientStructure<Self, &Self, false> {
        EuclideanRemainderQuotientStructure::new_ring(self, modulus)
    }

    fn into_quotient_ring(
        self,
        modulus: Self::Set,
    ) -> EuclideanRemainderQuotientStructure<Self, Self, false> {
        EuclideanRemainderQuotientStructure::new_ring(self, modulus)
    }
}
impl<Ring: EuclideanDomainSignature> RingToQuotientRingSignature for Ring {}

pub trait RingToQuotientFieldSignature:
    EuclideanDomainSignature + FactoringMonoidSignature
{
    fn quotient_field(
        &self,
        modulus: Self::Set,
    ) -> Result<EuclideanRemainderQuotientStructure<Self, &Self, true>, ()> {
        EuclideanRemainderQuotientStructure::new_field(self, modulus)
    }

    fn into_quotient_field(
        self,
        modulus: Self::Set,
    ) -> Result<EuclideanRemainderQuotientStructure<Self, Self, true>, ()> {
        EuclideanRemainderQuotientStructure::new_field(self, modulus)
    }

    fn quotient_field_unchecked(
        &self,
        modulus: Self::Set,
    ) -> EuclideanRemainderQuotientStructure<Self, &Self, true> {
        EuclideanRemainderQuotientStructure::new_field_unchecked(self, modulus)
    }

    fn into_quotient_field_unchecked(
        self,
        modulus: Self::Set,
    ) -> EuclideanRemainderQuotientStructure<Self, Self, true> {
        EuclideanRemainderQuotientStructure::new_field_unchecked(self, modulus)
    }
}
impl<Ring: EuclideanDomainSignature + FactoringMonoidSignature> RingToQuotientFieldSignature
    for Ring
{
}

impl<RS: EuclideanDomainSignature, RSB: BorrowedStructure<RS>, const IS_FIELD: bool> PartialEq
    for EuclideanRemainderQuotientStructure<RS, RSB, IS_FIELD>
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
    for EuclideanRemainderQuotientStructure<RS, RSB, IS_FIELD>
{
}

impl<RS: EuclideanDomainSignature, RSB: BorrowedStructure<RS>, const IS_FIELD: bool> Signature
    for EuclideanRemainderQuotientStructure<RS, RSB, IS_FIELD>
{
}

impl<RS: EuclideanDomainSignature, RSB: BorrowedStructure<RS>, const IS_FIELD: bool> SetSignature
    for EuclideanRemainderQuotientStructure<RS, RSB, IS_FIELD>
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
> ToStringSignature for EuclideanRemainderQuotientStructure<RS, RSB, IS_FIELD>
{
    fn to_string(&self, elem: &Self::Set) -> String {
        self.ring().to_string(elem)
    }
}

impl<RS: EuclideanDomainSignature, RSB: BorrowedStructure<RS>, const IS_FIELD: bool> EqSignature
    for EuclideanRemainderQuotientStructure<RS, RSB, IS_FIELD>
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
    RinglikeSpecializationSignature for EuclideanRemainderQuotientStructure<RS, RSB, IS_FIELD>
{
    fn try_ring_restructure(&self) -> Option<impl EqSignature<Set = Self::Set> + RingSignature> {
        Some(self.clone())
    }
}

impl<RS: EuclideanDomainSignature, RSB: BorrowedStructure<RS>, const IS_FIELD: bool> ZeroSignature
    for EuclideanRemainderQuotientStructure<RS, RSB, IS_FIELD>
{
    fn zero(&self) -> Self::Set {
        self.ring().zero()
    }
}

impl<RS: EuclideanDomainSignature, RSB: BorrowedStructure<RS>, const IS_FIELD: bool>
    AdditionSignature for EuclideanRemainderQuotientStructure<RS, RSB, IS_FIELD>
{
    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.ring().rem(&self.ring().add(a, b), &self.modulus)
    }
}

impl<RS: EuclideanDomainSignature, RSB: BorrowedStructure<RS>, const IS_FIELD: bool>
    CancellativeAdditionSignature for EuclideanRemainderQuotientStructure<RS, RSB, IS_FIELD>
{
    fn try_sub(&self, a: &Self::Set, b: &Self::Set) -> Option<Self::Set> {
        Some(self.sub(a, b))
    }
}

impl<RS: EuclideanDomainSignature, RSB: BorrowedStructure<RS>, const IS_FIELD: bool>
    TryNegateSignature for EuclideanRemainderQuotientStructure<RS, RSB, IS_FIELD>
{
    fn try_neg(&self, a: &Self::Set) -> Option<Self::Set> {
        Some(self.neg(a))
    }
}

impl<RS: EuclideanDomainSignature, RSB: BorrowedStructure<RS>, const IS_FIELD: bool>
    AdditiveMonoidSignature for EuclideanRemainderQuotientStructure<RS, RSB, IS_FIELD>
{
}

impl<RS: EuclideanDomainSignature, RSB: BorrowedStructure<RS>, const IS_FIELD: bool>
    AdditiveGroupSignature for EuclideanRemainderQuotientStructure<RS, RSB, IS_FIELD>
{
    fn neg(&self, a: &Self::Set) -> Self::Set {
        self.ring().neg(a)
    }

    fn sub(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.ring().sub(a, b)
    }
}

impl<RS: EuclideanDomainSignature, RSB: BorrowedStructure<RS>, const IS_FIELD: bool> OneSignature
    for EuclideanRemainderQuotientStructure<RS, RSB, IS_FIELD>
{
    fn one(&self) -> Self::Set {
        self.ring().one()
    }
}

impl<RS: EuclideanDomainSignature, RSB: BorrowedStructure<RS>, const IS_FIELD: bool>
    MultiplicationSignature for EuclideanRemainderQuotientStructure<RS, RSB, IS_FIELD>
{
    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.ring().rem(&self.ring().mul(a, b), &self.modulus)
    }
}

impl<RS: EuclideanDomainSignature, RSB: BorrowedStructure<RS>, const IS_FIELD: bool>
    MultiplicativeMonoidSignature for EuclideanRemainderQuotientStructure<RS, RSB, IS_FIELD>
{
}

impl<RS: EuclideanDomainSignature, RSB: BorrowedStructure<RS>, const IS_FIELD: bool>
    SemiRingSignature for EuclideanRemainderQuotientStructure<RS, RSB, IS_FIELD>
{
}

impl<RS: EuclideanDomainSignature, RSB: BorrowedStructure<RS>, const IS_FIELD: bool> RingSignature
    for EuclideanRemainderQuotientStructure<RS, RSB, IS_FIELD>
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
> MultiplicativeMonoidUnitsSignature for EuclideanRemainderQuotientStructure<RS, RSB, IS_FIELD>
{
    fn try_inv(&self, x: &Self::Set) -> Option<Self::Set> {
        if self.is_zero(x) {
            None
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
                Some(self.reduce(a))
            } else {
                None
            }
        }
    }
}

impl<RS: EuclideanDomainSignature + FavoriteAssociateSignature, RSB: BorrowedStructure<RS>>
    MultiplicativeIntegralMonoidSignature for EuclideanRemainderQuotientStructure<RS, RSB, true>
{
    fn try_div(&self, top: &Self::Set, bot: &Self::Set) -> Option<Self::Set> {
        Some(self.mul(top, &self.try_inv(bot)?))
    }
}

impl<RS: EuclideanDomainSignature + FavoriteAssociateSignature, RSB: BorrowedStructure<RS>>
    IntegralDomainSignature for EuclideanRemainderQuotientStructure<RS, RSB, true>
{
}

impl<RS: EuclideanDomainSignature + FavoriteAssociateSignature, RSB: BorrowedStructure<RS>>
    FieldSignature for EuclideanRemainderQuotientStructure<RS, RSB, true>
{
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_nzq::*;

    #[test]
    fn test() {
        let mod5 =
            EuclideanRemainderQuotientStructure::new_field(Integer::structure(), Integer::from(5))
                .unwrap();

        assert!(mod5.equal(
            &mod5.try_div(&Integer::from(3), &Integer::from(2)).unwrap(),
            &Integer::from(9)
        ));
    }
}
