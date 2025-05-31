use crate::structure::*;
use algebraeon_sets::structure::*;
use std::borrow::Borrow;

#[derive(Debug, Clone)]
pub struct QuotientStructure<RS: EuclideanDomainSignature, const IS_FIELD: bool> {
    ring: RS,
    modulus: RS::Set,
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> QuotientStructure<RS, IS_FIELD> {
    pub fn new_unchecked(ring: RS, modulus: RS::Set) -> Self {
        assert!(!ring.is_zero(&modulus));
        Self { ring, modulus }
    }

    pub fn ring(&self) -> &RS {
        &self.ring
    }

    pub fn modulus(&self) -> &RS::Set {
        &self.modulus
    }

    pub fn reduce(&self, a: impl Borrow<RS::Set>) -> RS::Set {
        self.ring.rem(a.borrow(), &self.modulus)
    }
}

impl<RS: EuclideanDomainSignature> QuotientStructure<RS, false> {
    pub fn new_ring(ring: RS, modulus: RS::Set) -> Self {
        Self::new_unchecked(ring, modulus)
    }
}

impl<RS: EuclideanDomainSignature> QuotientStructure<RS, true> {
    pub fn new_field_unchecked(ring: RS, modulus: RS::Set) -> Self {
        Self::new_unchecked(ring, modulus)
    }
}

impl<RS: EuclideanDomainSignature + FactorableSignature> QuotientStructure<RS, true> {
    pub fn new_field(ring: RS, modulus: RS::Set) -> Self {
        #[allow(clippy::manual_assert)]
        #[cfg(debug_assertions)]
        if !ring.is_irreducible(&modulus) {
            panic!(
                "The modulus must be irreducible to form a quotient field. Got {:?}.",
                modulus
            );
        }
        Self::new_unchecked(ring, modulus)
    }
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> PartialEq
    for QuotientStructure<RS, IS_FIELD>
{
    fn eq(&self, other: &Self) -> bool {
        if self.ring == other.ring {
            debug_assert_eq!(
                self.ring.equal(&self.modulus, &other.modulus),
                other.ring.equal(&self.modulus, &other.modulus)
            );
            self.ring.equal(&self.modulus, &other.modulus)
        } else {
            false
        }
    }
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> Eq for QuotientStructure<RS, IS_FIELD> {}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> Signature
    for QuotientStructure<RS, IS_FIELD>
{
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> SetSignature
    for QuotientStructure<RS, IS_FIELD>
{
    type Set = RS::Set;

    fn is_element(&self, _x: &Self::Set) -> bool {
        true
    }
}

impl<RS: EuclideanDomainSignature + ToStringSignature, const IS_FIELD: bool> ToStringSignature
    for QuotientStructure<RS, IS_FIELD>
{
    fn to_string(&self, elem: &Self::Set) -> String {
        self.ring.to_string(elem)
    }
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> EqSignature
    for QuotientStructure<RS, IS_FIELD>
{
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        self.ring.is_zero(
            &self
                .ring
                .rem(&self.ring.add(a, &self.ring.neg(b)), &self.modulus),
        )
    }
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> SemiRingSignature
    for QuotientStructure<RS, IS_FIELD>
{
    fn zero(&self) -> Self::Set {
        self.ring.zero()
    }

    fn one(&self) -> Self::Set {
        self.ring.one()
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.ring.rem(&self.ring.add(a, b), &self.modulus)
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.ring.rem(&self.ring.mul(a, b), &self.modulus)
    }
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> RingSignature
    for QuotientStructure<RS, IS_FIELD>
{
    fn neg(&self, a: &Self::Set) -> Self::Set {
        self.ring.neg(a)
    }
}

impl<RS: EuclideanDomainSignature + FavoriteAssociateSignature, const IS_FIELD: bool>
    SemiRingUnitsSignature for QuotientStructure<RS, IS_FIELD>
{
    fn inv(&self, x: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        if self.is_zero(x) {
            Err(RingDivisionError::DivideByZero)
        } else {
            let (g, a, b) = self.ring.euclidean_xgcd(x.clone(), self.modulus.clone());
            debug_assert!(
                self.ring.equal(
                    &g,
                    &self
                        .ring
                        .add(&self.ring.mul(&a, x), &self.ring.mul(&b, &self.modulus))
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

impl<RS: EuclideanDomainSignature + FavoriteAssociateSignature> IntegralDomainSignature
    for QuotientStructure<RS, true>
{
    fn div(&self, top: &Self::Set, bot: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        match self.inv(bot) {
            Ok(bot_inv) => Ok(self.mul(top, &bot_inv)),
            Err(err) => Err(err),
        }
    }
}

impl<RS: EuclideanDomainSignature + FavoriteAssociateSignature> FieldSignature
    for QuotientStructure<RS, true>
{
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_nzq::*;

    #[test]
    fn test() {
        let mod5 = QuotientStructure::new_field(Integer::structure(), Integer::from(5));

        assert!(mod5.equal(
            &mod5.div(&Integer::from(3), &Integer::from(2)).unwrap(),
            &Integer::from(9)
        ));
    }
}
