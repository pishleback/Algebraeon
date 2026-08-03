use crate::structure::*;
use algebraeon_structures::*;
use std::{borrow::Cow, rc::Rc};

/// A quotient of a Euclidean domain by a non-zero element.
#[derive(Debug, Clone)]
pub struct EuclideanRemainderQuotientStructure<RS: EuclideanDomainSignature, const IS_FIELD: bool> {
    ring: Rc<RS>,
    modulus: RS::Elem,
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool>
    EuclideanRemainderQuotientStructure<RS, IS_FIELD>
{
    pub fn new_unchecked(ring: Rc<RS>, modulus: RS::Elem) -> Rc<Self> {
        debug_assert!(!ring.is_zero(&modulus));
        Self { ring, modulus }.into()
    }

    pub fn ring(&self) -> &Rc<RS> {
        &self.ring
    }

    // not a unique reduction
    pub fn reduce(&self, x: &RS::Elem) -> RS::Elem {
        self.ring().rem(x, &self.modulus)
    }
}

impl<RS: EuclideanDomainSignature> EuclideanRemainderQuotientStructure<RS, false> {
    fn try_new_ring(ring: Rc<RS>, modulus: RS::Elem) -> Option<Rc<Self>> {
        if ring.is_zero(&modulus) {
            None
        } else {
            Some(Self::new_unchecked(ring, modulus))
        }
    }
}

impl<RS: EuclideanDomainSignature + FactoringMonoidSignature>
    EuclideanRemainderQuotientStructure<RS, true>
{
    fn new_field_unchecked(ring: Rc<RS>, modulus: RS::Elem) -> Rc<Self> {
        debug_assert!(ring.is_irreducible(&modulus));
        Self::new_unchecked(ring, modulus)
    }

    fn try_new_field(ring: Rc<RS>, modulus: RS::Elem) -> Option<Rc<Self>> {
        if ring.is_zero(&modulus) || !ring.is_irreducible(&modulus) {
            None
        } else {
            Some(Self::new_unchecked(ring, modulus))
        }
    }
}

pub trait RingToQuotientRingSignature: EuclideanDomainSignature {
    fn euclidean_quotient_ring(
        self: &Rc<Self>,
        modulus: Self::Elem,
    ) -> Option<Rc<EuclideanRemainderQuotientStructure<Self, false>>> {
        EuclideanRemainderQuotientStructure::try_new_ring(self.clone(), modulus)
    }
}
impl<Ring: EuclideanDomainSignature> RingToQuotientRingSignature for Ring {}

pub trait RingToQuotientFieldSignature:
    EuclideanDomainSignature + FactoringMonoidSignature
{
    fn quotient_field(
        self: &Rc<Self>,
        modulus: Self::Elem,
    ) -> Option<Rc<EuclideanRemainderQuotientStructure<Self, true>>> {
        EuclideanRemainderQuotientStructure::try_new_field(self.clone(), modulus)
    }

    fn quotient_field_unchecked(
        self: &Rc<Self>,
        modulus: Self::Elem,
    ) -> Rc<EuclideanRemainderQuotientStructure<Self, true>> {
        EuclideanRemainderQuotientStructure::new_field_unchecked(self.clone(), modulus)
    }
}
impl<Ring: EuclideanDomainSignature + FactoringMonoidSignature> RingToQuotientFieldSignature
    for Ring
{
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> PartialEq
    for EuclideanRemainderQuotientStructure<RS, IS_FIELD>
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

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> Eq
    for EuclideanRemainderQuotientStructure<RS, IS_FIELD>
{
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> Signature
    for EuclideanRemainderQuotientStructure<RS, IS_FIELD>
{
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> SetSignature
    for EuclideanRemainderQuotientStructure<RS, IS_FIELD>
{
    type Elem = RS::Elem;

    fn validate_element(self: &Rc<Self>, _x: &Self::Elem) -> Result<(), String> {
        Ok(())
    }
}

impl<RS: EuclideanDomainSignature + ToStringSignature, const IS_FIELD: bool> ToStringSignature
    for EuclideanRemainderQuotientStructure<RS, IS_FIELD>
{
    fn to_string(self: &Rc<Self>, elem: &Self::Elem) -> String {
        self.ring().to_string(elem)
    }
}

impl<RS: EuclideanDomainSignature + FavoriteAssociateSignature, const IS_FIELD: bool>
    QuotientSetSignature<RS> for EuclideanRemainderQuotientStructure<RS, IS_FIELD>
{
    fn pre_quotient_set(self: &Rc<Self>) -> Rc<RS> {
        self.ring().clone()
    }

    fn project(self: &Rc<Self>, x: <RS as SetSignature>::Elem) -> Self::Elem {
        x
    }

    fn project_ref(self: &Rc<Self>, x: &<RS as SetSignature>::Elem) -> Self::Elem {
        x.clone()
    }

    fn unproject(self: &Rc<Self>, x: Self::Elem) -> <RS as SetSignature>::Elem {
        x
    }

    fn unproject_ref(self: &Rc<Self>, x: &Self::Elem) -> <RS as SetSignature>::Elem {
        x.clone()
    }
}

impl<RS: EuclideanDomainSignature + FavoriteAssociateSignature, const IS_FIELD: bool>
    QuotientRingSignature<RS> for EuclideanRemainderQuotientStructure<RS, IS_FIELD>
{
}

impl<RS: EuclideanDomainSignature + FavoriteAssociateSignature, const IS_FIELD: bool>
    QuotientRingGetPrincipalIdealSignature<RS>
    for EuclideanRemainderQuotientStructure<RS, IS_FIELD>
{
    fn modulus<'a>(&'a self) -> Cow<'a, RS::Elem> {
        Cow::Borrowed(&self.modulus)
    }
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> EqSignature
    for EuclideanRemainderQuotientStructure<RS, IS_FIELD>
{
    fn equal(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> bool {
        self.ring().is_zero(
            &self
                .ring()
                .rem(&self.ring().add(a, &self.ring().neg(b)), &self.modulus),
        )
    }
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> RinglikeSpecializationSignature
    for EuclideanRemainderQuotientStructure<RS, IS_FIELD>
{
    fn try_ring_restructure(
        self: Rc<Self>,
    ) -> Option<Rc<impl EqSignature<Elem = Self::Elem> + RingSignature>> {
        Some(self)
    }
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> ZeroSignature
    for EuclideanRemainderQuotientStructure<RS, IS_FIELD>
{
    fn zero(self: &Rc<Self>) -> Self::Elem {
        self.ring().zero()
    }
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> AdditionSignature
    for EuclideanRemainderQuotientStructure<RS, IS_FIELD>
{
    fn add(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.ring().rem(&self.ring().add(a, b), &self.modulus)
    }
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> CancellativeAdditionSignature
    for EuclideanRemainderQuotientStructure<RS, IS_FIELD>
{
    fn try_sub(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.sub(a, b))
    }
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> TryNegateSignature
    for EuclideanRemainderQuotientStructure<RS, IS_FIELD>
{
    fn try_neg(self: &Rc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.neg(a))
    }
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> AdditiveMonoidSignature
    for EuclideanRemainderQuotientStructure<RS, IS_FIELD>
{
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> AdditiveGroupSignature
    for EuclideanRemainderQuotientStructure<RS, IS_FIELD>
{
    fn neg(self: &Rc<Self>, a: &Self::Elem) -> Self::Elem {
        self.ring().neg(a)
    }

    fn sub(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.ring().sub(a, b)
    }
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> OneSignature
    for EuclideanRemainderQuotientStructure<RS, IS_FIELD>
{
    fn one(self: &Rc<Self>) -> Self::Elem {
        self.ring().one()
    }
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> MultiplicationSignature
    for EuclideanRemainderQuotientStructure<RS, IS_FIELD>
{
    fn mul(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.ring().rem(&self.ring().mul(a, b), &self.modulus)
    }
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> CommutativeMultiplicationSignature
    for EuclideanRemainderQuotientStructure<RS, IS_FIELD>
{
}

impl<RS: EuclideanDomainSignature + FavoriteAssociateSignature, const IS_FIELD: bool>
    TryReciprocalSignature for EuclideanRemainderQuotientStructure<RS, IS_FIELD>
{
    fn try_reciprocal(self: &Rc<Self>, x: &Self::Elem) -> Option<Self::Elem> {
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
                Some(self.reduce(&a))
            } else {
                None
            }
        }
    }
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> MultiplicativeMonoidSignature
    for EuclideanRemainderQuotientStructure<RS, IS_FIELD>
{
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> MultiplicativeAbsorptionMonoidSignature
    for EuclideanRemainderQuotientStructure<RS, IS_FIELD>
{
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> LeftDistributiveMultiplicationOverAddition
    for EuclideanRemainderQuotientStructure<RS, IS_FIELD>
{
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> RightDistributiveMultiplicationOverAddition
    for EuclideanRemainderQuotientStructure<RS, IS_FIELD>
{
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> SemiRingSignature
    for EuclideanRemainderQuotientStructure<RS, IS_FIELD>
{
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> RingSignature
    for EuclideanRemainderQuotientStructure<RS, IS_FIELD>
{
    fn is_reduced(self: &Rc<Self>) -> Result<bool, String> {
        if IS_FIELD {
            return Ok(true);
        }
        Err("not implemented yet: is_reduced for quotient rings that are not fields".to_string())
    }
}

impl<RS: EuclideanDomainSignature + FavoriteAssociateSignature> CancellativeMultiplicationSignature
    for EuclideanRemainderQuotientStructure<RS, true>
{
    fn try_divide(self: &Rc<Self>, top: &Self::Elem, bot: &Self::Elem) -> Option<Self::Elem> {
        Some(self.mul(top, &self.try_reciprocal(bot)?))
    }
}

impl<RS: EuclideanDomainSignature + FavoriteAssociateSignature>
    MultiplicativeIntegralMonoidSignature for EuclideanRemainderQuotientStructure<RS, true>
{
}

impl<RS: EuclideanDomainSignature + FavoriteAssociateSignature> IntegralDomainSignature
    for EuclideanRemainderQuotientStructure<RS, true>
{
}

impl<RS: EuclideanDomainSignature + FavoriteAssociateSignature> FieldSignature
    for EuclideanRemainderQuotientStructure<RS, true>
{
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test() {
        let mod5 = EuclideanRemainderQuotientStructure::try_new_field(
            Integer::structure(),
            Integer::from(5),
        )
        .unwrap();

        assert!(
            mod5.equal(
                &mod5
                    .try_divide(&Integer::from(3), &Integer::from(2))
                    .unwrap(),
                &Integer::from(9)
            )
        );
    }
}
