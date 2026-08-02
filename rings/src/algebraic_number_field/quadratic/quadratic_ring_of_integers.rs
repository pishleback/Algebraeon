use std::sync::Arc;

use crate::{
    algebraic_number_field::{
        AlgebraicIntegerRingSignature, AlgebraicNumberFieldSignature, QuadraticNumberFieldElement,
        QuadraticNumberFieldStructure,
    },
    structure::{
        AdditionSignature, AdditiveGroupSignature, AdditiveMonoidSignature,
        CancellativeAdditionSignature, CancellativeMultiplicationSignature, CharZeroRingSignature,
        CharacteristicSignature, CommutativeMultiplicationSignature, DedekindDomainSignature,
        IntegralDomainSignature, LeftDistributiveMultiplicationOverAddition,
        MetaCharZeroRingSignature, MultiplicationSignature,
        MultiplicativeAbsorptionMonoidSignature, MultiplicativeIntegralMonoidSignature,
        MultiplicativeMonoidSignature, OneSignature, RightDistributiveMultiplicationOverAddition,
        RingSignature, RinglikeSpecializationSignature, SemiModuleSignature, SemiRingSignature,
        TryNegateSignature, TryReciprocalSignature, ZeroSignature,
    },
};
use algebraeon_structures::*;

#[derive(Debug, Clone)]
pub struct QuadraticRingOfIntegersStructure<D: BorrowedElem<Integer>> {
    // A squarefree integer
    qanf: Arc<QuadraticNumberFieldStructure<D>>,
}

impl<D: BorrowedElem<Integer>> PartialEq for QuadraticRingOfIntegersStructure<D> {
    fn eq(&self, other: &Self) -> bool {
        self.qanf == other.qanf
    }
}

impl<D: BorrowedElem<Integer>> Eq for QuadraticRingOfIntegersStructure<D> {}

impl QuadraticRingOfIntegersStructure<Integer> {
    /// Given a squarefree integer `d` other than 1, return the quadratic ring of integers for the number field `QQ[sqrt(d)]`.
    ///
    /// Returns an `Err` if `d` is not squarefree or is 1
    pub fn new(d: Integer) -> Result<Arc<Self>, String> {
        Ok(Self {
            qanf: QuadraticNumberFieldStructure::new(d)?,
        }
        .into())
    }
}

impl<D: BorrowedElem<Integer>> QuadraticRingOfIntegersStructure<D> {
    pub fn new_unchecked(d: D) -> Arc<Self> {
        Self {
            qanf: QuadraticNumberFieldStructure::new_unchecked(d),
        }
        .into()
    }
}

impl<D: BorrowedElem<Integer>> QuadraticRingOfIntegersStructure<D> {
    pub fn d(&self) -> &Integer {
        self.qanf.d()
    }

    pub fn anf(&self) -> Arc<QuadraticNumberFieldStructure<&Integer>> {
        QuadraticNumberFieldStructure::new_unchecked(self.d())
    }
}

impl<D: BorrowedElem<Integer>> Signature for QuadraticRingOfIntegersStructure<D> {}

impl<D: BorrowedElem<Integer>> SetSignature for QuadraticRingOfIntegersStructure<D> {
    type Elem = QuadraticNumberFieldElement;

    fn validate_element(self: &Arc<Self>, a: &Self::Elem) -> Result<(), String> {
        if self.anf().is_algebraic_integer(a) {
            Ok(())
        } else {
            Err("It is not an algebraic integer".to_string())
        }
    }
}

impl<D: BorrowedElem<Integer>> EqSignature for QuadraticRingOfIntegersStructure<D> {
    fn equal(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> bool {
        self.anf().equal(a, b)
    }
}

impl<D: BorrowedElem<Integer>> RinglikeSpecializationSignature
    for QuadraticRingOfIntegersStructure<D>
{
    fn try_ring_restructure(
        self: Arc<Self>,
    ) -> Option<Arc<impl EqSignature<Elem = Self::Elem> + RingSignature>> {
        Some(self)
    }

    fn try_char_zero_ring_restructure(
        self: Arc<Self>,
    ) -> Option<Arc<impl EqSignature<Elem = Self::Elem> + CharZeroRingSignature>> {
        Some(self)
    }
}

impl<D: BorrowedElem<Integer>> ZeroSignature for QuadraticRingOfIntegersStructure<D> {
    fn zero(self: &Arc<Self>) -> Self::Elem {
        self.anf().zero()
    }
}

impl<D: BorrowedElem<Integer>> AdditionSignature for QuadraticRingOfIntegersStructure<D> {
    fn add(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.anf().add(a, b)
    }
}

impl<D: BorrowedElem<Integer>> CancellativeAdditionSignature
    for QuadraticRingOfIntegersStructure<D>
{
    fn try_sub(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.sub(a, b))
    }
}

impl<D: BorrowedElem<Integer>> TryNegateSignature for QuadraticRingOfIntegersStructure<D> {
    fn try_neg(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.neg(a))
    }
}

impl<D: BorrowedElem<Integer>> AdditiveMonoidSignature for QuadraticRingOfIntegersStructure<D> {}

impl<D: BorrowedElem<Integer>> AdditiveGroupSignature for QuadraticRingOfIntegersStructure<D> {
    fn neg(self: &Arc<Self>, a: &Self::Elem) -> Self::Elem {
        self.anf().neg(a)
    }
}

impl<D: BorrowedElem<Integer>> OneSignature for QuadraticRingOfIntegersStructure<D> {
    fn one(self: &Arc<Self>) -> Self::Elem {
        self.anf().one()
    }
}

impl<D: BorrowedElem<Integer>> MultiplicationSignature for QuadraticRingOfIntegersStructure<D> {
    fn mul(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.anf().mul(a, b)
    }
}

impl<D: BorrowedElem<Integer>> CommutativeMultiplicationSignature
    for QuadraticRingOfIntegersStructure<D>
{
}

impl<D: BorrowedElem<Integer>> MultiplicativeMonoidSignature
    for QuadraticRingOfIntegersStructure<D>
{
}

impl<D: BorrowedElem<Integer>> MultiplicativeAbsorptionMonoidSignature
    for QuadraticRingOfIntegersStructure<D>
{
}

impl<D: BorrowedElem<Integer>> LeftDistributiveMultiplicationOverAddition
    for QuadraticRingOfIntegersStructure<D>
{
}

impl<D: BorrowedElem<Integer>> RightDistributiveMultiplicationOverAddition
    for QuadraticRingOfIntegersStructure<D>
{
}

impl<D: BorrowedElem<Integer>> SemiRingSignature for QuadraticRingOfIntegersStructure<D> {}

impl<D: BorrowedElem<Integer>> RingSignature for QuadraticRingOfIntegersStructure<D> {}

impl<D: BorrowedElem<Integer>> TryReciprocalSignature for QuadraticRingOfIntegersStructure<D> {
    fn try_reciprocal(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        let b = self.anf().try_reciprocal(a)?;
        if self.anf().is_algebraic_integer(&b) {
            Some(b)
        } else {
            None
        }
    }
}

impl<D: BorrowedElem<Integer>> CancellativeMultiplicationSignature
    for QuadraticRingOfIntegersStructure<D>
{
    fn try_divide(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        let d = self.anf().try_divide(a, b)?;
        if self.anf().is_algebraic_integer(&d) {
            Some(d)
        } else {
            None
        }
    }
}

impl<D: BorrowedElem<Integer>> MultiplicativeIntegralMonoidSignature
    for QuadraticRingOfIntegersStructure<D>
{
}

impl<D: BorrowedElem<Integer>> IntegralDomainSignature for QuadraticRingOfIntegersStructure<D> {}

impl<D: BorrowedElem<Integer>> CharacteristicSignature for QuadraticRingOfIntegersStructure<D> {
    fn characteristic(self: &Arc<Self>) -> Natural {
        Natural::ZERO
    }
}

impl<D: BorrowedElem<Integer>> CharZeroRingSignature for QuadraticRingOfIntegersStructure<D> {
    fn try_to_int(self: &Arc<Self>, x: &Self::Elem) -> Option<Integer> {
        self.anf().try_to_int(x)
    }
}

impl<D: BorrowedElem<Integer>> DedekindDomainSignature for QuadraticRingOfIntegersStructure<D> {}

impl<D: BorrowedElem<Integer>> AlgebraicIntegerRingSignature<QuadraticNumberFieldStructure<D>>
    for QuadraticRingOfIntegersStructure<D>
{
    fn anf(self: &Arc<Self>) -> Arc<QuadraticNumberFieldStructure<D>> {
        self.qanf.clone()
    }

    fn try_from_anf(
        self: &Arc<Self>,
        y: &QuadraticNumberFieldElement,
    ) -> Option<QuadraticNumberFieldElement> {
        let d_mod_4 = self.d() % Integer::from(4);
        if d_mod_4 == Integer::from(1) {
            // if d = 1 mod 4 then the ring of integers is {a + b * (1/2 + 1/2 sqrt(d)) : a, b in ZZ}
            let y2 = self.anf().scalar_mul(y, &Rational::TWO);
            if let Some(y2_rational_part_int) = y2.rational_part.try_to_int()
                && let Some(y2_algebraic_part_int) = y2.algebraic_part.try_to_int()
                && (y2_rational_part_int + y2_algebraic_part_int) % Integer::TWO == Integer::ZERO
            {
                Some(y.clone())
            } else {
                None
            }
        } else if d_mod_4 == Integer::from(2) || d_mod_4 == Integer::from(3) {
            // if d = 2 or 3 mod 4 then the ring of integers is {a + b * sqrt(d) : a, b in ZZ}
            if y.rational_part.is_integer() && y.algebraic_part.is_integer() {
                Some(y.clone())
            } else {
                None
            }
        } else {
            unreachable!()
        }
    }

    fn to_anf(self: &Arc<Self>, x: &QuadraticNumberFieldElement) -> QuadraticNumberFieldElement {
        x.clone()
    }

    fn integral_basis(self: &Arc<Self>) -> Vec<Self::Elem> {
        todo!()
    }
}
