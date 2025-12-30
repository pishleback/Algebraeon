use crate::{
    algebraic_number_field::{
        AlgebraicIntegerRingSignature, AlgebraicNumberFieldSignature, QuadraticNumberFieldElement,
        QuadraticNumberFieldStructure,
    },
    structure::{
        AdditiveGroupSignature, AdditiveMonoidSignature, CharZeroRingSignature,
        CharacteristicSignature, DedekindDomainSignature, IntegralDomainSignature,
        RingDivisionError, RingSignature, SemiRingSignature, SemiRingUnitsSignature,
    },
};
use algebraeon_nzq::{Integer, Natural};
use algebraeon_sets::structure::{BorrowedSet, EqSignature, SetSignature, Signature};

#[derive(Debug, Clone)]
pub struct QuadraticRingOfIntegersStructure<D: BorrowedSet<Integer>> {
    // A squarefree integer
    qanf: QuadraticNumberFieldStructure<D>,
}

impl<D: BorrowedSet<Integer>> PartialEq for QuadraticRingOfIntegersStructure<D> {
    fn eq(&self, other: &Self) -> bool {
        self.qanf == other.qanf
    }
}

impl<D: BorrowedSet<Integer>> Eq for QuadraticRingOfIntegersStructure<D> {}

impl QuadraticRingOfIntegersStructure<Integer> {
    /// Given a squarefree integer `d` other than 1, return the quadratic ring of integers for the number field `QQ[sqrt(d)]`.
    ///
    /// Returns an `Err` if `d` is not squarefree.
    pub fn new(d: Integer) -> Result<Self, ()> {
        Ok(Self {
            qanf: QuadraticNumberFieldStructure::new(d)?,
        })
    }
}

impl<D: BorrowedSet<Integer>> QuadraticRingOfIntegersStructure<D> {
    pub fn new_unchecked(d: D) -> Self {
        Self {
            qanf: QuadraticNumberFieldStructure::new_unchecked(d),
        }
    }
}

impl<D: BorrowedSet<Integer>> QuadraticRingOfIntegersStructure<D> {
    pub fn d(&self) -> &Integer {
        self.qanf.d()
    }

    pub fn anf<'d>(&'d self) -> QuadraticNumberFieldStructure<&'d Integer> {
        QuadraticNumberFieldStructure::new_unchecked(self.d())
    }
}

impl<D: BorrowedSet<Integer>> Signature for QuadraticRingOfIntegersStructure<D> {}

impl<D: BorrowedSet<Integer>> SetSignature for QuadraticRingOfIntegersStructure<D> {
    type Set = QuadraticNumberFieldElement;

    fn is_element(&self, a: &Self::Set) -> Result<(), String> {
        if self.anf().is_algebraic_integer(a) {
            Ok(())
        } else {
            Err("It is not an algebraic integer".to_string())
        }
    }
}

impl<D: BorrowedSet<Integer>> EqSignature for QuadraticRingOfIntegersStructure<D> {
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        self.anf().equal(a, b)
    }
}

impl<D: BorrowedSet<Integer>> AdditiveMonoidSignature for QuadraticRingOfIntegersStructure<D> {
    fn zero(&self) -> Self::Set {
        self.anf().zero()
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.anf().add(a, b)
    }
}

impl<D: BorrowedSet<Integer>> AdditiveGroupSignature for QuadraticRingOfIntegersStructure<D> {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        self.anf().neg(a)
    }
}

impl<D: BorrowedSet<Integer>> SemiRingSignature for QuadraticRingOfIntegersStructure<D> {
    fn one(&self) -> Self::Set {
        self.anf().one()
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.anf().mul(a, b)
    }
}

impl<D: BorrowedSet<Integer>> RingSignature for QuadraticRingOfIntegersStructure<D> {}

impl<D: BorrowedSet<Integer>> SemiRingUnitsSignature for QuadraticRingOfIntegersStructure<D> {
    fn inv(&self, a: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        let b = self.anf().inv(a)?;
        if self.anf().is_algebraic_integer(&b) {
            Ok(b)
        } else {
            Err(RingDivisionError::NotDivisible)
        }
    }
}

impl<D: BorrowedSet<Integer>> IntegralDomainSignature for QuadraticRingOfIntegersStructure<D> {
    fn div(&self, a: &Self::Set, b: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        let d = self.anf().div(a, b)?;
        if self.anf().is_algebraic_integer(&d) {
            Ok(d)
        } else {
            Err(RingDivisionError::NotDivisible)
        }
    }
}

impl<D: BorrowedSet<Integer>> CharacteristicSignature for QuadraticRingOfIntegersStructure<D> {
    fn characteristic(&self) -> Natural {
        Natural::ZERO
    }
}

impl<D: BorrowedSet<Integer>> CharZeroRingSignature for QuadraticRingOfIntegersStructure<D> {
    fn try_to_int(&self, x: &Self::Set) -> Option<Integer> {
        self.anf().try_to_int(x)
    }
}

impl<D: BorrowedSet<Integer>> DedekindDomainSignature for QuadraticRingOfIntegersStructure<D> {}

impl<D: BorrowedSet<Integer>> AlgebraicIntegerRingSignature<QuadraticNumberFieldStructure<D>>
    for QuadraticRingOfIntegersStructure<D>
{
    fn anf(&self) -> &QuadraticNumberFieldStructure<D> {
        &self.qanf
    }
}
