use crate::{
    algebraic_number_field::{
        AlgebraicIntegerRingSignature, AlgebraicNumberFieldSignature, QuadraticNumberFieldElement,
        QuadraticNumberFieldStructure,
    },
    structure::{
        AdditiveGroupSignature, AdditiveMonoidSignature, CharZeroRingSignature,
        CharacteristicSignature, DedekindDomainSignature, IntegralDomainSignature,
        MetaCharZeroRing, MultiplicativeMonoidSignature, MultiplicativeMonoidUnitsSignature,
        RingSignature, SemiModuleSignature, SemiRingSignature, SetWithZeroSignature,
    },
};
use algebraeon_nzq::{Integer, Natural, Rational};
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

impl<D: BorrowedSet<Integer>> SetWithZeroSignature for QuadraticRingOfIntegersStructure<D> {
    fn zero(&self) -> Self::Set {
        self.anf().zero()
    }
}

impl<D: BorrowedSet<Integer>> AdditiveMonoidSignature for QuadraticRingOfIntegersStructure<D> {
    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.anf().add(a, b)
    }

    fn try_neg(&self, a: &Self::Set) -> Option<Self::Set> {
        Some(self.neg(a))
    }

    fn try_sub(&self, a: &Self::Set, b: &Self::Set) -> Option<Self::Set> {
        Some(self.sub(a, b))
    }
}

impl<D: BorrowedSet<Integer>> AdditiveGroupSignature for QuadraticRingOfIntegersStructure<D> {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        self.anf().neg(a)
    }
}

impl<D: BorrowedSet<Integer>> MultiplicativeMonoidSignature
    for QuadraticRingOfIntegersStructure<D>
{
    fn one(&self) -> Self::Set {
        self.anf().one()
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.anf().mul(a, b)
    }
}

impl<D: BorrowedSet<Integer>> SemiRingSignature for QuadraticRingOfIntegersStructure<D> {}

impl<D: BorrowedSet<Integer>> RingSignature for QuadraticRingOfIntegersStructure<D> {}

impl<D: BorrowedSet<Integer>> MultiplicativeMonoidUnitsSignature
    for QuadraticRingOfIntegersStructure<D>
{
    fn try_inv(&self, a: &Self::Set) -> Option<Self::Set> {
        let b = self.anf().try_inv(a)?;
        if self.anf().is_algebraic_integer(&b) {
            Some(b)
        } else {
            None
        }
    }
}

impl<D: BorrowedSet<Integer>> IntegralDomainSignature for QuadraticRingOfIntegersStructure<D> {
    fn try_div(&self, a: &Self::Set, b: &Self::Set) -> Option<Self::Set> {
        let d = self.anf().try_div(a, b)?;
        if self.anf().is_algebraic_integer(&d) {
            Some(d)
        } else {
            None
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

    fn try_from_anf(&self, y: &QuadraticNumberFieldElement) -> Option<QuadraticNumberFieldElement> {
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

    fn to_anf(&self, x: &QuadraticNumberFieldElement) -> QuadraticNumberFieldElement {
        x.clone()
    }
    fn integral_basis(&self) -> Vec<Self::Set> {
        todo!()
    }
}
