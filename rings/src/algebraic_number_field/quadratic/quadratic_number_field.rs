use crate::algebraic_number_field::QuadraticRingOfIntegersStructure;
use crate::algebraic_number_field::structure::{
    AlgebraicIntegerRingInAlgebraicNumberFieldSignature,
    RingOfIntegersToAlgebraicNumberFieldInclusion,
};
use crate::structure::{
    AdditiveMonoidEqSignature, FreeModuleSignature, MetaAdditiveMonoidEq, MetaCharZeroRing,
    MetaFactorableSignature, PrincipalRationalSubfieldInclusion, RingDivisionError,
    RingHomomorphismRangeModuleStructure, SemiModuleSignature,
};
use crate::{
    algebraic_number_field::structure::AlgebraicNumberFieldSignature,
    structure::{
        AdditiveGroupSignature, AdditiveMonoidSignature, CharZeroFieldSignature,
        CharZeroRingSignature, CharacteristicSignature, FieldSignature, IntegralDomainSignature,
        RingSignature, SemiRingSignature, SemiRingUnitsSignature,
    },
};
use algebraeon_nzq::{Integer, Natural, Rational, RationalCanonicalStructure};
use algebraeon_sets::structure::{BorrowedStructure, MetaType, Morphism};
use algebraeon_sets::structure::{
    CanonicalStructure, CountableSetSignature, EqSignature, FiniteSetSignature, SetSignature,
    Signature,
};
use std::borrow::Cow;

#[derive(Debug, Clone, Copy, PartialEq, Eq, CanonicalStructure)]
#[canonical_structure(eq)]
pub enum QuadraticNumberFieldBasis {
    Rational,
    Algebraic,
}

impl CountableSetSignature for QuadraticNumberFieldBasisCanonicalStructure {
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Set> + Clone {
        vec![
            QuadraticNumberFieldBasis::Rational,
            QuadraticNumberFieldBasis::Algebraic,
        ]
        .into_iter()
    }
}

impl FiniteSetSignature for QuadraticNumberFieldBasisCanonicalStructure {
    fn size(&self) -> usize {
        2
    }
}

#[derive(Debug, Clone)]
pub struct QuadraticNumberFieldElement {
    /// Represents `rational_part + algebraic_part * sqrt(d)`
    pub rational_part: Rational,
    pub algebraic_part: Rational,
}

impl QuadraticNumberFieldElement {
    const ZERO: Self = Self {
        rational_part: Rational::ZERO,
        algebraic_part: Rational::ZERO,
    };
    const ONE: Self = Self {
        rational_part: Rational::ONE,
        algebraic_part: Rational::ZERO,
    };
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QuadraticNumberFieldStructure<D: BorrowedStructure<Integer>> {
    // A squarefree integer
    d: D,
}

impl QuadraticNumberFieldStructure<Integer> {
    /// Given a squarefree integer `d`, return the quadratic number field `QQ[sqrt(d)]`.
    ///
    /// Returns an `Err` if `d` is not squarefree.
    pub fn new(d: Integer) -> Result<Self, ()> {
        if d.is_squarefree() {
            Ok(Self { d })
        } else {
            Err(())
        }
    }
}

impl<D: BorrowedStructure<Integer>> QuadraticNumberFieldStructure<D> {
    pub fn new_unchecked(d: D) -> Self {
        debug_assert!(d.borrow().is_squarefree());
        Self { d }
    }
}

impl<D: BorrowedStructure<Integer>> QuadraticNumberFieldStructure<D> {
    pub fn d(&self) -> &Integer {
        self.d.borrow()
    }

    pub fn roi<'d>(&'d self) -> QuadraticRingOfIntegersStructure<&'d Integer> {
        QuadraticRingOfIntegersStructure::new_unchecked(self.d())
    }
}

impl<D: BorrowedStructure<Integer>> Signature for QuadraticNumberFieldStructure<D> {}

impl<D: BorrowedStructure<Integer>> SetSignature for QuadraticNumberFieldStructure<D> {
    type Set = QuadraticNumberFieldElement;

    fn is_element(&self, _: &Self::Set) -> Result<(), String> {
        Ok(())
    }
}

impl<D: BorrowedStructure<Integer>> EqSignature for QuadraticNumberFieldStructure<D> {
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        a.rational_part == b.rational_part && a.algebraic_part == b.algebraic_part
    }
}

impl<D: BorrowedStructure<Integer>> AdditiveMonoidSignature for QuadraticNumberFieldStructure<D> {
    fn zero(&self) -> Self::Set {
        Self::Set::ZERO
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        QuadraticNumberFieldElement {
            rational_part: &a.rational_part + &b.rational_part,
            algebraic_part: &a.algebraic_part + &b.algebraic_part,
        }
    }
}

impl<D: BorrowedStructure<Integer>> AdditiveGroupSignature for QuadraticNumberFieldStructure<D> {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        QuadraticNumberFieldElement {
            rational_part: -&a.rational_part,
            algebraic_part: -&a.algebraic_part,
        }
    }
}

impl<D: BorrowedStructure<Integer>> SemiRingSignature for QuadraticNumberFieldStructure<D> {
    fn one(&self) -> Self::Set {
        Self::Set::ONE
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        // (x + y sqrd(d))(z + w sqrt(d)) = (xz + dyw) + (xw + yz) sqrt(d)
        QuadraticNumberFieldElement {
            rational_part: &a.rational_part * &b.rational_part
                + Rational::from(self.d()) * &a.algebraic_part * &b.algebraic_part,
            algebraic_part: &a.rational_part * &b.algebraic_part
                + &a.algebraic_part * &b.rational_part,
        }
    }
}

impl<D: BorrowedStructure<Integer>> RingSignature for QuadraticNumberFieldStructure<D> {}

impl<D: BorrowedStructure<Integer>> SemiRingUnitsSignature for QuadraticNumberFieldStructure<D> {
    fn inv(&self, a: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        // (x + y sqrt(d))^{-1} = (a - b sqrt(d)) / (x^2 + dy^2)
        debug_assert!(!self.d().is_zero()); // it's squarefree in particular non-zero
        let d = &a.rational_part * &a.rational_part
            + Rational::from(self.d()) * &a.algebraic_part * &a.algebraic_part;
        debug_assert_eq!(d == Rational::ZERO, self.is_zero(a));
        if d == Rational::ZERO {
            return Err(RingDivisionError::DivideByZero);
        }
        Ok(QuadraticNumberFieldElement {
            rational_part: &a.rational_part / &d,
            algebraic_part: -&a.algebraic_part / d,
        })
    }
}

impl<D: BorrowedStructure<Integer>> IntegralDomainSignature for QuadraticNumberFieldStructure<D> {
    fn div(&self, a: &Self::Set, b: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        Ok(self.mul(a, &self.inv(b)?))
    }
}

impl<D: BorrowedStructure<Integer>> CharacteristicSignature for QuadraticNumberFieldStructure<D> {
    fn characteristic(&self) -> Natural {
        Natural::ZERO
    }
}

impl<D: BorrowedStructure<Integer>> CharZeroRingSignature for QuadraticNumberFieldStructure<D> {
    fn try_to_int(&self, x: &Self::Set) -> Option<Integer> {
        if x.algebraic_part == Rational::ZERO {
            x.rational_part.try_to_int()
        } else {
            None
        }
    }
}

impl<D: BorrowedStructure<Integer>> FieldSignature for QuadraticNumberFieldStructure<D> {}

impl<D: BorrowedStructure<Integer>> CharZeroFieldSignature for QuadraticNumberFieldStructure<D> {
    fn try_to_rat(&self, x: &Self::Set) -> Option<Rational> {
        if x.algebraic_part == Rational::ZERO {
            Some(x.rational_part.clone())
        } else {
            None
        }
    }
}

impl<D: BorrowedStructure<Integer>> SemiModuleSignature<RationalCanonicalStructure>
    for QuadraticNumberFieldStructure<D>
{
    fn ring(&self) -> &RationalCanonicalStructure {
        Rational::structure_ref()
    }

    fn scalar_mul(&self, a: &Self::Set, x: &Rational) -> Self::Set {
        QuadraticNumberFieldElement {
            rational_part: x * &a.rational_part,
            algebraic_part: x * &a.algebraic_part,
        }
    }
}

impl<'h, D: BorrowedStructure<Integer>, B: BorrowedStructure<QuadraticNumberFieldStructure<D>>>
    FreeModuleSignature<RationalCanonicalStructure>
    for RingHomomorphismRangeModuleStructure<
        'h,
        RationalCanonicalStructure,
        QuadraticNumberFieldStructure<D>,
        PrincipalRationalSubfieldInclusion<QuadraticNumberFieldStructure<D>, B>,
    >
{
    type Basis = QuadraticNumberFieldBasisCanonicalStructure;

    fn basis_set(&self) -> impl std::borrow::Borrow<Self::Basis> {
        QuadraticNumberFieldBasis::structure()
    }

    fn to_component<'a>(
        &self,
        b: &QuadraticNumberFieldBasis,
        v: &'a QuadraticNumberFieldElement,
    ) -> std::borrow::Cow<'a, Rational> {
        match b {
            QuadraticNumberFieldBasis::Rational => Cow::Borrowed(&v.rational_part),
            QuadraticNumberFieldBasis::Algebraic => Cow::Borrowed(&v.algebraic_part),
        }
    }

    fn from_component(&self, b: &QuadraticNumberFieldBasis, r: &Rational) -> Self::Set {
        match b {
            QuadraticNumberFieldBasis::Rational => QuadraticNumberFieldElement {
                rational_part: r.clone(),
                algebraic_part: Rational::ZERO,
            },
            QuadraticNumberFieldBasis::Algebraic => QuadraticNumberFieldElement {
                rational_part: Rational::ZERO,
                algebraic_part: r.clone(),
            },
        }
    }
}

impl<
    D: BorrowedStructure<Integer>,
    RB: BorrowedStructure<QuadraticRingOfIntegersStructure<D>>,
    KB: BorrowedStructure<QuadraticNumberFieldStructure<D>>,
>
    RingOfIntegersToAlgebraicNumberFieldInclusion<
        QuadraticRingOfIntegersStructure<D>,
        RB,
        QuadraticNumberFieldStructure<D>,
        KB,
    >
{
    pub fn d(&self) -> &Integer {
        let d = self.anf().d();
        debug_assert_eq!(d, self.roi().d());
        d
    }
}

impl<
    D: BorrowedStructure<Integer>,
    RB: BorrowedStructure<QuadraticRingOfIntegersStructure<D>>,
    KB: BorrowedStructure<QuadraticNumberFieldStructure<D>>,
> AlgebraicIntegerRingInAlgebraicNumberFieldSignature
    for RingOfIntegersToAlgebraicNumberFieldInclusion<
        QuadraticRingOfIntegersStructure<D>,
        RB,
        QuadraticNumberFieldStructure<D>,
        KB,
    >
{
    type AlgebraicNumberField = QuadraticNumberFieldStructure<D>;
    type RingOfIntegers = QuadraticRingOfIntegersStructure<D>;

    fn discriminant(&self) -> Integer {
        let d_mod_4 = self.d() % Integer::from(4);
        if d_mod_4 == Integer::from(1) {
            self.d().clone()
        } else if d_mod_4 == Integer::from(2) || d_mod_4 == Integer::from(3) {
            Integer::from(4) * self.d()
        } else {
            unreachable!()
        }
    }

    fn try_anf_to_roi(
        &self,
        y: &QuadraticNumberFieldElement,
    ) -> Option<QuadraticNumberFieldElement> {
        let d_mod_4 = self.d() % Integer::from(4);
        if d_mod_4 == Integer::from(1) {
            // if d = 1 mod 4 then the ring of integers is {a + b * (1/2 + 1/2 sqrt(d)) : a, b in ZZ}
            let y2 = self.range().scalar_mul(y, &Rational::TWO);
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

    fn roi_to_anf(&self, x: &QuadraticNumberFieldElement) -> QuadraticNumberFieldElement {
        x.clone()
    }
}

impl<D: BorrowedStructure<Integer>> AlgebraicNumberFieldSignature
    for QuadraticNumberFieldStructure<D>
{
    type Basis = QuadraticNumberFieldBasisCanonicalStructure;
    type RingOfIntegers = QuadraticRingOfIntegersStructure<D>;
    type RationalInclusion<B: BorrowedStructure<Self>> =
        PrincipalRationalSubfieldInclusion<Self, B>;

    fn roi(&self) -> Self::RingOfIntegers {
        QuadraticRingOfIntegersStructure::new_unchecked(self.d.clone())
    }

    fn finite_dimensional_rational_extension<'a>(&'a self) -> Self::RationalInclusion<&'a Self> {
        PrincipalRationalSubfieldInclusion::new(self)
    }

    fn into_finite_dimensional_rational_extension(self) -> Self::RationalInclusion<Self> {
        PrincipalRationalSubfieldInclusion::new(self)
    }

    fn is_algebraic_integer(&self, a: &Self::Set) -> bool {
        self.ring_of_integers_extension()
            .try_anf_to_roi(a)
            .is_some()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

    // ZZ[i]
    #[test]
    fn qanf_neg_1() {
        let anf = QuadraticNumberFieldStructure::new(Integer::from(-1)).unwrap();

        let a = QuadraticNumberFieldElement {
            rational_part: Rational::from(1),
            algebraic_part: Rational::from(4),
        };
        let b = QuadraticNumberFieldElement {
            rational_part: Rational::from(2),
            algebraic_part: Rational::from(3),
        };

        assert!(anf.equal(
            &anf.add(&a, &b),
            &QuadraticNumberFieldElement {
                rational_part: Rational::from(3),
                algebraic_part: Rational::from(7)
            }
        ));

        assert!(anf.equal(
            &anf.neg(&a),
            &QuadraticNumberFieldElement {
                rational_part: Rational::from(-1),
                algebraic_part: Rational::from(-4)
            }
        ));

        assert!(anf.equal(
            &anf.mul(&a, &b),
            &QuadraticNumberFieldElement {
                rational_part: Rational::from(-10),
                algebraic_part: Rational::from(11)
            }
        ));

        assert_eq!(
            anf.clone().into_ring_of_integers_extension().discriminant(),
            Integer::from(-4)
        );

        assert!(anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from(0),
            algebraic_part: Rational::from(0)
        }));

        assert!(anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from(2),
            algebraic_part: Rational::from(3)
        }));

        assert!(!anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from_str("1/2").unwrap(),
            algebraic_part: Rational::from_str("1/2").unwrap(),
        }));

        assert!(!anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from_str("7/2").unwrap(),
            algebraic_part: Rational::from_str("-3/2").unwrap(),
        }));

        assert!(!anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from_str("1").unwrap(),
            algebraic_part: Rational::from_str("1/2").unwrap(),
        }));

        assert!(!anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from_str("7/2").unwrap(),
            algebraic_part: Rational::from_str("-6").unwrap(),
        }));
    }

    // ZZ[sqrt(2)]
    #[test]
    fn qanf_pos_2() {
        let anf = QuadraticNumberFieldStructure::new(Integer::from(2)).unwrap();

        let a = QuadraticNumberFieldElement {
            rational_part: Rational::from(1),
            algebraic_part: Rational::from(4),
        };
        let b = QuadraticNumberFieldElement {
            rational_part: Rational::from(2),
            algebraic_part: Rational::from(3),
        };

        assert!(anf.equal(
            &anf.add(&a, &b),
            &QuadraticNumberFieldElement {
                rational_part: Rational::from(3),
                algebraic_part: Rational::from(7)
            }
        ));

        assert!(anf.equal(
            &anf.neg(&a),
            &QuadraticNumberFieldElement {
                rational_part: Rational::from(-1),
                algebraic_part: Rational::from(-4)
            }
        ));

        assert!(anf.equal(
            &anf.mul(&a, &b),
            &QuadraticNumberFieldElement {
                rational_part: Rational::from(26),
                algebraic_part: Rational::from(11)
            }
        ));

        assert_eq!(
            anf.clone().into_ring_of_integers_extension().discriminant(),
            Integer::from(8)
        );

        assert!(anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from(0),
            algebraic_part: Rational::from(0)
        }));

        assert!(anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from(2),
            algebraic_part: Rational::from(3)
        }));

        assert!(!anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from_str("1/2").unwrap(),
            algebraic_part: Rational::from_str("1/2").unwrap(),
        }));

        assert!(!anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from_str("7/2").unwrap(),
            algebraic_part: Rational::from_str("-3/2").unwrap(),
        }));

        assert!(!anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from_str("1").unwrap(),
            algebraic_part: Rational::from_str("1/2").unwrap(),
        }));

        assert!(!anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from_str("7/2").unwrap(),
            algebraic_part: Rational::from_str("-6").unwrap(),
        }));
    }

    // ZZ[sqrt(3)]
    #[test]
    fn qanf_pos_3() {
        let anf = QuadraticNumberFieldStructure::new(Integer::from(3)).unwrap();

        let a = QuadraticNumberFieldElement {
            rational_part: Rational::from(1),
            algebraic_part: Rational::from(4),
        };
        let b = QuadraticNumberFieldElement {
            rational_part: Rational::from(2),
            algebraic_part: Rational::from(3),
        };

        assert!(anf.equal(
            &anf.mul(&a, &b),
            &QuadraticNumberFieldElement {
                rational_part: Rational::from(38),
                algebraic_part: Rational::from(11)
            }
        ));

        assert_eq!(
            anf.clone().into_ring_of_integers_extension().discriminant(),
            Integer::from(12)
        );

        assert!(anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from(0),
            algebraic_part: Rational::from(0)
        }));

        assert!(anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from(2),
            algebraic_part: Rational::from(3)
        }));

        assert!(!anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from_str("1/2").unwrap(),
            algebraic_part: Rational::from_str("1/2").unwrap(),
        }));

        assert!(!anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from_str("7/2").unwrap(),
            algebraic_part: Rational::from_str("-3/2").unwrap(),
        }));

        assert!(!anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from_str("1").unwrap(),
            algebraic_part: Rational::from_str("1/2").unwrap(),
        }));

        assert!(!anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from_str("7/2").unwrap(),
            algebraic_part: Rational::from_str("-6").unwrap(),
        }));
    }

    // ZZ[1/2 + 1/2 sqrt(5)]
    #[test]
    fn qanf_pos_5() {
        let anf = QuadraticNumberFieldStructure::new(Integer::from(5)).unwrap();

        let a = QuadraticNumberFieldElement {
            rational_part: Rational::from(1),
            algebraic_part: Rational::from(4),
        };
        let b = QuadraticNumberFieldElement {
            rational_part: Rational::from(2),
            algebraic_part: Rational::from(3),
        };

        assert!(anf.equal(
            &anf.mul(&a, &b),
            &QuadraticNumberFieldElement {
                rational_part: Rational::from(62),
                algebraic_part: Rational::from(11)
            }
        ));

        assert_eq!(
            anf.clone().into_ring_of_integers_extension().discriminant(),
            Integer::from(5)
        );

        assert!(anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from(0),
            algebraic_part: Rational::from(0)
        }));

        assert!(anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from(2),
            algebraic_part: Rational::from(3)
        }));

        assert!(anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from_str("1/2").unwrap(),
            algebraic_part: Rational::from_str("1/2").unwrap(),
        }));

        assert!(anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from_str("7/2").unwrap(),
            algebraic_part: Rational::from_str("-3/2").unwrap(),
        }));

        assert!(!anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from_str("1").unwrap(),
            algebraic_part: Rational::from_str("1/2").unwrap(),
        }));

        assert!(!anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from_str("7/2").unwrap(),
            algebraic_part: Rational::from_str("-6").unwrap(),
        }));
    }

    // ZZ[1/2 + 1/2 sqrt(-3)]
    #[test]
    fn qanf_neg_3() {
        let anf = QuadraticNumberFieldStructure::new(Integer::from(-3)).unwrap();

        let a = QuadraticNumberFieldElement {
            rational_part: Rational::from(1),
            algebraic_part: Rational::from(4),
        };
        let b = QuadraticNumberFieldElement {
            rational_part: Rational::from(2),
            algebraic_part: Rational::from(3),
        };

        assert!(anf.equal(
            &anf.mul(&a, &b),
            &QuadraticNumberFieldElement {
                rational_part: Rational::from(-34),
                algebraic_part: Rational::from(11)
            }
        ));

        assert_eq!(
            anf.clone().into_ring_of_integers_extension().discriminant(),
            Integer::from(-3)
        );

        assert!(anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from(0),
            algebraic_part: Rational::from(0)
        }));

        assert!(anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from(2),
            algebraic_part: Rational::from(3)
        }));

        assert!(anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from_str("1/2").unwrap(),
            algebraic_part: Rational::from_str("1/2").unwrap(),
        }));

        assert!(anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from_str("7/2").unwrap(),
            algebraic_part: Rational::from_str("-3/2").unwrap(),
        }));

        assert!(!anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from_str("1").unwrap(),
            algebraic_part: Rational::from_str("1/2").unwrap(),
        }));

        assert!(!anf.is_algebraic_integer(&QuadraticNumberFieldElement {
            rational_part: Rational::from_str("7/2").unwrap(),
            algebraic_part: Rational::from_str("-6").unwrap(),
        }));
    }
}
