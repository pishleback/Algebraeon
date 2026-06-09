use crate::algebraic_number_field::{
    AlgebraicIntegerRingSignature, AlgebraicNumberFieldSignature, QuadraticRingOfIntegersStructure,
};
use crate::structure::{
    AdditionSignature, AdditiveGroupSignature, AdditiveMonoidSignature,
    CancellativeAdditionSignature, CancellativeMultiplicationSignature, CharZeroFieldSignature,
    CharZeroRingSignature, CharacteristicSignature, CommutativeMultiplicationSignature,
    FieldSignature, IntegralDomainSignature, LeftDistributiveMultiplicationOverAddition,
    MetaFactoringMonoidNaturalExponent, MetaZeroEqSignature, MultiplicationSignature,
    MultiplicativeAbsorptionMonoidSignature, MultiplicativeIntegralMonoidSignature,
    MultiplicativeMonoidSignature, OneSignature, RightDistributiveMultiplicationOverAddition,
    RingSignature, RinglikeSpecializationSignature, SemiRingSignature, TryNegateSignature,
    TryReciprocalSignature, ZeroSignature,
};
use crate::structure::{
    FreeModuleSignature, MetaCharZeroRingSignature, PrincipalRationalMap,
    RingHomomorphismRangeModuleStructure, SemiModuleSignature, ZeroEqSignature,
};
use algebraeon_nzq::{Integer, Natural, Rational, RationalCanonicalStructure};
use algebraeon_sets::structure::{
    CanonicalStructure, CountableSetSignature, EqSignature, FiniteSetSignature,
};
use algebraeon_sets::structure::{InjectiveFunction};
use algebraeon_structures::*;
use std::borrow::Cow;

#[derive(Debug, Clone, Copy, PartialEq, Eq, CanonicalStructure)]
#[canonical_structure(eq)]
pub enum QuadraticNumberFieldBasis {
    Rational,
    Algebraic,
}

impl CountableSetSignature for QuadraticNumberFieldBasisCanonicalStructure {
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> + Clone {
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

#[derive(Debug, Clone)]
pub struct QuadraticNumberFieldStructure<D: BorrowedElem<Integer>> {
    // A squarefree integer. Not the discriminant.
    d: D,
}

impl<D: BorrowedElem<Integer>> PartialEq for QuadraticNumberFieldStructure<D> {
    fn eq(&self, other: &Self) -> bool {
        self.d.borrow() == other.d.borrow()
    }
}

impl<D: BorrowedElem<Integer>> Eq for QuadraticNumberFieldStructure<D> {}

impl QuadraticNumberFieldStructure<Integer> {
    /// Given a squarefree integer `d` other than 1, return the quadratic number field `QQ[sqrt(d)]`.
    ///
    /// Returns an `Err` if `d` is not squarefree or is 1
    pub fn new(d: Integer) -> Result<Self, String> {
        if d != Integer::ONE && d.is_squarefree() {
            Ok(Self { d })
        } else {
            Err(format!("`{d}` is not square-free or is equal to 1"))
        }
    }
}

impl<D: BorrowedElem<Integer>> QuadraticNumberFieldStructure<D> {
    pub fn new_unchecked(d: D) -> Self {
        debug_assert!(d.borrow() != &Integer::ONE && d.borrow().is_squarefree());
        Self { d }
    }
}

impl<D: BorrowedElem<Integer>> QuadraticNumberFieldStructure<D> {
    pub fn d(&self) -> &Integer {
        self.d.borrow()
    }

    pub fn roi(&self) -> QuadraticRingOfIntegersStructure<&Integer> {
        QuadraticRingOfIntegersStructure::new_unchecked(self.d())
    }
}

impl<D: BorrowedElem<Integer>> Signature for QuadraticNumberFieldStructure<D> {}

impl<D: BorrowedElem<Integer>> SetSignature for QuadraticNumberFieldStructure<D> {
    type Elem = QuadraticNumberFieldElement;

    fn validate_element(&self, _: &Self::Elem) -> Result<(), String> {
        Ok(())
    }
}

impl<D: BorrowedElem<Integer>> EqSignature for QuadraticNumberFieldStructure<D> {
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        a.rational_part == b.rational_part && a.algebraic_part == b.algebraic_part
    }
}

impl<D: BorrowedElem<Integer>> RinglikeSpecializationSignature
    for QuadraticNumberFieldStructure<D>
{
    fn try_ring_restructure(&self) -> Option<impl EqSignature<Elem = Self::Elem> + RingSignature> {
        Some(self.clone())
    }

    fn try_char_zero_ring_restructure(
        &self,
    ) -> Option<impl EqSignature<Elem = Self::Elem> + CharZeroRingSignature> {
        Some(self.clone())
    }
}

impl<D: BorrowedElem<Integer>> ZeroSignature for QuadraticNumberFieldStructure<D> {
    fn zero(&self) -> Self::Elem {
        Self::Elem::ZERO
    }
}

impl<D: BorrowedElem<Integer>> AdditionSignature for QuadraticNumberFieldStructure<D> {
    fn add(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        QuadraticNumberFieldElement {
            rational_part: &a.rational_part + &b.rational_part,
            algebraic_part: &a.algebraic_part + &b.algebraic_part,
        }
    }
}

impl<D: BorrowedElem<Integer>> CancellativeAdditionSignature for QuadraticNumberFieldStructure<D> {
    fn try_sub(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.sub(a, b))
    }
}

impl<D: BorrowedElem<Integer>> TryNegateSignature for QuadraticNumberFieldStructure<D> {
    fn try_neg(&self, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.neg(a))
    }
}

impl<D: BorrowedElem<Integer>> AdditiveMonoidSignature for QuadraticNumberFieldStructure<D> {}

impl<D: BorrowedElem<Integer>> AdditiveGroupSignature for QuadraticNumberFieldStructure<D> {
    fn neg(&self, a: &Self::Elem) -> Self::Elem {
        QuadraticNumberFieldElement {
            rational_part: -&a.rational_part,
            algebraic_part: -&a.algebraic_part,
        }
    }
}

impl<D: BorrowedElem<Integer>> OneSignature for QuadraticNumberFieldStructure<D> {
    fn one(&self) -> Self::Elem {
        Self::Elem::ONE
    }
}

impl<D: BorrowedElem<Integer>> MultiplicationSignature for QuadraticNumberFieldStructure<D> {
    fn mul(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        // (x + y sqrd(d))(z + w sqrt(d)) = (xz + dyw) + (xw + yz) sqrt(d)
        QuadraticNumberFieldElement {
            rational_part: &a.rational_part * &b.rational_part
                + Rational::from(self.d()) * &a.algebraic_part * &b.algebraic_part,
            algebraic_part: &a.rational_part * &b.algebraic_part
                + &a.algebraic_part * &b.rational_part,
        }
    }
}

impl<D: BorrowedElem<Integer>> CommutativeMultiplicationSignature
    for QuadraticNumberFieldStructure<D>
{
}

impl<D: BorrowedElem<Integer>> MultiplicativeMonoidSignature for QuadraticNumberFieldStructure<D> {}

impl<D: BorrowedElem<Integer>> MultiplicativeAbsorptionMonoidSignature
    for QuadraticNumberFieldStructure<D>
{
}

impl<D: BorrowedElem<Integer>> LeftDistributiveMultiplicationOverAddition
    for QuadraticNumberFieldStructure<D>
{
}

impl<D: BorrowedElem<Integer>> RightDistributiveMultiplicationOverAddition
    for QuadraticNumberFieldStructure<D>
{
}

impl<D: BorrowedElem<Integer>> SemiRingSignature for QuadraticNumberFieldStructure<D> {}

impl<D: BorrowedElem<Integer>> RingSignature for QuadraticNumberFieldStructure<D> {}

impl<D: BorrowedElem<Integer>> TryReciprocalSignature for QuadraticNumberFieldStructure<D> {
    fn try_reciprocal(&self, a: &Self::Elem) -> Option<Self::Elem> {
        // (x + y sqrt(d))^{-1} = (a - b sqrt(d)) / (x^2 + dy^2)
        debug_assert!(!self.d().is_zero()); // it's squarefree in particular non-zero
        let d = &a.rational_part * &a.rational_part
            + Rational::from(self.d()) * &a.algebraic_part * &a.algebraic_part;
        debug_assert_eq!(d == Rational::ZERO, self.is_zero(a));
        if d == Rational::ZERO {
            return None;
        }
        Some(QuadraticNumberFieldElement {
            rational_part: &a.rational_part / &d,
            algebraic_part: -&a.algebraic_part / d,
        })
    }
}

impl<D: BorrowedElem<Integer>> CancellativeMultiplicationSignature
    for QuadraticNumberFieldStructure<D>
{
    fn try_divide(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.mul(a, &self.try_reciprocal(b)?))
    }
}

impl<D: BorrowedElem<Integer>> MultiplicativeIntegralMonoidSignature
    for QuadraticNumberFieldStructure<D>
{
}

impl<D: BorrowedElem<Integer>> IntegralDomainSignature for QuadraticNumberFieldStructure<D> {}

impl<D: BorrowedElem<Integer>> CharacteristicSignature for QuadraticNumberFieldStructure<D> {
    fn characteristic(&self) -> Natural {
        Natural::ZERO
    }
}

impl<D: BorrowedElem<Integer>> CharZeroRingSignature for QuadraticNumberFieldStructure<D> {
    fn try_to_int(&self, x: &Self::Elem) -> Option<Integer> {
        if x.algebraic_part == Rational::ZERO {
            x.rational_part.try_to_int()
        } else {
            None
        }
    }
}

impl<D: BorrowedElem<Integer>> FieldSignature for QuadraticNumberFieldStructure<D> {}

impl<D: BorrowedElem<Integer>> CharZeroFieldSignature for QuadraticNumberFieldStructure<D> {
    fn try_to_rat(&self, x: &Self::Elem) -> Option<Rational> {
        if x.algebraic_part == Rational::ZERO {
            Some(x.rational_part.clone())
        } else {
            None
        }
    }
}

impl<D: BorrowedElem<Integer>> SemiModuleSignature<RationalCanonicalStructure>
    for QuadraticNumberFieldStructure<D>
{
    fn ring(&self) -> &RationalCanonicalStructure {
        Rational::structure_ref()
    }

    fn scalar_mul(&self, a: &Self::Elem, x: &Rational) -> Self::Elem {
        QuadraticNumberFieldElement {
            rational_part: x * &a.rational_part,
            algebraic_part: x * &a.algebraic_part,
        }
    }
}

impl<'h, D: BorrowedElem<Integer>, B: BorrowedStructure<QuadraticNumberFieldStructure<D>>>
    FreeModuleSignature<RationalCanonicalStructure>
    for RingHomomorphismRangeModuleStructure<
        'h,
        RationalCanonicalStructure,
        QuadraticNumberFieldStructure<D>,
        PrincipalRationalMap<QuadraticNumberFieldStructure<D>, B>,
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
    ) -> Cow<'a, Rational> {
        match b {
            QuadraticNumberFieldBasis::Rational => Cow::Borrowed(&v.rational_part),
            QuadraticNumberFieldBasis::Algebraic => Cow::Borrowed(&v.algebraic_part),
        }
    }

    fn from_component(&self, b: &QuadraticNumberFieldBasis, r: &Rational) -> Self::Elem {
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

// impl<D: BorrowedSet<Integer>, RB: BorrowedStructure<QuadraticRingOfIntegersStructure<D>>>
//     RingOfIntegersToAlgebraicNumberFieldInclusion<
//         QuadraticNumberFieldStructure<D>,
//         QuadraticRingOfIntegersStructure<D>,
//         RB,
//     >
// {
//     pub fn d(&self) -> &Integer {
//         let d = self.anf().d();
//         debug_assert_eq!(d, self.roi().d());
//         d
//     }
// }

impl<D: BorrowedElem<Integer>> QuadraticNumberFieldStructure<D> {
    pub fn into_ring_of_integers(self) -> QuadraticRingOfIntegersStructure<D> {
        QuadraticRingOfIntegersStructure::new_unchecked(self.d)
    }

    pub fn ring_of_integers(&self) -> QuadraticRingOfIntegersStructure<&Integer> {
        QuadraticRingOfIntegersStructure::new_unchecked(self.d.borrow())
    }
}

impl<D: BorrowedElem<Integer>> AlgebraicNumberFieldSignature for QuadraticNumberFieldStructure<D> {
    type Basis = QuadraticNumberFieldBasisCanonicalStructure;
    type RationalInclusion<B: BorrowedStructure<Self>> = PrincipalRationalMap<Self, B>;

    fn inbound_finite_dimensional_rational_extension(&self) -> Self::RationalInclusion<&Self> {
        PrincipalRationalMap::new(self)
    }

    fn into_inbound_finite_dimensional_rational_extension(self) -> Self::RationalInclusion<Self> {
        PrincipalRationalMap::new(self)
    }

    fn generator(&self) -> Self::Elem {
        QuadraticNumberFieldElement {
            rational_part: Rational::ZERO,
            algebraic_part: Rational::ONE,
        }
    }

    fn discriminant(&self) -> Integer {
        let d_mod_4 = self.d() % Integer::from(4);
        if d_mod_4 == Integer::from(1) {
            // d if d = 1 (mod 4)
            self.d().clone()
        } else if d_mod_4 == Integer::from(2) || d_mod_4 == Integer::from(3) {
            // 4d if d = 2, 3 (mod 4)
            Integer::from(4) * self.d()
        } else {
            // d != 0 (mod 4) since d is squarefree
            unreachable!()
        }
    }

    fn integral_basis(&self) -> Vec<Self::Elem> {
        todo!()
    }

    fn is_algebraic_integer(&self, a: &Self::Elem) -> bool {
        self.ring_of_integers()
            .outbound_roi_to_anf_inclusion()
            .try_preimage(a)
            .is_some()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

    // ZZ[1]
    #[test]
    fn qanf_pos_1() {
        let anf = QuadraticNumberFieldStructure::new(Integer::from(1));
        debug_assert!(anf.is_err());
    }

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

        assert_eq!(anf.discriminant(), Integer::from(-4));

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

        assert_eq!(anf.discriminant(), Integer::from(8));

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

        assert_eq!(anf.discriminant(), Integer::from(12));

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

        assert_eq!(anf.discriminant(), Integer::from(5));

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

        assert_eq!(anf.discriminant(), Integer::from(-3));

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
