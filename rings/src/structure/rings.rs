use super::*;
use crate::polynomial::*;
use algebraeon_macros::{signature_meta_trait, skip_meta};
use algebraeon_nzq::{Integer, Natural, NaturalCanonicalStructure, Rational, traits::*};
use algebraeon_sets::structure::*;
use algebraeon_structures::*;
use std::{borrow::Borrow, fmt::Debug};

mod unconstructable_universal_structure {
    use crate::structure::{
        AdditionSignature, AdditiveGroupSignature, AdditiveMonoidSignature,
        CancellativeAdditionSignature, CharZeroRingSignature, CharacteristicSignature,
        CommutativeMultiplicationSignature, LeftDistributiveMultiplicationOverAddition,
        MultiplicationSignature, MultiplicativeAbsorptionMonoidSignature,
        MultiplicativeMonoidSignature, OneSignature, RightDistributiveMultiplicationOverAddition,
        RingSignature, RinglikeSpecializationSignature, SemiRingSignature, TryNegateSignature,
        TryReciprocalSignature, ZeroSignature,
    };
    use algebraeon_sets::structure::EqSignature;
    use algebraeon_structures::*;
    use std::fmt::Debug;
    use std::marker::PhantomData;

    pub struct UnconstructableStructure<Set> {
        _set: PhantomData<Set>,
    }

    unsafe impl<Set> Send for UnconstructableStructure<Set> {}

    unsafe impl<Set> Sync for UnconstructableStructure<Set> {}

    impl<Set> Debug for UnconstructableStructure<Set> {
        fn fmt(&self, _f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            unreachable!()
        }
    }

    impl<Set> Clone for UnconstructableStructure<Set> {
        fn clone(&self) -> Self {
            unreachable!()
        }
    }

    impl<Set> PartialEq for UnconstructableStructure<Set> {
        fn eq(&self, _other: &Self) -> bool {
            unreachable!()
        }
    }

    impl<Set> Eq for UnconstructableStructure<Set> {}

    impl<Set> Signature for UnconstructableStructure<Set> {}

    impl<Set: Debug + Clone + Send + Sync> SetSignature for UnconstructableStructure<Set> {
        type Elem = Set;

        fn validate_element(&self, _x: &Self::Elem) -> Result<(), String> {
            unreachable!()
        }
    }

    impl<Set: Debug + Clone + Send + Sync> EqSignature for UnconstructableStructure<Set> {
        fn equal(&self, _a: &Self::Elem, _b: &Self::Elem) -> bool {
            unreachable!()
        }
    }

    impl<Set: Debug + Clone + Send + Sync> RinglikeSpecializationSignature
        for UnconstructableStructure<Set>
    {
    }

    impl<Set: Debug + Clone + Send + Sync> ZeroSignature for UnconstructableStructure<Set> {
        fn zero(&self) -> Self::Elem {
            unreachable!()
        }
    }

    impl<Set: Debug + Clone + Send + Sync> AdditionSignature for UnconstructableStructure<Set> {
        fn add(&self, _a: &Self::Elem, _b: &Self::Elem) -> Self::Elem {
            unreachable!()
        }
    }

    impl<Set: Debug + Clone + Send + Sync> CancellativeAdditionSignature
        for UnconstructableStructure<Set>
    {
        fn try_sub(&self, _a: &Self::Elem, _b: &Self::Elem) -> Option<Self::Elem> {
            unreachable!()
        }
    }

    impl<Set: Debug + Clone + Send + Sync> TryNegateSignature for UnconstructableStructure<Set> {
        fn try_neg(&self, _a: &Self::Elem) -> Option<Self::Elem> {
            unreachable!()
        }
    }

    impl<Set: Debug + Clone + Send + Sync> AdditiveMonoidSignature for UnconstructableStructure<Set> {}

    impl<Set: Debug + Clone + Send + Sync> AdditiveGroupSignature for UnconstructableStructure<Set> {
        fn neg(&self, _a: &Self::Elem) -> Self::Elem {
            unreachable!()
        }
    }

    impl<Set: Debug + Clone + Send + Sync> OneSignature for UnconstructableStructure<Set> {
        fn one(&self) -> Self::Elem {
            unreachable!()
        }
    }

    impl<Set: Debug + Clone + Send + Sync> MultiplicationSignature for UnconstructableStructure<Set> {
        fn mul(&self, _a: &Self::Elem, _b: &Self::Elem) -> Self::Elem {
            unreachable!()
        }
    }

    impl<Set: Debug + Clone + Send + Sync> CommutativeMultiplicationSignature
        for UnconstructableStructure<Set>
    {
    }

    impl<Set: Debug + Clone + Send + Sync> TryReciprocalSignature for UnconstructableStructure<Set> {
        fn try_reciprocal(&self, _a: &Self::Elem) -> Option<Self::Elem> {
            unreachable!()
        }
    }

    impl<Set: Debug + Clone + Send + Sync> MultiplicativeMonoidSignature
        for UnconstructableStructure<Set>
    {
    }

    impl<Set: Debug + Clone + Send + Sync> MultiplicativeAbsorptionMonoidSignature
        for UnconstructableStructure<Set>
    {
    }

    impl<Set: Debug + Clone + Send + Sync> LeftDistributiveMultiplicationOverAddition
        for UnconstructableStructure<Set>
    {
    }

    impl<Set: Debug + Clone + Send + Sync> RightDistributiveMultiplicationOverAddition
        for UnconstructableStructure<Set>
    {
    }

    impl<Set: Debug + Clone + Send + Sync> SemiRingSignature for UnconstructableStructure<Set> {}

    impl<Set: Debug + Clone + Send + Sync> CharacteristicSignature for UnconstructableStructure<Set> {
        fn characteristic(&self) -> algebraeon_nzq::Natural {
            unreachable!()
        }
    }

    impl<Set: Debug + Clone + Send + Sync> RingSignature for UnconstructableStructure<Set> {}

    impl<Set: Debug + Clone + Send + Sync> CharZeroRingSignature for UnconstructableStructure<Set> {
        fn try_to_int(&self, _x: &Self::Elem) -> Option<algebraeon_nzq::Integer> {
            unreachable!()
        }
    }
}

/*
All methods are allowed to return None, and if they do it should only affect speed of algorithms, not correctness.
Where methods do not return None, the resulting structure may be used for optimizations only.

The methods currently require cloning some structures when sucessful.
That's because it's currently prohibitively messy to allow returning by value or by reference depending on context.
Stabilisation of Rust features such as impl aliases may make this more feasible in future.
*/
pub trait RinglikeSpecializationSignature: SetSignature + ToOwned<Owned = Self> {
    /*
    Used by:
     - Polynomial rings to determine whether the karatsuba is usable.
     */
    fn try_ring_restructure(&self) -> Option<impl EqSignature<Elem = Self::Elem> + RingSignature> {
        Option::<unconstructable_universal_structure::UnconstructableStructure<Self::Elem>>::None
    }

    /*
    Used by:
     - Formatting polynomials as strings: If the set of coefficients has this structure then it's possible to call .try_to_int(..) which can allow for nicer formatting at integer coefficients.
     */
    fn try_char_zero_ring_restructure(
        &self,
    ) -> Option<impl EqSignature<Elem = Self::Elem> + CharZeroRingSignature> {
        Option::<unconstructable_universal_structure::UnconstructableStructure<Self::Elem>>::None
    }
}

/// A set with a special element `0`.
#[signature_meta_trait]
pub trait ZeroSignature: RinglikeSpecializationSignature {
    fn zero(&self) -> Self::Elem;
}

/// A set with an associative commutative binary operation of addition.
#[signature_meta_trait]
pub trait AdditionSignature: RinglikeSpecializationSignature {
    fn add(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem;

    fn add_mut(&self, a: &mut Self::Elem, b: &Self::Elem) {
        *a = self.add(a, b);
    }
}

/// When `a + x` = `a + y` implies `x` = `y` for all `a`, `x`, `y`.
#[signature_meta_trait]
pub trait CancellativeAdditionSignature: AdditionSignature {
    /// Return the unique `x` such that `a` = `b + x`, or `None` if no such `x` exists.
    fn try_sub(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem>;
}

#[signature_meta_trait]
pub trait TryNegateSignature: ZeroSignature + AdditionSignature {
    fn try_neg(&self, a: &Self::Elem) -> Option<Self::Elem>;
}

#[signature_meta_trait]
pub trait AdditiveMonoidSignature: ZeroSignature + AdditionSignature + TryNegateSignature {
    fn sum(&self, vals: &[impl Borrow<Self::Elem>]) -> Self::Elem {
        let mut sum = self.zero();
        for val in vals {
            self.add_mut(&mut sum, val.borrow());
        }
        sum
    }
}

#[signature_meta_trait]
pub trait AdditiveGroupSignature: AdditiveMonoidSignature + CancellativeAdditionSignature {
    fn neg(&self, a: &Self::Elem) -> Self::Elem;

    fn sub(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.add(a, &self.neg(b))
    }

    fn sub_mut(&self, a: &mut Self::Elem, b: &Self::Elem) {
        *a = self.sub(a, b);
    }
}

#[signature_meta_trait]
pub trait ZeroEqSignature: ZeroSignature + EqSignature {
    fn is_zero(&self, a: &Self::Elem) -> bool {
        self.equal(a, &self.zero())
    }
}
impl<R: ZeroSignature + EqSignature> ZeroEqSignature for R {}

/// A set with a special element `1`.
#[signature_meta_trait]
pub trait OneSignature: RinglikeSpecializationSignature {
    fn one(&self) -> Self::Elem;
}

/// A set with an associative binary opperation `*`.
#[signature_meta_trait]
pub trait MultiplicationSignature: RinglikeSpecializationSignature {
    fn mul(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem;

    fn mul_mut(&self, a: &mut Self::Elem, b: &Self::Elem) {
        *a = self.mul(a, b);
    }
}

/// When `*` is commutative.
#[signature_meta_trait]
pub trait CommutativeMultiplicationSignature: MultiplicationSignature {}

/// When `1 * a` = `a * 1` = `a` for all `a`.
#[signature_meta_trait]
pub trait MultiplicativeMonoidSignature: OneSignature + MultiplicationSignature {
    fn product(&self, vals: &[impl Borrow<Self::Elem>]) -> Self::Elem {
        let mut prod = self.one();
        for val in vals {
            self.mul_mut(&mut prod, val.borrow());
        }
        prod
    }

    fn nat_pow(&self, a: &Self::Elem, n: &Natural) -> Self::Elem {
        if *n == Natural::ZERO {
            self.one()
        } else if *n == Natural::ONE {
            a.clone()
        } else {
            debug_assert!(*n >= Natural::TWO);
            let bits: Vec<_> = n.bits().collect();
            let mut pows = vec![a.clone()];
            while pows.len() < bits.len() {
                pows.push(self.mul(pows.last().unwrap(), pows.last().unwrap()));
            }
            let count = bits.len();
            debug_assert_eq!(count, pows.len());
            let mut ans = self.one();
            for i in 0..count {
                if bits[i] {
                    self.mul_mut(&mut ans, &pows[i]);
                }
            }
            ans
        }
    }
}

#[signature_meta_trait]
pub trait TryLeftReciprocalSignature: OneSignature + MultiplicationSignature {
    /// `x` such that `x*a`=`1` or `None` if no such `x` exists.
    fn try_left_reciprocal(&self, a: &Self::Elem) -> Option<Self::Elem>;
}

#[signature_meta_trait]
pub trait TryRightReciprocalSignature: OneSignature + MultiplicationSignature {
    /// `x` such that `a*x`=`1` or `None` if no such `x` exists.
    fn try_right_reciprocal(&self, a: &Self::Elem) -> Option<Self::Elem>;
}

#[signature_meta_trait]
pub trait TryReciprocalSignature: OneSignature + MultiplicationSignature {
    /// `b` such that `a*b`=`1` and `b*a`=`1` or `None` if no such `b` exists.
    fn try_reciprocal(&self, a: &Self::Elem) -> Option<Self::Elem>;

    fn is_unit(&self, a: &Self::Elem) -> bool {
        self.try_reciprocal(a).is_some()
    }

    #[skip_meta]
    fn units(&self) -> MultiplicativeMonoidUnitsStructure<Self, &Self> {
        MultiplicativeMonoidUnitsStructure::new(self)
    }

    #[skip_meta]
    fn into_units(self) -> MultiplicativeMonoidUnitsStructure<Self, Self> {
        MultiplicativeMonoidUnitsStructure::new(self)
    }
}

impl<S: TryReciprocalSignature + CommutativeMultiplicationSignature> TryLeftReciprocalSignature
    for S
{
    fn try_left_reciprocal(&self, a: &Self::Elem) -> Option<Self::Elem> {
        self.try_reciprocal(a)
    }
}

impl<S: TryReciprocalSignature + CommutativeMultiplicationSignature> TryRightReciprocalSignature
    for S
{
    fn try_right_reciprocal(&self, a: &Self::Elem) -> Option<Self::Elem> {
        self.try_reciprocal(a)
    }
}

#[signature_meta_trait]
pub trait MultiplicativeMonoidTryInverseSignature:
    MultiplicativeMonoidSignature + TryReciprocalSignature
{
    fn try_int_pow(&self, a: &Self::Elem, n: &Integer) -> Option<Self::Elem> {
        if *n == Integer::ZERO {
            Some(self.one())
        } else if *n > Integer::ZERO {
            Some(self.nat_pow(a, &Abs::abs(n)))
        } else {
            Some(self.nat_pow(&self.try_reciprocal(a)?, &Abs::abs(n)))
        }
    }
}
impl<S: MultiplicativeMonoidSignature + TryReciprocalSignature>
    MultiplicativeMonoidTryInverseSignature for S
{
}

#[signature_meta_trait]
pub trait MultiplicativeMonoidSquareOpsSignature: MultiplicativeMonoidSignature {
    fn is_square(&self, a: &Self::Elem) -> bool {
        self.sqrt_if_square(a).is_some()
    }

    fn sqrt_if_square(&self, a: &Self::Elem) -> Option<Self::Elem>;
}

/// 0 is such that `a*0` = `0*a` = `0` for all a in the monoid.
/// such an element is unqiue if it exists.
#[signature_meta_trait]
pub trait MultiplicativeAbsorptionMonoidSignature:
    MultiplicativeMonoidSignature + ZeroSignature
{
}

/// When `a*(b + c) = a*b + a*c`.
#[signature_meta_trait]
pub trait LeftDistributiveMultiplicationOverAddition:
    AdditionSignature + MultiplicationSignature
{
}

/// When `(a + b)*c = a*c + b*c`.
#[signature_meta_trait]
pub trait RightDistributiveMultiplicationOverAddition:
    AdditionSignature + MultiplicationSignature
{
}

/// When `a * x` = `a * y` implies `x` = `y` for all `a`, `x`, `y`.
#[signature_meta_trait]
pub trait LeftCancellativeMultiplicationSignature: MultiplicationSignature {
    /// Try to find `x` such that `a` = `b * x`.
    fn try_left_divide(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem>;
}

/// When `x * a` = `y * a` implies `x` = `y` for all `a`, `x`, `y`.
#[signature_meta_trait]
pub trait RightCancellativeMultiplicationSignature: MultiplicationSignature {
    /// Try to find `x` such that `a` = `x * b`.
    fn try_right_divide(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem>;
}

#[signature_meta_trait]
pub trait CancellativeMultiplicationSignature:
    CommutativeMultiplicationSignature
    + LeftCancellativeMultiplicationSignature
    + RightCancellativeMultiplicationSignature
{
    fn try_divide(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem>;

    /// return true iff a is divisible by b
    fn divisible(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        self.try_divide(a, b).is_some()
    }
}

impl<S: CancellativeMultiplicationSignature> LeftCancellativeMultiplicationSignature for S {
    fn try_left_divide(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        self.try_divide(a, b)
    }
}

impl<S: CancellativeMultiplicationSignature> RightCancellativeMultiplicationSignature for S {
    fn try_right_divide(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        self.try_divide(a, b)
    }
}

#[signature_meta_trait]
pub trait AreAssociateMultiplicationSignature:
    CancellativeMultiplicationSignature + ZeroEqSignature
{
    fn are_associate(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        if self.equal(a, &self.zero()) && self.equal(b, &self.zero()) {
            true
        } else {
            self.try_divide(a, b).is_some() && self.try_divide(b, a).is_some()
        }
    }
}
impl<S: CancellativeMultiplicationSignature + ZeroEqSignature> AreAssociateMultiplicationSignature
    for S
{
}

#[signature_meta_trait]
pub trait MultiplicativeIntegralMonoidSignature:
    MultiplicativeAbsorptionMonoidSignature
    + TryLeftReciprocalSignature
    + TryRightReciprocalSignature
    + LeftCancellativeMultiplicationSignature
    + RightCancellativeMultiplicationSignature
{
}

// semi-rings need not have cancellative addition
// but reciprocals are unique when they exist because all reciprocals are two-sided due to commutativity
#[signature_meta_trait]
pub trait SemiRingSignature:
    AdditiveMonoidSignature
    + TryNegateSignature
    + MultiplicativeAbsorptionMonoidSignature
    + CommutativeMultiplicationSignature
    + LeftDistributiveMultiplicationOverAddition
    + RightDistributiveMultiplicationOverAddition
{
    fn from_nat(&self, x: impl Into<Natural>) -> Self::Elem {
        let x = x.into();
        if x == Natural::ZERO {
            self.zero()
        } else if x == Natural::ONE {
            self.one()
        } else {
            let two = self.add(&self.one(), &self.one());
            debug_assert!(x >= Natural::TWO);
            let bits: Vec<bool> = x.bits().collect();
            let mut ans = self.zero();
            let mut v = self.one();
            for b in bits {
                if b {
                    self.add_mut(&mut ans, &v);
                }
                self.mul_mut(&mut v, &two);
            }
            ans
        }
    }
}

#[signature_meta_trait]
pub trait SemiRingEqSignature: SemiRingSignature + EqSignature {}
impl<R: SemiRingSignature + EqSignature> SemiRingEqSignature for R {}

#[signature_meta_trait]
pub trait RingUnitsSignature: RingSignature + TryReciprocalSignature {}
impl<Ring: RingSignature + TryReciprocalSignature> RingUnitsSignature for Ring {}

#[signature_meta_trait]
pub trait RingSignature: SemiRingSignature + AdditiveGroupSignature {
    /// Determine whether the ring is reduced.
    ///
    /// Returns `Ok(true)` if the ring is reduced, `Ok(false)` if it is not,
    /// and `Err` when the implementation cannot decide.
    fn is_reduced(&self) -> Result<bool, String> {
        Err("unable to decide whether the ring is reduced".to_string())
    }

    fn bracket(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.sub(&self.mul(a, b), &self.mul(b, a))
    }

    fn from_int(&self, x: impl Into<Integer>) -> Self::Elem {
        let x = x.into();
        if x < Integer::ZERO {
            self.neg(&self.from_int(-x))
        } else {
            self.from_nat(x.abs())
        }
    }

    #[skip_meta]
    fn inbound_principal_integer_map(&self) -> PrincipalIntegerMap<Self, &Self> {
        PrincipalIntegerMap::new(self)
    }

    #[skip_meta]
    fn into_inbound_principal_integer_map(self) -> PrincipalIntegerMap<Self, Self> {
        PrincipalIntegerMap::new(self)
    }
}

#[signature_meta_trait]
pub trait RingEqSignature: RingSignature + EqSignature {}
impl<R: RingSignature + EqSignature> RingEqSignature for R {}

#[signature_meta_trait]
pub trait IntegralDomainSignature:
    RingSignature
    + TryReciprocalSignature
    + MultiplicativeIntegralMonoidSignature
    + CancellativeMultiplicationSignature
    + EqSignature
{
    fn try_from_rat(&self, x: &Rational) -> Option<Self::Elem> {
        let n = Fraction::numerator(x);
        let d = Fraction::denominator(x);
        debug_assert!(!d.is_zero());
        self.try_divide(&self.from_int(n), &self.from_nat(d))
    }
}

#[signature_meta_trait]
pub trait FavoriteAssociateSignature: TryReciprocalSignature + EqSignature {
    //For associate class of elements, choose a unique representative
    //write self=unit*assoc and return (unit, assoc)
    //0 is required to return (1, 0)
    //every unit u is required to return (u, 1) i.e. 1 is the favorite associate of every unit
    //it seems to happen that the product of favorite associates is another favorite associate. Should this be a requirement?

    fn factor_fav_assoc(&self, a: &Self::Elem) -> (Self::Elem, Self::Elem);
    fn fav_assoc(&self, a: &Self::Elem) -> Self::Elem {
        self.factor_fav_assoc(a).1
    }
    fn is_fav_assoc(&self, a: &Self::Elem) -> bool {
        let (_u, b) = self.factor_fav_assoc(a);
        self.equal(a, &b)
    }
}

#[signature_meta_trait]
pub trait CharacteristicSignature: SemiRingSignature {
    fn characteristic(&self) -> Natural;
}

#[signature_meta_trait]
pub trait OrderedRingSignature: IntegralDomainSignature + OrdSignature {
    fn abs(&self, a: &Self::Elem) -> Self::Elem {
        match self.cmp(a, &self.zero()) {
            std::cmp::Ordering::Less => self.neg(a),
            std::cmp::Ordering::Equal => self.zero(),
            std::cmp::Ordering::Greater => a.clone(),
        }
    }
}

#[signature_meta_trait]
pub trait FiniteUnitsSignature: RingSignature {
    fn all_units(&self) -> Vec<Self::Elem>;

    fn all_units_and_zero(&self) -> Vec<Self::Elem> {
        let mut elems = vec![self.zero()];
        elems.append(&mut self.all_units());
        elems
    }
}

impl<R: RingSignature + TryReciprocalSignature> FiniteUnitsSignature for R
where
    for<'a> MultiplicativeMonoidUnitsStructure<R, &'a R>: FiniteSetSignature<Elem = R::Elem>,
{
    fn all_units(&self) -> Vec<Self::Elem> {
        self.units().list_all_elements()
    }
}

#[signature_meta_trait]
pub trait GreatestCommonDivisorSignature:
    FavoriteAssociateSignature + IntegralDomainSignature
{
    //any gcds should be the standard associate representative
    //euclidean_gcd can be used to implement this
    fn gcd(&self, x: &Self::Elem, y: &Self::Elem) -> Self::Elem;
    fn gcd_list(&self, elems: Vec<impl Borrow<Self::Elem>>) -> Self::Elem {
        let mut gcd = self.zero();
        for x in elems {
            gcd = self.gcd(&gcd, x.borrow());
        }
        gcd
    }
    fn lcm(&self, x: &Self::Elem, y: &Self::Elem) -> Self::Elem {
        if self.is_zero(x) && self.is_zero(y) {
            self.zero()
        } else {
            self.try_divide(&self.mul(x, y), &self.gcd(x, y)).unwrap()
        }
    }
    fn lcm_list(&self, elems: Vec<impl Borrow<Self::Elem>>) -> Self::Elem {
        let mut lcm = self.one();
        for x in elems {
            lcm = self.lcm(&lcm, x.borrow());
        }
        lcm
    }
}

#[signature_meta_trait]
pub trait BezoutDomainSignature: GreatestCommonDivisorSignature {
    //any gcds should be the standard associate representative
    fn xgcd(&self, a: &Self::Elem, b: &Self::Elem) -> (Self::Elem, Self::Elem, Self::Elem); //(g, x, y) s.t. g = ax + by
    fn xgcd_list(&self, elems: Vec<&Self::Elem>) -> (Self::Elem, Vec<Self::Elem>) {
        // println!("{:?}", elems);
        match elems.len() {
            0 => (self.zero(), vec![]),
            1 => {
                let (unit, assoc) = self.factor_fav_assoc(elems[0]);
                (assoc, vec![self.try_reciprocal(&unit).unwrap()])
            }
            2 => {
                let (g, x, y) = self.xgcd(elems[0], elems[1]);
                (g, vec![x, y])
            }
            n => {
                let k = n / 2;
                let (g1, coeffs1) = self.xgcd_list((0..k).map(|i| elems[i]).collect());
                let (g2, coeffs2) = self.xgcd_list((k..n).map(|i| elems[i]).collect());
                let (g, x, y) = self.xgcd(&g1, &g2);
                let mut coeffs = vec![];
                for c in coeffs1 {
                    coeffs.push(self.mul(&x, &c));
                }
                for c in coeffs2 {
                    coeffs.push(self.mul(&y, &c));
                }
                (g, coeffs)
            }
        }
    }
}

#[signature_meta_trait]
pub trait EuclideanDivisionSignature: SemiRingEqSignature {
    /// None for 0 and Some(norm) for everything else
    fn norm(&self, elem: &Self::Elem) -> Option<Natural>;

    /// None if b is 0 and Some((q, r)) such that a=bq+r and r=0 or norm(r) < norm(b)
    fn quorem(&self, a: &Self::Elem, b: &Self::Elem) -> Option<(Self::Elem, Self::Elem)>;

    fn quo(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        self.quorem(a, b).map(|(q, _r)| q)
    }

    fn rem(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        if self.is_zero(b) {
            a.clone()
        } else {
            let (_q, r) = self.quorem(a, b).unwrap();
            r
        }
    }

    fn euclidean_gcd(&self, mut x: Self::Elem, mut y: Self::Elem) -> Self::Elem
    where
        Self: FavoriteAssociateSignature,
    {
        //Euclidean algorithm
        while !self.is_zero(&y) {
            let r = self.rem(&x, &y);
            (x, y) = (y, r);
        }
        let (_unit, assoc) = self.factor_fav_assoc(&x);
        assoc
    }

    fn euclidean_xgcd(
        &self,
        mut x: Self::Elem,
        mut y: Self::Elem,
    ) -> (Self::Elem, Self::Elem, Self::Elem)
    where
        Self: FavoriteAssociateSignature + IntegralDomainSignature,
    {
        let orig_x = x.clone();
        let orig_y = y.clone();

        let mut pa = self.one();
        let mut a = self.zero();
        let mut pb = self.zero();
        let mut b = self.one();

        while !self.is_zero(&y) {
            let (q, r) = self.quorem(&x, &y).unwrap();
            let new_a = self.add(&pa, &self.neg(&self.mul(&q, &a)));
            (a, pa) = (new_a, a);
            let new_b = self.add(&pb, &self.neg(&self.mul(&q, &b)));
            (b, pb) = (new_b, b);
            (x, y) = (y, r);
        }
        let (unit, ass_x) = self.factor_fav_assoc(&x);
        // g = u*g_ass
        // g = xa+by
        // xa+by=u*g_ass
        debug_assert!(self.is_unit(&unit));
        let (g, a, b) = (
            ass_x,
            self.try_divide(&pa, &unit).unwrap(),
            self.try_divide(&pb, &unit).unwrap(),
        );
        // println!("{:?} = {:?} * {:?} + {:?} * {:?}", g, a, orig_x, b, orig_y);
        debug_assert!(self.equal(
            &self.add(&self.mul(&a, &orig_x), &self.mul(&b, &orig_y)),
            &g
        ));
        debug_assert!(self.equal(&g, &self.euclidean_gcd(orig_x, orig_y)));
        (g, a, b)
    }
}

#[signature_meta_trait]
pub trait EuclideanDomainSignature: EuclideanDivisionSignature + IntegralDomainSignature {}
impl<Ring: EuclideanDivisionSignature + IntegralDomainSignature> EuclideanDomainSignature for Ring {}

#[signature_meta_trait]
pub trait InfiniteSignature: RinglikeSpecializationSignature {
    fn generate_distinct_elements(&self) -> Box<dyn Iterator<Item = Self::Elem>>;
}

#[signature_meta_trait]
pub trait FieldSignature: IntegralDomainSignature {
    fn from_rat(&self, x: &Rational) -> Self::Elem {
        self.try_from_rat(x).unwrap()
    }
}

impl<FS: FieldSignature> FavoriteAssociateSignature for FS {
    fn factor_fav_assoc(&self, a: &Self::Elem) -> (Self::Elem, Self::Elem) {
        if self.is_zero(a) {
            (self.one(), self.zero())
        } else {
            (a.clone(), self.one())
        }
    }
}

impl<FS: FieldSignature> EuclideanDivisionSignature for FS {
    fn norm(&self, elem: &Self::Elem) -> Option<Natural> {
        if self.is_zero(elem) {
            None
        } else {
            Some(Natural::from(1u8))
        }
    }

    fn quorem(&self, a: &Self::Elem, b: &Self::Elem) -> Option<(Self::Elem, Self::Elem)> {
        if self.is_zero(b) {
            None
        } else {
            Some((self.try_divide(a, b).unwrap(), self.zero()))
        }
    }
}

impl<FS: FieldSignature> GreatestCommonDivisorSignature for FS {
    fn gcd(&self, x: &Self::Elem, y: &Self::Elem) -> Self::Elem {
        self.euclidean_gcd(x.clone(), y.clone())
    }
}

impl<FS: FieldSignature> BezoutDomainSignature for FS {
    fn xgcd(&self, x: &Self::Elem, y: &Self::Elem) -> (Self::Elem, Self::Elem, Self::Elem) {
        self.euclidean_xgcd(x.clone(), y.clone())
    }
}

/// When `.characteristic()` always returns 0
#[signature_meta_trait]
pub trait CharZeroRingSignature: RingSignature + CharacteristicSignature {
    fn try_to_int(&self, x: &Self::Elem) -> Option<Integer>;
}

impl<RS: CharZeroRingSignature + 'static> InfiniteSignature for RS {
    fn generate_distinct_elements(&self) -> Box<dyn Iterator<Item = <Self as SetSignature>::Elem>> {
        struct IntegerIterator<RS: CharZeroRingSignature> {
            ring: RS,
            next: Integer,
        }

        impl<RS: CharZeroRingSignature> Iterator for IntegerIterator<RS> {
            type Item = RS::Elem;

            fn next(&mut self) -> Option<Self::Item> {
                let next = self.next.clone();
                if Integer::ZERO < next {
                    self.next = -self.next.clone();
                } else {
                    self.next = Integer::from(1) - self.next.clone();
                }
                Some(self.ring.from_int(next))
            }
        }

        Box::new(IntegerIterator {
            ring: self.clone(),
            next: Integer::from(0),
        })
    }
}

#[signature_meta_trait]
pub trait CharZeroFieldSignature: FieldSignature + CharZeroRingSignature {
    fn try_to_rat(&self, x: &Self::Elem) -> Option<Rational>;

    #[skip_meta]
    fn inbound_principal_rational_map(&self) -> PrincipalRationalMap<Self, &Self> {
        PrincipalRationalMap::new(self)
    }
    #[skip_meta]
    fn into_inbound_principal_rational_map(self) -> PrincipalRationalMap<Self, Self> {
        PrincipalRationalMap::new(self)
    }
}

#[signature_meta_trait]
pub trait FiniteFieldSignature: FieldSignature + FiniteUnitsSignature + FiniteSetSignature {
    // Return (p, k) where p is a prime and |F| = p^k
    fn characteristic_and_power(&self) -> (Natural, Natural);
}

//is a subset of the complex numbers
#[signature_meta_trait]
pub trait ComplexSubsetSignature: RinglikeSpecializationSignature {
    fn as_f32_real_and_imaginary_parts(&self, z: &Self::Elem) -> (f32, f32);
    fn as_f64_real_and_imaginary_parts(&self, z: &Self::Elem) -> (f64, f64);
}

//is a subset of the real numbers
#[signature_meta_trait]
pub trait RealSubsetSignature: ComplexSubsetSignature {
    fn as_f64(&self, x: &Self::Elem) -> f64 {
        let (r, i) = self.as_f64_real_and_imaginary_parts(x);
        debug_assert_eq!(i, 0.0);
        r
    }

    fn as_f32(&self, x: &Self::Elem) -> f32 {
        let (r, i) = self.as_f32_real_and_imaginary_parts(x);
        debug_assert_eq!(i, 0.0);
        r
    }
}

#[signature_meta_trait]
pub trait RealRoundingSignature: RealSubsetSignature {
    fn floor(&self, x: &Self::Elem) -> Integer; //round down
    fn ceil(&self, x: &Self::Elem) -> Integer; //round up
    fn round(&self, x: &Self::Elem) -> Integer; //round closets, either direction is fine if mid way
}

#[signature_meta_trait]
pub trait RealFromFloatSignature: RealSubsetSignature {
    fn from_f64_approx(&self, x: f64) -> Self::Elem;
    fn from_f32_approx(&self, x: f32) -> Self::Elem {
        self.from_f64_approx(f64::from(x))
    }
}

#[signature_meta_trait]
pub trait ComplexConjugateSignature: ComplexSubsetSignature {
    fn conjugate(&self, x: &Self::Elem) -> Self::Elem;
}
impl<RS: RealSubsetSignature> ComplexConjugateSignature for RS {
    fn conjugate(&self, x: &Self::Elem) -> Self::Elem {
        x.clone()
    }
}

#[signature_meta_trait]
pub trait PositiveRealNthRootSignature: ComplexSubsetSignature {
    //if x is a non-negative real number, return the nth root of x
    //may also return Ok for other well-defined values such as for 1st root of any x and 0th root of any non-zero x, but is not required to
    fn nth_root(&self, x: &Self::Elem, n: usize) -> Result<Self::Elem, ()>;
    fn square_root(&self, x: &Self::Elem) -> Result<Self::Elem, ()> {
        self.nth_root(x, 2)
    }
    fn cube_root(&self, x: &Self::Elem) -> Result<Self::Elem, ()> {
        self.nth_root(x, 3)
    }
}

// TODO: Move this sort of struture to the field inclusion homomorphism
#[signature_meta_trait]
pub trait AlgebraicClosureSignature: FieldSignature
where
    //TODO: can this allow polynomial structures taking a reference to the base field rather than an instance?
    PolynomialStructure<Self::BFS, Self::BFS>: FactoringMonoidSignature<FactoredExponent = NaturalCanonicalStructure>
        + SetSignature<Elem = Polynomial<<Self::BFS as SetSignature>::Elem>>,
{
    type BFS: FieldSignature; //base field structure

    fn base_field(&self) -> Self::BFS;

    fn base_field_inclusion(&self, x: &<Self::BFS as SetSignature>::Elem) -> Self::Elem;

    //return None for the zero polynomial
    fn all_roots_list(
        &self,
        poly: &Polynomial<<Self::BFS as SetSignature>::Elem>,
    ) -> Option<Vec<Self::Elem>>;

    fn all_roots_unique(
        &self,
        poly: &Polynomial<<Self::BFS as SetSignature>::Elem>,
    ) -> Option<Vec<Self::Elem>> {
        let base_field_poly = self.base_field().into_polynomials();
        self.all_roots_list(
            &base_field_poly
                .factorizations()
                .expand_squarefree(&base_field_poly.factor(poly)),
        )
    }

    fn all_roots_powers(
        &self,
        poly: &Polynomial<<Self::BFS as SetSignature>::Elem>,
    ) -> Option<Vec<(Self::Elem, usize)>> {
        let mut root_powers = vec![];
        let base_field_poly = self.base_field().into_polynomials();
        for (factor, k) in base_field_poly.factor(poly).into_powers()? {
            for root in self.all_roots_list(&factor).unwrap() {
                root_powers.push((root, (&k).try_into().unwrap()));
            }
        }
        Some(root_powers)
    }
}

/// The free ring of rank 0 is the integers
/// The free ring of rank 1 is the polynomial ring over the integers
/// The free ring of rank n is the multipolynomial ring over the integers
#[signature_meta_trait]
pub trait FreeRingSignature: RingSignature {
    type Generator: Clone + Debug + PartialEq + Eq + std::hash::Hash + Send + Sync;

    fn free_generators(&self) -> std::collections::HashSet<Self::Generator>;
    fn free_rank(&self) -> usize {
        self.free_generators().len()
    }
}
