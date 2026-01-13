use super::*;
use crate::polynomial::*;
use algebraeon_macros::{signature_meta_trait, skip_meta};
use algebraeon_nzq::{Integer, Natural, NaturalCanonicalStructure, Rational, traits::*};
use algebraeon_sets::structure::*;
use std::{borrow::Borrow, fmt::Debug};

mod unconstructable_universal_structure {
    use algebraeon_sets::structure::{EqSignature, SetSignature, Signature};
    use std::fmt::Debug;
    use std::marker::PhantomData;

    use crate::structure::{
        AdditiveGroupSignature, AdditiveMonoidSignature, CancellativeAdditiveMonoidSignature,
        CharZeroRingSignature, CharacteristicSignature, MultiplicativeMonoidSignature,
        RingSignature, RinglikeSpecializationSignature, SemiRingSignature, SetWithZeroSignature,
    };

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
        type Set = Set;

        fn is_element(&self, _x: &Self::Set) -> Result<(), String> {
            unreachable!()
        }
    }

    impl<Set: Debug + Clone + Send + Sync> EqSignature for UnconstructableStructure<Set> {
        fn equal(&self, _a: &Self::Set, _b: &Self::Set) -> bool {
            unreachable!()
        }
    }

    impl<Set: Debug + Clone + Send + Sync> RinglikeSpecializationSignature
        for UnconstructableStructure<Set>
    {
    }

    impl<Set: Debug + Clone + Send + Sync> SetWithZeroSignature for UnconstructableStructure<Set> {
        fn zero(&self) -> Self::Set {
            unreachable!()
        }
    }

    impl<Set: Debug + Clone + Send + Sync> AdditiveMonoidSignature for UnconstructableStructure<Set> {
        fn add(&self, _a: &Self::Set, _b: &Self::Set) -> Self::Set {
            unreachable!()
        }

        fn try_neg(&self, _a: &Self::Set) -> Option<Self::Set> {
            unreachable!()
        }
    }

    impl<Set: Debug + Clone + Send + Sync> CancellativeAdditiveMonoidSignature
        for UnconstructableStructure<Set>
    {
        fn try_sub(&self, _a: &Self::Set, _b: &Self::Set) -> Option<Self::Set> {
            unreachable!()
        }
    }

    impl<Set: Debug + Clone + Send + Sync> AdditiveGroupSignature for UnconstructableStructure<Set> {
        fn neg(&self, _a: &Self::Set) -> Self::Set {
            unreachable!()
        }
    }

    impl<Set: Debug + Clone + Send + Sync> MultiplicativeMonoidSignature
        for UnconstructableStructure<Set>
    {
        fn one(&self) -> Self::Set {
            unreachable!()
        }

        fn mul(&self, _a: &Self::Set, _b: &Self::Set) -> Self::Set {
            unreachable!()
        }
    }

    impl<Set: Debug + Clone + Send + Sync> SemiRingSignature for UnconstructableStructure<Set> {}

    impl<Set: Debug + Clone + Send + Sync> CharacteristicSignature for UnconstructableStructure<Set> {
        fn characteristic(&self) -> algebraeon_nzq::Natural {
            unreachable!()
        }
    }

    impl<Set: Debug + Clone + Send + Sync> RingSignature for UnconstructableStructure<Set> {}

    impl<Set: Debug + Clone + Send + Sync> CharZeroRingSignature for UnconstructableStructure<Set> {
        fn try_to_int(&self, _x: &Self::Set) -> Option<algebraeon_nzq::Integer> {
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
    fn try_ring_restructure(&self) -> Option<impl EqSignature<Set = Self::Set> + RingSignature> {
        Option::<unconstructable_universal_structure::UnconstructableStructure<Self::Set>>::None
    }

    /*
    Used by:
     - Formatting polynomials as strings: If the set of coefficients has this structure then it's possible to call .try_to_int(..) which can allow for nicer formatting at integer coefficients.
     */
    fn try_char_zero_ring_restructure(
        &self,
    ) -> Option<impl EqSignature<Set = Self::Set> + CharZeroRingSignature> {
        Option::<unconstructable_universal_structure::UnconstructableStructure<Self::Set>>::None
    }
}

#[signature_meta_trait]
pub trait SetWithZeroSignature: RinglikeSpecializationSignature {
    fn zero(&self) -> Self::Set;
}

#[signature_meta_trait]
pub trait AdditiveMonoidSignature: SetWithZeroSignature {
    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set;

    fn add_mut(&self, a: &mut Self::Set, b: &Self::Set) {
        *a = self.add(a, b);
    }

    fn try_neg(&self, a: &Self::Set) -> Option<Self::Set>;

    fn sum(&self, vals: Vec<impl Borrow<Self::Set>>) -> Self::Set {
        let mut sum = self.zero();
        for val in vals {
            self.add_mut(&mut sum, val.borrow());
        }
        sum
    }
}

#[signature_meta_trait]
pub trait CancellativeAdditiveMonoidSignature: AdditiveMonoidSignature {
    fn try_sub(&self, a: &Self::Set, b: &Self::Set) -> Option<Self::Set>;
}

#[signature_meta_trait]
pub trait SetWithZeroAndEqSignature: SetWithZeroSignature + EqSignature {
    fn is_zero(&self, a: &Self::Set) -> bool {
        self.equal(a, &self.zero())
    }
}
impl<R: SetWithZeroSignature + EqSignature> SetWithZeroAndEqSignature for R {}

#[signature_meta_trait]
pub trait MultiplicativeMonoidSignature: RinglikeSpecializationSignature {
    fn one(&self) -> Self::Set;

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set;

    fn mul_mut(&self, a: &mut Self::Set, b: &Self::Set) {
        *a = self.mul(a, b);
    }

    fn product(&self, vals: Vec<impl Borrow<Self::Set>>) -> Self::Set {
        let mut prod = self.one();
        for val in vals {
            self.mul_mut(&mut prod, val.borrow());
        }
        prod
    }

    fn nat_pow(&self, a: &Self::Set, n: &Natural) -> Self::Set {
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
pub trait MultiplicativeMonoidSquareOpsSignature: RinglikeSpecializationSignature {
    fn is_square(&self, a: &Self::Set) -> bool {
        self.sqrt_if_square(a).is_some()
    }

    fn sqrt_if_square(&self, a: &Self::Set) -> Option<Self::Set>;
}

/// 0 is such that 0*a=0 for all a in the monoid.
/// such an element is unqiue if it exists.
#[signature_meta_trait]
pub trait MultiplicativeMonoidWithZeroSignature:
    MultiplicativeMonoidSignature + SetWithZeroSignature
{
}
impl<R: MultiplicativeMonoidSignature + SetWithZeroSignature> MultiplicativeMonoidWithZeroSignature
    for R
{
}

#[signature_meta_trait]
pub trait MultiplicativeMonoidUnitsSignature: MultiplicativeMonoidSignature {
    /// b such that a*b=1 and b*a=1
    /// Err(DivideByZero) if b is zero
    /// Err(NotDivisible) if no such b exists
    fn try_inv(&self, a: &Self::Set) -> Option<Self::Set>;

    fn is_unit(&self, a: &Self::Set) -> bool {
        self.try_inv(a).is_some()
    }

    fn try_int_pow(&self, a: &Self::Set, n: &Integer) -> Option<Self::Set> {
        if *n == Integer::ZERO {
            Some(self.one())
        } else if *n > Integer::ZERO {
            Some(self.nat_pow(a, &Abs::abs(n)))
        } else {
            Some(self.nat_pow(&self.try_inv(a)?, &Abs::abs(n)))
        }
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

#[signature_meta_trait]
pub trait MultiplicativeIntegralMonoidSignature:
    MultiplicativeMonoidUnitsSignature + SetWithZeroSignature + EqSignature
{
    fn try_div(&self, a: &Self::Set, b: &Self::Set) -> Option<Self::Set>;

    /// return true iff a is divisible by b
    fn divisible(&self, a: &Self::Set, b: &Self::Set) -> bool {
        self.try_div(a, b).is_some()
    }
    fn are_associate(&self, a: &Self::Set, b: &Self::Set) -> bool {
        if self.equal(a, &self.zero()) && self.equal(b, &self.zero()) {
            true
        } else {
            self.try_div(a, b).is_some() && self.try_div(b, a).is_some()
        }
    }
}

#[signature_meta_trait]
pub trait MultiplicativeGroupSignature: MultiplicativeMonoidSignature {
    fn inv(&self, a: &Self::Set) -> Self::Set;

    fn int_pow(&self, a: &Self::Set, n: &Integer) -> Self::Set {
        if *n == Integer::ZERO {
            self.one()
        } else if *n > Integer::ZERO {
            self.nat_pow(a, &Abs::abs(n))
        } else {
            self.nat_pow(&self.inv(a), &Abs::abs(n))
        }
    }
}

#[signature_meta_trait]
pub trait FavoriteAssociateSignature: MultiplicativeMonoidUnitsSignature + EqSignature {
    //For associate class of elements, choose a unique representative
    //write self=unit*assoc and return (unit, assoc)
    //0 is required to return (1, 0)
    //every unit u is required to return (u, 1) i.e. 1 is the favorite associate of every unit
    //it seems to happen that the product of favorite associates is another favorite associate. Should this be a requirement?

    fn factor_fav_assoc(&self, a: &Self::Set) -> (Self::Set, Self::Set);
    fn fav_assoc(&self, a: &Self::Set) -> Self::Set {
        self.factor_fav_assoc(a).1
    }
    fn is_fav_assoc(&self, a: &Self::Set) -> bool {
        let (_u, b) = self.factor_fav_assoc(a);
        self.equal(a, &b)
    }
}

#[signature_meta_trait]
pub trait SemiRingSignature:
    AdditiveMonoidSignature + MultiplicativeMonoidWithZeroSignature
{
    fn from_nat(&self, x: impl Into<Natural>) -> Self::Set {
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
pub trait CharacteristicSignature: SemiRingSignature {
    fn characteristic(&self) -> Natural;
}

#[signature_meta_trait]
pub trait RingUnitsSignature: RingSignature + MultiplicativeMonoidUnitsSignature {}
impl<Ring: RingSignature + MultiplicativeMonoidUnitsSignature> RingUnitsSignature for Ring {}

#[signature_meta_trait]
pub trait AdditiveGroupSignature: CancellativeAdditiveMonoidSignature {
    fn neg(&self, a: &Self::Set) -> Self::Set;

    fn sub(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.add(a, &self.neg(b))
    }
}

#[signature_meta_trait]
pub trait RingSignature: SemiRingSignature + AdditiveGroupSignature {
    /// Determine whether the ring is reduced.
    ///
    /// Returns `Ok(true)` if the ring is reduced, `Ok(false)` if it is not,
    /// and `Err` when the implementation cannot decide.
    fn is_reduced(&self) -> Result<bool, String> {
        Err("unable to decide whether the ring is reduced".to_string())
    }

    fn bracket(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.sub(&self.mul(a, b), &self.mul(b, a))
    }

    fn from_int(&self, x: impl Into<Integer>) -> Self::Set {
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
    RingUnitsSignature + MultiplicativeIntegralMonoidSignature + EqSignature
{
    fn try_from_rat(&self, x: &Rational) -> Option<Self::Set> {
        let n = Fraction::numerator(x);
        let d = Fraction::denominator(x);
        debug_assert!(!d.is_zero());
        self.try_div(&self.from_int(n), &self.from_nat(d))
    }
}

#[signature_meta_trait]
pub trait OrderedRingSignature: IntegralDomainSignature {
    // <= satisfying translation invariance and multiplication by positive scalar
    fn ring_cmp(&self, a: &Self::Set, b: &Self::Set) -> std::cmp::Ordering;
    fn abs(&self, a: &Self::Set) -> Self::Set {
        match self.ring_cmp(a, &self.zero()) {
            std::cmp::Ordering::Less => self.neg(a),
            std::cmp::Ordering::Equal => self.zero(),
            std::cmp::Ordering::Greater => a.clone(),
        }
    }
}

#[signature_meta_trait]
pub trait FiniteUnitsSignature: RingSignature {
    fn all_units(&self) -> Vec<Self::Set>;

    fn all_units_and_zero(&self) -> Vec<Self::Set> {
        let mut elems = vec![self.zero()];
        elems.append(&mut self.all_units());
        elems
    }
}

impl<R: RingSignature + MultiplicativeMonoidUnitsSignature> FiniteUnitsSignature for R
where
    for<'a> MultiplicativeMonoidUnitsStructure<R, &'a R>: FiniteSetSignature<Set = R::Set>,
{
    fn all_units(&self) -> Vec<Self::Set> {
        self.units().list_all_elements()
    }
}

#[signature_meta_trait]
pub trait GreatestCommonDivisorSignature:
    FavoriteAssociateSignature + IntegralDomainSignature
{
    //any gcds should be the standard associate representative
    //euclidean_gcd can be used to implement this
    fn gcd(&self, x: &Self::Set, y: &Self::Set) -> Self::Set;
    fn gcd_list(&self, elems: Vec<impl Borrow<Self::Set>>) -> Self::Set {
        let mut gcd = self.zero();
        for x in elems {
            gcd = self.gcd(&gcd, x.borrow());
        }
        gcd
    }
    fn lcm(&self, x: &Self::Set, y: &Self::Set) -> Self::Set {
        if self.is_zero(x) && self.is_zero(y) {
            self.zero()
        } else {
            self.try_div(&self.mul(x, y), &self.gcd(x, y)).unwrap()
        }
    }
    fn lcm_list(&self, elems: Vec<impl Borrow<Self::Set>>) -> Self::Set {
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
    fn xgcd(&self, a: &Self::Set, b: &Self::Set) -> (Self::Set, Self::Set, Self::Set); //(g, x, y) s.t. g = ax + by
    fn xgcd_list(&self, elems: Vec<&Self::Set>) -> (Self::Set, Vec<Self::Set>) {
        // println!("{:?}", elems);
        match elems.len() {
            0 => (self.zero(), vec![]),
            1 => {
                let (unit, assoc) = self.factor_fav_assoc(elems[0]);
                (assoc, vec![self.try_inv(&unit).unwrap()])
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
    fn norm(&self, elem: &Self::Set) -> Option<Natural>;

    /// None if b is 0 and Some((q, r)) such that a=bq+r and r=0 or norm(r) < norm(b)
    fn quorem(&self, a: &Self::Set, b: &Self::Set) -> Option<(Self::Set, Self::Set)>;

    fn quo(&self, a: &Self::Set, b: &Self::Set) -> Option<Self::Set> {
        self.quorem(a, b).map(|(q, _r)| q)
    }

    fn rem(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        if self.is_zero(b) {
            a.clone()
        } else {
            let (_q, r) = self.quorem(a, b).unwrap();
            r
        }
    }

    fn euclidean_gcd(&self, mut x: Self::Set, mut y: Self::Set) -> Self::Set
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
        mut x: Self::Set,
        mut y: Self::Set,
    ) -> (Self::Set, Self::Set, Self::Set)
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
            self.try_div(&pa, &unit).unwrap(),
            self.try_div(&pb, &unit).unwrap(),
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
    fn generate_distinct_elements(&self) -> Box<dyn Iterator<Item = Self::Set>>;
}

#[signature_meta_trait]
pub trait FieldSignature: IntegralDomainSignature {
    fn from_rat(&self, x: &Rational) -> Self::Set {
        self.try_from_rat(x).unwrap()
    }
}

impl<FS: FieldSignature> FavoriteAssociateSignature for FS {
    fn factor_fav_assoc(&self, a: &Self::Set) -> (Self::Set, Self::Set) {
        if self.is_zero(a) {
            (self.one(), self.zero())
        } else {
            (a.clone(), self.one())
        }
    }
}

impl<FS: FieldSignature> EuclideanDivisionSignature for FS {
    fn norm(&self, elem: &Self::Set) -> Option<Natural> {
        if self.is_zero(elem) {
            None
        } else {
            Some(Natural::from(1u8))
        }
    }

    fn quorem(&self, a: &Self::Set, b: &Self::Set) -> Option<(Self::Set, Self::Set)> {
        if self.is_zero(b) {
            None
        } else {
            Some((self.try_div(a, b).unwrap(), self.zero()))
        }
    }
}

impl<FS: FieldSignature> GreatestCommonDivisorSignature for FS {
    fn gcd(&self, x: &Self::Set, y: &Self::Set) -> Self::Set {
        self.euclidean_gcd(x.clone(), y.clone())
    }
}

impl<FS: FieldSignature> BezoutDomainSignature for FS {
    fn xgcd(&self, x: &Self::Set, y: &Self::Set) -> (Self::Set, Self::Set, Self::Set) {
        self.euclidean_xgcd(x.clone(), y.clone())
    }
}

/// When `.characteristic()` always returns 0
#[signature_meta_trait]
pub trait CharZeroRingSignature: RingSignature + CharacteristicSignature {
    fn try_to_int(&self, x: &Self::Set) -> Option<Integer>;
}

impl<RS: CharZeroRingSignature + 'static> InfiniteSignature for RS {
    fn generate_distinct_elements(&self) -> Box<dyn Iterator<Item = <Self as SetSignature>::Set>> {
        struct IntegerIterator<RS: CharZeroRingSignature> {
            ring: RS,
            next: Integer,
        }

        impl<RS: CharZeroRingSignature> Iterator for IntegerIterator<RS> {
            type Item = RS::Set;

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
    fn try_to_rat(&self, x: &Self::Set) -> Option<Rational>;

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
    fn as_f32_real_and_imaginary_parts(&self, z: &Self::Set) -> (f32, f32);
    fn as_f64_real_and_imaginary_parts(&self, z: &Self::Set) -> (f64, f64);
}

//is a subset of the real numbers
#[signature_meta_trait]
pub trait RealSubsetSignature: ComplexSubsetSignature {
    fn as_f64(&self, x: &Self::Set) -> f64 {
        let (r, i) = self.as_f64_real_and_imaginary_parts(x);
        debug_assert_eq!(i, 0.0);
        r
    }
    fn as_f32(&self, x: &Self::Set) -> f32 {
        let (r, i) = self.as_f32_real_and_imaginary_parts(x);
        debug_assert_eq!(i, 0.0);
        r
    }
}

#[signature_meta_trait]
pub trait RealRoundingSignature: RealSubsetSignature {
    fn floor(&self, x: &Self::Set) -> Integer; //round down
    fn ceil(&self, x: &Self::Set) -> Integer; //round up
    fn round(&self, x: &Self::Set) -> Integer; //round closets, either direction is fine if mid way
}

#[signature_meta_trait]
pub trait RealFromFloatSignature: RealSubsetSignature {
    fn from_f64_approx(&self, x: f64) -> Self::Set;
    fn from_f32_approx(&self, x: f32) -> Self::Set {
        self.from_f64_approx(f64::from(x))
    }
}

#[signature_meta_trait]
pub trait ComplexConjugateSignature: RinglikeSpecializationSignature {
    fn conjugate(&self, x: &Self::Set) -> Self::Set;
}

impl<RS: RealSubsetSignature> ComplexConjugateSignature for RS {
    fn conjugate(&self, x: &Self::Set) -> Self::Set {
        x.clone()
    }
}

#[signature_meta_trait]
pub trait PositiveRealNthRootSignature: ComplexSubsetSignature {
    //if x is a non-negative real number, return the nth root of x
    //may also return Ok for other well-defined values such as for 1st root of any x and 0th root of any non-zero x, but is not required to
    fn nth_root(&self, x: &Self::Set, n: usize) -> Result<Self::Set, ()>;
    fn square_root(&self, x: &Self::Set) -> Result<Self::Set, ()> {
        self.nth_root(x, 2)
    }
    fn cube_root(&self, x: &Self::Set) -> Result<Self::Set, ()> {
        self.nth_root(x, 3)
    }
}

// TODO: Move this sort of struture to the field inclusion homomorphism
// #[signature_meta_trait]
pub trait AlgebraicClosureSignature: FieldSignature
where
    //TODO: can this allow polynomial structures taking a reference to the base field rather than an instance?
    PolynomialStructure<Self::BFS, Self::BFS>: FactoringMonoidSignature<FactoredExponent = NaturalCanonicalStructure>
        + SetSignature<Set = Polynomial<<Self::BFS as SetSignature>::Set>>,
{
    type BFS: FieldSignature; //base field structure

    fn base_field(&self) -> Self::BFS;

    fn base_field_inclusion(&self, x: &<Self::BFS as SetSignature>::Set) -> Self::Set;

    //return None for the zero polynomial
    fn all_roots_list(
        &self,
        poly: &Polynomial<<Self::BFS as SetSignature>::Set>,
    ) -> Option<Vec<Self::Set>>;

    fn all_roots_unique(
        &self,
        poly: &Polynomial<<Self::BFS as SetSignature>::Set>,
    ) -> Option<Vec<Self::Set>> {
        let base_field_poly = self.base_field().into_polynomials();
        self.all_roots_list(
            &base_field_poly
                .factorizations()
                .expand_squarefree(&base_field_poly.factor(poly)),
        )
    }

    fn all_roots_powers(
        &self,
        poly: &Polynomial<<Self::BFS as SetSignature>::Set>,
    ) -> Option<Vec<(Self::Set, usize)>> {
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
pub trait FreeRingSignature: RingSignature {
    type Generator: Clone + Debug + PartialEq + Eq + std::hash::Hash + Send + Sync;

    fn free_generators(&self) -> std::collections::HashSet<Self::Generator>;
    fn free_rank(&self) -> usize {
        self.free_generators().len()
    }
}
