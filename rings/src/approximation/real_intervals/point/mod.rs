use crate::{
    approximation::{
        rational_interval::RationalInterval,
        real_intervals::{RationalPoint, Subset, SubsetsStructure},
    },
    continued_fraction::{SimpleContinuedFraction, ToSimpleContinuedFractionSignature},
    structure::{
        AdditionSignature, AdditiveGroupSignature, AdditiveMonoidSignature,
        CancellativeAdditionSignature, CommutativeMultiplicationSignature, ComplexSubsetSignature,
        LeftDistributiveMultiplicationOverAddition, MetaRealRoundingSignature,
        MetaRealSubsetSignature, MetaTryReciprocalSignature, MultiplicationSignature,
        MultiplicativeAbsorptionMonoidSignature, MultiplicativeMonoidSignature, OneSignature,
        RealRoundingSignature, RealSubsetSignature, RightDistributiveMultiplicationOverAddition,
        RingSignature, RinglikeSpecializationSignature, SemiRingSignature, TryNegateSignature,
        TryReciprocalSignature, ZeroSignature,
    },
};
use algebraeon_macros::CanonicalStructure;
use algebraeon_sets::approximations::ApproximatePointsSignature;
use algebraeon_structures::*;
use std::{
    fmt::Debug,
    sync::{Arc, Mutex, MutexGuard},
};

pub trait RealApproximatePointInterface: Debug + Send + Sync {
    fn rational_interval_neighbourhood(&self) -> Subset;
    fn length(&self) -> Rational {
        self.rational_interval_neighbourhood().length()
    }
    fn refine(&mut self);
    fn refine_to_length(&mut self, length: &Rational) {
        while self.length() > *length {
            self.refine();
        }
    }
}

#[derive(Debug, Clone, CanonicalStructure)]
pub struct RealApproximatePoint {
    repr: Arc<Mutex<dyn RealApproximatePointInterface>>,
}

impl RealApproximatePoint {
    pub fn new<R: RealApproximatePointInterface + 'static>(repr: R) -> Self {
        Self {
            repr: Arc::new(Mutex::new(repr)),
        }
    }

    pub fn from_rat(rat: Rational) -> Self {
        Self::new(rational::RationalPoint { x: rat })
    }

    pub fn from_continued_fraction(cf: impl SimpleContinuedFraction + 'static) -> Self {
        Self::new(continued_fraction::SimpleContinuedFractionPoint::from(cf))
    }

    pub fn lock(&self) -> MutexGuard<'_, dyn RealApproximatePointInterface + 'static> {
        self.repr.lock().unwrap()
    }
}

impl ApproximatePointsSignature for RealApproximatePointCanonicalStructure {
    type Precision = RationalCanonicalStructure;
    type OpenSubsetsStructure = SubsetsStructure;

    fn open_neighbourhood(&self, approx_point: &Self::Elem) -> Subset {
        approx_point.lock().rational_interval_neighbourhood()
    }

    fn precision(&self, approx_point: &Self::Elem) -> Rational {
        approx_point.lock().length()
    }

    fn refine(&self, approx_point: &mut Self::Elem) {
        approx_point.lock().refine();
    }

    fn refine_to(&self, approx_point: &mut Self::Elem, length: &Rational) {
        approx_point.lock().refine_to_length(length);
    }
}

#[derive(Debug)]
struct AddPoints {
    first: RealApproximatePoint,
    second: RealApproximatePoint,
}
impl RealApproximatePointInterface for AddPoints {
    fn rational_interval_neighbourhood(&self) -> Subset {
        let first_nbd = self.first.lock().rational_interval_neighbourhood();
        let second_nbd = self.second.lock().rational_interval_neighbourhood();
        match (first_nbd, second_nbd) {
            (Subset::Singleton(first), Subset::Singleton(second)) => {
                Subset::Singleton(first + second)
            }
            (Subset::Singleton(rational), Subset::Interval(interval))
            | (Subset::Interval(interval), Subset::Singleton(rational)) => Subset::Interval(
                RationalInterval::new_unchecked(&rational + interval.a(), rational + interval.b()),
            ),
            (Subset::Interval(first), Subset::Interval(second)) => Subset::Interval(
                RationalInterval::new_unchecked(first.a() + second.a(), first.b() + second.b()),
            ),
        }
    }

    fn length(&self) -> Rational {
        self.first.lock().length() + self.second.lock().length()
    }

    fn refine(&mut self) {
        let first_length = self.first.lock().length();
        let second_length = self.second.lock().length();
        if first_length >= second_length {
            self.first.lock().refine();
        } else {
            self.second.lock().refine();
        }
    }

    fn refine_to_length(&mut self, length: &Rational) {
        let half_length = length * Rational::ONE_HALF;
        self.first.lock().refine_to_length(&half_length);
        self.second.lock().refine_to_length(&half_length);
    }
}

#[derive(Debug)]
struct NegPoint {
    pt: RealApproximatePoint,
}
impl RealApproximatePointInterface for NegPoint {
    fn rational_interval_neighbourhood(&self) -> Subset {
        match self.pt.lock().rational_interval_neighbourhood() {
            Subset::Singleton(rational) => Subset::Singleton(-rational),
            Subset::Interval(interval) => Subset::Interval(RationalInterval::new_unchecked(
                -interval.b(),
                -interval.a(),
            )),
        }
    }

    fn length(&self) -> Rational {
        self.pt.lock().length()
    }

    fn refine(&mut self) {
        self.pt.lock().refine();
    }

    fn refine_to_length(&mut self, length: &Rational) {
        self.pt.lock().refine_to_length(length)
    }
}

#[derive(Debug)]
struct MulPoints {
    first: RealApproximatePoint,
    second: RealApproximatePoint,
}
impl RealApproximatePointInterface for MulPoints {
    fn rational_interval_neighbourhood(&self) -> Subset {
        let first_nbd = self.first.lock().rational_interval_neighbourhood();
        let second_nbd = self.second.lock().rational_interval_neighbourhood();
        match (first_nbd, second_nbd) {
            (Subset::Singleton(first), Subset::Singleton(second)) => {
                Subset::Singleton(first * second)
            }
            (Subset::Singleton(rational), Subset::Interval(interval))
            | (Subset::Interval(interval), Subset::Singleton(rational)) => {
                match rational.cmp(&Rational::ZERO) {
                    std::cmp::Ordering::Less => Subset::Interval(RationalInterval::new_unchecked(
                        &rational * interval.b(),
                        rational * interval.a(),
                    )),
                    std::cmp::Ordering::Equal => Subset::Singleton(Rational::ZERO),
                    std::cmp::Ordering::Greater => {
                        Subset::Interval(RationalInterval::new_unchecked(
                            &rational * interval.a(),
                            rational * interval.b(),
                        ))
                    }
                }
            }
            (Subset::Interval(first), Subset::Interval(second)) => {
                Subset::Interval(RationalInterval::new_unchecked(
                    std::cmp::min(first.a() * second.b(), first.b() * second.a()),
                    std::cmp::max(first.a() * second.a(), first.b() * second.b()),
                ))
            }
        }
    }

    fn length(&self) -> Rational {
        self.rational_interval_neighbourhood().length()
    }

    fn refine(&mut self) {
        let first_length = self.first.lock().length();
        let second_length = self.second.lock().length();
        if first_length >= second_length {
            self.first.lock().refine();
        } else {
            self.second.lock().refine();
        }
    }
}

impl RinglikeSpecializationSignature for RealApproximatePointCanonicalStructure {}

impl ZeroSignature for RealApproximatePointCanonicalStructure {
    fn zero(&self) -> Self::Elem {
        RealApproximatePoint::new(rational::RationalPoint { x: Rational::ZERO })
    }
}

impl AdditionSignature for RealApproximatePointCanonicalStructure {
    fn add(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        RealApproximatePoint::new(AddPoints {
            first: a.clone(),
            second: b.clone(),
        })
    }
}

impl CancellativeAdditionSignature for RealApproximatePointCanonicalStructure {
    fn try_sub(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.sub(a, b))
    }
}

impl TryNegateSignature for RealApproximatePointCanonicalStructure {
    fn try_neg(&self, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.neg(a))
    }
}

impl AdditiveMonoidSignature for RealApproximatePointCanonicalStructure {}

impl AdditiveGroupSignature for RealApproximatePointCanonicalStructure {
    fn neg(&self, a: &Self::Elem) -> Self::Elem {
        RealApproximatePoint::new(NegPoint { pt: a.clone() })
    }
}

impl OneSignature for RealApproximatePointCanonicalStructure {
    fn one(&self) -> Self::Elem {
        RealApproximatePoint::new(rational::RationalPoint { x: Rational::ONE })
    }
}

impl MultiplicationSignature for RealApproximatePointCanonicalStructure {
    fn mul(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        RealApproximatePoint::new(MulPoints {
            first: a.clone(),
            second: b.clone(),
        })
    }
}

impl CommutativeMultiplicationSignature for RealApproximatePointCanonicalStructure {}

impl MultiplicativeMonoidSignature for RealApproximatePointCanonicalStructure {}

impl MultiplicativeAbsorptionMonoidSignature for RealApproximatePointCanonicalStructure {}

impl LeftDistributiveMultiplicationOverAddition for RealApproximatePointCanonicalStructure {}

impl RightDistributiveMultiplicationOverAddition for RealApproximatePointCanonicalStructure {}

impl SemiRingSignature for RealApproximatePointCanonicalStructure {}

impl RingSignature for RealApproximatePointCanonicalStructure {}

impl ToSimpleContinuedFractionSignature for RealApproximatePointCanonicalStructure {}

#[derive(Debug)]
struct InvPoint {
    pt: RealApproximatePoint,
}
impl RealApproximatePointInterface for InvPoint {
    fn rational_interval_neighbourhood(&self) -> Subset {
        loop {
            let nbd = self.pt.lock().rational_interval_neighbourhood();
            match nbd {
                Subset::Singleton(rational) => {
                    return Subset::Singleton(rational.try_reciprocal().expect(
                        "\
Inverse called on an approximate value which later turned out to be exactly 0.",
                    ));
                }
                Subset::Interval(interval) => {
                    let a = interval.a();
                    let b = interval.b();
                    match (a.cmp(&Rational::ZERO), b.cmp(&Rational::ZERO)) {
                        (std::cmp::Ordering::Less, std::cmp::Ordering::Less)
                        | (std::cmp::Ordering::Greater, std::cmp::Ordering::Greater) => {
                            return Subset::Interval(RationalInterval::new_unchecked(
                                b.try_reciprocal().unwrap(),
                                a.try_reciprocal().unwrap(),
                            ));
                        }
                        _ => {
                            self.pt.lock().refine();
                        }
                    }
                }
            }
        }
    }

    fn refine(&mut self) {
        self.pt.lock().refine();
    }
}

impl TryReciprocalSignature for RealApproximatePointCanonicalStructure {
    /// # Warning
    /// May fail to halt if the input is zero.
    fn try_reciprocal(&self, a: &Self::Elem) -> Option<Self::Elem> {
        let nbd = a.lock().rational_interval_neighbourhood();
        match nbd {
            Subset::Singleton(rational) => Some(RealApproximatePoint::new(RationalPoint {
                x: rational.try_reciprocal()?,
            })),
            Subset::Interval(_) => Some(RealApproximatePoint::new(InvPoint { pt: a.clone() })),
        }
    }
}

impl ComplexSubsetSignature for RealApproximatePointCanonicalStructure {
    fn as_f32_real_and_imaginary_parts(&self, z: &Self::Elem) -> (f32, f32) {
        loop {
            let nbd = z.lock().rational_interval_neighbourhood();
            match nbd {
                Subset::Singleton(rational) => {
                    return (rational.as_f32(), 0.0);
                }
                Subset::Interval(interval) => {
                    let a = interval.a().as_f32();
                    let b = interval.b().as_f32();
                    if a == b {
                        return (a, 0.0);
                    }
                    z.lock().refine();
                }
            }
        }
    }

    fn as_f64_real_and_imaginary_parts(&self, z: &Self::Elem) -> (f64, f64) {
        loop {
            let nbd = z.lock().rational_interval_neighbourhood();
            match nbd {
                Subset::Singleton(rational) => {
                    return (rational.as_f64(), 0.0);
                }
                Subset::Interval(interval) => {
                    let a = interval.a().as_f64();
                    let b = interval.b().as_f64();
                    if a == b {
                        return (a, 0.0);
                    }
                    z.lock().refine();
                }
            }
        }
    }
}

impl RealSubsetSignature for RealApproximatePointCanonicalStructure {}

impl RealRoundingSignature for RealApproximatePointCanonicalStructure {
    /// # Warning
    /// May fail to halt on integer inputs.
    fn floor(&self, x: &Self::Elem) -> Integer {
        loop {
            let nbd = x.lock().rational_interval_neighbourhood();
            match nbd {
                Subset::Singleton(rational) => {
                    return rational.floor();
                }
                Subset::Interval(interval) => {
                    let a = algebraeon_structures::Floor::floor(interval.a());
                    let b = algebraeon_structures::Floor::floor(interval.b());
                    if a == b {
                        return a;
                    } else {
                        self.refine(&mut x.clone());
                    }
                }
            }
        }
    }

    /// # Warning!
    /// May fail to halt on integer inputs.
    fn ceil(&self, x: &Self::Elem) -> Integer {
        loop {
            let nbd = x.lock().rational_interval_neighbourhood();
            match nbd {
                Subset::Singleton(rational) => {
                    return rational.ceil();
                }
                Subset::Interval(interval) => {
                    let a = algebraeon_structures::Ceil::ceil(interval.a());
                    let b = algebraeon_structures::Ceil::ceil(interval.b());
                    if a == b {
                        return a;
                    } else {
                        self.refine(&mut x.clone());
                    }
                }
            }
        }
    }

    fn round(&self, x: &Self::Elem) -> Integer {
        self.floor(&self.add(
            x,
            &RealApproximatePoint::new(RationalPoint {
                x: Rational::ONE_HALF,
            }),
        ))
    }
}

pub mod continued_fraction;
pub mod pi;
pub mod rational;
