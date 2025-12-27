use crate::{
    approximation::{
        rational_interval::RationalInterval,
        real_intervals::{RationalPoint, Subset, SubsetsStructure},
    },
    continued_fraction::SimpleContinuedFraction,
    structure::{
        AdditiveGroupSignature, AdditiveMonoidSignature, ComplexSubsetSignature, MetaRealRounding,
        MetaRealSubset, MetaSemiRingUnitsSignature, RealRoundingSignature, RealSubsetSignature,
        RingSignature, SemiRingSignature, SemiRingUnitsSignature,
    },
};
use algebraeon_nzq::{Integer, Rational, RationalCanonicalStructure, traits::Floor};
use algebraeon_sets::{
    approximations::ApproximatePointsSignature,
    structure::{CanonicalStructure, MetaType, SetSignature, Signature},
};
use std::{
    fmt::Debug,
    sync::{Arc, Mutex, MutexGuard},
};

pub trait PointInterface: Debug + Send + Sync {
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
pub struct Point {
    repr: Arc<Mutex<dyn PointInterface>>,
}

impl Point {
    pub fn new<R: PointInterface + 'static>(repr: R) -> Self {
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

    pub fn lock(&self) -> MutexGuard<'_, dyn PointInterface + 'static> {
        self.repr.lock().unwrap()
    }
}

impl ApproximatePointsSignature for PointCanonicalStructure {
    type Precision = RationalCanonicalStructure;
    type OpenSubsetsStructure = SubsetsStructure;

    fn open_neighbourhood(&self, approx_point: &Self::Set) -> Subset {
        approx_point.lock().rational_interval_neighbourhood()
    }

    fn precision(&self, approx_point: &Self::Set) -> Rational {
        approx_point.lock().length()
    }

    fn refine(&self, approx_point: &mut Self::Set) {
        approx_point.lock().refine();
    }

    fn refine_to(&self, approx_point: &mut Self::Set, length: &Rational) {
        approx_point.lock().refine_to_length(length);
    }
}

#[derive(Debug)]
pub struct AddPoints {
    first: Point,
    second: Point,
}
impl PointInterface for AddPoints {
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
pub struct NegPoint {
    pt: Point,
}
impl PointInterface for NegPoint {
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
pub struct MulPoints {
    first: Point,
    second: Point,
}
impl PointInterface for MulPoints {
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

impl AdditiveMonoidSignature for PointCanonicalStructure {
    fn zero(&self) -> Self::Set {
        Point::new(rational::RationalPoint { x: Rational::ZERO })
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        Point::new(AddPoints {
            first: a.clone(),
            second: b.clone(),
        })
    }
}

impl AdditiveGroupSignature for PointCanonicalStructure {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        Point::new(NegPoint { pt: a.clone() })
    }
}

impl SemiRingSignature for PointCanonicalStructure {
    fn one(&self) -> Self::Set {
        Point::new(rational::RationalPoint { x: Rational::ONE })
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        Point::new(MulPoints {
            first: a.clone(),
            second: b.clone(),
        })
    }
}

impl RingSignature for PointCanonicalStructure {}

#[derive(Debug)]
pub struct InvPoint {
    pt: Point,
}
impl PointInterface for InvPoint {
    fn rational_interval_neighbourhood(&self) -> Subset {
        loop {
            let nbd = self.pt.lock().rational_interval_neighbourhood();
            match nbd {
                Subset::Singleton(rational) => {
                    return Subset::Singleton(rational.inv().expect(
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
                                b.inv().unwrap(),
                                a.inv().unwrap(),
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

impl SemiRingUnitsSignature for PointCanonicalStructure {
    /// # Warning
    /// May fail to halt if the input is zero.
    fn inv(&self, a: &Self::Set) -> Result<Self::Set, crate::structure::RingDivisionError> {
        let nbd = a.lock().rational_interval_neighbourhood();
        match nbd {
            Subset::Singleton(rational) => Ok(Point::new(RationalPoint { x: rational.inv()? })),
            Subset::Interval(_) => Ok(Point::new(InvPoint { pt: a.clone() })),
        }
    }
}

impl ComplexSubsetSignature for PointCanonicalStructure {
    fn as_f32_real_and_imaginary_parts(&self, z: &Self::Set) -> (f32, f32) {
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

    fn as_f64_real_and_imaginary_parts(&self, z: &Self::Set) -> (f64, f64) {
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

impl RealSubsetSignature for PointCanonicalStructure {}

impl RealRoundingSignature for PointCanonicalStructure {
    /// # Warning
    /// May fail to halt on integer inputs.
    fn floor(&self, x: &Self::Set) -> Integer {
        loop {
            let nbd = x.lock().rational_interval_neighbourhood();
            match nbd {
                Subset::Singleton(rational) => {
                    return rational.floor();
                }
                Subset::Interval(interval) => {
                    let a = algebraeon_nzq::traits::Floor::floor(interval.a());
                    let b = algebraeon_nzq::traits::Floor::floor(interval.b());
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
    fn ceil(&self, x: &Self::Set) -> Integer {
        loop {
            let nbd = x.lock().rational_interval_neighbourhood();
            match nbd {
                Subset::Singleton(rational) => {
                    return rational.ceil();
                }
                Subset::Interval(interval) => {
                    let a = algebraeon_nzq::traits::Ceil::ceil(interval.a());
                    let b = algebraeon_nzq::traits::Ceil::ceil(interval.b());
                    if a == b {
                        return a;
                    } else {
                        self.refine(&mut x.clone());
                    }
                }
            }
        }
    }

    fn round(&self, x: &Self::Set) -> Integer {
        self.floor(&self.add(
            x,
            &Point::new(RationalPoint {
                x: Rational::ONE_HALF,
            }),
        ))
    }
}

pub mod continued_fraction;
pub mod pi;
pub mod rational;
