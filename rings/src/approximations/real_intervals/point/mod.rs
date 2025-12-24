use algebraeon_nzq::{Integer, Rational, RationalCanonicalStructure, traits::Floor};
use algebraeon_sets::{
    approximations::ApproximatePointsSignature,
    structure::{SetSignature, Signature},
};
use std::{
    fmt::Debug,
    sync::{Arc, Mutex},
};

use crate::{
    approximations::{
        rational_interval::RationalInterval,
        real_intervals::{RationalPoint, Subset, SubsetsStructure},
    },
    structure::{
        AdditiveGroupSignature, AdditiveMonoidSignature, ComplexSubsetSignature, MetaRealRounding,
        MetaRealSubset, MetaSemiRingUnitsSignature, RealRoundingSignature, RealSubsetSignature,
        RingSignature, SemiRingSignature, SemiRingUnitsSignature,
    },
};

pub trait Point: Debug + Send + Sync {
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

pub fn point<P: Point + 'static>(p: P) -> <PointsStructure as SetSignature>::Set {
    Arc::new(Mutex::new(p))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PointsStructure {}

pub fn points() -> PointsStructure {
    PointsStructure {}
}

impl Signature for PointsStructure {}

impl SetSignature for PointsStructure {
    type Set = Arc<Mutex<dyn Point>>;

    fn is_element(&self, _: &Self::Set) -> Result<(), String> {
        Ok(())
    }
}

impl ApproximatePointsSignature for PointsStructure {
    type Precision = RationalCanonicalStructure;
    type OpenSubsetsStructure = SubsetsStructure;

    fn open_neighbourhood(&self, approx_point: &Self::Set) -> Subset {
        approx_point
            .lock()
            .unwrap()
            .rational_interval_neighbourhood()
    }

    fn precision(&self, approx_point: &Self::Set) -> Rational {
        approx_point.lock().unwrap().length()
    }

    fn refine(&self, approx_point: &mut Self::Set) {
        approx_point.lock().unwrap().refine();
    }

    fn refine_to(&self, approx_point: &mut Self::Set, length: &Rational) {
        approx_point.lock().unwrap().refine_to_length(length);
    }
}

#[derive(Debug)]
pub struct AddPoints {
    first: Arc<Mutex<dyn Point>>,
    second: Arc<Mutex<dyn Point>>,
}
impl Point for AddPoints {
    fn rational_interval_neighbourhood(&self) -> Subset {
        let first_nbd = self.first.lock().unwrap().rational_interval_neighbourhood();
        let second_nbd = self
            .second
            .lock()
            .unwrap()
            .rational_interval_neighbourhood();
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
        self.first.lock().unwrap().length() + self.second.lock().unwrap().length()
    }

    fn refine(&mut self) {
        let first_length = self.first.lock().unwrap().length();
        let second_length = self.second.lock().unwrap().length();
        if first_length >= second_length {
            self.first.lock().unwrap().refine();
        } else {
            self.second.lock().unwrap().refine();
        }
    }

    fn refine_to_length(&mut self, length: &Rational) {
        let half_length = length * Rational::ONE_HALF;
        self.first.lock().unwrap().refine_to_length(&half_length);
        self.second.lock().unwrap().refine_to_length(&half_length);
    }
}

#[derive(Debug)]
pub struct NegPoint {
    pt: Arc<Mutex<dyn Point>>,
}
impl Point for NegPoint {
    fn rational_interval_neighbourhood(&self) -> Subset {
        match self.pt.lock().unwrap().rational_interval_neighbourhood() {
            Subset::Singleton(rational) => Subset::Singleton(-rational),
            Subset::Interval(interval) => Subset::Interval(RationalInterval::new_unchecked(
                -interval.b(),
                -interval.a(),
            )),
        }
    }

    fn length(&self) -> Rational {
        self.pt.lock().unwrap().length()
    }

    fn refine(&mut self) {
        self.pt.lock().unwrap().refine();
    }

    fn refine_to_length(&mut self, length: &Rational) {
        self.pt.lock().unwrap().refine_to_length(length)
    }
}

#[derive(Debug)]
pub struct MulPoints {
    first: Arc<Mutex<dyn Point>>,
    second: Arc<Mutex<dyn Point>>,
}
impl Point for MulPoints {
    fn rational_interval_neighbourhood(&self) -> Subset {
        let first_nbd = self.first.lock().unwrap().rational_interval_neighbourhood();
        let second_nbd = self
            .second
            .lock()
            .unwrap()
            .rational_interval_neighbourhood();
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
        let first_length = self.first.lock().unwrap().length();
        let second_length = self.second.lock().unwrap().length();
        if first_length >= second_length {
            self.first.lock().unwrap().refine();
        } else {
            self.second.lock().unwrap().refine();
        }
    }
}

impl AdditiveMonoidSignature for PointsStructure {
    fn zero(&self) -> Self::Set {
        point(rational::RationalPoint { x: Rational::ZERO })
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        point(AddPoints {
            first: a.clone(),
            second: b.clone(),
        })
    }
}

impl AdditiveGroupSignature for PointsStructure {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        point(NegPoint { pt: a.clone() })
    }
}

impl SemiRingSignature for PointsStructure {
    fn one(&self) -> Self::Set {
        point(rational::RationalPoint { x: Rational::ONE })
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        point(MulPoints {
            first: a.clone(),
            second: b.clone(),
        })
    }
}

impl RingSignature for PointsStructure {}

#[derive(Debug)]
pub struct InvPoint {
    pt: Arc<Mutex<dyn Point>>,
}
impl Point for InvPoint {
    fn rational_interval_neighbourhood(&self) -> Subset {
        loop {
            let nbd = self.pt.lock().unwrap().rational_interval_neighbourhood();
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
                            self.pt.lock().unwrap().refine();
                        }
                    }
                }
            }
        }
    }

    fn refine(&mut self) {
        self.pt.lock().unwrap().refine();
    }
}

impl SemiRingUnitsSignature for PointsStructure {
    /// # Warning
    /// Mail fail to halt if the input is zero.
    fn inv(&self, a: &Self::Set) -> Result<Self::Set, crate::structure::RingDivisionError> {
        let nbd = a.lock().unwrap().rational_interval_neighbourhood();
        match nbd {
            Subset::Singleton(rational) => Ok(point(RationalPoint { x: rational.inv()? })),
            Subset::Interval(_) => Ok(point(InvPoint { pt: a.clone() })),
        }
    }
}

impl ComplexSubsetSignature for PointsStructure {
    fn as_f32_real_and_imaginary_parts(&self, z: &Self::Set) -> (f32, f32) {
        loop {
            let nbd = z.lock().unwrap().rational_interval_neighbourhood();
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
                    z.lock().unwrap().refine();
                }
            }
        }
    }

    fn as_f64_real_and_imaginary_parts(&self, z: &Self::Set) -> (f64, f64) {
        loop {
            let nbd = z.lock().unwrap().rational_interval_neighbourhood();
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
                    z.lock().unwrap().refine();
                }
            }
        }
    }
}

impl RealSubsetSignature for PointsStructure {}

impl RealRoundingSignature for PointsStructure {
    /// # Warning
    /// May fail to halt on integer inputs.
    fn floor(&self, x: &Self::Set) -> Integer {
        loop {
            let nbd = x.lock().unwrap().rational_interval_neighbourhood();
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
            let nbd = x.lock().unwrap().rational_interval_neighbourhood();
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
            &point(RationalPoint {
                x: Rational::ONE_HALF,
            }),
        ))
    }
}

pub mod continued_fraction;
pub mod pi;
pub mod rational;
