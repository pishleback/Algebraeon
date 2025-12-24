use super::open_interval::OpenRationalInterval;
use crate::{
    coontinued_fractions::SimpleContinuedFraction,
    structure::{
        AdditiveGroupSignature, AdditiveMonoidSignature, ComplexSubsetSignature, MetaRealSubset,
        RealSubsetSignature, RingSignature, SemiRingSignature,
    },
};
use algebraeon_nzq::{Integer, Rational, RationalCanonicalStructure};
use algebraeon_sets::{
    approximations::{ApproximatePointsSignature, SubsetsSignature},
    structure::{SetSignature, Signature},
};
use std::{
    fmt::Debug,
    sync::{Arc, Mutex},
};

#[derive(Debug, Clone)]
pub enum Subset {
    Singleton(Rational),
    Interval(OpenRationalInterval),
}

impl Subset {
    pub fn length(&self) -> Rational {
        match self {
            Subset::Singleton(_) => Rational::ZERO,
            Subset::Interval(interval) => interval.b() - interval.a(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SubsetsStructure {}

pub fn subsets() -> SubsetsStructure {
    SubsetsStructure {}
}

impl Signature for SubsetsStructure {}

impl SetSignature for SubsetsStructure {
    type Set = Subset;

    fn is_element(&self, _: &Self::Set) -> Result<(), String> {
        Ok(())
    }
}

impl SubsetsSignature for SubsetsStructure {}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PointsStructure {}

pub fn points() -> PointsStructure {
    PointsStructure {}
}

pub fn point<P: Point + 'static>(p: P) -> <PointsStructure as SetSignature>::Set {
    Arc::new(Mutex::new(p))
}

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
            | (Subset::Interval(interval), Subset::Singleton(rational)) => {
                Subset::Interval(OpenRationalInterval::new_unchecked(
                    &rational + interval.a(),
                    rational + interval.b(),
                ))
            }
            (Subset::Interval(first), Subset::Interval(second)) => Subset::Interval(
                OpenRationalInterval::new_unchecked(first.a() + second.a(), first.b() + second.b()),
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
            Subset::Interval(interval) => Subset::Interval(OpenRationalInterval::new_unchecked(
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
                    std::cmp::Ordering::Less => {
                        Subset::Interval(OpenRationalInterval::new_unchecked(
                            &rational * interval.b(),
                            rational * interval.a(),
                        ))
                    }
                    std::cmp::Ordering::Equal => Subset::Singleton(Rational::ZERO),
                    std::cmp::Ordering::Greater => {
                        Subset::Interval(OpenRationalInterval::new_unchecked(
                            &rational * interval.a(),
                            rational * interval.b(),
                        ))
                    }
                }
            }
            (Subset::Interval(first), Subset::Interval(second)) => {
                Subset::Interval(OpenRationalInterval::new_unchecked(
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
        point(RationalPoint { x: Rational::ZERO })
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
        point(RationalPoint { x: Rational::ONE })
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        point(MulPoints {
            first: a.clone(),
            second: b.clone(),
        })
    }
}

impl RingSignature for PointsStructure {}

impl ComplexSubsetSignature for PointsStructure {
    fn as_f32_real_and_imaginary_parts(&self, z: &Self::Set) -> (f32, f32) {
        let mut z = z.lock().unwrap();
        z.refine_to_length(&Rational::from_integers(
            Integer::ONE,
            Integer::from(100000000i64),
        ));
        let rat = match z.rational_interval_neighbourhood() {
            Subset::Singleton(rational) => rational,
            Subset::Interval(interval) => Rational::ONE_HALF * (interval.a() + interval.b()),
        };
        (rat.as_f32(), 0.0)
    }

    fn as_f64_real_and_imaginary_parts(&self, z: &Self::Set) -> (f64, f64) {
        let mut z = z.lock().unwrap();
        z.refine_to_length(&Rational::from_integers(
            Integer::ONE,
            Integer::from(10000000000000000i64),
        ));
        let rat = match z.rational_interval_neighbourhood() {
            Subset::Singleton(rational) => rational,
            Subset::Interval(interval) => Rational::ONE_HALF * (interval.a() + interval.b()),
        };
        (rat.as_f64(), 0.0)
    }
}

impl RealSubsetSignature for PointsStructure {}

#[derive(Debug, Clone)]
pub struct RationalPoint {
    pub x: Rational,
}

impl Point for RationalPoint {
    fn rational_interval_neighbourhood(&self) -> Subset {
        Subset::Singleton(self.x.clone())
    }

    fn length(&self) -> Rational {
        Rational::ZERO
    }

    fn refine(&mut self) {
        // Nothing needs to be done
    }
}

#[derive(Debug, Clone)]
pub struct SimpleContinuedFractionPoint<SCFG: SimpleContinuedFraction> {
    coeffs: SCFG,
    n: usize,
    // pn / pn is the latest rational approximation
    pn: Integer,
    qn: Integer,
    // This is the previous pn / qn
    pn_prev: Integer,
    qn_prev: Integer,
    // The sequence of rationals pn/qn aternates around its limit
    // so pn/qn and pn_prev/qn_prev define a shortening open interval containing the limit
}

impl<SCFG: SimpleContinuedFraction> From<SCFG> for SimpleContinuedFractionPoint<SCFG> {
    fn from(mut coeffs: SCFG) -> Self {
        let a1 = coeffs.next();
        let a2 = coeffs.next();
        // p0/q0  =  1 / 0
        // p1/q1  =  a1 / 1
        // p2/q2  =  a2*p1+p0 / a2*q1 + q0  =  a2*a1+1 / a2
        let p1 = a1.clone();
        let q1 = Integer::ONE;
        let p2 = &a2 * a1 + Integer::ONE;
        let q2 = a2;
        Self {
            coeffs,
            n: 2,
            pn: p2,
            qn: q2,
            pn_prev: p1,
            qn_prev: q1,
        }
    }
}

impl<SCFG: SimpleContinuedFraction> SimpleContinuedFractionPoint<SCFG> {
    pub fn a(&self) -> Rational {
        if !self.n.is_multiple_of(2) {
            Rational::from_integers(&self.pn, &self.qn)
        } else {
            Rational::from_integers(&self.pn_prev, &self.qn_prev)
        }
    }

    pub fn b(&self) -> Rational {
        if self.n.is_multiple_of(2) {
            Rational::from_integers(&self.pn, &self.qn)
        } else {
            Rational::from_integers(&self.pn_prev, &self.qn_prev)
        }
    }

    pub fn refine(&mut self) {
        self.n += 1;
        let an = self.coeffs.next();
        // p_n / q_n  =  a_n*p_{n-1}+p_{n-2} / a_n*q_{n-1}+q_{n-2}
        let pn = &an * &self.pn + &self.pn_prev;
        let qn = an * &self.qn + &self.qn_prev;
        std::mem::swap(&mut self.pn_prev, &mut self.pn);
        std::mem::swap(&mut self.qn_prev, &mut self.qn);
        self.pn = pn;
        self.qn = qn;
    }
}

impl<SCFG: SimpleContinuedFraction> Point for SimpleContinuedFractionPoint<SCFG> {
    fn rational_interval_neighbourhood(&self) -> Subset {
        Subset::Interval(OpenRationalInterval::new_unchecked(self.a(), self.b()))
    }

    fn refine(&mut self) {
        self.refine();
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use algebraeon_nzq::Natural;

    use crate::coontinued_fractions::eulers_constant;

    use super::*;

    #[test]
    fn test_eulers_constant() {
        let mut e = SimpleContinuedFractionPoint::from(eulers_constant());

        // 2 < e < 3
        if let Subset::Interval(interval) = e.rational_interval_neighbourhood() {
            assert_eq!(interval.a(), &Integer::from(2));
            assert_eq!(interval.b(), &Integer::from(3));
        } else {
            panic!()
        }

        e.refine();

        // 8/3 < e < 3
        if let Subset::Interval(interval) = e.rational_interval_neighbourhood() {
            assert_eq!(interval.a(), &Rational::from_str("8/3").unwrap());
            assert_eq!(interval.b(), &Integer::from(3));
        } else {
            panic!()
        }

        e.refine();

        // 8/3 < e < 11/4
        if let Subset::Interval(interval) = e.rational_interval_neighbourhood() {
            assert_eq!(interval.a(), &Rational::from_str("8/3").unwrap());
            assert_eq!(interval.b(), &Rational::from_str("11/4").unwrap());
        } else {
            panic!()
        }

        assert_eq!(points().as_f32(&point(e.clone())), std::f32::consts::E);
        assert_eq!(points().as_f64(&point(e)), std::f64::consts::E);
    }

    #[test]
    fn test_ring_opps() {
        let e = point(SimpleContinuedFractionPoint::from(eulers_constant()));
        // compute e+(-e)^3
        let v = points().add(
            &e,
            &points().nat_pow(&points().neg(&e), &Natural::from(3u32)),
        );
        let f = points().as_f64(&v);
        debug_assert_eq!(f, -17.367_255_094_728_623);
    }
}
