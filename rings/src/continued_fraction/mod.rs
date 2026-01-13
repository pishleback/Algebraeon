use crate::structure::{FieldSignature, RealRoundingSignature, RingUnitsSignature};
use algebraeon_macros::signature_meta_trait;
use algebraeon_nzq::{Integer, Rational};
use algebraeon_sets::structure::MetaType;
use std::{
    borrow::Cow,
    fmt::Debug,
    sync::{Arc, Mutex},
};

pub trait SimpleContinuedFraction: Debug + Clone + Send + Sync {
    /// Return the nth continued fraction coefficient
    /// n=0 must return Some(any Integer)
    /// n>0 must return Some(any Integer >= 1) or None
    /// Calls with the same n must return the same Integer.
    fn coeff(&'_ self, n: usize) -> Option<Cow<'_, Integer>>;

    fn iter(&'_ self) -> impl Iterator<Item = Cow<'_, Integer>> {
        (0usize..)
            .map(|n| self.coeff(n))
            .take_while(|c| c.is_some())
            .map(|c| c.unwrap())
    }

    fn rational_approximations(self) -> RationalApproximations<Self>
    where
        Self: Sized,
    {
        RationalApproximations::new(self)
    }
}

#[derive(Debug, Clone)]
pub struct RationalApproximations<SCF: SimpleContinuedFraction> {
    pqs: Arc<Mutex<Vec<(Integer, Integer)>>>,
    scf: SCF,
}

impl<SCF: SimpleContinuedFraction> RationalApproximations<SCF> {
    pub fn new(scf: SCF) -> Self {
        let a0 = scf.coeff(0).unwrap();
        let a1 = scf.coeff(1).unwrap();

        let q0 = Integer::from(1);
        let p1 = a1.as_ref() * a0.as_ref() + Integer::ONE;
        let p0 = a0.into_owned();
        let q1 = a1.into_owned();

        Self {
            pqs: Arc::new(Mutex::new(vec![(p0, q0), (p1, q1)])),
            scf,
        }
    }

    pub fn get_rat(&self, n: usize) -> Rational {
        let mut pqs = self.pqs.lock().unwrap();
        while pqs.len() < n + 1 {
            let a0 = self.scf.coeff(pqs.len()).unwrap();
            let (p1, q1) = &pqs[pqs.len() - 1];
            let (p2, q2) = &pqs[pqs.len() - 2];
            let p0 = a0.as_ref() * p1 + p2;
            let q0 = a0.as_ref() * q1 + q2;
            pqs.push((p0, q0));
        }
        Rational::from_integers(&pqs[n].0, &pqs[n].1)
    }
}

#[derive(Debug, Clone)]
pub struct RationalSimpleContinuedFraction {
    coeffs: Vec<Integer>,
}

impl SimpleContinuedFraction for RationalSimpleContinuedFraction {
    fn coeff(&'_ self, n: usize) -> Option<Cow<'_, Integer>> {
        Some(Cow::Borrowed(self.coeffs.get(n)?))
    }
}

#[derive(Debug, Clone)]
pub struct PeriodicSimpleContinuedFraction {
    initial: Vec<Integer>,
    repeats: Vec<Integer>,
}

impl PeriodicSimpleContinuedFraction {
    pub fn new(initial: Vec<Integer>, repeats: Vec<Integer>) -> Result<Self, ()> {
        #[allow(clippy::needless_range_loop)]
        for i in 1..initial.len() {
            if initial[i] <= Integer::ZERO {
                return Err(());
            }
        }
        if repeats.is_empty() {
            return Err(());
        }
        #[allow(clippy::needless_range_loop)]
        for i in 0..repeats.len() {
            if repeats[i] <= Integer::ZERO {
                return Err(());
            }
        }
        Ok(Self { initial, repeats })
    }

    pub fn new_unchecked(initial: Vec<Integer>, repeats: Vec<Integer>) -> Self {
        #[cfg(debug_assertions)]
        Self::new(initial.clone(), repeats.clone()).unwrap();
        Self { initial, repeats }
    }
}

impl SimpleContinuedFraction for PeriodicSimpleContinuedFraction {
    fn coeff(&'_ self, n: usize) -> Option<Cow<'_, Integer>> {
        Some(Cow::Borrowed(if let Some(c) = self.initial.get(n) {
            c
        } else {
            &self.repeats[(n - self.initial.len()) % self.repeats.len()]
        }))
    }
}

pub trait IrrationalSimpleContinuedFractionGenerator: Debug + Send + Sync {
    /// First value returned can be any integer
    /// All subsequent values must be >= 1
    fn next(&mut self) -> Integer;
    // fn into_iter(mut self) -> impl Iterator<Item = Integer>
    // where
    //     Self: 'static,
    // {
    //     (0usize..).map(move |_| self.next())
    // }
    fn into_continued_fraction(self) -> IrrationalSimpleContinuedFraction
    where
        Self: Sized + 'static,
    {
        self.into()
    }
}

#[derive(Debug)]
struct IrrationalSimpleContinuedFractionCache {
    coeffs: Vec<Integer>,
    coeff_gen: Box<dyn IrrationalSimpleContinuedFractionGenerator>,
}

#[derive(Debug, Clone)]
pub struct IrrationalSimpleContinuedFraction {
    cache: Arc<Mutex<IrrationalSimpleContinuedFractionCache>>,
}

impl<G: IrrationalSimpleContinuedFractionGenerator + 'static> From<G>
    for IrrationalSimpleContinuedFraction
{
    fn from(g: G) -> Self {
        Self {
            cache: Arc::new(Mutex::new(IrrationalSimpleContinuedFractionCache {
                coeffs: vec![],
                coeff_gen: Box::new(g),
            })),
        }
    }
}

impl SimpleContinuedFraction for IrrationalSimpleContinuedFraction {
    fn coeff(&'_ self, n: usize) -> Option<Cow<'_, Integer>> {
        let mut cache = self.cache.lock().unwrap();
        while cache.coeffs.len() < n + 1 {
            let c = cache.coeff_gen.next();
            cache.coeffs.push(c);
        }
        Some(Cow::Owned(cache.coeffs[n].clone()))
    }
}

mod eulers_constant;

pub fn eulers_constant() -> IrrationalSimpleContinuedFraction {
    eulers_constant::EulersConstantSimpleContinuedFractionGenerator::new().into()
}

#[derive(Debug)]
struct SimpleContinuedFractionFromRealStructureCache<R: ToSimpleContinuedFractionSignature> {
    value: Option<R::Set>, // terminated iff None
    coeffs: Vec<Integer>,
}

#[derive(Debug, Clone)]
pub struct SimpleContinuedFractionFromRealStructure<R: ToSimpleContinuedFractionSignature> {
    ring: R,
    cache: Arc<Mutex<SimpleContinuedFractionFromRealStructureCache<R>>>,
}

impl<R: ToSimpleContinuedFractionSignature> SimpleContinuedFractionFromRealStructure<R> {
    fn new(ring: R, value: R::Set) -> Self {
        Self {
            ring,
            cache: Arc::new(Mutex::new(SimpleContinuedFractionFromRealStructureCache {
                value: Some(value),
                coeffs: vec![],
            })),
        }
    }
}

impl<R: ToSimpleContinuedFractionSignature> SimpleContinuedFraction
    for SimpleContinuedFractionFromRealStructure<R>
{
    fn coeff(&'_ self, n: usize) -> Option<Cow<'_, Integer>> {
        let mut cache = self.cache.lock().unwrap();
        while cache.coeffs.len() < n + 1
            && let Some(cache_value) = cache.value.as_ref()
        {
            let c = self.ring.floor(cache_value);
            if let Some(value) = self
                .ring
                .try_inv(&self.ring.sub(cache_value, &self.ring.from_int(&c)))
            {
                cache.value = Some(value);
            } else {
                cache.value = None;
            }
            cache.coeffs.push(c);
        }
        Some(Cow::Owned(cache.coeffs.get(n)?.clone()))
    }
}

/// Implementing this trait is only valid if self.try_inv only returns None when given 0
#[signature_meta_trait]
pub trait ToSimpleContinuedFractionSignature: RealRoundingSignature + RingUnitsSignature {
    /// # Warning
    /// `value` should be irrational.
    /// If `value` is not irrational then calls to `.next()` on the resulting `SimpleContinuedFraction` may panic.
    fn simple_continued_fraction(
        &self,
        value: Self::Set,
    ) -> SimpleContinuedFractionFromRealStructure<Self> {
        SimpleContinuedFractionFromRealStructure::new(self.clone(), value)
    }
}
impl<R: RealRoundingSignature + FieldSignature + RingUnitsSignature>
    ToSimpleContinuedFractionSignature for R
{
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use super::*;

    fn rat_cf(rat: &'static str) -> Vec<Integer> {
        Rational::from_str(rat)
            .unwrap()
            .simple_continued_fraction()
            .iter()
            .map(|x| x.as_ref().clone())
            .collect::<Vec<_>>()
    }

    #[test]
    fn test_rational() {
        assert_eq!(rat_cf("-2"), vec![Integer::from(-2)]);
        assert_eq!(rat_cf("-1"), vec![Integer::from(-1)]);
        assert_eq!(rat_cf("0"), vec![Integer::from(0)]);
        assert_eq!(rat_cf("1"), vec![Integer::from(1)]);
        assert_eq!(rat_cf("2"), vec![Integer::from(2)]);

        assert_eq!(
            rat_cf("9/10"),
            vec![Integer::from(0), Integer::from(1), Integer::from(9)]
        );

        assert_eq!(rat_cf("-9/10"), vec![Integer::from(-1), Integer::from(10)]);

        assert_eq!(
            rat_cf("-5678/1234"),
            vec![
                Integer::from(-5),
                Integer::from(2),
                Integer::from(1),
                Integer::from(1),
                Integer::from(30),
                Integer::from(4)
            ]
        );
    }

    #[test]
    fn test_periodic() {
        assert_eq!(
            PeriodicSimpleContinuedFraction::new(
                vec![Integer::from(1), Integer::from(2)],
                vec![Integer::from(3), Integer::from(4), Integer::from(5)],
            )
            .unwrap()
            .iter()
            .take(10)
            .map(|x| x.as_ref().clone())
            .collect::<Vec<_>>(),
            vec![
                Integer::from(1),
                Integer::from(2),
                Integer::from(3),
                Integer::from(4),
                Integer::from(5),
                Integer::from(3),
                Integer::from(4),
                Integer::from(5),
                Integer::from(3),
                Integer::from(4)
            ]
        );

        assert_eq!(
            PeriodicSimpleContinuedFraction::new(
                vec![],
                vec![Integer::from(3), Integer::from(4), Integer::from(5)]
            )
            .unwrap()
            .iter()
            .take(7)
            .map(|x| x.as_ref().clone())
            .collect::<Vec<_>>(),
            vec![
                Integer::from(3),
                Integer::from(4),
                Integer::from(5),
                Integer::from(3),
                Integer::from(4),
                Integer::from(5),
                Integer::from(3)
            ]
        );

        assert!(PeriodicSimpleContinuedFraction::new(vec![], vec![]).is_err());
        assert!(PeriodicSimpleContinuedFraction::new(vec![], vec![Integer::from(-1)]).is_err());
        assert!(
            PeriodicSimpleContinuedFraction::new(
                vec![Integer::from(1), Integer::from(-1)],
                vec![Integer::from(1)]
            )
            .is_err()
        );
    }
}
