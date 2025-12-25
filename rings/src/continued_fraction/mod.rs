use algebraeon_nzq::{Integer, Rational};
use algebraeon_sets::structure::MetaType;
use std::{
    borrow::Cow,
    fmt::Debug,
    sync::{Arc, Mutex},
};

use crate::structure::{RealRoundingSignature, RingDivisionError, RingUnitsSignature};

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

impl IrrationalSimpleContinuedFraction {
    pub fn new<G: IrrationalSimpleContinuedFractionGenerator + 'static>(g: G) -> Self {
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
    IrrationalSimpleContinuedFraction::new(
        eulers_constant::EulersConstantSimpleContinuedFractionGenerator::new(),
    )
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
            match self
                .ring
                .inv(&self.ring.sub(cache_value, &self.ring.from_int(&c)))
            {
                Ok(value) => {
                    cache.value = Some(value);
                }
                Err(RingDivisionError::DivideByZero) => {
                    cache.value = None;
                }
                Err(RingDivisionError::NotDivisible) => unreachable!(),
            }
            cache.coeffs.push(c);
        }
        Some(Cow::Owned(cache.coeffs.get(n)?.clone()))
    }
}

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
impl<R: RealRoundingSignature + RingUnitsSignature> ToSimpleContinuedFractionSignature for R {}

pub trait MetaToSimpleContinuedFraction: MetaType
where
    Self::Signature: ToSimpleContinuedFractionSignature,
{
    /// # Warning
    /// `value` should be irrational.
    /// If `value` is not irrational then calls to `.next()` on the resulting `SimpleContinuedFraction` may panic.
    fn simple_continued_fraction(
        self,
    ) -> SimpleContinuedFractionFromRealStructure<Self::Signature> {
        Self::structure().simple_continued_fraction(self)
    }
}
impl<R: MetaType> MetaToSimpleContinuedFraction for R where
    Self::Signature: ToSimpleContinuedFractionSignature
{
}
