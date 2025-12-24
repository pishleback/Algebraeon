use algebraeon_nzq::Integer;
use std::fmt::Debug;

use crate::structure::{RealRoundingSignature, RingUnitsSignature};

pub trait SimpleContinuedFraction: Debug + Clone + Send + Sync {
    /// First value returned can be any integer
    /// All subsequent values must be >= 1
    fn next(&mut self) -> Integer;
    fn into_iter(mut self) -> impl Iterator<Item = Integer>
    where
        Self: 'static,
    {
        (0usize..).map(move |_| self.next())
    }
}

mod eulers_constant;

pub fn eulers_constant() -> impl SimpleContinuedFraction {
    eulers_constant::EulersConstantSimpleContinuedFraction::new()
}

#[derive(Debug, Clone)]
pub struct SimpleContinuedFractionFromRealStructure<R: ToSimpleContinuedFractionSignature> {
    ring: R,
    value: R::Set,
}

impl<R: ToSimpleContinuedFractionSignature> SimpleContinuedFraction
    for SimpleContinuedFractionFromRealStructure<R>
{
    fn next(&mut self) -> Integer {
        let c = self.ring.floor(&self.value);
        self.value = self
            .ring
            .inv(&self.ring.sub(&self.value, &self.ring.from_int(&c)))
            .unwrap();
        c
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
        SimpleContinuedFractionFromRealStructure {
            ring: self.clone(),
            value,
        }
    }
}
impl<R: RealRoundingSignature + RingUnitsSignature> ToSimpleContinuedFractionSignature for R {}
