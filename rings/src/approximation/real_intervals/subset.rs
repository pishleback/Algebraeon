use crate::approximation::rational_interval::RationalInterval;
use algebraeon_nzq::Rational;
use algebraeon_sets::{
    approximations::SubsetsSignature,
    structure::{SetSignature, Signature},
};

#[derive(Debug, Clone)]
pub enum Subset {
    Singleton(Rational),
    Interval(RationalInterval),
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
