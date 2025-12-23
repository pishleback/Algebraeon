use algebraeon_nzq::{Rational, RationalCanonicalStructure};
use algebraeon_sets::{
    approximations::{ApproximatePointsSignature, SubsetsSignature},
    structure::{SetSignature, Signature},
};

use super::open_interval::OpenRationalInterval;
use std::{
    fmt::Debug,
    sync::{Arc, Mutex},
};

use super::*;

#[derive(Debug, Clone)]
pub enum Subset {
    Singleton(Rational),
    Interval(OpenRationalInterval),
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

pub trait Point: Debug + Send + Sync {
    fn rational_interval_neighbourhood(&self) -> Subset;
    fn length(&self) -> Rational;
    fn refine_to_length(&mut self, length: &Rational);
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

    fn refine_to(&self, approx_point: Self::Set, length: &Rational) -> Self::Set {
        approx_point.lock().unwrap().refine_to_length(length);
        approx_point
    }
}

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

    fn refine_to_length(&mut self, _: &Rational) {
        // Nothing needs to be done
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test() {
        let reals = points();

        println!("{:?}", reals);
    }
}
