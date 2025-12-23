use algebraeon_nzq::{Rational, RationalCanonicalStructure};

use crate::open_interval::OpenRationalInterval;
use std::{
    fmt::Debug,
    sync::{Arc, Mutex},
};

use super::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RealLineIsolatingIntervalsStructure {}

impl Signature for RealLineIsolatingIntervalsStructure {}

impl SetSignature for RealLineIsolatingIntervalsStructure {
    type Set = OpenRationalInterval;

    fn is_element(&self, _: &Self::Set) -> Result<(), String> {
        Ok(())
    }
}

impl TopologicalSpaceOpenSubsetsSignature for RealLineIsolatingIntervalsStructure {}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RealApproximatePointsStructure {}

pub trait RealApproximatePoint: Debug + Send + Sync {
    fn rational_interval_neighbourhood(&self) -> OpenRationalInterval;
    fn length(&self) -> Rational;
    fn refine_to_length(&mut self, length: &Rational);
}

impl Signature for RealApproximatePointsStructure {}

impl SetSignature for RealApproximatePointsStructure {
    type Set = Arc<Mutex<dyn RealApproximatePoint>>;

    fn is_element(&self, _: &Self::Set) -> Result<(), String> {
        Ok(())
    }
}

impl TopologicalSpaceApproximatePointsSignature for RealApproximatePointsStructure {
    type Precision = RationalCanonicalStructure;
    type OpenSubsetsStructure = RealLineIsolatingIntervalsStructure;

    fn open_neighbourhood(&self, approx_point: &Self::Set) -> OpenRationalInterval {
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
