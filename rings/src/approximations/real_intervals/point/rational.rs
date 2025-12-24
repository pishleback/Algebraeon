use crate::approximations::real_intervals::{Point, Subset};
use algebraeon_nzq::Rational;
use std::fmt::Debug;

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
