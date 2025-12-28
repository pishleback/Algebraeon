use crate::{
    approximation::{
        rational_interval::RationalInterval,
        real_intervals::{RealApproximatePointInterface, Subset},
    },
    continued_fraction::{RationalApproximations, SimpleContinuedFraction},
};
use algebraeon_nzq::Rational;
use std::fmt::Debug;

#[derive(Debug, Clone)]
pub struct SimpleContinuedFractionPoint<SCF: SimpleContinuedFraction> {
    n: usize,
    scf: RationalApproximations<SCF>,
}

impl<SCF: SimpleContinuedFraction> From<SCF> for SimpleContinuedFractionPoint<SCF> {
    fn from(scf: SCF) -> Self {
        Self {
            scf: scf.rational_approximations(),
            n: 0,
        }
    }
}

impl<SCF: SimpleContinuedFraction> SimpleContinuedFractionPoint<SCF> {
    pub fn a(&self) -> Rational {
        if self.n.is_multiple_of(2) {
            self.scf.get_rat(self.n)
        } else {
            self.scf.get_rat(self.n + 1)
        }
    }

    pub fn b(&self) -> Rational {
        if self.n.is_multiple_of(2) {
            self.scf.get_rat(self.n + 1)
        } else {
            self.scf.get_rat(self.n)
        }
    }

    pub fn refine(&mut self) {
        self.n += 1;
    }
}

impl<SCF: SimpleContinuedFraction> RealApproximatePointInterface for SimpleContinuedFractionPoint<SCF> {
    fn rational_interval_neighbourhood(&self) -> Subset {
        Subset::Interval(RationalInterval::new_unchecked(self.a(), self.b()))
    }

    fn refine(&mut self) {
        self.refine();
    }
}
