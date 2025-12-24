use crate::{
    approximation::{
        rational_interval::RationalInterval,
        real_intervals::{Point, Subset},
    },
    continued_fraction::SimpleContinuedFraction,
};
use algebraeon_nzq::{Integer, Rational};
use std::fmt::Debug;

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
        Subset::Interval(RationalInterval::new_unchecked(self.a(), self.b()))
    }

    fn refine(&mut self) {
        self.refine();
    }
}
