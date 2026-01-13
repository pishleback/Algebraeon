/*
source: https://math.stackexchange.com/questions/2761351/arbitrarily-accurate-rational-approximations-of-pi-with-explicit-error-bound

    \pi = sum_{k=0}^{\infty} r(k)/16^k
    where r(k) = 4/(8k+1) - 2/(8k+4) - 1/(8k+5) - 1/(8k+6)

The series converges quickly to \pi from below and the tail of the series is bounded by
    sum_{k=n+1}^{\infty} r(k)/16^k   <   106 / 16^{n-1}*12285

*/

use algebraeon_nzq::{Integer, Natural, Rational};

use crate::{
    approximation::{
        rational_interval::RationalInterval,
        real_intervals::{RealApproximatePointInterface, Subset},
    },
    structure::MetaMultiplicativeMonoidSignature,
};

fn r(k: Natural) -> Rational {
    let k8 = Integer::from(8) * Integer::from(k);
    Rational::from_integers(Integer::from(4), &k8 + Integer::from(1))
        - Rational::from_integers(Integer::from(2), &k8 + Integer::from(4))
        - Rational::from_integers(Integer::from(1), &k8 + Integer::from(5))
        - Rational::from_integers(Integer::from(1), k8 + Integer::from(6))
}

#[derive(Debug, Clone)]
pub struct Pi {
    n: usize,
    t: Rational, // sum_{k=0}^n r(k)/16^k
}

impl Pi {
    pub fn new() -> Self {
        Self {
            n: 0,
            t: r(Natural::ZERO),
        }
    }

    fn tail_bound(&self) -> Rational {
        if self.n == 0 {
            Rational::from_integers(Integer::from(1), Integer::from(100))
        } else {
            Rational::from_integers(
                Integer::from(106),
                (Integer::from(16).nat_pow(&Natural::from(self.n - 1))) * Integer::from(12285),
            )
        }
    }
}

impl RealApproximatePointInterface for Pi {
    fn rational_interval_neighbourhood(&self) -> Subset {
        Subset::Interval(RationalInterval::new_unchecked(
            self.t.clone(),
            &self.t + self.tail_bound(),
        ))
    }

    fn length(&self) -> Rational {
        self.tail_bound()
    }

    fn refine(&mut self) {
        self.n += 1;
        let k = Natural::from(self.n);
        let d = Rational::from(16).nat_pow(&k);
        self.t += r(k) / d;
    }
}
