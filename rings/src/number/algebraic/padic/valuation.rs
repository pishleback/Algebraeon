use super::*;
use algebraeon_nzq::traits::{Abs, DivMod};
use std::ops::{Add, Mul, Sub};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Valuation {
    Infinity,
    Finite(Integer),
}
impl PartialOrd for Valuation {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some({
            match (self, other) {
                (Valuation::Infinity, Valuation::Infinity) => std::cmp::Ordering::Equal,
                (Valuation::Infinity, Valuation::Finite(_)) => std::cmp::Ordering::Greater,
                (Valuation::Finite(_), Valuation::Infinity) => std::cmp::Ordering::Less,
                (Valuation::Finite(finite_self), Valuation::Finite(finite_other)) => {
                    finite_self.cmp(finite_other)
                }
            }
        })
    }
}
impl Ord for Valuation {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other).unwrap()
    }
}
impl Valuation {
    pub fn unwrap_int(self) -> Integer {
        match self {
            Valuation::Infinity => panic!("unwrap_int() called on an infinite value"),
            Valuation::Finite(v) => v,
        }
    }

    pub fn unwrap_nat(self) -> Natural {
        match self {
            Valuation::Infinity => panic!("unwrap_nat() called on an infinite value"),
            Valuation::Finite(v) => {
                if v < Integer::ZERO {
                    panic!("unwrap_nat() called on a negative finite valuation");
                } else {
                    v.abs()
                }
            }
        }
    }
}

impl Add<&Valuation> for Valuation {
    type Output = Valuation;
    fn add(self, other: &Valuation) -> Self::Output {
        match (self, other) {
            (Valuation::Infinity, Valuation::Infinity) => Valuation::Infinity,
            (Valuation::Infinity, Valuation::Finite(_)) => Valuation::Infinity,
            (Valuation::Finite(_), Valuation::Infinity) => Valuation::Infinity,
            (Valuation::Finite(a), Valuation::Finite(b)) => Valuation::Finite(a + b),
        }
    }
}
impl Add<Valuation> for &Valuation {
    type Output = Valuation;
    fn add(self, other: Valuation) -> Self::Output {
        match (self, other) {
            (Valuation::Infinity, Valuation::Infinity) => Valuation::Infinity,
            (Valuation::Infinity, Valuation::Finite(_)) => Valuation::Infinity,
            (Valuation::Finite(_), Valuation::Infinity) => Valuation::Infinity,
            (Valuation::Finite(a), Valuation::Finite(b)) => Valuation::Finite(a + b),
        }
    }
}
impl Add<Valuation> for Valuation {
    type Output = Valuation;
    fn add(self, other: Valuation) -> Self::Output {
        match (self, other) {
            (Valuation::Infinity, Valuation::Infinity) => Valuation::Infinity,
            (Valuation::Infinity, Valuation::Finite(_)) => Valuation::Infinity,
            (Valuation::Finite(_), Valuation::Infinity) => Valuation::Infinity,
            (Valuation::Finite(a), Valuation::Finite(b)) => Valuation::Finite(a + b),
        }
    }
}
impl Add<&Valuation> for &Valuation {
    type Output = Valuation;
    fn add(self, other: &Valuation) -> Self::Output {
        match (self, other) {
            (Valuation::Infinity, Valuation::Infinity) => Valuation::Infinity,
            (Valuation::Infinity, Valuation::Finite(_)) => Valuation::Infinity,
            (Valuation::Finite(_), Valuation::Infinity) => Valuation::Infinity,
            (Valuation::Finite(a), Valuation::Finite(b)) => Valuation::Finite(a + b),
        }
    }
}

/*
Infinite - Infinite is undefined since it coresponds in valuation land to zero/zero
Infinite - Finite is Infinite since it coresponds in valuation land to zero/finite
Finite - Infinite is undefined since it coresponds in valuation land to finite/zero
*/
impl Sub<&Valuation> for Valuation {
    type Output = Valuation;
    fn sub(self, other: &Valuation) -> Self::Output {
        match (self, other) {
            (_, Valuation::Infinity) => {
                panic!("Can't subtract valuation Infinity")
            }
            (Valuation::Infinity, Valuation::Finite(_)) => Valuation::Infinity,
            (Valuation::Finite(a), Valuation::Finite(b)) => Valuation::Finite(a - b),
        }
    }
}
impl Sub<Valuation> for &Valuation {
    type Output = Valuation;
    fn sub(self, other: Valuation) -> Self::Output {
        match (self, other) {
            (_, Valuation::Infinity) => {
                panic!("Can't subtract valuation Infinity")
            }
            (Valuation::Infinity, Valuation::Finite(_)) => Valuation::Infinity,
            (Valuation::Finite(a), Valuation::Finite(b)) => Valuation::Finite(a - b),
        }
    }
}
impl Sub<Valuation> for Valuation {
    type Output = Valuation;
    fn sub(self, other: Valuation) -> Self::Output {
        match (self, other) {
            (_, Valuation::Infinity) => {
                panic!("Can't subtract valuation Infinity")
            }
            (Valuation::Infinity, Valuation::Finite(_)) => Valuation::Infinity,
            (Valuation::Finite(a), Valuation::Finite(b)) => Valuation::Finite(a - b),
        }
    }
}
impl Sub<&Valuation> for &Valuation {
    type Output = Valuation;
    fn sub(self, other: &Valuation) -> Self::Output {
        match (self, other) {
            (_, Valuation::Infinity) => {
                panic!("Can't subtract valuation Infinity")
            }
            (Valuation::Infinity, Valuation::Finite(_)) => Valuation::Infinity,
            (Valuation::Finite(a), Valuation::Finite(b)) => Valuation::Finite(a - b),
        }
    }
}

/*
Finite * Infinity = Infinity is justified by thinking about the p-adic norm instead where non-zero * zero = zero
Infinity * Infinity = Infinity is justified in the same way by the p-adic norm since zero * zero = zero
*/
impl Mul<&Valuation> for Valuation {
    type Output = Valuation;
    fn mul(self, other: &Valuation) -> Self::Output {
        match (self, other) {
            (Valuation::Infinity, Valuation::Infinity) => Valuation::Infinity,
            (Valuation::Infinity, Valuation::Finite(_)) => Valuation::Infinity,
            (Valuation::Finite(_), Valuation::Infinity) => Valuation::Infinity,
            (Valuation::Finite(a), Valuation::Finite(b)) => Valuation::Finite(a * b),
        }
    }
}
impl Mul<Valuation> for &Valuation {
    type Output = Valuation;
    fn mul(self, other: Valuation) -> Self::Output {
        match (self, other) {
            (Valuation::Infinity, Valuation::Infinity) => Valuation::Infinity,
            (Valuation::Infinity, Valuation::Finite(_)) => Valuation::Infinity,
            (Valuation::Finite(_), Valuation::Infinity) => Valuation::Infinity,
            (Valuation::Finite(a), Valuation::Finite(b)) => Valuation::Finite(a * b),
        }
    }
}
impl Mul<Valuation> for Valuation {
    type Output = Valuation;
    fn mul(self, other: Valuation) -> Self::Output {
        match (self, other) {
            (Valuation::Infinity, Valuation::Infinity) => Valuation::Infinity,
            (Valuation::Infinity, Valuation::Finite(_)) => Valuation::Infinity,
            (Valuation::Finite(_), Valuation::Infinity) => Valuation::Infinity,
            (Valuation::Finite(a), Valuation::Finite(b)) => Valuation::Finite(a * b),
        }
    }
}
impl Mul<&Valuation> for &Valuation {
    type Output = Valuation;
    fn mul(self, other: &Valuation) -> Self::Output {
        match (self, other) {
            (Valuation::Infinity, Valuation::Infinity) => Valuation::Infinity,
            (Valuation::Infinity, Valuation::Finite(_)) => Valuation::Infinity,
            (Valuation::Finite(_), Valuation::Infinity) => Valuation::Infinity,
            (Valuation::Finite(a), Valuation::Finite(b)) => Valuation::Finite(a * b),
        }
    }
}

pub fn padic_nat_valuation(p: &Natural, mut n: Natural) -> Valuation {
    debug_assert!(is_prime(p));
    if n == Natural::ZERO {
        Valuation::Infinity
    } else {
        let mut k = 0usize;
        let mut r;
        loop {
            (n, r) = n.div_mod(p);
            if r == Natural::ZERO {
                k += 1;
                continue;
            } else {
                break;
            }
        }
        Valuation::Finite(Integer::from(k))
    }
}

pub fn padic_int_valuation(p: &Natural, n: Integer) -> Valuation {
    padic_nat_valuation(p, n.abs())
}

pub fn padic_rat_valuation(p: &Natural, r: Rational) -> Valuation {
    let (n, d) = r.into_abs_numerator_and_denominator();
    match padic_nat_valuation(p, n) {
        Valuation::Infinity => Valuation::Infinity,
        Valuation::Finite(vn) => Valuation::Finite(vn - padic_nat_valuation(p, d).unwrap_int()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_padic_valuation() {
        assert_eq!(
            padic_int_valuation(&Natural::from(2u32), Integer::from(0)),
            Valuation::Infinity
        );
        assert_eq!(
            padic_int_valuation(&Natural::from(2u32), Integer::from(5)),
            Valuation::Finite(Integer::from(0))
        );
        assert_eq!(
            padic_int_valuation(&Natural::from(2u32), Integer::from(-12)),
            Valuation::Finite(Integer::from(2))
        );
        assert_eq!(
            padic_int_valuation(&Natural::from(2u32), Integer::from(256)),
            Valuation::Finite(Integer::from(8))
        );

        assert_eq!(
            padic_int_valuation(&Natural::from(7u32), Integer::from(0)),
            Valuation::Infinity
        );
        assert_eq!(
            padic_int_valuation(&Natural::from(7u32), Integer::from(-98)),
            Valuation::Finite(Integer::from(2))
        );
        assert_eq!(
            padic_int_valuation(&Natural::from(7u32), Integer::from(42)),
            Valuation::Finite(Integer::from(1))
        );
    }
}
