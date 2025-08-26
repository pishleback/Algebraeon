use crate::{
    natural::factorization::primes::is_prime,
    structure::{MetaFactorableSignature, QuotientStructure, SemiRingSignature},
};
use algebraeon_nzq::{Integer, Natural, traits::Abs};
use algebraeon_sets::structure::MetaType;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QuadraticSymbolValue {
    Zero,
    Pos,
    Neg,
}

impl std::ops::Mul for QuadraticSymbolValue {
    type Output = Self;

    fn mul(self, other: Self) -> Self::Output {
        match (self, other) {
            (Self::Zero, _) | (_, Self::Zero) => Self::Zero,
            (Self::Pos, Self::Pos) | (Self::Neg, Self::Neg) => Self::Pos,
            (Self::Pos, Self::Neg) | (Self::Neg, Self::Pos) => Self::Neg,
        }
    }
}

impl QuadraticSymbolValue {
    pub fn nat_pow(&self, k: &Natural) -> Self {
        match self {
            Self::Zero => Self::Zero,
            Self::Pos => Self::Pos,
            Self::Neg => {
                if k % Natural::TWO == Natural::ZERO {
                    Self::Pos
                } else {
                    Self::Neg
                }
            }
        }
    }
}

impl TryFrom<&Integer> for QuadraticSymbolValue {
    type Error = ();
    fn try_from(value: &Integer) -> Result<Self, Self::Error> {
        if *value == Integer::ZERO {
            Ok(Self::Zero)
        } else if *value == Integer::ONE {
            Ok(Self::Pos)
        } else if *value == -Integer::ONE {
            Ok(Self::Neg)
        } else {
            Err(())
        }
    }
}

impl TryFrom<Integer> for QuadraticSymbolValue {
    type Error = ();
    fn try_from(value: Integer) -> Result<Self, Self::Error> {
        Self::try_from(&value)
    }
}

impl Into<Integer> for QuadraticSymbolValue {
    fn into(self) -> Integer {
        match self {
            QuadraticSymbolValue::Zero => Integer::ZERO,
            QuadraticSymbolValue::Pos => Integer::ONE,
            QuadraticSymbolValue::Neg => -Integer::ONE,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LegendreSymbolError {
    BottomNotOddPrime,
}

pub fn legendre_symbol(
    a: &Integer,
    p: &Natural,
) -> Result<QuadraticSymbolValue, LegendreSymbolError> {
    if p % Natural::TWO == Natural::ZERO || !is_prime(p) {
        Err(LegendreSymbolError::BottomNotOddPrime)
    } else {
        let mod_p = QuotientStructure::new_field_unchecked(Integer::structure(), Integer::from(p));
        let v = mod_p.reduce(mod_p.nat_pow(a, &((p - Natural::ONE) / Natural::TWO)));
        if v == Integer::ZERO {
            Ok(QuadraticSymbolValue::Zero)
        } else if v == Integer::ONE {
            Ok(QuadraticSymbolValue::Pos)
        } else {
            debug_assert_eq!(v, Integer::from(p) - Integer::ONE);
            Ok(QuadraticSymbolValue::Neg)
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum JacobiSymbolError {
    BottomEven,
}

// Jacobi symbol computation without factorization
// Based on the algorithm from: https://en.wikipedia.org/wiki/Jacobi_symbol#Implementation_in_C++
pub fn jacobi_symbol(a: &Integer, n: &Natural) -> Result<QuadraticSymbolValue, JacobiSymbolError> {
    if n % Natural::TWO == Natural::ZERO {
        return Err(JacobiSymbolError::BottomEven);
    }

    if n == &Natural::ONE {
        return Ok(QuadraticSymbolValue::Pos);
    }

    let mut a = a.clone();
    let mut n = n.clone();
    let mut result = QuadraticSymbolValue::Pos;

    // Handle negative numbers: (-1/n) = (-1)^((n-1)/2)
    if a < Integer::ZERO {
        a = -a;
        if &n % Natural::from(4u32) == Natural::from(3u32) {
            result = result * QuadraticSymbolValue::Neg;
        }
    }

    loop {
        a = a % Integer::from(&n);

        if a == Integer::ZERO {
            return Ok(QuadraticSymbolValue::Zero);
        }

        // Extract powers of 2 from a
        while &a % Integer::TWO == Integer::ZERO {
            a = &a / Integer::TWO;

            let n_mod_8 = &n % Natural::from(8u32);
            if n_mod_8 == Natural::from(3u32) || n_mod_8 == Natural::from(5u32) {
                result = result * QuadraticSymbolValue::Neg;
            }
        }

        if a == Integer::ONE {
            return Ok(result);
        }

        let a_nat = Natural::try_from(&a).unwrap();

        if &a_nat % Natural::from(4u32) == Natural::from(3u32)
            && &n % Natural::from(4u32) == Natural::from(3u32)
        {
            result = result * QuadraticSymbolValue::Neg;
        }

        // Swap a and n for next iteration
        let temp = a_nat;
        a = Integer::from(&n);
        n = temp;
    }
}

pub fn kronecker_symbol(a: &Integer, n: &Integer) -> QuadraticSymbolValue {
    match n.factor() {
        None => {
            // n == 0
            if a.abs() == Natural::ONE {
                // (1/0) = (-1/0) = 1
                QuadraticSymbolValue::Pos
            } else {
                // (a/0) = 0 for all a != +/-1
                QuadraticSymbolValue::Zero
            }
        }
        Some(n) => {
            let (u, powers) = n.into_unit_and_powers();
            let mut val = if u == Integer::ONE {
                QuadraticSymbolValue::Pos
            } else {
                debug_assert_eq!(u, -Integer::ONE);
                if a >= &Integer::ZERO {
                    QuadraticSymbolValue::Pos
                } else {
                    QuadraticSymbolValue::Neg
                }
            };
            for (p, k) in powers {
                val = val * {
                    if p == Natural::TWO {
                        if a % Integer::TWO == Integer::ZERO {
                            QuadraticSymbolValue::Zero
                        } else {
                            let a_mod_8 = a % Integer::from(8);
                            if a_mod_8 == Integer::from(1) || a_mod_8 == Integer::from(7) {
                                QuadraticSymbolValue::Pos
                            } else {
                                debug_assert!(
                                    a_mod_8 == Integer::from(3) || a_mod_8 == Integer::from(5)
                                );
                                QuadraticSymbolValue::Neg
                            }
                        }
                    } else {
                        legendre_symbol(a, &p.abs()).unwrap()
                    }
                }
                .nat_pow(&k);
            }
            val
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_legendre_symbol() {
        assert_eq!(
            legendre_symbol(&0.into(), &0u32.into()),
            Err(LegendreSymbolError::BottomNotOddPrime)
        );
        assert_eq!(
            legendre_symbol(&0.into(), &1u32.into()),
            Err(LegendreSymbolError::BottomNotOddPrime)
        );
        assert_eq!(
            legendre_symbol(&0.into(), &2u32.into()),
            Err(LegendreSymbolError::BottomNotOddPrime)
        );
        assert_eq!(
            legendre_symbol(&0.into(), &3u32.into()),
            Ok(QuadraticSymbolValue::Zero)
        );
        assert_eq!(
            legendre_symbol(&0.into(), &4u32.into()),
            Err(LegendreSymbolError::BottomNotOddPrime)
        );

        assert_eq!(
            legendre_symbol(&1.into(), &0u32.into()),
            Err(LegendreSymbolError::BottomNotOddPrime)
        );
        assert_eq!(
            legendre_symbol(&1.into(), &1u32.into()),
            Err(LegendreSymbolError::BottomNotOddPrime)
        );
        assert_eq!(
            legendre_symbol(&1.into(), &2u32.into()),
            Err(LegendreSymbolError::BottomNotOddPrime)
        );
        assert_eq!(
            legendre_symbol(&1.into(), &3u32.into()),
            Ok(QuadraticSymbolValue::Pos)
        );
        assert_eq!(
            legendre_symbol(&1.into(), &4u32.into()),
            Err(LegendreSymbolError::BottomNotOddPrime)
        );

        assert_eq!(
            legendre_symbol(&0.into(), &5u32.into()),
            Ok(QuadraticSymbolValue::Zero)
        );
        assert_eq!(
            legendre_symbol(&1.into(), &5u32.into()),
            Ok(QuadraticSymbolValue::Pos)
        );
        assert_eq!(
            legendre_symbol(&2.into(), &5u32.into()),
            Ok(QuadraticSymbolValue::Neg)
        );
        assert_eq!(
            legendre_symbol(&3.into(), &5u32.into()),
            Ok(QuadraticSymbolValue::Neg)
        );
        assert_eq!(
            legendre_symbol(&4.into(), &5u32.into()),
            Ok(QuadraticSymbolValue::Pos)
        );

        assert_eq!(
            legendre_symbol(&0.into(), &7u32.into()),
            Ok(QuadraticSymbolValue::Zero)
        );
        assert_eq!(
            legendre_symbol(&1u32.into(), &7u32.into()),
            Ok(QuadraticSymbolValue::Pos)
        );
        assert_eq!(
            legendre_symbol(&2.into(), &7u32.into()),
            Ok(QuadraticSymbolValue::Pos)
        );
        assert_eq!(
            legendre_symbol(&3.into(), &7u32.into()),
            Ok(QuadraticSymbolValue::Neg)
        );
        assert_eq!(
            legendre_symbol(&4.into(), &7u32.into()),
            Ok(QuadraticSymbolValue::Pos)
        );
        assert_eq!(
            legendre_symbol(&5.into(), &7u32.into()),
            Ok(QuadraticSymbolValue::Neg)
        );
        assert_eq!(
            legendre_symbol(&6.into(), &7u32.into()),
            Ok(QuadraticSymbolValue::Neg)
        );
    }

    #[test]
    fn test_jacobi_symbol() {
        assert_eq!(
            jacobi_symbol(&0.into(), &0u32.into()),
            Err(JacobiSymbolError::BottomEven)
        );
        assert_eq!(
            jacobi_symbol(&0.into(), &1u32.into()),
            Ok(QuadraticSymbolValue::Pos)
        );
        assert_eq!(
            jacobi_symbol(&0.into(), &2u32.into()),
            Err(JacobiSymbolError::BottomEven)
        );
        assert_eq!(
            jacobi_symbol(&0.into(), &6u32.into()),
            Err(JacobiSymbolError::BottomEven)
        );
        assert_eq!(
            jacobi_symbol(&3.into(), &9u32.into()),
            Ok(QuadraticSymbolValue::Zero)
        );
        for p in [3u32, 5, 7] {
            for a in 0..10 {
                let (a, p) = (Integer::from(a), Natural::from(p));
                let ls = legendre_symbol(&a, &p);
                let js = jacobi_symbol(&a, &p);
                assert_eq!(ls.unwrap(), js.unwrap());
            }
        }
    }

    #[test]
    fn test_kronecker_symbol() {
        let vals: Vec<Vec<i32>> = vec![
            vec![0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0],
            vec![
                -1, 0, -1, 1, 1, -1, 0, 1, -1, -1, 1, 0, 1, -1, -1, 1, 0, 1, -1,
            ],
            vec![0, -1, 0, -1, 0, -1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
            vec![
                0, -1, 1, 0, -1, 1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0,
            ],
            vec![0, 1, 0, 1, 0, -1, 0, 1, 0, -1, 0, -1, 0, 1, 0, 1, 0, -1, 0],
            vec![
                -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            ],
            vec![0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            vec![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            vec![0, -1, 0, -1, 0, 1, 0, 1, 0, -1, 0, -1, 0, 1, 0, 1, 0, -1, 0],
            vec![
                0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0,
            ],
            vec![0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
            vec![
                1, 0, 1, -1, -1, 1, 0, 1, -1, -1, 1, 0, 1, -1, -1, 1, 0, 1, -1,
            ],
            vec![0, -1, 0, 0, 0, -1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0],
            vec![1, 1, -1, 1, -1, -1, 0, 1, 1, -1, 1, -1, -1, 0, 1, 1, -1, 1],
            vec![0, -1, 0, -1, 0, 1, 0, 1, 0, -1, 0, -1, 0, 1, 0, 1, 0, -1, 0],
            vec![0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0],
            vec![0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, -1, 0, 1, 0, -1, 0],
            vec![
                1, -1, -1, -1, 1, -1, 0, 1, -1, 1, 1, 1, -1, -1, -1, 1, -1, 0, 1,
            ],
            vec![0, 1, 0, 0, 0, -1, 0, 1, 0, 0, 0, -1, 0, 1, 0, 0, 0, -1, 0],
        ];
        for (n, row) in vals.into_iter().enumerate() {
            for (a, val) in row.iter().enumerate() {
                if val == &2 {
                    continue;
                }
                let (a, n) = (Integer::from(a as isize - 6), Integer::from(n as isize - 6));
                let v: Integer = kronecker_symbol(&a, &n).into();
                println!("({}/{})={}", a, n, v);
                assert_eq!(v, Integer::from(*val));
            }
        }
    }
}
