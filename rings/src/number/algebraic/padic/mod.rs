use malachite_base::num::arithmetic::traits::{DivMod, UnsignedAbs};
use malachite_base::num::basic::traits::{One, Two, Zero};
use malachite_nz::{integer::Integer, natural::Natural};
use malachite_q::Rational;

use crate::{number::natural::primes::*, polynomial::polynomial::*, structure::structure::*};

// Some algorithms here on p-adic root isolation can be found in
// Sturm, Thomas & Weispfenning, Volker. (2004). P-adic Root Isolation. Revista de la Real Academia de Ciencias Exactas, Físicas y Naturales. Serie A, Matemáticas.
// https://www.researchgate.net/profile/Thomas-Sturm-2/publication/2925550_P-adic_Root_Isolation/links/580b8c0708aeef1bfeeb5db8/P-adic-Root-Isolation.pdf?origin=scientificContributions

pub mod valuation;
use valuation::*;

mod isolate;

#[derive(Debug, Clone)]
enum PAdicAlgebraicRootHenselLiftable {
    No,
    Yes {
        // f'(a) where a is the approximate root OR equivelently the lifted root.
        dpoly_valuation: Integer,
    },
}

#[derive(Debug, Clone)]
pub struct PAdicAlgebraicRoot {
    // A prime number
    p: Natural,
    // An irreducible polynomial of degree >= 2
    poly: Polynomial<Integer>,
    // A p-adic isolating ball containing exactly one root of the polynomial
    approx: isolate::PAdicRationalBall,
}

impl PAdicAlgebraicRoot {
    fn new(p: Natural, poly: Polynomial<Integer>, approx: isolate::PAdicRationalBall) -> Self {
        debug_assert!(is_prime(&p));
        debug_assert!(poly.is_irreducible());
        Self { p, poly, approx }
    }

    fn refine(&mut self, ndigits: &Integer) {
        while self.approx.ndigits() < ndigits {
            self.approx = isolate::refine(
                &self.p,
                &self.poly,
                &self.approx,
                self.approx.ndigits() + Integer::ONE,
            );
            // verify that the root was lifted correctly
            debug_assert_eq!(
                PAdicRational::from_rational(
                    self.p.clone(),
                    self.poly
                        .apply_map(|coeff| Rational::from(coeff))
                        .evaluate(self.approx.center())
                )
                .truncate(self.approx.ndigits())
                .rational_value(),
                Rational::ZERO
            );
        }
    }
}

#[derive(Debug, Clone)]
pub struct PAdicRational {
    // A prime number
    p: Natural,
    rat: Rational,
}

impl PAdicRational {
    pub fn from_rational(p: Natural, rat: Rational) -> Self {
        debug_assert!(is_prime(&p));
        Self { p, rat }
    }
}

#[derive(Debug, Clone)]
pub enum PAdicAlgebraic {
    Rational(PAdicRational),
    Algebraic(PAdicAlgebraicRoot),
}

impl PAdicAlgebraic {
    pub fn from_rational(p: Natural, rat: Rational) -> Self {
        Self::Rational(PAdicRational::from_rational(p, rat))
    }

    pub fn p(&self) -> &Natural {
        match self {
            PAdicAlgebraic::Rational(x) => &x.p,
            PAdicAlgebraic::Algebraic(x) => &x.p,
        }
    }
}

impl std::fmt::Display for PAdicAlgebraic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let p = self.p();
        let n = Integer::from(if p < &Natural::from(10u32) { 6 } else { 3 });
        write!(f, "{}", self.clone().truncate(&n).string_repr())
    }
}

pub mod truncation {
    use super::*;

    // Represent value * p^shift with 0 <= value < p^digits
    #[derive(Debug, Clone)]
    pub enum Truncated {
        Zero {
            p: Natural,
        },
        NonZero {
            p: Natural,
            value: Natural, // non-zero mod p i.e. valuation 0
            shift: Integer,
            num_digits: Natural, // >= 1
        },
    }

    impl Truncated {
        pub fn digits(&self) -> Option<(Vec<Natural>, Integer)> {
            match self {
                Truncated::Zero { .. } => None,
                Truncated::NonZero {
                    p,
                    value,
                    shift,
                    num_digits,
                } => Some({
                    let mut k = Natural::ZERO;
                    let mut digits = vec![];
                    let mut v = value.clone();
                    while &k < num_digits {
                        let r;
                        (v, r) = v.div_mod(p);
                        digits.push(r);
                        k += Natural::ONE;
                    }
                    (digits, shift.clone())
                }),
            }
        }

        pub fn rational_value(&self) -> Rational {
            match self {
                Truncated::Zero { .. } => Rational::ZERO,
                Truncated::NonZero {
                    p, value, shift, ..
                } => Rational::from(value) * Rational::from(p).int_pow(shift).unwrap(),
            }
        }

        pub fn string_repr(&self) -> String {
            let p = match self {
                Truncated::Zero { p } => p,
                Truncated::NonZero { p, .. } => p,
            };
            match self.digits() {
                None => "0".into(),
                Some((digits, mut shift)) => {
                    use std::fmt::Write;
                    let seps = p >= &Natural::from(10u32);
                    let mut rev_digits = digits.into_iter().rev().collect::<Vec<_>>();
                    while shift > 0 {
                        rev_digits.push(Natural::ZERO);
                        shift -= Integer::ONE;
                    }
                    debug_assert!(shift <= 0);
                    let shift = (-shift).unsigned_abs();
                    let mut s = String::new();
                    write!(&mut s, "...").unwrap();
                    for (i, d) in rev_digits.into_iter().rev().enumerate().rev() {
                        write!(&mut s, "{}", d).unwrap();
                        if i != 0 {
                            if seps {
                                if i == shift {
                                    write!(&mut s, ";").unwrap();
                                } else {
                                    write!(&mut s, ",").unwrap();
                                }
                            } else {
                                if i == shift {
                                    write!(&mut s, ".").unwrap();
                                }
                            }
                        }
                    }
                    s
                }
            }
        }
    }

    impl PAdicRational {
        pub fn truncate(&self, cutoffv: &Integer) -> Truncated {
            match padic_rat_valuation(&self.p, self.rat.clone()) {
                Valuation::Finite(shift) => {
                    let shifted_rat =
                        &self.rat * Rational::from(&self.p).int_pow(&-&shift).unwrap();
                    let (n, d) = (shifted_rat.numerator(), shifted_rat.denominator());
                    debug_assert_eq!(padic_int_valuation(&self.p, n.clone()).unwrap_nat(), 0);
                    debug_assert_eq!(padic_int_valuation(&self.p, d.clone()).unwrap_nat(), 0);
                    if cutoffv <= &shift {
                        Truncated::Zero { p: self.p.clone() }
                    } else {
                        let num_digits = cutoffv - &shift;
                        debug_assert!(num_digits >= 1);
                        let num_digits = num_digits.unsigned_abs();
                        let pn = Integer::from(&self.p).nat_pow(&num_digits); // p^{num_digits}
                        let (g, _, d_inv) = Integer::xgcd(&pn, &d);
                        debug_assert_eq!(g, Integer::ONE);
                        let value = Integer::rem(&(n * d_inv), &pn);
                        debug_assert!(value > 0);
                        let value = value.unsigned_abs();
                        Truncated::NonZero {
                            p: self.p.clone(),
                            value,
                            shift,
                            num_digits,
                        }
                    }
                }
                Valuation::Infinity => Truncated::Zero { p: self.p.clone() },
            }
        }
    }

    impl PAdicAlgebraicRoot {
        pub fn truncate(&mut self, modulus_pow: &Integer) -> Truncated {
            self.refine(modulus_pow);
            let rat = self.approx.center();
            PAdicRational {
                p: self.p.clone(),
                rat: rat.clone(),
            }
            .truncate(modulus_pow)
        }
    }

    impl PAdicAlgebraic {
        // Truncate modulo p^cutoffv
        // e.g.
        //  cutoffv=0 is modulo 1
        //  cutoffv=1 is modulo p
        pub fn truncate(&mut self, cutoffv: &Integer) -> Truncated {
            match self {
                PAdicAlgebraic::Rational(a) => a.truncate(cutoffv),
                PAdicAlgebraic::Algebraic(a) => a.truncate(cutoffv),
            }
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn test_padic_digits() {
            assert_eq!(
                PAdicAlgebraic::from_rational(Natural::from(2u32), Rational::ZERO,)
                    .truncate(&6.into())
                    .digits(),
                None
            );

            assert_eq!(
                PAdicAlgebraic::from_rational(
                    Natural::from(2u32),
                    Rational::from_integers(Integer::from(9), Integer::from(1)),
                )
                .truncate(&0.into())
                .digits(),
                None
            );

            assert_eq!(
                PAdicAlgebraic::from_rational(
                    Natural::from(2u32),
                    Rational::from_integers(Integer::from(9), Integer::from(1)),
                )
                .truncate(&6.into())
                .digits(),
                Some((
                    vec![
                        Natural::from(1u32),
                        Natural::from(0u32),
                        Natural::from(0u32),
                        Natural::from(1u32),
                        Natural::from(0u32),
                        Natural::from(0u32),
                    ],
                    Integer::from(0)
                ))
            );

            assert_eq!(
                PAdicAlgebraic::from_rational(
                    Natural::from(2u32),
                    Rational::from_integers(Integer::from(10), Integer::from(1)),
                )
                .truncate(&6.into())
                .digits(),
                Some((
                    vec![
                        Natural::from(1u32),
                        Natural::from(0u32),
                        Natural::from(1u32),
                        Natural::from(0u32),
                        Natural::from(0u32),
                    ],
                    Integer::from(1)
                ))
            );

            assert_eq!(
                PAdicAlgebraic::from_rational(
                    Natural::from(2u32),
                    Rational::from_integers(Integer::from(31), Integer::from(36)),
                )
                .truncate(&10.into())
                .digits(),
                Some((
                    vec![
                        Natural::from(1u32),
                        Natural::from(1u32),
                        Natural::from(1u32),
                        Natural::from(0u32),
                        Natural::from(0u32),
                        Natural::from(1u32),
                        Natural::from(1u32),
                        Natural::from(1u32),
                        Natural::from(0u32),
                        Natural::from(0u32),
                        Natural::from(0u32),
                        Natural::from(1u32),
                    ],
                    Integer::from(-2)
                ))
            );
        }
    }
}

impl Polynomial<Integer> {
    fn all_padic_roots_irreducible(&self, p: &Natural) -> Vec<PAdicAlgebraic> {
        debug_assert!(self.is_irreducible());
        let n = self.degree().unwrap();
        if n == 1 {
            // Rational root
            let a = self.coeff(1);
            let b = self.coeff(0);
            // f(x) = ax + b
            // root = -b/a
            vec![PAdicAlgebraic::Rational(PAdicRational {
                p: p.clone(),
                rat: -Rational::from_integers(b, a),
            })]
        } else {
            isolate::isolate(p, self)
                .into_iter()
                .map(|root| {
                    PAdicAlgebraic::Algebraic(PAdicAlgebraicRoot::new(
                        p.clone(),
                        self.clone(),
                        root,
                    ))
                })
                .collect()
        }
    }

    pub fn all_padic_roots(&self, p: &Natural) -> Vec<PAdicAlgebraic> {
        debug_assert!(is_prime(p));
        assert_ne!(self, &Self::zero());
        let factors = self.factor().unwrap();
        let mut roots = vec![];
        for (factor, k) in factors.factors() {
            for root in factor.all_padic_roots_irreducible(p) {
                let mut i = Natural::from(0u8);
                while &i < k {
                    roots.push(root.clone());
                    i += Natural::from(1u8);
                }
            }
        }
        roots
    }
}

#[cfg(test)]
mod tests {
    use crate::structure::elements::*;

    use super::*;

    #[test]
    fn test_all_padic_roots_example1() {
        let x = Polynomial::<Integer>::var().into_ergonomic();
        let f = (16 * x.pow(2) - 17).into_verbose();
        println!("f = {}", f);
        let roots = f.all_padic_roots(&Natural::from(2u32));
        assert_eq!(roots.len(), 2);
        for root in roots {
            println!("{}", root);
        }
    }

    #[test]
    fn test_all_padic_roots_example2() {
        let x = Polynomial::<Integer>::var().into_ergonomic();
        let f = (x.pow(6) - 2).into_verbose();
        println!("f = {}", f);
        assert_eq!(f.all_padic_roots(&Natural::from(2u32)).len(), 0);
        assert_eq!(f.all_padic_roots(&Natural::from(7u32)).len(), 0);
        assert_eq!(f.all_padic_roots(&Natural::from(727u32)).len(), 6);
        for root in f.all_padic_roots(&Natural::from(727u32)) {
            println!("{}", root);
        }
    }

    #[test]
    fn test_all_padic_roots_example3() {
        let x = Polynomial::<Integer>::var().into_ergonomic();
        let f = (3 * x - 2).into_verbose();
        println!("{:?}", f);
        assert_eq!(f.all_padic_roots(&Natural::from(2u32)).len(), 1);
        assert_eq!(f.all_padic_roots(&Natural::from(3u32)).len(), 1);
        assert_eq!(f.all_padic_roots(&Natural::from(5u32)).len(), 1);
    }
}
