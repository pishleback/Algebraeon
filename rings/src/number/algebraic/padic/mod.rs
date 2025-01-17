use malachite_base::num::arithmetic::traits::{DivMod, UnsignedAbs};
use malachite_base::num::basic::traits::{One, Two, Zero};
use malachite_nz::{integer::Integer, natural::Natural};
use malachite_q::Rational;

use crate::{number::natural::primes::*, polynomial::polynomial::*, structure::structure::*};

pub mod valuation;
use valuation::*;

mod isolate;

#[derive(Debug, Clone)]
pub struct IsolatingBall {
    p: Natural,
    c: Rational,
    v: Valuation,
}

impl IsolatingBall {
    pub fn overlap(x: &Self, y: &Self) -> bool {
        let p = &x.p;
        debug_assert_eq!(p, &y.p);
        debug_assert!(is_prime(p));
        let vdiff = padic_rat_valuation(p, &x.c - &y.c);
        vdiff >= x.v && vdiff >= y.v
    }
}

#[derive(Debug, Clone)]
pub struct PAdicAlgebraicRoot {
    // A prime number
    p: Natural,
    // An irreducible primitive fav-assoc polynomial of degree >= 2
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

    pub fn refine(&mut self, ndigits: &Integer) {
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

    pub fn isolating_ball(&self) -> IsolatingBall {
        IsolatingBall {
            p: self.p.clone(),
            c: self.approx.center().clone(),
            v: Valuation::Finite(self.approx.ndigits().clone()),
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

    pub fn isolating_ball(&self) -> IsolatingBall {
        IsolatingBall {
            p: self.p.clone(),
            c: self.rat.clone(),
            v: Valuation::Infinity,
        }
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

    pub fn isolating_ball(&self) -> IsolatingBall {
        match self {
            PAdicAlgebraic::Rational(x) => x.isolating_ball(),
            PAdicAlgebraic::Algebraic(x) => x.isolating_ball(),
        }
    }

    pub fn refine(&mut self, ndigits: &Integer) {
        match self {
            PAdicAlgebraic::Rational(_) => {}
            PAdicAlgebraic::Algebraic(x) => x.refine(ndigits),
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

pub mod structure {
    use crate::number::algebraic::poly_tools::root_sum_poly;

    use super::*;
    use algebraeon_sets::structure::*;

    impl PAdicRational {
        fn equal(a: &Self, b: &Self) -> bool {
            debug_assert_eq!(a.p, b.p);
            a.rat == b.rat
        }

        fn add(a: &Self, b: &Self) -> Self {
            let p = &a.p;
            debug_assert_eq!(p, &b.p);
            Self {
                p: p.clone(),
                rat: &a.rat + &b.rat,
            }
        }

        fn neg(self) -> Self {
            Self {
                p: self.p,
                rat: -self.rat,
            }
        }
    }

    impl PAdicAlgebraicRoot {
        fn equal_mut(a: &mut Self, b: &mut Self) -> bool {
            debug_assert_eq!(a.p, b.p);
            if a.poly != b.poly {
                return false;
            }
            let ndigits = Integer::min(a.approx.ndigits().clone(), b.approx.ndigits().clone());
            a.truncate(&ndigits).rational_value() == b.truncate(&ndigits).rational_value()
        }

        fn neg(mut self) -> Self {
            self.poly = Polynomial::compose(
                &self.poly,
                &Polynomial::from_coeffs(vec![Integer::from(0), Integer::from(-1)]),
            );
            self.approx = self.approx.neg();
            self
        }

        fn add_rat(&self, rat: &PAdicRational) -> Self {
            let p = self.p.clone();
            debug_assert_eq!(p, rat.p);
            let poly = Polynomial::compose(
                &self.poly.apply_map(|c| Rational::from(c)),
                &Polynomial::from_coeffs(vec![-&rat.rat, Rational::ONE]),
            )
            .primitive_part_fof();
            let approx = self.approx.clone().add_rat(&rat.rat);
            Self { p, poly, approx }
        }

        fn add_mut(a: &mut Self, b: &mut Self) -> PAdicAlgebraic {
            let p = a.p.clone();
            debug_assert_eq!(p, b.p);
            let mut candidates = root_sum_poly(&a.poly, &b.poly)
                .primitive_squarefree_part()
                .all_padic_roots(&p);
            let mut k = Integer::ZERO;
            while candidates.len() > 1 {
                a.refine(&k);
                b.refine(&k);
                let aball = a.isolating_ball();
                let bball = b.isolating_ball();
                let sball = IsolatingBall {
                    p: p.clone(),
                    c: aball.c + bball.c,
                    v: std::cmp::min(aball.v.clone(), bball.v.clone()),
                };
                candidates = candidates
                    .into_iter()
                    .filter_map(|mut root| {
                        root.refine(&k);
                        let rball = root.isolating_ball();
                        match IsolatingBall::overlap(&rball, &sball) {
                            true => Some(root),
                            false => None,
                        }
                    })
                    .collect();
                k += Integer::ONE;
            }
            debug_assert_eq!(candidates.len(), 1);
            candidates.into_iter().next().unwrap()
        }
    }

    impl PAdicAlgebraic {
        fn neg(self) -> Self {
            match self {
                PAdicAlgebraic::Rational(x) => PAdicAlgebraic::Rational(x.neg()),
                PAdicAlgebraic::Algebraic(x) => PAdicAlgebraic::Algebraic(x.neg()),
            }
        }
    }

    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct PAdicAlgebraicStructure {
        p: Natural,
    }

    impl Structure for PAdicAlgebraicStructure {
        type Set = PAdicAlgebraic;
    }

    impl PAdicAlgebraicStructure {
        pub fn new(p: Natural) -> Self {
            if !is_prime(&p) {
                panic!("{} is not prime", p)
            }
            Self { p }
        }
    }

    impl PAdicAlgebraicStructure {
        fn check_is_element(&self, a: &<Self as Structure>::Set) {
            #[cfg(debug_assertions)]
            if &self.p != a.p() {
                panic!(
                    "{}-adic structure cannot use {}-adic elements",
                    self.p,
                    a.p()
                );
            }
        }
    }

    impl PartialEqStructure for PAdicAlgebraicStructure {
        fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
            self.check_is_element(a);
            self.check_is_element(b);
            match (a, b) {
                (PAdicAlgebraic::Rational(a), PAdicAlgebraic::Rational(b)) => {
                    PAdicRational::equal(a, b)
                }
                (PAdicAlgebraic::Rational(_), PAdicAlgebraic::Algebraic(_)) => false,
                (PAdicAlgebraic::Algebraic(_), PAdicAlgebraic::Rational(_)) => false,
                (PAdicAlgebraic::Algebraic(a), PAdicAlgebraic::Algebraic(b)) => {
                    PAdicAlgebraicRoot::equal_mut(&mut a.clone(), &mut b.clone())
                }
            }
        }
    }

    impl EqStructure for PAdicAlgebraicStructure {}

    impl SemiRingStructure for PAdicAlgebraicStructure {
        fn zero(&self) -> Self::Set {
            PAdicAlgebraic::Rational(PAdicRational {
                p: self.p.clone(),
                rat: Rational::ZERO,
            })
        }

        fn one(&self) -> Self::Set {
            PAdicAlgebraic::Rational(PAdicRational {
                p: self.p.clone(),
                rat: Rational::ONE,
            })
        }

        fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
            self.check_is_element(a);
            self.check_is_element(b);
            match (a, b) {
                (PAdicAlgebraic::Rational(a), PAdicAlgebraic::Rational(b)) => {
                    PAdicAlgebraic::Rational(PAdicRational::add(a, b))
                }
                (PAdicAlgebraic::Rational(a), PAdicAlgebraic::Algebraic(b)) => {
                    PAdicAlgebraic::Algebraic(b.add_rat(a))
                }
                (PAdicAlgebraic::Algebraic(a), PAdicAlgebraic::Rational(b)) => {
                    PAdicAlgebraic::Algebraic(a.add_rat(b))
                }
                (PAdicAlgebraic::Algebraic(a), PAdicAlgebraic::Algebraic(b)) => {
                    PAdicAlgebraicRoot::add_mut(&mut a.clone(), &mut b.clone())
                }
            }
        }

        fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
            self.check_is_element(a);
            self.check_is_element(b);
            match (a, b) {
                (PAdicAlgebraic::Rational(a), PAdicAlgebraic::Rational(b)) => {
                    PAdicAlgebraic::Rational(PAdicRational {
                        p: self.p.clone(),
                        rat: &a.rat * &b.rat,
                    })
                }
                (PAdicAlgebraic::Rational(a), PAdicAlgebraic::Algebraic(b)) => todo!(),
                (PAdicAlgebraic::Algebraic(a), PAdicAlgebraic::Rational(b)) => todo!(),
                (PAdicAlgebraic::Algebraic(a), PAdicAlgebraic::Algebraic(b)) => todo!(),
            }
        }
    }

    impl RingStructure for PAdicAlgebraicStructure {
        fn neg(&self, a: &Self::Set) -> Self::Set {
            self.check_is_element(a);
            a.clone().neg()
        }
    }

    impl IntegralDomainStructure for PAdicAlgebraicStructure {
        fn div(&self, a: &Self::Set, b: &Self::Set) -> Result<Self::Set, RingDivisionError> {
            self.check_is_element(a);
            self.check_is_element(b);
            match (a, b) {
                (PAdicAlgebraic::Rational(a), PAdicAlgebraic::Rational(b)) => {
                    Ok(PAdicAlgebraic::Rational(PAdicRational {
                        p: self.p.clone(),
                        rat: Rational::div(&a.rat, &b.rat)?,
                    }))
                }
                (PAdicAlgebraic::Rational(a), PAdicAlgebraic::Algebraic(b)) => todo!(),
                (PAdicAlgebraic::Algebraic(a), PAdicAlgebraic::Rational(b)) => todo!(),
                (PAdicAlgebraic::Algebraic(a), PAdicAlgebraic::Algebraic(b)) => todo!(),
            }
        }

        fn from_rat(&self, x: &Rational) -> Option<Self::Set> {
            Some(PAdicAlgebraic::Rational(PAdicRational {
                p: self.p.clone(),
                rat: x.clone(),
            }))
        }
    }

    #[cfg(test)]
    mod tests {
        use crate::structure::elements::*;

        use super::*;

        #[test]
        fn test_padic_field_opps() {
            let ring = PAdicAlgebraicStructure::new(Natural::from(5u32));
            let x = Polynomial::<Integer>::var().into_ergonomic();

            let a = {
                let f = (x.pow(3) - 3 * x.pow(2) - x.pow(1) + 1).into_verbose();
                let r = f.all_padic_roots(&Natural::from(5u32));
                assert_eq!(r.len(), 1);
                r.into_iter().next().unwrap()
            };

            let b = {
                let f = (x.pow(4) + x.pow(2) - 2 * x.pow(1) - 1).into_verbose();
                let r = f.all_padic_roots(&Natural::from(5u32));
                assert_eq!(r.len(), 1);
                r.into_iter().next().unwrap()
            };

            let c = {
                let f = (x.pow(5) + x.pow(2) + 2 * x.pow(1) + 1).into_verbose();
                let r = f.all_padic_roots(&Natural::from(5u32));
                assert_eq!(r.len(), 1);
                r.into_iter().next().unwrap()
            };

            let d = PAdicAlgebraic::from_rational(
                Natural::from(5u32),
                Rational::from_integers(Integer::from(2), Integer::from(7)),
            );

            println!("a = {}", a);
            println!("b = {}", b);
            println!("c = {}", c);
            println!("d = {}", d);

            println!("-a = {}", ring.neg(&a));
            println!("-b = {}", ring.neg(&b));
            println!("-c = {}", ring.neg(&c));
            println!("-d = {}", ring.neg(&d));

            println!("a+b = {}", ring.add(&a, &b));
            println!("a+c = {}", ring.add(&a, &c));
            println!("d+b = {}", ring.add(&d, &b));
            println!("d+c = {}", ring.add(&d, &c));
            println!("c+c = {}", ring.add(&c, &c));

            /*

            let a = {
                let f = (x.pow(3) - 3 * x.pow(2) - x.pow(1) + 1).into_verbose();
                let r = f.all_padic_roots(&Natural::from(5u32));
                assert_eq!(r.len(), 1);
                r.into_iter().next().unwrap().shift_by(-1)
            };

            let b = {
                let f = (x.pow(4) + x.pow(2) - 2 * x.pow(1) - 1).into_verbose();
                let r = f.all_padic_roots(&Natural::from(5u32));
                assert_eq!(r.len(), 1);
                r.into_iter().next().unwrap().shift_by(-1)
            };

            let c = {
                let f = (x.pow(5) + x.pow(2) + 2 * x.pow(1) + 1).into_verbose();
                let r = f.all_padic_roots(&Natural::from(5u32));
                assert_eq!(r.len(), 1);
                r.into_iter().next().unwrap().shift_by(-1)
            };

            let d = a.clone().shift_by(4);

            let x = ring
                .from_rat(&Rational::from_integers(
                    Integer::from(2),
                    Integer::from(3 * 125),
                ))
                .unwrap();

            println!("a = {}", a);
            println!("b = {}", b);
            println!("c = {}", c);
            println!("d = {}", d);
            println!("x = {}", x);

            println!("-a = {}", ring.neg(&a));
            debug_assert_eq!(
                ring.neg(&a).reduce_modulo_valuation(5).digits(),
                (
                    vec![
                        Natural::from(3u8),
                        Natural::from(0u8),
                        Natural::from(2u8),
                        Natural::from(3u8),
                        Natural::from(3u8),
                        Natural::from(1u8)
                    ],
                    -1
                )
            );
            debug_assert_eq!(a.valuation(), Some(-1));

            println!("-b = {}", ring.neg(&b));
            debug_assert_eq!(
                ring.neg(&b).reduce_modulo_valuation(5).digits(),
                (
                    vec![
                        Natural::from(3u8),
                        Natural::from(1u8),
                        Natural::from(3u8),
                        Natural::from(2u8),
                        Natural::from(4u8),
                        Natural::from(3u8)
                    ],
                    -1
                )
            );

            println!("-x = {}", ring.neg(&x));
            debug_assert_eq!(
                ring.neg(&x).reduce_modulo_valuation(3).digits(),
                (
                    vec![
                        Natural::from(1u8),
                        Natural::from(3u8),
                        Natural::from(1u8),
                        Natural::from(3u8),
                        Natural::from(1u8),
                        Natural::from(3u8),
                    ],
                    -3
                )
            );

            println!("a+x = {}", ring.add(&a, &x));
            debug_assert_eq!(
                ring.add(&a, &x).reduce_modulo_valuation(6).digits(),
                (
                    vec![
                        Natural::from(4u8),
                        Natural::from(1u8),
                        Natural::from(0u8),
                        Natural::from(1u8),
                        Natural::from(1u8),
                        Natural::from(3u8),
                        Natural::from(4u8),
                        Natural::from(4u8),
                        Natural::from(1u8),
                    ],
                    -3
                )
            );

            println!("b+x = {}", ring.add(&b, &x));
            debug_assert_eq!(
                ring.add(&b, &x).reduce_modulo_valuation(6).digits(),
                (
                    vec![
                        Natural::from(4u8),
                        Natural::from(1u8),
                        Natural::from(0u8),
                        Natural::from(0u8),
                        Natural::from(0u8),
                        Natural::from(4u8),
                        Natural::from(3u8),
                        Natural::from(2u8),
                        Natural::from(2u8),
                    ],
                    -3
                )
            );

            println!("c+x = {}", ring.add(&c, &x));
            debug_assert_eq!(
                ring.add(&c, &x).reduce_modulo_valuation(6).digits(),
                (
                    vec![
                        Natural::from(4u8),
                        Natural::from(1u8),
                        Natural::from(4u8),
                        Natural::from(2u8),
                        Natural::from(1u8),
                        Natural::from(1u8),
                        Natural::from(0u8),
                        Natural::from(2u8),
                        Natural::from(0u8),
                    ],
                    -3
                )
            );

            println!("a+b = {}", ring.add(&a, &b));
            debug_assert_eq!(
                ring.add(&a, &b).reduce_modulo_valuation(6).digits(),
                (
                    vec![
                        Natural::from(4u8),
                        Natural::from(2u8),
                        Natural::from(4u8),
                        Natural::from(3u8),
                        Natural::from(1u8),
                        Natural::from(4u8),
                        Natural::from(2u8),
                    ],
                    -1
                )
            );

            println!("c+d = {}", ring.add(&c, &d));
            // debug_assert_eq!(
            //     ring.add(&a, &b).reduce_modulo_valuation(6).digits(),
            //     (
            //         vec![
            //             Natural::from(4u8),
            //             Natural::from(2u8),
            //             Natural::from(4u8),
            //             Natural::from(3u8),
            //             Natural::from(1u8),
            //             Natural::from(4u8),
            //             Natural::from(2u8),
            //         ],
            //         -1
            //     )
            // );
            */
        }
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
