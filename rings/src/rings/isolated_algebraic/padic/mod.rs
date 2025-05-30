use crate::{
    polynomial::*,
    rings::{natural::factorization::primes::*, valuation::*},
    structure::*,
};
use algebraeon_nzq::*;
use algebraeon_sets::structure::MetaType;
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
    use algebraeon_nzq::traits::{Abs, DivMod, Fraction};

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
                Truncated::Zero { p } | Truncated::NonZero { p, .. } => p,
            };
            match self.digits() {
                None => "0".into(),
                Some((digits, mut shift)) => {
                    use std::fmt::Write;
                    let seps = p >= &Natural::from(10u32);
                    let mut rev_digits = digits.into_iter().rev().collect::<Vec<_>>();
                    while shift > Integer::ZERO {
                        rev_digits.push(Natural::ZERO);
                        shift -= Integer::ONE;
                    }
                    debug_assert!(shift <= Integer::ZERO);
                    let shift = (-shift).abs();
                    let mut s = String::new();
                    write!(&mut s, "...").unwrap();
                    for (i, d) in rev_digits.into_iter().rev().enumerate().rev() {
                        write!(&mut s, "{}", d).unwrap();
                        #[allow(clippy::collapsible_else_if)]
                        if i != 0 {
                            if seps {
                                if Integer::from(i) == shift {
                                    write!(&mut s, ";").unwrap();
                                } else {
                                    write!(&mut s, ",").unwrap();
                                }
                            } else {
                                if Integer::from(i) == shift {
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
                    let (n, d) = shifted_rat.numerator_and_denominator();
                    let d = Integer::from(d);
                    debug_assert_eq!(
                        padic_int_valuation(&self.p, n.clone()).unwrap_nat(),
                        Natural::ZERO
                    );
                    debug_assert_eq!(
                        padic_int_valuation(&self.p, d.clone()).unwrap_nat(),
                        Natural::ZERO
                    );
                    if cutoffv <= &shift {
                        Truncated::Zero { p: self.p.clone() }
                    } else {
                        let num_digits = cutoffv - &shift;
                        debug_assert!(num_digits >= Integer::ONE);
                        let num_digits = num_digits.abs();
                        let pn = Integer::from(&self.p).nat_pow(&num_digits); // p^{num_digits}
                        let (g, _, d_inv) = Integer::xgcd(&pn, &d);
                        debug_assert_eq!(g, Integer::ONE);
                        let value = Integer::rem(&(n * d_inv), &pn);
                        debug_assert!(value > Integer::ZERO);
                        let value = value.abs();
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
        for (factor, k) in Polynomial::<Integer>::structure()
            .factorizations()
            .to_powers(&factors)
        {
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
    use crate::rings::isolated_algebraic::poly_tools::{
        root_product_poly, root_rat_mul_poly, root_sum_poly,
    };

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

        fn inv(self) -> Result<Self, RingDivisionError> {
            Ok(Self {
                p: self.p,
                rat: Rational::inv(&self.rat)?,
            })
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
            /*
            Let x be the root approximated by a: |x - a| <= v
            Let b be the rational value
            Then x+b is approxmated by a+b:
                |(x+b) - (a+b)| = |x-a| <= v
            */
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
            /*
            Let x be the root approximated by a: |x - a| <= v
            Let y be the root approximated by b: |y - b| <= w
            Then x+y is approximated by a+b:
                |(x + y) - (a + b)|
                = |(x - a) + (y - b)|
                <= max(|x - a|, |y - b|)
                <= max(v, w)
            */
            let p = a.p.clone();
            debug_assert_eq!(p, b.p);
            let rsp = root_sum_poly(&a.poly, &b.poly);
            let rsppsqfp = rsp.primitive_squarefree_part();
            let mut candidates = rsppsqfp.all_padic_roots(&p);
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

        fn mul_rat(&self, rat: &PAdicRational) -> PAdicAlgebraic {
            /*
            Let x be the root approximated by a: |x - a| <= v
            Let b be the rational value
            Then xb is approxmated by ab:
                |(xb) - (ab)|
                <= |x-a|.|b|
                = v.|b|
            */
            let p = self.p.clone();
            debug_assert_eq!(p, rat.p);
            if rat.rat == Rational::ZERO {
                PAdicAlgebraic::Rational(PAdicRational {
                    p,
                    rat: Rational::ZERO,
                })
            } else {
                let poly = root_rat_mul_poly(self.poly.clone(), &rat.rat);
                let approx = self.approx.clone().mul_rat(&p, &rat.rat);
                PAdicAlgebraic::Algebraic(Self { p, poly, approx })
            }
        }

        fn mul_mut(a: &mut Self, b: &mut Self) -> PAdicAlgebraic {
            /*
            Let x be the root approximated by a: |x - a| <= v
            Let y be the root approximated by b: |y - b| <= w
            Then xy is approximated by ab:
                |xy - ab|
                = |xy - xb + xb - ab|
                = |x(y - b) + b(x - a)|
                <= max(|x(y - b)|, |b(x - a)|)
                = max(|x|.|y-b|, |b|.|x - a|)
                <= max(|x|.w, |b|.v)
                <= |b|.v
                By a symmetric argument it is also <= |a|.w, so
                |xy - ab| <= max(|b|.v, |a|.w)

            */
            let p = a.p.clone();
            debug_assert_eq!(p, b.p);
            let mut candidates = root_product_poly(&a.poly, &b.poly)
                .primitive_squarefree_part()
                .all_padic_roots(&p);
            let mut k = Integer::ZERO;
            while candidates.len() > 1 {
                a.refine(&k);
                b.refine(&k);
                let aball = a.isolating_ball();
                let bball = b.isolating_ball();
                let v = std::cmp::min(
                    padic_rat_valuation(&p, bball.c.clone()) * aball.v.clone(),
                    padic_rat_valuation(&p, aball.c.clone()) * bball.v.clone(),
                );
                let sball = IsolatingBall {
                    p: p.clone(),
                    c: aball.c * bball.c,
                    v,
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

        fn inv_mut(&mut self) -> Result<PAdicAlgebraic, RingDivisionError> {
            /*
            Let x be the root approximated by a: |x-a| <= v
            Then 1/x is approximated by 1/a:
                |1/x - 1/a|
                = |(a - x) / xa|
                = |a-x| / |x|.|a|
                <= v / |x - a + a|.|a|
                <= v / max(|x - a|, |a|).|a|
                <= v / max(v, |a|).|a|
            */
            let inv_poly = Polynomial::from_coeffs(
                self.poly.clone().into_coeffs().into_iter().rev().collect(),
            )
            .fav_assoc();
            debug_assert!(inv_poly.is_irreducible());
            let p = self.p.clone();
            let mut candidates = inv_poly.all_padic_roots(&p);
            let mut k = Integer::ZERO;
            while candidates.len() > 1 {
                self.refine(&k);
                let sball = self.isolating_ball();
                if !sball.c.is_zero() {
                    let vc = padic_rat_valuation(&p, sball.c.clone());
                    #[cfg(debug_assertions)]
                    vc.clone().unwrap_int();
                    let iball = IsolatingBall {
                        p: p.clone(),
                        c: sball.c.inv().unwrap(),
                        v: &sball.v - std::cmp::min(&sball.v, &vc) - vc,
                    };
                    candidates = candidates
                        .into_iter()
                        .filter_map(|mut root| {
                            root.refine(&k);
                            let rball = root.isolating_ball();
                            match IsolatingBall::overlap(&rball, &iball) {
                                true => Some(root),
                                false => None,
                            }
                        })
                        .collect();
                }
                k += Integer::ONE;
            }
            debug_assert_eq!(candidates.len(), 1);
            Ok(candidates.into_iter().next().unwrap())
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

    impl Signature for PAdicAlgebraicStructure {}

    impl SetSignature for PAdicAlgebraicStructure {
        type Set = PAdicAlgebraic;

        fn is_element(&self, x: &Self::Set) -> bool {
            &self.p == x.p()
        }
    }

    impl PAdicAlgebraicStructure {
        pub fn new(p: Natural) -> Self {
            if !is_prime(&p) {
                panic!("{} is not prime", p)
            }
            Self { p }
        }
    }

    impl EqSignature for PAdicAlgebraicStructure {
        fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
            debug_assert!(self.is_element(a));
            debug_assert!(self.is_element(b));
            match (a, b) {
                (PAdicAlgebraic::Rational(a), PAdicAlgebraic::Rational(b)) => {
                    PAdicRational::equal(a, b)
                }
                (PAdicAlgebraic::Rational(_), PAdicAlgebraic::Algebraic(_))
                | (PAdicAlgebraic::Algebraic(_), PAdicAlgebraic::Rational(_))
                => false,
                (PAdicAlgebraic::Algebraic(a), PAdicAlgebraic::Algebraic(b)) => {
                    PAdicAlgebraicRoot::equal_mut(&mut a.clone(), &mut b.clone())
                }
            }
        }
    }

    impl SemiRingSignature for PAdicAlgebraicStructure {
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
            debug_assert!(self.is_element(a));
            debug_assert!(self.is_element(b));
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
            debug_assert!(self.is_element(a));
            debug_assert!(self.is_element(b));
            match (a, b) {
                (PAdicAlgebraic::Rational(a), PAdicAlgebraic::Rational(b)) => {
                    PAdicAlgebraic::Rational(PAdicRational {
                        p: self.p.clone(),
                        rat: &a.rat * &b.rat,
                    })
                }
                (PAdicAlgebraic::Rational(a), PAdicAlgebraic::Algebraic(b)) => b.mul_rat(a),
                (PAdicAlgebraic::Algebraic(a), PAdicAlgebraic::Rational(b)) => a.mul_rat(b),
                (PAdicAlgebraic::Algebraic(a), PAdicAlgebraic::Algebraic(b)) => {
                    PAdicAlgebraicRoot::mul_mut(&mut a.clone(), &mut b.clone())
                }
            }
        }
    }

    impl CharacteristicSignature for PAdicAlgebraicStructure {
        fn characteristic(&self) -> Natural {
            Natural::ZERO
        }
    }

    impl RingSignature for PAdicAlgebraicStructure {
        fn neg(&self, a: &Self::Set) -> Self::Set {
            debug_assert!(self.is_element(a));
            a.clone().neg()
        }
    }

    impl SemiRingUnitsSignature for PAdicAlgebraicStructure {
        fn inv(&self, a: &PAdicAlgebraic) -> Result<PAdicAlgebraic, RingDivisionError> {
            debug_assert!(self.is_element(a));
            match a {
                PAdicAlgebraic::Rational(a) => Ok(PAdicAlgebraic::Rational(a.clone().inv()?)),
                PAdicAlgebraic::Algebraic(a) => Ok(a.clone().inv_mut()?),
            }
        }
    }

    impl IntegralDomainSignature for PAdicAlgebraicStructure {
        fn div(&self, a: &Self::Set, b: &Self::Set) -> Result<Self::Set, RingDivisionError> {
            debug_assert!(self.is_element(a));
            debug_assert!(self.is_element(b));
            #[allow(clippy::single_match)]
            match (a, b) {
                (PAdicAlgebraic::Rational(a), PAdicAlgebraic::Rational(b)) => {
                    return Ok(PAdicAlgebraic::Rational(PAdicRational {
                        p: self.p.clone(),
                        rat: Rational::div(&a.rat, &b.rat)?,
                    }));
                }
                _ => {}
            }
            Ok(self.mul(a, &self.inv(b)?))
        }

        fn from_rat(&self, x: &Rational) -> Option<Self::Set> {
            Some(PAdicAlgebraic::Rational(PAdicRational {
                p: self.p.clone(),
                rat: x.clone(),
            }))
        }
    }

    impl FieldSignature for PAdicAlgebraicStructure {}

    #[cfg(test)]
    mod tests {
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

            let e = PAdicAlgebraic::from_rational(
                Natural::from(5u32),
                Rational::from_integers(Integer::from(-2), Integer::from(1)),
            );

            println!("a = {}", a);
            assert_eq!(
                a.clone().truncate(&6.into()).digits(),
                Some((
                    vec![2, 4, 2, 1, 1, 3]
                        .into_iter()
                        .map(|d| Natural::from(d as u8))
                        .collect(),
                    0.into()
                ))
            );
            println!("b = {}", b);
            assert_eq!(
                b.clone().truncate(&6.into()).digits(),
                Some((
                    vec![2, 3, 1, 2, 0, 1]
                        .into_iter()
                        .map(|d| Natural::from(d as u8))
                        .collect(),
                    0.into()
                ))
            );
            println!("c = {}", c);
            assert_eq!(
                c.clone().truncate(&6.into()).digits(),
                Some((
                    vec![1, 1, 3, 4, 1, 0]
                        .into_iter()
                        .map(|d| Natural::from(d as u8))
                        .collect(),
                    0.into()
                ))
            );
            println!("d = {}", d);
            assert_eq!(
                d.clone().truncate(&6.into()).digits(),
                Some((
                    vec![1, 2, 1, 4, 2, 3]
                        .into_iter()
                        .map(|d| Natural::from(d as u8))
                        .collect(),
                    0.into()
                ))
            );
            println!("e = {}", e);
            assert_eq!(
                e.clone().truncate(&6.into()).digits(),
                Some((
                    vec![3, 4, 4, 4, 4, 4]
                        .into_iter()
                        .map(|d| Natural::from(d as u8))
                        .collect(),
                    0.into()
                ))
            );

            println!("-a = {}", ring.neg(&a));
            assert_eq!(
                ring.neg(&a).truncate(&6.into()).digits(),
                Some((
                    vec![3, 0, 2, 3, 3, 1]
                        .into_iter()
                        .map(|d| Natural::from(d as u8))
                        .collect(),
                    0.into()
                ))
            );
            println!("-b = {}", ring.neg(&b));
            assert_eq!(
                ring.neg(&b).truncate(&6.into()).digits(),
                Some((
                    vec![3, 1, 3, 2, 4, 3]
                        .into_iter()
                        .map(|d| Natural::from(d as u8))
                        .collect(),
                    0.into()
                ))
            );
            println!("-c = {}", ring.neg(&c));
            assert_eq!(
                ring.neg(&c).truncate(&6.into()).digits(),
                Some((
                    vec![4, 3, 1, 0, 3, 4]
                        .into_iter()
                        .map(|d| Natural::from(d as u8))
                        .collect(),
                    0.into()
                ))
            );
            println!("-d = {}", ring.neg(&d));
            assert_eq!(
                ring.neg(&d).truncate(&6.into()).digits(),
                Some((
                    vec![4, 2, 3, 0, 2, 1]
                        .into_iter()
                        .map(|d| Natural::from(d as u8))
                        .collect(),
                    0.into()
                ))
            );
            println!("-e = {}", ring.neg(&e));
            assert_eq!(
                ring.neg(&e).truncate(&6.into()).digits(),
                Some((
                    vec![2, 0, 0, 0, 0, 0]
                        .into_iter()
                        .map(|d| Natural::from(d as u8))
                        .collect(),
                    0.into()
                ))
            );

            println!("a+b = {}", ring.add(&a, &b));
            assert_eq!(
                ring.add(&a, &b).truncate(&6.into()).digits(),
                Some((
                    vec![4, 2, 4, 3, 1, 4]
                        .into_iter()
                        .map(|d| Natural::from(d as u8))
                        .collect(),
                    0.into()
                ))
            );
            println!("a+c = {}", ring.add(&a, &c));
            assert_eq!(
                ring.add(&a, &c).truncate(&6.into()).digits(),
                Some((
                    vec![3, 0, 1, 1, 3, 3]
                        .into_iter()
                        .map(|d| Natural::from(d as u8))
                        .collect(),
                    0.into()
                ))
            );
            println!("d+b = {}", ring.add(&d, &b));
            assert_eq!(
                ring.add(&d, &b).truncate(&6.into()).digits(),
                Some((
                    vec![3, 0, 3, 1, 3, 4]
                        .into_iter()
                        .map(|d| Natural::from(d as u8))
                        .collect(),
                    0.into()
                ))
            );
            println!("d+c = {}", ring.add(&d, &c));
            assert_eq!(
                ring.add(&d, &c).truncate(&6.into()).digits(),
                Some((
                    vec![2, 3, 4, 3, 4, 3]
                        .into_iter()
                        .map(|d| Natural::from(d as u8))
                        .collect(),
                    0.into()
                ))
            );
            println!("c+c = {}", ring.add(&c, &c));
            assert_eq!(
                ring.add(&c, &c).truncate(&6.into()).digits(),
                Some((
                    vec![2, 2, 1, 4, 3, 0]
                        .into_iter()
                        .map(|d| Natural::from(d as u8))
                        .collect(),
                    0.into()
                ))
            );

            println!("a*b = {}", ring.mul(&a, &b));
            assert_eq!(
                ring.mul(&a, &b).truncate(&6.into()).digits(),
                Some((
                    vec![4, 4, 0, 0, 4, 4]
                        .into_iter()
                        .map(|d| Natural::from(d as u8))
                        .collect(),
                    0.into()
                ))
            );
            println!("a*c = {}", ring.mul(&a, &c));
            assert_eq!(
                ring.mul(&a, &c).truncate(&6.into()).digits(),
                Some((
                    vec![2, 1, 3, 0, 1, 0]
                        .into_iter()
                        .map(|d| Natural::from(d as u8))
                        .collect(),
                    0.into()
                ))
            );
            println!("d*b = {}", ring.mul(&d, &b));
            assert_eq!(
                ring.mul(&d, &b).truncate(&6.into()).digits(),
                Some((
                    vec![2, 2, 0, 2, 4, 3]
                        .into_iter()
                        .map(|d| Natural::from(d as u8))
                        .collect(),
                    0.into()
                ))
            );
            println!("d*c = {}", ring.mul(&d, &c));
            assert_eq!(
                ring.mul(&d, &c).truncate(&6.into()).digits(),
                Some((
                    vec![1, 3, 1, 1, 1, 2]
                        .into_iter()
                        .map(|d| Natural::from(d as u8))
                        .collect(),
                    0.into()
                ))
            );
            println!("c*c = {}", ring.mul(&c, &c));
            assert_eq!(
                ring.mul(&c, &c).truncate(&6.into()).digits(),
                Some((
                    vec![1, 2, 2, 0, 2, 0]
                        .into_iter()
                        .map(|d| Natural::from(d as u8))
                        .collect(),
                    0.into()
                ))
            );
            println!("a*e = {}", ring.mul(&a, &e));
            assert_eq!(
                ring.mul(&a, &e).truncate(&6.into()).digits(),
                Some((
                    vec![1, 1, 4, 1, 2, 3]
                        .into_iter()
                        .map(|d| Natural::from(d as u8))
                        .collect(),
                    0.into()
                ))
            );

            println!("a^-1 = {}", ring.inv(&a).unwrap());
            assert_eq!(
                ring.inv(&a).unwrap().truncate(&6.into()).digits(),
                Some((
                    vec![3, 1, 1, 4, 2, 1]
                        .into_iter()
                        .map(|d| Natural::from(d as u8))
                        .collect(),
                    0.into()
                ))
            );
            println!("b^-1 = {}", ring.inv(&b).unwrap());
            assert_eq!(
                ring.inv(&b).unwrap().truncate(&6.into()).digits(),
                Some((
                    vec![3, 0, 0, 4, 0, 0]
                        .into_iter()
                        .map(|d| Natural::from(d as u8))
                        .collect(),
                    0.into()
                ))
            );
            println!("c^-1 = {}", ring.inv(&c).unwrap());
            assert_eq!(
                ring.inv(&c).unwrap().truncate(&6.into()).digits(),
                Some((
                    vec![1, 4, 2, 0, 3, 4]
                        .into_iter()
                        .map(|d| Natural::from(d as u8))
                        .collect(),
                    0.into()
                ))
            );
            println!("d^-1 = {}", ring.inv(&d).unwrap());
            assert_eq!(
                ring.inv(&d).unwrap().truncate(&6.into()).digits(),
                Some((
                    vec![1, 3, 2, 2, 2, 2]
                        .into_iter()
                        .map(|d| Natural::from(d as u8))
                        .collect(),
                    0.into()
                ))
            );
            println!("e^-1 = {}", ring.inv(&e).unwrap());
            assert_eq!(
                ring.inv(&e).unwrap().truncate(&6.into()).digits(),
                Some((
                    vec![2, 2, 2, 2, 2, 2]
                        .into_iter()
                        .map(|d| Natural::from(d as u8))
                        .collect(),
                    0.into()
                ))
            );
        }
    }
}

#[cfg(test)]
mod tests {
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
