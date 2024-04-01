use std::collections::HashSet;
use std::f32::RADIX;
use std::fmt::Display;
use std::ops::Mul;
use std::rc::Rc;
use std::str::FromStr;

use itertools::Itertools;
use malachite_base::num::arithmetic::traits::NegAssign;
use malachite_base::num::basic::traits::One;
use malachite_base::num::basic::traits::OneHalf;
use malachite_base::num::basic::traits::Two;
use malachite_base::num::basic::traits::Zero;
use malachite_nz::integer::Integer;
use malachite_nz::natural::Natural;
use malachite_q::arithmetic::traits::SimplestRationalInInterval;
use malachite_q::Rational;
use rayon::collections::binary_heap::Iter;

use crate::rings::number::algebraic::bisection_gen::RationalSimpleBetweenGenerator;
use crate::rings::number::natural::*;
use crate::rings::polynomial::multipoly::*;
use crate::rings::polynomial::polynomial::*;
use crate::rings::ring_structure::cannonical::*;
use crate::rings::ring_structure::structure::*;
use crate::rings::structure::*;

fn rat_to_string(a: Rational) -> String {
    if a == 0 {
        return "0".into();
    }
    let neg = a < Rational::from(0);
    let (mant, exp): (f64, _) = a
        .sci_mantissa_and_exponent_with_rounding(
            malachite_base::rounding_modes::RoundingMode::Nearest,
        )
        .unwrap();
    let mut b = (2.0 as f64).powf(exp as f64) * mant;
    if neg {
        b = -b;
    }
    b = (1000.0 * b).round() / 1000.0;
    b.to_string()
}

fn root_sum_poly(p: &Polynomial<Integer>, q: &Polynomial<Integer>) -> Polynomial<Integer> {
    let x = Variable::new(String::from("x"));
    let z = Variable::new(String::from("z"));

    let p = p.apply_map(|c| MultiPolynomial::constant(c.clone()));
    let q = q.apply_map(|c| MultiPolynomial::constant(c.clone()));
    let r = q
        .evaluate(&MultiPolynomial::add(
            &MultiPolynomial::var(z.clone()),
            &MultiPolynomial::neg(&MultiPolynomial::var(x.clone())),
        ))
        .expand(&x);

    let root_sum_poly = Polynomial::resultant(&p, &r)
        .expand(&z)
        .apply_map(|c| MultiPolynomial::as_constant(c).unwrap());
    root_sum_poly.primitive_squarefree_part()
}

fn root_prod_poly(p: &Polynomial<Integer>, q: &Polynomial<Integer>) -> Polynomial<Integer> {
    let x = Variable::new(String::from("x"));
    let t = Variable::new(String::from("t"));

    let p = p.apply_map(|c| MultiPolynomial::constant(c.clone()));
    let q = q
        .apply_map(|c| MultiPolynomial::constant(c.clone()))
        .evaluate(&MultiPolynomial::var(x.clone()));
    let r = q.homogenize(&t).expand(&t);
    //x ** q.degree() * q(t * x ** -1)

    let root_prod_poly = Polynomial::resultant(&p, &r)
        .expand(&x)
        .apply_map(|c| MultiPolynomial::as_constant(c).unwrap());
    root_prod_poly.primitive_squarefree_part()
}

fn root_pos_rat_mul_poly(poly: Polynomial<Integer>, rat: &Rational) -> Polynomial<Integer> {
    debug_assert!(rat > &Rational::ZERO);
    debug_assert!(poly.is_irreducible());
    //we are multiplying by a so need to replace f(x) with f(x/a)
    //e.g. f(x) = x-1 and multiply root by 3 then replace f(x) with
    //f(x/3) = 3/x-1 = x-3
    //e.g. f(x) = 1 + x + x^2 replace it with f(d/n * x) = 1 + d/n x + d^2/n^2 x^2 = n^2 + ndx + d^2 x

    // println!("poly = {}", poly);
    // println!("rat = {}", rat);

    let rat_mul_poly = Polynomial::from_coeffs({
        let degree = poly.degree().unwrap();
        let (n, d) = (Rational::numerator(rat), Rational::denominator(rat));
        let mut n_pows = vec![Integer::from(1)];
        let mut d_pows = vec![Integer::from(1)];

        {
            let mut n_pow = n.clone();
            let mut d_pow = d.clone();
            for _i in 0..degree {
                n_pows.push(n_pow.clone());
                d_pows.push(d_pow.clone());
                n_pow *= &n;
                d_pow *= &d;
            }
        }

        // println!("n_pows = {:?}", n_pows);
        // println!("d_pows = {:?}", d_pows);

        debug_assert_eq!(n_pows.len(), degree + 1);
        debug_assert_eq!(d_pows.len(), degree + 1);

        let coeffs = poly
            .into_coeffs()
            .iter()
            .enumerate()
            .map(|(i, c)| &d_pows[i] * &n_pows[degree - i] * c)
            .collect();
        coeffs
    })
    .primitive_part()
    .unwrap();

    // println!("rat_mul_poly = {}", rat_mul_poly);
    // println!("rat_mul_poly = {}", rat_mul_poly.factor().unwrap());

    debug_assert!(rat_mul_poly.is_irreducible());

    rat_mul_poly
}

fn unique_linear_root(poly: &Polynomial<Integer>) -> Rational {
    debug_assert_eq!(poly.degree().unwrap(), 1);
    -Rational::from_integers(poly.coeff(0), poly.coeff(1))
}

fn evaluate_at_rational(poly: &Polynomial<Integer>, val: &Rational) -> Rational {
    poly.apply_map(|x| Rational::from(x)).evaluate(&val)
}

fn bisect_box(
    poly: &Polynomial<Integer>,
    n: usize,
    a: &Rational,
    b: &Rational,
    c: &Rational,
    d: &Rational,
) -> (
    (usize, Rational, Rational, Rational, Rational),
    (usize, Rational, Rational, Rational, Rational),
) {
    let ba = b - a;
    let dc = d - c;
    if ba >= dc {
        //bisect (a, b)
        for m in RationalSimpleBetweenGenerator::new_between_middle(
            a,
            b,
            &Rational::from_integers(Integer::from(1), Integer::from(3)),
        ) {
            match (
                poly.count_complex_roots(&a, &m, &c, &d),
                poly.count_complex_roots(&m, &b, &c, &d),
            ) {
                (Some(n1), Some(n2)) => {
                    debug_assert_eq!(n1 + n2, n);
                    return (
                        (n1, a.clone(), m.clone(), c.clone(), d.clone()),
                        (n2, m, b.clone(), c.clone(), d.clone()),
                    );
                }
                _ => {}
            }
        }
    } else {
        //bisect (c, d)
        for m in RationalSimpleBetweenGenerator::new_between_middle(
            c,
            d,
            &Rational::from_integers(Integer::from(1), Integer::from(3)),
        ) {
            match (
                poly.count_complex_roots(&a, &b, &c, &m),
                poly.count_complex_roots(&a, &b, &m, &d),
            ) {
                (Some(n1), Some(n2)) => {
                    debug_assert_eq!(n1 + n2, n);
                    return (
                        (n1, a.clone(), b.clone(), c.clone(), m.clone()),
                        (n2, a.clone(), b.clone(), m, d.clone()),
                    );
                }
                _ => {}
            }
        }
    }
    unreachable!()

    // let pgen = NaturalPrimeGenerator::new();
    // for p in pgen {
    //     let mut x = Natural::from(1u8);
    //     while x < p {
    //         {
    //             let f = Rational::from_naturals_ref(&x, &p);

    //             let ba = b - a;
    //             let dc = d - c;

    //             let ((a1, b1, c1, d1), (a2, b2, c2, d2)) = {
    //                 if ba >= dc {
    //                     let m = a + f * ba;
    //                     (
    //                         (a.clone(), m.clone(), c.clone(), d.clone()),
    //                         (m, b.clone(), c.clone(), d.clone()),
    //                     )
    //                 } else {
    //                     let m = c + f * dc;
    //                     (
    //                         (a.clone(), b.clone(), c.clone(), m.clone()),
    //                         (a.clone(), b.clone(), m, d.clone()),
    //                     )
    //                 }
    //             };

    //             match (
    //                 poly.count_complex_roots(&a1, &b1, &c1, &d1),
    //                 poly.count_complex_roots(&a2, &b2, &c2, &d2),
    //             ) {
    //                 (Some(n1), Some(n2)) => {
    //                     debug_assert_eq!(n1 + n2, n);
    //                     return ((n1, a1, b1, c1, d1), (n2, a2, b2, c2, d2));
    //                 }
    //                 _ => {}
    //             }
    //             x += Natural::from(1u8);
    //         }
    //     }
    // }
    // unreachable!();
}

#[derive(Debug, Clone)]
enum LowerBound {
    Inf,
    Finite(Rational),
}

#[derive(Debug, Clone)]
enum UpperBound {
    Inf,
    Finite(Rational),
}

impl LowerBound {
    pub fn neg(self) -> UpperBound {
        match self {
            LowerBound::Inf => UpperBound::Inf,
            LowerBound::Finite(a) => UpperBound::Finite(-a),
        }
    }
}

impl UpperBound {
    pub fn neg(self) -> LowerBound {
        match self {
            UpperBound::Inf => LowerBound::Inf,
            UpperBound::Finite(a) => LowerBound::Finite(-a),
        }
    }
}

impl std::hash::Hash for LowerBound {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        match self {
            LowerBound::Inf => {}
            LowerBound::Finite(x) => {
                x.hash(state);
            }
        }
    }
}

impl std::hash::Hash for UpperBound {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        match self {
            UpperBound::Inf => {}
            UpperBound::Finite(x) => {
                x.hash(state);
            }
        }
    }
}

impl PartialEq<UpperBound> for LowerBound {
    fn eq(&self, other: &UpperBound) -> bool {
        match (self, other) {
            (LowerBound::Finite(a), UpperBound::Finite(b)) => a == b,
            _ => false,
        }
    }
}

impl PartialOrd<UpperBound> for LowerBound {
    fn partial_cmp(&self, other: &UpperBound) -> Option<std::cmp::Ordering> {
        match (self, other) {
            (LowerBound::Finite(a), UpperBound::Finite(b)) => a.partial_cmp(b),
            _ => Some(std::cmp::Ordering::Less),
        }
    }
}

impl PartialEq<Rational> for LowerBound {
    fn eq(&self, b: &Rational) -> bool {
        match self {
            LowerBound::Finite(a) => a == b,
            _ => false,
        }
    }
}

impl PartialOrd<Rational> for LowerBound {
    fn partial_cmp(&self, b: &Rational) -> Option<std::cmp::Ordering> {
        match self {
            LowerBound::Finite(a) => a.partial_cmp(b),
            _ => Some(std::cmp::Ordering::Less),
        }
    }
}

impl PartialEq<UpperBound> for Rational {
    fn eq(&self, other: &UpperBound) -> bool {
        match other {
            UpperBound::Finite(b) => self == b,
            _ => false,
        }
    }
}

impl PartialOrd<UpperBound> for Rational {
    fn partial_cmp(&self, other: &UpperBound) -> Option<std::cmp::Ordering> {
        match other {
            UpperBound::Finite(b) => self.partial_cmp(b),
            _ => Some(std::cmp::Ordering::Less),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum Interleave {
    First,
    Second,
}

#[derive(Debug, Clone)]
enum SquarefreePolyRealRootInterval {
    Rational(Rational),
    //lower bound, upper bound, increasing
    //increasing = false : decreasing i.e. poly(a) > poly(b), true : increasing i.e. poly(a) < poly(b)
    Real(Rational, Rational, bool),
}

#[derive(Debug, Clone)]
struct SquarefreePolyRealRoots {
    poly_sqfr: Polynomial<Integer>,
    //an ordered list of isolating intervals for the squarefree polynomial
    //e.g. if r represents a real algebraic number and | represents a rational root
    //        (      r    )      |  ( r     )   |   |   (        r   )
    //note: it is allowed that some r might actually be rational but not known to be
    intervals: Vec<SquarefreePolyRealRootInterval>,
}

impl SquarefreePolyRealRoots {
    pub fn check_invariants(&self) -> Result<(), &'static str> {
        //poly should be squarefree
        if self
            .poly_sqfr
            .clone()
            .primitive_squarefree_part()
            .degree()
            .unwrap()
            != self.poly_sqfr.degree().unwrap()
        {
            return Err("poly should be squarefree");
        }

        //check the isolating intervals
        if self.intervals.len() != 0 {
            for i in 0..self.intervals.len() - 1 {
                let int1 = &self.intervals[i];
                let int2 = &self.intervals[i + 1];
                match (int1, int2) {
                    (
                        SquarefreePolyRealRootInterval::Rational(a),
                        SquarefreePolyRealRootInterval::Rational(x),
                    ) => {
                        if !(a < x) {
                            return Err("interval values should be strictly increasing");
                        }
                    }
                    (
                        SquarefreePolyRealRootInterval::Rational(a),
                        SquarefreePolyRealRootInterval::Real(x, y, _),
                    ) => {
                        if !(a < x) {
                            return Err("interval values should be strictly increasing");
                        }
                        if !(x < y) {
                            return Err("interval values should be strictly increasing");
                        }
                    }
                    (
                        SquarefreePolyRealRootInterval::Real(a, b, _),
                        SquarefreePolyRealRootInterval::Rational(x),
                    ) => {
                        if !(a < b) {
                            return Err("interval values should be strictly increasing");
                        }
                        if !(b < x) {
                            return Err("interval values should be strictly increasing");
                        }
                    }
                    (
                        SquarefreePolyRealRootInterval::Real(a, b, _),
                        SquarefreePolyRealRootInterval::Real(x, y, _),
                    ) => {
                        if !(a < b) {
                            return Err("interval values should be strictly increasing");
                        }
                        if !(b <= x) {
                            return Err("interval values should be increasing");
                        }
                        if !(x < y) {
                            return Err("interval values should be strictly increasing");
                        }
                    }
                }
            }
        }

        for interval in self.intervals.iter() {
            match interval {
                SquarefreePolyRealRootInterval::Rational(a) => {
                    if evaluate_at_rational(&self.poly_sqfr, a) != Rational::from(0) {
                        return Err("poly should be zero at a rational root");
                    }
                }
                SquarefreePolyRealRootInterval::Real(a, b, incr) => {
                    let at_a = evaluate_at_rational(&self.poly_sqfr, a);
                    let at_b = evaluate_at_rational(&self.poly_sqfr, b);

                    if at_a == Rational::from(0) || at_b == Rational::from(0) {
                        return Err("poly should not be zero at boundary of isolating interval");
                    }

                    if (at_a > Rational::from(0)) == (at_b > Rational::from(0)) {
                        return Err("sign of poly should be different at a and at b");
                    }

                    if *incr {
                        if !((at_a < Rational::from(0)) && (at_b > Rational::from(0))) {
                            return Err("sign of poly should go from neg to pos here");
                        }
                    } else {
                        if !((at_a > Rational::from(0)) && (at_b < Rational::from(0))) {
                            return Err("sign of poly should go from pos to neg here");
                        }
                    }
                }
            }
        }

        Ok(())
    }

    fn get_wide_interval(&self, idx: usize) -> (LowerBound, UpperBound) {
        assert!(idx < self.intervals.len());

        let wide_a = {
            if idx == 0 {
                LowerBound::Inf
            } else {
                LowerBound::Finite({
                    match &self.intervals[idx - 1] {
                        SquarefreePolyRealRootInterval::Rational(a) => a.clone(),
                        SquarefreePolyRealRootInterval::Real(_, prev_b, _) => prev_b.clone(),
                    }
                })
            }
        };
        let wide_b = {
            if idx == self.intervals.len() - 1 {
                UpperBound::Inf
            } else {
                UpperBound::Finite({
                    match &self.intervals[idx + 1] {
                        SquarefreePolyRealRootInterval::Rational(a) => a.clone(),
                        SquarefreePolyRealRootInterval::Real(prev_a, _, _) => prev_a.clone(),
                    }
                })
            }
        };
        debug_assert!(wide_a.clone() < wide_b.clone());
        (wide_a, wide_b)
    }

    fn refine(&mut self, idx: usize) {
        assert!(idx < self.intervals.len());

        match &mut self.intervals[idx] {
            SquarefreePolyRealRootInterval::Rational(_a) => {}
            SquarefreePolyRealRootInterval::Real(a, b, dir) => {
                let m = RationalSimpleBetweenGenerator::new_between_middle(
                    a,
                    b,
                    &Rational::from_integers(Integer::from(1), Integer::from(3)),
                )
                .next()
                .unwrap();

                match evaluate_at_rational(&self.poly_sqfr, &m).cmp(&Rational::from(0)) {
                    std::cmp::Ordering::Less => match dir {
                        true => {
                            *a = m;
                        }
                        false => {
                            *b = m;
                        }
                    },
                    std::cmp::Ordering::Equal => {
                        self.intervals[idx] = SquarefreePolyRealRootInterval::Rational(m);
                    }
                    std::cmp::Ordering::Greater => match dir {
                        true => {
                            *b = m;
                        }
                        false => {
                            *a = m;
                        }
                    },
                }
            }
        }
    }

    fn refine_all(&mut self) {
        for idx in 0..self.intervals.len() {
            self.refine(idx);
        }
    }

    fn get_real_root(&self, idx: usize) -> RealAlgebraic {
        // println!("get_real_root of {:?} at {}", self, idx);
        // println!("{}", self.poly_sqfr);
        assert!(idx < self.intervals.len());
        match &self.intervals[idx] {
            SquarefreePolyRealRootInterval::Rational(rat) => RealAlgebraic::Rational(rat.clone()),
            SquarefreePolyRealRootInterval::Real(a, b, dir) => {
                let (unit, factors) = self.poly_sqfr.factor().unwrap().unit_and_factors();
                for (factor, k) in factors.into_iter() {
                    // println!("factor = {}", factor);
                    debug_assert_eq!(k, Natural::ONE); //square free
                    let deg = factor.degree().unwrap();
                    debug_assert_ne!(deg, 0);
                    let at_a = evaluate_at_rational(&factor, a);
                    let at_b = evaluate_at_rational(&factor, b);
                    debug_assert_ne!(at_a, Rational::ZERO);
                    debug_assert_ne!(at_b, Rational::ZERO);
                    let sign_a = at_a >= Rational::ZERO;
                    let sign_b = at_b >= Rational::ZERO;
                    // println!("at_a = {}", at_a);
                    // println!("at_b = {}", at_b);
                    if deg == 1 {
                        if sign_a != sign_b {
                            return RealAlgebraic::Rational(unique_linear_root(&factor));
                        }
                    } else {
                        if sign_a && !sign_b {
                            let (wide_a, wide_b) = self.get_wide_interval(idx);
                            return RealAlgebraic::Real(RealAlgebraicRoot {
                                poly: factor,
                                tight_a: a.clone(),
                                tight_b: b.clone(),
                                wide_a: wide_a,
                                wide_b: wide_b,
                                dir: false,
                            });
                        } else if !sign_a && sign_b {
                            let (wide_a, wide_b) = self.get_wide_interval(idx);
                            return RealAlgebraic::Real(RealAlgebraicRoot {
                                poly: factor,
                                tight_a: a.clone(),
                                tight_b: b.clone(),
                                wide_a: wide_a,
                                wide_b: wide_b,
                                dir: true,
                            });
                        }
                    }
                    debug_assert_eq!(sign_a, sign_b);
                }
                unreachable!()
            }
        }
    }

    fn to_real_roots(self) -> Vec<RealAlgebraic> {
        debug_assert!(self.poly_sqfr.is_irreducible());
        let deg = self.poly_sqfr.degree().unwrap();
        if deg == 0 {
            vec![]
        } else if deg == 1 {
            if self.intervals.len() == 0 {
                vec![]
            } else if self.intervals.len() == 1 {
                match self.intervals.into_iter().next().unwrap() {
                    SquarefreePolyRealRootInterval::Rational(a) => {
                        vec![RealAlgebraic::Rational(a)]
                    }
                    SquarefreePolyRealRootInterval::Real(_, _, _) => {
                        // panic!("degree 1 polynomial should have had rational root found exactly");
                        vec![RealAlgebraic::Rational(unique_linear_root(&self.poly_sqfr))]
                    }
                }
            } else {
                panic!();
            }
        } else {
            let mut roots = vec![];
            for (idx, interval) in self.intervals.iter().enumerate() {
                roots.push({
                    match interval {
                        SquarefreePolyRealRootInterval::Rational(a) => {
                            RealAlgebraic::Rational(a.clone())
                        }
                        SquarefreePolyRealRootInterval::Real(tight_a, tight_b, dir) => {
                            let (wide_a, wide_b) = self.get_wide_interval(idx);
                            RealAlgebraic::Real(RealAlgebraicRoot {
                                poly: self.poly_sqfr.clone(),
                                tight_a: tight_a.clone(),
                                tight_b: tight_b.clone(),
                                wide_a,
                                wide_b,
                                dir: *dir,
                            })
                        }
                    }
                });
            }
            roots
        }
    }

    //separate the isolating intervals of the roots in roots1 and roots2
    //return Err if a root in roots1 and a root in roots2 are equal and thus cant be separated
    //ends of real root intervals should not equal rational roots in the other one
    fn separate(roots1: &mut Self, roots2: &mut Self) -> Result<Vec<(Interleave, usize)>, ()> {
        let poly_gcd_sqfr = Polynomial::subresultant_gcd(&roots1.poly_sqfr, &roots2.poly_sqfr);
        let (_, poly_gcd_sqfr) = poly_gcd_sqfr.factor_primitive().unwrap();

        //compute which of the roots are equals to some root from the other one
        let is_gcdroot_1: Vec<_> = roots1
            .intervals
            .iter()
            .map(|root| match root {
                SquarefreePolyRealRootInterval::Rational(x) => {
                    evaluate_at_rational(&poly_gcd_sqfr, x) == Rational::from(0)
                }
                SquarefreePolyRealRootInterval::Real(a, b, _dir) => {
                    debug_assert_ne!(evaluate_at_rational(&poly_gcd_sqfr, a), Rational::from(0));
                    debug_assert_ne!(evaluate_at_rational(&poly_gcd_sqfr, b), Rational::from(0));
                    (evaluate_at_rational(&poly_gcd_sqfr, a) > Rational::from(0))
                        != (evaluate_at_rational(&poly_gcd_sqfr, b) > Rational::from(0))
                }
            })
            .collect();
        let is_gcdroot_2: Vec<_> = roots2
            .intervals
            .iter()
            .map(|root| match root {
                SquarefreePolyRealRootInterval::Rational(x) => {
                    evaluate_at_rational(&poly_gcd_sqfr, x) == Rational::from(0)
                }
                SquarefreePolyRealRootInterval::Real(a, b, _dir) => {
                    debug_assert_ne!(evaluate_at_rational(&poly_gcd_sqfr, a), Rational::from(0));
                    debug_assert_ne!(evaluate_at_rational(&poly_gcd_sqfr, b), Rational::from(0));
                    (evaluate_at_rational(&poly_gcd_sqfr, a) > Rational::from(0))
                        != (evaluate_at_rational(&poly_gcd_sqfr, b) > Rational::from(0))
                }
            })
            .collect();

        //do the separation
        let mut all_roots = vec![];

        let mut idx1 = 0;
        let mut idx2 = 0;
        while idx1 < roots1.intervals.len() && idx2 < roots2.intervals.len() {
            let (wide1_a, wide1_b) = roots1.get_wide_interval(idx1);
            let (wide2_a, wide2_b) = roots2.get_wide_interval(idx2);

            let root1 = &roots1.intervals[idx1];
            let root2 = &roots2.intervals[idx2];

            //check if the roots are equal
            if is_gcdroot_1[idx1] && is_gcdroot_2[idx2] {
                match root1 {
                    SquarefreePolyRealRootInterval::Rational(x1) => {
                        if &wide2_a < x1 && x1 < &wide2_b {
                            return Err(());
                        }
                    }
                    SquarefreePolyRealRootInterval::Real(a1, b1, _dir1) => {
                        if &wide2_a < a1 && b1 < &wide2_b {
                            return Err(());
                        }
                    }
                }
                match root2 {
                    SquarefreePolyRealRootInterval::Rational(x2) => {
                        if &wide1_a < x2 && x2 < &wide1_b {
                            return Err(());
                        }
                    }
                    SquarefreePolyRealRootInterval::Real(a2, b2, _dir2) => {
                        if &wide1_a < a2 && b2 < &wide1_b {
                            return Err(());
                        }
                    }
                }
            }

            //check if one is bigger than the other
            match (root1, root2) {
                (
                    SquarefreePolyRealRootInterval::Rational(x1),
                    SquarefreePolyRealRootInterval::Rational(x2),
                ) => match x1.cmp(x2) {
                    std::cmp::Ordering::Less => {
                        all_roots.push((Interleave::First, idx1));
                        idx1 += 1;
                        continue;
                    }
                    std::cmp::Ordering::Equal => panic!(),
                    std::cmp::Ordering::Greater => {
                        all_roots.push((Interleave::Second, idx2));
                        idx2 += 1;
                        continue;
                    }
                },
                (
                    SquarefreePolyRealRootInterval::Rational(x1),
                    SquarefreePolyRealRootInterval::Real(a2, b2, _dir2),
                ) => {
                    if x1 < a2 {
                        all_roots.push((Interleave::First, idx1));
                        idx1 += 1;
                        continue;
                    }
                    if b2 < x1 {
                        all_roots.push((Interleave::Second, idx2));
                        idx2 += 1;
                        continue;
                    }
                }
                (
                    SquarefreePolyRealRootInterval::Real(a1, b1, _dir1),
                    SquarefreePolyRealRootInterval::Rational(x2),
                ) => {
                    if x2 < a1 {
                        all_roots.push((Interleave::Second, idx2));
                        idx2 += 1;
                        continue;
                    }
                    if b1 < x2 {
                        all_roots.push((Interleave::First, idx1));
                        idx1 += 1;
                        continue;
                    }
                }
                (
                    SquarefreePolyRealRootInterval::Real(a1, b1, _dir1),
                    SquarefreePolyRealRootInterval::Real(a2, b2, _dir2),
                ) => {
                    if b2 < a1 {
                        all_roots.push((Interleave::Second, idx2));
                        idx2 += 1;
                        continue;
                    }
                    if b1 < a2 {
                        all_roots.push((Interleave::First, idx1));
                        idx1 += 1;
                        continue;
                    }
                }
            }

            //refine and try again
            roots1.refine(idx1);
            roots2.refine(idx2);
        }

        debug_assert!(idx1 == roots1.intervals.len() || idx2 == roots2.intervals.len());

        while idx1 < roots1.intervals.len() {
            all_roots.push((Interleave::First, idx1));
            idx1 += 1;
        }

        while idx2 < roots2.intervals.len() {
            all_roots.push((Interleave::Second, idx2));
            idx2 += 1;
        }

        for r1 in &roots1.intervals {
            for r2 in &roots2.intervals {
                match (r1, r2) {
                    (
                        SquarefreePolyRealRootInterval::Rational(a),
                        SquarefreePolyRealRootInterval::Rational(x),
                    ) => {
                        debug_assert!(a != x);
                    }
                    (
                        SquarefreePolyRealRootInterval::Rational(a),
                        SquarefreePolyRealRootInterval::Real(x, y, _),
                    ) => {
                        debug_assert!(a < x || y < a);
                    }
                    (
                        SquarefreePolyRealRootInterval::Real(a, b, _),
                        SquarefreePolyRealRootInterval::Rational(x),
                    ) => {
                        debug_assert!(x < a || b < x);
                    }
                    (
                        SquarefreePolyRealRootInterval::Real(a, b, _),
                        SquarefreePolyRealRootInterval::Real(x, y, _),
                    ) => {
                        debug_assert!(b < x || y < a);
                    }
                }
            }
        }

        debug_assert_eq!(
            all_roots.len(),
            roots1.intervals.len() + roots2.intervals.len()
        );

        Ok(all_roots)
    }
}

impl Polynomial<Integer> {
    fn sign_variations(&self) -> usize {
        //https://en.wikipedia.org/wiki/Descartes'_rule_of_signs
        //equals the number of strictly positive real roots modulo 2
        //and number of positive real roots is less than this number
        let nonzero_coeffs = self
            .coeffs()
            .into_iter()
            .filter(|c| c != &&Integer::zero())
            .collect_vec();
        let mut v = 0;
        for i in 0..nonzero_coeffs.len() - 1 {
            if (nonzero_coeffs[i] < &0) != (nonzero_coeffs[i + 1] < &0) {
                v += 1;
            }
        }
        v
    }

    //Collins and Akritas algorithm %https://en.wikipedia.org/wiki/Real-root_isolation
    fn isolate_real_roots_by_collin_akritas(&self) -> Vec<(Natural, usize, bool)> {
        //input: p(x), a square-free polynomial, such that p(0) p(1) ≠ 0, for which the roots in the interval [0, 1] are searched
        //output: a list of triples (c, k, h) representing isolating intervals of the form [c/2^k, (c+h)/2^k]
        debug_assert_ne!(self.evaluate(&Integer::zero()), Integer::zero());
        debug_assert_ne!(self.evaluate(&Integer::one()), Integer::zero());
        debug_assert_eq!(
            self.clone().primitive_squarefree_part().degree().unwrap(),
            self.degree().unwrap()
        );

        let mut l = vec![(Natural::from(0u8), 0, self.clone())];
        let mut isol = vec![];
        while l.len() != 0 {
            let (c, k, mut q) = l.pop().unwrap();
            if q.evaluate(&Integer::from(0)) == Integer::from(0) {
                //q = q/x
                q = Self::div(&q, &Self::var()).unwrap();
                isol.push((c.clone(), k.clone(), false)); //rational root
            }
            let v = Self::compose(
                &q.reversed(),
                &Self::from_coeffs(vec![Integer::from(1), Integer::from(1)]),
            )
            .sign_variations();
            if v == 1 {
                isol.push((c, k, true)); //root
            } else if v > 1 {
                //bisect
                //q_small(x) = 2^n q(x/2)
                let q_small = q.apply_map_with_powers(|(i, coeff)| {
                    coeff * Integer::from(2) << (q.degree().unwrap() - i)
                });
                l.push((
                    (c.clone() << 1) + Natural::from(1u8),
                    k + 1,
                    Self::compose(
                        &q_small,
                        &Self::from_coeffs(vec![Integer::from(1), Integer::from(1)]),
                    ),
                ));
                l.push((c << 1, k + 1, q_small));
            }
        }
        isol
    }

    //isolate all real roots of a squarefree (no repeated roots) polynomial between a and b
    fn real_roots_squarefree(
        self,
        opt_a: Option<&Rational>,
        opt_b: Option<&Rational>,
        include_a: bool,
        include_b: bool,
    ) -> SquarefreePolyRealRoots {
        assert_ne!(self, Self::zero());
        //poly should be squarefree
        debug_assert_eq!(
            self.clone().primitive_squarefree_part().degree().unwrap(),
            self.degree().unwrap()
        );

        match (opt_a, opt_b) {
            (Some(a), Some(b)) => {
                assert!(a < b);
            }
            _ => {}
        }

        let d = self.degree().unwrap();
        if d == 0 {
            //constant polynomial has no roots
            SquarefreePolyRealRoots {
                poly_sqfr: self,
                intervals: vec![],
            }
        } else if d == 1 {
            //poly = a+bx
            //root = -a/b
            let root = unique_linear_root(&self);

            if {
                match opt_a {
                    Some(a) => match a.cmp(&root) {
                        std::cmp::Ordering::Less => true,
                        std::cmp::Ordering::Equal => include_a,
                        std::cmp::Ordering::Greater => false,
                    },
                    None => true,
                }
            } && {
                match opt_b {
                    Some(b) => match b.cmp(&root) {
                        std::cmp::Ordering::Greater => true,
                        std::cmp::Ordering::Equal => include_b,
                        std::cmp::Ordering::Less => false,
                    },
                    None => true,
                }
            } {
                SquarefreePolyRealRoots {
                    poly_sqfr: self,
                    intervals: vec![SquarefreePolyRealRootInterval::Rational(root)],
                }
            } else {
                SquarefreePolyRealRoots {
                    poly_sqfr: self,
                    intervals: vec![],
                }
            }
        } else {
            if opt_a.is_none() || opt_b.is_none() {
                //compute a bound M on the absolute value of any root
                //m = (Cauchy's bound + 1) https://captainblack.wordpress.com/2009/03/08/cauchys-upper-bound-for-the-roots-of-a-polynomial/
                let m = Rational::from(2)
                    + Rational::from_naturals(
                        itertools::max((0..d).map(|i| self.coeff(i).unsigned_abs_ref().clone()))
                            .unwrap(),
                        self.coeff(d).unsigned_abs_ref().clone(),
                    );

                debug_assert!(m > Rational::ZERO);

                return match opt_a {
                    Some(a_val) => match opt_b {
                        Some(_b_val) => panic!(),
                        None => {
                            self.real_roots_squarefree(Some(a_val), Some(&m), include_a, include_b)
                        }
                    },
                    None => match opt_b {
                        Some(b_val) => {
                            self.real_roots_squarefree(Some(&-m), Some(b_val), include_a, include_b)
                        }
                        None => {
                            let neg_m = -m.clone();
                            self.real_roots_squarefree(Some(&neg_m), Some(&m), include_a, include_b)
                        }
                    },
                };
            }
            let (a, b) = (opt_a.unwrap(), opt_b.unwrap());
            debug_assert!(a < b);

            //deal with end roots
            let mut poly_no_endroots = self.clone();
            let mut intervals = vec![];
            if evaluate_at_rational(&self, a) == Rational::from(0) {
                poly_no_endroots = Self::div(
                    &poly_no_endroots,
                    &Polynomial::from_coeffs(vec![
                        -Rational::numerator(a),
                        Rational::denominator(a),
                    ]),
                )
                .unwrap();
                if include_a {
                    intervals.push(SquarefreePolyRealRootInterval::Rational(a.clone()));
                }
            }
            let mut do_add_b = false;
            if evaluate_at_rational(&self, b) == Rational::from(0) {
                poly_no_endroots = Self::div(
                    &poly_no_endroots,
                    &Polynomial::from_coeffs(vec![
                        -Rational::numerator(b),
                        Rational::denominator(b),
                    ]),
                )
                .unwrap();
                if include_b {
                    do_add_b = true;
                }
            }

            debug_assert_ne!(
                evaluate_at_rational(&poly_no_endroots, a),
                Rational::from(0)
            );
            debug_assert_ne!(
                evaluate_at_rational(&poly_no_endroots, b),
                Rational::from(0)
            );

            //apply a transformation to p so that its roots in (a, b) are moved to roots in (0, 1)
            let (_, trans_poly) = Polynomial::compose(
                &poly_no_endroots.apply_map(|c| Rational::from(c)),
                &Polynomial::from_coeffs(vec![a.clone(), b.clone() - a.clone()]),
            )
            .factor_primitive_fof();

            'interval_loop: for (c, k, h) in trans_poly.isolate_real_roots_by_collin_akritas() {
                let d = Natural::from(1u8) << k;
                let mut interval_a = (b - a) * Rational::from_naturals(c.clone(), d.clone()) + a;
                if h {
                    let mut interval_b =
                        (b - a) * Rational::from_naturals(&c + Natural::from(1u8), d.clone()) + a;

                    //at the moment, interval_a and interval_b might be rational roots
                    //we need to strink them a little bit if so
                    if evaluate_at_rational(&self, &interval_a) == Rational::from(0)
                        || evaluate_at_rational(&self, &interval_b) == Rational::from(0)
                    {
                        let interval_m = (&interval_a + &interval_b) / Rational::from(2);
                        let mut shrunk_inerval_a = (&interval_a + &interval_m) / Rational::from(2);
                        let mut shrunk_inerval_b = (&interval_m + &interval_b) / Rational::from(2);
                        loop {
                            let shrunk_a_sign = evaluate_at_rational(&self, &shrunk_inerval_a)
                                .cmp(&Rational::from(0));
                            let shrunk_b_sign = evaluate_at_rational(&self, &shrunk_inerval_b)
                                .cmp(&Rational::from(0));
                            if shrunk_a_sign.is_eq() {
                                intervals.push(SquarefreePolyRealRootInterval::Rational(
                                    shrunk_inerval_a,
                                ));
                                continue 'interval_loop;
                            } else if shrunk_b_sign.is_eq() {
                                intervals.push(SquarefreePolyRealRootInterval::Rational(
                                    shrunk_inerval_b,
                                ));
                                continue 'interval_loop;
                            } else if shrunk_a_sign.is_ge() != shrunk_b_sign.is_ge() {
                                break;
                            }
                            shrunk_inerval_a =
                                (&interval_a + &shrunk_inerval_a) / Rational::from(2);
                            shrunk_inerval_b =
                                (&shrunk_inerval_b + &interval_b) / Rational::from(2);
                        }
                        interval_a = shrunk_inerval_a;
                        interval_b = shrunk_inerval_b;
                    }

                    let sign_a = evaluate_at_rational(&self, &interval_a) > Rational::from(0);
                    let sign_b = evaluate_at_rational(&self, &interval_b) > Rational::from(0);
                    debug_assert_ne!(evaluate_at_rational(&self, &interval_a), Rational::from(0));
                    debug_assert_ne!(evaluate_at_rational(&self, &interval_b), Rational::from(0));
                    debug_assert_ne!(sign_a, sign_b);
                    intervals.push(SquarefreePolyRealRootInterval::Real(
                        interval_a, interval_b, sign_b,
                    ));
                } else {
                    intervals.push(SquarefreePolyRealRootInterval::Rational(interval_a));
                }
            }

            if do_add_b {
                intervals.push(SquarefreePolyRealRootInterval::Rational(b.clone()));
            }

            let roots = SquarefreePolyRealRoots {
                poly_sqfr: self,
                intervals,
            };
            #[cfg(debug_assertions)]
            roots.check_invariants().unwrap();
            roots
        }
    }

    //isolate all real roots of the irreducible poly in the open interval (a, b)
    fn all_real_roots_squarefree(&self) -> SquarefreePolyRealRoots {
        self.clone().real_roots_squarefree(None, None, false, false)
    }

    //isolate all real roots of the irreducible poly in the open interval (a, b)
    fn real_roots_irreducible(
        &self,
        opt_a: Option<&Rational>,
        opt_b: Option<&Rational>,
        include_a: bool,
        include_b: bool,
    ) -> Vec<RealAlgebraic> {
        debug_assert!(self.is_irreducible());

        self.clone()
            .real_roots_squarefree(opt_a, opt_b, include_a, include_b)
            .to_real_roots()
    }

    //get the real roots with multiplicity of poly
    pub fn real_roots(
        &self,
        a: Option<&Rational>,
        b: Option<&Rational>,
        include_a: bool,
        include_b: bool,
    ) -> Vec<RealAlgebraic> {
        assert_ne!(self, &Self::zero());
        let factors = self.factor().unwrap();
        let mut roots = vec![];
        for (factor, k) in factors.factors() {
            for root in factor.real_roots_irreducible(a, b, include_a, include_b) {
                let mut i = Natural::from(0u8);
                while &i < k {
                    roots.push(root.clone());
                    i += Natural::from(1u8);
                }
            }
        }
        roots
    }

    pub fn all_real_roots(&self) -> Vec<RealAlgebraic> {
        self.real_roots(None, None, false, false)
    }

    fn at_fixed_re_or_im_impl<const RE_OR_IM: bool>(
        &self,
        a: &Rational,
    ) -> (Polynomial<Integer>, Polynomial<Integer>) {
        //find real and imag polys of
        //poly(a + xi) if RE_OR_IM = false
        //poly(x + ai) if RE_OR_IM = true
        //up to rational multiples (its the roots we care about)
        match self.degree() {
            Some(n) => {
                let (a_numer, a_denom) = (Rational::numerator(a), Rational::denominator(a));
                //multiply everything by a_d^n so that everything is integers

                //compute 1, a, a^2, a^3, ..., a^n (after multiplying everything by a_d)
                // a_d^n(a_n/a_d)^k = a_n^k a_d^{n-k}
                let mut a_numer_pow = vec![Integer::from(1)];
                let mut a_denom_pow = vec![Integer::from(1)];
                for k in 1..n + 1 {
                    a_numer_pow.push(&a_numer * &a_numer_pow[k - 1]);
                    a_denom_pow.push(&a_denom * &a_denom_pow[k - 1]);
                }
                let mut a_pow = vec![];
                for k in 0..n + 1 {
                    a_pow.push(&a_numer_pow[k] * &a_denom_pow[n - k]);
                }

                let mut re = Vec::with_capacity(n + 1);
                let mut im = Vec::with_capacity(n + 1);
                for _ in 0..n + 1 {
                    re.push(Integer::from(0));
                    im.push(Integer::from(0));
                }
                let mut n_choose = vec![Integer::from(1)];
                for n in 0..n + 1 {
                    if n == 0 {
                        debug_assert_eq!(n_choose, vec![Integer::from(1)]);
                    } else if n == 1 {
                        debug_assert_eq!(n_choose, vec![Integer::from(1), Integer::from(1)]);
                    } else if n == 2 {
                        debug_assert_eq!(
                            n_choose,
                            vec![Integer::from(1), Integer::from(2), Integer::from(1)]
                        );
                    } else if n == 3 {
                        debug_assert_eq!(
                            n_choose,
                            vec![
                                Integer::from(1),
                                Integer::from(3),
                                Integer::from(3),
                                Integer::from(1)
                            ]
                        );
                    }

                    //if fixed real add
                    //(a + xi)^n = \sum_{k=0,1,...,n} \binom{n}{k} a^{n-k} (xi)^k
                    //           = \sum_{k=0,1,...,n} \binom{n}{k} a^{n-k} x^k i^k
                    //           = \sum_{k=0,1,...,n} {
                    //               k = 0 mod 4        + \binom{n}{k} a^{n-k} x^k
                    //               k = 1 mod 4        + \binom{n}{k} a^{n-k} x^k i
                    //               k = 2 mod 4        - \binom{n}{k} a^{n-k} x^k
                    //               k = 3 mod 4        - \binom{n}{k} a^{n-k} x^k i
                    //                                }
                    //
                    //if fixed imag add
                    //(a + xi)^n = \sum_{k=0,1,...,n} \binom{n}{k} a^{n-k} (xi)^k
                    //           = \sum_{k=0,1,...,n} \binom{n}{k} a^{n-k} x^k i^k
                    //           = \sum_{k=0,1,...,n} {
                    //               k = 0 mod 4        + \binom{n}{k} a^{n-k} x^k
                    //               k = 1 mod 4        + \binom{n}{k} a^{n-k} x^k i
                    //               k = 2 mod 4        - \binom{n}{k} a^{n-k} x^k
                    //               k = 3 mod 4        - \binom{n}{k} a^{n-k} x^k i
                    //
                    if self.coeff(n) != Integer::from(0) {
                        let mut k = 0;
                        loop {
                            //k = 0 mod 4
                            re[{
                                match RE_OR_IM {
                                    false => k,
                                    true => n - k,
                                }
                            }] += self.coeff(n)
                                * &n_choose[k]
                                * &a_pow[{
                                    match RE_OR_IM {
                                        false => n - k,
                                        true => k,
                                    }
                                }];
                            if k == n {
                                break;
                            }
                            k += 1;
                            //k = 1 mod 4
                            im[{
                                match RE_OR_IM {
                                    false => k,
                                    true => n - k,
                                }
                            }] += self.coeff(n)
                                * &n_choose[k]
                                * &a_pow[{
                                    match RE_OR_IM {
                                        false => n - k,
                                        true => k,
                                    }
                                }];
                            if k == n {
                                break;
                            }
                            k += 1;
                            //k = 2 mod 4
                            re[{
                                match RE_OR_IM {
                                    false => k,
                                    true => n - k,
                                }
                            }] -= self.coeff(n)
                                * &n_choose[k]
                                * &a_pow[{
                                    match RE_OR_IM {
                                        false => n - k,
                                        true => k,
                                    }
                                }];
                            if k == n {
                                break;
                            }
                            k += 1;
                            //k = 3 mod 4
                            im[{
                                match RE_OR_IM {
                                    false => k,
                                    true => n - k,
                                }
                            }] -= self.coeff(n)
                                * &n_choose[k]
                                * &a_pow[{
                                    match RE_OR_IM {
                                        false => n - k,
                                        true => k,
                                    }
                                }];
                            if k == n {
                                break;
                            }
                            k += 1;
                        }
                    }
                    //update n choose k
                    //e.g. for n=3 do
                    //[1, 3, 3, 1]
                    //[1, 3, 3, 1, 1]
                    //[1, 3, 3, 4, 1]
                    //[1, 3, 6, 4, 1]
                    //[1, 4, 6, 4, 1]
                    n_choose.push(Integer::from(1));
                    for i in (1..n + 1).rev() {
                        n_choose[i] = &n_choose[i] + &n_choose[i - 1];
                    }
                }
                (Polynomial::from_coeffs(re), Polynomial::from_coeffs(im))
            }
            None => (Self::zero(), Self::zero()),
        }
    }

    fn at_fixed_re(&self, a: &Rational) -> (Polynomial<Integer>, Polynomial<Integer>) {
        self.at_fixed_re_or_im_impl::<false>(a)
    }

    fn at_fixed_im(&self, a: &Rational) -> (Polynomial<Integer>, Polynomial<Integer>) {
        self.at_fixed_re_or_im_impl::<true>(a)
    }

    //count how many complex roots are in the box a < re < b, c < im < d
    //or return None if there is a root on the boundary
    pub fn count_complex_roots(
        &self,
        a: &Rational,
        b: &Rational,
        c: &Rational,
        d: &Rational,
    ) -> Option<usize> {
        assert!(a < b);
        assert!(c < d);

        //the idea is to compute the winding number of the path around the boundary of the box
        //this is done by computing where the value of the polynomial crosses the real and imaginary axes as the input traces the path
        //the crossing points and their order is done using the exact total ordering of real polynomial roots
        let (a_vert_re, a_vert_im) = self.at_fixed_re(a);
        let (b_vert_re, b_vert_im) = self.at_fixed_re(b);
        let (c_horz_re, c_horz_im) = self.at_fixed_im(c);
        let (d_horz_re, d_horz_im) = self.at_fixed_im(d);

        // println!("poly = {} abcd = {} {} {} {}", &self, a, b, c, d);
        // println!(
        //     "a_vert_re = {}, a_vert_im = {}",
        //     Polynomial::to_string(&a_vert_re),
        //     Polynomial::to_string(&a_vert_im)
        // );
        // println!(
        //     "b_vert_re = {}, b_vert_im = {}",
        //     Polynomial::to_string(&b_vert_re),
        //     Polynomial::to_string(&b_vert_im)
        // );
        // println!(
        //     "c_horz_re = {}, c_horz_im = {}",
        //     Polynomial::to_string(&c_horz_re),
        //     Polynomial::to_string(&c_horz_im)
        // );
        // println!(
        //     "d_horz_re = {}, d_horz_im = {}",
        //     Polynomial::to_string(&d_horz_re),
        //     Polynomial::to_string(&d_horz_im)
        // );

        // //checks will fail - the real and imaginary parts are only up to scalar multiples
        // debug_assert_eq!(
        //     evaluate_at_rational(&a_vert_re, c),
        //     evaluate_at_rational(&c_horz_re, a)
        // );
        // debug_assert_eq!(
        //     evaluate_at_rational(&a_vert_re, d),
        //     evaluate_at_rational(&d_horz_re, a)
        // );
        // debug_assert_eq!(
        //     evaluate_at_rational(&b_vert_re, c),
        //     evaluate_at_rational(&c_horz_re, b)
        // );
        // debug_assert_eq!(
        //     evaluate_at_rational(&b_vert_re, d),
        //     evaluate_at_rational(&d_horz_re, b)
        // );
        // debug_assert_eq!(
        //     evaluate_at_rational(&a_vert_im, c),
        //     evaluate_at_rational(&c_horz_im, a)
        // );
        // debug_assert_eq!(
        //     evaluate_at_rational(&a_vert_im, d),
        //     evaluate_at_rational(&d_horz_im, a)
        // );
        // debug_assert_eq!(
        //     evaluate_at_rational(&b_vert_im, c),
        //     evaluate_at_rational(&c_horz_im, b)
        // );
        // debug_assert_eq!(
        //     evaluate_at_rational(&b_vert_im, d),
        //     evaluate_at_rational(&d_horz_im, b)
        // );

        //compute squarefree versions because when only care about the roots without multiplicity
        let a_vert_re_sqfr = a_vert_re.clone().primitive_squarefree_part();
        let a_vert_im_sqfr = a_vert_im.clone().primitive_squarefree_part();
        let b_vert_re_sqfr = b_vert_re.clone().primitive_squarefree_part();
        let b_vert_im_sqfr = b_vert_im.clone().primitive_squarefree_part();
        let c_horz_re_sqfr = c_horz_re.clone().primitive_squarefree_part();
        let c_horz_im_sqfr = c_horz_im.clone().primitive_squarefree_part();
        let d_horz_re_sqfr = d_horz_re.clone().primitive_squarefree_part();
        let d_horz_im_sqfr = d_horz_im.clone().primitive_squarefree_part();

        //trace an anticlockwise path around the box and create a list of crossings which encode what happens to the value of the polynomial
        #[derive(Debug)]
        enum Crossing {
            PosRe,
            PosIm,
            NegRe,
            NegIm,
        }

        fn crossings<const REVERSE: bool>(
            re: &Polynomial<Integer>,
            mut re_sqfr: Polynomial<Integer>,
            im: &Polynomial<Integer>,
            mut im_sqfr: Polynomial<Integer>,
            s: &Rational,
            t: &Rational,
        ) -> Option<Vec<Crossing>> {
            // println!(
            //     "REVERSE={} re={}, re_sqfr={}, im={}, im_sqfr={}",
            //     REVERSE,
            //     Polynomial::to_string(&re),
            //     Polynomial::to_string(&re_sqfr),
            //     Polynomial::to_string(&im),
            //     Polynomial::to_string(&im_sqfr)
            // );
            debug_assert_eq!(re == &Polynomial::zero(), re_sqfr == Polynomial::zero());
            debug_assert_eq!(im == &Polynomial::zero(), im_sqfr == Polynomial::zero());
            //because if the real and imaginary part are both constant at 0 then poly has infinitely many complex zeros which is not possible
            debug_assert!(re_sqfr != Polynomial::zero() || im_sqfr != Polynomial::zero());
            if re_sqfr == Polynomial::zero() {
                //the image is doing a path confied to the imaginary axis
                let roots_im = im.real_roots(Some(s), Some(t), true, true);
                if roots_im.len() == 0 {
                    //the image stays once side of the real axis
                    let val = evaluate_at_rational(im, s);
                    debug_assert_eq!(val > 0, evaluate_at_rational(im, t) > 0);
                    if val > 0 {
                        Some(vec![Crossing::PosIm]) //this whole line segment is a positive imaginary crossing
                    } else {
                        Some(vec![Crossing::NegIm]) //this whole line segment is a negative imaginary crossing
                    }
                } else {
                    //the image crosses the real axis and hence passes through 0
                    None
                }
            } else if im_sqfr == Polynomial::zero() {
                //the image is doing a path confied to the real axis
                let roots_re = re.real_roots(Some(s), Some(t), true, true);
                if roots_re.len() == 0 {
                    //the image stays one side of the imaginary axis
                    let val = evaluate_at_rational(re, s);
                    debug_assert_eq!(val > 0, evaluate_at_rational(re, t) > 0);
                    if val > 0 {
                        Some(vec![Crossing::PosRe]) //this whole line segment is a positive real crossing
                    } else {
                        Some(vec![Crossing::NegRe]) //this whole line segment is a negative real crossing
                    }
                } else {
                    //the image crosses the imaginary axis and hence passes through 0
                    None
                }
            } else {
                //want to isolate roots of squarefree polynomials without factoring
                //get ordered real roots in some structure
                //get ordered imag roots in some structure
                //do a merge sort pass to interleave the real and imaginary roots in the correct order
                //    if a real root equals an imaginary root then there is a root on the boundary
                //for each root of one type, compute the sign of the other part when evaluated at the root

                let mut crossings = vec![];

                //check the value of the real and imaginary part at the vertex at the start of this path
                let v = {
                    match REVERSE {
                        false => s,
                        true => t,
                    }
                };
                match (
                    evaluate_at_rational(re, v).cmp(&Rational::from(0)),
                    evaluate_at_rational(im, v).cmp(&Rational::from(0)),
                ) {
                    (std::cmp::Ordering::Equal, std::cmp::Ordering::Equal) => {
                        //the polynomial is zero at vertex v
                        return None;
                    }
                    (std::cmp::Ordering::Equal, std::cmp::Ordering::Less) => {
                        crossings.push(Crossing::NegIm);
                    }
                    (std::cmp::Ordering::Equal, std::cmp::Ordering::Greater) => {
                        crossings.push(Crossing::PosIm);
                    }
                    (std::cmp::Ordering::Less, std::cmp::Ordering::Equal) => {
                        crossings.push(Crossing::NegRe);
                    }
                    (std::cmp::Ordering::Greater, std::cmp::Ordering::Equal) => {
                        crossings.push(Crossing::PosRe);
                    }
                    (_, _) => {}
                }

                if evaluate_at_rational(&re_sqfr, s) == Rational::from(0) {
                    re_sqfr = Polynomial::div(
                        &re_sqfr,
                        &Polynomial::from_coeffs(vec![
                            -Rational::numerator(s),
                            Rational::denominator(s),
                        ]),
                    )
                    .unwrap();
                }

                if evaluate_at_rational(&re_sqfr, t) == Rational::from(0) {
                    re_sqfr = Polynomial::div(
                        &re_sqfr,
                        &Polynomial::from_coeffs(vec![
                            -Rational::numerator(t),
                            Rational::denominator(t),
                        ]),
                    )
                    .unwrap();
                }

                if evaluate_at_rational(&im_sqfr, s) == Rational::from(0) {
                    im_sqfr = Polynomial::div(
                        &im_sqfr,
                        &Polynomial::from_coeffs(vec![
                            -Rational::numerator(s),
                            Rational::denominator(s),
                        ]),
                    )
                    .unwrap();
                }
                if evaluate_at_rational(&im_sqfr, t) == Rational::from(0) {
                    im_sqfr = Polynomial::div(
                        &im_sqfr,
                        &Polynomial::from_coeffs(vec![
                            -Rational::numerator(t),
                            Rational::denominator(t),
                        ]),
                    )
                    .unwrap();
                }
                debug_assert_ne!(evaluate_at_rational(&re_sqfr, s), Rational::from(0));
                debug_assert_ne!(evaluate_at_rational(&re_sqfr, t), Rational::from(0));
                debug_assert_ne!(evaluate_at_rational(&im_sqfr, s), Rational::from(0));
                debug_assert_ne!(evaluate_at_rational(&im_sqfr, t), Rational::from(0));

                let mut re_roots =
                    re_sqfr
                        .clone()
                        .real_roots_squarefree(Some(s), Some(t), REVERSE, !REVERSE);
                let mut im_roots =
                    im_sqfr
                        .clone()
                        .real_roots_squarefree(Some(s), Some(t), REVERSE, !REVERSE);

                // println!("re_roots = {:?}", re_roots);
                // println!("im_roots = {:?}", im_roots);

                debug_assert!(re_roots.check_invariants().is_ok());
                debug_assert!(im_roots.check_invariants().is_ok());

                match SquarefreePolyRealRoots::separate(&mut re_roots, &mut im_roots) {
                    Ok(all_roots) => {
                        //the isolating intervals for re_roots and im_roots no longer overlap
                        //we can use this to our advantage...

                        for (interleave, root_idx) in {
                            match REVERSE {
                                false => all_roots,
                                true => all_roots.into_iter().rev().collect(),
                            }
                        } {
                            // println!("interleave = {:?} root_idx = {:?}", interleave, root_idx);
                            match interleave {
                                Interleave::First => {
                                    //a real root
                                    loop {
                                        let re_root = &re_roots.intervals[root_idx];
                                        match re_root {
                                            SquarefreePolyRealRootInterval::Rational(x) => {
                                                match evaluate_at_rational(&im, x)
                                                    .cmp(&Rational::from(0))
                                                {
                                                    std::cmp::Ordering::Less => {
                                                        crossings.push(Crossing::NegIm);
                                                        break;
                                                    }
                                                    std::cmp::Ordering::Equal => panic!(),
                                                    std::cmp::Ordering::Greater => {
                                                        crossings.push(Crossing::PosIm);
                                                        break;
                                                    }
                                                }
                                            }
                                            SquarefreePolyRealRootInterval::Real(a, b, _) => {
                                                match evaluate_at_rational(&im, a)
                                                    .cmp(&Rational::from(0))
                                                {
                                                    std::cmp::Ordering::Less => {
                                                        crossings.push(Crossing::NegIm);
                                                        break;
                                                    }
                                                    std::cmp::Ordering::Equal => {
                                                        //need to refine
                                                    }
                                                    std::cmp::Ordering::Greater => {
                                                        crossings.push(Crossing::PosIm);
                                                        break;
                                                    }
                                                }
                                                match evaluate_at_rational(&im, b)
                                                    .cmp(&Rational::from(0))
                                                {
                                                    std::cmp::Ordering::Less => {
                                                        crossings.push(Crossing::NegIm);
                                                        break;
                                                    }
                                                    std::cmp::Ordering::Equal => {
                                                        //need to refine
                                                    }
                                                    std::cmp::Ordering::Greater => {
                                                        crossings.push(Crossing::PosIm);
                                                        break;
                                                    }
                                                }
                                            }
                                        }
                                        re_roots.refine(root_idx);
                                    }
                                }
                                Interleave::Second => {
                                    //an imaginary root
                                    loop {
                                        let im_root = &im_roots.intervals[root_idx];
                                        match im_root {
                                            SquarefreePolyRealRootInterval::Rational(x) => {
                                                match evaluate_at_rational(&re, x)
                                                    .cmp(&Rational::from(0))
                                                {
                                                    std::cmp::Ordering::Less => {
                                                        crossings.push(Crossing::NegRe);
                                                        break;
                                                    }
                                                    std::cmp::Ordering::Equal => panic!(),
                                                    std::cmp::Ordering::Greater => {
                                                        crossings.push(Crossing::PosRe);
                                                        break;
                                                    }
                                                }
                                            }
                                            SquarefreePolyRealRootInterval::Real(a, b, _) => {
                                                match evaluate_at_rational(&re, a)
                                                    .cmp(&Rational::from(0))
                                                {
                                                    std::cmp::Ordering::Less => {
                                                        crossings.push(Crossing::NegRe);
                                                        break;
                                                    }
                                                    std::cmp::Ordering::Equal => {
                                                        //need to refine
                                                    }
                                                    std::cmp::Ordering::Greater => {
                                                        crossings.push(Crossing::PosRe);
                                                        break;
                                                    }
                                                }
                                                match evaluate_at_rational(&re, b)
                                                    .cmp(&Rational::from(0))
                                                {
                                                    std::cmp::Ordering::Less => {
                                                        crossings.push(Crossing::NegRe);
                                                        break;
                                                    }
                                                    std::cmp::Ordering::Equal => {
                                                        //need to refine
                                                    }
                                                    std::cmp::Ordering::Greater => {
                                                        crossings.push(Crossing::PosRe);
                                                        break;
                                                    }
                                                }
                                            }
                                        }
                                        im_roots.refine(root_idx);
                                    }
                                }
                            }
                        }
                        Some(crossings)
                    }
                    Err(()) => None,
                }
            }
        }

        /*
             a           b

         d   +-----------+
             |           |
             |           |
             |           |
             |           |
             |           |
         c   +-----------+

        */

        // println!("c = {:?}", crossings::<false>(&c_horz_re, c_horz_re_sqfr.clone(), &c_horz_im, c_horz_im_sqfr.clone(), a, b));
        // println!("b = {:?}", crossings::<false>(&b_vert_re, b_vert_re_sqfr.clone(), &b_vert_im, b_vert_im_sqfr.clone(), c, d));
        // println!("d = {:?}", crossings::<true>(&d_horz_re, d_horz_re_sqfr.clone(), &d_horz_im, d_horz_im_sqfr.clone(), a, b));
        // println!("a = {:?}", crossings::<true>(&a_vert_re, a_vert_re_sqfr.clone(), &a_vert_im, a_vert_im_sqfr.clone(), c, d));

        let mut winding = vec![];
        for cr in vec![
            crossings::<false>(&c_horz_re, c_horz_re_sqfr, &c_horz_im, c_horz_im_sqfr, a, b),
            crossings::<false>(&b_vert_re, b_vert_re_sqfr, &b_vert_im, b_vert_im_sqfr, c, d),
            crossings::<true>(&d_horz_re, d_horz_re_sqfr, &d_horz_im, d_horz_im_sqfr, a, b),
            crossings::<true>(&a_vert_re, a_vert_re_sqfr, &a_vert_im, a_vert_im_sqfr, c, d),
        ] {
            match cr {
                Some(mut w) => winding.append(&mut w),
                None => {
                    return None;
                }
            }
        }

        // println!("winding = {:?}", winding);

        //compute the winding number = number of roots
        if winding.len() <= 0 {
            Some(0)
        } else {
            fn axis_pair_to_num_offset(ax1: &Crossing, ax2: &Crossing) -> isize {
                match (ax1, ax2) {
                    (Crossing::PosRe, Crossing::PosRe) => 0,
                    (Crossing::PosRe, Crossing::PosIm) => 1,
                    (Crossing::PosRe, Crossing::NegRe) => panic!(),
                    (Crossing::PosRe, Crossing::NegIm) => -1,
                    (Crossing::PosIm, Crossing::PosRe) => -1,
                    (Crossing::PosIm, Crossing::PosIm) => 0,
                    (Crossing::PosIm, Crossing::NegRe) => 1,
                    (Crossing::PosIm, Crossing::NegIm) => panic!(),
                    (Crossing::NegRe, Crossing::PosRe) => panic!(),
                    (Crossing::NegRe, Crossing::PosIm) => -1,
                    (Crossing::NegRe, Crossing::NegRe) => 0,
                    (Crossing::NegRe, Crossing::NegIm) => 1,
                    (Crossing::NegIm, Crossing::PosRe) => 1,
                    (Crossing::NegIm, Crossing::PosIm) => panic!(),
                    (Crossing::NegIm, Crossing::NegRe) => -1,
                    (Crossing::NegIm, Crossing::NegIm) => 0,
                }
            }

            let mut num: isize = 0;
            num += axis_pair_to_num_offset(&winding[winding.len() - 1], &winding[0]);
            for i in 0..winding.len() - 1 {
                num += axis_pair_to_num_offset(&winding[i], &winding[i + 1]);
            }

            if num < 0 {
                panic!("winding should always be overall anti-clockwise");
            }
            let num = num as usize;
            match num % 4 {
                0 => Some(num / 4),
                _ => panic!("invalid remainder modulo four"),
            }
        }
    }

    fn uhp_complex_roots_irreducible_impl(
        &self,
        num_real_roots: usize,
    ) -> Vec<ComplexAlgebraicRoot> {
        debug_assert!(self.is_irreducible());
        debug_assert_eq!(self.all_real_roots().len(), num_real_roots);
        let deg = self.degree().unwrap();
        debug_assert!(num_real_roots <= deg);
        if num_real_roots == deg {
            vec![]
        } else {
            //search the upper half plane for the complete roots with positive imaginary part
            debug_assert_eq!((deg - num_real_roots) % 2, 0);
            let target_uhp_num = (deg - num_real_roots) / 2;

            let mut a = Rational::from(-1);
            let mut b = Rational::from(1);
            let mut c = Rational::from_signeds(1, 2);
            let mut d = Rational::from(2);

            loop {
                match self.count_complex_roots(&a, &b, &c, &d) {
                    Some(n) => {
                        debug_assert!(n <= target_uhp_num);
                        if n == target_uhp_num {
                            break;
                        }
                    }
                    None => {
                        //boundary root
                    }
                }
                a *= Rational::from(2);
                b *= Rational::from(2);
                c *= Rational::from_signeds(1, 2);
                d *= Rational::from(2);
            }

            fn bisect(
                poly: &Polynomial<Integer>,
                n: usize,
                a: &Rational,
                b: &Rational,
                c: &Rational,
                d: &Rational,
            ) -> Vec<ComplexAlgebraicRoot> {
                debug_assert!(a < b);
                debug_assert!(c < d);
                debug_assert_eq!(poly.count_complex_roots(&a, &b, &c, &d).unwrap(), n);
                if n == 0 {
                    vec![]
                } else if n == 1 {
                    vec![ComplexAlgebraicRoot {
                        poly: poly.clone(),
                        tight_a: a.clone(),
                        tight_b: b.clone(),
                        tight_c: c.clone(),
                        tight_d: d.clone(),
                    }]
                } else {
                    let ((n1, a1, b1, c1, d1), (n2, a2, b2, c2, d2)) =
                        bisect_box(poly, n, a, b, c, d);

                    let mut roots = bisect(poly, n1, &a1, &b1, &c1, &d1);
                    roots.append(&mut bisect(poly, n2, &a2, &b2, &c2, &d2));
                    return roots;
                }
            }

            bisect(self, target_uhp_num, &a, &b, &c, &d)
        }
    }

    fn lhp_complex_roots_irreducible_impl(
        &self,
        num_real_roots: usize,
    ) -> Vec<ComplexAlgebraicRoot> {
        self.uhp_complex_roots_irreducible_impl(num_real_roots)
            .into_iter()
            .map(|root| root.conj())
            .collect()
    }

    pub fn all_complex_roots_irreducible(&self) -> Vec<ComplexAlgebraic> {
        debug_assert!(self.is_irreducible());
        let deg = self.degree().unwrap();

        let mut all_roots = vec![];
        for real_root in self.all_real_roots() {
            all_roots.push(ComplexAlgebraic::Real(real_root));
        }
        let num_real_roots = all_roots.len();

        debug_assert!(num_real_roots <= deg);
        if num_real_roots == deg {
            return all_roots;
        }

        for complex_root in self.uhp_complex_roots_irreducible_impl(num_real_roots) {
            all_roots.push(ComplexAlgebraic::Complex(complex_root.clone().conj()));
            all_roots.push(ComplexAlgebraic::Complex(complex_root));
        }

        debug_assert_eq!(all_roots.len(), deg);
        all_roots
    }

    pub fn all_complex_roots(&self) -> Vec<ComplexAlgebraic> {
        assert_ne!(self, &Self::zero());
        let factors = self.factor().unwrap();
        let mut roots = vec![];
        for (factor, k) in factors.factors() {
            for root in factor.all_complex_roots_irreducible() {
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

#[derive(Debug, Clone)]
pub struct RealAlgebraicRoot {
    poly: Polynomial<Integer>, //a primitive irreducible polynomial of degree >= 2 with a unique real root between a and b
    //an arbitrarily small interval containing the root. May be mutated
    tight_a: Rational, //tight lower bound
    tight_b: Rational, //tight upper bound
    //a heuristically large interval containing the root. Should not shrink
    wide_a: LowerBound, //wide lower bound. None means -inf
    wide_b: UpperBound, //wide upper bound. None means +inf
    //false : decreasing i.e. poly(a) > poly(b), true : increasing i.e. poly(a) < poly(b)
    dir: bool,
}

impl Display for RealAlgebraicRoot {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.poly.num_coeffs() == 3 {
            //quadratic
            let a = self.poly.coeff(2);
            let b = self.poly.coeff(1);
            let c = self.poly.coeff(0);
            debug_assert!(a > Integer::ZERO);

            let d = &b * &b - Integer::from(4) * &a * &c;
            let mut d_sq = Integer::ONE;
            let mut d_sqfreee = Integer::ONE;
            let (d_sign, d_factors) = d.factor().unwrap().unit_and_factors();
            for (d_factor, k) in d_factors {
                d_sq *= d_factor.nat_pow(&(&k / Natural::TWO));
                if k % Natural::TWO == Natural::ONE {
                    d_sqfreee *= d_factor;
                }
            }
            debug_assert_eq!(d_sign, Integer::ONE); //because we are a real number
            debug_assert_eq!(d, &d_sqfreee * &d_sq * &d_sq);

            let two_a = Integer::TWO * a;

            let x = Rational::from_integers(-b, two_a.clone());
            let y = Rational::from_integers(d_sq, two_a);
            debug_assert!(y > Rational::ZERO);
            let r = d_sqfreee;

            let mut tight_a_abs = self.tight_a.clone();
            if tight_a_abs < Rational::ZERO {
                tight_a_abs = -tight_a_abs;
            }

            let mut tight_b_abs = self.tight_b.clone();
            if tight_b_abs < Rational::ZERO {
                tight_b_abs = -tight_b_abs;
            }

            let sign = tight_a_abs < tight_b_abs;
            if x == Rational::ZERO {
                if y == Rational::ONE {
                    write!(
                        f,
                        "{}√{}",
                        match sign {
                            true => "",
                            false => "-",
                        },
                        r
                    );
                } else {
                    write!(
                        f,
                        "{}{}√{}",
                        match sign {
                            true => "",
                            false => "-",
                        },
                        y,
                        r
                    );
                }
            } else {
                if y == Rational::ONE {
                    write!(
                        f,
                        "{}{}√{}",
                        x,
                        match sign {
                            true => "+",
                            false => "-",
                        },
                        r
                    );
                } else {
                    write!(
                        f,
                        "{}{}{}√{}",
                        x,
                        match sign {
                            true => "+",
                            false => "-",
                        },
                        y,
                        r
                    );
                }
            }
        } else {
            let mut root = self.clone();
            root.refine_to_accuracy(&Rational::from_integers(
                Integer::from(1),
                Integer::from(100),
            ));
            let m = (&root.tight_a + &root.tight_b) / Rational::TWO;

            write!(f, "≈");
            write!(f, "{}", rat_to_string(m));
            // write!(f, "±");
            // write!(f, "{}", rat_to_string(self.accuracy() / Rational::TWO));
        }
        Ok(())
    }
}

impl RealAlgebraicRoot {
    pub fn check_invariants(&self) -> Result<(), &'static str> {
        if !(self.tight_a < self.tight_b) {
            return Err("tight a should be strictly less than b");
        }
        if !(self.wide_a.clone() < self.wide_b.clone()) {
            return Err("wide a should be strictly less than b");
        }
        if self.poly
            != self
                .poly
                .clone()
                .primitive_squarefree_part()
                .factor_fav_assoc()
                .1
        {
            return Err("poly should be primitive and favoriate associate");
        }
        if !self.poly.is_irreducible() {
            return Err("poly should be irreducible");
        }
        if self.poly.degree().unwrap() < 2 {
            return Err("poly should have degree at least 2");
        }
        let at_a = self.evaluate(&self.tight_a);
        let at_b = self.evaluate(&self.tight_b);
        assert_ne!(at_a, Rational::from(0));
        assert_ne!(at_b, Rational::from(0));
        let sign_a = &at_a > &Rational::from(0);
        let sign_b = &at_b > &Rational::from(0);
        if sign_a == sign_b {
            return Err("sign at a and b should be different");
        }
        if self.dir != (sign_a == false) {
            return Err("dir is incorrect");
        }
        Ok(())
    }

    fn new_wide_bounds(poly: Polynomial<Integer>, wide_a: Rational, wide_b: Rational) -> Self {
        let dir = poly.apply_map(|x| Rational::from(x)).evaluate(&wide_a) < Rational::from(0);
        let x = Self {
            poly,
            tight_a: wide_a.clone(),
            tight_b: wide_b.clone(),
            wide_a: LowerBound::Finite(wide_a),
            wide_b: UpperBound::Finite(wide_b),
            dir,
        };
        debug_assert!(x.check_invariants().is_ok());
        x
    }

    fn evaluate(&self, val: &Rational) -> Rational {
        evaluate_at_rational(&self.poly, val)
    }

    pub fn accuracy(&self) -> Rational {
        &self.tight_b - &self.tight_a
    }

    pub fn refine(&mut self) {
        let m = RationalSimpleBetweenGenerator::new_between_middle(
            &self.tight_a,
            &self.tight_b,
            &Rational::from_integers(Integer::from(1), Integer::from(3)),
        )
        .next()
        .unwrap();
        let m_sign = self.evaluate(&m) > Rational::from(0);
        match self.dir == m_sign {
            true => {
                self.tight_b = m;
            }
            false => {
                self.tight_a = m;
            }
        }
    }

    pub fn refine_to_accuracy(&mut self, accuracy: &Rational) {
        while &self.accuracy() > accuracy {
            self.refine();
        }
    }

    pub fn cmp_mut(&mut self, other: &mut Self) -> std::cmp::Ordering {
        let polys_are_eq = self.poly == other.poly; //polys should be irreducible primitive fav-assoc so this is valid
        loop {
            //test for equality: If the tight bounds on one are within the wide bounds of the other
            if polys_are_eq {
                if other.wide_a <= self.tight_a && self.tight_b <= other.wide_b {
                    return std::cmp::Ordering::Equal;
                }
                if self.wide_a <= other.tight_a && other.tight_b <= self.wide_b {
                    return std::cmp::Ordering::Equal;
                }
            }

            //test for inequality: If the tight bounds are disjoint
            if self.tight_b <= other.tight_a {
                return std::cmp::Ordering::Less;
            }
            if other.tight_b <= self.tight_a {
                return std::cmp::Ordering::Greater;
            }

            //refine
            self.refine();
            other.refine();
        }
    }

    pub fn cmp_rat_mut(&mut self, other: &Rational) -> std::cmp::Ordering {
        loop {
            // println!("cmp_rat_mut {:?}", self);

            //test for inequality: other is outside the tight bounds
            if &self.tight_b <= other {
                return std::cmp::Ordering::Less;
            }
            if other <= &self.tight_a {
                return std::cmp::Ordering::Greater;
            }

            // println!("refine");

            //refine
            self.refine();
        }
    }

    fn neg_mut(&mut self) {
        let (unit, fav_assoc) = Polynomial::compose(
            &self.poly,
            &Polynomial::from_coeffs(vec![Integer::from(0), Integer::from(-1)]),
        )
        .factor_fav_assoc();
        if unit == Polynomial::one() {
            self.poly = fav_assoc;
            self.dir = !self.dir;
        } else if unit == Polynomial::neg(&Polynomial::one()) {
            self.poly = fav_assoc;
        } else {
            panic!();
        }
        (self.tight_a, self.tight_b) = (-self.tight_b.clone(), -self.tight_a.clone());
        (self.wide_a, self.wide_b) = (self.wide_b.clone().neg(), self.wide_a.clone().neg());
    }

    fn neg(mut self) -> Self {
        self.neg_mut();
        self
    }

    pub fn min_poly(&self) -> Polynomial<Rational> {
        self.poly.apply_map(|c| Rational::from(c)).fav_assoc()
    }
}

#[derive(Debug, Clone)]
pub struct ComplexAlgebraicRoot {
    tight_a: Rational, //tight lower bound for the real part
    tight_b: Rational, //tight upper bound for the real part
    tight_c: Rational, //tight lower bound for the imaginary part
    tight_d: Rational, //tight upper bound for the imaginary part

    poly: Polynomial<Integer>, //a primitive irreducible polynomial of degree >= 2 with a unique non-real complex root in the box defined by (a, b, c, d)
}

impl Display for ComplexAlgebraicRoot {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.poly.num_coeffs() == 3 {
            //quadratic
            let a = self.poly.coeff(2);
            let b = self.poly.coeff(1);
            let c = self.poly.coeff(0);
            debug_assert!(a > Integer::ZERO);

            let d = &b * &b - Integer::from(4) * &a * &c;
            let mut d_sq = Integer::ONE;
            let mut d_sqfreee = Integer::ONE;
            let (d_sign, d_factors) = d.factor().unwrap().unit_and_factors();
            for (d_factor, k) in d_factors {
                d_sq *= d_factor.nat_pow(&(&k / Natural::TWO));
                if k % Natural::TWO == Natural::ONE {
                    d_sqfreee *= d_factor;
                }
            }
            debug_assert_eq!(d_sign, -Integer::ONE); //because we are a real number
            debug_assert_eq!(-d, &d_sqfreee * &d_sq * &d_sq);

            let two_a = Integer::TWO * a;

            let x = Rational::from_integers(-b, two_a.clone());
            let y = Rational::from_integers(d_sq, two_a);
            debug_assert!(y > Rational::ZERO);
            let r = d_sqfreee;
            let r_str = {
                if r == Integer::ONE {
                    String::from("i")
                } else {
                    String::from("i√") + &r.to_string()
                }
            };

            let mut tight_c_abs = self.tight_c.clone();
            if tight_c_abs < Rational::ZERO {
                tight_c_abs = -tight_c_abs;
            }

            let mut tight_d_abs = self.tight_d.clone();
            if tight_d_abs < Rational::ZERO {
                tight_d_abs = -tight_d_abs;
            }

            let sign = tight_c_abs < tight_d_abs;
            if x == Rational::ZERO {
                if y == Rational::ONE {
                    write!(
                        f,
                        "{}{}",
                        match sign {
                            true => "",
                            false => "-",
                        },
                        r_str
                    );
                } else {
                    write!(
                        f,
                        "{}{}{}",
                        match sign {
                            true => "",
                            false => "-",
                        },
                        y,
                        r_str
                    );
                }
            } else {
                if y == Rational::ONE {
                    write!(
                        f,
                        "{}{}{}",
                        x,
                        match sign {
                            true => "+",
                            false => "-",
                        },
                        r_str
                    );
                } else {
                    write!(
                        f,
                        "{}{}{}{}",
                        x,
                        match sign {
                            true => "+",
                            false => "-",
                        },
                        y,
                        r_str
                    );
                }
            }
        } else {
            let mut root = self.clone();
            root.refine_to_accuracy(&Rational::from_integers(
                Integer::from(1),
                Integer::from(100),
            ));

            let m_re = (&root.tight_a + &root.tight_b) / Rational::TWO;
            let m_im = (&root.tight_c + &root.tight_d) / Rational::TWO;

            write!(f, "≈");
            write!(f, "{}", rat_to_string(m_re));
            // write!(f, "±");
            // write!(f, "{}", rat_to_string(self.accuracy_re() / Rational::TWO));
            if m_im >= 0 {
                write!(f, "+");
            }
            write!(f, "{}", rat_to_string(m_im));
            // write!(f, "±");
            // write!(f, "{}", rat_to_string(self.accuracy_im() / Rational::TWO));
            write!(f, "i");
        }
        Ok(())
    }
}

impl ComplexAlgebraicRoot {
    pub fn check_invariants(&self) -> Result<(), &'static str> {
        if !(self.tight_a < self.tight_b) {
            return Err("tight a should be strictly less than b");
        }
        if !(self.tight_c < self.tight_d) {
            return Err("tight c should be strictly less than d");
        }
        // if !(self.wide_a < self.wide_b) {
        //     return Err("wide a should be strictly less than b");
        // }
        // if !(self.wide_c < self.wide_d) {
        //     return Err("wide c should be strictly less than d");
        // }

        if !self.poly.is_irreducible() {
            return Err("Isolated complex root minimal polynomial should be irreducible");
        }

        if self.poly.degree().unwrap() < 2 {
            return Err("Isolated complex root minimal polynomial should have degree at least 2");
        }

        match self.poly.count_complex_roots(
            &self.tight_a,
            &self.tight_b,
            &self.tight_c,
            &self.tight_d,
        ) {
            Some(1) => {}
            Some(_) => {
                return Err("Isolated complex root must exactly 1 root with none on the boundary");
            }
            None => {
                return Err(
                    "Isolated complex root must contain exactly 1 root with none on the boundary",
                );
            }
        }

        let real_roots_in_box = if self.tight_c < Rational::ZERO && Rational::ZERO < self.tight_d {
            self.poly
                .all_real_roots()
                .into_iter()
                .filter(|x| {
                    &RealAlgebraic::Rational(self.tight_a.clone()) < x
                        && x < &RealAlgebraic::Rational(self.tight_b.clone())
                })
                .collect_vec()
        } else {
            vec![]
        };
        if !real_roots_in_box.is_empty() {
            return Err("Isolated complex root must not be a real root");
        }

        Ok(())
    }

    pub fn accuracy_re(&self) -> Rational {
        &self.tight_b - &self.tight_a
    }

    pub fn accuracy_im(&self) -> Rational {
        &self.tight_d - &self.tight_c
    }

    fn conj(mut self) -> Self {
        (self.tight_c, self.tight_d) = (-self.tight_d, -self.tight_c);
        self
    }

    pub fn neg_mut(&mut self) {
        self.poly = Polynomial::compose(
            &self.poly,
            &Polynomial::from_coeffs(vec![Integer::from(0), Integer::from(-1)]),
        )
        .fav_assoc();
        (self.tight_a, self.tight_b) = (-self.tight_b.clone(), -self.tight_a.clone());
        (self.tight_c, self.tight_d) = (-self.tight_d.clone(), -self.tight_c.clone());
    }

    pub fn neg(mut self) -> Self {
        self.neg_mut();
        self
    }

    pub fn refine(&mut self) {
        let ((n1, a1, b1, c1, d1), (n2, a2, b2, c2, d2)) = bisect_box(
            &self.poly,
            1,
            &self.tight_a,
            &self.tight_b,
            &self.tight_c,
            &self.tight_d,
        );

        match (n1, n2) {
            (1, 0) => {
                self.tight_a = a1;
                self.tight_b = b1;
                self.tight_c = c1;
                self.tight_d = d1;
            }
            (0, 1) => {
                self.tight_a = a2;
                self.tight_b = b2;
                self.tight_c = c2;
                self.tight_d = d2;
            }
            _ => {
                unreachable!();
            }
        }

        #[cfg(debug_assertions)]
        self.check_invariants().unwrap();
    }

    pub fn refine_to_accuracy(&mut self, accuracy: &Rational) {
        while &self.accuracy_re() > accuracy || &self.accuracy_im() > accuracy {
            self.refine();
        }
    }

    pub fn eq_mut(&mut self, other: &mut Self) -> bool {
        let poly = &self.poly;
        //polys should be irreducible primitive fav-assoc so this is valid
        if poly == &other.poly {
            let mut uhp_roots = poly.uhp_complex_roots_irreducible_impl(
                poly.real_roots_irreducible(None, None, false, false).len(),
            );
            //want to identify which of the uhp_roots or their conjugates self and other are
            //use false/true for uph/lhp and use usize for the index in uhp_roots
            let mut identify = |root: &mut Self| -> (bool, usize) {
                let mut possible = vec![];
                for i in 0..uhp_roots.len() {
                    possible.push((false, i));
                    possible.push((true, i));
                }

                while possible.len() > 1 {
                    possible = possible
                        .into_iter()
                        .filter(|(conj, idx)| {
                            let (possible_root_tight_a, possible_root_tight_b) =
                                (&uhp_roots[*idx].tight_a, &uhp_roots[*idx].tight_b);
                            let (possible_root_tight_c, possible_root_tight_d) = match conj {
                                false => (
                                    uhp_roots[*idx].tight_c.clone(),
                                    uhp_roots[*idx].tight_d.clone(),
                                ),
                                true => (-&uhp_roots[*idx].tight_d, -&uhp_roots[*idx].tight_c),
                            };
                            //does the region for root overlap with region for possible_tight_root?
                            possible_root_tight_a < &root.tight_b
                                && &root.tight_a < possible_root_tight_b
                                && possible_root_tight_c < root.tight_d
                                && root.tight_c < possible_root_tight_d
                        })
                        .collect();

                    root.refine();
                    for uhp_root in &mut uhp_roots {
                        uhp_root.refine();
                    }
                }

                debug_assert_ne!(possible.len(), 0);
                possible[0]
            };
            identify(self) == identify(other)
        } else {
            false
        }
    }

    pub fn min_poly(&self) -> Polynomial<Rational> {
        self.poly.apply_map(|c| Rational::from(c)).fav_assoc()
    }

    // let polys_are_eq = self.poly == other.poly; //polys should be irreducible primitive fav-assoc so this is valid
    // loop {
    //     //test for equality: If the tight bounds on one are within the wide bounds of the other
    //     if polys_are_eq {
    //         if other.wide_a <= self.tight_a && self.tight_b <= other.wide_b {
    //             return std::cmp::Ordering::Equal;
    //         }
    //         if self.wide_a <= other.tight_a && other.tight_b <= self.wide_b {
    //             return std::cmp::Ordering::Equal;
    //         }
    //     }

    //     //test for inequality: If the tight bounds are disjoint
    //     if self.tight_b <= other.tight_a {
    //         return std::cmp::Ordering::Less;
    //     }
    //     if other.tight_b <= self.tight_a {
    //         return std::cmp::Ordering::Greater;
    //     }

    //     //refine
    //     self.refine();
    //     other.refine();
    // }
}

#[derive(Debug, Clone)]
pub enum RealAlgebraic {
    Rational(Rational),
    Real(RealAlgebraicRoot),
}

impl RealAlgebraic {
    pub fn check_invariants(&self) -> Result<(), &'static str> {
        match self {
            RealAlgebraic::Rational(_x) => {}
            RealAlgebraic::Real(x) => match x.check_invariants() {
                Ok(()) => {}
                Err(e) => {
                    return Err(e);
                }
            },
        }
        Ok(())
    }

    pub fn cmp_mut(&mut self, other: &mut Self) -> std::cmp::Ordering {
        {
            match self {
                RealAlgebraic::Rational(self_rep) => match other {
                    RealAlgebraic::Rational(other_rep) => self_rep.cmp(&other_rep),
                    RealAlgebraic::Real(other_rep) => other_rep.cmp_rat_mut(self_rep).reverse(),
                },
                RealAlgebraic::Real(self_rep) => match other {
                    RealAlgebraic::Rational(other_rep) => self_rep.cmp_rat_mut(other_rep),
                    RealAlgebraic::Real(other_rep) => self_rep.cmp_mut(other_rep),
                },
            }
        }
    }

    pub fn min_poly(&self) -> Polynomial<Rational> {
        match self {
            RealAlgebraic::Rational(rat) => Polynomial::from_coeffs(vec![-rat, Rational::ONE]),
            RealAlgebraic::Real(real_root) => real_root.min_poly(),
        }
    }
}

impl PositiveRealNthRootStructure for CannonicalStructure<RealAlgebraic> {
    fn nth_root(&self, x: &Self::Set, n: usize) -> Result<Self::Set, ()> {
        if n == 0 {
            panic!()
        } else if n == 1 {
            Ok(x.clone())
        } else {
            match x.cmp(&mut RealAlgebraic::zero()) {
                std::cmp::Ordering::Less => Err(()),
                std::cmp::Ordering::Equal => Ok(self.zero()),
                std::cmp::Ordering::Greater => {
                    let poly = match x {
                        RealAlgebraic::Rational(rat) => {
                            Polynomial::from_coeffs(vec![-rat.numerator(), rat.denominator()])
                        }
                        RealAlgebraic::Real(real) => real.poly.clone(),
                    };
                    let mut coeffs = vec![];
                    for (i, c) in poly.into_coeffs().into_iter().enumerate() {
                        if i != 0 {
                            for _j in 0..(n - 1) {
                                coeffs.push(Integer::ZERO);
                            }
                        }
                        coeffs.push(c)
                    }
                    let nthroot_poly = Polynomial::from_coeffs(coeffs).primitive_squarefree_part();
                    // println!("nthroot_poly = {:?}", nthroot_poly);
                    let possible_nthroots = nthroot_poly.all_real_roots_squarefree();
                    //search through the intervals starting from the largest.
                    //if it is positive, check that the nth power of the bounds surrounds self
                    //otherwise, the first one that does not have both bounds positive must be the nth root

                    // println!("{:?}", possible_nthroots);

                    for (idx, interval) in possible_nthroots.intervals.iter().enumerate().rev() {
                        match interval {
                            SquarefreePolyRealRootInterval::Rational(rat) => {
                                debug_assert!(rat > &Rational::ZERO);
                                match x {
                                    RealAlgebraic::Rational(self_rat) => {
                                        debug_assert_eq!(&rat.nat_pow(&Natural::from(n)), self_rat);
                                    }
                                    RealAlgebraic::Real(_) => debug_assert!(false),
                                }
                                return Ok(possible_nthroots.get_real_root(idx));
                            }
                            SquarefreePolyRealRootInterval::Real(a, b, dir) => {
                                debug_assert!(a < b);
                                if &Rational::ZERO < a {
                                    let a_pow = a.nat_pow(&Natural::from(n));
                                    let b_pow = b.nat_pow(&Natural::from(n));
                                    if &RealAlgebraic::Rational(a_pow) <= x
                                        && x <= &RealAlgebraic::Rational(b_pow)
                                    {
                                        return Ok(possible_nthroots.get_real_root(idx));
                                    }
                                } else {
                                    return Ok(possible_nthroots.get_real_root(idx));
                                }
                            }
                        }
                    }

                    unreachable!()
                }
            }
        }
    }
}

impl PartialEq for RealAlgebraic {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other).is_eq()
    }
}

impl Eq for RealAlgebraic {}

impl PartialOrd for RealAlgebraic {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.clone().cmp_mut(&mut other.clone()))
    }
}

impl Ord for RealAlgebraic {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other).unwrap()
    }
}

#[derive(Debug, Clone)]
// #[derive(Debug, Clone, PartialEq, Eq)] //todo: this
pub enum ComplexAlgebraic {
    Real(RealAlgebraic),
    Complex(ComplexAlgebraicRoot),
}

impl ComplexAlgebraic {
    pub fn i() -> Self {
        let i = ComplexAlgebraic::Complex(ComplexAlgebraicRoot {
            tight_a: Rational::from_integers(Integer::from(-1), Integer::from(1)),
            tight_b: Rational::from_integers(Integer::from(1), Integer::from(1)),
            tight_c: Rational::from_integers(Integer::from(0), Integer::from(1)),
            tight_d: Rational::from_integers(Integer::from(2), Integer::from(1)),
            poly: Polynomial::from_coeffs(vec![Integer::ONE, Integer::ZERO, Integer::ONE]),
        });
        #[cfg(debug_assertions)]
        i.check_invariants().unwrap();
        i
    }
}

impl ComplexAlgebraic {
    pub fn check_invariants(&self) -> Result<(), &'static str> {
        match self {
            ComplexAlgebraic::Real(x) => match x.check_invariants() {
                Ok(()) => {}
                Err(e) => {
                    return Err(e);
                }
            },
            ComplexAlgebraic::Complex(x) => match x.check_invariants() {
                Ok(()) => {}
                Err(e) => {
                    return Err(e);
                }
            },
        }
        Ok(())
    }

    pub fn eq_mut(&mut self, other: &mut Self) -> bool {
        match (self, other) {
            (ComplexAlgebraic::Real(a), ComplexAlgebraic::Real(b)) => {
                RealAlgebraic::cmp_mut(a, b).is_eq()
            }
            (ComplexAlgebraic::Real(a), ComplexAlgebraic::Complex(b)) => false,
            (ComplexAlgebraic::Complex(a), ComplexAlgebraic::Real(b)) => false,
            (ComplexAlgebraic::Complex(a), ComplexAlgebraic::Complex(b)) => a.eq_mut(b),
        }
    }
}

impl StructuredType for RealAlgebraic {
    type Structure = CannonicalStructure<RealAlgebraic>;

    fn structure() -> Rc<Self::Structure> {
        CannonicalStructure::new().into()
    }
}

impl Display for RealAlgebraic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RealAlgebraic::Rational(a) => write!(f, "{}", a),
            RealAlgebraic::Real(a) => write!(f, "{}", a),
        }
    }
}

impl EqualityStructure for CannonicalStructure<RealAlgebraic> {
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        a == b
    }
}

impl RingStructure for CannonicalStructure<RealAlgebraic> {
    fn zero(&self) -> Self::Set {
        RealAlgebraic::Rational(Rational::from(0))
    }

    fn one(&self) -> Self::Set {
        RealAlgebraic::Rational(Rational::from(1))
    }

    fn neg(&self, a: &Self::Set) -> Self::Set {
        match a {
            RealAlgebraic::Rational(a) => RealAlgebraic::Rational(-a),
            RealAlgebraic::Real(root) => RealAlgebraic::Real(root.clone().neg()),
        }
    }

    fn add(&self, alg1: &Self::Set, alg2: &Self::Set) -> Self::Set {
        // println!("add {:?} {:?}", alg1, alg2);

        fn add_rat(mut elem: RealAlgebraicRoot, rat: &Rational) -> RealAlgebraicRoot {
            elem.tight_a += rat;
            elem.tight_b += rat;
            match &elem.wide_a {
                LowerBound::Inf => {}
                LowerBound::Finite(a) => {
                    elem.wide_a = LowerBound::Finite(a + rat);
                }
            }
            match &elem.wide_b {
                UpperBound::Inf => {}
                UpperBound::Finite(b) => {
                    elem.wide_b = UpperBound::Finite(b + rat);
                }
            }

            //compose with x - n/d = dx - n
            elem.poly = Polynomial::compose(
                &elem.poly.apply_map(|c| Rational::from(c)),
                &Polynomial::from_coeffs(vec![-rat, Rational::ONE]),
            )
            .primitive_part_fof();

            #[cfg(debug_assertions)]
            elem.check_invariants().unwrap();
            elem
        }

        match (alg1, alg2) {
            (RealAlgebraic::Rational(rat1), RealAlgebraic::Rational(rat2)) => {
                RealAlgebraic::Rational(rat1 + rat2)
            }
            (RealAlgebraic::Rational(rat1), RealAlgebraic::Real(alg2)) => {
                RealAlgebraic::Real(add_rat(alg2.clone(), rat1))
            }
            (RealAlgebraic::Real(alg1), RealAlgebraic::Rational(rat2)) => {
                RealAlgebraic::Real(add_rat(alg1.clone(), rat2))
            }
            (RealAlgebraic::Real(alg1), RealAlgebraic::Real(alg2)) => {
                let mut alg1 = alg1.clone();
                let mut alg2 = alg2.clone();

                let factored_rsp = root_sum_poly(&alg1.poly, &alg2.poly).factor().unwrap();
                let polys: Vec<_> = factored_rsp.factors().iter().map(|(f, _k)| f).collect();
                //the sum of alg1 and alg2 is exactly one root of exactly one of the irreducible polynomials in polys
                //the task now is to refine alg1 and alg2 until the root is identified

                let mut root_groups: Vec<_> = polys
                    .into_iter()
                    .map(|p| p.all_real_roots_squarefree())
                    .collect();

                //store indicies of possible roots
                let mut possible = std::collections::HashSet::new();
                for i in 0..root_groups.len() {
                    for j in 0..root_groups[i].intervals.len() {
                        possible.insert((i, j));
                    }
                }

                while possible.len() > 1 {
                    let ans_tight_a = &alg1.tight_a + &alg2.tight_a;
                    let ans_tight_b = &alg1.tight_b + &alg2.tight_b;
                    //filter out roots which dont overlap with the known range for the sum root
                    possible = possible
                        .into_iter()
                        .filter(|(i, j)| match &root_groups[*i].intervals[*j] {
                            SquarefreePolyRealRootInterval::Rational(x) => {
                                &ans_tight_a < x && x < &ans_tight_b
                            }
                            SquarefreePolyRealRootInterval::Real(ta, tb, _dir) => {
                                ta < &ans_tight_b && &ans_tight_a < tb
                            }
                        })
                        .collect();

                    alg1.refine();
                    alg2.refine();
                    for (i, j) in &possible {
                        root_groups[*i].refine(*j);
                    }
                }
                assert_eq!(possible.len(), 1);
                let (i, j) = possible.into_iter().next().unwrap();
                let ans = root_groups
                    .into_iter()
                    .nth(i)
                    .unwrap()
                    .to_real_roots()
                    .into_iter()
                    .nth(j)
                    .unwrap();
                #[cfg(debug_assertions)]
                ans.check_invariants().unwrap();
                ans
            }
        }
    }

    fn mul(&self, elem1: &Self::Set, elem2: &Self::Set) -> Self::Set {
        match elem1.cmp(&self.zero()) {
            std::cmp::Ordering::Less => {
                return self.neg(&self.mul(&self.neg(elem1), elem2));
            }
            std::cmp::Ordering::Equal => {
                return self.zero();
            }
            std::cmp::Ordering::Greater => {}
        }

        match elem2.cmp(&self.zero()) {
            std::cmp::Ordering::Less => {
                return self.neg(&self.mul(elem1, &self.neg(elem2)));
            }
            std::cmp::Ordering::Equal => {
                return self.zero();
            }
            std::cmp::Ordering::Greater => {}
        }

        debug_assert!(elem1 > &self.zero());
        debug_assert!(elem2 > &self.zero());

        fn mul_pos_rat(mut elem: RealAlgebraicRoot, rat: &Rational) -> RealAlgebraicRoot {
            debug_assert!(rat > &Rational::from(0));
            elem.tight_a *= rat;
            elem.tight_b *= rat;
            match &elem.wide_a {
                LowerBound::Inf => {}
                LowerBound::Finite(a) => {
                    elem.wide_a = LowerBound::Finite(a * rat);
                }
            }
            match &elem.wide_b {
                UpperBound::Inf => {}
                UpperBound::Finite(b) => {
                    elem.wide_b = UpperBound::Finite(b * rat);
                }
            }
            elem.poly = root_pos_rat_mul_poly(elem.poly, rat);
            #[cfg(debug_assertions)]
            elem.check_invariants().unwrap();
            elem
        }

        match (elem1, elem2) {
            (RealAlgebraic::Rational(rat1), RealAlgebraic::Rational(rat2)) => {
                RealAlgebraic::Rational(rat1 * rat2)
            }
            (RealAlgebraic::Rational(rat1), RealAlgebraic::Real(alg2)) => {
                RealAlgebraic::Real(mul_pos_rat(alg2.clone(), rat1))
            }
            (RealAlgebraic::Real(alg1), RealAlgebraic::Rational(rat2)) => {
                RealAlgebraic::Real(mul_pos_rat(alg1.clone(), rat2))
            }
            (RealAlgebraic::Real(alg1), RealAlgebraic::Real(alg2)) => {
                let mut alg1 = alg1.clone();
                let mut alg2 = alg2.clone();

                let factored_rpp = root_prod_poly(&alg1.poly, &alg2.poly).factor().unwrap();
                let polys: Vec<_> = factored_rpp.factors().iter().map(|(f, _k)| f).collect();
                //the sum of alg1 and alg2 is exactly one root of exactly one of the irreducible polynomials in polys
                //the task now is to refine alg1 and alg2 until the root is identified

                let mut root_groups: Vec<_> = polys
                    .into_iter()
                    .map(|p| p.all_real_roots_squarefree())
                    .collect();

                //store indicies of possible roots
                let mut possible = std::collections::HashSet::new();
                for i in 0..root_groups.len() {
                    for j in 0..root_groups[i].intervals.len() {
                        possible.insert((i, j));
                    }
                }

                while possible.len() > 1 {
                    let ans_tight_a = &alg1.tight_a * &alg2.tight_a;
                    let ans_tight_b = &alg1.tight_b * &alg2.tight_b;
                    //filter out roots which dont overlap with the known range for the sum root
                    possible = possible
                        .into_iter()
                        .filter(|(i, j)| match &root_groups[*i].intervals[*j] {
                            SquarefreePolyRealRootInterval::Rational(x) => {
                                &ans_tight_a < x && x < &ans_tight_b
                            }
                            SquarefreePolyRealRootInterval::Real(ta, tb, _dir) => {
                                ta < &ans_tight_b && &ans_tight_a < tb
                            }
                        })
                        .collect();

                    alg1.refine();
                    alg2.refine();
                    for (i, j) in &possible {
                        root_groups[*i].refine(*j);
                    }
                }
                assert_eq!(possible.len(), 1);
                let (i, j) = possible.into_iter().next().unwrap();
                let ans = root_groups
                    .into_iter()
                    .nth(i)
                    .unwrap()
                    .to_real_roots()
                    .into_iter()
                    .nth(j)
                    .unwrap();
                #[cfg(debug_assertions)]
                ans.check_invariants().unwrap();
                ans
            }
        }
    }
}

impl IntegralDomainStructure for CannonicalStructure<RealAlgebraic> {
    fn div(&self, a: &Self::Set, b: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        match self.inv(b) {
            Ok(b_inv) => Ok(self.mul(a, &b_inv)),
            Err(err) => Err(err),
        }
    }

    fn inv(&self, a: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        let mut a = a.clone();
        match RealAlgebraic::cmp_mut(&mut a, &mut self.zero()) {
            std::cmp::Ordering::Less => match self.inv(&self.neg(&a)) {
                Ok(neg_elem_inv) => Ok(self.neg(&neg_elem_inv)),
                Err(err) => Err(err),
            },
            std::cmp::Ordering::Equal => Err(RingDivisionError::DivideByZero),
            std::cmp::Ordering::Greater => match a {
                RealAlgebraic::Rational(x) => {
                    Ok(RealAlgebraic::Rational(Rational::inv(&x).unwrap()))
                }
                RealAlgebraic::Real(mut root) => {
                    debug_assert!(root.tight_a >= Rational::from(0));
                    while root.tight_a == Rational::from(0) {
                        root.refine();
                    }
                    debug_assert!(Rational::from(0) < root.tight_a);
                    (root.tight_a, root.tight_b) = (
                        Rational::inv(&root.tight_b).unwrap(),
                        Rational::inv(&root.tight_a).unwrap(),
                    );
                    (root.wide_a, root.wide_b) = (
                        {
                            match root.wide_b {
                                UpperBound::Inf => LowerBound::Finite(Rational::from(0)),
                                UpperBound::Finite(x) => match Rational::inv(&x) {
                                    Ok(x_inv) => LowerBound::Finite(x_inv),
                                    Err(RingDivisionError::DivideByZero) => panic!("wide upper bound of strictly positive root should be strictly positive i.e. non-zero"),
                                    Err(_) => panic!()
                                },
                            }
                        },
                        {
                            match root.wide_a {
                                LowerBound::Inf => UpperBound::Inf,
                                LowerBound::Finite(x) => match x.cmp(&Rational::from(0)) {
                                    std::cmp::Ordering::Less => UpperBound::Inf,
                                    std::cmp::Ordering::Equal => UpperBound::Inf,
                                    std::cmp::Ordering::Greater => {
                                        UpperBound::Finite(Rational::inv(&x).unwrap())
                                    }
                                },
                            }
                        },
                    );
                    let (unit, fav_assoc) = Polynomial::from_coeffs(
                        root.poly.into_coeffs().into_iter().rev().collect(),
                    )
                    .factor_fav_assoc();
                    if unit == Polynomial::one() {
                        root.poly = fav_assoc;
                        root.dir = !root.dir;
                    } else if unit == Polynomial::neg(&Polynomial::one()) {
                        root.poly = fav_assoc;
                    } else {
                        panic!();
                    }
                    let ans = RealAlgebraic::Real(root);
                    #[cfg(debug_assertions)]
                    ans.check_invariants().unwrap();
                    Ok(ans)
                }
            },
        }
    }
}

impl FieldStructure for CannonicalStructure<RealAlgebraic> {}

impl CharZeroStructure for CannonicalStructure<RealAlgebraic> {}

impl ComplexSubsetStructure for CannonicalStructure<RealAlgebraic> {}

impl RealSubsetStructure for CannonicalStructure<RealAlgebraic> {}

impl RealToFloatStructure for CannonicalStructure<RealAlgebraic> {
    fn as_f64(&self, x: &Self::Set) -> f64 {
        todo!()
    }
}

impl RealRoundingStructure for CannonicalStructure<RealAlgebraic> {
    fn floor(&self, x: &Self::Set) -> Integer {
        todo!()
    }
    fn ceil(&self, x: &Self::Set) -> Integer {
        todo!()
    }
    fn round(&self, x: &Self::Set) -> Integer {
        todo!()
    }
}

impl PartialEq for ComplexAlgebraic {
    fn eq(&self, other: &Self) -> bool {
        Self::eq_mut(&mut self.clone(), &mut other.clone())
    }
}

impl Eq for ComplexAlgebraic {}

//box gen must yield boxes which converge to some root of the polynomial
fn identify_complex_root(
    poly: Polynomial<Integer>,
    mut box_gen: impl Iterator<Item = (Rational, Rational, Rational, Rational)>,
) -> ComplexAlgebraic {
    let poly = poly.primitive_squarefree_part();
    //find the irreducible factor poly which contains the root being converged to
    let (mut a, mut b, mut c, mut d) = box_gen.next().unwrap();

    let irr_poly = {
        let (unit, factors) = poly.factor().unwrap().unit_and_factors();
        let irr_polys = factors.into_iter().map(|(f, k)| f).collect_vec();
        let mut possible_irr_poly_idxs: HashSet<_> = (0..irr_polys.len()).collect();
        loop {
            debug_assert!(!possible_irr_poly_idxs.is_empty());
            possible_irr_poly_idxs = possible_irr_poly_idxs
                .into_iter()
                .filter(|idx| irr_polys[*idx].count_complex_roots(&a, &b, &c, &d) != Some(0))
                .collect();
            if possible_irr_poly_idxs.len() == 1 {
                break;
            }
            (a, b, c, d) = box_gen.next().unwrap();
        }
        debug_assert_eq!(possible_irr_poly_idxs.len(), 1);
        irr_polys
            .into_iter()
            .nth(possible_irr_poly_idxs.into_iter().next().unwrap())
            .unwrap()
    };

    let mut roots = irr_poly.all_complex_roots();
    let mut possible_roots: HashSet<_> = (0..roots.len()).collect();
    loop {
        debug_assert!(!possible_roots.is_empty());
        possible_roots = possible_roots
            .into_iter()
            .filter(|idx| match &roots[*idx] {
                ComplexAlgebraic::Real(RealAlgebraic::Rational(root)) => {
                    &a < root && root < &b && c < 0 && 0 < d
                }
                ComplexAlgebraic::Real(RealAlgebraic::Real(root)) => {
                    a < root.tight_b && root.tight_a < b && c < 0 && 0 < d
                }
                ComplexAlgebraic::Complex(root) => {
                    a < root.tight_b && root.tight_a < b && c < root.tight_d && root.tight_c < d
                }
            })
            .collect();
        if possible_roots.len() == 1 {
            break;
        }
        (a, b, c, d) = box_gen.next().unwrap();
        for idx in &possible_roots {
            match &mut roots[*idx] {
                ComplexAlgebraic::Real(RealAlgebraic::Rational(root)) => {}
                ComplexAlgebraic::Real(RealAlgebraic::Real(root)) => {
                    root.refine();
                }
                ComplexAlgebraic::Complex(root) => {
                    root.refine();
                }
            }
        }
    }
    debug_assert_eq!(possible_roots.len(), 1);
    let ans = roots
        .into_iter()
        .nth(possible_roots.into_iter().next().unwrap())
        .unwrap();
    #[cfg(debug_assertions)]
    ans.check_invariants().unwrap();
    ans
}

impl StructuredType for ComplexAlgebraic {
    type Structure = CannonicalStructure<ComplexAlgebraic>;

    fn structure() -> Rc<Self::Structure> {
        CannonicalStructure::new().into()
    }
}

impl Display for ComplexAlgebraic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ComplexAlgebraic::Real(a) => write!(f, "{}", a),
            ComplexAlgebraic::Complex(a) => write!(f, "{}", a),
        }
    }
}

impl EqualityStructure for CannonicalStructure<ComplexAlgebraic> {
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        a == b
    }
}

impl RingStructure for CannonicalStructure<ComplexAlgebraic> {
    fn zero(&self) -> Self::Set {
        ComplexAlgebraic::Real(RealAlgebraic::Rational(Rational::zero()))
    }

    fn one(&self) -> Self::Set {
        ComplexAlgebraic::Real(RealAlgebraic::Rational(Rational::one()))
    }

    fn neg(&self, a: &Self::Set) -> Self::Set {
        match a {
            ComplexAlgebraic::Real(root) => ComplexAlgebraic::Real(RealAlgebraic::neg(root)),
            ComplexAlgebraic::Complex(root) => ComplexAlgebraic::Complex(root.clone().neg()),
        }
    }

    fn add(&self, alg1: &Self::Set, alg2: &Self::Set) -> Self::Set {
        // println!("add {:?} {:?}", alg1, alg2);
        // alg1.check_invariants().unwrap();
        // alg2.check_invariants().unwrap();

        fn add_real(mut cpx: ComplexAlgebraicRoot, real: RealAlgebraic) -> ComplexAlgebraic {
            match real {
                RealAlgebraic::Rational(rat) => {
                    cpx.tight_a += &rat;
                    cpx.tight_b += &rat;
                    //compose with x - n/d = dx - n
                    cpx.poly = Polynomial::compose(
                        &cpx.poly.apply_map(|c| Rational::from(c)),
                        &Polynomial::from_coeffs(vec![-rat, Rational::ONE]),
                    )
                    .primitive_part_fof();

                    let cpx = ComplexAlgebraic::Complex(cpx);

                    #[cfg(debug_assertions)]
                    cpx.check_invariants().unwrap();
                    cpx
                }
                RealAlgebraic::Real(mut real) => {
                    identify_complex_root(
                        root_sum_poly(&cpx.poly, &real.poly),
                        (0..).map(|i| {
                            if i != 0 {
                                cpx.refine();
                                real.refine();
                            }
                            let ans_tight_a = &cpx.tight_a + &real.tight_a;
                            let ans_tight_b = &cpx.tight_b + &real.tight_b;
                            let ans_tight_c = cpx.tight_c.clone();
                            let ans_tight_d = cpx.tight_d.clone();
                            (ans_tight_a, ans_tight_b, ans_tight_c, ans_tight_d)
                        }),
                    )

                    // loop {
                    //     let ans_tight_a = &cpx.tight_a + &real.tight_a;
                    //     let ans_tight_b = &cpx.tight_b + &real.tight_b;
                    //     let ans_tight_c = cpx.tight_c.clone();
                    //     let ans_tight_d = cpx.tight_d.clone();

                    //     match sum_poly.count_complex_roots(
                    //         &ans_tight_a,
                    //         &ans_tight_b,
                    //         &ans_tight_c,
                    //         &ans_tight_d,
                    //     ) {
                    //         Some(count) => {
                    //             if count == 0 {
                    //                 unreachable!();
                    //             } else if count == 1 {
                    //                 for (irr_poly, k) in
                    //                     sum_poly.factor().unwrap().unit_and_factors().1
                    //                 {
                    //                     debug_assert_eq!(k, Natural::ONE);
                    //                     if irr_poly
                    //                         .count_complex_roots(
                    //                             &ans_tight_a,
                    //                             &ans_tight_b,
                    //                             &ans_tight_c,
                    //                             &ans_tight_d,
                    //                         )
                    //                         .unwrap()
                    //                         == 1
                    //                     {
                    //                         let ans = ComplexAlgebraicRoot {
                    //                             tight_a: ans_tight_a,
                    //                             tight_b: ans_tight_b,
                    //                             tight_c: ans_tight_c,
                    //                             tight_d: ans_tight_d,
                    //                             poly: irr_poly,
                    //                         };
                    //                         #[cfg(debug_assertions)]
                    //                         ans.check_invariants().unwrap();
                    //                         return ans;
                    //                     }
                    //                 }
                    //             }
                    //         }
                    //         None => {}
                    //     }

                    //     cpx.refine();
                    //     real.refine();
                    // }
                }
            }
        }

        match (alg1, alg2) {
            (ComplexAlgebraic::Real(real1), ComplexAlgebraic::Real(real2)) => {
                ComplexAlgebraic::Real(RealAlgebraic::add(real1, real2))
            }
            (ComplexAlgebraic::Real(real1), ComplexAlgebraic::Complex(cpx2)) => {
                add_real(cpx2.clone(), real1.clone())
            }
            (ComplexAlgebraic::Complex(cpx1), ComplexAlgebraic::Real(real2)) => {
                add_real(cpx1.clone(), real2.clone())
            }
            (ComplexAlgebraic::Complex(cpx1), ComplexAlgebraic::Complex(cpx2)) => {
                let mut cpx1 = cpx1.clone();
                let mut cpx2 = cpx2.clone();

                identify_complex_root(
                    root_sum_poly(&cpx1.poly, &cpx2.poly),
                    (0..).map(|i| {
                        if i != 0 {
                            cpx1.refine();
                            cpx2.refine();
                        }
                        let ans_tight_a = &cpx1.tight_a + &cpx2.tight_a;
                        let ans_tight_b = &cpx1.tight_b + &cpx2.tight_b;
                        let ans_tight_c = &cpx1.tight_c + &cpx2.tight_c;
                        let ans_tight_d = &cpx1.tight_d + &cpx2.tight_d;
                        (ans_tight_a, ans_tight_b, ans_tight_c, ans_tight_d)
                    }),
                )
            }
        }
    }

    fn mul(&self, alg1: &Self::Set, alg2: &Self::Set) -> Self::Set {
        // println!("mul {:?} {:?}", alg1, alg2);
        // alg1.check_invariants().unwrap();
        // alg2.check_invariants().unwrap();

        fn mul_real(mut cpx: ComplexAlgebraicRoot, real: RealAlgebraic) -> ComplexAlgebraic {
            match real {
                RealAlgebraic::Rational(rat) => match rat.cmp(&Rational::ZERO) {
                    std::cmp::Ordering::Equal => ComplexAlgebraic::zero(),
                    std::cmp::Ordering::Less => mul_real(cpx.neg(), RealAlgebraic::Rational(-rat)),
                    std::cmp::Ordering::Greater => {
                        debug_assert!(rat > Rational::ZERO);
                        cpx.tight_a *= &rat;
                        cpx.tight_b *= &rat;
                        cpx.tight_c *= &rat;
                        cpx.tight_d *= &rat;
                        cpx.poly = root_pos_rat_mul_poly(cpx.poly, &rat);
                        #[cfg(debug_assertions)]
                        cpx.check_invariants().is_ok();
                        ComplexAlgebraic::Complex(cpx)
                    }
                },
                RealAlgebraic::Real(mut real) => {
                    identify_complex_root(
                        root_prod_poly(&cpx.poly, &real.poly),
                        (0..).map(|i| {
                            if i != 0 {
                                cpx.refine();
                                real.refine();
                            }
                            let mut pts_re = vec![];
                            let mut pts_im = vec![];
                            for (re, im) in [
                                (&cpx.tight_a, &cpx.tight_c),
                                (&cpx.tight_a, &cpx.tight_d),
                                (&cpx.tight_b, &cpx.tight_c),
                                (&cpx.tight_b, &cpx.tight_d),
                            ] {
                                for t in [&real.tight_a, &real.tight_b] {
                                    pts_re.push(t * re);
                                    pts_im.push(t * im);
                                }
                            }

                            let ans_tight_a = pts_re.iter().min().unwrap().clone();
                            let ans_tight_b = pts_re.into_iter().max().unwrap();
                            let ans_tight_c = pts_im.iter().min().unwrap().clone();
                            let ans_tight_d = pts_im.into_iter().max().unwrap();

                            //allow the bounds to expand a little bit so that the bounds are simpler
                            let diff_re = &ans_tight_b - &ans_tight_a;
                            let diff_im = &ans_tight_d - &ans_tight_c;
                            let ans_tight_a = Rational::simplest_rational_in_closed_interval(
                                &(&ans_tight_a - (&diff_re * Rational::from_str("1/10").unwrap())),
                                &ans_tight_a,
                            );
                            let ans_tight_b = Rational::simplest_rational_in_closed_interval(
                                &ans_tight_b,
                                &(&ans_tight_b + (&diff_re * Rational::from_str("1/10").unwrap())),
                            );
                            let ans_tight_c = Rational::simplest_rational_in_closed_interval(
                                &(&ans_tight_c - (&diff_im * Rational::from_str("1/10").unwrap())),
                                &ans_tight_c,
                            );
                            let ans_tight_d = Rational::simplest_rational_in_closed_interval(
                                &ans_tight_d,
                                &(&ans_tight_d + (&diff_im * Rational::from_str("1/10").unwrap())),
                            );
                            (ans_tight_a, ans_tight_b, ans_tight_c, ans_tight_d)
                        }),
                    )
                }
            }
        }

        match (alg1, alg2) {
            (ComplexAlgebraic::Real(real1), ComplexAlgebraic::Real(real2)) => {
                ComplexAlgebraic::Real(RealAlgebraic::mul(real1, real2))
            }
            (ComplexAlgebraic::Real(real1), ComplexAlgebraic::Complex(cpx2)) => {
                mul_real(cpx2.clone(), real1.clone())
            }
            (ComplexAlgebraic::Complex(cpx1), ComplexAlgebraic::Real(real2)) => {
                mul_real(cpx1.clone(), real2.clone())
            }
            (ComplexAlgebraic::Complex(cpx1), ComplexAlgebraic::Complex(cpx2)) => {
                let mut cpx1 = cpx1.clone();
                let mut cpx2 = cpx2.clone();

                identify_complex_root(
                    root_prod_poly(&cpx1.poly, &cpx2.poly),
                    (0..).map(|i| {
                        if i != 0 {
                            cpx1.refine();
                            cpx2.refine();
                        }

                        let mut pts_re = vec![];
                        let mut pts_im = vec![];
                        for (re1, im1) in [
                            (&cpx1.tight_a, &cpx1.tight_c),
                            (&cpx1.tight_a, &cpx1.tight_d),
                            (&cpx1.tight_b, &cpx1.tight_c),
                            (&cpx1.tight_b, &cpx1.tight_d),
                        ] {
                            for (re2, im2) in [
                                (&cpx2.tight_a, &cpx2.tight_c),
                                (&cpx2.tight_a, &cpx2.tight_d),
                                (&cpx2.tight_b, &cpx2.tight_c),
                                (&cpx2.tight_b, &cpx2.tight_d),
                            ] {
                                pts_re.push(re1 * re2 - im1 * im2);
                                pts_im.push(re1 * im2 + im1 * re2);
                            }
                        }

                        let ans_tight_a = pts_re.iter().min().unwrap().clone();
                        let ans_tight_b = pts_re.into_iter().max().unwrap();
                        let ans_tight_c = pts_im.iter().min().unwrap().clone();
                        let ans_tight_d = pts_im.into_iter().max().unwrap();

                        //allow the bounds to expand a little bit so that the bounds are simpler
                        let diff_re = &ans_tight_b - &ans_tight_a;
                        let diff_im = &ans_tight_d - &ans_tight_c;
                        let ans_tight_a = Rational::simplest_rational_in_closed_interval(
                            &(&ans_tight_a - (&diff_re * Rational::from_str("1/10").unwrap())),
                            &ans_tight_a,
                        );
                        let ans_tight_b = Rational::simplest_rational_in_closed_interval(
                            &ans_tight_b,
                            &(&ans_tight_b + (&diff_re * Rational::from_str("1/10").unwrap())),
                        );
                        let ans_tight_c = Rational::simplest_rational_in_closed_interval(
                            &(&ans_tight_c - (&diff_im * Rational::from_str("1/10").unwrap())),
                            &ans_tight_c,
                        );
                        let ans_tight_d = Rational::simplest_rational_in_closed_interval(
                            &ans_tight_d,
                            &(&ans_tight_d + (&diff_im * Rational::from_str("1/10").unwrap())),
                        );

                        (ans_tight_a, ans_tight_b, ans_tight_c, ans_tight_d)
                    }),
                )

                // loop {

                //     // println!("cpx1 = {:?}", cpx1);
                //     // println!("cpx2 = {:?}", cpx2);
                //     // println!(
                //     //     "abcd = {} {} {} {}",
                //     //     ans_tight_a, ans_tight_b, ans_tight_c, ans_tight_d
                //     // );

                //     match prod_poly.count_complex_roots(
                //         &ans_tight_a,
                //         &ans_tight_b,
                //         &ans_tight_c,
                //         &ans_tight_d,
                //     ) {
                //         Some(count) => {
                //             if count == 0 {
                //                 unreachable!();
                //             } else if count == 1 {
                //                 for (irr_poly, k) in
                //                     prod_poly.factor().unwrap().unit_and_factors().1
                //                 {
                //                     debug_assert_eq!(k, Natural::ONE);
                //                     if irr_poly
                //                         .count_complex_roots(
                //                             &ans_tight_a,
                //                             &ans_tight_b,
                //                             &ans_tight_c,
                //                             &ans_tight_d,
                //                         )
                //                         .unwrap()
                //                         == 1
                //                     {
                //                         let ans;
                //                         if irr_poly.degree().unwrap() == 1 {
                //                             ans = ComplexAlgebraic::Real(RealAlgebraic::Rational(
                //                                 unique_linear_root(&irr_poly),
                //                             ));
                //                         } else {
                //                             ans = ComplexAlgebraic::Complex(ComplexAlgebraicRoot {
                //                                 tight_a: ans_tight_a,
                //                                 tight_b: ans_tight_b,
                //                                 tight_c: ans_tight_c,
                //                                 tight_d: ans_tight_d,
                //                                 poly: irr_poly,
                //                             });
                //                         }
                //                         #[cfg(debug_assertions)]
                //                         ans.check_invariants().unwrap();
                //                         return ans;
                //                     }
                //                 }
                //             }
                //         }
                //         None => {}
                //     }

                //     cpx1.refine();
                //     cpx2.refine();
                // }
            }
        }
    }
}

impl IntegralDomainStructure for CannonicalStructure<ComplexAlgebraic> {
    fn inv(&self, a: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        // println!("inv {:?}", a);
        // a.check_invariants().unwrap();

        match a {
            ComplexAlgebraic::Real(a) => Ok(ComplexAlgebraic::Real(a.inv()?)),
            ComplexAlgebraic::Complex(a) => {
                let mut root = a.clone();

                let inv_poly = Polynomial::from_coeffs(
                    root.poly.clone().into_coeffs().into_iter().rev().collect(),
                )
                .fav_assoc();
                debug_assert!(inv_poly.is_irreducible());

                // println!("root = {} = {:?}", root, root);
                // println!("inv_poly = {}", inv_poly);
                // for root in inv_poly.all_complex_roots() {
                //     println!("root = {}", root);
                // }

                // z = the root represented by self
                // a = the center of the approximation
                // eps = the radius of the approximation

                //result used to get region which contains the reciprocal
                // if for some 0 < lambda < 1 we have |z - a| < eps <= (1 - lambda) |a| then |1/z - 1/a| < eps / lambda |a|^2

                loop {
                    //estimate eps such that |z - a| < eps
                    let eps = ((&root.tight_b - &root.tight_a) + (&root.tight_d - &root.tight_c))
                        / Rational::TWO;

                    let w_re = (&root.tight_a + &root.tight_b) / Rational::TWO;
                    let w_im = (&root.tight_c + &root.tight_d) / Rational::TWO;
                    let w_mag_sq = &w_re * &w_re + &w_im * &w_im;

                    // println!(
                    //     "w = {} + {} i",
                    //     rat_to_string(w_re.clone()),
                    //     rat_to_string(w_im.clone())
                    // );

                    //refine until eps < |a|
                    if &eps * &eps > w_mag_sq {
                        root.refine();
                        continue;
                    }

                    // find 0 < lambda < 1 such that eps < (1 - lambda) * |a|
                    let mut lambda = Rational::ONE_HALF;
                    while &eps * &eps >= (Rational::ONE - &lambda) * &w_mag_sq {
                        lambda = lambda * &Rational::ONE_HALF;
                    }
                    let lambda = lambda;
                    // println!("lambda = {}", lambda);

                    // |1/z - 1/a| < eps / lambda |a|^2

                    let w_recip_re = w_re / &w_mag_sq;
                    let w_recip_im = -w_im / &w_mag_sq;
                    let delta = &eps / (lambda * w_mag_sq);

                    let inv_tight_a = &w_recip_re - &delta;
                    let inv_tight_b = &w_recip_re + &delta;
                    let inv_tight_c = &w_recip_im - &delta;
                    let inv_tight_d = &w_recip_im + &delta;

                    // println!(
                    //     "w_inv = {} + {} i",
                    //     rat_to_string(w_recip_re),
                    //     rat_to_string(w_recip_im)
                    // );
                    // println!(
                    //     "eps = {}  delta = {}",
                    //     rat_to_string(eps),
                    //     rat_to_string(delta)
                    // );

                    match inv_poly.count_complex_roots(
                        &inv_tight_a,
                        &inv_tight_b,
                        &inv_tight_c,
                        &inv_tight_d,
                    ) {
                        Some(count) => {
                            if count == 0 {
                                unreachable!();
                            } else if count == 1 {
                                let ans = ComplexAlgebraic::Complex(ComplexAlgebraicRoot {
                                    tight_a: inv_tight_a,
                                    tight_b: inv_tight_b,
                                    tight_c: inv_tight_c,
                                    tight_d: inv_tight_d,
                                    poly: inv_poly,
                                });
                                #[cfg(debug_assertions)]
                                ans.check_invariants().unwrap();
                                return Ok(ans);
                            }
                        }
                        None => {}
                    }

                    root.refine();
                }
            }
        }
    }

    fn div(&self, a: &Self::Set, b: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        Ok(self.mul(a, &self.inv(b)?))
    }
}

impl FieldStructure for CannonicalStructure<ComplexAlgebraic> {}

impl CharZeroStructure for CannonicalStructure<ComplexAlgebraic> {}

impl ComplexSubsetStructure for CannonicalStructure<ComplexAlgebraic> {}

impl ComplexConjugateStructure for CannonicalStructure<ComplexAlgebraic> {
    fn conjugate(&self, x: &Self::Set) -> Self::Set {
        match x {
            ComplexAlgebraic::Real(x) => ComplexAlgebraic::Real(x.clone()),
            ComplexAlgebraic::Complex(x) => ComplexAlgebraic::Complex(x.clone().conj()),
        }
    }
}

impl PositiveRealNthRootStructure for CannonicalStructure<ComplexAlgebraic> {
    fn nth_root(&self, x: &Self::Set, n: usize) -> Result<Self::Set, ()> {
        match x {
            ComplexAlgebraic::Real(x) => Ok(ComplexAlgebraic::Real(x.nth_root(n)?)),
            ComplexAlgebraic::Complex(x) => Err(()),
        }
    }
}

impl ComplexAlgebraic {
    pub fn min_poly(&self) -> Polynomial<Rational> {
        match self {
            ComplexAlgebraic::Real(real_algebraic) => real_algebraic.min_poly(),
            ComplexAlgebraic::Complex(complex_root) => complex_root.min_poly(),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::rings::polynomial::polynomial::Polynomial;

    use super::*;

    #[test]
    fn test_root_sum_poly() {
        for (f, g, exp) in vec![
            (
                Polynomial::from_coeffs(vec![Integer::from(0)]),
                Polynomial::from_coeffs(vec![Integer::from(0)]),
                Polynomial::from_coeffs(vec![Integer::from(0)]),
            ),
            (
                Polynomial::from_coeffs(vec![Integer::from(0)]),
                Polynomial::from_coeffs(vec![Integer::from(1)]),
                Polynomial::from_coeffs(vec![Integer::from(1)]),
            ),
            (
                Polynomial::from_coeffs(vec![Integer::from(1)]),
                Polynomial::from_coeffs(vec![Integer::from(1)]),
                Polynomial::from_coeffs(vec![Integer::from(1)]),
            ),
            (
                Polynomial::from_coeffs(vec![Integer::from(-3), Integer::from(1)]),
                Polynomial::from_coeffs(vec![Integer::from(-5), Integer::from(1)]),
                Polynomial::from_coeffs(vec![Integer::from(-8), Integer::from(1)]),
            ),
            (
                Polynomial::from_coeffs(vec![Integer::from(1)]),
                Polynomial::from_coeffs(vec![Integer::from(-7), Integer::from(1)]),
                Polynomial::from_coeffs(vec![Integer::from(1)]),
            ),
            (
                Polynomial::from_coeffs(vec![Integer::from(-1), Integer::from(2)]),
                Polynomial::from_coeffs(vec![Integer::from(-1), Integer::from(3)]),
                Polynomial::from_coeffs(vec![Integer::from(-5), Integer::from(6)]),
            ),
            (
                Polynomial::from_coeffs(vec![
                    Integer::from(-1),
                    Integer::from(-2),
                    Integer::from(1),
                ]),
                Polynomial::from_coeffs(vec![
                    Integer::from(-2),
                    Integer::from(0),
                    Integer::from(1),
                ]),
                Polynomial::from_coeffs(vec![
                    Integer::from(-7),
                    Integer::from(5),
                    Integer::from(3),
                    Integer::from(-1),
                ]),
            ),
        ] {
            println!();
            let rsp = root_sum_poly(&f, &g);
            println!("f = {}", Polynomial::to_string(&f));
            println!("g = {}", Polynomial::to_string(&g));
            println!(
                "exp = {}    exp_factored = {:?}",
                Polynomial::to_string(&exp),
                exp.factorize_by_kroneckers_method()
            );
            println!(
                "rsp = {}    rsp_factored = {:?}",
                Polynomial::to_string(&rsp),
                rsp.factorize_by_kroneckers_method()
            );
            assert!(Polynomial::are_associate(&exp, &rsp));
        }
    }

    #[test]
    fn test_root_prod_poly() {
        for (f, g, exp) in vec![
            (
                Polynomial::from_coeffs(vec![Integer::from(0)]),
                Polynomial::from_coeffs(vec![Integer::from(0)]),
                Polynomial::from_coeffs(vec![Integer::from(0)]),
            ),
            (
                Polynomial::from_coeffs(vec![Integer::from(0)]),
                Polynomial::from_coeffs(vec![Integer::from(1)]),
                Polynomial::from_coeffs(vec![Integer::from(1)]),
            ),
            (
                Polynomial::from_coeffs(vec![Integer::from(1)]),
                Polynomial::from_coeffs(vec![Integer::from(1)]),
                Polynomial::from_coeffs(vec![Integer::from(1)]),
            ),
            (
                Polynomial::from_coeffs(vec![Integer::from(-3), Integer::from(1)]),
                Polynomial::from_coeffs(vec![Integer::from(-5), Integer::from(1)]),
                Polynomial::from_coeffs(vec![Integer::from(-15), Integer::from(1)]),
            ),
            (
                Polynomial::from_coeffs(vec![Integer::from(1)]),
                Polynomial::from_coeffs(vec![Integer::from(-7), Integer::from(1)]),
                Polynomial::from_coeffs(vec![Integer::from(1)]),
            ),
            (
                Polynomial::from_coeffs(vec![Integer::from(-1), Integer::from(2)]),
                Polynomial::from_coeffs(vec![Integer::from(-1), Integer::from(3)]),
                Polynomial::from_coeffs(vec![Integer::from(-1), Integer::from(6)]),
            ),
            (
                Polynomial::from_coeffs(vec![
                    Integer::from(-1),
                    Integer::from(-2),
                    Integer::from(1),
                ]),
                Polynomial::from_coeffs(vec![
                    Integer::from(-2),
                    Integer::from(0),
                    Integer::from(1),
                ]),
                Polynomial::from_coeffs(vec![
                    Integer::from(4),
                    Integer::from(0),
                    Integer::from(-12),
                    Integer::from(0),
                    Integer::from(1),
                ]),
            ),
            (
                Polynomial::from_coeffs(vec![
                    Integer::from(-2),
                    Integer::from(0),
                    Integer::from(1),
                ]),
                Polynomial::from_coeffs(vec![
                    Integer::from(-2),
                    Integer::from(0),
                    Integer::from(1),
                ]),
                Polynomial::from_coeffs(vec![
                    Integer::from(-4),
                    Integer::from(0),
                    Integer::from(1),
                ]),
            ),
        ] {
            println!();
            let rpp = root_prod_poly(&f, &g);
            println!("f = {}", Polynomial::to_string(&f));
            println!("g = {}", Polynomial::to_string(&g));
            println!(
                "exp = {}    exp_factored = {:?}",
                Polynomial::to_string(&exp),
                exp.factorize_by_kroneckers_method()
            );
            println!(
                "rpp = {}    rpp_factored = {:?}",
                Polynomial::to_string(&rpp),
                rpp.factorize_by_kroneckers_method()
            );
            assert!(Polynomial::are_associate(&exp, &rpp));
        }
    }

    #[test]
    fn test_squarefree_polynomial_real_root_isolation() {
        let f = Polynomial::product(vec![
            &Polynomial::from_coeffs(vec![
                Integer::from(-2),
                Integer::from(-4),
                Integer::from(-2),
            ]),
            &Polynomial::from_coeffs(vec![Integer::from(6), Integer::from(0), Integer::from(-3)]),
            &Polynomial::from_coeffs(vec![Integer::from(1), Integer::from(-3), Integer::from(1)]),
            &Polynomial::from_coeffs(vec![
                Integer::from(2),
                Integer::from(-3),
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
                Integer::from(1),
            ]),
            &Polynomial::from_coeffs(vec![
                Integer::from(1),
                Integer::from(-3),
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
                Integer::from(1),
            ]),
            &Polynomial::from_coeffs(vec![
                Integer::from(-1),
                Integer::from(12),
                Integer::from(-4),
                Integer::from(-15),
                Integer::from(5),
                Integer::from(3),
                Integer::from(-1),
            ]),
        ]);
        let f = f.primitive_squarefree_part();
        //f is a squarefree polynomial with lots of roots
        println!("f = {:?}", f);
        let intervals = Polynomial::real_roots_squarefree(f, None, None, false, false);
        println!("intervals = {:?}", &intervals);
        intervals.check_invariants().unwrap();

        let f =
            Polynomial::from_coeffs(vec![Integer::from(1), Integer::from(-3), Integer::from(1)]);
        println!("f = {:?}", f);
        let mut intervals = Polynomial::real_roots_squarefree(f, None, None, false, false);
        intervals.check_invariants().unwrap();
        intervals.clone().to_real_roots();
        for root in intervals.clone().to_real_roots() {
            println!("root = {:?}", root);
            root.check_invariants().unwrap();
        }
        println!("refine");
        for _i in 0..10 {
            intervals.refine_all();
        }
        for root in intervals.clone().to_real_roots() {
            println!("root = {:?}", root);
            root.check_invariants().unwrap();
        }

        let f = Polynomial::product(vec![
            &Polynomial::from_coeffs(vec![Integer::from(-1), Integer::from(1)]),
            &Polynomial::from_coeffs(vec![Integer::from(-2), Integer::from(1)]),
            &Polynomial::from_coeffs(vec![Integer::from(-3), Integer::from(1)]),
            &Polynomial::from_coeffs(vec![Integer::from(-4), Integer::from(1)]),
        ]);
        assert_eq!(
            f.clone()
                .real_roots_squarefree(
                    Some(&Rational::from(1)),
                    Some(&Rational::from(4)),
                    false,
                    false
                )
                .intervals
                .len(),
            2
        );
        assert_eq!(
            f.clone()
                .real_roots_squarefree(
                    Some(&Rational::from(1)),
                    Some(&Rational::from(4)),
                    true,
                    false
                )
                .intervals
                .len(),
            3
        );
        assert_eq!(
            f.clone()
                .real_roots_squarefree(
                    Some(&Rational::from(1)),
                    Some(&Rational::from(4)),
                    false,
                    true
                )
                .intervals
                .len(),
            3
        );
        assert_eq!(
            f.clone()
                .real_roots_squarefree(
                    Some(&Rational::from(1)),
                    Some(&Rational::from(4)),
                    true,
                    true
                )
                .intervals
                .len(),
            4
        );
    }

    #[test]
    fn test_real_roots_squarefree() {
        //a test where real_roots_squarefree find rational roots while refining the bounds on a real root such that no rational root lies on the end points
        let poly = Polynomial::from_coeffs(vec![
            Integer::from(0),
            Integer::from(9),
            Integer::from(0),
            Integer::from(-4),
        ]);
        let roots = poly.real_roots_squarefree(
            Some(&Rational::from(-2)),
            Some(&Rational::from(2)),
            false,
            true,
        );
        println!("{:?}", roots);
        roots.check_invariants().unwrap();
    }

    #[test]
    fn test_real_root_irreducible_count() {
        assert_eq!(
            Polynomial::from_coeffs(vec![
                Integer::from(3),
                Integer::from(-3),
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
                Integer::from(1)
            ])
            .real_roots_irreducible(None, None, false, false)
            .len(),
            1
        );
        assert_eq!(
            Polynomial::from_coeffs(vec![
                Integer::from(1),
                Integer::from(-3),
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
                Integer::from(1)
            ])
            .real_roots_irreducible(None, None, false, false)
            .len(),
            3
        );
    }

    #[test]
    fn test_real_algebraic_ordering() {
        let mut all_roots = vec![];
        for f in vec![
            Polynomial::from_coeffs(vec![
                Integer::from(-2),
                Integer::from(-4),
                Integer::from(-2),
            ]),
            Polynomial::from_coeffs(vec![Integer::from(6), Integer::from(0), Integer::from(-3)]),
            Polynomial::from_coeffs(vec![Integer::from(1), Integer::from(-3), Integer::from(1)]),
            Polynomial::from_coeffs(vec![
                Integer::from(2),
                Integer::from(-3),
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
                Integer::from(1),
            ]),
            Polynomial::from_coeffs(vec![
                Integer::from(1),
                Integer::from(-3),
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
                Integer::from(1),
            ]),
            Polynomial::from_coeffs(vec![
                Integer::from(-1),
                Integer::from(12),
                Integer::from(-4),
                Integer::from(-15),
                Integer::from(5),
                Integer::from(3),
                Integer::from(-1),
            ]),
        ] {
            for root in Polynomial::real_roots(&f, None, None, false, false) {
                all_roots.push(root.clone());
            }
        }

        all_roots.sort();

        for mut root in &mut all_roots {
            root.check_invariants().unwrap();
            match &mut root {
                RealAlgebraic::Rational(_a) => {}
                RealAlgebraic::Real(a) => {
                    a.refine_to_accuracy(&Rational::from_signeds(1, i64::MAX))
                }
            }
            println!("    {} {:?}", root.to_string(), root);
        }

        let mut all_roots_sorted_by_lower_tight_bound = all_roots.clone();
        all_roots_sorted_by_lower_tight_bound.sort_by_key(|root| match root {
            RealAlgebraic::Rational(a) => a.clone(),
            RealAlgebraic::Real(r) => r.tight_a.clone(),
        });
        assert_eq!(all_roots, all_roots_sorted_by_lower_tight_bound);
    }

    #[test]
    fn test_at_fixed_re_and_im() {
        let f = Polynomial::from_coeffs(vec![
            Integer::from(-1),
            Integer::from(3),
            Integer::from(0),
            Integer::from(1),
        ]);

        println!("f = {}", Polynomial::to_string(&f));

        let (vert_re_f, vert_im_f) = Polynomial::at_fixed_re(&f, &Rational::TWO);
        println!("re = {}", Polynomial::to_string(&vert_re_f));
        println!("im = {}", Polynomial::to_string(&vert_im_f));
        // f(z) = z^3 + 3z - 1
        // f(2 + xi) = (2 + xi)^3 + 3(2 + xi) - 1
        //           = 8 + 12xi - 6x^2 - x^3i + 6 + 3xi - 1
        //           = 13 + 15ix - 6x^2 - ix^3
        debug_assert_eq!(
            vert_re_f,
            Polynomial::from_coeffs(vec![Integer::from(13), Integer::from(0), Integer::from(-6)])
        );
        debug_assert_eq!(
            vert_im_f,
            Polynomial::from_coeffs(vec![
                Integer::from(0),
                Integer::from(15),
                Integer::from(0),
                Integer::from(-1)
            ])
        );

        let (vert_re_f, vert_im_f) = Polynomial::at_fixed_re(&f, &Rational::from_signeds(1, 2));
        println!("re = {}", Polynomial::to_string(&vert_re_f));
        println!("im = {}", Polynomial::to_string(&vert_im_f));
        // f(z) = z^3 + 3z - 1
        // f(1/2 + xi) = 5 + 30ix - 12x^2 - 8ix^3
        debug_assert_eq!(
            vert_re_f,
            Polynomial::from_coeffs(vec![Integer::from(5), Integer::from(0), Integer::from(-12)])
        );
        debug_assert_eq!(
            vert_im_f,
            Polynomial::from_coeffs(vec![
                Integer::from(0),
                Integer::from(30),
                Integer::from(0),
                Integer::from(-8)
            ])
        );

        let (vert_re_f, vert_im_f) = Polynomial::at_fixed_im(&f, &Rational::TWO);
        println!("re = {}", Polynomial::to_string(&vert_re_f));
        println!("im = {}", Polynomial::to_string(&vert_im_f));
        // f(z) = z^3 + 3z - 1
        // f(x + 2i) = -1 -2i -9x + 6ix^2 + x^3
        debug_assert_eq!(
            &vert_re_f,
            &Polynomial::from_coeffs(vec![
                Integer::from(-1),
                Integer::from(-9),
                Integer::from(0),
                Integer::from(1)
            ])
        );
        debug_assert_eq!(
            &vert_im_f,
            &Polynomial::from_coeffs(vec![Integer::from(-2), Integer::from(0), Integer::from(6),])
        );

        let (vert_re_f, vert_im_f) = Polynomial::at_fixed_im(&f, &Rational::from_signeds(1, 2));
        println!("re = {}", Polynomial::to_string(&vert_re_f));
        println!("im = {}", Polynomial::to_string(&vert_im_f));
        // f(z) = z^3 + 3z - 1
        // f(x + 1/2i) = -8 +11i + 18x + 12ix^2 + 8x^3
        debug_assert_eq!(
            vert_re_f,
            Polynomial::from_coeffs(vec![
                Integer::from(-8),
                Integer::from(18),
                Integer::from(0),
                Integer::from(8)
            ])
        );
        debug_assert_eq!(
            vert_im_f,
            Polynomial::from_coeffs(vec![Integer::from(11), Integer::from(0), Integer::from(12)])
        );
    }

    #[test]
    fn test_count_complex_roots() {
        //cyclotomic polynomials in a box of sidelength 4
        for k in 1..19 {
            let f = Polynomial::add(
                &Polynomial::var_pow(k),
                &Polynomial::neg(&Polynomial::one()),
            );
            let n = f
                .count_complex_roots(
                    &Rational::from(-2),
                    &Rational::from(2),
                    &Rational::from(-2),
                    &Rational::from(2),
                )
                .unwrap();
            assert_eq!(n, k);
        }

        //cyclotomic polynomials in a box with a boundary root iff k=0 mod 2
        for k in 1..19 {
            let f = Polynomial::add(
                &Polynomial::var_pow(k),
                &Polynomial::neg(&Polynomial::one()),
            );
            let n = Polynomial::count_complex_roots(
                &f,
                &Rational::from(-1),
                &Rational::from(3),
                &Rational::from(-3),
                &Rational::from(3),
            );
            if k % 2 == 0 {
                assert!(n.is_none());
            } else {
                assert_eq!(n.unwrap(), k);
            }
        }

        //cyclotomic polynomials in a box with a boundary root iff k=0 mod 4
        for k in 1..19 {
            let f = Polynomial::add(
                &Polynomial::var_pow(k),
                &Polynomial::neg(&Polynomial::one()),
            );
            let n = Polynomial::count_complex_roots(
                &f,
                &Rational::from(-2),
                &Rational::from(2),
                &Rational::from(-1),
                &Rational::from(1),
            );
            if k % 4 == 0 {
                assert!(n.is_none());
            } else {
                assert_eq!(n.unwrap(), k);
            }
        }

        //other test cases
        assert_eq!(
            Some(1),
            Polynomial::count_complex_roots(
                &Polynomial::from_coeffs(vec![
                    Integer::from(2),
                    Integer::from(-8),
                    Integer::from(1),
                    Integer::from(-4),
                    Integer::from(0),
                    Integer::from(1),
                ]),
                &Rational::from(-1),
                &Rational::from(1),
                &Rational::from(1),
                &Rational::from(2),
            )
        );

        assert_eq!(
            Some(3),
            Polynomial::count_complex_roots(
                &Polynomial::from_coeffs(vec![
                    Integer::from(2),
                    Integer::from(-8),
                    Integer::from(1),
                    Integer::from(-4),
                    Integer::from(0),
                    Integer::from(1),
                ]),
                &Rational::from(-3),
                &Rational::from(3),
                &Rational::from(-1),
                &Rational::from(1),
            )
        );

        //polynomial with roots 2+3i, 2-3i and counting box with 2+3i as a vertex
        assert_eq!(
            None,
            Polynomial::count_complex_roots(
                &Polynomial::from_coeffs(vec![
                    Integer::from(13),
                    Integer::from(-4),
                    Integer::from(1),
                ]),
                &Rational::from(2),
                &Rational::from(3),
                &Rational::from(3),
                &Rational::from(4),
            )
        );

        //x^2-x+1
        let f =
            Polynomial::from_coeffs(vec![Integer::from(1), Integer::from(-1), Integer::from(1)]);
        let n = f
            .count_complex_roots(
                &Rational::from(-1),
                &Rational::from(1),
                &Rational::from(-1),
                &Rational::from(1),
            )
            .unwrap();
        assert_eq!(n, 2);
    }

    #[test]
    fn test_count_complex_roots_big_values() {
        {
            let a = Rational::from_str("667/19382").unwrap(); //0.034413373232896505
            let b = Rational::from_str("754/4405").unwrap(); //0.17116912599318956
            let c = Rational::from_str("899/9691").unwrap(); //0.09276648436693839
            let d = Rational::from_str("3161/19382").unwrap(); //0.16308946445155298

            //0.951343405 -6.707838852x + 27.141574009x^2 = 0
            // x = 0.123571 - 0.140646 i
            // x = 0.123571 + 0.140646 i

            let poly = Polynomial::from_coeffs(vec![
                Integer::from_str("951343405").unwrap(),
                Integer::from_str("-6707838852").unwrap(),
                Integer::from_str("27141574009").unwrap(),
            ]);
            assert_eq!(poly.count_complex_roots(&a, &b, &c, &d), Some(1));
        }

        // {
        //     let a = Rational::from_str("-163099/329494").unwrap(); //-0.49499839147298585
        //     let b = Rational::from_str("-26827/74885").unwrap(); //-0.3582426387126928
        //     let c = Rational::from_str("899/9691").unwrap(); //0.09276648436693839
        //     let d = Rational::from_str("3161/19382").unwrap(); //0.16308946445155298

        //     //2139048288466 + 8191288386270x + 7843914888601x^2 = 0
        //     //same as
        //     //2.139048288466 + 8.191288386270x + 7.843914888601x^2 = 0
        //     //roots are -0.52 +/- a tiny complex bit

        //     let poly = Polynomial::from_coeffs(vec![
        //         Integer::from_str("2139048288466").unwrap(),
        //         Integer::from_str("8191288386270").unwrap(),
        //         Integer::from_str("7843914888601").unwrap(),
        //     ]);
        //     assert_eq!(poly.count_complex_roots(&a, &b, &c, &d), Some(1));
        // }
    }

    #[test]
    fn test_real_neg() {
        {
            let f = Polynomial::from_coeffs(vec![
                Integer::from(-2),
                Integer::from(0),
                Integer::from(1),
            ]);
            let roots = Polynomial::all_real_roots(&f);

            assert_eq!(roots.len(), 2);
            let a = &roots[0];
            let b = &roots[1];

            let a_neg = RealAlgebraic::neg(a);
            let b_neg = RealAlgebraic::neg(b);

            a_neg.check_invariants().unwrap();
            b_neg.check_invariants().unwrap();

            println!("a = {}", a.to_string());
            println!("b = {}", b.to_string());
            println!("a_neg = {}", a_neg.to_string());
            println!("b_neg = {}", b_neg.to_string());

            assert_ne!(a, b);
            assert_eq!(a, &b_neg);
            assert_eq!(b, &a_neg);
        }
        {
            let f = Polynomial::from_coeffs(vec![
                Integer::from(-1),
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
                Integer::from(3),
                Integer::from(1),
            ]);
            let roots = Polynomial::all_real_roots(&f);

            assert_eq!(roots.len(), 3);
            for root in roots {
                RealAlgebraic::neg(&root).check_invariants().unwrap();
            }
        }
        {
            //example where f(g(x)) is not primitive even though f and g are
            let f = Polynomial::from_coeffs(vec![
                Integer::from(-4),
                Integer::from(-1),
                Integer::from(1),
            ]);
            let roots = Polynomial::all_real_roots(&f);
            for root in roots {
                let root2 = RealAlgebraic::add(
                    &root,
                    &RealAlgebraic::from_rat(&Rational::from_signeds(1, 2)).unwrap(),
                );
                root2.check_invariants().unwrap();
            }
        }
    }

    #[test]
    fn test_real_add() {
        let f =
            Polynomial::from_coeffs(vec![Integer::from(-2), Integer::from(0), Integer::from(3)]);
        let roots = Polynomial::all_real_roots(&f);
        let a = RealAlgebraic::sum(roots.iter().collect());
        assert_eq!(a, RealAlgebraic::zero());

        let f = Polynomial::from_coeffs(vec![
            Integer::from(-7),
            Integer::from(0),
            Integer::from(100),
        ]);
        let roots = Polynomial::all_real_roots(&f);
        let a = RealAlgebraic::sum(roots.iter().collect());
        assert_eq!(a, RealAlgebraic::zero());

        let f = Polynomial::from_coeffs(vec![
            Integer::from(-100),
            Integer::from(0),
            Integer::from(7),
        ]);
        let roots = Polynomial::all_real_roots(&f);
        let a = RealAlgebraic::sum(roots.iter().collect());
        assert_eq!(a, RealAlgebraic::zero());
    }

    #[test]
    fn test_real_mul() {
        let f = Polynomial::from_coeffs(vec![
            Integer::from(-100),
            Integer::from(0),
            Integer::from(7),
        ]);
        // (x-a)(x-b) = x^2 - 100/7
        // so ab=-100/7
        let roots = Polynomial::all_real_roots(&f);
        let a = RealAlgebraic::product(roots.iter().collect());
        assert_eq!(
            a,
            RealAlgebraic::from_rat(&Rational::from_signeds(-100, 7)).unwrap()
        );
    }

    #[test]
    fn test_real_nth_root() {
        let x = &Polynomial::<Integer>::var().into_ring();
        let f = ((4 * x.pow(5) - 12 * x.pow(3) + 8 * x + 1)
            * (x + 1)
            * (x)
            * (x - 1)
            * (x - 2)
            * (x - 3)
            * (x - 4)
            * (x - 5)
            * (x - 144)
            * (x.pow(2) - 3))
            .into_set();
        let n = 2;

        for root in f.all_real_roots() {
            println!();
            println!("root = {}", root);
            match root.nth_root(n) {
                Ok(nth_root) => {
                    println!("YES {}-root = {}", n, nth_root);
                    debug_assert!(RealAlgebraic::zero() <= root);
                    debug_assert!(RealAlgebraic::zero() <= nth_root);
                    debug_assert_eq!(nth_root.nat_pow(&Natural::from(n)), root);
                }
                Err(()) => {
                    println!("NO {}-root", n);
                    debug_assert!(RealAlgebraic::zero() > root);
                }
            }
        }
    }

    #[test]
    fn test_all_complex_roots() {
        let f = Polynomial::from_coeffs(vec![
            Integer::from(-1),
            Integer::from(-1),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(1),
        ]);
        let roots = Polynomial::all_complex_roots(&f);
        assert_eq!(roots.len(), Polynomial::degree(&f).unwrap());
        for root in &roots {
            root.check_invariants().unwrap();
        }
    }

    #[test]
    fn test_complex_equal() {
        let x = &Polynomial::<Integer>::var().into_ring();
        let f = (x.pow(5) - x + 1).into_set();
        assert!(f.is_irreducible());
        let mut count = 0;
        for a in f.all_complex_roots() {
            for b in f.all_complex_roots() {
                if a == b {
                    count += 1;
                }
            }
        }
        assert_eq!(count, f.degree().unwrap());
    }

    #[test]
    fn test_complex_root_sum() {
        let f = Polynomial::from_coeffs(vec![
            Integer::from(1),
            Integer::from(0),
            Integer::from(0),
            Integer::from(1),
        ]);
        let roots = f.all_complex_roots();
        let s = ComplexAlgebraic::sum(roots.iter().collect());
        println!("{:?}", s);
        assert_eq!(s, ComplexAlgebraic::zero());

        let f = Polynomial::from_coeffs(vec![
            Integer::from(7),
            Integer::from(-3),
            Integer::from(42),
            Integer::from(9),
        ]);
        println!("f = {}", f);
        let roots = f.all_complex_roots();
        for root in &roots {
            println!("root = {}", root);
        }
        let s = ComplexAlgebraic::sum(roots.iter().collect());
        println!("root sum = {:?}", s);
        assert_eq!(
            s,
            ComplexAlgebraic::Real(RealAlgebraic::Rational(Rational::from_integers(
                Integer::from(-14),
                Integer::from(3)
            )))
        );
    }

    #[test]
    fn test_complex_mul() {
        let f = Polynomial::from_coeffs(vec![
            Integer::from(1),
            Integer::from(0),
            Integer::from(0),
            Integer::from(1),
        ]);
        let roots = f.all_complex_roots();
        let s = ComplexAlgebraic::product(roots.iter().collect());
        println!("{:?}", s);
        assert_eq!(s, ComplexAlgebraic::one().neg());

        let f = Polynomial::from_coeffs(vec![
            Integer::from(7),
            Integer::from(-3),
            Integer::from(42),
            Integer::from(9),
        ]);
        println!("f = {}", f);
        let roots = f.all_complex_roots();
        for root in &roots {
            println!("root = {}", root);
        }
        let s = ComplexAlgebraic::product(roots.iter().collect());
        println!("root prod = {:?}", s);
        assert_eq!(
            s,
            ComplexAlgebraic::Real(RealAlgebraic::Rational(Rational::from_integers(
                Integer::from(-7),
                Integer::from(9)
            )))
        );
    }

    #[test]
    fn test_complex_inv() {
        let x = &Polynomial::<Integer>::var().into_ring();
        let f = (x.pow(4) - x + 1).into_set();

        for root in f.all_complex_roots() {
            assert_eq!(
                ComplexAlgebraic::mul(&root.inv().unwrap(), &root),
                ComplexAlgebraic::one()
            );
        }
    }
}
