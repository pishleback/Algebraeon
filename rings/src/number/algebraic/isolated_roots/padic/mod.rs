use std::rc::Rc;

use algebraeon_sets::structure::MetaType;
use malachite_base::num::arithmetic::traits::UnsignedAbs;
use malachite_base::num::{
    arithmetic::traits::{DivMod, Pow},
    basic::traits::{One, Zero},
};
use malachite_nz::{integer::Integer, natural::Natural};
use malachite_q::Rational;

use crate::{
    number::natural::primes::is_prime,
    polynomial::polynomial::*,
    ring_structure::{quotient::QuotientStructure, structure::*},
};

fn pos_int_to_nat(x: Integer) -> Natural {
    debug_assert!(x >= Integer::ZERO);
    x.unsigned_abs()
}

fn padic_valuation(p: &Natural, mut n: Integer) -> Option<usize> {
    debug_assert!(is_prime(p));
    let p = Integer::from(p);
    if n == Natural::ZERO {
        None
    } else {
        let mut k = 0;
        let mut r;
        loop {
            (n, r) = n.div_mod(&p);
            if r == Natural::ZERO {
                k += 1;
                continue;
            } else {
                break;
            }
        }
        Some(k)
    }
}

fn padic_digits(p: &Natural, mut n: Natural) -> Vec<Natural> {
    debug_assert!(is_prime(p));
    let mut digits = vec![];
    let mut r;
    while n != 0 {
        (n, r) = n.div_mod(p);
        digits.push(r);
    }
    digits
}

#[derive(Debug, Clone)]
struct PAdicIntegerAlgebraicRoot {
    p: Natural,                 // a prime number
    poly: Polynomial<Integer>,  // a primitive irreducible degree >= 2 polynomial
    dpoly: Polynomial<Integer>, // the derivative of poly
    dpoly_valuation: usize, // f'(a) where a is the approximate root OR equivelently the lifted root.
    approx_root: Integer,   // an approximation to the root represented by this struct modulo p^k
    k: usize,
    // Furthermore approx_root must satisfy unique lifting to p-adic integers. By hensels lemma, this is the case whenever
    // |f(a)|_p < |f'(a)|^2
    // where f = poly and a = approx root
    // The true root b will agree with the approximate root a in the first dpoly_valuation+1 digits since |a-b| < |f'(a)| = |f'(b)|
}

impl PAdicIntegerAlgebraicRoot {
    fn modulus(&self) -> Natural {
        self.p.clone().pow(self.k as u64)
    }

    fn check(&self) -> Result<(), &'static str> {
        let pk = self.modulus();
        if !is_prime(&self.p) {
            return Err("p not prime");
        }
        match self.poly.degree() {
            Some(d) => {
                if d <= 1 {
                    return Err("deg(poly) <= 1");
                }
            }
            None => {
                return Err("poly = 0");
            }
        }
        if !Polynomial::is_irreducible(&self.poly) {
            return Err("poly is not irreducible");
        }
        if self.poly != self.poly.clone().primitive_part().unwrap() {
            return Err("f is not primitive");
        }
        if self.dpoly != self.poly.clone().derivative() {
            return Err("dpoly is not the derivative of poly");
        }
        if self.approx_root >= pk {
            return Err("approx root >= p^k");
        }

        let poly_mod_pk = PolynomialStructure::new(
            QuotientStructure::new_ring(Integer::structure(), Integer::from(pk)).into(),
        );

        let vfa = padic_valuation(&self.p, poly_mod_pk.evaluate(&self.poly, &self.approx_root));
        let vdfa = padic_valuation(
            &self.p,
            poly_mod_pk.evaluate(&self.dpoly, &self.approx_root),
        );
        match (vfa, vdfa) {
            (None, None) => {
                return Err("f(a) = f'(a) = 0");
            }
            (None, Some(_)) => {}
            (Some(_), None) => {
                return Err("f(a) != 0 and f'(a) = 0");
            }
            (Some(poly_val), Some(dpoly_val)) => {
                if !(poly_val > 2 * dpoly_val) {
                    return Err("|f(a)|_p < |f'(a)|^2 does not hold");
                }
            }
        }
        if vdfa.unwrap() != self.dpoly_valuation {
            return Err("v(f'(a)) does not match true value");
        }
        Ok(())
    }

    fn refine(&mut self) {
        self.k += 1;
        // p^{k+1}
        let pk1 = self.modulus();
        // Z/p^{k+1}Z
        let mod_pk1 = Rc::new(QuotientStructure::new_ring(
            Integer::structure(),
            Integer::from(&pk1),
        ));
        // Z/p^{k+1}Z[x]
        let poly_mod_pk1 = PolynomialStructure::new(mod_pk1.clone());

        // Update approximate root by Newtons method:
        // a <- a - f(a)/f'(a)
        let mut fa = poly_mod_pk1.evaluate(&self.poly, &self.approx_root); // f(a)
        let mut dfa = poly_mod_pk1.evaluate(&self.dpoly, &self.approx_root); // f'(a)
        let m = Integer::from(self.p.clone().pow(self.dpoly_valuation as u64)); // m is a common divisor of the top and bottom of f(a)/f'(a)
        (fa, dfa) = (fa / &m, dfa / &m);
        // now dfa != 0 mod p so we can find an inverse modulo p^{k+1}
        let (g, _, inv_dfa) = Integer::xgcd(&Integer::from(pk1), &dfa);
        debug_assert_eq!(g, Integer::ONE);
        self.approx_root =
            mod_pk1.add(&self.approx_root, &mod_pk1.neg(&mod_pk1.mul(&fa, &inv_dfa)));
    }

    fn refine_to_valuation(&mut self, k: usize) {
        while self.correct_approx_valuation() < k {
            self.refine();
        }
    }

    fn correct_approx_valuation(&self) -> usize {
        let pk = self.p.clone().pow(self.k as u64);
        let mod_pk = Rc::new(QuotientStructure::new_ring(
            Integer::structure(),
            Integer::from(&pk),
        ));
        let poly_mod_pk = PolynomialStructure::new(mod_pk.clone());
        // a = a-f(a)/f'(a) each iteration
        // v(f(a)) increases by at least v(f(a))-2v(f'(a)) each refinement
        // so the first v(f(a))-v(f'(a)) digits are correct
        let fa = poly_mod_pk.evaluate(&self.poly, &self.approx_root);
        let vfa = padic_valuation(&self.p, fa).unwrap_or(self.k);
        vfa - self.dpoly_valuation
    }

    fn reduce_modulo_valuation(&mut self, k: usize) -> Natural {
        self.refine_to_valuation(k);
        let pk = Integer::from(self.p.clone().pow(k as u64));
        pos_int_to_nat(Integer::rem(&self.approx_root, &pk))
    }

    fn equal_mut(&mut self, other: &mut Self) -> bool {
        let p = &self.p;
        if p != &other.p {
            return false;
        }
        if self.poly != other.poly {
            return false;
        }
        let dpoly_valuation = self.dpoly_valuation;
        if dpoly_valuation != other.dpoly_valuation {
            return false;
        }
        self.reduce_modulo_valuation(dpoly_valuation + 1)
            == other.reduce_modulo_valuation(dpoly_valuation + 1)
    }

    fn rightshift(&mut self) -> Option<Self> {
        if self.reduce_modulo_valuation(1) != Integer::ZERO {
            None
        } else {
            let (mul, poly) = self
                .poly
                .apply_map_with_powers(|(power, coeff)| {
                    coeff * Integer::from(&self.p).nat_pow(&Natural::from(power))
                })
                .factor_primitive()
                .unwrap();
            let dpoly = self.poly.clone().derivative();
            let approx_root = &self.approx_root / Integer::from(&self.p);
            let dpoly_valuation = self.dpoly_valuation + 1 - padic_valuation(&self.p, mul).unwrap();
            Some(Self {
                p: self.p.clone(),
                poly,
                dpoly,
                dpoly_valuation,
                approx_root,
                k: self.k,
            })
        }
    }

    /// Divide by the largest power of p possible, so that self != 0 mod p, and return the power
    fn fully_rightshift(&mut self) -> (Self, usize) {
        let mut rshifted = self.clone();
        let mut k = 0;
        loop {
            match rshifted.rightshift() {
                Some(new_rshifted) => {
                    rshifted = new_rshifted;
                    k += 1;
                }
                None => {
                    return (rshifted, k);
                }
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct PAdicAlgebraicRoot {
    // Multiply int_root by p^k to get the p-adic root represented by this struct
    root: PAdicIntegerAlgebraicRoot, // should be non-zero modulo p
    k: isize,                        // how much to left shift by
}

impl From<PAdicIntegerAlgebraicRoot> for PAdicAlgebraicRoot {
    fn from(mut value: PAdicIntegerAlgebraicRoot) -> Self {
        let (root, k) = value.fully_rightshift();
        PAdicAlgebraicRoot {
            root,
            k: k as isize,
        }
    }
}

impl PAdicAlgebraicRoot {
    fn shift_by(&mut self, k: isize) {
        self.k += k;
    }
}

impl PAdicAlgebraicRoot {
    fn check(&mut self) -> Result<(), &'static str> {
        if self.root.reduce_modulo_valuation(1) == Integer::ZERO {
            return Err("self.root = 0 mod p");
        }
        Ok(())
    }

    fn reduce_modulo_valuation<'a>(&'a mut self, k: isize) -> PAdicDigits<'a> {
        if k < self.k {
            PAdicDigits {
                p: &self.root.p,
                value: Natural::ZERO,
                shift: 0,
            }
        } else {
            let num_digits = (k - self.k) as usize;
            let value = self.root.reduce_modulo_valuation(num_digits);
            PAdicDigits {
                p: &self.root.p,
                value,
                shift: self.k,
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct PAdicRational {
    p: Natural, // a prime number
    rat: Rational,
}
impl PAdicRational {
    fn reduce_modulo_valuation<'a>(&'a self, k: isize) -> PAdicDigits<'a> {
        if self.rat == Rational::ZERO {
            PAdicDigits {
                p: &self.p,
                value: Natural::ZERO,
                shift: 0,
            }
        } else {
            let (n, d) = (self.rat.numerator(), self.rat.denominator());
            let a_uint = padic_valuation(&self.p, n.clone()).unwrap();
            let b_uint = padic_valuation(&self.p, d.clone()).unwrap();
            let (a, b) = (a_uint as isize, b_uint as isize);
            let shift = a - b;
            if k < shift {
                PAdicDigits {
                    p: &self.p,
                    value: Natural::ZERO,
                    shift: 0,
                }
            } else {
                let num_digits = (k - shift) as usize;
                let pa = Integer::from(&self.p).nat_pow(&Natural::from(a_uint)); // p^a
                let pb = Integer::from(&self.p).nat_pow(&Natural::from(b_uint)); // p^b
                let pn = Integer::from(&self.p).nat_pow(&Natural::from(num_digits)); // p^{num_digits}
                let (n, d) = (n / &pa, d / &pb);
                let (g, _, d_inv) = Integer::xgcd(&pn, &d);
                debug_assert_eq!(g, Integer::ONE);
                let value = pos_int_to_nat(Integer::rem(&(n * d_inv), &pn));
                PAdicDigits {
                    p: &self.p,
                    value,
                    shift,
                }
            }
        }
    }
}

#[derive(Debug, Clone)]
struct PAdicDigits<'a> {
    p: &'a Natural,
    value: Natural,
    shift: isize,
}
impl<'a> PAdicDigits<'a> {
    fn digits(&mut self) -> (Vec<Natural>, isize) {
        (padic_digits(self.p, self.value.clone()), self.shift)
    }
}

/// Store an algebraic p-adic number
#[derive(Debug, Clone)]
pub enum PAdicAlgebraic {
    Rational(PAdicRational),
    Algebraic(PAdicAlgebraicRoot),
}
impl From<PAdicIntegerAlgebraicRoot> for PAdicAlgebraic {
    fn from(value: PAdicIntegerAlgebraicRoot) -> Self {
        PAdicAlgebraic::Algebraic(value.into())
    }
}
impl From<PAdicAlgebraicRoot> for PAdicAlgebraic {
    fn from(value: PAdicAlgebraicRoot) -> Self {
        PAdicAlgebraic::Algebraic(value)
    }
}
impl PAdicAlgebraic {
    fn reduce_modulo_valuation<'a>(&'a mut self, k: isize) -> PAdicDigits<'a> {
        match self {
            PAdicAlgebraic::Rational(rational) => rational.reduce_modulo_valuation(k),
            PAdicAlgebraic::Algebraic(root) => root.reduce_modulo_valuation(k),
        }
    }
}

impl Polynomial<Integer> {
    fn all_padic_roots_irreducible(&self, p: &Natural) -> Vec<PAdicAlgebraic> {
        debug_assert!(is_prime(p));
        debug_assert!(self.is_irreducible());

        let f = self.clone();
        let d = f.degree().unwrap();

        debug_assert!(d > 0);
        if d == 1 {
            // Rational root
            let a = f.coeff(1);
            let b = f.coeff(0);
            // f(x) = ax + b
            // root = -b/a
            vec![PAdicAlgebraic::Rational(PAdicRational {
                p: p.clone(),
                rat: -Rational::from_integers(b, a),
            })]
        } else {
            // Algebraic root

            // Apply f(x) -> f(x/p) until the leading coefficient is not divisible by p i.e. is a p-adic integer.
            // Call the resulting polynomial g(x)
            let shift = padic_valuation(p, f.leading_coeff().unwrap()).unwrap();
            let g = f
                .apply_map_with_powers(|(power, coeff)| {
                    coeff * Integer::from(p).nat_pow(&Natural::from(shift * (d - power)))
                })
                .primitive_part()
                .unwrap();
            let dg = g.clone().derivative();

            // Manually lift roots until they uniquely lift to the p-adic integers
            let mut k = 0;
            let mut pk = Natural::ONE;
            let mut ununique_liftable_roots = vec![Natural::ZERO];
            let mut unique_liftable_roots = vec![];
            while !ununique_liftable_roots.is_empty() {
                let mut lifted_roots = vec![];
                let mod_pk1 =
                    QuotientStructure::new_ring(Integer::structure(), Integer::from(&pk * p));
                let poly_mod_pk1 = PolynomialStructure::new(mod_pk1.into());
                for root in ununique_liftable_roots {
                    let mut offset = Natural::ZERO;
                    while offset < *p {
                        let possible_lifted_root = &root + &pk * &offset;
                        let g_eval =
                            poly_mod_pk1.evaluate(&g, &Integer::from(&possible_lifted_root));
                        if g_eval == Integer::ZERO {
                            let g_valuation = k + 1;
                            let dg_eval =
                                poly_mod_pk1.evaluate(&dg, &Integer::from(&possible_lifted_root));
                            let dg_valuation = padic_valuation(p, dg_eval).unwrap_or(k + 1);
                            if g_valuation > 2 * dg_valuation {
                                let mut padic_int_root = PAdicIntegerAlgebraicRoot {
                                    p: p.clone(),
                                    poly: g.clone(),
                                    dpoly: dg.clone(),
                                    dpoly_valuation: dg_valuation,
                                    approx_root: Integer::from(possible_lifted_root),
                                    k: k + 1,
                                };
                                #[cfg(debug_assertions)]
                                padic_int_root.check().unwrap();
                                if !unique_liftable_roots
                                    .iter_mut()
                                    .any(|unique_liftable_root| {
                                        PAdicIntegerAlgebraicRoot::equal_mut(
                                            unique_liftable_root,
                                            &mut padic_int_root,
                                        )
                                    })
                                {
                                    unique_liftable_roots.push(padic_int_root);
                                }
                            } else {
                                lifted_roots.push(possible_lifted_root);
                            }
                        }
                        offset += Natural::ONE;
                    }
                }
                ununique_liftable_roots = lifted_roots;
                k += 1;
                pk *= p;
            }

            unique_liftable_roots
                .into_iter()
                .map(|root| {
                    {
                        let mut root = PAdicAlgebraicRoot::from(root);
                        root.shift_by(-(shift as isize));
                        root
                    }
                    .into()
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
    use crate::ring_structure::elements::IntoErgonomic;

    use super::*;

    #[test]
    fn test_padic_valuation() {
        assert_eq!(
            padic_valuation(&Natural::from(2u32), Integer::from(0)),
            None
        );
        assert_eq!(
            padic_valuation(&Natural::from(2u32), Integer::from(5)),
            Some(0)
        );
        assert_eq!(
            padic_valuation(&Natural::from(2u32), Integer::from(-12)),
            Some(2)
        );
        assert_eq!(
            padic_valuation(&Natural::from(2u32), Integer::from(256)),
            Some(8)
        );

        assert_eq!(
            padic_valuation(&Natural::from(7u32), Integer::from(0)),
            None
        );
        assert_eq!(
            padic_valuation(&Natural::from(7u32), Integer::from(-98)),
            Some(2)
        );
        assert_eq!(
            padic_valuation(&Natural::from(7u32), Integer::from(42)),
            Some(1)
        );
    }

    #[test]
    fn test_rational_padic_reduction() {}

    #[test]
    fn test_valid_padic_root_and_refine() {
        // Set up a 7-adic square-root of 2
        let poly =
            Polynomial::from_coeffs(vec![Integer::from(-2), Integer::from(0), Integer::from(1)]);
        let mut root = PAdicIntegerAlgebraicRoot {
            p: Natural::from(7u32),
            poly: poly.clone(),
            dpoly: poly.derivative(),
            dpoly_valuation: 0,
            approx_root: Integer::from(3),
            k: 1,
        };
        // Test refining it
        root.check().unwrap();
        root.refine();
        root.check().unwrap();
        root.refine();
        root.check().unwrap();
        root.refine();
        root.check().unwrap();
        root.refine();
        root.check().unwrap();
        println!("{:?}", root);

        debug_assert_eq!(root.reduce_modulo_valuation(0), Integer::from(0));
        debug_assert_eq!(root.reduce_modulo_valuation(1), Integer::from(3));
        debug_assert_eq!(root.reduce_modulo_valuation(2), Integer::from(10));
        debug_assert_eq!(root.reduce_modulo_valuation(3), Integer::from(108));
        debug_assert_eq!(
            root.reduce_modulo_valuation(20),
            Integer::from(75182500718243698u64)
        );
    }

    #[test]
    fn test_padic_root_rightshift() {
        let poly =
            Polynomial::from_coeffs(vec![Integer::from(-98), Integer::from(0), Integer::from(1)]);
        let mut root = PAdicIntegerAlgebraicRoot {
            p: Natural::from(7u32),
            poly: poly.clone(),
            dpoly: poly.derivative(),
            dpoly_valuation: 1,
            approx_root: Integer::from(21),
            k: 2,
        };

        let x = root.reduce_modulo_valuation(10);
        println!("{:?}", padic_digits(&root.p, x));
        root.check().unwrap();
        assert_eq!(root.reduce_modulo_valuation(0), Integer::from(0));
        assert_eq!(root.reduce_modulo_valuation(1), Integer::from(0));
        assert_eq!(root.reduce_modulo_valuation(2), Integer::from(21));
        assert_eq!(root.reduce_modulo_valuation(3), Integer::from(70));
        assert_eq!(root.reduce_modulo_valuation(4), Integer::from(756));
        root.check().unwrap();
        assert_eq!(
            PAdicAlgebraic::from(root.clone())
                .reduce_modulo_valuation(10)
                .digits(),
            (
                vec![3, 1, 2, 6, 1, 2, 1, 2, 4]
                    .into_iter()
                    .map(|x| (x as u32).into())
                    .collect(),
                1
            )
        );

        let (mut rshift_root, k) = root.clone().fully_rightshift();
        assert_eq!(k, 1);

        let x = rshift_root.reduce_modulo_valuation(10);
        println!("{:?}", padic_digits(&rshift_root.p, x));
        rshift_root.check().unwrap();
        assert_eq!(rshift_root.reduce_modulo_valuation(0), Integer::from(0));
        assert_eq!(rshift_root.reduce_modulo_valuation(1), Integer::from(3));
        assert_eq!(rshift_root.reduce_modulo_valuation(2), Integer::from(10));
        assert_eq!(rshift_root.reduce_modulo_valuation(3), Integer::from(108));
        assert_eq!(rshift_root.reduce_modulo_valuation(4), Integer::from(2166));
        rshift_root.check().unwrap();

        assert!(rshift_root.rightshift().is_none());
    }

    #[test]
    fn test_padic_digits() {
        println!(
            "{:?}",
            PAdicAlgebraic::Rational(PAdicRational {
                p: Natural::from(2u32),
                rat: Rational::from_integers(Integer::from(-1), Integer::from(1)),
            })
            .reduce_modulo_valuation(100)
            .digits()
        );

        assert_eq!(
            PAdicAlgebraic::Rational(PAdicRational {
                p: Natural::from(2u32),
                rat: Rational::from_integers(Integer::from(9), Integer::from(1)),
            })
            .reduce_modulo_valuation(30)
            .digits(),
            (
                vec![
                    Natural::from(1u32),
                    Natural::from(0u32),
                    Natural::from(0u32),
                    Natural::from(1u32)
                ],
                0
            )
        );

        assert_eq!(
            PAdicAlgebraic::Rational(PAdicRational {
                p: Natural::from(2u32),
                rat: Rational::from_integers(Integer::from(10), Integer::from(1)),
            })
            .reduce_modulo_valuation(30)
            .digits(),
            (
                vec![
                    Natural::from(1u32),
                    Natural::from(0u32),
                    Natural::from(1u32)
                ],
                1
            )
        );

        assert_eq!(
            PAdicAlgebraic::Rational(PAdicRational {
                p: Natural::from(2u32),
                rat: Rational::from_integers(Integer::from(31), Integer::from(36)),
            })
            .reduce_modulo_valuation(10)
            .digits(),
            (
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
                -2
            )
        );
    }

    #[test]
    fn test_all_padic_roots_example1() {
        let x = Polynomial::<Integer>::var().into_ergonomic();
        let f = (16 * x.pow(2) - 17).into_verbose();
        println!("{:?}", f);
        let roots = f.all_padic_roots(&Natural::from(2u32));
        assert_eq!(roots.len(), 2);
        for mut root in roots {
            println!("{:?}", root.reduce_modulo_valuation(10).digits());
        }
    }

    #[test]
    fn test_all_padic_roots_example2() {
        let x = Polynomial::<Integer>::var().into_ergonomic();
        let f = (x.pow(6) - 2).into_verbose();
        println!("{:?}", f);
        assert_eq!(f.all_padic_roots(&Natural::from(2u32)).len(), 0);
        assert_eq!(f.all_padic_roots(&Natural::from(7u32)).len(), 0);
        assert_eq!(f.all_padic_roots(&Natural::from(727u32)).len(), 6);
    }

    #[test]
    fn test_all_padic_roots_example3() {
        let x = Polynomial::<Integer>::var().into_ergonomic();
        let f = (3 * x - 2).into_verbose();
        println!("{:?}", f);
        for mut root in f.all_padic_roots(&Natural::from(7u32)) {
            println!("{:?}", root.reduce_modulo_valuation(10).digits());
        }
        assert_eq!(f.all_padic_roots(&Natural::from(2u32)).len(), 1);
        assert_eq!(f.all_padic_roots(&Natural::from(3u32)).len(), 1);
        assert_eq!(f.all_padic_roots(&Natural::from(5u32)).len(), 1);
    }
}
