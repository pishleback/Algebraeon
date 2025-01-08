use std::rc::Rc;

use algebraeon_sets::structure::MetaType;
use malachite_base::num::{
    arithmetic::traits::{DivMod, Pow},
    basic::traits::{One, Zero},
};
use malachite_nz::{integer::Integer, natural::Natural};
use malachite_q::Rational;

use super::natural::primes::is_prime;
use crate::{
    polynomial::polynomial::*,
    ring_structure::{quotient::QuotientStructure, structure::*},
};

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

#[derive(Debug, Clone)]
struct PAdicIntegerAlgebraicRoot {
    p: Natural,                 // a prime number
    poly: Polynomial<Integer>,  // an irreducible degree >= 2 polynomial
    dpoly: Polynomial<Integer>, // the derivative of poly
    approx_root: Integer, // an approximation to the root represented by this struct modulo p^k
    k: usize,
    // Furthermore approx_root must satisfy unique lifting to p-adic integers. By hensels lemma, this is the case whenever
    // |f(a)|_p < |f'(a)|^2
    // where f = poly and a = approx root
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
        if self.dpoly != self.poly.clone().derivative() {
            return Err("dpoly is not the derivative of poly");
        }
        if self.approx_root >= pk {
            return Err("approx root >= p^k");
        }

        let poly_mod_pk = PolynomialStructure::new(
            QuotientStructure::new_ring(Integer::structure(), Integer::from(pk)).into(),
        );

        match (
            padic_valuation(&self.p, poly_mod_pk.evaluate(&self.poly, &self.approx_root)),
            padic_valuation(
                &self.p,
                poly_mod_pk.evaluate(&self.dpoly, &self.approx_root),
            ),
        ) {
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
        let (g, _, dpoly_at_approx_root_inv) = Integer::xgcd(
            &Integer::from(pk1),
            &poly_mod_pk1.evaluate(&self.dpoly, &self.approx_root),
        );
        debug_assert_eq!(g, Integer::ONE);
        self.approx_root = mod_pk1.add(
            &self.approx_root,
            &mod_pk1.neg(&mod_pk1.mul(
                &poly_mod_pk1.evaluate(&self.poly, &self.approx_root),
                &dpoly_at_approx_root_inv,
            )),
        );
    }
}

#[derive(Debug, Clone)]
pub struct PAdicAlgebraicRoot {
    // Multiply int_root by p^k to get the p-adic root represented by this struct
    root: PAdicIntegerAlgebraicRoot, // should be non-zero modulo p
    k: Integer,
}

#[derive(Debug, Clone)]
pub enum PAdicAlgebraic {
    Rational(Rational),
    Algebraic(PAdicAlgebraicRoot),
}

#[cfg(test)]
mod tests {
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
    fn test_valid_padic_root_and_refine() {
        let poly =
            Polynomial::from_coeffs(vec![Integer::from(-2), Integer::from(0), Integer::from(1)]);
        let mut root = PAdicIntegerAlgebraicRoot {
            p: Natural::from(7u32),
            poly: poly.clone(),
            dpoly: poly.derivative(),
            approx_root: Integer::from(3),
            k: 1,
        };
        assert!(root.check().is_ok());
        root.refine();
        assert!(root.check().is_ok());
        root.refine();
        assert!(root.check().is_ok());
        root.refine();
        assert!(root.check().is_ok());
        root.refine();
        assert!(root.check().is_ok());
        println!("{:?}", root);
    }
}
