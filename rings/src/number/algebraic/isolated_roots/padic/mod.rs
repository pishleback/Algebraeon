use malachite_base::num::arithmetic::traits::DivMod;
use malachite_base::num::basic::traits::Zero;
use malachite_nz::{integer::Integer, natural::Natural};
use malachite_q::Rational;

use crate::{number::natural::primes::*, polynomial::polynomial::*, structure::structure::*};

// Some algorithms here on p-adic root isolation can be found in
// Sturm, Thomas & Weispfenning, Volker. (2004). P-adic Root Isolation. Revista de la Real Academia de Ciencias Exactas, Físicas y Naturales. Serie A, Matemáticas.
// https://www.researchgate.net/profile/Thomas-Sturm-2/publication/2925550_P-adic_Root_Isolation/links/580b8c0708aeef1bfeeb5db8/P-adic-Root-Isolation.pdf?origin=scientificContributions

fn padic_int_valuation(p: &Natural, mut n: Integer) -> Option<Natural> {
    debug_assert!(is_prime(p));
    let p = Integer::from(p);
    if n == Natural::ZERO {
        None
    } else {
        let mut k = 0usize;
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
        Some(Natural::from(k))
    }
}

mod balancable_pairs {
    use std::collections::HashSet;

    use super::*;

    /// A balancable pair of a polynomial consists of two of the monomial terms satisfying some conditions:
    /// For a fixed choice of prime $p$, a balancable pair of a polynomial $$f(x) = f_nx^n + ... + f_2x^2 + f_1x + f_0$$ is
    /// a pair of indicies $0 \le i < j \le n$ such that
    /// - $f_i \ne 0$
    /// - $f_j \ne 0$
    /// - $\frac{v_p(f_j) - v_p(f_i)}{j - i}$ is an integer, called the balancing value
    #[derive(Debug)]
    pub struct BalancablePair<'a> {
        p: Natural,
        f: &'a Polynomial<Integer>,
        n: usize, // n = deg(f)
        i: usize,
        j: usize,
        vfi: Natural, // v_p(f_i)
        vfj: Natural, // v_p(f_j)
        bv: Integer,  // the balancing value
    }

    impl<'a> BalancablePair<'a> {
        /// A balancable pair is critical if
        /// $$\frac{jv_p(f_i) - iv_p(f_j)}{j-i} = v_p(f_i) + i\frac{v_p(f_i) - v_p(f_j)}{j-i} = \min_{k = 1, \dots, n} v_p(f_k) + k \frac{v_p(f_i) - v_p(f_j)}{j - i}$$
        /// If this is the case then the balancing value is called a critical value for $f$.
        /// If $\alpha \in \mathbb{Q}_p$ is such that $f(\alpha) = 0$ then $v_p(\alpha)$ is a critical value for $f$.
        pub fn is_critical(&self) -> bool {
            let mut min = None;
            for k in 1..(self.n + 1) {
                match padic_int_valuation(&self.p, self.f.coeff(k)) {
                    Some(vfk) => {
                        let val = Integer::from(vfk) + Integer::from(k) * &self.bv;
                        if min.is_none() {
                            min = Some(val);
                        } else {
                            if &val < min.as_ref().unwrap() {
                                min = Some(val);
                            }
                        }
                    }
                    None => {
                        // vfk = inf
                    }
                }
            }
            min.unwrap() == Integer::from(&self.vfi) + Integer::from(self.i) * &self.bv
        }

        pub fn balancing_value(&self) -> &Integer {
            &self.bv
        }
    }

    impl Polynomial<Integer> {
        pub fn balancable_pairs<'a>(&'a self, p: &Natural) -> Vec<BalancablePair<'a>> {
            assert!(!self.is_zero());
            debug_assert!(is_prime(&p));
            let mut bps = vec![];
            let n = self.degree().unwrap();
            let coeff_valuations = (0..(n + 1))
                .map(|k| padic_int_valuation(p, self.coeff(k)))
                .collect::<Vec<_>>();
            for i in 0..(n + 1) {
                for j in (i + 1)..(n + 1) {
                    match (&coeff_valuations[i], &coeff_valuations[j]) {
                        (Some(vfi), Some(vfj)) => {
                            match Integer::div(
                                &(Integer::from(vfi) - Integer::from(vfj)),
                                &Integer::from(j - i),
                            ) {
                                Ok(bv) => bps.push(BalancablePair {
                                    p: p.clone(),
                                    f: &self,
                                    n,
                                    i,
                                    j,
                                    vfi: vfi.clone(),
                                    vfj: vfj.clone(),
                                    bv,
                                }),
                                Err(_) => {
                                    //\frac{v_p(f_j)-v_p(f_i)}{j-i} is not an inteer
                                }
                            }
                        }
                        _ => {
                            // f_i=0 or f_j=0
                        }
                    }
                }
            }
            bps
        }

        pub fn critical_values(&self, p: &Natural) -> HashSet<Integer> {
            let mut cvs = HashSet::new();
            for bp in self.balancable_pairs(p) {
                if bp.is_critical() {
                    cvs.insert(bp.balancing_value().clone());
                }
            }
            cvs
        }
    }

    #[cfg(test)]
    mod tests {
        use crate::structure::elements::*;

        use super::*;

        #[test]
        fn test_balancable_pair() {
            let x = Polynomial::<Integer>::var().into_ergonomic();
            let f = (81 * x.pow(4) - 6 * &x + 5).into_verbose();
            println!("f = {}", f);
            {
                let p = Natural::from(2u32);
                let mut bps = f.balancable_pairs(&p).into_iter();

                let bp = bps.next().unwrap();
                assert_eq!(bp.p, p);
                assert_eq!(bp.n, 4);
                assert_eq!(bp.i, 0);
                assert_eq!(bp.j, 1);
                assert_eq!(bp.vfi, Natural::from(0u32));
                assert_eq!(bp.vfj, Natural::from(1u32));
                assert_eq!(bp.bv, Integer::from(-1));
                assert!(!bp.is_critical());

                let bp = bps.next().unwrap();
                assert_eq!(bp.p, p);
                assert_eq!(bp.n, 4);
                assert_eq!(bp.i, 0);
                assert_eq!(bp.j, 4);
                assert_eq!(bp.vfi, Natural::from(0u32));
                assert_eq!(bp.vfj, Natural::from(0u32));
                assert_eq!(bp.bv, Integer::from(0));
                assert!(bp.is_critical());

                assert!(bps.next().is_none());

                debug_assert_eq!(f.critical_values(&p), HashSet::from([Integer::from(0)]));
            }
            {
                let p = Natural::from(3u32);
                let mut bps = f.balancable_pairs(&p).into_iter();

                let bp = bps.next().unwrap();
                assert_eq!(bp.p, p);
                assert_eq!(bp.n, 4);
                assert_eq!(bp.i, 0);
                assert_eq!(bp.j, 1);
                assert_eq!(bp.vfi, Natural::from(0u32));
                assert_eq!(bp.vfj, Natural::from(1u32));
                assert_eq!(bp.bv, Integer::from(-1));
                assert!(bp.is_critical());

                let bp = bps.next().unwrap();
                assert_eq!(bp.p, p);
                assert_eq!(bp.n, 4);
                assert_eq!(bp.i, 0);
                assert_eq!(bp.j, 4);
                assert_eq!(bp.vfi, Natural::from(0u32));
                assert_eq!(bp.vfj, Natural::from(4u32));
                assert_eq!(bp.bv, Integer::from(-1));
                assert!(bp.is_critical());

                let bp = bps.next().unwrap();
                assert_eq!(bp.p, p);
                assert_eq!(bp.n, 4);
                assert_eq!(bp.i, 1);
                assert_eq!(bp.j, 4);
                assert_eq!(bp.vfi, Natural::from(1u32));
                assert_eq!(bp.vfj, Natural::from(4u32));
                assert_eq!(bp.bv, Integer::from(-1));
                assert!(bp.is_critical());

                assert!(bps.next().is_none());

                debug_assert_eq!(f.critical_values(&p), HashSet::from([Integer::from(-1)]));
            }
            {
                let p = Natural::from(5u32);
                let mut bps = f.balancable_pairs(&p).into_iter();

                let bp = bps.next().unwrap();
                assert_eq!(bp.p, p);
                assert_eq!(bp.n, 4);
                assert_eq!(bp.i, 0);
                assert_eq!(bp.j, 1);
                assert_eq!(bp.vfi, Natural::from(1u32));
                assert_eq!(bp.vfj, Natural::from(0u32));
                assert_eq!(bp.bv, Integer::from(1));
                assert!(bp.is_critical());

                let bp = bps.next().unwrap();
                assert_eq!(bp.p, p);
                assert_eq!(bp.n, 4);
                assert_eq!(bp.i, 1);
                assert_eq!(bp.j, 4);
                assert_eq!(bp.vfi, Natural::from(0u32));
                assert_eq!(bp.vfj, Natural::from(0u32));
                assert_eq!(bp.bv, Integer::from(0));
                assert!(bp.is_critical());

                assert!(bps.next().is_none());

                debug_assert_eq!(
                    f.critical_values(&p),
                    HashSet::from([Integer::from(0), Integer::from(1)])
                );
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct PAdicBall {
    a: Rational,
    v: Integer,
}

impl Polynomial<Integer> {
    fn all_padic_roots_irreducible(&self, p: &Natural) -> Vec<PAdicBall> {
        todo!()
    }

    pub fn all_padic_roots(&self, p: &Natural) -> Vec<PAdicBall> {
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
    use super::*;

    #[test]
    fn test_padic_valuation() {
        assert_eq!(
            padic_int_valuation(&Natural::from(2u32), Integer::from(0)),
            None
        );
        assert_eq!(
            padic_int_valuation(&Natural::from(2u32), Integer::from(5)),
            Some(Natural::from(0u32))
        );
        assert_eq!(
            padic_int_valuation(&Natural::from(2u32), Integer::from(-12)),
            Some(Natural::from(2u32))
        );
        assert_eq!(
            padic_int_valuation(&Natural::from(2u32), Integer::from(256)),
            Some(Natural::from(8u32))
        );

        assert_eq!(
            padic_int_valuation(&Natural::from(7u32), Integer::from(0)),
            None
        );
        assert_eq!(
            padic_int_valuation(&Natural::from(7u32), Integer::from(-98)),
            Some(Natural::from(2u32))
        );
        assert_eq!(
            padic_int_valuation(&Natural::from(7u32), Integer::from(42)),
            Some(Natural::from(1u32))
        );
    }
}
