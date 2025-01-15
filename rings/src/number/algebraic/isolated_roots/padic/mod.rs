use malachite_base::num::arithmetic::traits::{DivMod, UnsignedAbs};
use malachite_base::num::basic::traits::{One, Two, Zero};
use malachite_nz::{integer::Integer, natural::Natural};
use malachite_q::Rational;

use crate::{number::natural::primes::*, polynomial::polynomial::*, structure::structure::*};

// Some algorithms here on p-adic root isolation can be found in
// Sturm, Thomas & Weispfenning, Volker. (2004). P-adic Root Isolation. Revista de la Real Academia de Ciencias Exactas, Físicas y Naturales. Serie A, Matemáticas.
// https://www.researchgate.net/profile/Thomas-Sturm-2/publication/2925550_P-adic_Root_Isolation/links/580b8c0708aeef1bfeeb5db8/P-adic-Root-Isolation.pdf?origin=scientificContributions

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
                    v.unsigned_abs()
                }
            }
        }
    }
}

fn padic_int_valuation(p: &Natural, mut n: Integer) -> Valuation {
    debug_assert!(is_prime(p));
    let p = Integer::from(p);
    if n == Natural::ZERO {
        Valuation::Infinity
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
        Valuation::Finite(Integer::from(k))
    }
}

mod balancable_pairs {
    use std::collections::HashSet;

    use super::*;

    /// A balancable pair of a polynomial consists of two of the monomial terms satisfying some conditions.
    ///
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
        vfi: Integer, // v_p(f_i)
        vfj: Integer, // v_p(f_j)
        bv: Integer,  // the balancing value
    }

    impl<'a> BalancablePair<'a> {
        /// A balancable pair is critical if
        /// $$\frac{jv_p(f_i) - iv_p(f_j)}{j-i} = v_p(f_i) + i\frac{v_p(f_i) - v_p(f_j)}{j-i} = \min_{k = 1, \dots, n} v_p(f_k) + k \frac{v_p(f_i) - v_p(f_j)}{j - i}$$
        /// If this is the case then the balancing value is called a critical value for $f$.
        /// If $\alpha \in \mathbb{Q}_p$ is such that $f(\alpha) = 0$ then $v_p(\alpha)$ is a critical value for $f$.
        pub fn is_critical(&self) -> bool {
            let min = (1..(self.n + 1))
                .filter_map(|k| match padic_int_valuation(&self.p, self.f.coeff(k)) {
                    Valuation::Infinity => None,
                    Valuation::Finite(vfk) => {
                        Some(Integer::from(vfk) + Integer::from(k) * &self.bv)
                    }
                })
                .min()
                .unwrap();
            min == self.crossmul_balancing_value()
        }

        /// $\frac{v_p(f_j) - v_p(f_i)}{j - i}$
        pub fn balancing_value(&self) -> &Integer {
            &self.bv
        }

        /// $\frac{iv_p(f_j) - jv_p(f_i)}{j - i}$
        pub fn crossmul_balancing_value(&self) -> Integer {
            &self.vfi + Integer::from(self.i) * &self.bv
        }

        /// The polynomial
        /// $$f_{i,j}(x) = p^{-\frac{iv_p(f_j) - jv_p(f_i)}{j - i}} \sum_{k=0}^n f_k \left(p^{\frac{v_p(f_j) - v_p(f_i)}{j - i}} x\right)^k$$
        /// It is such that for $\alpha \in \mathbb{Q}_p$
        /// - $f_{i,j}(\alpha) = 0$ if and only if $f\left(p^{\frac{v_p(f_j) - v_p(f_i)}{j - i}} \alpha\right) = 0$
        /// - $f_{i,j}^{\prime}(\alpha) = 0$ if and only if $f^{\prime}\left(p^{\frac{v_p(f_j) - v_p(f_i)}{j - i}} \alpha\right) = 0$
        pub fn normalization(&self) -> Polynomial<Integer> {
            debug_assert!(self.is_critical());
            let cmbv = self.crossmul_balancing_value();
            Polynomial::from_coeffs(
                (0..(self.n + 1))
                    .map(|k| {
                        let p_pow = self.balancing_value() * Integer::from(k) - &cmbv;
                        // compute f_k*p^p_pow
                        if p_pow >= Integer::ZERO {
                            let p_pow = p_pow.unsigned_abs();
                            self.f.coeff(k) * Integer::from(self.p.nat_pow(&p_pow))
                        } else {
                            let neg_p_pow = (-p_pow).unsigned_abs();
                            Integer::div(
                                &self.f.coeff(k),
                                &Integer::from(self.p.nat_pow(&neg_p_pow)),
                            )
                            .unwrap()
                        }
                    })
                    .collect(),
            )
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
                        (Valuation::Finite(vfi), Valuation::Finite(vfj)) => {
                            match Integer::div(&(vfi - vfj), &Integer::from(j - i)) {
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
                debug_assert_eq!(
                    bp.normalization(),
                    (10125 * x.pow(4) - 6 * &x + 1).into_verbose()
                );

                let bp = bps.next().unwrap();
                assert_eq!(bp.p, p);
                assert_eq!(bp.n, 4);
                assert_eq!(bp.i, 1);
                assert_eq!(bp.j, 4);
                assert_eq!(bp.vfi, Natural::from(0u32));
                assert_eq!(bp.vfj, Natural::from(0u32));
                assert_eq!(bp.bv, Integer::from(0));
                assert!(bp.is_critical());
                debug_assert_eq!(
                    bp.normalization(),
                    (81 * x.pow(4) - 6 * &x + 5).into_verbose()
                );

                assert!(bps.next().is_none());

                debug_assert_eq!(
                    f.critical_values(&p),
                    HashSet::from([Integer::from(0), Integer::from(1)])
                );
            }
        }
    }
}

mod isolate {
    use super::*;

    /// Represent all p-adic integers which differ from `center` by a valuation of at most `valuation_radius`
    #[derive(Debug)]
    pub struct PAdicIntegerVal0Ball {
        center: Integer,
        valuation_radius: Natural,
    }

    /// A brute-force algorithm to isolate all roots of f of valuation 0
    /// Return ([r1, ...], alpha) such that each ri is alpha-close to exactly one valuation 0 root of f and vise-versa
    pub fn isolatebf0(p: &Natural, f: &Polynomial<Integer>) -> Vec<PAdicIntegerVal0Ball> {
        debug_assert!(!f.is_zero());
        debug_assert!(f.is_squarefree());
        debug_assert!(is_prime(p));
        let n = f.degree().unwrap();
        if n == 0 {
            vec![]
        } else {
            let mut roots = vec![];
            // disc(f) != 0 since f is squarefree
            let alpha = padic_int_valuation(p, f.clone().discriminant().unwrap()).unwrap_nat();
            let two_alpha = Natural::TWO * &alpha;
            let mut i = Natural::ONE;
            let max_i = p.nat_pow(&(&two_alpha + Natural::ONE));
            while i < max_i {
                if &i % p != 0 {
                    if padic_int_valuation(p, f.evaluate(&Integer::from(&i)))
                        > Valuation::Finite(Integer::from(&two_alpha))
                    {
                        let int_i = Integer::from(&i);
                        if !roots.iter().any(|alt_i| {
                            padic_int_valuation(p, &int_i - alt_i)
                                > Valuation::Finite(Integer::from(&alpha))
                        }) {
                            roots.push(int_i)
                        }
                    }
                }
                i += Natural::ONE;
            }
            roots
                .into_iter()
                .map(|i| PAdicIntegerVal0Ball {
                    center: i,
                    valuation_radius: alpha.clone(),
                })
                .collect()
        }
    }

    /// Return ([..., (ri, bi), ...], alpha) such that each ri is bi-close to exactly one valuation 0 root of f and vise-versa
    pub fn isolate0(p: &Natural, f: &Polynomial<Integer>) -> Vec<PAdicIntegerVal0Ball> {
        debug_assert!(!f.is_zero());
        debug_assert!(f.is_squarefree());
        debug_assert!(is_prime(p));
        let n = f.degree().unwrap();
        if n == 0 {
            vec![]
        } else {
            let mut roots = vec![];
            // disc(f) != 0 since f is squarefree
            let alpha = padic_int_valuation(p, f.clone().discriminant().unwrap()).unwrap_nat();
            let mut i = Natural::ONE;
            while &i < p {
                // i = 1, ..., p-1
                roots.append(&mut isorefine(p, f, &alpha, &i, &Natural::ZERO));
                i += Natural::ONE;
            }
            roots
        }
    }
    fn isorefine(
        p: &Natural,
        f: &Polynomial<Integer>,
        alpha: &Natural,
        i: &Natural,
        beta: &Natural,
    ) -> Vec<PAdicIntegerVal0Ball> {
        debug_assert!(!f.is_zero());
        debug_assert!(f.is_squarefree());
        debug_assert!(is_prime(p));
        if padic_int_valuation(p, f.clone().derivative().evaluate(&Integer::from(i)))
            <= Valuation::Finite(Integer::from(beta))
        {
            return isorefine1(p, f, alpha, i, beta);
        }
        let mut roots = vec![];
        if beta < alpha {
            let beta_plus_one = beta + Natural::ONE;
            let mut k = Natural::ZERO;
            while &k < p {
                // k = 0, ..., p-1
                roots.append(&mut isorefine(
                    p,
                    f,
                    alpha,
                    &(i + &k * p.nat_pow(&beta_plus_one)),
                    &beta_plus_one,
                ));
                k += Natural::ONE;
            }
        }
        roots
    }
    fn isorefine1(
        p: &Natural,
        f: &Polynomial<Integer>,
        alpha: &Natural,
        i: &Natural,
        beta: &Natural,
    ) -> Vec<PAdicIntegerVal0Ball> {
        debug_assert!(!f.is_zero());
        debug_assert!(f.is_squarefree());
        debug_assert!(is_prime(p));
        let val_f_at_i = padic_int_valuation(p, f.evaluate(&Integer::from(i)));
        if !(Valuation::Finite(Integer::from(beta)) < val_f_at_i) {
            return vec![];
        }
        if Valuation::Finite(Integer::from(Natural::TWO * beta)) < val_f_at_i {
            return vec![PAdicIntegerVal0Ball {
                center: Integer::from(i),
                valuation_radius: beta.clone(),
            }];
        }
        let mut roots = vec![];
        if beta < alpha {
            let beta_plus_one = beta + Natural::ONE;
            let mut k = Natural::ONE;
            while &k < p {
                // k = 1, ..., p-1
                roots.append(&mut isorefine1(
                    p,
                    f,
                    alpha,
                    &(i + &k * p.nat_pow(&beta_plus_one)),
                    &beta_plus_one,
                ));
                k += Natural::ONE;
            }
        }
        roots
    }

    #[cfg(test)]
    mod tests {
        use crate::structure::elements::*;

        use super::*;

        #[test]
        fn test_isolatebf0() {
            let x = Polynomial::<Integer>::var().into_ergonomic();
            let f = (x.pow(2) - 1).into_verbose();
            let p = Natural::from(2u32);
            println!("f = {}", f);
            println!("{:?}", isolatebf0(&p, &f));
            assert_eq!(isolatebf0(&p, &f).len(), 2);
        }
        #[test]
        fn test_isolate0() {
            let x = Polynomial::<Integer>::var().into_ergonomic();
            let f = (x.pow(2) - 1).into_verbose();
            let p = Natural::from(2u32);
            println!("f = {}", f);
            println!("{:?}", isolate0(&p, &f));
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
