use super::*;
use algebraeon_sets::structure::MetaType;
use std::collections::HashSet;

// Some algorithms here on p-adic root isolation can be found in
// Sturm, Thomas & Weispfenning, Volker. (2004). P-adic Root Isolation. Revista de la Real Academia de Ciencias Exactas, Físicas y Naturales. Serie A, Matemáticas.
// https://www.researchgate.net/profile/Thomas-Sturm-2/publication/2925550_P-adic_Root_Isolation/links/580b8c0708aeef1bfeeb5db8/P-adic-Root-Isolation.pdf?origin=scientificContributions

mod balancable_pairs {
    use std::collections::HashSet;

    use algebraeon_nzq::traits::Abs;

    use super::*;

    /// A balancable pair of a polynomial consists of two of the monomial terms satisfying some conditions.
    ///
    /// For a fixed choice of prime $p$, a balancable pair of a polynomial $$f(x) = f_nx^n + ... + f_2x^2 + f_1x + f_0$$ is
    /// a pair of indices $0 \le i < j \le n$ such that
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
            let min = (0..=self.n)
                .filter_map(
                    |k| match padic_int_valuation(&self.p, self.f.coeff(k).into_owned()) {
                        Valuation::Infinity => None,
                        Valuation::Finite(vfk) => Some(vfk + Integer::from(k) * &self.bv),
                    },
                )
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
        /// $$f_{i,j}(x) = p^{-\frac{iv_p(f_i) - jv_p(f_j)}{j - i}} \sum_{k=0}^n f_k \left(p^{\frac{v_p(f_i) - v_p(f_j)}{j - i}} x\right)^k$$
        /// It is such that for $\alpha \in \mathbb{Q}_p$
        /// - $f_{i,j}(\alpha) = 0$ if and only if $f\left(p^{\frac{v_p(f_j) - v_p(f_i)}{j - i}} \alpha\right) = 0$
        /// - $f_{i,j}^{\prime}(\alpha) = 0$ if and only if $f^{\prime}\left(p^{\frac{v_p(f_j) - v_p(f_i)}{j - i}} \alpha\right) = 0$
        pub fn normalization(&self) -> Polynomial<Integer> {
            debug_assert!(self.f.is_irreducible());
            debug_assert!(self.is_critical());
            let cmbv = self.crossmul_balancing_value();
            Polynomial::from_coeffs(
                (0..=self.n)
                    .map(|k| {
                        let p_pow = self.balancing_value() * Integer::from(k) - &cmbv;
                        // compute f_k*p^p_pow
                        if p_pow >= Integer::ZERO {
                            let p_pow = p_pow.abs();
                            self.f.coeff(k).as_ref() * Integer::from(self.p.nat_pow(&p_pow))
                        } else {
                            let neg_p_pow = (-p_pow).abs();
                            Integer::try_div(
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
            debug_assert!(p.is_prime());
            let mut bps = vec![];
            let n = self.degree().unwrap();
            let coeff_valuations = (0..=n)
                .map(|k| padic_int_valuation(p, self.coeff(k).into_owned()))
                .collect::<Vec<_>>();
            for i in 0..=n {
                for j in (i + 1)..=n {
                    #[allow(clippy::single_match_else)]
                    match (&coeff_valuations[i], &coeff_valuations[j]) {
                        (Valuation::Finite(vfi), Valuation::Finite(vfj)) => {
                            match Integer::try_div(&(vfi - vfj), &Integer::from(j - i)) {
                                Some(bv) => bps.push(BalancablePair {
                                    p: p.clone(),
                                    f: self,
                                    n,
                                    i,
                                    j,
                                    vfi: vfi.clone(),
                                    vfj: vfj.clone(),
                                    bv,
                                }),
                                None => {
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

// Represent all p-adic numbers which differ from s by a valuation up to v
// In terms of digits, a is correct in the first v digits
#[derive(Debug, Clone)]
pub struct PAdicRationalBall {
    a: Rational,
    // Note: this value is one larger than what's in the paper
    // In the paper, v is the valuation corresponding the radius of an _open_ ball
    // Here v, is the (signed) number of correct p-adic digits i.e. the valuation corresponding the radius of a _closed_ ball
    v: Integer,
}

impl PAdicRationalBall {
    pub fn center(&self) -> &Rational {
        &self.a
    }

    pub fn ndigits(&self) -> &Integer {
        &self.v
    }

    pub fn neg(self) -> Self {
        Self {
            a: -self.a,
            v: self.v,
        }
    }

    pub fn add_rat(self, offset: &Rational) -> Self {
        Self {
            a: self.a + offset,
            v: self.v,
        }
    }

    pub fn mul_rat(self, p: &Natural, scale: &Rational) -> Self {
        debug_assert_ne!(scale, &Rational::ZERO);
        Self {
            a: self.a * scale,
            v: self.v + padic_rat_valuation(p, scale.clone()).unwrap_int(),
        }
    }
}

/// A brute-force algorithm to isolate all roots of f of valuation 0
/// Return ([r1, ...], alpha) such that each ri is alpha-close to exactly one valuation 0 root of f and vise-versa
fn isolatebf0(p: &Natural, f: &Polynomial<Integer>) -> Vec<PAdicRationalBall> {
    debug_assert!(!f.is_zero());
    debug_assert!(f.is_squarefree());
    debug_assert!(p.is_prime());
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
            #[allow(clippy::collapsible_if)]
            if &i % p != Natural::ZERO {
                if padic_int_valuation(p, f.evaluate(&Integer::from(&i)))
                    > Valuation::Finite(Integer::from(&two_alpha))
                {
                    let int_i = Integer::from(&i);
                    if !roots.iter().any(|alt_i| {
                        padic_int_valuation(p, &int_i - alt_i)
                            > Valuation::Finite(Integer::from(&alpha))
                    }) {
                        roots.push(int_i);
                    }
                }
            }
            i += Natural::ONE;
        }
        roots
            .into_iter()
            .map(|i| PAdicRationalBall {
                a: Rational::from(i),
                v: Integer::from(&alpha) + Integer::ONE,
            })
            .collect()
    }
}

/// Return ([..., (ri, bi), ...], alpha) such that each ri is bi-close to exactly one valuation 0 root of f and vise-versa.
/// Computes one p-adic digit at a time.
fn isolate0(p: &Natural, f: &Polynomial<Integer>) -> Vec<PAdicRationalBall> {
    debug_assert!(!f.is_zero());
    debug_assert!(f.is_squarefree());
    debug_assert!(p.is_prime());
    let n = f.degree().unwrap();
    if n == 0 {
        vec![]
    } else {
        let df = f.clone().derivative();
        let mut roots = vec![];
        // disc(f) != 0 since f is squarefree
        let alpha = padic_int_valuation(p, f.clone().discriminant().unwrap()).unwrap_nat();
        let mut i = Natural::ONE;
        while &i < p {
            // i = 1, ..., p-1
            roots.append(&mut isorefine(p, f, &df, &alpha, &i, &Natural::ONE));
            i += Natural::ONE;
        }
        roots
    }
}

fn isorefine(
    p: &Natural,
    f: &Polynomial<Integer>,
    df: &Polynomial<Integer>,
    alpha: &Natural,
    i: &Natural,
    beta: &Natural,
) -> Vec<PAdicRationalBall> {
    debug_assert!(!f.is_zero());
    debug_assert!(f.is_squarefree());
    debug_assert!(p.is_prime());
    debug_assert_eq!(&f.clone().derivative(), df);
    let p_tothe_beta = p.nat_pow(beta);
    let vdfi = padic_int_valuation(
        p,
        Integer::structure()
            .into_quotient_ring(Integer::from(&p_tothe_beta))
            .into_polynomials()
            .evaluate(df, &Integer::from(i)),
    );
    if vdfi < Valuation::Finite(Integer::from(beta)) {
        return isorefine1(p, f, alpha, i, beta);
    }
    let mut roots = vec![];
    if beta <= alpha {
        let beta_plus_one = beta + Natural::ONE;
        let mut k = Natural::ZERO;
        while &k < p {
            // k = 0, ..., p-1
            roots.append(&mut isorefine(
                p,
                f,
                df,
                alpha,
                &(i + &k * &p_tothe_beta),
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
) -> Vec<PAdicRationalBall> {
    debug_assert!(!f.is_zero());
    debug_assert!(f.is_squarefree());
    debug_assert!(p.is_prime());
    let val_f_at_i = padic_int_valuation(p, f.evaluate(&Integer::from(i)));
    if Valuation::Finite(Integer::from(beta)) > val_f_at_i {
        return vec![];
    }
    if Valuation::Finite(Integer::from(Natural::TWO * beta) - Integer::TWO) < val_f_at_i {
        return vec![PAdicRationalBall {
            a: Rational::from(i),
            v: Integer::from(beta),
        }];
    }
    let mut roots = vec![];
    if beta <= alpha {
        let beta_plus_one = beta + Natural::ONE;
        let mut k = Natural::ONE;
        while &k < p {
            // k = 1, ..., p-1
            roots.append(&mut isorefine1(
                p,
                f,
                alpha,
                &(i + &k * p.nat_pow(beta)),
                &beta_plus_one,
            ));
            k += Natural::ONE;
        }
    }
    roots
}

/// Refine an isolated root of valuation 0
fn refine0(
    p: &Natural,
    f: &Polynomial<Integer>,
    mut r: PAdicRationalBall,
    target_beta: &Integer,
) -> PAdicRationalBall {
    debug_assert!(r.ndigits() >= &Integer::ZERO);
    let mut do_hensel_lift = false;
    while &r.v < target_beta {
        let fa = f
            .apply_map(|coeff| Rational::from(coeff))
            .evaluate(r.center());
        let dfa = f
            .apply_map(|coeff| Rational::from(coeff))
            .derivative()
            .evaluate(r.center());

        do_hensel_lift = {
            do_hensel_lift || {
                let vfa = padic_rat_valuation(p, fa.clone());
                let vdfa = padic_rat_valuation(p, dfa.clone());
                vfa > match vdfa {
                    Valuation::Infinity => Valuation::Infinity,
                    Valuation::Finite(vdfa_finite) => Valuation::Finite(vdfa_finite * Integer::TWO),
                }
            }
        };

        if do_hensel_lift {
            let new_a = r.a - fa / dfa;
            let new_v = r.v + Integer::ONE;
            let new_a = PAdicRational::from_rational(p.clone(), new_a)
                .truncate(&new_v)
                .rational_value();

            // Hensel lift
            r = PAdicRationalBall { a: new_a, v: new_v };
        } else {
            // Non Hensel lift
            let rv_plus_one = &r.v + Integer::ONE;
            r = refine0_impl(p, f, r, &rv_plus_one).unwrap();
        }
    }
    r
}

fn refine0_impl(
    p: &Natural,
    f: &Polynomial<Integer>,
    r: PAdicRationalBall,
    target_beta: &Integer,
) -> Option<PAdicRationalBall> {
    debug_assert!(r.ndigits() >= &Integer::ZERO);
    let PAdicRationalBall { a: c, v: beta } = &r;
    let vfc = padic_rat_valuation(p, f.apply_map(|coeff| Rational::from(coeff)).evaluate(c));
    if Valuation::Finite(beta.clone()) > vfc {
        return None;
    }
    if target_beta < beta && Valuation::Finite(Integer::TWO * target_beta - Integer::TWO) < vfc {
        return Some(r);
    }
    let beta_plus_one = beta + Integer::ONE;
    let mut k = Natural::ZERO;
    while &k < p {
        // k = 0, ..., p-1
        if let Some(refined_r) = refine0_impl(
            p,
            f,
            PAdicRationalBall {
                a: c + Rational::from(&k) * Rational::from(p).try_int_pow(beta).unwrap(),
                v: beta_plus_one.clone(),
            },
            target_beta,
        ) {
            return Some(refined_r);
        }
        k += Natural::ONE;
    }
    None
}

pub fn isolate(p: &Natural, f: &Polynomial<Integer>) -> Vec<PAdicRationalBall> {
    debug_assert!(p.is_prime());
    debug_assert!(!f.is_zero());
    debug_assert!(f.is_squarefree());
    debug_assert!(f.degree().unwrap() >= 2);
    // compute a set containing one balanced pair for each critical value of f
    let critical_balanced_pairs = {
        let mut critical_balanced_pairs = vec![];
        let mut critical_values = HashSet::new();
        for balancable_pair in f.balancable_pairs(p) {
            if balancable_pair.is_critical() {
                let critical_value = balancable_pair.balancing_value();
                if !critical_values.contains(critical_value) {
                    critical_values.insert(critical_value.clone());
                    critical_balanced_pairs.push(balancable_pair);
                }
            }
        }
        critical_balanced_pairs
    };
    let mut roots = vec![];
    for critical_balanced_pair in critical_balanced_pairs {
        let kappa = critical_balanced_pair.balancing_value();
        for PAdicRationalBall { a: c, v: beta } in
            isolate0(p, &critical_balanced_pair.normalization())
        {
            roots.push(PAdicRationalBall {
                a: Rational::from(p).try_int_pow(kappa).unwrap() * c,
                v: beta + kappa,
            });
        }
    }
    roots
}

pub fn refine(
    p: &Natural,
    f: &Polynomial<Integer>,
    r: &PAdicRationalBall,
    target_gamma: Integer,
) -> PAdicRationalBall {
    let PAdicRationalBall { a: d, v: gamma } = &r;
    // find a balancable pair for f with critical value v(d)
    let vd = padic_rat_valuation(p, d.clone()).unwrap_int(); // d != 0 since it should come from a finite shift of something of valuation 0
    let bp = (|| {
        for bp in f.balancable_pairs(p) {
            if bp.is_critical() && vd == bp.balancing_value().clone() {
                return bp;
            }
        }
        unreachable!()
    })();
    let g = bp.normalization();
    let c = d * Rational::from(p).try_int_pow(&(-&vd)).unwrap();
    let beta = gamma - &vd;
    let target_beta = &target_gamma - &vd;
    let PAdicRationalBall { a: refined_c, v: _ } =
        refine0(p, &g, PAdicRationalBall { a: c, v: beta }, &target_beta);
    let refined_d = refined_c * Rational::from(p).try_int_pow(&vd).unwrap();
    PAdicRationalBall {
        a: refined_d,
        v: target_gamma,
    }
}

#[cfg(test)]
mod tests {
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
        let f = (x.pow(2) - 17).into_verbose();
        let p = Natural::from(2u32);
        println!("f = {}", f);
        println!("{:?}", isolate0(&p, &f));
        assert_eq!(isolate0(&p, &f).len(), 2);
    }

    #[test]
    fn test_isolate() {
        let x = Polynomial::<Integer>::var().into_ergonomic();
        let f = (256 * x.pow(2) - 17).into_verbose();
        let p = Natural::from(2u32);
        println!("f = {}", f);
        println!("{:?}", isolate(&p, &f));
        assert_eq!(isolate(&p, &f).len(), 2);
    }

    #[test]
    fn test_refine() {
        let x = Polynomial::<Integer>::var().into_ergonomic();
        let f = (256 * x.pow(2) - 17).into_verbose();
        let p = Natural::from(2u32);
        println!("f = {}", f);
        for r in isolate(&p, &f) {
            println!("r = {:?}", r);
            let s = refine(&p, &f, &r, Integer::from(100));
            println!("s = {:?}", s);
        }
    }
}
