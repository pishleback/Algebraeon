use super::{Polynomial, polynomial_structure::*};
use crate::structure::*;
use algebraeon_structures::*;

#[derive(Debug, Clone)]
enum HenselFactorizationNodeCases<
    const LIFTED_BEZOUT_COEFFS: bool,
    RS: EuclideanDomainSignature + GreatestCommonDivisorSignature + FactoringMonoidSignature,
> {
    Leaf,
    Branch {
        // let alpha be the leading coefficient of `h`. Then
        // alpha is invertible modulo `i` (thus also any power of `i`)
        // `h = alpha*f*g mod i^n`
        // `af + bg = 1 mod i`, or if LIFTED_BEZOUT_COEFFS = true, `af + bg = 1 mod i^n`
        f_factorization: Box<HenselFactorizationNode<LIFTED_BEZOUT_COEFFS, RS>>, // defined modulo i^n
        g_factorization: Box<HenselFactorizationNode<LIFTED_BEZOUT_COEFFS, RS>>, // defined modulo i^n
        a: Polynomial<RS::Elem>,                                                 // defined modulo i
        b: Polynomial<RS::Elem>,                                                 // defined modulo i
    },
}

#[derive(Debug, Clone)]
struct HenselFactorizationNode<
    const LIFTED_BEZOUT_COEFFS: bool,
    RS: EuclideanDomainSignature + GreatestCommonDivisorSignature + FactoringMonoidSignature,
> {
    h: Polynomial<RS::Elem>,
    factorization: HenselFactorizationNodeCases<LIFTED_BEZOUT_COEFFS, RS>,
}

/// A factorization of a polynomial `h` as a product `h = alpha*f_1*f_2*...*f_k mod i^n` for some `k >= 1` where `f_1, f_2, ..., f_n` are monic and `alpha` is the leading coefficient of `h` and is non-zero modulo `i`
#[derive(Debug, Clone)]
pub struct HenselFactorization<
    const LIFTED_BEZOUT_COEFFS: bool,
    RS: EuclideanDomainSignature + GreatestCommonDivisorSignature + FactoringMonoidSignature,
> {
    ring: RS,
    i: RS::Elem,
    n: Natural,
    factorization: HenselFactorizationNode<LIFTED_BEZOUT_COEFFS, RS>, //defined absolutely and factored modulo i^n
}

impl<
    const LIFTED_BEZOUT_COEFFS: bool,
    RS: EuclideanDomainSignature + GreatestCommonDivisorSignature + FactoringMonoidSignature,
> HenselFactorization<LIFTED_BEZOUT_COEFFS, RS>
{
    pub fn bezout_coeff_modulus_base(&self) -> &RS::Elem {
        &self.i
    }

    pub fn bezout_coeff_modulus_power(&self) -> Natural {
        if LIFTED_BEZOUT_COEFFS {
            self.n.clone()
        } else {
            Natural::ONE
        }
    }

    pub fn factorization_modulus_base(&self) -> &RS::Elem {
        &self.i
    }

    pub fn factorization_modulus_power(&self) -> &Natural {
        &self.n
    }
}

impl<
    const LIFTED_BEZOUT_COEFFS: bool,
    RS: EuclideanDomainSignature + GreatestCommonDivisorSignature + FactoringMonoidSignature,
> HenselFactorizationNodeCases<LIFTED_BEZOUT_COEFFS, RS>
{
    #[allow(unused)]
    fn check(
        &self,
        ring: &RS,
        h: &Polynomial<RS::Elem>,
        i: &RS::Elem,
        n: &Natural,
    ) -> Result<(), &'static str> {
        match self {
            HenselFactorizationNodeCases::Leaf => {}
            HenselFactorizationNodeCases::Branch {
                f_factorization,
                g_factorization,
                a,
                b,
            } => {
                f_factorization.check(ring, i, n)?;
                g_factorization.check(ring, i, n)?;

                let poly_ring = ring.polynomials();
                let ring_mod_i = ring.euclidean_quotient_ring(i.clone()).unwrap();
                let poly_ring_mod_i = ring_mod_i.polynomials();
                let poly_ring_mod_i_tothe_n = ring
                    .euclidean_quotient_ring(ring.nat_pow(i, n))
                    .unwrap()
                    .into_polynomials();

                //af + bg = 1 mod i
                if !poly_ring_mod_i.is_zero(&poly_ring_mod_i.sum(&[
                    poly_ring_mod_i.mul(a, &f_factorization.h),
                    poly_ring_mod_i.mul(b, &g_factorization.h),
                    poly_ring_mod_i.neg(&poly_ring_mod_i.one()),
                ])) {
                    return Err("af + bg != 1 mod i");
                }

                //af + bg = 1 mod i^n
                #[allow(clippy::collapsible_if)]
                if LIFTED_BEZOUT_COEFFS {
                    if !poly_ring_mod_i_tothe_n.is_zero(&poly_ring_mod_i_tothe_n.sum(&[
                        poly_ring_mod_i_tothe_n.mul(a, &f_factorization.h),
                        poly_ring_mod_i_tothe_n.mul(b, &g_factorization.h),
                        poly_ring_mod_i_tothe_n.neg(&poly_ring_mod_i_tothe_n.one()),
                    ])) {
                        return Err("af + bg != 1 mod i^n");
                    }
                }

                // //deg(a) < deg(g) and deg(b) < deg(f)
                // if poly_ring.degree(a).unwrap() >= poly_ring.degree(&g_factorization.h).unwrap() {
                //     return Err("deg(a) >= deg(g)");
                // }
                // if poly_ring.degree(b).unwrap() >= poly_ring.degree(&f_factorization.h).unwrap() {
                //     return Err("deg(b) >= deg(f)");
                // }

                //h = alpha*f*g mod i^n
                if !poly_ring_mod_i_tothe_n.equal(
                    h,
                    &poly_ring_mod_i_tothe_n.product(&[
                        &Polynomial::constant(poly_ring.leading_coeff(h).unwrap().clone()),
                        &f_factorization.h,
                        &g_factorization.h,
                    ]),
                ) {
                    return Err("h != alpha*f*g mod i^n");
                }
            }
        }
        Ok(())
    }

    fn new_split(
        ring: &RS,
        p: &RS::Elem,
        n: &Natural,
        first_fs: Vec<&Polynomial<RS::Elem>>,
        second_fs: Vec<&Polynomial<RS::Elem>>,
    ) -> Self {
        let poly_ring = ring.polynomials();
        let poly_ring_mod_p = ring.quotient_field_unchecked(p.clone()).into_polynomials();

        //first_h and second_h are defined modulo p^n
        let first_h = poly_ring
            .product(&first_fs)
            .apply_map(|c| ring.rem(c, &ring.nat_pow(p, n)));
        let second_h = poly_ring
            .product(&second_fs)
            .apply_map(|c| ring.rem(c, &ring.nat_pow(p, n)));

        //find a, b such that af + bg = 1 mod i
        let (u, a, b) = poly_ring_mod_p.xgcd(&first_h, &second_h);
        assert!(
            poly_ring_mod_p.equal(&u, &poly_ring_mod_p.one()),
            "Factors should be coprime modulo i"
        );

        Self::Branch {
            f_factorization: Box::new(HenselFactorizationNode::new(ring, p, n, first_h, first_fs)),
            g_factorization: Box::new(HenselFactorizationNode::new(
                ring, p, n, second_h, second_fs,
            )),
            a,
            b,
        }
    }

    /// Build a branch node from two coprime lists of modular factors when the base
    /// factorization is known only modulo the prime `p` (i.e. as elements of the
    /// field `F_p`).
    ///
    /// Forms the two subproducts `first_h ≡ ∏ first_fs` and `second_h ≡ ∏ second_fs`
    /// modulo `p^n`, computes Bézout cofactors `a, b` with `a·first_h + b·second_h
    /// ≡ 1 (mod p)` in `F_p[x]`, lifts those cofactors back to ring elements, and
    /// recurses on each half. This mirrors [`new_split`](Self::new_split) but takes
    /// an explicit field `F_p` instead of assuming the ring carries a field-quotient
    /// view, so it works with backends such as Montgomery form.
    fn new_split_mod_field<
        FieldModP: QuotientRingGetPrincipalIdealSignature<RS> + FieldSignature,
    >(
        ring: &RS,
        field_mod_p: &FieldModP,
        n: &Natural,
        first_fs: Vec<&Polynomial<RS::Elem>>,
        second_fs: Vec<&Polynomial<RS::Elem>>,
    ) -> Self {
        let poly_ring = ring.polynomials();
        let poly_ring_mod_p = field_mod_p.polynomials();
        let p = field_mod_p.modulus();

        let first_h = poly_ring
            .product(&first_fs)
            .apply_map(|c| ring.rem(c, &ring.nat_pow(p.as_ref(), n)));
        let second_h = poly_ring
            .product(&second_fs)
            .apply_map(|c| ring.rem(c, &ring.nat_pow(p.as_ref(), n)));
        let first_h_mod_p = first_h.apply_map(|c| field_mod_p.project_ref(c));
        let second_h_mod_p = second_h.apply_map(|c| field_mod_p.project_ref(c));
        let (u, a, b) = poly_ring_mod_p.xgcd(&first_h_mod_p, &second_h_mod_p);
        assert!(
            poly_ring_mod_p.equal(&u, &poly_ring_mod_p.one()),
            "Factors should be coprime modulo p"
        );

        Self::Branch {
            f_factorization: Box::new(HenselFactorizationNode::new_mod_field(
                ring,
                field_mod_p,
                n,
                first_h,
                first_fs,
            )),
            g_factorization: Box::new(HenselFactorizationNode::new_mod_field(
                ring,
                field_mod_p,
                n,
                second_h,
                second_fs,
            )),
            a: a.apply_map(|c| field_mod_p.unproject_ref(c)),
            b: b.apply_map(|c| field_mod_p.unproject_ref(c)),
        }
    }

    fn factor_list<'a>(&'a self, h: &'a Polynomial<RS::Elem>) -> Vec<&'a Polynomial<RS::Elem>> {
        match self {
            HenselFactorizationNodeCases::Leaf => vec![h],
            HenselFactorizationNodeCases::Branch {
                f_factorization,
                g_factorization,
                a: _,
                b: _,
            } => {
                let mut fs = vec![];
                fs.extend(f_factorization.factor_list());
                fs.extend(g_factorization.factor_list());
                fs
            }
        }
    }
}

impl<RS: EuclideanDomainSignature + GreatestCommonDivisorSignature + FactoringMonoidSignature>
    HenselFactorizationNodeCases<true, RS>
{
    fn dont_lift_bezout_coeffs(self) -> HenselFactorizationNodeCases<false, RS> {
        match self {
            HenselFactorizationNodeCases::Leaf => HenselFactorizationNodeCases::Leaf,
            HenselFactorizationNodeCases::Branch {
                f_factorization,
                g_factorization,
                a,
                b,
            } => HenselFactorizationNodeCases::Branch {
                f_factorization: Box::new(f_factorization.dont_lift_bezout_coeffs()),
                g_factorization: Box::new(g_factorization.dont_lift_bezout_coeffs()),
                a,
                b,
            },
        }
    }
}

/// Remainder of `dividend` on division by the monic polynomial `monic_divisor`,
/// reducing every coefficient modulo `modulus` at each step.
///
/// Computing this remainder by exact division over the integers makes the
/// intermediate quotient/remainder coefficients grow geometrically with the
/// degree (for a degree-`d` divisor with large coefficients they can reach
/// ~`d * bits(modulus)` bits), which dominates the Hensel lift for high-degree
/// inputs. Reducing modulo `modulus` after every multiply keeps every value
/// bounded, and the lift only ever needs the result modulo `modulus` anyway.
fn rem_by_monic_mod<RS: EuclideanDomainSignature>(
    ring: &RS,
    dividend: &Polynomial<RS::Elem>,
    monic_divisor: &Polynomial<RS::Elem>,
    modulus: &RS::Elem,
) -> Polynomial<RS::Elem> {
    let poly_ring = ring.polynomials();
    let div_deg = poly_ring.degree(monic_divisor).unwrap();
    let Some(dividend_deg) = poly_ring.degree(dividend) else {
        return poly_ring.zero();
    };
    if dividend_deg < div_deg {
        return dividend.apply_map(|c| ring.rem(c, modulus));
    }

    let mut rem: Vec<RS::Elem> = (0..=dividend_deg)
        .map(|j| poly_ring.coeff(dividend, j).as_ref().clone())
        .collect();
    let div_coeffs: Vec<RS::Elem> = (0..div_deg)
        .map(|j| ring.rem(poly_ring.coeff(monic_divisor, j).as_ref(), modulus))
        .collect();

    for i in (div_deg..=dividend_deg).rev() {
        // `monic_divisor` is monic, so the next quotient coefficient is just the
        // current leading remainder coefficient. Reduce only this pivot: the lower
        // coefficients each accumulate at most `div_deg` subtractions of size
        // `modulus^2` before they themselves become a pivot, so they stay bounded
        // without a reduction on every inner step (reductions are the expensive op).
        let q = ring.rem(&rem[i], modulus);
        rem[i] = ring.zero();
        if ring.is_zero(&q) {
            continue;
        }
        let base = i - div_deg;
        for j in 0..div_deg {
            let prod = ring.mul(&q, &div_coeffs[j]);
            rem[base + j] = ring.add(&rem[base + j], &ring.neg(&prod));
        }
    }

    rem.truncate(div_deg);
    Polynomial::from_coeffs(rem.iter().map(|c| ring.rem(c, modulus)).collect())
}

/// https://en.wikipedia.org/wiki/Hensel%27s_lemma
///
/// Input is such that:
///  - `h = fg mod i^n`
///  - `af + bg = 1 mod i`
///
/// Output (df, dg, f', g') such that:
///  - `h = (f + df)(g + dg) mod i^{n+1}`
///  - `a(f + df) + b(g + dg) = 1 mod i`
///  - `f' = f + df mod i^{n+1}`
///  - `g' = g + dg mod i^{n+1}`
#[allow(clippy::too_many_arguments)]
fn compute_lift_factors<
    RS: EuclideanDomainSignature + GreatestCommonDivisorSignature + FactoringMonoidSignature,
>(
    ring: &RS,
    i: &RS::Elem,
    n: &Natural,
    a: &Polynomial<RS::Elem>,
    b: &Polynomial<RS::Elem>,
    f: &Polynomial<RS::Elem>,
    g: &Polynomial<RS::Elem>,
    h: &Polynomial<RS::Elem>,
) -> (
    Polynomial<RS::Elem>,
    Polynomial<RS::Elem>,
    Polynomial<RS::Elem>,
    Polynomial<RS::Elem>,
) {
    let poly_ring = ring.polynomials();
    let modulus = ring.nat_pow(i, &(n + Natural::ONE));

    let alpha = poly_ring.leading_coeff(h).unwrap();
    let (gcd, beta, gamma) = ring.euclidean_xgcd(alpha.clone(), i.clone());
    debug_assert!(ring.equal(&gcd, &ring.one()));
    drop(gcd);
    drop(gamma);

    let ring_mod_i = ring.euclidean_quotient_ring(i.clone()).unwrap();
    debug_assert!(ring_mod_i.equal(alpha, poly_ring.leading_coeff(h).unwrap()));

    let delta_h = poly_ring
        .add(
            h,
            &poly_ring.neg(&poly_ring.product(&[&Polynomial::constant(alpha.clone()), f, g])),
        )
        .apply_map(|c| ring.rem(c, &modulus));

    //found delta_h such that
    //delta_h = h - alpha*f*g mod i^n+1
    debug_assert!(poly_ring.degree(&delta_h).unwrap_or(0) < poly_ring.degree(h).unwrap());

    // rg = (a * delta_h) mod g, rf = (b * delta_h) mod f, all modulo i^{n+1}.
    // Only the remainders are needed; the quotients are discarded. Performing the
    // division modulo i^{n+1} avoids catastrophic integer-coefficient blow-up.
    let rg = rem_by_monic_mod(ring, &poly_ring.mul(a, &delta_h), g, &modulus);
    let rf = rem_by_monic_mod(ring, &poly_ring.mul(b, &delta_h), f, &modulus);

    let delta_f = poly_ring
        .mul(&Polynomial::constant(beta.clone()), &rf)
        .apply_map(|c| ring.rem(c, &modulus));
    let delta_g = poly_ring
        .mul(&Polynomial::constant(beta.clone()), &rg)
        .apply_map(|c| ring.rem(c, &modulus));

    let lifted_f = poly_ring
        .add(&delta_f, f)
        .apply_map(|c| ring.rem(c, &modulus));
    let lifted_g = poly_ring
        .add(&delta_g, g)
        .apply_map(|c| ring.rem(c, &modulus));

    (delta_f, delta_g, lifted_f, lifted_g)
}

impl<RS: EuclideanDomainSignature + GreatestCommonDivisorSignature + FactoringMonoidSignature>
    HenselFactorizationNodeCases<false, RS>
{
    fn linear_lift(&mut self, ring: &RS, i: &RS::Elem, n: &Natural, h: &Polynomial<RS::Elem>) {
        match self {
            HenselFactorizationNodeCases::Leaf => {}
            HenselFactorizationNodeCases::Branch {
                f_factorization,
                g_factorization,
                a,
                b,
            } => {
                let f = &f_factorization.h;
                let g = &g_factorization.h;
                let (_, _, lifted_f, lifted_g) = compute_lift_factors(ring, i, n, a, b, f, g, h);
                f_factorization.h = lifted_f;
                g_factorization.h = lifted_g;

                f_factorization.linear_lift(ring, i, n);
                g_factorization.linear_lift(ring, i, n);
            }
        }

        #[cfg(debug_assertions)]
        self.check(ring, h, i, &(n + Natural::ONE)).unwrap();
    }
}

impl<RS: EuclideanDomainSignature + GreatestCommonDivisorSignature + FactoringMonoidSignature>
    HenselFactorizationNodeCases<true, RS>
{
    fn quadratic_lift(&mut self, ring: &RS, i: &RS::Elem, n: &Natural, h: &Polynomial<RS::Elem>) {
        match self {
            HenselFactorizationNodeCases::Leaf => {}
            HenselFactorizationNodeCases::Branch {
                f_factorization,
                g_factorization,
                a,
                b,
            } => {
                let pring_mod_i2n = ring
                    .euclidean_quotient_ring(ring.nat_pow(i, &(n * Natural::TWO)))
                    .unwrap()
                    .into_polynomials();

                let f = &f_factorization.h;
                let g = &g_factorization.h;
                let (delta_f, delta_g, lifted_f, lifted_g) =
                    compute_lift_factors(ring, &ring.nat_pow(i, n), &Natural::ONE, a, b, f, g, h);

                // beta = af + bg - 1 mod i^n
                let beta = pring_mod_i2n.sum(&[
                    pring_mod_i2n.mul(a, f),
                    pring_mod_i2n.mul(b, g),
                    pring_mod_i2n.neg(&pring_mod_i2n.one()),
                ]);

                // big_delta = beta + a * delta_f + b * delta_g mod i^n
                let big_delta = pring_mod_i2n.sum(&[
                    beta,
                    pring_mod_i2n.mul(a, &delta_f),
                    pring_mod_i2n.mul(b, &delta_g),
                ]);

                // a * lifted_f + b * lifted_g = 1 + big_delta
                debug_assert!(pring_mod_i2n.equal(
                    &pring_mod_i2n.add(
                        &pring_mod_i2n.mul(a, &lifted_f),
                        &pring_mod_i2n.mul(b, &lifted_g)
                    ),
                    &pring_mod_i2n.add(&big_delta, &pring_mod_i2n.one())
                ));

                // delta_a = -a * big_delta
                // delta_b = -b * big_delta
                let delta_a = pring_mod_i2n.neg(&pring_mod_i2n.mul(a, &big_delta));
                let delta_b = pring_mod_i2n.neg(&pring_mod_i2n.mul(b, &big_delta));

                f_factorization.h = lifted_f;
                g_factorization.h = lifted_g;
                *a = pring_mod_i2n.add(a, &delta_a);
                *b = pring_mod_i2n.add(b, &delta_b);

                f_factorization.quadratic_lift(ring, i, n);
                g_factorization.quadratic_lift(ring, i, n);
            }
        }

        #[cfg(debug_assertions)]
        self.check(ring, h, i, &(n * Natural::TWO)).unwrap();
    }
}

impl<
    const LIFTED_BEZOUT_COEFFS: bool,
    RS: EuclideanDomainSignature + GreatestCommonDivisorSignature + FactoringMonoidSignature,
> HenselFactorizationNode<LIFTED_BEZOUT_COEFFS, RS>
{
    #[allow(unused)]
    fn check(&self, ring: &RS, i: &RS::Elem, n: &Natural) -> Result<(), &'static str> {
        // let poly_ring = PolynomialStructure::new(ring.clone().into());
        // if !poly_ring.is_monic(&self.h) {
        //     return Err("h is not monic");
        // }
        self.factorization.check(ring, &self.h, i, n)?;
        Ok(())
    }

    fn new(
        ring: &RS,
        p: &RS::Elem,
        n: &Natural,
        h: Polynomial<RS::Elem>,
        mut fs: Vec<&Polynomial<RS::Elem>>,
    ) -> Self {
        debug_assert!(!fs.is_empty());
        match fs.len() {
            0 => panic!(),
            1 => Self {
                h,
                factorization: HenselFactorizationNodeCases::Leaf,
            },
            fs_len => {
                debug_assert!(fs_len >= 2);
                let second_fs = fs.split_off(fs_len / 2);
                let first_fs = fs;
                debug_assert!(!first_fs.is_empty());
                debug_assert!(!second_fs.is_empty());
                debug_assert_eq!(first_fs.len() + second_fs.len(), fs_len);

                //find an inverse beta to alpha modulo p
                let alpha = ring.polynomials().leading_coeff(&h).unwrap();
                let (g, _beta, _gamma) = ring.euclidean_xgcd(alpha.clone(), p.clone());
                debug_assert!(ring.equal(&g, &ring.one()));

                Self {
                    h,
                    factorization: HenselFactorizationNodeCases::new_split(
                        ring, p, n, first_fs, second_fs,
                    ),
                }
            }
        }
    }

    /// Recursively build a product-tree node for the factors `fs`, each known
    /// modulo the prime field `F_p`, by splitting the list in half at every level.
    /// A single factor becomes a leaf; otherwise the two halves are combined by
    /// [`new_split_mod_field`](HenselFactorizationNodeCases::new_split_mod_field).
    /// This is the field-quotient counterpart of [`new`](Self::new).
    fn new_mod_field<FieldModP: QuotientRingGetPrincipalIdealSignature<RS> + FieldSignature>(
        ring: &RS,
        field_mod_p: &FieldModP,
        n: &Natural,
        h: Polynomial<RS::Elem>,
        mut fs: Vec<&Polynomial<RS::Elem>>,
    ) -> Self {
        debug_assert!(!fs.is_empty());
        match fs.len() {
            0 => panic!(),
            1 => Self {
                h,
                factorization: HenselFactorizationNodeCases::Leaf,
            },
            fs_len => {
                let second_fs = fs.split_off(fs_len / 2);
                let first_fs = fs;
                Self {
                    h,
                    factorization: HenselFactorizationNodeCases::new_split_mod_field(
                        ring,
                        field_mod_p,
                        n,
                        first_fs,
                        second_fs,
                    ),
                }
            }
        }
    }

    fn factor_list(&self) -> Vec<&Polynomial<RS::Elem>> {
        self.factorization.factor_list(&self.h)
    }
}

impl<RS: EuclideanDomainSignature + GreatestCommonDivisorSignature + FactoringMonoidSignature>
    HenselFactorizationNode<true, RS>
{
    fn dont_lift_bezout_coeffs(self) -> HenselFactorizationNode<false, RS> {
        HenselFactorizationNode {
            h: self.h,
            factorization: self.factorization.dont_lift_bezout_coeffs(),
        }
    }
}

impl<RS: EuclideanDomainSignature + GreatestCommonDivisorSignature + FactoringMonoidSignature>
    HenselFactorizationNode<false, RS>
{
    fn linear_lift(&mut self, ring: &RS, i: &RS::Elem, n: &Natural) {
        self.factorization.linear_lift(ring, i, n, &self.h);
    }
}

impl<RS: EuclideanDomainSignature + GreatestCommonDivisorSignature + FactoringMonoidSignature>
    HenselFactorizationNode<true, RS>
{
    fn quadratic_lift(&mut self, ring: &RS, i: &RS::Elem, n: &Natural) {
        self.factorization.quadratic_lift(ring, i, n, &self.h);
    }
}

impl<
    const LIFTED_BEZOUT_COEFFS: bool,
    RS: EuclideanDomainSignature + GreatestCommonDivisorSignature + FactoringMonoidSignature,
> HenselFactorization<LIFTED_BEZOUT_COEFFS, RS>
{
    #[allow(unused)]
    fn check(&self) -> Result<(), &'static str> {
        self.factorization.check(&self.ring, &self.i, &self.n)
    }

    pub fn new(
        ring: RS,
        p: RS::Elem,
        n: Natural,
        h: Polynomial<RS::Elem>,
        fs: Vec<Polynomial<RS::Elem>>,
    ) -> Self {
        let poly_ring = ring.polynomials();

        debug_assert!(ring.is_irreducible(&p));

        // h and all fs monic
        assert!(!fs.is_empty());
        // debug_assert!(poly_ring.is_monic(&h));
        for f in &fs {
            debug_assert!(poly_ring.is_monic(f));
        }
        // h = product of fs modulo i^n
        let poly_ring_mod_p_tothe_n = ring
            .euclidean_quotient_ring(ring.nat_pow(&p, &n))
            .unwrap()
            .into_polynomials();
        let alpha = poly_ring.leading_coeff(&h).unwrap();
        debug_assert!(poly_ring_mod_p_tothe_n.equal(
            &h,
            &poly_ring_mod_p_tothe_n.mul(
                &Polynomial::constant(alpha.clone()),
                &poly_ring_mod_p_tothe_n.product(&fs.iter().collect::<Vec<_>>())
            )
        ));
        // fs are coprime mod i - checked when computing bezout coefficients
        let ans_ring = ring.clone();
        let factors = HenselFactorizationNode::new(&ring, &p, &n, h, fs.iter().collect());
        let ans = Self {
            ring: ans_ring,
            i: p,
            n,
            factorization: factors,
        };
        #[cfg(debug_assertions)]
        ans.check().unwrap();
        ans
    }

    pub fn modulus(&self) -> RS::Elem {
        self.ring.nat_pow(&self.i, &self.n)
    }

    //return the lifted factors in order
    pub fn factors(&self) -> Vec<&Polynomial<RS::Elem>> {
        self.factorization.factor_list()
    }
}

impl<RS: EuclideanDomainSignature + GreatestCommonDivisorSignature + FactoringMonoidSignature>
    HenselFactorization<true, RS>
{
    /// Build a Hensel factorization of `h` from its monic irreducible factors `fs`
    /// taken modulo the prime field `F_p`, so the initial modulus is `p^1`.
    ///
    /// The factors must satisfy `h ≡ lc(h) · ∏ fs (mod p)` and be pairwise coprime
    /// in `F_p[x]`; the Bézout cofactors of the product tree are computed in `F_p`.
    /// The result can then be quadratically (or linearly) lifted to any `p^k`. This
    /// is the product-tree counterpart of [`new`](Self::new), which instead takes a
    /// base factorization already given modulo `p^n`.
    pub fn new_from_mod_field_factors<
        FieldModP: QuotientRingGetPrincipalIdealSignature<RS> + FieldSignature,
    >(
        ring: RS,
        field_mod_p: FieldModP,
        h: Polynomial<RS::Elem>,
        fs: Vec<Polynomial<RS::Elem>>,
    ) -> Self {
        let p = field_mod_p.modulus().into_owned();
        let n = Natural::ONE;
        let factors =
            HenselFactorizationNode::new_mod_field(&ring, &field_mod_p, &n, h, fs.iter().collect());
        let ans = Self {
            ring,
            i: p,
            n,
            factorization: factors,
        };
        #[cfg(debug_assertions)]
        ans.check().unwrap();
        ans
    }

    pub fn dont_lift_bezout_coeffs(self) -> HenselFactorization<false, RS> {
        // When LIFTED_BEZOUT_COEFFS = true, the pair (self.i, self.n) means the factorization is modulo i^n and the bezout coeffs are modulo i^n
        // When LIFTED_BEZOUT_COEFFS = false, the pair (self.i, self.n) means the factorization is modulo i^n and the bezout coeffs are modulo i
        // So when LIFTED_BEZOUT_COEFFS goes from true -> false it's not logically incorrect to send (i, n) -> (i, n), but it's better to send (i, n) -> (i^n, 1)
        let i = self.ring.nat_pow(&self.i, &self.n);
        let n = Natural::ONE;
        HenselFactorization {
            ring: self.ring,
            i,
            n,
            factorization: self.factorization.dont_lift_bezout_coeffs(),
        }
    }
}

impl<RS: EuclideanDomainSignature + GreatestCommonDivisorSignature + FactoringMonoidSignature>
    HenselFactorization<false, RS>
{
    /// update the Hensel factorization sending (i, i^n) -> (i, i^{n+1})
    pub fn linear_lift(&mut self) {
        self.factorization.linear_lift(&self.ring, &self.i, &self.n);
        self.n += Natural::ONE;
    }
}

impl<RS: EuclideanDomainSignature + GreatestCommonDivisorSignature + FactoringMonoidSignature>
    HenselFactorization<true, RS>
{
    /// update the Hensel factorization sending (i, i^n) -> (i, i^{2n})
    pub fn quadratic_lift(&mut self) {
        self.factorization
            .quadratic_lift(&self.ring, &self.i, &self.n);
        self.n *= Natural::TWO;
    }
}

impl<
    RS: FactoringMonoidSignature<FactoredExponent = NaturalCanonicalStructure>
        + EuclideanDomainSignature
        + GreatestCommonDivisorSignature,
    RSB: BorrowedStructure<RS>,
    RSQB: BorrowedStructure<EuclideanRemainderQuotientStructure<RS, RSB, true>>,
    RSQPB: BorrowedStructure<
        PolynomialStructure<EuclideanRemainderQuotientStructure<RS, RSB, true>, RSQB>,
    >,
    NB: BorrowedStructure<NaturalCanonicalStructure>,
>
    FactoringStructure<
        PolynomialStructure<EuclideanRemainderQuotientStructure<RS, RSB, true>, RSQB>,
        RSQPB,
        NaturalCanonicalStructure,
        NB,
    >
where
    PolynomialStructure<EuclideanRemainderQuotientStructure<RS, RSB, true>, RSQB>:
        SetSignature<Elem = Polynomial<RS::Elem>>
            + FactoringMonoidSignature<FactoredExponent = NaturalCanonicalStructure>,
{
    /// If the polynomial is squarefree return a hensel factorization, otherwise return None
    pub fn into_hensel_factorization(
        &self,
        a: NonZeroFactored<Polynomial<RS::Elem>, Natural>,
        h: Polynomial<RS::Elem>,
    ) -> Option<HenselFactorization<true, RS>>
    where
        RS: EuclideanDomainSignature + GreatestCommonDivisorSignature,
    {
        let poly_ring_mod = self.objects().clone();
        let ring_mod = poly_ring_mod.coeff_ring();
        let ring = ring_mod.ring().clone();
        // let poly_ring: PolynomialStructure<RS> = PolynomialStructure::new(ring.clone());

        let mut fs = vec![];
        let (_unit, factors) = a.into_unit_and_powers();
        for (factor, power) in factors {
            if power == Natural::ONE {
                fs.push(factor);
            } else {
                return None;
            }
        }

        let hensel_factorization =
            HenselFactorization::new(ring, ring_mod.modulus().into_owned(), Natural::ONE, h, fs);
        #[cfg(debug_assertions)]
        hensel_factorization.check().unwrap();
        Some(hensel_factorization)
    }

    /// If the polynomial is squarefree return a hensel factorization, otherwise return None
    pub fn into_hensel_factorization_unchecked(
        &self,
        a: NonZeroFactored<Polynomial<RS::Elem>, Natural>,
        h: Polynomial<RS::Elem>,
    ) -> HenselFactorization<true, RS>
    where
        RS: EuclideanDomainSignature + GreatestCommonDivisorSignature,
    {
        let poly_ring_mod = self.objects().clone();
        let ring_mod = poly_ring_mod.coeff_ring();
        let ring = ring_mod.ring().clone();

        let mut fs = vec![];
        let (_unit, factors) = a.into_unit_and_powers();
        for (factor, power) in factors {
            debug_assert!(power == Natural::ONE);
            fs.push(factor);
        }

        let hensel_factorization =
            HenselFactorization::new(ring, ring_mod.modulus().into_owned(), Natural::ONE, h, fs);
        #[cfg(debug_assertions)]
        hensel_factorization.check().unwrap();
        hensel_factorization
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hensel_factorization_example1_mod5() {
        let f1 =
            Polynomial::from_coeffs(vec![Integer::from(1), Integer::from(-1), Integer::from(1)]);
        let f2 = Polynomial::from_coeffs(vec![
            Integer::from(4),
            Integer::from(3),
            Integer::from(2),
            Integer::from(1),
            Integer::from(1),
        ]);
        let f3 = Polynomial::from_coeffs(vec![
            Integer::from(4),
            Integer::from(3),
            Integer::from(2),
            Integer::from(1),
        ]);
        let h = Polynomial::product(&[&Polynomial::constant(Integer::from(8)), &f1, &f2, &f3]);
        let h = h.apply_map(|c| Integer::rem(c, &Integer::from(5)));
        //h = prod fs mod 5
        //fs are coprime
        //h and fs are monic

        //set up bezout coefficients for hensel lifting the factorization modulo 25, 125, ...
        let mut hensel_fact = HenselFactorization::<true, _>::new(
            Integer::structure(),
            Integer::from(5),
            Natural::from(1u8),
            h.clone(),
            vec![f1, f2, f3],
        )
        .dont_lift_bezout_coeffs();
        hensel_fact.check().unwrap();
        println!("5^1: {:?}", hensel_fact.factors());
        for i in 2..20 {
            hensel_fact.linear_lift();
            hensel_fact.check().unwrap();
            println!("5^{}: {:?}", i, hensel_fact.factors());
            let lifted_product = Polynomial::mul(
                &Polynomial::constant(Polynomial::leading_coeff(&h).unwrap().clone()),
                &Polynomial::product(&hensel_fact.factors()),
            )
            .apply_map(|c| Integer::rem(c, &hensel_fact.modulus()));
            println!("{:?} {:?}", lifted_product, h);
            assert_eq!(lifted_product, h);
        }
    }

    #[test]
    fn rem_by_monic_mod_matches_exact_division() {
        let ring = Integer::structure();
        // x^3 + 2x + 5 = x * (x^2 + 1) + (x + 5), so the remainder is x + 5.
        let dividend = Polynomial::from_coeffs(vec![
            Integer::from(5),
            Integer::from(2),
            Integer::from(0),
            Integer::from(1),
        ]);
        let monic_divisor =
            Polynomial::from_coeffs(vec![Integer::from(1), Integer::from(0), Integer::from(1)]);
        let modulus = Integer::from(1000);
        let rem = rem_by_monic_mod(&ring, &dividend, &monic_divisor, &modulus);
        assert_eq!(
            rem,
            Polynomial::from_coeffs(vec![Integer::from(5), Integer::from(1)])
        );

        // x^2 mod (x - 1) = 1, with coefficients reduced mod 7.
        let rem2 = rem_by_monic_mod(
            &ring,
            &Polynomial::from_coeffs(vec![Integer::from(0), Integer::from(0), Integer::from(1)]),
            &Polynomial::from_coeffs(vec![Integer::from(-1), Integer::from(1)]),
            &Integer::from(7),
        );
        assert_eq!(rem2, Polynomial::constant(Integer::from(1)));

        // A dividend of smaller degree is returned unchanged (mod the modulus).
        let small = Polynomial::from_coeffs(vec![Integer::from(3), Integer::from(2)]);
        let rem3 = rem_by_monic_mod(&ring, &small, &monic_divisor, &modulus);
        assert_eq!(rem3, small);
    }

    #[test]
    fn new_from_mod_field_factors_lifts_to_higher_powers() {
        // h = (x^2 + 2)(x + 1)(x + 4) = x^4 + x^2 + 3 modulo 5, with three pairwise
        // coprime monic factors over F_5.
        let field_mod_5 = Integer::structure().into_quotient_field_unchecked(Integer::from(5));
        let f1 =
            Polynomial::from_coeffs(vec![Integer::from(2), Integer::from(0), Integer::from(1)]);
        let f2 = Polynomial::from_coeffs(vec![Integer::from(1), Integer::from(1)]);
        let f3 = Polynomial::from_coeffs(vec![Integer::from(4), Integer::from(1)]);
        let h = Polynomial::from_coeffs(vec![
            Integer::from(3),
            Integer::from(0),
            Integer::from(1),
            Integer::from(0),
            Integer::from(1),
        ]);

        let mut hensel_fact = HenselFactorization::<true, _>::new_from_mod_field_factors(
            Integer::structure(),
            field_mod_5,
            h.clone(),
            vec![f1, f2, f3],
        );
        assert_eq!(hensel_fact.factors().len(), 3);
        assert_eq!(hensel_fact.modulus(), Integer::from(5));

        // Each quadratic lift keeps the product congruent to h modulo the new power.
        for _ in 0..4 {
            hensel_fact.quadratic_lift();
            hensel_fact.check().unwrap();
            let lifted_product = Polynomial::product(&hensel_fact.factors())
                .apply_map(|c| Integer::rem(c, &hensel_fact.modulus()));
            assert_eq!(lifted_product, h);
        }
    }
}
