use super::polynomial::*;
use crate::structure::*;
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;

#[derive(Debug, Clone)]
enum HenselProduct<
    const LIFTED_BEZOUT_COEFFS: bool,
    RS: EuclideanDivisionStructure + GreatestCommonDivisorStructure + FactorableStructure,
> {
    Leaf,
    Branch {
        //let alpha be the leading coefficient of h. Then
        //alpha is invertible modulo i (thus also any power of i)
        //h = alpha*f*g mod i^n
        //af + bg = 1 mod i
        f_factorization: Box<HenselFactorizationImpl<LIFTED_BEZOUT_COEFFS, RS>>, //defined modulo i^n
        g_factorization: Box<HenselFactorizationImpl<LIFTED_BEZOUT_COEFFS, RS>>, //defined modulo i^n
        a: Polynomial<RS::Set>,                                                  //defined modulo i
        b: Polynomial<RS::Set>,                                                  //defined modulo i
    },
}

#[derive(Debug, Clone)]
struct HenselFactorizationImpl<
    const LIFTED_BEZOUT_COEFFS: bool,
    RS: EuclideanDivisionStructure + GreatestCommonDivisorStructure + FactorableStructure,
> {
    h: Polynomial<RS::Set>,
    factorization: HenselProduct<LIFTED_BEZOUT_COEFFS, RS>,
}

//represent a factorization of a polynomial h(x) as a product alpha f_1(x)f_2(x)...f_k(x) modulo i^n where f_1(x), f_2(x), ..., f_n(x) are monic and alpha is the leading coefficient of h and is non-zero modulo i
#[derive(Debug, Clone)]
pub struct HenselFactorization<
    const LIFTED_BEZOUT_COEFFS: bool,
    RS: EuclideanDivisionStructure + GreatestCommonDivisorStructure + FactorableStructure,
> {
    ring: RS,
    i: RS::Set,
    n: Natural,
    factorization: HenselFactorizationImpl<LIFTED_BEZOUT_COEFFS, RS>, //defined absolutely and factored modulo i^n
}

impl<
    const LIFTED_BEZOUT_COEFFS: bool,
    RS: EuclideanDivisionStructure + GreatestCommonDivisorStructure + FactorableStructure,
> HenselProduct<LIFTED_BEZOUT_COEFFS, RS>
{
    fn check(
        &self,
        ring: &RS,
        h: &Polynomial<RS::Set>,
        i: &RS::Set,
        n: &Natural,
    ) -> Result<(), &'static str> {
        match self {
            HenselProduct::Leaf => {}
            HenselProduct::Branch {
                f_factorization,
                g_factorization,
                a,
                b,
            } => {
                f_factorization.check(ring, i, n)?;
                g_factorization.check(ring, i, n)?;

                let poly_ring = PolynomialStructure::new(ring.clone());
                let ring_mod_i = QuotientStructure::new_ring(ring.clone(), i.clone());
                let poly_ring_mod_i = PolynomialStructure::new(ring_mod_i.clone());
                let poly_ring_mod_i_tothe_n = PolynomialStructure::new(
                    QuotientStructure::new_ring(ring.clone(), ring.nat_pow(i, n)),
                );

                //af + bg = 1 mod i
                if !poly_ring_mod_i.is_zero(&poly_ring_mod_i.sum(vec![
                    poly_ring_mod_i.mul(a, &f_factorization.h),
                    poly_ring_mod_i.mul(b, &g_factorization.h),
                    poly_ring_mod_i.neg(&poly_ring_mod_i.one()),
                ])) {
                    return Err("af + bg != 1 mod i");
                }

                //af + bg = 1 mod i^n
                if LIFTED_BEZOUT_COEFFS {
                    if !poly_ring_mod_i_tothe_n.is_zero(&poly_ring_mod_i_tothe_n.sum(vec![
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
                    &poly_ring_mod_i_tothe_n.product(vec![
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
        p: &RS::Set,
        n: &Natural,
        first_fs: Vec<&Polynomial<RS::Set>>,
        second_fs: Vec<&Polynomial<RS::Set>>,
    ) -> Self {
        let poly_ring = PolynomialStructure::new(ring.clone());
        let poly_ring_mod_p =
            PolynomialStructure::new(QuotientStructure::new_field(ring.clone(), p.clone()));

        //first_h and second_h are defined modulo p^n
        let first_h = poly_ring
            .product(first_fs.clone())
            .apply_map(|c| ring.rem(c, &ring.nat_pow(p, n)));
        let second_h = poly_ring
            .product(second_fs.clone())
            .apply_map(|c| ring.rem(c, &ring.nat_pow(p, n)));

        //find a, b such that af + bg = 1 mod i
        let (u, a, b) = poly_ring_mod_p.xgcd(&first_h, &second_h);
        if !poly_ring_mod_p.equal(&u, &poly_ring_mod_p.one()) {
            panic!("Factors should be coprime modulo i");
        }

        Self::Branch {
            f_factorization: Box::new(HenselFactorizationImpl::new(ring, p, n, first_h, first_fs)),
            g_factorization: Box::new(HenselFactorizationImpl::new(
                ring, p, n, second_h, second_fs,
            )),
            a,
            b,
        }
    }

    fn factor_list<'a>(&'a self, h: &'a Polynomial<RS::Set>) -> Vec<&'a Polynomial<RS::Set>> {
        match self {
            HenselProduct::Leaf => vec![h],
            HenselProduct::Branch {
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

impl<RS: EuclideanDivisionStructure + GreatestCommonDivisorStructure + FactorableStructure>
    HenselProduct<true, RS>
{
    fn dont_lift_bezout_coeffs(self) -> HenselProduct<false, RS> {
        match self {
            HenselProduct::Leaf => HenselProduct::Leaf,
            HenselProduct::Branch {
                f_factorization,
                g_factorization,
                a,
                b,
            } => HenselProduct::Branch {
                f_factorization: Box::new(f_factorization.dont_lift_bezout_coeffs()),
                g_factorization: Box::new(g_factorization.dont_lift_bezout_coeffs()),
                a,
                b,
            },
        }
    }
}

fn compute_lift_factors<
    RS: EuclideanDivisionStructure + GreatestCommonDivisorStructure + FactorableStructure,
>(
    ring: &RS,
    i: &RS::Set,
    n: &Natural,
    a: &Polynomial<RS::Set>,
    b: &Polynomial<RS::Set>,
    f: &Polynomial<RS::Set>,
    g: &Polynomial<RS::Set>,
    h: &Polynomial<RS::Set>,
) -> (
    Polynomial<RS::Set>,
    Polynomial<RS::Set>,
    Polynomial<RS::Set>,
    Polynomial<RS::Set>,
) {
    let poly_ring = PolynomialStructure::new(ring.clone());

    let alpha = poly_ring.leading_coeff(h).unwrap();
    let (gcd, beta, gamma) = ring.euclidean_xgcd(alpha.clone(), i.clone());
    debug_assert!(ring.equal(&gcd, &ring.one()));
    drop(gcd);
    drop(gamma);

    let ring_mod_i = QuotientStructure::new_ring(ring.clone(), i.clone());
    debug_assert!(ring_mod_i.equal(alpha, poly_ring.leading_coeff(h).unwrap()));

    let delta_h = poly_ring
        .add(
            h,
            &poly_ring.neg(&poly_ring.product(vec![&Polynomial::constant(alpha.clone()), f, g])),
        )
        .apply_map(|c| ring.rem(c, &ring.nat_pow(i, &(n + Natural::ONE))));

    //found delta_h such that
    //delta_h = h - alpha*f*g mod i^n+1
    let poly_ring_mod_i_tothe_nplusone = PolynomialStructure::new(QuotientStructure::new_ring(
        ring.clone(),
        ring.nat_pow(i, &(n + Natural::ONE)),
    ));
    debug_assert!(poly_ring_mod_i_tothe_nplusone.equal(
        &delta_h,
        &poly_ring_mod_i_tothe_nplusone.add(
            h,
            &poly_ring_mod_i_tothe_nplusone.neg(&poly_ring_mod_i_tothe_nplusone.product(vec![
                &Polynomial::constant(alpha.clone()),
                f,
                g,
            ]))
        ),
    ));
    let poly_ring = PolynomialStructure::new(ring.clone());

    debug_assert!(poly_ring.degree(&delta_h).unwrap_or(0) < poly_ring.degree(h).unwrap());

    //(qg, rg) = quorem(a * delta_h, g)
    //(qf, rf) = quorem(b * delta_h, f)
    let (qg, rg) = poly_ring
        .try_quorem(&poly_ring.mul(a, &delta_h), g)
        .unwrap();
    let (qf, rf) = poly_ring
        .try_quorem(&poly_ring.mul(b, &delta_h), f)
        .unwrap();

    //qf + qg = 0 mod i^{n+1}
    debug_assert!(
        poly_ring_mod_i_tothe_nplusone.is_zero(&poly_ring_mod_i_tothe_nplusone.add(&qf, &qg))
    );

    let delta_f = poly_ring.mul(&Polynomial::constant(beta.clone()), &rf);
    let delta_g = poly_ring.mul(&Polynomial::constant(beta.clone()), &rg);

    let lifted_f = poly_ring
        .add(&delta_f, f)
        .apply_map(|c| ring.rem(c, &ring.nat_pow(i, &(n + Natural::ONE))));
    let lifted_g = poly_ring
        .add(&delta_g, g)
        .apply_map(|c| ring.rem(c, &ring.nat_pow(i, &(n + Natural::ONE))));

    (delta_f, delta_g, lifted_f, lifted_g)
}

impl<RS: EuclideanDivisionStructure + GreatestCommonDivisorStructure + FactorableStructure>
    HenselProduct<false, RS>
{
    fn linear_lift(&mut self, ring: &RS, i: &RS::Set, n: &Natural, h: &Polynomial<RS::Set>) {
        match self {
            HenselProduct::Leaf => {}
            HenselProduct::Branch {
                f_factorization,
                g_factorization,
                a,
                b,
            } => {
                let f = &f_factorization.h;
                let g = &g_factorization.h;
                let (_, _, lifted_f, lifted_g) = compute_lift_factors(ring, i, n, a, b, &f, &g, h);
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

impl<RS: EuclideanDivisionStructure + GreatestCommonDivisorStructure + FactorableStructure>
    HenselProduct<true, RS>
{
    fn quadratic_lift(&mut self, ring: &RS, i: &RS::Set, n: &Natural, h: &Polynomial<RS::Set>) {
        match self {
            HenselProduct::Leaf => {}
            HenselProduct::Branch {
                f_factorization,
                g_factorization,
                a,
                b,
            } => {
                let pring_mod_i2n = PolynomialStructure::new(QuotientStructure::new_ring(
                    ring.clone(),
                    ring.nat_pow(i, &(n * Natural::TWO)),
                ));

                let f = &f_factorization.h;
                let g = &g_factorization.h;
                let (delta_f, delta_g, lifted_f, lifted_g) =
                    compute_lift_factors(ring, &ring.nat_pow(i, n), &Natural::ONE, a, b, &f, &g, h);

                // beta = af + bg - 1 mod i^n
                let beta = pring_mod_i2n.sum(vec![
                    pring_mod_i2n.mul(a, f),
                    pring_mod_i2n.mul(b, g),
                    pring_mod_i2n.neg(&pring_mod_i2n.one()),
                ]);

                // big_delta = beta + a * delta_f + b * delta_g mod i^n
                let big_delta = pring_mod_i2n.sum(vec![
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
    RS: EuclideanDivisionStructure + GreatestCommonDivisorStructure + FactorableStructure,
> HenselFactorizationImpl<LIFTED_BEZOUT_COEFFS, RS>
{
    fn check(&self, ring: &RS, i: &RS::Set, n: &Natural) -> Result<(), &'static str> {
        // let poly_ring = PolynomialStructure::new(ring.clone().into());
        // if !poly_ring.is_monic(&self.h) {
        //     return Err("h is not monic");
        // }
        self.factorization.check(ring, &self.h, i, n)?;
        Ok(())
    }

    fn new(
        ring: &RS,
        p: &RS::Set,
        n: &Natural,
        h: Polynomial<RS::Set>,
        mut fs: Vec<&Polynomial<RS::Set>>,
    ) -> Self {
        debug_assert!(fs.len() >= 1);
        match fs.len() {
            0 => panic!(),
            1 => Self {
                h,
                factorization: HenselProduct::Leaf,
            },
            fs_len => {
                debug_assert!(fs_len >= 2);
                let second_fs = fs.split_off(fs_len / 2);
                let first_fs = fs;
                debug_assert!(first_fs.len() >= 1);
                debug_assert!(second_fs.len() >= 1);
                debug_assert_eq!(first_fs.len() + second_fs.len(), fs_len);

                //find an inverse beta to alpha modulo p
                let alpha = PolynomialStructure::new(ring.clone())
                    .leading_coeff(&h)
                    .unwrap();
                let (g, _beta, _gamma) = ring.euclidean_xgcd(alpha.clone(), p.clone());
                debug_assert!(ring.equal(&g, &ring.one()));

                Self {
                    h,
                    factorization: HenselProduct::new_split(ring, p, n, first_fs, second_fs),
                }
            }
        }
    }

    fn factor_list(&self) -> Vec<&Polynomial<RS::Set>> {
        self.factorization.factor_list(&self.h)
    }
}

impl<RS: EuclideanDivisionStructure + GreatestCommonDivisorStructure + FactorableStructure>
    HenselFactorizationImpl<true, RS>
{
    fn dont_lift_bezout_coeffs(self) -> HenselFactorizationImpl<false, RS> {
        HenselFactorizationImpl {
            h: self.h,
            factorization: self.factorization.dont_lift_bezout_coeffs(),
        }
    }
}

impl<RS: EuclideanDivisionStructure + GreatestCommonDivisorStructure + FactorableStructure>
    HenselFactorizationImpl<false, RS>
{
    fn linear_lift(&mut self, ring: &RS, i: &RS::Set, n: &Natural) {
        self.factorization.linear_lift(ring, i, n, &self.h);
    }
}

impl<RS: EuclideanDivisionStructure + GreatestCommonDivisorStructure + FactorableStructure>
    HenselFactorizationImpl<true, RS>
{
    fn quadratic_lift(&mut self, ring: &RS, i: &RS::Set, n: &Natural) {
        self.factorization.quadratic_lift(ring, i, n, &self.h);
    }
}

impl<
    const LIFTED_BEZOUT_COEFFS: bool,
    RS: EuclideanDivisionStructure + GreatestCommonDivisorStructure + FactorableStructure,
> HenselFactorization<LIFTED_BEZOUT_COEFFS, RS>
{
    fn check(&self) -> Result<(), &'static str> {
        self.factorization.check(&self.ring, &self.i, &self.n)
    }

    pub fn new(
        ring: RS,
        p: RS::Set,
        n: Natural,
        h: Polynomial<RS::Set>,
        fs: Vec<Polynomial<RS::Set>>,
    ) -> Self {
        let poly_ring = PolynomialStructure::new(ring.clone());

        debug_assert!(ring.is_irreducible(&p));

        // h and all fs monic
        assert!(fs.len() >= 1);
        // debug_assert!(poly_ring.is_monic(&h));
        for f in fs.iter() {
            debug_assert!(poly_ring.is_monic(f));
        }
        // h = product of fs modulo i^n
        let poly_ring_mod_p_tothe_n = PolynomialStructure::new(QuotientStructure::new_ring(
            ring.clone(),
            ring.nat_pow(&p, &n),
        ));
        let alpha = poly_ring.leading_coeff(&h).unwrap();
        debug_assert!(poly_ring_mod_p_tothe_n.equal(
            &h,
            &poly_ring_mod_p_tothe_n.mul(
                &Polynomial::constant(alpha.clone()),
                &poly_ring_mod_p_tothe_n.product(fs.iter().collect())
            )
        ));
        // fs are coprime mod i - checked when computing bezout coefficients
        let ans_ring = ring.clone();
        let factors = HenselFactorizationImpl::new(&ring, &p, &n, h, fs.iter().collect());
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

    pub fn modolus(&self) -> RS::Set {
        self.ring.nat_pow(&self.i, &self.n)
    }

    //return the lifted factors in order
    pub fn factors(&self) -> Vec<&Polynomial<RS::Set>> {
        self.factorization.factor_list()
    }
}

impl<RS: EuclideanDivisionStructure + GreatestCommonDivisorStructure + FactorableStructure>
    HenselFactorization<true, RS>
{
    pub fn dont_lift_bezout_coeffs(self) -> HenselFactorization<false, RS> {
        HenselFactorization {
            ring: self.ring,
            i: self.i,
            n: self.n,
            factorization: self.factorization.dont_lift_bezout_coeffs(),
        }
    }
}

impl<RS: EuclideanDivisionStructure + GreatestCommonDivisorStructure + FactorableStructure>
    HenselFactorization<false, RS>
{
    pub fn linear_lift(&mut self) {
        self.factorization.linear_lift(&self.ring, &self.i, &self.n);
        self.n += Natural::ONE;
    }
}

impl<RS: EuclideanDivisionStructure + GreatestCommonDivisorStructure + FactorableStructure>
    HenselFactorization<true, RS>
{
    pub fn quadratic_lift(&mut self) {
        self.factorization
            .quadratic_lift(&self.ring, &self.i, &self.n);
        self.n *= Natural::TWO;
    }
}

impl<RS: FactorableStructure + EuclideanDivisionStructure + GreatestCommonDivisorStructure>
    Factored<PolynomialStructure<QuotientStructure<RS, true>>>
where
    PolynomialStructure<QuotientStructure<RS, true>>:
        SetStructure<Set = Polynomial<RS::Set>> + FactorableStructure,
{
    /// If the polynomial is squarefree return a hensel factorization, otherwise return None
    pub fn into_hensel_factorization(
        self,
        h: Polynomial<RS::Set>,
    ) -> Option<HenselFactorization<true, RS>>
    where
        RS: EuclideanDivisionStructure + GreatestCommonDivisorStructure,
    {
        let poly_ring_mod = self.ring().clone();
        let ring_mod = poly_ring_mod.coeff_ring();
        let ring = ring_mod.ring();
        // let poly_ring: PolynomialStructure<RS> = PolynomialStructure::new(ring.clone());

        let mut fs = vec![];
        let (_unit, factors) = self.unit_and_factors();
        for (factor, power) in factors {
            if power == Natural::ONE {
                fs.push(factor);
            } else {
                return None;
            }
        }

        let hensel_factorization =
            HenselFactorization::new(ring, ring_mod.modulus().clone(), Natural::ONE, h, fs);
        #[cfg(debug_assertions)]
        hensel_factorization.check().unwrap();
        Some(hensel_factorization)
    }
}

#[cfg(test)]
mod tests {
    use algebraeon_nzq::*;

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
        let h = Polynomial::product(vec![&Polynomial::constant(Integer::from(8)), &f1, &f2, &f3]);
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
                &Polynomial::constant(Polynomial::leading_coeff(&h).unwrap()),
                &Polynomial::product(hensel_fact.factors()),
            )
            .apply_map(|c| Integer::rem(c, &hensel_fact.modolus()));
            println!("{:?} {:?}", lifted_product, h);
            assert_eq!(lifted_product, h);
        }
    }
}
