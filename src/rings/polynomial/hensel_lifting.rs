use std::borrow::Borrow;
use std::rc::Rc;

use malachite_base::num::basic::traits::One;
use malachite_nz::natural::Natural;

use super::super::super::structure::*;
use super::super::ring_structure::quotient::*;
use super::super::ring_structure::structure::*;
use super::polynomial::*;

// impl<RS: EuclideanDivisionStructure + GreatestCommonDivisorStructure + UniqueFactorizationStructure> PolynomialStructure<RS> {
//     fn polynomial_reduce_modulo(&self,
//         poly: &Polynomial<RS::Set>,
//         i: impl Borrow<RS::Set>,
//         n: impl Borrow<Natural>,
//     ) -> Polynomial<EuclideanQuotient<RS>> {
//         let i = i.borrow();
//         let n = n.borrow();
//         poly.apply_map_ref(|c| UniversalEuclideanQuotient::<false, Ring>::new(c.clone(), i.nat_pow(n)))
//     }
// }

#[derive(Debug, Clone)]
enum HenselProduct<
    RS: EuclideanDivisionStructure + GreatestCommonDivisorStructure + UniqueFactorizationStructure,
> {
    Leaf,
    Branch {
        //h = fg mod i^n
        //af + bg = 1 mod i
        f_factorization: Box<HenselFactorizationImpl<RS>>, //defined modulo i^n
        g_factorization: Box<HenselFactorizationImpl<RS>>, //defined modulo i^n
        a: Polynomial<RS::Set>,                            //defined modulo i
        b: Polynomial<RS::Set>,                            //defined modulo i
    },
}

#[derive(Debug, Clone)]
struct HenselFactorizationImpl<
    RS: EuclideanDivisionStructure + GreatestCommonDivisorStructure + UniqueFactorizationStructure,
> {
    h: Polynomial<RS::Set>,
    factorization: HenselProduct<RS>,
}

//represent a factorization of a monic polynomial h(x) into coprime monic polynomials f_1(x)f_2(x)...f_k(x) modulo i^n
#[derive(Debug, Clone)]
pub struct HenselFactorization<
    RS: EuclideanDivisionStructure + GreatestCommonDivisorStructure + UniqueFactorizationStructure,
> {
    ring: Rc<RS>,
    i: RS::Set,
    n: Natural,
    factors_impl: HenselFactorizationImpl<RS>, //defined absolutely and factored modulo i^n
}

impl<
        RS: EuclideanDivisionStructure + GreatestCommonDivisorStructure + UniqueFactorizationStructure,
    > HenselProduct<RS>
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

                //af + bg = 1 mod i
                let poly_ring_mod_i = PolynomialStructure::new(
                    QuotientStructure::new_ring(ring.clone().into(), i.clone()).into(),
                );
                if !poly_ring_mod_i.is_zero(&poly_ring_mod_i.sum(vec![
                    poly_ring_mod_i.mul(a, &f_factorization.h),
                    poly_ring_mod_i.mul(b, &g_factorization.h),
                    poly_ring_mod_i.neg(&poly_ring_mod_i.one()),
                ])) {
                    return Err("af + bg != 1 mod i");
                }
                //h = fg mod i^n
                let poly_ring_mod_i_tothe_n = PolynomialStructure::new(
                    QuotientStructure::new_ring(ring.clone().into(), ring.nat_pow(i, n)).into(),
                );
                if !poly_ring_mod_i_tothe_n.is_zero(
                    &poly_ring_mod_i_tothe_n.add(
                        &poly_ring_mod_i_tothe_n.neg(
                            &poly_ring_mod_i_tothe_n.mul(&f_factorization.h, &g_factorization.h),
                        ),
                        h,
                    ),
                ) {
                    return Err("h != fg mod i^n");
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
        let poly_ring = PolynomialStructure::new(ring.clone().into());
        let poly_ring_mod_p = PolynomialStructure::new(
            QuotientStructure::new_field(ring.clone().into(), p.clone()).into(),
        );

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
                a,
                b,
            } => {
                let mut fs = vec![];
                fs.extend(f_factorization.factor_list());
                fs.extend(g_factorization.factor_list());
                fs
            }
        }
    }

    fn linear_lift(&mut self, ring: &RS, i: &RS::Set, n: &Natural, h: &Polynomial<RS::Set>) {
        let poly_ring = PolynomialStructure::new(ring.clone().into());

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

                // println!("{:?}", self);
                //  h = fg mod i^n  =>  h - fg = 0 mod i^n
                let h_minus_fg = poly_ring.add(&poly_ring.neg(&poly_ring.mul(f, g)), h);
                let delta_h_over_i_tothe_n = h_minus_fg
                    .apply_map(|c| ring.rem(&ring.quo(c, &ring.nat_pow(i, n)).unwrap(), i));
                // println!("delta_h_over_i_tothe_n = {:?}", &delta_h_over_i_tothe_n);

                let delta_h = poly_ring.mul(
                    &delta_h_over_i_tothe_n,
                    &Polynomial::constant(ring.nat_pow(i, n)),
                );

                //found delta_h such that
                //delta_h = h - fg mod i^n+1
                let poly_ring_mod_i_tothe_nplusone = PolynomialStructure::new(
                    QuotientStructure::new_ring(
                        ring.clone().into(),
                        ring.nat_pow(i, &(n + Natural::ONE)),
                    )
                    .into(),
                );

                debug_assert!(poly_ring_mod_i_tothe_nplusone.equal(
                    &poly_ring_mod_i_tothe_nplusone.add(
                        &poly_ring_mod_i_tothe_nplusone
                            .neg(&poly_ring_mod_i_tothe_nplusone.mul(f, g)),
                        h
                    ),
                    &delta_h
                ));

                // println!("delta_h = {:?}", &delta_h);

                //(qg, rg) = quorem(a * delta_h, g)
                //(qf, rf) = quorem(b * delta_h, f)
                let (qg, rg) = poly_ring
                    .try_quorem(&poly_ring.mul(a, &delta_h), g)
                    .unwrap();
                let (qf, rf) = poly_ring
                    .try_quorem(&poly_ring.mul(b, &delta_h), f)
                    .unwrap();

                // println!("qg = {:?}  rg = {:?}", qg, rg);
                // println!("qf = {:?}  rf = {:?}", qf, rf);

                //qf + qg = 0 mod i^{n+1}
                debug_assert!(poly_ring_mod_i_tothe_nplusone
                    .is_zero(&poly_ring_mod_i_tothe_nplusone.add(&qf, &qg),));

                let lifted_f = poly_ring
                    .add(&rf, f)
                    .apply_map(|c| ring.rem(c, &ring.nat_pow(i, &(n + Natural::ONE))));
                let lifted_g = poly_ring
                    .add(&rg, g)
                    .apply_map(|c| ring.rem(c, &ring.nat_pow(i, &(n + Natural::ONE))));

                f_factorization.h = lifted_f;
                g_factorization.h = lifted_g;

                f_factorization.linear_lift(ring, i, n);
                g_factorization.linear_lift(ring, i, n);
            }
        }
    }
}

impl<
        RS: EuclideanDivisionStructure + GreatestCommonDivisorStructure + UniqueFactorizationStructure,
    > HenselFactorizationImpl<RS>
{
    fn check(&self, ring: &RS, i: &RS::Set, n: &Natural) -> Result<(), &'static str> {
        let poly_ring = PolynomialStructure::new(ring.clone().into());

        if !poly_ring.is_monic(&self.h) {
            return Err("h is not monic");
        }
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
                // println!("{:?}, {:?}", first_fs, second_fs);
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

    fn linear_lift(&mut self, ring: &RS, i: &RS::Set, n: &Natural) {
        self.factorization.linear_lift(ring, i, n, &self.h);
    }
}

impl<
        RS: EuclideanDivisionStructure + GreatestCommonDivisorStructure + UniqueFactorizationStructure,
    > HenselFactorization<RS>
{
    fn check(&self) -> Result<(), &'static str> {
        self.factors_impl
            .check(self.ring.as_ref(), &self.i, &self.n)
    }

    pub fn new(
        ring: Rc<RS>,
        p: RS::Set,
        n: Natural,
        h: Polynomial<RS::Set>,
        fs: Vec<Polynomial<RS::Set>>,
    ) -> Self {
        let poly_ring = PolynomialStructure::new(ring.clone());

        debug_assert!(ring.is_irreducible(&p));

        // h and all fs monic
        assert!(fs.len() >= 1);
        debug_assert!(poly_ring.is_monic(&h));
        for f in fs.iter() {
            debug_assert!(poly_ring.is_monic(f));
        }
        // h = product of fs modulo i^n
        let poly_ring_mod_p_tothe_n = PolynomialStructure::new(
            QuotientStructure::new_ring(ring.clone(), ring.nat_pow(&p, &n)).into(),
        );
        debug_assert!(poly_ring_mod_p_tothe_n
            .equal(&h, &poly_ring_mod_p_tothe_n.product(fs.iter().collect())));
        // fs are coprime mod i - checked when computing bezout coefficients
        let factors = HenselFactorizationImpl::new(ring.as_ref(), &p, &n, h, fs.iter().collect());
        let ans = Self {
            ring,
            i: p,
            n,
            factors_impl: factors,
        };
        debug_assert!(ans.check().is_ok());
        ans
    }

    pub fn modolus(&self) -> RS::Set {
        self.ring.nat_pow(&self.i, &self.n)
    }

    //return the lifted factors in order
    pub fn factors(&self) -> Vec<&Polynomial<RS::Set>> {
        self.factors_impl.factor_list()
    }

    pub fn linear_lift(&mut self) {
        self.factors_impl
            .linear_lift(self.ring.as_ref(), &self.i, &self.n);
        self.n += Natural::ONE;
    }
}

#[cfg(test)]
mod tests {
    use malachite_nz::integer::Integer;

    use super::super::super::ring_structure::cannonical::*;

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
        let h = Polynomial::product(vec![&f1, &f2, &f3]);
        let h = h.apply_map(|c| Integer::rem(c, &Integer::from(5)));
        //h = prod fs mod 5
        //fs are coprime
        //h and fs are monic

        //set up bezout coefficients for hensel lifting the factorization modulo 25, 125, ...
        let mut hensel_fact = HenselFactorization::new(
            Integer::structure(),
            Integer::from(5),
            Natural::from(1u8),
            h.clone(),
            vec![f1, f2, f3],
        );
        hensel_fact.check().unwrap();
        println!("5^1: {:?}", hensel_fact.factors());
        for i in 2..20 {
            hensel_fact.linear_lift();
            hensel_fact.check().unwrap();
            println!("5^{}: {:?}", i, hensel_fact.factors());
            let lifted_product = Polynomial::product(hensel_fact.factors())
                .apply_map(|c| Integer::rem(c, &hensel_fact.modolus()));
            assert_eq!(lifted_product, h);
        }
    }
}
