use std::borrow::Borrow;

use malachite_base::num::basic::traits::One;
use malachite_nz::natural::Natural;

use super::poly::*;
use super::ring::*;

fn polynomial_reduce_modulo<Ring: EuclideanDomain + UniqueFactorizationDomain>(
    poly: &Polynomial<Ring>,
    i: impl Borrow<Ring>,
    n: impl Borrow<Natural>,
) -> Polynomial<EuclideanQuotient<false, Ring>> {
    let i = i.borrow();
    let n = n.borrow();
    poly.apply_map(|c| EuclideanQuotient::<false, Ring>::new(c.clone(), i.nat_pow(n)))
}

#[derive(Debug, Clone)]
enum HenselProduct<Ring: EuclideanDomain + UniqueFactorizationDomain> {
    Leaf,
    Branch {
        //h = fg mod i^n
        //af + bg = 1 mod i
        f_factorization: Box<HenselFactorizationImpl<Ring>>, //defined modulo i^n
        g_factorization: Box<HenselFactorizationImpl<Ring>>, //defined modulo i^n
        a: Polynomial<Ring>,                                 //defined modulo i
        b: Polynomial<Ring>,                                 //defined modulo i
    },
}

impl<Ring: EuclideanDomain + UniqueFactorizationDomain> HenselProduct<Ring> {
    fn check(&self, h: &Polynomial<Ring>, i: &Ring, n: &Natural) -> Result<(), &'static str> {
        match self {
            HenselProduct::Leaf => {}
            HenselProduct::Branch {
                f_factorization,
                g_factorization,
                a,
                b,
            } => {
                f_factorization.check(i, n)?;
                g_factorization.check(i, n)?;
                //af + bg = 1 mod i
                if polynomial_reduce_modulo(
                    &Polynomial::sum(vec![
                        &Polynomial::mul_refs(a, &f_factorization.h),
                        &Polynomial::mul_refs(b, &g_factorization.h),
                        &Polynomial::one().neg(),
                    ]),
                    i,
                    &Natural::from(1u8),
                ) != Polynomial::zero()
                {
                    return Err("af + bg != 1 mod i");
                }
                //h = fg mod i^n
                if polynomial_reduce_modulo(
                    &Polynomial::add_ref(
                        Polynomial::mul_refs(&f_factorization.h, &g_factorization.h).neg(),
                        h,
                    ),
                    i,
                    n,
                ) != Polynomial::zero()
                {
                    return Err("h != fg mod i^n");
                }
            }
        }
        Ok(())
    }

    fn new_split(
        p: &Ring,
        n: &Natural,
        first_fs: Vec<&Polynomial<Ring>>,
        second_fs: Vec<&Polynomial<Ring>>,
    ) -> Self {
        let first_h = Polynomial::product(first_fs.clone())
            .apply_map(|c| Ring::rem_lref(c, p.nat_pow(n)).unwrap());
        let second_h = Polynomial::product(second_fs.clone())
            .apply_map(|c| Ring::rem_lref(c, p.nat_pow(n)).unwrap());

        let (u, a, b) = Polynomial::xgcd(
            first_h.apply_map(|c| EuclideanQuotient::<true, Ring>::new(c.clone(), p.clone())),
            second_h.apply_map(|c| EuclideanQuotient::<true, Ring>::new(c.clone(), p.clone())),
        );
        if u != Polynomial::one() {
            panic!("Factors should be coprime modulo i");
        }
        let a = a.apply_map(|c| c.clone().lift());
        let b = b.apply_map(|c| c.clone().lift());

        Self::Branch {
            f_factorization: Box::new(HenselFactorizationImpl::new(p, n, first_h, first_fs)),
            g_factorization: Box::new(HenselFactorizationImpl::new(p, n, second_h, second_fs)),
            a,
            b,
        }
    }

    fn factor_list<'a>(&'a self, h: &'a Polynomial<Ring>) -> Vec<&'a Polynomial<Ring>> {
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

    fn linear_lift(&mut self, i: &Ring, n: &Natural, h: &Polynomial<Ring>) {
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
                let h_minus_fg = Polynomial::add_ref(Polynomial::mul_refs(f, g).neg(), h);
                let delta_h_over_i_tothe_n = h_minus_fg.apply_map(|c| {
                    Ring::rem_rref(Ring::quo_lref(c, i.nat_pow(n)).unwrap(), i).unwrap()
                });
                // println!("delta_h_over_i_tothe_n = {:?}", &delta_h_over_i_tothe_n);

                let delta_h =
                    Polynomial::mul(delta_h_over_i_tothe_n, Polynomial::constant(i.nat_pow(n)));

                //found delta_h such that
                //delta_h = h - fg mod i^n+1
                debug_assert_eq!(
                    polynomial_reduce_modulo(
                        &Polynomial::add_ref(Polynomial::mul_refs(f, g).neg(), h),
                        i,
                        n + Natural::from(1u8)
                    ),
                    polynomial_reduce_modulo(&delta_h, i, n + Natural::from(1u8))
                );

                // println!("delta_h = {:?}", &delta_h);

                //(qg, rg) = quorem(a * delta_h, g)
                //(qf, rf) = quorem(b * delta_h, f)
                let (qg, rg) =
                    Polynomial::try_quorem_rref(Polynomial::mul_refs(a, &delta_h), g).unwrap();
                let (qf, rf) =
                    Polynomial::try_quorem_rref(Polynomial::mul_refs(b, &delta_h), f).unwrap();

                // println!("qg = {:?}  rg = {:?}", qg, rg);
                // println!("qf = {:?}  rf = {:?}", qf, rf);

                //qf + qg = 0 mod i^{n+1}
                debug_assert_eq!(
                    polynomial_reduce_modulo(
                        &Polynomial::add_refs(&qf, &qg),
                        i,
                        n + Natural::from(1u8)
                    ),
                    Polynomial::zero()
                );

                let lifted_f = Polynomial::add_ref(rf, f)
                    .apply_map(|c| Ring::rem_lref(c, i.nat_pow(n + Natural::from(1u8))).unwrap());
                let lifted_g = Polynomial::add_ref(rg, g)
                    .apply_map(|c| Ring::rem_lref(c, i.nat_pow(n + Natural::from(1u8))).unwrap());

                f_factorization.h = lifted_f;
                g_factorization.h = lifted_g;

                f_factorization.linear_lift(i, n);
                g_factorization.linear_lift(i, n);
            }
        }
    }
}

#[derive(Debug, Clone)]
struct HenselFactorizationImpl<Ring: EuclideanDomain + UniqueFactorizationDomain> {
    h: Polynomial<Ring>,
    factorization: HenselProduct<Ring>,
}

impl<Ring: EuclideanDomain + UniqueFactorizationDomain> HenselFactorizationImpl<Ring> {
    fn check(&self, i: &Ring, n: &Natural) -> Result<(), &'static str> {
        if !self.h.is_monic() {
            return Err("h is not monic");
        }
        self.factorization.check(&self.h, i, n)?;
        Ok(())
    }

    fn new(p: &Ring, n: &Natural, h: Polynomial<Ring>, mut fs: Vec<&Polynomial<Ring>>) -> Self {
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
                    factorization: HenselProduct::new_split(p, n, first_fs, second_fs),
                }
            }
        }
    }

    fn factor_list(&self) -> Vec<&Polynomial<Ring>> {
        self.factorization.factor_list(&self.h)
    }

    fn linear_lift(&mut self, i: &Ring, n: &Natural) {
        self.factorization.linear_lift(i, n, &self.h);
    }
}

//represent a factorization of a monic polynomial h(x) into coprime monic polynomials f_1(x)f_2(x)...f_k(x) modulo i^n
#[derive(Debug, Clone)]
pub struct HenselFactorization<Ring: EuclideanDomain + UniqueFactorizationDomain> {
    i: Ring,
    n: Natural,
    factors_impl: HenselFactorizationImpl<Ring>, //defined absolutely and factored modulo i^n
}

impl<Ring: EuclideanDomain + UniqueFactorizationDomain> HenselFactorization<Ring> {
    fn check(&self) -> Result<(), &'static str> {
        self.factors_impl.check(&self.i, &self.n)
    }

    pub fn new(p: Ring, n: Natural, h: Polynomial<Ring>, fs: Vec<Polynomial<Ring>>) -> Self {
        debug_assert!(p.is_irreducible().unwrap());

        // h and all fs monic
        assert!(fs.len() >= 1);
        debug_assert!(h.is_monic());
        for f in fs.iter() {
            debug_assert!(f.is_monic());
        }
        // h = product of fs modulo i^n
        debug_assert_eq!(
            polynomial_reduce_modulo(&h, &p, &n),
            polynomial_reduce_modulo(&Polynomial::product(fs.iter().collect()), &p, &n)
        );
        // fs are coprime mod i - checked when computing bezout coefficients
        let factors = HenselFactorizationImpl::new(&p, &n, h, fs.iter().collect());
        let ans = Self {
            i: p,
            n,
            factors_impl: factors,
        };
        debug_assert!(ans.check().is_ok());
        ans
    }

    pub fn modolus(&self) -> Ring {
        self.i.nat_pow(&self.n)
    }

    //return the lifted factors in order
    pub fn factors(&self) -> Vec<&Polynomial<Ring>> {
        self.factors_impl.factor_list()
    }

    pub fn linear_lift(&mut self) {
        self.factors_impl.linear_lift(&self.i, &self.n);
        self.n += Natural::ONE;
    }
}

#[cfg(test)]
mod tests {
    use malachite_nz::integer::Integer;

    use crate::Modulo;

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
        let h = h.apply_map(|c| Integer::rem_lref(c, Integer::from(5)).unwrap());
        //h = prod fs mod 5
        //fs are coprime
        //h and fs are monic

        //set up bezout coefficients for hensel lifting the factorization modulo 25, 125, ...
        let mut hensel_fact = HenselFactorization::new(
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
                .apply_map(|c| Integer::rem_lref(c, hensel_fact.modolus()).unwrap());
            assert_eq!(lifted_product, h);
        }
    }
}
