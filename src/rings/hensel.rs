use malachite_nz::natural::Natural;

use super::poly::*;
use super::ring::*;

fn polynomial_reduce_modulo<Ring: EuclideanDomain + UniqueFactorizationDomain>(
    poly: &Polynomial<Ring>,
    i: &Ring,
    n: &Natural,
) -> Polynomial<EuclideanQuotient<false, Ring>> {
    poly.apply_map(|c| EuclideanQuotient::<false, Ring>::new(c.clone(), i.nat_pow(n)))
}

#[derive(Debug, Clone)]
struct HenselPairFactorizationOld<Ring: EuclideanDomain + UniqueFactorizationDomain> {
    //represent a factorization of h(x) modulo i^n as f(x)g(x) with f and g coprime
    /*
        h, f, g monic
        af + bg = 1 mod i
        h = fg mod i^n
    */
    i: Ring,
    n: Natural,
    h: Polynomial<Ring>,
    f: Polynomial<Ring>,
    g: Polynomial<Ring>,
    a: Polynomial<Ring>,
    b: Polynomial<Ring>,
}

impl<Ring: EuclideanDomain + UniqueFactorizationDomain> HenselPairFactorizationOld<Ring> {
    fn new(
        i: Ring,
        h: Polynomial<Ring>,
        f: Polynomial<Ring>,
        g: Polynomial<Ring>,
    ) -> Result<Self, &'static str> {
        if polynomial_reduce_modulo(
            &Polynomial::add_ref(Polynomial::mul_refs(&f, &g).neg(), &h),
            &i,
            &Natural::from(1u8),
        ) != Polynomial::zero()
        {
            return Err("h != fg mod i");
        }

        let f_mod = f.apply_map(|c| EuclideanQuotient::<true, Ring>::new(c.clone(), i.clone()));
        let g_mod = g.apply_map(|c| EuclideanQuotient::<true, Ring>::new(c.clone(), i.clone()));
        let (u, a, b) = Polynomial::xgcd(f_mod, g_mod);
        if u != Polynomial::one() {
            return Err("f and g should be coprime modulo i");
        }
        let a = a.apply_map(|c| c.clone().lift());
        let b = b.apply_map(|c| c.clone().lift());

        let hensel_factorization = Self {
            i,
            n: Natural::from(1u8),
            h,
            f,
            g,
            a,
            b,
        };
        debug_assert!(hensel_factorization.check().is_ok());
        Ok(hensel_factorization)
    }
}

impl<Ring: EuclideanDomain + UniqueFactorizationDomain> HenselPairFactorizationOld<Ring> {
    fn check(&self) -> Result<(), &'static str> {
        if self.n == 0 {
            return Err("n == 0");
        }
        if !self.h.is_monic() {
            return Err("h should be monic");
        }
        if !self.f.is_monic() {
            return Err("f should be monic");
        }
        if !self.g.is_monic() {
            return Err("g should be monic");
        }
        //af + bg = 1 mod i
        if polynomial_reduce_modulo(
            &Polynomial::sum(vec![
                &Polynomial::mul_refs(&self.a, &self.f),
                &Polynomial::mul_refs(&self.b, &self.g),
                &Polynomial::one().neg(),
            ]),
            &self.i,
            &Natural::from(1u8),
        ) != Polynomial::zero()
        {
            return Err("af + bg != 1 mod i");
        }
        //h = fg mod i^n
        if polynomial_reduce_modulo(
            &Polynomial::add_ref(Polynomial::mul_refs(&self.f, &self.g).neg(), &self.h),
            &self.i,
            &self.n,
        ) != Polynomial::zero()
        {
            return Err("h != fg mod i^n");
        }
        Ok(())
    }

    fn linear_lift(&self) -> Self {
        // println!("{:?}", self);
        //  h = fg mod i^n  =>  h - fg = 0 mod i^n
        let h_minus_fg = Polynomial::add_ref(Polynomial::mul_refs(&self.f, &self.g).neg(), &self.h);
        let delta_h_over_i_tothe_n = h_minus_fg.apply_map(|c| {
            Ring::rem_rref(Ring::quo_lref(c, self.i.nat_pow(&self.n)).unwrap(), &self.i).unwrap()
        });
        // println!("delta_h_over_i_tothe_n = {:?}", &delta_h_over_i_tothe_n);

        let delta_h = Polynomial::mul(
            delta_h_over_i_tothe_n,
            Polynomial::constant(self.i.nat_pow(&self.n)),
        );

        //found delta_h such that
        //delta_h = h - fg mod i^n+1
        debug_assert_eq!(
            polynomial_reduce_modulo(
                &Polynomial::add_ref(Polynomial::mul_refs(&self.f, &self.g).neg(), &self.h),
                &self.i,
                &(&self.n + Natural::from(1u8))
            ),
            polynomial_reduce_modulo(&delta_h, &self.i, &(&self.n + Natural::from(1u8)))
        );

        // println!("delta_h = {:?}", &delta_h);

        //(qg, rg) = quorem(a * delta_h, g)
        //(qf, rf) = quorem(b * delta_h, f)
        let (qg, rg) =
            Polynomial::try_quorem_rref(Polynomial::mul_refs(&self.a, &delta_h), &self.g).unwrap();
        let (qf, rf) =
            Polynomial::try_quorem_rref(Polynomial::mul_refs(&self.b, &delta_h), &self.f).unwrap();

        // println!("qg = {:?}  rg = {:?}", qg, rg);
        // println!("qf = {:?}  rf = {:?}", qf, rf);

        //qf + qg = 0 mod i^{n+1}
        debug_assert_eq!(
            polynomial_reduce_modulo(
                &Polynomial::add_refs(&qf, &qg),
                &self.i,
                &(&self.n + Natural::from(1u8))
            ),
            Polynomial::zero()
        );

        let hf_lifted = Self {
            i: self.i.clone(),
            n: &self.n + Natural::from(1u8),
            h: self.h.clone(),
            f: Polynomial::add_ref(rf, &self.f),
            g: Polynomial::add_ref(rg, &self.g),
            a: self.a.clone(),
            b: self.b.clone(),
        };
        debug_assert!(hf_lifted.check().is_ok());
        hf_lifted
    }
}

#[derive(Debug, Clone)]
enum HenselProduct<Ring: EuclideanDomain + UniqueFactorizationDomain> {
    Irreducible,
    Split {
        //h = fg mod i^n
        //af + bg = 1 mod i
        f: Box<HenselFactorizationImpl<Ring>>, //defined modulo i^n
        g: Box<HenselFactorizationImpl<Ring>>, //defined modulo i^n
        a: Polynomial<Ring>,
        b: Polynomial<Ring>,
    },
}

impl<Ring: EuclideanDomain + UniqueFactorizationDomain> HenselProduct<Ring> {
    fn check(&self, h: &Polynomial<Ring>, i: &Ring, n: &Natural) -> Result<(), &'static str> {
        match self {
            HenselProduct::Irreducible => {}
            HenselProduct::Split { f, g, a, b } => {
                f.check(i, n)?;
                g.check(i, n)?;
                //af + bg = 1 mod i
                if polynomial_reduce_modulo(
                    &Polynomial::sum(vec![
                        &Polynomial::mul_refs(a, &f.h),
                        &Polynomial::mul_refs(b, &g.h),
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
                    &Polynomial::add_ref(Polynomial::mul_refs(&f.h, &g.h).neg(), h),
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
        i: &Ring,
        n: &Natural,
        first_fs: Vec<&Polynomial<Ring>>,
        second_fs: Vec<&Polynomial<Ring>>,
    ) -> Self {
        let first_h = Polynomial::product(first_fs.clone())
            .apply_map(|c| Ring::rem_lref(c, i.nat_pow(n)).unwrap());
        let second_h = Polynomial::product(second_fs.clone())
            .apply_map(|c| Ring::rem_lref(c, i.nat_pow(n)).unwrap());

        let (u, a, b) = Polynomial::xgcd(
            first_h.apply_map(|c| EuclideanQuotient::<true, Ring>::new(c.clone(), i.clone())),
            second_h.apply_map(|c| EuclideanQuotient::<true, Ring>::new(c.clone(), i.clone())),
        );
        if u != Polynomial::one() {
            panic!("Factors should be coprime modulo i");
        }
        let a = a.apply_map(|c| c.clone().lift());
        let b = b.apply_map(|c| c.clone().lift());

        Self::Split {
            f: Box::new(HenselFactorizationImpl::new(i, n, first_h, first_fs)),
            g: Box::new(HenselFactorizationImpl::new(i, n, second_h, second_fs)),
            a,
            b,
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

    fn new(i: &Ring, n: &Natural, h: Polynomial<Ring>, mut fs: Vec<&Polynomial<Ring>>) -> Self {
        debug_assert!(fs.len() >= 1);
        match fs.len() {
            0 => panic!(),
            1 => Self {
                h,
                factorization: HenselProduct::Irreducible,
            },
            fs_len => {
                debug_assert!(fs_len >= 2);
                let second_fs = fs.split_off(fs_len / 2);
                let first_fs = fs;
                debug_assert!(first_fs.len() >= 1);
                debug_assert!(second_fs.len() >= 1);
                debug_assert_eq!(first_fs.len() + second_fs.len(), fs_len);
                println!("{:?}, {:?}", first_fs, second_fs);
                Self {
                    h,
                    factorization: HenselProduct::new_split(i, n, first_fs, second_fs),
                }
            }
        }
    }
}

//represent a factorization of a monic polynomial h(x) into coprime monic polynomials f_1(x)f_2(x)...f_k(x) modulo i^n
#[derive(Debug, Clone)]
pub struct HenselFactorization<Ring: EuclideanDomain + UniqueFactorizationDomain> {
    i: Ring,
    n: Natural,
    factors: HenselFactorizationImpl<Ring>, //defined absolutely and factored modulo i^n
}

impl<Ring: EuclideanDomain + UniqueFactorizationDomain> HenselFactorization<Ring> {
    fn check(&self) -> Result<(), &'static str> {
        self.factors.check(&self.i, &self.n)
    }

    fn new(i: Ring, n: Natural, h: Polynomial<Ring>, fs: Vec<Polynomial<Ring>>) -> Self {
        // h and all fs monic
        assert!(fs.len() >= 1);
        debug_assert!(h.is_monic());
        for f in fs.iter() {
            debug_assert!(f.is_monic());
        }
        // h = product of fs modulo i^n
        debug_assert_eq!(
            polynomial_reduce_modulo(&h, &i, &n),
            polynomial_reduce_modulo(&Polynomial::product(fs.iter().collect()), &i, &n)
        );
        // fs are coprime mod i - checked when computing bezout coefficients
        let factors = HenselFactorizationImpl::new(&i, &n, h, fs.iter().collect());
        let ans = Self { i, n, factors };
        debug_assert!(ans.check().is_ok());
        ans
    }

    fn linear_lift(&mut self) {
        todo!()
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
        let hensel_fact_1 =
            HenselFactorization::new(Integer::from(5), Natural::from(1u8), h, vec![f1, f2, f3]);
        hensel_fact_1.check().unwrap();
        println!("{:#?}", hensel_fact_1);
    }

    #[test]
    fn test_hensel_factorization_example1_mod5_old() {
        //(4 + 19x + 18x^2 + 12x^3 + 7x^4 + 5x^5 + x^6) = (1 + 4x + x^2)(4 + 3x + 2x^2 + x^3 + x^4)
        //(4 + 4x + 3x^2 + 2x^3 + 2x^4 + 0x^5 + x^6) = (1 + 4x + x^2)(4 + 3x + 2x^2 + x^3 + x^4) mod 5

        let h = Polynomial::from_coeffs(vec![
            Integer::from(4 + 15),
            Integer::from(4 + 15),
            Integer::from(3 + 15),
            Integer::from(2 + 15),
            Integer::from(2 + 15),
            Integer::from(0 + 15),
            Integer::from(1),
        ]);
        let f = Polynomial::from_coeffs(vec![
            Integer::from(1 + 15),
            Integer::from(4 + 15),
            Integer::from(1),
        ]);
        let g = Polynomial::from_coeffs(vec![
            Integer::from(4 + 15),
            Integer::from(3 + 15),
            Integer::from(2 + 15),
            Integer::from(1 + 15),
            Integer::from(1),
        ]);

        let h_fact_1 = HenselPairFactorizationOld::new(Integer::from(5), h, f, g).unwrap();
        assert!(h_fact_1.check().is_ok());
        println!("{:?}", h_fact_1);

        let h_fact_2 = h_fact_1.linear_lift();
        assert!(h_fact_2.check().is_ok());
        println!("{:?}", h_fact_2);

        let h_fact_3 = h_fact_2.linear_lift();
        assert!(h_fact_3.check().is_ok());
        println!("{:?}", h_fact_3);

        let h_fact_4 = h_fact_3.linear_lift();
        assert!(h_fact_4.check().is_ok());
        println!("{:?}", h_fact_4);

        let h_fact_5 = h_fact_4.linear_lift();
        assert!(h_fact_5.check().is_ok());
        println!("{:?}", h_fact_5);
    }

    #[test]
    fn test_hensel_factorization_example2_mod5_old() {
        //(4 + 19x + 18x^2 + 12x^3 + 7x^4 + 5x^5 + x^6) = (1 + 4x + x^2)(4 + 3x + 2x^2 + x^3 + x^4)
        //(4 + 19x + 18x^2 + 12x^3 + 7x^4 + 5x^5 + x^6) = (1 + 4x + x^2)(4 + 3x + 2x^2 + x^3 + x^4) mod 5

        //in this exactly we are lifting a factorization which works in Z, so the lifted factors are identical

        let h = Polynomial::from_coeffs(vec![
            Integer::from(4),
            Integer::from(19),
            Integer::from(18),
            Integer::from(12),
            Integer::from(7),
            Integer::from(5),
            Integer::from(1),
        ]);
        let f = Polynomial::from_coeffs(vec![Integer::from(1), Integer::from(4), Integer::from(1)]);
        let g = Polynomial::from_coeffs(vec![
            Integer::from(4),
            Integer::from(3),
            Integer::from(2),
            Integer::from(1),
            Integer::from(1),
        ]);

        let h_fact_1 = HenselPairFactorizationOld::new(Integer::from(5), h, f, g).unwrap();
        assert!(h_fact_1.check().is_ok());
        println!("{:?}", h_fact_1);

        let h_fact_2 = h_fact_1.linear_lift();
        assert!(h_fact_2.check().is_ok());
        println!("{:?}", h_fact_2);

        let h_fact_3 = h_fact_2.linear_lift();
        assert!(h_fact_3.check().is_ok());
        println!("{:?}", h_fact_3);

        let h_fact_4 = h_fact_3.linear_lift();
        assert!(h_fact_4.check().is_ok());
        println!("{:?}", h_fact_4);
    }

    // #[test]
    // fn test_hensel_factorization_not_coprime() {
    //     //h(x) = f(x)h(x) = (1 + 4x + x^2)(1 + 4x + x^2) -> error because not coprime

    //     let f = Polynomial::from_coeffs(vec![Integer::from(1), Integer::from(4), Integer::from(1)]);
    //     let g = Polynomial::from_coeffs(vec![Integer::from(1), Integer::from(4), Integer::from(1)]);
    //     assert!(HenselFactorization::new(Integer::from(5), f, g).is_err());
    // }

    // #[test]
    // fn test_hensel_linear_lift() {
    //     //h(x) = f(x)h(x) = (1 + 4x + x^2)(4 + 3x + 2x^2 + x^3 + x^4)

    //     let f =
    //         Polynomial::from_coeffs(vec![Integer::from(1), Integer::from(-1), Integer::from(1)]);
    //     let g = Polynomial::from_coeffs(vec![
    //         Integer::from(4),
    //         Integer::from(3),
    //         Integer::from(2),
    //         Integer::from(1),
    //         Integer::from(1),
    //     ]);
    //     let h_fact_1 = HenselFactorization::new(Integer::from(5), f, g).unwrap();
    //     assert!(h_fact_1.check().is_ok());
    //     println!("{:?}", h_fact_1);

    //     let h_fact_2 = h_fact_1.linear_lift();
    //     println!("{:?}", h_fact_2);
    // }
}
