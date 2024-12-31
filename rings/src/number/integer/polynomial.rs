use malachite_nz::integer::Integer;

use crate::{
    number::integer::*,
    number::natural::{functions::*, primes::*},
    polynomial::polynomial::*,
    ring_structure::quotient::*,
};

impl GreatestCommonDivisorStructure for PolynomialStructure<CannonicalStructure<Integer>> {
    fn gcd(&self, x: &Self::Set, y: &Self::Set) -> Self::Set {
        self.gcd_by_primitive_subresultant(x.clone(), y.clone())
    }
}

impl UniqueFactorizationStructure for PolynomialStructure<CannonicalStructure<Integer>> {
    fn factor(&self, p: &Self::Set) -> Option<Factored<Self>> {
        // self.factorize_by_kroneckers_method(p)
        Polynomial::factorize_by_zassenhaus_algorithm(p.clone())
    }
}

impl Polynomial<Integer> {
    fn find_factor_primitive_sqfree_by_zassenhaus_algorithm(
        self,
    ) -> Option<(Polynomial<Integer>, Polynomial<Integer>)> {
        // println!("{}", self);

        //ways to improve this:
        //1) use LLL to find the right combination of modular factors
        //2) don't use factor by find factor, just do it all in one go and partition the modular factors into true factors

        //This might be of use when I do the LLL bit, especially 10 and 11: http://people.csail.mit.edu/madhu/FT98/
        let f = self;
        let f_deg = f.degree().unwrap();
        debug_assert_ne!(f_deg, 0);

        //https://en.wikipedia.org/wiki/Landau-Mignotte_bound
        //find a bound B such that the coefficients of any factors of B have absolute value bounded by B

        //use B = (n choose floor(n/2)) * (sum of absolute value of coeffs of f)

        let mut factor_coeff_bound = Natural::ZERO;
        for c in f.coeffs() {
            factor_coeff_bound += c.unsigned_abs();
        }
        factor_coeff_bound *= choose_usize(f_deg, f_deg / 2);
        //for the lifted factors we have
        // 0 <= coeff < hensel_factorization_f_over_p.modolus()
        //so replacing the coefficient with a different representative if necessary we have
        // -floor(hensel_factorization_f_over_p.modolus()/2) <= coeff <= floor(hensel_factorization_f_over_p.modolus()/2)
        // i.e. |coeff| <= floor(hensel_factorization_f_over_p.modolus()/2)
        // we know that coefficients of any integer factor of f satisfies
        // |coeff| <= factor_coeff_bound
        //so we want
        // factor_coeff_bound <= floor(hensel_factorization_f_over_p.modolus()/2)
        // i.e. 2 * factor_coeff_bound <= hensel_factorization_f_over_p.modolus()
        let target_modolus = Natural::TWO * factor_coeff_bound;

        // println!("zassenhaus: {}", f);
        if f_deg == 1 {
            None
        } else {
            let prime_gen = PrimeGenerator::new();
            for p in prime_gen {
                let mod_p = QuotientStructure::new_field(Integer::structure(), Integer::from(&p));
                let poly_mod_p = PolynomialStructure::new(mod_p.into());

                if poly_mod_p.degree(&f).unwrap() == f_deg {
                    let facotred_f_mod_p = poly_mod_p.factor(&f).unwrap();
                    // println!("f mod {} = {:?}", p, facotred_f_mod_p);
                    match facotred_f_mod_p.into_hensel_factorization(f.clone()) {
                        Some(mut hensel_factorization_f_over_p) => {
                            // println!("good prime = {}", p);
                            // println!("{:?}", hensel_factorization_f_over_p.factors());

                            while hensel_factorization_f_over_p.modolus() < target_modolus {
                                hensel_factorization_f_over_p.linear_lift();
                            }

                            let modulus = hensel_factorization_f_over_p.modolus();

                            let lifted_factors = hensel_factorization_f_over_p.factors();

                            // println!("{:?}", hensel_factorization_f_over_p.modolus());
                            // println!("{:?}", lifted_factors);
                            // println!("lifted_factors.len() = {:?}", lifted_factors.len());

                            for subset in (0..lifted_factors.len())
                                .map(|_i| vec![false, true])
                                .multi_cartesian_product()
                            {
                                let possible_factor = Polynomial::mul(
                                    &Polynomial::constant(f.leading_coeff().unwrap()),
                                    &Polynomial::product(
                                        (0..lifted_factors.len())
                                            .filter(|i| subset[*i])
                                            .map(|i| lifted_factors[i])
                                            .collect(),
                                    ),
                                )
                                .apply_map(|c| {
                                    let c = Integer::rem(c, &modulus);
                                    if c > Integer::quo(&modulus, &Integer::TWO).unwrap() {
                                        c - &modulus
                                    } else {
                                        c.clone()
                                    }
                                })
                                .primitive_part() //factoring f(x) = 49x^2-10000 had possible_factor = 49x-700, which is only a factor over the rationals and not over the integers unless we take the primitive part which is 7x-100, soo this seems to make sense though I cant properly justify it right now.
                                .unwrap();

                                // println!("possible_factor = {:?}", possible_factor);

                                // println!("{:?}", Polynomial::div(&f, &possible_factor));

                                if possible_factor.degree().unwrap() != 0
                                    && possible_factor.degree().unwrap() != f_deg
                                {
                                    match Polynomial::div(&f, &possible_factor) {
                                        Ok(other_factor) => {
                                            return Some((possible_factor, other_factor))
                                        }
                                        Err(RingDivisionError::NotDivisible) => {}
                                        Err(RingDivisionError::DivideByZero) => unreachable!(),
                                    }
                                }
                            }
                            return None;
                        }
                        None => {}
                    }
                }
            }
            unreachable!("Because there are infinitely many primes")
        }
    }

    pub fn factorize_by_zassenhaus_algorithm(
        self,
    ) -> Option<Factored<PolynomialStructure<CannonicalStructure<Integer>>>> {
        if self == Self::zero() {
            None
        } else {
            Some(
                Self::structure().factorize_by_primitive_sqfree_factorize_by_yuns_algorithm(
                    self,
                    &|f| {
                        factorize_by_find_factor(&Self::structure(), f, &|f| {
                            Self::find_factor_primitive_sqfree_by_zassenhaus_algorithm(f)
                        })
                    },
                ),
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::elements::*;

    #[test]
    fn test_zassenhaus_against_kroneckers() {
        let x = &Polynomial::<Integer>::var().into_ring();

        let f = ((2 * x.pow(3) + 6 * x.pow(2) - 4) * (6 * x.pow(5) + 7 * x.pow(4) - 4)).into_set();
        println!("{}", f);
        // println!("{}", f.clone().factorize_by_kroneckers_method().unwrap());
        // println!("{}", f.clone().factorize_by_zassenhaus_algorithm().unwrap());
        assert!(Factored::equal(
            &f.clone().factorize_by_kroneckers_method().unwrap(),
            &f.clone().factorize_by_zassenhaus_algorithm().unwrap()
        ));

        let f = (49 * x.pow(2) - 10000).into_set();
        println!("{}", f);
        // println!("{}", f.clone().factorize_by_kroneckers_method().unwrap());
        // println!("{}", f.clone().factorize_by_zassenhaus_algorithm().unwrap());
        assert!(Factored::equal(
            &f.clone().factorize_by_kroneckers_method().unwrap(),
            &f.clone().factorize_by_zassenhaus_algorithm().unwrap()
        ));
    }
}
