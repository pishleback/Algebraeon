use super::factor;
use crate::{
    rings::natural::factorization::primes::is_prime,
    structure::{Factored, FactoredSignature, SemiRingSignature},
};
use algebraeon_nzq::{Natural, NaturalCanonicalStructure, gcd, traits::ModPow};
use algebraeon_sets::structure::MetaType;
use itertools::Itertools;
use std::{borrow::Borrow, collections::HashMap};

impl FactoredSignature<FactoredNatural> for NaturalCanonicalStructure {
    type PrimeObject = Natural;

    type Object = Natural;

    fn object_divides(&self, a: &Self::Object, b: &Self::Object) -> bool {
        b % a == Natural::ZERO
    }

    fn object_is_prime(&self, object: &Self::PrimeObject) -> bool {
        is_prime(object)
    }

    fn prime_to_object(&self, prime: Self::PrimeObject) -> Self::Object {
        prime
    }

    fn object_product(&self, objects: Vec<&Self::Object>) -> Self::Object {
        self.product(objects)
    }
}

#[derive(Debug, Clone)]
pub struct FactoredNatural {
    primes: HashMap<Natural, Natural>,
}

impl Factored for FactoredNatural {
    type Structure = NaturalCanonicalStructure;

    fn factored_structure<'a>(&'a self) -> impl 'a + Borrow<Self::Structure> {
        Natural::structure()
    }

    fn from_factor_powers_impl(
        _structure: Self::Structure,
        factor_powers: Vec<(Natural, Natural)>,
    ) -> Self {
        Self {
            primes: factor_powers.into_iter().collect(),
        }
    }

    fn factor_powers(&self) -> Vec<(&Natural, &Natural)> {
        self.primes.iter().collect()
    }

    fn into_factor_powers(self) -> Vec<(Natural, Natural)> {
        self.primes.into_iter().collect()
    }

    fn expanded(&self) -> Natural {
        let mut t = Natural::ONE;
        for (p, k) in &self.primes {
            t *= p.pow(k);
        }
        t
    }

    fn mul(mut a: Self, b: Self) -> Self {
        for (p, k) in b.into_factor_powers() {
            *a.primes.entry(p).or_insert(Natural::ZERO) += k;
        }
        a
    }
}

impl std::fmt::Display for FactoredNatural {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.primes.is_empty() {
            write!(f, "1")?;
        } else {
            for (i, (p, k)) in self
                .primes
                .iter()
                .sorted_by_cached_key(|(p, _k)| (*p).clone())
                .enumerate()
            {
                if i != 0 {
                    write!(f, " × ")?;
                }
                write!(f, "{}", p)?;
                if k != &Natural::ONE {
                    write!(f, "^")?;
                    write!(f, "{}", k)?;
                }
            }
        }
        Ok(())
    }
}

impl FactoredNatural {
    pub fn mul_prime(&mut self, p: Natural) {
        debug_assert!(is_prime(&p));
        *self.primes.entry(p).or_insert(Natural::ZERO) += Natural::ONE;
    }

    pub fn euler_totient(&self) -> Natural {
        let mut t = Natural::ONE;
        for (p, k) in &self.primes {
            t *= (p - &Natural::ONE) * p.pow(&(k - &Natural::ONE));
        }
        t
    }

    pub fn distinct_prime_factors(&self) -> Vec<&Natural> {
        let mut primes = vec![];
        for (p, _) in &self.primes {
            primes.push(p);
        }
        primes
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IsPrimitiveRootResult {
    NonUnit,
    No,
    Yes,
}
impl FactoredNatural {
    /// Return whether x is a primitive root modulo the value represented by self
    pub fn is_primitive_root(&self, x: &Natural) -> IsPrimitiveRootResult {
        let n_factored = self;
        let n = n_factored.expanded();
        if gcd(x.clone(), n.clone()) != Natural::ONE {
            IsPrimitiveRootResult::NonUnit
        } else {
            let phi_n = n_factored.euler_totient();
            let x_mod_n = x % &n;
            for p in factor(phi_n.clone()).unwrap().distinct_prime_factors() {
                if (&x_mod_n).mod_pow(&phi_n / p, &n) == Natural::ONE {
                    return IsPrimitiveRootResult::No;
                }
            }
            IsPrimitiveRootResult::Yes
        }
    }
}
