use super::factor;
use crate::{
    linear::ordered_set_free_module::FreeModuleOverOrderedSetStructure,
    rings::natural::factorization::primes::is_prime,
    structure::{AdditiveMonoidSignature, FactoredSignature, SemiRingSignature},
};
use algebraeon_nzq::{Natural, NaturalCanonicalStructure, gcd, traits::ModPow};
use algebraeon_sets::structure::*;
use itertools::Itertools;

pub trait NaturalCanonicalFactorizationStructure {
    fn factorizations(&self) -> NaturalFactorizationStructure {
        NaturalFactorizationStructure {}
    }
}
impl NaturalCanonicalFactorizationStructure for NaturalCanonicalStructure {}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NaturalFactorizationStructure {}

impl NaturalFactorizationStructure {
    fn powers_semimodule(
        &self,
    ) -> FreeModuleOverOrderedSetStructure<
        NaturalCanonicalStructure,
        NaturalCanonicalStructure,
        NaturalCanonicalStructure,
        NaturalCanonicalStructure,
    > {
        FreeModuleOverOrderedSetStructure::new(Natural::structure(), Natural::structure())
    }
}

impl Signature for NaturalFactorizationStructure {}

impl SetSignature for NaturalFactorizationStructure {
    type Set = Vec<(Natural, Natural)>;

    fn is_element(&self, x: &Self::Set) -> bool {
        if !self.powers_semimodule().is_element(x) {
            return false;
        }
        for (prime, _) in self.to_powers_unchecked(x) {
            if !is_prime(prime) {
                return false;
            }
        }
        true
    }
}

impl FactoredSignature for NaturalFactorizationStructure {
    type PrimeObject = Natural;
    type Object = Natural;

    fn object_divides(&self, a: &Self::Object, b: &Self::Object) -> bool {
        b % a == Natural::ZERO
    }

    fn try_object_is_prime(&self, object: &Self::PrimeObject) -> Option<bool> {
        Some(is_prime(object))
    }

    fn prime_into_object(&self, prime: Self::PrimeObject) -> Self::Object {
        prime
    }

    fn object_product(&self, objects: Vec<&Self::Object>) -> Self::Object {
        Natural::structure().product(objects)
    }

    fn new_powers_unchecked(&self, factor_powers: Vec<(Natural, Natural)>) -> Self::Set {
        factor_powers
    }

    fn to_powers_unchecked<'a>(&self, a: &'a Self::Set) -> Vec<(&'a Natural, &'a Natural)> {
        a.iter().map(|(p, k)| (p, k)).collect()
    }

    fn into_powers_unchecked(&self, a: Self::Set) -> Vec<(Natural, Natural)> {
        a
    }

    fn expanded(&self, a: &Self::Set) -> Natural {
        let mut t = Natural::ONE;
        for (p, k) in a {
            t *= p.pow(k);
        }
        t
    }

    fn mul(&self, a: Self::Set, b: Self::Set) -> Self::Set {
        self.powers_semimodule().add(&a, &b)
    }
}

impl ToStringSignature for NaturalFactorizationStructure {
    fn to_string(&self, elem: &Self::Set) -> String {
        use std::fmt::Write;
        let mut f = String::new();
        if elem.is_empty() {
            write!(f, "1").unwrap();
        } else {
            for (i, (p, k)) in elem
                .iter()
                .sorted_by_cached_key(|(p, _k)| (*p).clone())
                .enumerate()
            {
                if i != 0 {
                    write!(f, " × ").unwrap();
                }
                write!(f, "{}", p).unwrap();
                if k != &Natural::ONE {
                    write!(f, "^").unwrap();
                    write!(f, "{}", k).unwrap();
                }
            }
        }
        f
    }
}

impl NaturalFactorizationStructure {
    pub fn mul_prime(&self, f: &mut Vec<(Natural, Natural)>, p: Natural) {
        debug_assert!(is_prime(&p));
        *f = self.powers_semimodule().add(f, &vec![(p, Natural::ONE)]);
    }

    pub fn euler_totient(&self, f: &Vec<(Natural, Natural)>) -> Natural {
        let mut t = Natural::ONE;
        for (p, k) in f {
            t *= (p - &Natural::ONE) * p.pow(&(k - &Natural::ONE));
        }
        t
    }

    pub fn distinct_prime_factors<'a>(&self, f: &'a Vec<(Natural, Natural)>) -> Vec<&'a Natural> {
        let mut primes = vec![];
        for (p, _) in f {
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
impl NaturalFactorizationStructure {
    /// Return whether x is a primitive root modulo the factorized value
    pub fn is_primitive_root(
        &self,
        x: &Natural,
        n_factored: &Vec<(Natural, Natural)>,
    ) -> IsPrimitiveRootResult {
        let factorizations = Natural::structure().factorizations();
        let n = factorizations.expanded(n_factored);
        if gcd(x.clone(), n.clone()) != Natural::ONE {
            IsPrimitiveRootResult::NonUnit
        } else {
            let phi_n = factorizations.euler_totient(n_factored);
            let x_mod_n = x % &n;
            for p in factorizations.distinct_prime_factors(&factor(phi_n.clone()).unwrap()) {
                if (&x_mod_n).mod_pow(&phi_n / p, &n) == Natural::ONE {
                    return IsPrimitiveRootResult::No;
                }
            }
            IsPrimitiveRootResult::Yes
        }
    }
}
