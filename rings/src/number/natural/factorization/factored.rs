use algebraeon_nzq::Natural;
use itertools::Itertools;
use std::collections::HashMap;

use crate::number::natural::{
    factorization::primes::is_prime,
    functions::{gcd, pow},
};

use super::factor;

#[derive(Debug, Clone)]
pub struct FactoredNatural {
    primes: HashMap<Natural, usize>,
}

impl std::fmt::Display for FactoredNatural {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.primes.is_empty() {
            write!(f, "1")?;
        } else {
            for (i, (p, &k)) in self
                .primes
                .iter()
                .sorted_by_cached_key(|(p, _k)| (*p).clone())
                .enumerate()
            {
                if i != 0 {
                    write!(f, " × ")?;
                }
                write!(f, "{}", p)?;
                if k != 1 {
                    write!(f, "^")?;
                    write!(f, "{}", k)?;
                }
            }
        }
        Ok(())
    }
}

impl FactoredNatural {
    pub fn new_unchecked(primes: HashMap<Natural, usize>) -> Self {
        for (p, k) in &primes {
            debug_assert!(is_prime(p));
            debug_assert!(k > &0);
        }
        Self { primes }
    }

    pub fn one() -> Self {
        Self {
            primes: HashMap::new(),
        }
    }

    pub fn mul_prime(&mut self, p: Natural) {
        debug_assert!(is_prime(&p));
        *self.primes.entry(p).or_insert(0) += 1;
    }

    pub fn from_prime_unchecked(prime: Natural) -> Self {
        Self::new_unchecked(HashMap::from([(prime, 1)]))
    }

    pub fn powers(&self) -> &HashMap<Natural, usize> {
        &self.primes
    }

    pub fn into_powers(self) -> HashMap<Natural, usize> {
        self.primes
    }

    pub fn expand(&self) -> Natural {
        let mut t = Natural::ONE;
        for (p, k) in &self.primes {
            t *= pow(p, &(*k).into());
        }
        t
    }

    pub fn euler_totient(&self) -> Natural {
        let mut t = Natural::ONE;
        for (p, k) in &self.primes {
            t *= (p - &Natural::ONE) * pow(p, &(k - 1).into());
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
        let n = n_factored.expand();
        if gcd(x.clone(), n.clone()) != Natural::ONE {
            IsPrimitiveRootResult::NonUnit
        } else {
            let phi_n = n_factored.euler_totient();
            let x_mod_n = x % &n;
            for p in factor(phi_n.clone()).unwrap().distinct_prime_factors() {
                if x_mod_n.mod_pow_ref(&phi_n / p, &n) == Natural::ONE {
                    return IsPrimitiveRootResult::No;
                }
            }
            IsPrimitiveRootResult::Yes
        }
    }
}
