use std::collections::HashMap;

use malachite_nz::natural::Natural;

pub struct NaturalPrimeGenerator {
    n: Natural,
    primes: Vec<Natural>,
}

impl NaturalPrimeGenerator {
    pub fn new() -> Self {
        Self {
            n: Natural::from(2u8),
            primes: vec![],
        }
    }
}

impl Iterator for NaturalPrimeGenerator {
    type Item = Natural;

    fn next(&mut self) -> Option<Self::Item> {
        'next_loop: loop {
            //todo: only check primes up to sqrt n
            for p in &self.primes {
                if &self.n % p == 0 {
                    self.n += Natural::from(1u8);
                    continue 'next_loop;
                }
            }
            let next_p = self.n.clone();
            self.n += Natural::from(1u8);
            self.primes.push(next_p.clone());
            return Some(next_p);
        }
    }
}

pub fn factorial(n: Natural) -> Natural {
    let mut k = Natural::from(1u8);
    let mut i = Natural::from(1u8);
    while i <= n {
        k *= &i;
        i += Natural::from(1u8);
    }
    k
}

pub fn factor_nat(mut n: Natural) -> HashMap<Natural, Natural> {
    //TODO: more efficient implementations
    assert_ne!(n, 0);
    let mut fs = HashMap::new();
    let mut p = Natural::from(2u8);
    while n > 1 && p <= n {
        while &n % &p == 0 {
            *fs.entry(p.clone()).or_insert(Natural::from(0u8)) += Natural::from(1u8);
            n /= &p;
        }
        p += Natural::from(1u8);
    }
    fs
}