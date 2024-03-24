use std::{collections::HashMap, f32::NAN};

use malachite_base::num::basic::traits::{One, Zero};
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

pub fn choose(a: usize, b: usize) -> Natural {
    if b > a {
        Natural::ZERO
    } else {
        let mut vals = vec![Natural::ONE];
        for i in 0..a {
            vals = {
                let mut next_vals = vec![Natural::ONE];
                for i in 0..vals.len() - 1 {
                    next_vals.push(&vals[i] + &vals[i + 1]);
                }
                next_vals.push(Natural::ONE);
                next_vals
            };
        }
        vals.into_iter().nth(b).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_choose() {
        assert_eq!(choose(0, 0), Natural::from(1u32));

        assert_eq!(choose(1, 0), Natural::from(1u32));
        assert_eq!(choose(1, 1), Natural::from(1u32));

        assert_eq!(choose(2, 0), Natural::from(1u32));
        assert_eq!(choose(2, 1), Natural::from(2u32));
        assert_eq!(choose(2, 2), Natural::from(1u32));

        assert_eq!(choose(3, 0), Natural::from(1u32));
        assert_eq!(choose(3, 1), Natural::from(3u32));
        assert_eq!(choose(3, 2), Natural::from(3u32));
        assert_eq!(choose(3, 3), Natural::from(1u32));

        assert_eq!(choose(4, 0), Natural::from(1u32));
        assert_eq!(choose(4, 1), Natural::from(4u32));
        assert_eq!(choose(4, 2), Natural::from(6u32));
        assert_eq!(choose(4, 3), Natural::from(4u32));
        assert_eq!(choose(4, 4), Natural::from(1u32));

        assert_eq!(choose(3, 4), Natural::from(0u32));
    }
}
