use std::collections::HashMap;

use malachite_base::num::basic::traits::{One, Two, Zero};
use malachite_nz::natural::Natural;

pub fn nat_to_usize(n: &Natural) -> Result<usize, ()> {
    let limbs = n.to_limbs_asc();
    if limbs.len() == 0 {
        Ok(0)
    } else if limbs.len() == 1 {
        let n = limbs[0];
        if Natural::from(n) > Natural::from(usize::MAX) {
            Err(())
        } else {
            Ok(n as usize)
        }
    } else {
        Err(())
    }
}

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

pub fn factor_natural_by_try_primes(mut n: Natural) -> HashMap<Natural, Natural> {
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
        for _i in 0..a {
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

pub fn pow(x: &Natural, n: &Natural) -> Natural {
    use malachite_base::num::logic::traits::BitIterable;
    if *n == 0 {
        Natural::ONE
    } else if *n == 1 {
        x.clone()
    } else {
        debug_assert!(*n >= 2);
        let bits: Vec<_> = n.bits().collect();
        let mut pows = vec![x.clone()];
        while pows.len() < bits.len() {
            pows.push(pows.last().unwrap() * pows.last().unwrap());
        }
        let count = bits.len();
        debug_assert_eq!(count, pows.len());
        let mut ans = Natural::ONE;
        for i in 0..count {
            if bits[i] {
                ans *= &pows[i];
            }
        }
        ans
    }
}

/// Compute the floor of the nth root of a
pub fn nth_root(x: &Natural, n: &Natural) -> Natural {
    if n == &Natural::ZERO {
        panic!()
    } else if n == &Natural::ONE {
        x.clone()
    } else if x == &Natural::ZERO {
        Natural::ZERO
    } else {
        let mut a = Natural::ONE;
        let mut b = x.clone();
        while &a + &Natural::ONE < b {
            let m = (&a + &b) / Natural::TWO;
            if pow(&m, &n) <= *x {
                a = m;
            } else {
                b = m;
            }
        }
        a
    }
}

/// Return the number of bits needed to store n i.e. ceil(log2(n)) for all non-zero n
pub fn bitcount(n: &Natural) -> usize {
    use malachite_base::num::logic::traits::BitIterable;
    n.bits().len()
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IsPowerTestResult {
    Zero,
    One,
    Power(Natural, Natural),
    No,
}

pub fn is_power_test(n: &Natural) -> IsPowerTestResult {
    if *n == Natural::ZERO {
        IsPowerTestResult::Zero
    } else if *n == Natural::ONE {
        IsPowerTestResult::One
    } else {
        for k in 2..bitcount(n) {
            let k = Natural::from(k);
            let a = nth_root(n, &k);
            if *n == pow(&a, &k) {
                return IsPowerTestResult::Power(a, k);
            }
        }
        IsPowerTestResult::No
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PrimalityTestResult {
    Zero,
    One,
    Prime,
    Composite,
    Unknown,
}

// https://cr.yp.to/papers/aks.pdf
pub fn aks_primality_test(n: &Natural) -> PrimalityTestResult {
    if *n == Natural::ZERO {
        PrimalityTestResult::Zero
    } else if *n == Natural::ONE {
        PrimalityTestResult::One
    } else {
        PrimalityTestResult::Unknown
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use super::*;

    #[test]
    fn test_primality_test() {
        println!("{:?}", aks_primality_test(&Natural::from_str("6").unwrap()));
    }

    #[test]
    fn test_pow() {
        assert_eq!(
            pow(&Natural::from(0usize), &Natural::from(0usize)),
            Natural::from(1usize)
        );
        assert_eq!(
            pow(&Natural::from(0usize), &Natural::from(1usize)),
            Natural::from(0usize)
        );
        assert_eq!(
            pow(&Natural::from(0usize), &Natural::from(2usize)),
            Natural::from(0usize)
        );
        assert_eq!(
            pow(&Natural::from(0usize), &Natural::from(3usize)),
            Natural::from(0usize)
        );

        assert_eq!(
            pow(&Natural::from(1usize), &Natural::from(0usize)),
            Natural::from(1usize)
        );
        assert_eq!(
            pow(&Natural::from(1usize), &Natural::from(1usize)),
            Natural::from(1usize)
        );
        assert_eq!(
            pow(&Natural::from(1usize), &Natural::from(2usize)),
            Natural::from(1usize)
        );
        assert_eq!(
            pow(&Natural::from(1usize), &Natural::from(3usize)),
            Natural::from(1usize)
        );

        assert_eq!(
            pow(&Natural::from(2usize), &Natural::from(0usize)),
            Natural::from(1usize)
        );
        assert_eq!(
            pow(&Natural::from(2usize), &Natural::from(1usize)),
            Natural::from(2usize)
        );
        assert_eq!(
            pow(&Natural::from(2usize), &Natural::from(2usize)),
            Natural::from(4usize)
        );
        assert_eq!(
            pow(&Natural::from(2usize), &Natural::from(3usize)),
            Natural::from(8usize)
        );

        assert_eq!(
            pow(&Natural::from(3usize), &Natural::from(0usize)),
            Natural::from(1usize)
        );
        assert_eq!(
            pow(&Natural::from(3usize), &Natural::from(1usize)),
            Natural::from(3usize)
        );
        assert_eq!(
            pow(&Natural::from(3usize), &Natural::from(2usize)),
            Natural::from(9usize)
        );
        assert_eq!(
            pow(&Natural::from(3usize), &Natural::from(3usize)),
            Natural::from(27usize)
        );
    }

    #[test]
    fn test_nth_root() {
        for n in 1..10usize {
            for x in 0..10usize {
                let x = Natural::from(x);
                let n = Natural::from(n);
                let r = nth_root(&x, &n);
                println!("{}th root of {} is {}", &n, &x, &r);
                assert!(pow(&r, &n) <= x);
                assert!(x == 0 || pow(&(r + Natural::ONE), &n) > x);
            }
        }
    }

    #[test]
    fn test_is_power_test() {
        assert_eq!(
            is_power_test(&Natural::from(0usize)),
            IsPowerTestResult::Zero
        );
        assert_eq!(
            is_power_test(&Natural::from(1usize)),
            IsPowerTestResult::One
        );
        assert_eq!(is_power_test(&Natural::from(2usize)), IsPowerTestResult::No);
        assert_eq!(is_power_test(&Natural::from(3usize)), IsPowerTestResult::No);
        assert_eq!(
            is_power_test(&Natural::from(4usize)),
            IsPowerTestResult::Power(Natural::from(2usize), Natural::from(2usize))
        );
        assert_eq!(is_power_test(&Natural::from(5usize)), IsPowerTestResult::No);
        assert_eq!(is_power_test(&Natural::from(6usize)), IsPowerTestResult::No);
        assert_eq!(is_power_test(&Natural::from(7usize)), IsPowerTestResult::No);
        assert_eq!(
            is_power_test(&Natural::from(8usize)),
            IsPowerTestResult::Power(Natural::from(2usize), Natural::from(3usize))
        );
        assert_eq!(
            is_power_test(&Natural::from(9usize)),
            IsPowerTestResult::Power(Natural::from(3usize), Natural::from(2usize))
        );
        assert_eq!(
            is_power_test(&Natural::from(10usize)),
            IsPowerTestResult::No
        );

        println!("{:?}", is_power_test(&Natural::from(0usize)));
    }

    #[test]
    fn test_nat_to_usize() {
        assert_eq!(nat_to_usize(&Natural::from(0u8)).unwrap(), 0);
        assert_eq!(nat_to_usize(&Natural::from(1u8)).unwrap(), 1);
        assert_eq!(nat_to_usize(&Natural::from(2u8)).unwrap(), 2);
    }

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
