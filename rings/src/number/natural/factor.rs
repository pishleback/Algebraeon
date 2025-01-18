use super::functions::*;
use super::*;

#[derive(Debug, Clone)]
pub struct Factored {
    primes: HashMap<Natural, Natural>,
}

impl Factored {
    pub fn new_unchecked(primes: HashMap<Natural, Natural>) -> Self {
        for (_p, k) in &primes {
            // TODO
            // debug_assert!(is_prime(p));
            debug_assert!(*k > 0);
        }
        Self { primes }
    }

    pub fn from_prime_unchecked(prime: Natural) -> Self {
        Self::new_unchecked(HashMap::from([(prime, Natural::ONE)]))
    }

    pub fn powers(&self) -> &HashMap<Natural, Natural> {
        &self.primes
    }

    pub fn into_powers(self) -> HashMap<Natural, Natural> {
        self.primes
    }

    pub fn expand(&self) -> Natural {
        let mut t = Natural::ONE;
        for (p, k) in &self.primes {
            t *= pow(p, k);
        }
        t
    }

    pub fn euler_totient(&self) -> Natural {
        let mut t = Natural::ONE;
        for (p, k) in &self.primes {
            t *= (p - &Natural::ONE) * pow(p, &(k - &Natural::ONE));
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
impl Factored {
    /// Return whether x is a primitive root modulo the value represented by self
    pub fn is_primitive_root(&self, x: &Natural) -> IsPrimitiveRootResult {
        use malachite_base::num::arithmetic::traits::{Mod, ModPow};

        let n_factored = self;
        let n = n_factored.expand();
        if gcd(x.clone(), n.clone()) != Natural::ONE {
            IsPrimitiveRootResult::NonUnit
        } else {
            let phi_n = n_factored.euler_totient();
            let x_mod_n = x.mod_op(&n);
            for p in factor_by_try_divisors(phi_n.clone())
                .unwrap()
                .distinct_prime_factors()
            {
                if (&x_mod_n).mod_pow(&phi_n / p, &n) == Natural::ONE {
                    return IsPrimitiveRootResult::No;
                }
            }
            IsPrimitiveRootResult::Yes
        }
    }
}

pub fn factor_by_try_divisors(mut n: Natural) -> Option<Factored> {
    if n == Natural::ZERO {
        None
    } else {
        debug_assert_ne!(n, 0);
        let mut fs = HashMap::new();
        let mut d = Natural::TWO;
        loop {
            while &n % &d == 0 {
                *fs.entry(d.clone()).or_insert(Natural::from(0u8)) += Natural::from(1u8);
                n /= &d;
            }
            if n == 1 {
                break;
            }
            if &d * &d > n {
                fs.insert(n.clone(), Natural::ONE);
                break;
            }
            d += Natural::ONE;
        }
        Some(Factored::new_unchecked(fs))
    }
}

pub fn factor(n: &Natural) -> Option<Factored> {
    factor_by_try_divisors(n.clone())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_factor_natural_by_try_primes() {
        println!("{:?}", factor_by_try_divisors(Natural::from(12usize)));
        assert_eq!(
            factor_by_try_divisors(Natural::from(12usize))
                .unwrap()
                .into_powers(),
            HashMap::from([
                (Natural::from(2usize), Natural::from(2usize)),
                (Natural::from(3usize), Natural::from(1usize))
            ])
        );

        println!("{:?}", factor_by_try_divisors(Natural::from(100000001usize)));
    }

    #[test]
    fn test_euler_totient() {
        assert_eq!(
            factor_by_try_divisors(Natural::from(12usize))
                .unwrap()
                .euler_totient(),
            Natural::from(4usize)
        );
    }

    #[test]
    fn test_is_primitive_root() {
        assert_eq!(
            factor_by_try_divisors(Natural::from(761usize))
                .unwrap()
                .is_primitive_root(&Natural::from(0usize)),
            IsPrimitiveRootResult::NonUnit,
        );
        assert_eq!(
            factor_by_try_divisors(Natural::from(761usize))
                .unwrap()
                .is_primitive_root(&Natural::from(1usize)),
            IsPrimitiveRootResult::No,
        );
        assert_eq!(
            factor_by_try_divisors(Natural::from(761usize))
                .unwrap()
                .is_primitive_root(&Natural::from(2usize)),
            IsPrimitiveRootResult::No,
        );
        assert_eq!(
            factor_by_try_divisors(Natural::from(761usize))
                .unwrap()
                .is_primitive_root(&Natural::from(3usize)),
            IsPrimitiveRootResult::No,
        );
        assert_eq!(
            factor_by_try_divisors(Natural::from(761usize))
                .unwrap()
                .is_primitive_root(&Natural::from(4usize)),
            IsPrimitiveRootResult::No,
        );
        assert_eq!(
            factor_by_try_divisors(Natural::from(761usize))
                .unwrap()
                .is_primitive_root(&Natural::from(5usize)),
            IsPrimitiveRootResult::No,
        );
        assert_eq!(
            factor_by_try_divisors(Natural::from(761usize))
                .unwrap()
                .is_primitive_root(&Natural::from(6usize)),
            IsPrimitiveRootResult::Yes,
        );
    }
}
