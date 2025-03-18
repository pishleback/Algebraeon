use malachite_base::num::arithmetic::traits::AbsDiff;
use primes::is_prime;

use crate::polynomial::polynomial::Polynomial;

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
            debug_assert!(k > &Natural::ZERO);
        }
        Self { primes }
    }

    pub fn one() -> Self {
        Self {
            primes: HashMap::new(),
        }
    }

    fn mul_prime(&mut self, p: Natural) {
        *self.primes.entry(p).or_insert(Natural::from(0u8)) += Natural::from(1u8);
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
        let n_factored = self;
        let n = n_factored.expand();
        if gcd(x.clone(), n.clone()) != Natural::ONE {
            IsPrimitiveRootResult::NonUnit
        } else {
            let phi_n = n_factored.euler_totient();
            let x_mod_n = (&x).rem(&n);
            for p in factor(phi_n.clone()).unwrap().distinct_prime_factors() {
                if x_mod_n.mod_pow_ref(&phi_n / p, &n) == Natural::ONE {
                    return IsPrimitiveRootResult::No;
                }
            }
            IsPrimitiveRootResult::Yes
        }
    }
}

#[derive(Debug)]
struct Factorizer {
    prime_factors: Factored,
    to_factor: Vec<Natural>,
}

impl Factorizer {
    fn new(n: Natural) -> Self {
        Self {
            prime_factors: Factored::one(),
            to_factor: vec![n],
        }
    }

    fn new_with_trial_division(mut n: Natural, max_d: &Natural) -> Self {
        let mut prime_factors = Factored::one();
        let mut d = Natural::TWO;
        while &d * &d <= n && &d <= max_d {
            while &n % &d == Natural::ZERO {
                prime_factors.mul_prime(d.clone());
                n = &n / &d;
            }
            d += Natural::ONE;
        }
        Self {
            prime_factors: prime_factors,
            to_factor: vec![n],
        }
    }

    fn next_to_factor(&mut self) -> Option<Natural> {
        loop {
            match self.to_factor.last() {
                Some(n) => {
                    if n == &Natural::ONE {
                        self.to_factor.pop();
                    } else {
                        return Some(n.clone());
                    }
                }
                None => {
                    return None;
                }
            }
        }
    }

    fn found_factor(&mut self, d: Natural) {
        let n = self.to_factor.pop().unwrap();
        debug_assert_ne!(d, Natural::ZERO);
        debug_assert_ne!(d, n);
        debug_assert_eq!(&n % &d, Natural::ZERO);
        self.to_factor.push(&n / &d);
        self.to_factor.push(d);
    }

    fn found_prime_factor(&mut self, p: Natural) {
        debug_assert!(is_prime(&p));
        let n = self.to_factor.pop().unwrap();
        debug_assert_eq!(&n % &p, Natural::ZERO);
        if p != n {
            self.to_factor.push(n / &p);
        }
        self.prime_factors.mul_prime(p);
    }

    fn complete(self) -> Factored {
        debug_assert!(self.to_factor.is_empty());
        self.prime_factors
    }
}

pub fn factor(n: Natural) -> Option<Factored> {
    if n == Natural::ZERO {
        None
    } else {
        let mut f = Factorizer::new_with_trial_division(n, &Natural::from(10000u32));
        'MAIN: loop {
            match f.next_to_factor() {
                None => {
                    break;
                }
                Some(n) => {
                    debug_assert!(&n >= &Natural::TWO);
                    if n < Natural::from(1000000u32) {
                        // Trial division
                        let mut d = Natural::TWO;
                        while &d * &d <= n {
                            if &n % &d == Natural::ZERO {
                                f.found_prime_factor(d);
                                continue 'MAIN;
                            }
                            d += Natural::ONE;
                        }
                        f.found_prime_factor(n);
                    } else if is_prime(&n) {
                        f.found_prime_factor(n);
                    } else {
                        // Pollard-Rho

                        // g(x) = x^2 + 1
                        let g1 = Polynomial::<Natural>::from_coeffs(vec![
                            Natural::ONE,
                            Natural::ZERO,
                            Natural::ONE,
                        ]);
                        // g(g(x))
                        let g2 = Polynomial::compose(&g1, &g1);

                        'RHO_LOOP: for mut x in (2usize..).map(|n| Natural::from(n)) {
                            let mut y = x.clone();
                            loop {
                                x = g1.evaluate(&x) % &n;
                                y = g2.evaluate(&y) % &n;
                                let d = gcd(Natural::abs_diff(x.clone(), &y), n.clone());
                                if d > Natural::ONE {
                                    debug_assert!(d <= n);
                                    if d == n {
                                        continue 'RHO_LOOP;
                                    } else {
                                        f.found_factor(d);
                                        continue 'MAIN;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        Some(f.complete())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_factor() {
        println!("{:?}", factor(Natural::from(12usize)));
    }

    #[test]
    fn test_euler_totient() {
        assert_eq!(
            factor(Natural::from(12usize)).unwrap().euler_totient(),
            Natural::from(4usize)
        );
    }

    #[test]
    fn test_is_primitive_root() {
        assert_eq!(
            factor(Natural::from(761usize))
                .unwrap()
                .is_primitive_root(&Natural::from(0usize)),
            IsPrimitiveRootResult::NonUnit,
        );
        assert_eq!(
            factor(Natural::from(761usize))
                .unwrap()
                .is_primitive_root(&Natural::from(1usize)),
            IsPrimitiveRootResult::No,
        );
        assert_eq!(
            factor(Natural::from(761usize))
                .unwrap()
                .is_primitive_root(&Natural::from(2usize)),
            IsPrimitiveRootResult::No,
        );
        assert_eq!(
            factor(Natural::from(761usize))
                .unwrap()
                .is_primitive_root(&Natural::from(3usize)),
            IsPrimitiveRootResult::No,
        );
        assert_eq!(
            factor(Natural::from(761usize))
                .unwrap()
                .is_primitive_root(&Natural::from(4usize)),
            IsPrimitiveRootResult::No,
        );
        assert_eq!(
            factor(Natural::from(761usize))
                .unwrap()
                .is_primitive_root(&Natural::from(5usize)),
            IsPrimitiveRootResult::No,
        );
        assert_eq!(
            factor(Natural::from(761usize))
                .unwrap()
                .is_primitive_root(&Natural::from(6usize)),
            IsPrimitiveRootResult::Yes,
        );
    }
}
