use super::functions::*;
use super::*;
use crate::polynomial::polynomial::Polynomial;
use algebraeon_nzq::traits::AbsDiff;
use primes::is_prime;
use std::collections::HashMap;

mod ecm;

#[derive(Debug, Clone)]
pub struct Factored {
    primes: HashMap<Natural, usize>,
}

impl Factored {
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

    fn mul_prime(&mut self, p: Natural) {
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
impl Factored {
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

#[derive(Debug)]
enum Factor {
    Prime(Natural),
    Composite(Natural),
    StrictlyComposite(Natural),
}

#[derive(Debug, Clone)]
enum ToFactorType {
    Composite,
    StrictlyComposite,
}

#[derive(Debug, Clone)]
struct ToFactor {
    n: Natural,
    t: ToFactorType,
}
impl ToFactor {
    fn new_composite(n: Natural) -> Self {
        Self {
            n,
            t: ToFactorType::Composite,
        }
    }

    fn new_strictly_composite(n: Natural) -> Self {
        Self {
            n,
            t: ToFactorType::StrictlyComposite,
        }
    }
}

#[derive(Debug)]
struct Factorizer {
    prime_factors: Factored,
    to_factor: Vec<ToFactor>,
}

impl Factorizer {
    fn new(n: Natural) -> Self {
        debug_assert_ne!(n, Natural::ZERO);
        Self {
            prime_factors: Factored::one(),
            to_factor: vec![ToFactor::new_composite(n)],
        }
    }

    fn partially_factor_by_method(&mut self, algorithm: impl Fn(ToFactor) -> (Vec<Factor>, bool)) {
        let mut to_factor_now = self.to_factor.clone();
        self.to_factor = vec![];
        while !to_factor_now.is_empty() {
            let n = to_factor_now.pop().unwrap();
            #[cfg(debug_assertions)]
            let n_copy = n.n.clone();
            let mut prod = Natural::ONE;
            let (factors, terminate) = algorithm(n);
            for factor in factors {
                match factor {
                    Factor::Prime(p) => {
                        debug_assert_ne!(p, Natural::ONE);
                        debug_assert!(is_prime(&p));
                        prod *= &p;
                        self.prime_factors.mul_prime(p);
                    }
                    Factor::Composite(d) => {
                        debug_assert_ne!(d, Natural::ONE);
                        prod *= &d;
                        if terminate {
                            self.to_factor.push(ToFactor::new_composite(d));
                        } else {
                            to_factor_now.push(ToFactor::new_composite(d))
                        }
                    }
                    Factor::StrictlyComposite(d) => {
                        debug_assert_ne!(d, Natural::ONE);
                        prod *= &d;
                        if terminate {
                            self.to_factor.push(ToFactor::new_strictly_composite(d));
                        } else {
                            to_factor_now.push(ToFactor::new_strictly_composite(d))
                        }
                    }
                }
            }
            #[cfg(debug_assertions)]
            assert_eq!(n_copy, prod);
        }
    }

    fn complete(self) -> Factored {
        assert!(self.to_factor.is_empty());
        self.prime_factors
    }
}

fn trial_division(mut n: Natural, max_d: usize) -> Vec<Factor> {
    let mut factors = vec![];
    let mut d = 2usize;
    while Natural::from(d * d) <= n {
        debug_assert_ne!(Natural::from(d), n);
        if d > max_d {
            factors.push(Factor::Composite(n));
            return factors;
        }
        let d_nat = d.into();
        while &n % &d_nat == Natural::ZERO {
            factors.push(Factor::Prime(d_nat.clone()));
            n = &n / &d_nat;
        }
        if n == Natural::ONE {
            return factors;
        }
        d += 1;
    }
    factors.push(Factor::Prime(n));
    factors
}

fn pollard_rho(n: Natural, mut x: Natural, max_steps: usize) -> Vec<Factor> {
    debug_assert!(!is_prime(&n));

    println!("{}", n);

    // g(x) = x^2 + 1
    let g1 = Polynomial::<Natural>::from_coeffs(vec![Natural::ONE, Natural::ZERO, Natural::ONE]);
    // g(g(x))
    let g2 = Polynomial::compose(&g1, &g1);

    let mut y = x.clone();
    for _ in 0..max_steps {
        x = g1.evaluate(&x) % &n;
        y = g2.evaluate(&y) % &n;
        let d = gcd(Natural::abs_diff(x.clone(), &y), n.clone());
        if d > Natural::ONE {
            debug_assert!(d <= n);
            if d == n {
                return vec![Factor::StrictlyComposite(n)];
            } else {
                return vec![Factor::Composite(n / &d), Factor::Composite(d)];
            }
        }
    }
    return vec![Factor::StrictlyComposite(n)];
}

fn terminate_once_trivial(
    n: ToFactor,
    algorithm: impl Fn(ToFactor) -> Vec<Factor>,
) -> (Vec<Factor>, bool) {
    let factors = algorithm(n);
    let mut prime_count: usize = 0;
    let mut composite_count: usize = 0;
    for factor in &factors {
        match factor {
            Factor::Prime(_) => {
                prime_count += 1;
            }
            Factor::Composite(_) | Factor::StrictlyComposite(_) => {
                composite_count += 1;
            }
        }
    }
    (factors, composite_count == 1 && prime_count == 0)
}

fn exclude_prime_inputs(n: ToFactor, algorithm: impl Fn(Natural) -> Vec<Factor>) -> Vec<Factor> {
    match n.t {
        ToFactorType::Composite => {
            if is_prime(&n.n) {
                return vec![Factor::Prime(n.n)];
            }
        }
        ToFactorType::StrictlyComposite => {}
    }
    algorithm(n.n)
}

pub fn factor(n: Natural) -> Option<Factored> {
    if n == Natural::ZERO {
        None
    } else if n == Natural::ONE {
        Some(Factored::one())
    } else {
        let mut f = Factorizer::new(n);
        // Trial division
        f.partially_factor_by_method(|n| (trial_division(n.n, 1000000000000), true));

        // // Pollard-Rho
        // for x in [2u32, 3, 4] {
        //     f.partially_factor_by_method(|n| {
        //         terminate_once_trivial(n, |n| {
        //             exclude_prime_inputs(n, |n| pollard_rho(n, Natural::from(x), 10000))
        //         })
        //     });
        // }

        // ECM
        f.partially_factor_by_method(|n| {
            (
                exclude_prime_inputs(n, |n| ecm::factor_by_lenstra_elliptic_curve(n)),
                false,
            )
        });

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
