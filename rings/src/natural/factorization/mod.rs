use super::functions::*;
use super::*;
use crate::polynomial::{Polynomial, SemiRingToPolynomialSemiRingSignature};
use algebraeon_nzq::traits::AbsDiff;
pub use factored::*;
use primes::is_prime;

pub mod ecm;
pub mod factored;
pub mod primes;

#[derive(Debug, Clone)]
pub enum Factor {
    Prime(Natural),
    Composite(Natural),
    StrictlyComposite(Natural),
}

fn trivial_factor(n: ToFactor) -> Vec<Factor> {
    match n.t {
        ToFactorType::Composite => {
            if is_prime(&n.n) {
                vec![Factor::Prime(n.n)]
            } else {
                vec![Factor::Composite(n.n)]
            }
        }
        ToFactorType::StrictlyComposite => vec![Factor::StrictlyComposite(n.n)],
    }
}

pub fn trial_division(mut n: Natural, max_d: usize) -> Vec<Factor> {
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

pub fn pollard_rho(n: Natural, mut x: Natural, max_steps: usize) -> Vec<Factor> {
    debug_assert!(!is_prime(&n));

    let nat_polys = Natural::structure().into_polynomial_semiring();

    // g(x) = x^2 + 1
    let g1 = Polynomial::<Natural>::from_coeffs(vec![Natural::ONE, Natural::ZERO, Natural::ONE]);
    // g(g(x))
    let g2 = nat_polys.compose(&g1, &g1);

    let mut y = x.clone();
    for _ in 0..max_steps {
        x = nat_polys.evaluate(&g1, &x) % &n;
        y = nat_polys.evaluate(&g2, &y) % &n;
        let d = gcd(Natural::abs_diff(x.clone(), &y), n.clone());
        if d > Natural::ONE {
            debug_assert!(d <= n);
            return if d == n {
                vec![Factor::StrictlyComposite(n)]
            } else {
                vec![Factor::Composite(n / &d), Factor::Composite(d)]
            };
        }
    }
    vec![Factor::StrictlyComposite(n)]
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
    prime_factors: Vec<(Natural, Natural)>,
    to_factor: Vec<ToFactor>,
}

impl Factorizer {
    fn new(n: Natural) -> Self {
        debug_assert_ne!(n, Natural::ZERO);
        Self {
            prime_factors: Natural::structure().factorizations().new_trivial(),
            to_factor: vec![ToFactor::new_composite(n)],
        }
    }

    fn partially_factor_by_method(&mut self, algorithm: impl Fn(ToFactor) -> (Vec<Factor>, bool)) {
        let factorizations = Natural::structure().factorizations();
        let mut to_factor_now = self.to_factor.clone();
        self.to_factor = vec![];
        #[allow(clippy::manual_while_let_some)]
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
                        factorizations.mul_prime(&mut self.prime_factors, p);
                    }
                    Factor::Composite(d) => {
                        debug_assert_ne!(d, Natural::ONE);
                        prod *= &d;
                        if terminate {
                            self.to_factor.push(ToFactor::new_composite(d));
                        } else {
                            to_factor_now.push(ToFactor::new_composite(d));
                        }
                    }
                    Factor::StrictlyComposite(d) => {
                        debug_assert_ne!(d, Natural::ONE);
                        prod *= &d;
                        if terminate {
                            self.to_factor.push(ToFactor::new_strictly_composite(d));
                        } else {
                            to_factor_now.push(ToFactor::new_strictly_composite(d));
                        }
                    }
                }
            }
            #[cfg(debug_assertions)]
            assert_eq!(n_copy, prod);
        }
    }

    fn complete(self) -> Vec<(Natural, Natural)> {
        assert!(self.to_factor.is_empty());
        self.prime_factors
    }
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
    debug_assert!(!is_prime(&n.n));
    algorithm(n.n)
}

pub(super) fn factor(n: Natural) -> Option<Vec<(Natural, Natural)>> {
    if n == Natural::ZERO {
        None
    } else if n == Natural::ONE {
        Some(Natural::structure().factorizations().new_trivial())
    } else {
        let mut f = Factorizer::new(n);
        // Trial division
        f.partially_factor_by_method(|n| (trial_division(n.n, 100_000), true));

        // Pollard-Rho
        for x in [2u32, 3, 4] {
            f.partially_factor_by_method(|n| {
                terminate_once_trivial(n, |n| {
                    if n.n.bitcount() < 100 {
                        exclude_prime_inputs(n, |n| pollard_rho(n, Natural::from(x), 10000))
                    } else {
                        trivial_factor(n)
                    }
                })
            });
        }

        // ECM
        f.partially_factor_by_method(|n| {
            (
                exclude_prime_inputs(n, |n| {
                    let mut rng = algebraeon_nzq::Rng::new(0);
                    for fith_target_factor_digits in vec![2, 3, 2, 4, 2, 3].into_iter().cycle() {
                        if let Ok(d) = ecm::ecm_one_factor_target_digits(
                            &n,
                            fith_target_factor_digits,
                            &mut rng,
                        ) {
                            return vec![Factor::Composite(n / &d), Factor::Composite(d)];
                        }
                    }
                    unreachable!()
                }),
                false,
            )
        });

        Some(f.complete())
    }
}

#[cfg(test)]
mod tests {
    use crate::natural::factorization::factored::IsPrimitiveRootResult;

    use super::*;

    #[test]
    fn test_factor() {
        println!("{:?}", factor(Natural::from(12usize)));
    }

    #[test]
    fn test_euler_totient() {
        assert_eq!(
            Natural::structure()
                .factorizations()
                .euler_totient(&factor(Natural::from(12usize)).unwrap()),
            Natural::from(4usize)
        );
    }

    #[test]
    fn test_is_primitive_root() {
        let factorizations = Natural::structure().factorizations();
        assert_eq!(
            factorizations.is_primitive_root(
                &Natural::from(0usize),
                &factor(Natural::from(761usize)).unwrap(),
            ),
            IsPrimitiveRootResult::NonUnit,
        );
        assert_eq!(
            factorizations.is_primitive_root(
                &Natural::from(1usize),
                &factor(Natural::from(761usize)).unwrap(),
            ),
            IsPrimitiveRootResult::No,
        );
        assert_eq!(
            factorizations.is_primitive_root(
                &Natural::from(2usize),
                &factor(Natural::from(761usize)).unwrap(),
            ),
            IsPrimitiveRootResult::No,
        );
        assert_eq!(
            factorizations.is_primitive_root(
                &Natural::from(3usize),
                &factor(Natural::from(761usize)).unwrap(),
            ),
            IsPrimitiveRootResult::No,
        );
        assert_eq!(
            factorizations.is_primitive_root(
                &Natural::from(4usize),
                &factor(Natural::from(761usize)).unwrap(),
            ),
            IsPrimitiveRootResult::No,
        );
        assert_eq!(
            factorizations.is_primitive_root(
                &Natural::from(5usize),
                &factor(Natural::from(761usize)).unwrap(),
            ),
            IsPrimitiveRootResult::No,
        );
        assert_eq!(
            factorizations.is_primitive_root(
                &Natural::from(6usize),
                &factor(Natural::from(761usize)).unwrap(),
            ),
            IsPrimitiveRootResult::Yes,
        );
    }
}
