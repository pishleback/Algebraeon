use factor::Factored;
use malachite_base::num::arithmetic::traits::PowerOf2;

use malachite_base::num::arithmetic::traits::DivMod;
use malachite_base::num::arithmetic::traits::{Mod, ModPow};
use malachite_base::num::logic::traits::BitIterable;

use super::functions::*;
use super::*;

#[derive(Debug)]
pub struct PrimeGenerator {
    n: Natural,
    primes: Vec<Natural>,
}

impl PrimeGenerator {
    pub fn new() -> Self {
        Self {
            n: Natural::from(2u8),
            primes: vec![],
        }
    }
}

impl Iterator for PrimeGenerator {
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PrimalityTestResult {
    Zero,
    One,
    Prime,
    Composite,
}

pub fn naive_primality_test(n: &Natural) -> PrimalityTestResult {
    if *n == Natural::ZERO {
        PrimalityTestResult::Zero
    } else if *n == Natural::ONE {
        PrimalityTestResult::One
    } else {
        let mut prime_gen = PrimeGenerator::new();
        loop {
            let p = prime_gen.next().unwrap();
            match p.cmp(n) {
                std::cmp::Ordering::Less => {
                    continue;
                }
                std::cmp::Ordering::Equal => {
                    return PrimalityTestResult::Prime;
                }
                std::cmp::Ordering::Greater => {
                    return PrimalityTestResult::Composite;
                }
            }
        }
    }
}

// https://cr.yp.to/papers/aks.pdf
pub fn aks_primality_test(n: &Natural) -> PrimalityTestResult {
    match is_power_test(n) {
        IsPowerTestResult::Zero => PrimalityTestResult::Zero,
        IsPowerTestResult::One => PrimalityTestResult::One,
        IsPowerTestResult::Power(_, _) => PrimalityTestResult::Composite,
        IsPowerTestResult::No => {
            // println!("n = {:?}", n);

            // Select r0 >= 3
            // Use r0 ~ 0.01*log2(n)^2
            let r0 = {
                let aprox_log2_n = bitcount(n) - 1;
                let r0 = (&aprox_log2_n * &aprox_log2_n) / 100;
                if r0 < 3 {
                    3
                } else {
                    r0
                }
            };
            // println!("r0 = {:?}", r0);

            // Find the smallest prime r >= r0 such that n is a primitive root modulo r
            let mut prime_gen = PrimeGenerator::new();
            let mut r;
            loop {
                r = prime_gen.next().unwrap();
                if r < r0 {
                    continue;
                }
                match Factored::from_prime_unchecked(r.clone()).is_primitive_root(n) {
                    factor::IsPrimitiveRootResult::NonUnit => {
                        // n is divisible by r
                        match *n == r {
                            true => {
                                return PrimalityTestResult::Prime;
                            }
                            false => {
                                return PrimalityTestResult::Composite;
                            }
                        }
                    }
                    factor::IsPrimitiveRootResult::No => {}
                    factor::IsPrimitiveRootResult::Yes => {
                        break;
                    }
                }
            }
            // println!("r = {:?}", r);

            // Select d between 0 and phi(r)-1
            // Use d ~ 0.5*phi(r)
            let phi_r = &r - Natural::ONE;
            let d = &phi_r / Natural::TWO;
            // println!("d = {:?}", d);

            // Select i between 0 and d
            // Use i ~ 0.475*phi(r)
            let i = {
                let i = (Natural::from(475usize) * &phi_r) / Natural::from(1000usize);
                if i > d {
                    d.clone()
                } else {
                    i
                }
            };
            // println!("i = {:?}", i);

            // Select j between 0 and phi(r) - 1 - d
            // Use j ~ 0.475*phi(r)
            let j = {
                let j = (Natural::from(475usize) * &phi_r) / Natural::from(1000usize);
                let max = &phi_r - Natural::ONE - &d;
                if j > max {
                    max
                } else {
                    j
                }
            };
            // println!("j = {:?}", j);

            // Select s such that (2s choose i) * (d choose i) * (2s - i choose j) * (phi(r)-1-d choose j) >= n^ceil(sqrt(phi(r)/3))
            // Use the smallest such s found via binary search

            let s = {
                let rhs = pow(
                    n,
                    &nth_root_ceil(&(&phi_r / Natural::from(3usize)), &Natural::TWO),
                );
                // println!("rhs = {:?}", rhs);

                let lhs = |s: &Natural| -> Natural {
                    let two_s = Natural::TWO * s;
                    if i > two_s {
                        Natural::ZERO
                    } else {
                        choose(&two_s, &i)
                            * choose(&d, &i)
                            * choose(&two_s - &i, &j)
                            * choose(&phi_r - Natural::ONE - &d, &j)
                    }
                };

                let mut step_pow = 0u64;
                let mut s = Natural::ZERO;
                loop {
                    // Grow s until lhs(s) >= rhs
                    s += Natural::power_of_2(step_pow);
                    match lhs(&s).cmp(&rhs) {
                        std::cmp::Ordering::Less => {
                            step_pow += 1;
                            continue;
                        }
                        std::cmp::Ordering::Equal => {
                            break;
                        }
                        std::cmp::Ordering::Greater => {
                            step_pow -= 1;
                        }
                    }
                    // Shrink s until we find the smallest s such that lhs(s) >= rhs
                    loop {
                        match lhs(&(&s - Natural::power_of_2(step_pow))).cmp(&rhs) {
                            std::cmp::Ordering::Less => {
                                if step_pow == 0 {
                                    break;
                                }
                                step_pow -= 1;
                                continue;
                            }
                            std::cmp::Ordering::Equal => {
                                break;
                            }
                            std::cmp::Ordering::Greater => {
                                s -= Natural::power_of_2(step_pow);
                                if step_pow == 0 {
                                    break;
                                }
                                step_pow -= 1;
                                continue;
                            }
                        }
                    }
                    break;
                }
                #[cfg(debug_assertions)]
                {
                    assert!(lhs(&s) >= rhs);
                    if s > 1 {
                        assert!(lhs(&(&s - Natural::ONE)) < rhs);
                    }
                }
                s
            };
            // println!("s = {:?}", s);

            // Let S = {2, 3, ..., s, s + 1}
            let s_set = {
                let mut s_set = vec![];
                let mut b = Natural::TWO;
                while s_set.len() < s {
                    s_set.push(b.clone());
                    b += Natural::ONE;
                }
                s_set
            };
            // println!("S = {:?}", s_set);

            // For all b in S check whether gcd(n,b)=1. If not it is easy to tell if n is prime.
            for b in &s_set {
                let g = gcd(n.clone(), b.clone());
                if g != Natural::ONE {
                    match g == *n {
                        true => {
                            // b and thus also n=g is small, so we can do a naive test
                            return naive_primality_test(n);
                        }
                        false => {
                            return PrimalityTestResult::Composite;
                        }
                    }
                }
            }

            // For all b b' in S check whether gcd(n, bb'-1)=1
            for idx in 0..s_set.len() {
                for jdx in idx..s_set.len() {
                    let bi = &s_set[idx];
                    let bj = &s_set[jdx];
                    let g = gcd(n.clone(), bi * bj - Natural::ONE);
                    if g != Natural::ONE {
                        match g == *n {
                            true => {
                                // bi*bj-1 and thus also n=g is small, so we can do a naive test
                                return naive_primality_test(n);
                            }
                            false => {
                                return PrimalityTestResult::Composite;
                            }
                        }
                    }
                }
            }

            // For all distinct b b' in S check whether gcd(n, b - b')=1
            for idx in 0..s_set.len() {
                for jdx in (idx + 1)..s_set.len() {
                    let bi = &s_set[idx];
                    let bj = &s_set[jdx];
                    let g = gcd(n.clone(), bj - bi);
                    if g != Natural::ONE {
                        match g == *n {
                            true => {
                                // bj-bi and thus also n=g is small, so we can do a naive test
                                return naive_primality_test(n);
                            }
                            false => {
                                return PrimalityTestResult::Composite;
                            }
                        }
                    }
                }
            }

            // For all b in S check whether b^(n-1)=1 mod n
            for b in &s_set {
                if b.mod_op(n).mod_pow(n - Natural::ONE, n) != Natural::ONE {
                    return PrimalityTestResult::Composite;
                }
            }

            // For all b in S check whether (x+b)^n = x^n+b mod x^r-1, n. If not for some b then n is composite
            // In order to make use of fast integer algorithms, we encode polynomials as big integers thought of as r-length arrays of numbers each at most n^2
            // Thus each block of a polynomial bigint shall have length coeff_size := bitcount(r*n^2)

            let r_usize = nat_to_usize(&r).unwrap();
            let coeff_size = bitcount(n * n * &r);
            let coeff_mask = (Natural::ONE << coeff_size) - Natural::ONE;

            // Convert between polynomials and big integers by x <-> 2^k
            let zero_poly = || -> Vec<Natural> { vec![Natural::ZERO; r_usize] };
            let one_poly = || -> Vec<Natural> {
                let mut p = zero_poly();
                p[0] = Natural::ONE;
                p
            };

            let poly_to_bigint = |p: &Vec<Natural>| -> Natural {
                let mut p_bigint = Natural::ZERO;
                for (i, c) in p.iter().enumerate() {
                    p_bigint += c << (i * coeff_size);
                }
                p_bigint
            };
            #[cfg(debug_assertions)]
            let bigint_to_poly = |mut p_bigint: Natural| -> Vec<Natural> {
                let mut p = zero_poly();
                let coeff_modulus = Natural::ONE << coeff_size;
                let mut coeff;
                for i in 0..r_usize {
                    (p_bigint, coeff) = p_bigint.div_mod(&coeff_modulus);
                    p[i] = coeff.mod_op(n);
                }
                p
            };

            #[cfg(debug_assertions)]
            let mul_polys_naive = |a: &Vec<Natural>, b: &Vec<Natural>| -> Vec<Natural> {
                let mut s = zero_poly();
                for i in 0..r_usize {
                    for j in 0..r_usize {
                        s[(i + j) % r_usize] = (&s[(i + j) % r_usize] + &a[i] * &b[j]).mod_op(n);
                    }
                }
                s
            };

            // Reduce the coefficients of a modulo n
            let reduce_coeffs = |a: &Natural| -> Natural {
                let mut b = Natural::ZERO;
                for i in 0..r_usize {
                    let coeff =
                        ((&*a & (&coeff_mask << (i * coeff_size))) >> (i * coeff_size)).mod_op(n);
                    b |= coeff << (coeff_size * i);
                }
                b
            };

            let mul_polys = |a: &Natural, b: &Natural| -> Natural {
                // Use fast integer multiplication to compute the polynomial product
                let m = (Natural::ONE << (coeff_size * r_usize)) - Natural::ONE;
                let s = (a * b).mod_op(m);
                // println!("s = {:?}", s);
                let s = reduce_coeffs(&s);
                // println!("s = {:?}", s);
                #[cfg(debug_assertions)]
                {
                    // println!("a = {:?}", a);
                    // println!("b = {:?}", b);

                    // println!("a = {:?}", bigint_to_poly(a.clone()));
                    // println!("b = {:?}", bigint_to_poly(b.clone()));

                    assert_eq!(
                        s,
                        poly_to_bigint(&mul_polys_naive(
                            &bigint_to_poly(a.clone()),
                            &bigint_to_poly(b.clone())
                        ))
                    );
                }
                s
            };

            let square_poly = |poly: &Natural| -> Natural { mul_polys(poly, poly) };

            let nthpow_poly = |mut poly: Natural, k: &Natural| -> Natural {
                let mut s = poly_to_bigint(&one_poly());
                for k_bit in k.bits() {
                    if k_bit {
                        s = mul_polys(&s, &poly);
                    }
                    poly = square_poly(&poly);
                }
                s
            };

            for b in &s_set {
                // println!("b = {:?}", b);
                // (x + b)^n
                let lhs = {
                    let mut poly = zero_poly();
                    poly[0] = b.clone(); //this is ok since r >= 3
                    poly[1] = Natural::ONE;
                    nthpow_poly(poly_to_bigint(&poly), n)
                };

                // x^n + b
                let rhs = {
                    let mut poly = zero_poly();
                    poly[0] = b.clone(); //this is ok since r >= 3
                    poly[nat_to_usize(&n.mod_op(&r)).unwrap()] += Natural::ONE;
                    poly_to_bigint(&poly)
                };

                if lhs != rhs {
                    return PrimalityTestResult::Composite;
                }
            }

            // Otherwise n is prime
            PrimalityTestResult::Prime
        }
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use super::*;

    #[test]
    fn test_aks_primality_test() {
        assert_eq!(
            aks_primality_test(&Natural::from_str("0").unwrap()),
            PrimalityTestResult::Zero
        );
        assert_eq!(
            aks_primality_test(&Natural::from_str("1").unwrap()),
            PrimalityTestResult::One
        );
        assert_eq!(
            aks_primality_test(&Natural::from_str("2").unwrap()),
            PrimalityTestResult::Prime
        );
        assert_eq!(
            aks_primality_test(&Natural::from_str("3").unwrap()),
            PrimalityTestResult::Prime
        );
        assert_eq!(
            aks_primality_test(&Natural::from_str("4").unwrap()),
            PrimalityTestResult::Composite
        );
        assert_eq!(
            aks_primality_test(&Natural::from_str("17239").unwrap()),
            PrimalityTestResult::Prime
        );
        assert_eq!(
            aks_primality_test(&Natural::from_str("126487").unwrap()),
            PrimalityTestResult::Prime
        );
        assert_eq!(
            aks_primality_test(&Natural::from_str("198741").unwrap()),
            PrimalityTestResult::Composite
        );
    }
}
