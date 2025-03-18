use std::borrow::Borrow;

use primes::PrimeGenerator;

use super::*;

pub fn factorial(n: impl Borrow<Natural>) -> Natural {
    let mut k = Natural::from(1u8);
    let mut i = Natural::from(1u8);
    while i <= *n.borrow() {
        k *= &i;
        i += Natural::from(1u8);
    }
    k
}

pub fn choose_usize(a: usize, b: usize) -> Natural {
    // I think this is a bit slow
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

pub fn choose(a: impl Borrow<Natural>, b: impl Borrow<Natural>) -> Natural {
    if b.borrow() > a.borrow() {
        Natural::ZERO
    } else {
        let mut t = Natural::ONE;
        let mut m = a.borrow() - b.borrow();
        while &m < a.borrow() {
            m += Natural::ONE;
            t *= &m;
        }
        t / factorial(b.borrow())
    }

    // choose_usize(
    //     nat_to_usize(a.borrow()).unwrap(),
    //     nat_to_usize(b.borrow()).unwrap(),
    // )
}

pub fn pow(x: &Natural, n: &Natural) -> Natural {
    use malachite_base::num::logic::traits::BitIterable;
    if *n == Natural::ZERO {
        Natural::ONE
    } else if *n == Natural::ONE {
        x.clone()
    } else {
        debug_assert!(*n >= Natural::TWO);
        let bits: Vec<_> = n.to_malachite_ref().bits().collect();
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
pub fn nth_root_floor(x: &Natural, n: &Natural) -> Natural {
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

/// Compute the floor of the square root of x
pub fn sqrt_floor(x: &Natural) -> Natural {
    nth_root_floor(x, &Natural::TWO)
}

/// Compute the cil of the nth root of a
pub fn nth_root_ceil(x: &Natural, n: &Natural) -> Natural {
    let mut a = nth_root_floor(x, n);
    if pow(&a, n) != *x {
        a += Natural::ONE;
    }
    a
}

/// Compute the ceil of the square root of x
pub fn sqrt_ceil(x: &Natural) -> Natural {
    nth_root_ceil(x, &Natural::TWO)
}

/// Return the number of bits needed to store n i.e. ceil(log2(n)) for all non-zero n
pub fn bitcount(n: impl Borrow<Natural>) -> usize {
    use malachite_base::num::logic::traits::BitIterable;
    n.borrow().to_malachite_ref().bits().len()
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
        // The largest power n can possibly be is when n is a power of 2 and n = 2^{bitcount(n)+1}
        // So we only need to check n isn't a kth power up to k = bitcount(n)+1
        // We also only need to check for prime k
        let max_k = Natural::from(bitcount(n) + 1);
        for k in PrimeGenerator::new() {
            if k > max_k {
                return IsPowerTestResult::No;
            } else {
                let a = nth_root_floor(n, &k);
                if *n == pow(&a, &k) {
                    return IsPowerTestResult::Power(a, k);
                }
            }
        }
        unreachable!()
    }
}

pub fn gcd(mut x: Natural, mut y: Natural) -> Natural {
    while y != Natural::ZERO {
        let r = x.rem(&y);
        (x, y) = (y, r)
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_factorial() {
        assert_eq!(factorial(&Natural::from(0usize)), Natural::from(1usize));
        assert_eq!(factorial(&Natural::from(1usize)), Natural::from(1usize));
        assert_eq!(factorial(&Natural::from(2usize)), Natural::from(2usize));
        assert_eq!(factorial(&Natural::from(3usize)), Natural::from(6usize));
    }

    #[test]
    fn test_choose_usize() {
        assert_eq!(choose_usize(0, 0), Natural::from(1u32));

        assert_eq!(choose_usize(1, 0), Natural::from(1u32));
        assert_eq!(choose_usize(1, 1), Natural::from(1u32));

        assert_eq!(choose_usize(2, 0), Natural::from(1u32));
        assert_eq!(choose_usize(2, 1), Natural::from(2u32));
        assert_eq!(choose_usize(2, 2), Natural::from(1u32));

        assert_eq!(choose_usize(3, 0), Natural::from(1u32));
        assert_eq!(choose_usize(3, 1), Natural::from(3u32));
        assert_eq!(choose_usize(3, 2), Natural::from(3u32));
        assert_eq!(choose_usize(3, 3), Natural::from(1u32));

        assert_eq!(choose_usize(4, 0), Natural::from(1u32));
        assert_eq!(choose_usize(4, 1), Natural::from(4u32));
        assert_eq!(choose_usize(4, 2), Natural::from(6u32));
        assert_eq!(choose_usize(4, 3), Natural::from(4u32));
        assert_eq!(choose_usize(4, 4), Natural::from(1u32));

        assert_eq!(choose_usize(3, 4), Natural::from(0u32));
    }

    #[test]
    fn test_choose() {
        assert_eq!(
            choose(Natural::from(0usize), Natural::from(0usize)),
            Natural::from(1u32)
        );

        assert_eq!(
            choose(Natural::from(1usize), Natural::from(0usize)),
            Natural::from(1u32)
        );
        assert_eq!(
            choose(Natural::from(1usize), Natural::from(1usize)),
            Natural::from(1u32)
        );

        assert_eq!(
            choose(Natural::from(2usize), Natural::from(0usize)),
            Natural::from(1u32)
        );
        assert_eq!(
            choose(Natural::from(2usize), Natural::from(1usize)),
            Natural::from(2u32)
        );
        assert_eq!(
            choose(Natural::from(2usize), Natural::from(2usize)),
            Natural::from(1u32)
        );

        assert_eq!(
            choose(Natural::from(3usize), Natural::from(0usize)),
            Natural::from(1u32)
        );
        assert_eq!(
            choose(Natural::from(3usize), Natural::from(1usize)),
            Natural::from(3u32)
        );
        assert_eq!(
            choose(Natural::from(3usize), Natural::from(2usize)),
            Natural::from(3u32)
        );
        assert_eq!(
            choose(Natural::from(3usize), Natural::from(3usize)),
            Natural::from(1u32)
        );

        assert_eq!(
            choose(Natural::from(4usize), Natural::from(0usize)),
            Natural::from(1u32)
        );
        assert_eq!(
            choose(Natural::from(4usize), Natural::from(1usize)),
            Natural::from(4u32)
        );
        assert_eq!(
            choose(Natural::from(4usize), Natural::from(2usize)),
            Natural::from(6u32)
        );
        assert_eq!(
            choose(Natural::from(4usize), Natural::from(3usize)),
            Natural::from(4u32)
        );
        assert_eq!(
            choose(Natural::from(4usize), Natural::from(4usize)),
            Natural::from(1u32)
        );

        assert_eq!(
            choose(Natural::from(3usize), Natural::from(4usize)),
            Natural::from(0u32)
        );
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
                let r = nth_root_floor(&x, &n);
                println!("{}th root of {} is {}", &n, &x, &r);
                assert!(pow(&r, &n) <= x);
                assert!(x == Natural::ZERO || pow(&(r + Natural::ONE), &n) > x);
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
    fn test_gcd() {
        assert_eq!(
            gcd(Natural::from(0usize), Natural::from(0usize)),
            Natural::from(0usize)
        );
        assert_eq!(
            gcd(Natural::from(1usize), Natural::from(0usize)),
            Natural::from(1usize)
        );
        assert_eq!(
            gcd(Natural::from(0usize), Natural::from(1usize)),
            Natural::from(1usize)
        );
        assert_eq!(
            gcd(Natural::from(2usize), Natural::from(0usize)),
            Natural::from(2usize)
        );
        assert_eq!(
            gcd(Natural::from(0usize), Natural::from(2usize)),
            Natural::from(2usize)
        );

        assert_eq!(
            gcd(Natural::from(12usize), Natural::from(15usize)),
            Natural::from(3usize)
        );
    }
}
