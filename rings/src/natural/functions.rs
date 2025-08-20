use super::*;
use algebraeon_nzq::primes;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IsPowerTestResult {
    Zero,
    One,
    Power(Natural, usize),
    No,
}

/// Returns n=a^p where a is as large as possible and p is prime
pub fn is_power_test(n: &Natural) -> IsPowerTestResult {
    if *n == Natural::ZERO {
        IsPowerTestResult::Zero
    } else if *n == Natural::ONE {
        IsPowerTestResult::One
    } else {
        // The largest power n can possibly be is when n is a power of 2 and n = 2^{bitcount(n)+1}
        // So we only need to check n isn't a kth power up to k = bitcount(n)+1
        // We also only need to check for prime k
        let max_k = n.bitcount() + 1;
        for k in primes().take_while(|&k| k <= max_k) {
            let a = n.nth_root_floor(&k.into());
            if *n == a.pow(&k.into()) {
                return IsPowerTestResult::Power(a, k);
            }
        }
        IsPowerTestResult::No
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
            IsPowerTestResult::Power(Natural::from(2usize), 2)
        );
        assert_eq!(is_power_test(&Natural::from(5usize)), IsPowerTestResult::No);
        assert_eq!(is_power_test(&Natural::from(6usize)), IsPowerTestResult::No);
        assert_eq!(is_power_test(&Natural::from(7usize)), IsPowerTestResult::No);
        assert_eq!(
            is_power_test(&Natural::from(8usize)),
            IsPowerTestResult::Power(Natural::from(2usize), 3)
        );
        assert_eq!(
            is_power_test(&Natural::from(9usize)),
            IsPowerTestResult::Power(Natural::from(3usize), 2)
        );
        assert_eq!(
            is_power_test(&Natural::from(10usize)),
            IsPowerTestResult::No
        );

        println!("{:?}", is_power_test(&Natural::from(0usize)));
    }
}
