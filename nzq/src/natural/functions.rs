use super::*;

impl Natural {
    /// Factorial
    /// ```
    /// use algebraeon_nzq::Natural;
    /// assert_eq!(
    ///     Natural::from(0u32).factorial(),
    ///     Natural::from(1u32),
    /// );    
    /// assert_eq!(
    ///     Natural::from(1u32).factorial(),
    ///     Natural::from(1u32),
    /// );
    /// assert_eq!(
    ///     Natural::from(2u32).factorial(),
    ///     Natural::from(2u32),
    /// );
    /// assert_eq!(
    ///     Natural::from(5u32).factorial(),
    ///     Natural::from(120u32),
    /// );
    /// ```
    pub fn factorial(&self) -> Self {
        let mut i = Natural::ZERO;
        let mut f = Natural::ONE;
        while &i < self {
            i += Natural::ONE;
            f *= &i;
        }
        f
    }

    pub fn pow(&self, n: &Natural) -> Natural {
        let x = self;
        if *n == Natural::ZERO {
            Natural::ONE
        } else if *n == Natural::ONE {
            x.clone()
        } else {
            debug_assert!(*n >= Natural::TWO);
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
    pub fn nth_root_floor(&self, n: &Natural) -> Natural {
        let x = self;
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
                if m.pow(&n) <= *x {
                    a = m;
                } else {
                    b = m;
                }
            }
            a
        }
    }

    /// Compute the floor of the square root of x
    pub fn sqrt_floor(&self) -> Natural {
        self.nth_root_floor(&Natural::TWO)
    }

    /// Compute the cil of the nth root of a
    pub fn nth_root_ceil(&self, n: &Natural) -> Natural {
        let x = self;
        let mut a = x.nth_root_floor(n);
        if a.pow(n) != *x {
            a += Natural::ONE;
        }
        a
    }

    /// Compute the ceil of the square root of x
    pub fn sqrt_ceil(&self) -> Natural {
        self.nth_root_ceil(&Natural::TWO)
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
        t / b.borrow().factorial()
    }
}

pub fn gcd(mut x: Natural, mut y: Natural) -> Natural {
    while y != Natural::ZERO {
        let r = x % &y;
        (x, y) = (y, r)
    }
    x
}

pub fn lcm(x: Natural, y: Natural) -> Natural {
    (&x * &y) / gcd(x, y)
}

// fn stirling_partition_number(n: &Natural, x: &Natural) -> Natural {
//     todo!()
// }

// fn stirling_cycle_number(n: &Natural, x: &Natural) -> Natural {
//     todo!()
// }

#[cfg(test)]
mod tests {
    use super::*;

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
            Natural::from(0usize).pow(&Natural::from(0usize)),
            Natural::from(1usize)
        );
        assert_eq!(
            Natural::from(0usize).pow(&Natural::from(1usize)),
            Natural::from(0usize)
        );
        assert_eq!(
            Natural::from(0usize).pow(&Natural::from(2usize)),
            Natural::from(0usize)
        );
        assert_eq!(
            Natural::from(0usize).pow(&Natural::from(3usize)),
            Natural::from(0usize)
        );

        assert_eq!(
            Natural::from(1usize).pow(&Natural::from(0usize)),
            Natural::from(1usize)
        );
        assert_eq!(
            Natural::from(1usize).pow(&Natural::from(1usize)),
            Natural::from(1usize)
        );
        assert_eq!(
            Natural::from(1usize).pow(&Natural::from(2usize)),
            Natural::from(1usize)
        );
        assert_eq!(
            Natural::from(1usize).pow(&Natural::from(3usize)),
            Natural::from(1usize)
        );

        assert_eq!(
            Natural::from(2usize).pow(&Natural::from(0usize)),
            Natural::from(1usize)
        );
        assert_eq!(
            Natural::from(2usize).pow(&Natural::from(1usize)),
            Natural::from(2usize)
        );
        assert_eq!(
            Natural::from(2usize).pow(&Natural::from(2usize)),
            Natural::from(4usize)
        );
        assert_eq!(
            Natural::from(2usize).pow(&Natural::from(3usize)),
            Natural::from(8usize)
        );

        assert_eq!(
            Natural::from(3usize).pow(&Natural::from(0usize)),
            Natural::from(1usize)
        );
        assert_eq!(
            Natural::from(3usize).pow(&Natural::from(1usize)),
            Natural::from(3usize)
        );
        assert_eq!(
            Natural::from(3usize).pow(&Natural::from(2usize)),
            Natural::from(9usize)
        );
        assert_eq!(
            Natural::from(3usize).pow(&Natural::from(3usize)),
            Natural::from(27usize)
        );
    }

    #[test]
    fn test_nth_root() {
        for n in 1..10usize {
            for x in 0..10usize {
                let x = Natural::from(x);
                let n = Natural::from(n);
                let r = x.nth_root_floor(&n);
                println!("{}th root of {} is {}", &n, &x, &r);
                assert!(r.pow(&n) <= x);
                assert!(x == Natural::ZERO || (r + Natural::ONE).pow(&n) > x);
            }
        }
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
