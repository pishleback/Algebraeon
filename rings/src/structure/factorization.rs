use algebraeon_nzq::*;
use algebraeon_sets::structure::*;
use std::fmt::Debug;

pub trait FactoredSignature: SetSignature {
    /// A type used to hold objects before they are factored
    type Object: Clone + Debug;
    /// A type used to hold the prime objects.
    type PrimeObject: Clone + Debug;

    /// Same as `new_powers` but don't validate the input
    fn new_powers_unchecked(&self, factor_powers: Vec<(Self::PrimeObject, Natural)>) -> Self::Set;

    /// Construct a new factorization object given a list of powers of primes
    fn new_powers(&self, factor_powers: Vec<(Self::PrimeObject, Natural)>) -> Self::Set {
        #[cfg(debug_assertions)]
        {
            let mut ps: Vec<&Self::PrimeObject> = vec![];
            for (p, _) in &factor_powers {
                assert!(self.try_object_is_prime(p) != Some(false));
                assert!(!ps.iter().any(|q| self.prime_object_equivalent(p, q)));
                ps.push(p);
            }
        }
        self.new_powers_unchecked(
            factor_powers
                .into_iter()
                .filter(|(_, k)| k != &Natural::ZERO)
                .collect(),
        )
    }

    /// Return the trivial factorization
    fn new_trivial(&self) -> Self::Set {
        self.new_powers(vec![])
    }

    /// Construct a factorization from a prime object
    fn new_prime(&self, prime: Self::PrimeObject) -> Self::Set {
        debug_assert!(self.try_object_is_prime(&prime) != Some(false));
        self.new_powers(vec![(prime, Natural::ONE)])
    }

    /// Same as `to_powers` but don't validate the input
    fn to_powers_unchecked<'a>(
        &self,
        a: &'a Self::Set,
    ) -> Vec<(&'a Self::PrimeObject, &'a Natural)>;

    /// Given a factorization return a list of it's primes and their powers
    fn to_powers<'a>(&self, a: &'a Self::Set) -> Vec<(&'a Self::PrimeObject, &'a Natural)> {
        debug_assert!(self.is_element(a));
        self.to_powers_unchecked(a)
    }

    /// Same as `into_powers` but don't validate the input
    fn into_powers_unchecked(&self, a: Self::Set) -> Vec<(Self::PrimeObject, Natural)>;

    /// Consume a factorization and return a list of it's primes and their powers
    fn into_powers(&self, a: Self::Set) -> Vec<(Self::PrimeObject, Natural)> {
        debug_assert!(self.is_element(&a));
        self.into_powers_unchecked(a)
    }

    /// Given a factorization return a list of it's prime factors with duplicate entries according to the multiplicity
    fn to_primes<'a>(&self, a: &'a Self::Set) -> Vec<&'a Self::PrimeObject> {
        debug_assert!(self.is_element(a));
        let mut factors = vec![];
        for (p, k) in self.to_powers(a) {
            let k: usize = k.try_into().unwrap();
            for _ in 0..k {
                factors.push(p);
            }
        }
        factors
    }

    /// Consume a factorization and return a list of it's prime factors with duplicate entries according to the multiplicity
    fn into_primes(&self, a: Self::Set) -> Vec<Self::PrimeObject> {
        debug_assert!(self.is_element(&a));
        let mut factors = vec![];
        for (p, k) in self.into_powers(a) {
            let k: usize = k.try_into().unwrap();
            for _ in 0..k {
                factors.push(p.clone());
            }
        }
        factors
    }

    /// Given a factorization return a list of it's prime factors without multiplicity i.e. without repeats
    fn to_prime_support<'a>(&self, a: &'a Self::Set) -> Vec<&'a Self::PrimeObject> {
        debug_assert!(self.is_element(a));
        self.to_powers(a).into_iter().map(|(p, _)| p).collect()
    }

    /// Consume a factorization and return a list of it's prime factors without multiplicity i.e. without repeats
    fn into_prime_support(self, a: Self::Set) -> Vec<Self::PrimeObject> {
        debug_assert!(self.is_element(&a));
        self.into_powers(a).into_iter().map(|(p, _)| p).collect()
    }

    /// Return true if the factorization is trivial
    fn is_trivial(&self, a: &Self::Set) -> bool {
        debug_assert!(self.is_element(a));
        self.to_powers(a).into_iter().next().is_none()
    }

    /// Return true if the factorization is a single prime with multiplicity one
    fn is_prime(&self, a: &Self::Set) -> bool {
        debug_assert!(self.is_element(a));
        self.is_prime_unchecked(a)
    }

    /// Same as `is_prime` but don't validate inputs
    fn is_prime_unchecked(&self, a: &Self::Set) -> bool {
        let mut factor_powers = self.to_powers_unchecked(a).into_iter();
        match factor_powers.next() {
            Some((_, k)) => match factor_powers.next() {
                Some(_) => false,
                None => k == &Natural::ONE,
            },
            None => false,
        }
    }

    /// Return true if the factorization is square-free i.e. all multiplicities are zero or one
    fn is_squarefree(&self, a: &Self::Set) -> bool {
        debug_assert!(self.is_element(a));
        self.to_powers(a)
            .into_iter()
            .all(|(_, k)| k == &Natural::ONE)
    }

    /// The multiplicative identity object
    fn object_one(&self) -> Self::Object {
        self.object_product(vec![])
    }

    /// The product of two objects
    fn object_mul(&self, a: &Self::Object, b: &Self::Object) -> Self::Object {
        self.object_product(vec![a, b])
    }

    /// The product of objects
    fn object_product(&self, objects: Vec<&Self::Object>) -> Self::Object;

    /// Does a divide b?
    fn object_divides(&self, a: &Self::Object, b: &Self::Object) -> bool;

    /// Are objects equivelent with respect to divisibility?
    fn object_equivalent(&self, a: &Self::Object, b: &Self::Object) -> bool {
        self.object_divides(a, b) && self.object_divides(b, a)
    }

    /// Try to determine if an object is prime
    /// May be inconclusive but this is bad for debugging
    fn try_object_is_prime(&self, object: &Self::PrimeObject) -> Option<bool>;

    /// Are prime objects equivelent with respect to divisibility?
    fn prime_object_equivalent(&self, a: &Self::PrimeObject, b: &Self::PrimeObject) -> bool {
        self.object_equivalent(
            &self.prime_into_object(a.clone()),
            &self.prime_into_object(b.clone()),
        )
    }

    /// Convert a prime object into an object
    fn prime_into_object(&self, prime: Self::PrimeObject) -> Self::Object;

    /// Return the product of two factorizations
    fn mul(&self, a: Self::Set, b: Self::Set) -> Self::Set;

    /// Return a natural power of a factorization
    fn pow(self, a: Self::Set, n: &Natural) -> Self::Set {
        debug_assert!(self.is_element(&a));
        if *n == Natural::ZERO {
            self.new_trivial()
        } else if *n == Natural::ONE {
            a
        } else {
            debug_assert!(*n >= Natural::TWO);
            let bits: Vec<_> = n.bits().collect();
            let mut pows = vec![a.clone()];
            while pows.len() < bits.len() {
                pows.push(self.mul(pows.last().unwrap().clone(), pows.last().unwrap().clone()));
            }
            let count = bits.len();
            debug_assert_eq!(count, pows.len());
            let mut ans = self.new_trivial();
            for (i, pow) in pows.into_iter().enumerate() {
                if bits[i] {
                    ans = self.mul(ans, pow);
                }
            }
            ans
        }
    }

    /// Determine whether a divides b
    fn divides(&self, a: &Self::Set, b: &Self::Set) -> bool {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        for (pa, ka) in self.to_powers(a) {
            'EQUIV_PRIME_SEARCH: {
                for (pb, kb) in self.to_powers(b) {
                    if self.prime_object_equivalent(pa, pb) {
                        if ka > kb {
                            // `b` has the prime factor `pa` of `a` but its multiplicity is too small
                            return false;
                        }
                        break 'EQUIV_PRIME_SEARCH;
                    }
                }
                // `a` has a prime factor `pa` which `b` does not have
                return false;
            }
        }
        // every prime factor `pa` of `a` appears in `b` with at least the same multiplicity
        true
    }

    /// Are a and b equivelent with respect to division?
    fn equivalent(&self, a: &Self::Set, b: &Self::Set) -> bool {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        self.divides(a, b) && self.divides(b, a)
    }

    /// Expand a factorization into an object
    fn expanded(&self, a: &Self::Set) -> Self::Object;

    /// Expand a factorization into an object by seeing all multiplicities to one
    fn expanded_squarefree(&self, a: &Self::Set) -> Self::Object {
        debug_assert!(self.is_element(a));
        self.expanded(
            &self.new_powers(
                self.to_prime_support(a)
                    .into_iter()
                    .map(|p| (p.clone(), Natural::ONE))
                    .collect(),
            ),
        )
    }

    /// Return an iterator over all divisors of a factorization
    fn divisors<'a>(&'a self, a: &'a Self::Set) -> Box<dyn Iterator<Item = Self::Object> + 'a> {
        let factors = self.to_powers(a);
        if factors.is_empty() {
            Box::new(vec![self.object_one()].into_iter())
        } else {
            let mut factor_powers = vec![];
            for (p, k) in factors {
                let j = factor_powers.len();
                factor_powers.push(vec![]);
                let mut p_pow = self.object_one();
                let mut i = Natural::from(0u8);
                while &i <= k {
                    factor_powers[j].push(p_pow.clone());
                    p_pow = self.object_mul(&p_pow, &self.prime_into_object(p.clone()));
                    i += Natural::from(1u8);
                }
            }

            #[allow(clippy::redundant_closure_for_method_calls)]
            Box::new(
                itertools::Itertools::multi_cartesian_product(
                    factor_powers.into_iter().map(|p_pows| p_pows.into_iter()),
                )
                .map(move |prime_power_factors| {
                    self.object_product(prime_power_factors.iter().collect())
                        .clone()
                }),
            )
        }
    }

    /// The number of divisors of a factorization
    fn count_divisors(&self, a: &Self::Set) -> Natural {
        debug_assert!(self.is_element(a));
        let factors = self.to_powers(a);
        let mut count = Natural::from(1u8);
        for (_p, k) in factors {
            count *= k + Natural::ONE;
        }
        count
    }

    /// The greatest common divisor of two factorizations
    fn gcd(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        let factor_powers = self
            .to_powers(a)
            .into_iter()
            .filter_map(|(p, pk)| {
                for (q, qk) in self.to_powers(b) {
                    if self.prime_object_equivalent(p, q) {
                        return Some((p.clone(), std::cmp::min(pk, qk).clone()));
                    }
                }
                None
            })
            .collect();
        self.new_powers(factor_powers)
    }
}
