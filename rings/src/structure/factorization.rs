use algebraeon_nzq::*;
use algebraeon_sets::structure::*;
use std::fmt::Debug;

pub trait FactoredSignature: SetSignature {
    /// A type used to hold objects before they are factored
    type Object: Clone + Debug;
    /// A type used to hold the prime objects.
    type PrimeObject: Clone + Debug;

    fn from_factor_powers_impl(
        &self,
        factor_powers: Vec<(Self::PrimeObject, Natural)>,
    ) -> Self::Set;

    fn from_factor_powers(&self, factor_powers: Vec<(Self::PrimeObject, Natural)>) -> Self::Set {
        #[cfg(debug_assertions)]
        {
            let mut ps: Vec<&Self::PrimeObject> = vec![];
            for (p, _) in &factor_powers {
                assert!(self.object_is_prime(p));
                assert!(!ps.iter().any(|q| self.prime_object_equivalent(p, q)));
                ps.push(p);
            }
        }
        self.from_factor_powers_impl(
            factor_powers
                .into_iter()
                .filter(|(_, k)| k != &Natural::ZERO)
                .collect(),
        )
    }

    fn new_trivial(&self) -> Self::Set {
        self.from_factor_powers(vec![])
    }

    fn from_prime(&self, prime: Self::PrimeObject) -> Self::Set {
        debug_assert!(self.object_is_prime(&prime));
        self.from_factor_powers(vec![(prime, Natural::ONE)])
    }

    fn factor_powers<'a>(&self, a: &'a Self::Set) -> Vec<(&'a Self::PrimeObject, &'a Natural)>;

    fn into_factor_powers(&self, a: Self::Set) -> Vec<(Self::PrimeObject, Natural)>;

    fn factor_list<'a>(&self, a: &'a Self::Set) -> Vec<&'a Self::PrimeObject> {
        debug_assert!(self.is_element(a));
        let mut factors = vec![];
        for (p, k) in self.factor_powers(a) {
            let k: usize = k.try_into().unwrap();
            for _ in 0..k {
                factors.push(p);
            }
        }
        factors
    }

    fn into_factor_list(&self, a: Self::Set) -> Vec<Self::PrimeObject> {
        debug_assert!(self.is_element(&a));
        let mut factors = vec![];
        for (p, k) in self.into_factor_powers(a) {
            let k: usize = k.try_into().unwrap();
            for _ in 0..k {
                factors.push(p.clone())
            }
        }
        factors
    }

    fn squarefree_factor_list<'a>(&self, a: &'a Self::Set) -> Vec<&'a Self::PrimeObject> {
        debug_assert!(self.is_element(a));
        self.factor_powers(a).into_iter().map(|(p, _)| p).collect()
    }

    fn into_squarefree_factor_list(self, a: Self::Set) -> Vec<Self::PrimeObject> {
        debug_assert!(self.is_element(&a));
        self.into_factor_powers(a)
            .into_iter()
            .map(|(p, _)| p)
            .collect()
    }

    fn is_trivial(&self, a: &Self::Set) -> bool {
        debug_assert!(self.is_element(a));
        self.factor_powers(a).into_iter().next().is_none()
    }

    fn is_prime(&self, a: &Self::Set) -> bool {
        debug_assert!(self.is_element(a));
        let mut factor_powers = self.factor_powers(a).into_iter();
        match factor_powers.next() {
            Some((_, k)) => match factor_powers.next() {
                Some(_) => false,
                None => k == &Natural::ONE,
            },
            None => false,
        }
    }

    fn is_squarefree(&self, a: &Self::Set) -> bool {
        debug_assert!(self.is_element(a));
        self.factor_powers(a)
            .into_iter()
            .all(|(_, k)| k == &Natural::ONE)
    }

    /// return true iff a divides b
    fn object_divides(&self, a: &Self::Object, b: &Self::Object) -> bool;

    // not necessarily equal but equivalent wrt division
    fn object_equivalent(&self, a: &Self::Object, b: &Self::Object) -> bool {
        self.object_divides(a, b) && self.object_divides(b, a)
    }

    fn prime_object_equivalent(&self, a: &Self::PrimeObject, b: &Self::PrimeObject) -> bool {
        self.object_equivalent(
            &self.prime_to_object(a.clone()),
            &self.prime_to_object(b.clone()),
        )
    }

    /// return true if object actually is a prime object
    /// may return true if object doesn't represent a prime object, but this is not so good for debugging
    /// if returns false then object is definitely not a valid prime object
    fn object_is_prime(&self, object: &Self::PrimeObject) -> bool;

    fn prime_to_object(&self, prime: Self::PrimeObject) -> Self::Object;

    fn object_one(&self) -> Self::Object {
        self.object_product(vec![])
    }

    fn object_mul(&self, a: &Self::Object, b: &Self::Object) -> Self::Object {
        self.object_product(vec![a, b])
    }

    fn object_product(&self, objects: Vec<&Self::Object>) -> Self::Object;

    fn divides(&self, a: &Self::Set, b: &Self::Set) -> bool {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        for (pa, ka) in self.factor_powers(a) {
            'EQUIV_PRIME_SEARCH: {
                for (pb, kb) in self.factor_powers(b) {
                    if self.prime_object_equivalent(pa, pb) {
                        if !(ka <= kb) {
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

    fn equivalent(&self, a: &Self::Set, b: &Self::Set) -> bool {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        self.divides(a, b) && self.divides(b, a)
    }

    fn expanded(&self, a: &Self::Set) -> Self::Object;

    fn expanded_squarefree(&self, a: &Self::Set) -> Self::Object {
        debug_assert!(self.is_element(a));
        self.expanded(
            &self.from_factor_powers(
                self.squarefree_factor_list(a)
                    .into_iter()
                    .map(|p| (p.clone(), Natural::ONE))
                    .collect(),
            ),
        )
    }

    fn mul(&self, a: Self::Set, b: Self::Set) -> Self::Set;

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

    fn divisors<'a>(&'a self, a: &'a Self::Set) -> Box<dyn Iterator<Item = Self::Object> + 'a> {
        let factors = self.factor_powers(a);
        if factors.len() == 0 {
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
                    p_pow = self.object_mul(&p_pow, &self.prime_to_object(p.clone()));
                    i += Natural::from(1u8);
                }
            }

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

    fn count_divisors(&self, a: &Self::Set) -> Option<Natural> {
        debug_assert!(self.is_element(a));
        let factors = self.factor_powers(a);
        let mut count = Natural::from(1u8);
        for (_p, k) in factors {
            count *= k + Natural::ONE;
        }
        Some(count)
    }

    fn gcd(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        let factor_powers = self
            .factor_powers(a)
            .into_iter()
            .filter_map(|(p, pk)| {
                for (q, qk) in self.factor_powers(b) {
                    if self.prime_object_equivalent(&p, q) {
                        return Some((p.clone(), std::cmp::min(pk, qk).clone()));
                    }
                }
                None
            })
            .collect();
        self.from_factor_powers(factor_powers)
    }
}
