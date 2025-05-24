use algebraeon_nzq::*;
use algebraeon_sets::structure::*;
use std::borrow::Borrow;
use std::fmt::Debug;

pub trait FactoredSignature<F: Factored>: Signature {
    /// A type used to hold the prime objects.
    type PrimeObject: Clone + Debug;
    /// A type used to hold any object.
    type FactoredObject: Clone + Debug;

    /// return true iff a divides b
    fn object_divides(&self, a: &Self::FactoredObject, b: &Self::FactoredObject) -> bool;

    // not necessarily equal but equivalent wrt division
    fn object_equivalent(&self, a: &Self::FactoredObject, b: &Self::FactoredObject) -> bool {
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

    fn prime_to_object(&self, prime: Self::PrimeObject) -> Self::FactoredObject;

    fn object_one(&self) -> Self::FactoredObject {
        self.object_product(vec![])
    }

    fn object_mul(
        &self,
        a: &Self::FactoredObject,
        b: &Self::FactoredObject,
    ) -> Self::FactoredObject {
        self.object_product(vec![a, b])
    }

    fn object_product(&self, objects: Vec<&Self::FactoredObject>) -> Self::FactoredObject;
}

pub trait Factored: Debug + Clone + Sized {
    /// A structure for working with the above types.
    type Structure: FactoredSignature<Self>;

    fn factored_structure<'a>(&'a self) -> impl 'a + Borrow<Self::Structure>;

    fn from_factor_powers_impl(
        structure: Self::Structure,
        factor_powers: Vec<(
            <Self::Structure as FactoredSignature<Self>>::PrimeObject,
            Natural,
        )>,
    ) -> Self;

    fn from_factor_powers(
        structure: Self::Structure,
        factor_powers: Vec<(
            <Self::Structure as FactoredSignature<Self>>::PrimeObject,
            Natural,
        )>,
    ) -> Self {
        #[cfg(debug_assertions)]
        {
            let mut ps: Vec<&<Self::Structure as FactoredSignature<Self>>::PrimeObject> = vec![];
            for (p, _) in &factor_powers {
                assert!(structure.object_is_prime(p));
                assert!(!ps.iter().any(|q| structure.prime_object_equivalent(p, q)));
                ps.push(p);
            }
        }
        Self::from_factor_powers_impl(
            structure,
            factor_powers
                .into_iter()
                .filter(|(_, k)| k != &Natural::ZERO)
                .collect(),
        )
    }

    fn new_trivial(structure: Self::Structure) -> Self {
        Self::from_factor_powers(structure, vec![])
    }

    fn from_prime(
        structure: Self::Structure,
        prime: <Self::Structure as FactoredSignature<Self>>::PrimeObject,
    ) -> Self {
        debug_assert!(structure.object_is_prime(&prime));
        Self::from_factor_powers(structure, vec![(prime, Natural::ONE)])
    }

    fn factor_powers(
        &self,
    ) -> Vec<(
        &<Self::Structure as FactoredSignature<Self>>::PrimeObject,
        &Natural,
    )>;

    fn into_factor_powers(
        self,
    ) -> Vec<(
        <Self::Structure as FactoredSignature<Self>>::PrimeObject,
        Natural,
    )>;

    fn factor_list(&self) -> Vec<&<Self::Structure as FactoredSignature<Self>>::PrimeObject> {
        let mut factors = vec![];
        for (p, k) in self.factor_powers() {
            let k: usize = k.try_into().unwrap();
            for _ in 0..k {
                factors.push(p);
            }
        }
        factors
    }

    fn into_factor_list(self) -> Vec<<Self::Structure as FactoredSignature<Self>>::PrimeObject> {
        let mut factors = vec![];
        for (p, k) in self.into_factor_powers() {
            let k: usize = k.try_into().unwrap();
            for _ in 0..k {
                factors.push(p.clone())
            }
        }
        factors
    }

    fn squarefree_factor_list(
        &self,
    ) -> Vec<&<Self::Structure as FactoredSignature<Self>>::PrimeObject> {
        self.factor_powers().into_iter().map(|(p, _)| p).collect()
    }

    fn into_squarefree_factor_list(
        self,
    ) -> Vec<<Self::Structure as FactoredSignature<Self>>::PrimeObject> {
        self.into_factor_powers()
            .into_iter()
            .map(|(p, _)| p)
            .collect()
    }

    fn divides(a: &Self, b: &Self) -> bool {
        let structure =
            common_structure::<Self::Structure>(a.factored_structure(), b.factored_structure());
        for (pa, ka) in a.factor_powers() {
            'EQUIV_PRIME_SEARCH: {
                for (pb, kb) in b.factor_powers() {
                    if structure.prime_object_equivalent(pa, pb) {
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

    fn equivalent(a: &Self, b: &Self) -> bool {
        Self::divides(a, b) && Self::divides(b, a)
    }

    fn is_trivial(&self) -> bool {
        self.factor_powers().into_iter().next().is_none()
    }

    fn is_prime(&self) -> bool {
        let mut factor_powers = self.factor_powers().into_iter();
        match factor_powers.next() {
            Some((_, k)) => match factor_powers.next() {
                Some(_) => false,
                None => k == &Natural::ONE,
            },
            None => false,
        }
    }

    fn is_squarefree(&self) -> bool {
        self.factor_powers()
            .into_iter()
            .all(|(_, k)| k == &Natural::ONE)
    }

    fn expanded(&self) -> <Self::Structure as FactoredSignature<Self>>::FactoredObject;

    fn expanded_squarefree(&self) -> <Self::Structure as FactoredSignature<Self>>::FactoredObject {
        Self::from_factor_powers(
            self.factored_structure().borrow().clone(),
            self.squarefree_factor_list()
                .into_iter()
                .map(|p| (p.clone(), Natural::ONE))
                .collect(),
        )
        .expanded()
    }

    fn mul(a: Self, b: Self) -> Self;

    fn pow(self, n: &Natural) -> Self {
        let structure = self.factored_structure().borrow().clone();
        if *n == Natural::ZERO {
            Self::new_trivial(structure)
        } else if *n == Natural::ONE {
            self
        } else {
            debug_assert!(*n >= Natural::TWO);
            let bits: Vec<_> = n.bits().collect();
            let mut pows = vec![self.clone()];
            while pows.len() < bits.len() {
                pows.push(Self::mul(
                    pows.last().unwrap().clone(),
                    pows.last().unwrap().clone(),
                ));
            }
            let count = bits.len();
            debug_assert_eq!(count, pows.len());
            let mut ans = Self::new_trivial(structure);
            for (i, pow) in pows.into_iter().enumerate() {
                if bits[i] {
                    ans = Self::mul(ans, pow);
                }
            }
            ans
        }
    }

    fn divisors<'a>(
        &'a self,
    ) -> Box<dyn Iterator<Item = <Self::Structure as FactoredSignature<Self>>::FactoredObject> + 'a>
    {
        let structure = self.factored_structure();
        let factors = self.factor_powers();
        if factors.len() == 0 {
            Box::new(vec![structure.borrow().object_one()].into_iter())
        } else {
            let mut factor_powers = vec![];
            for (p, k) in factors {
                let j = factor_powers.len();
                factor_powers.push(vec![]);
                let mut p_pow = structure.borrow().object_one();
                let mut i = Natural::from(0u8);
                while &i <= k {
                    factor_powers[j].push(p_pow.clone());
                    p_pow = structure
                        .borrow()
                        .object_mul(&p_pow, &structure.borrow().prime_to_object(p.clone()));
                    i += Natural::from(1u8);
                }
            }

            Box::new(
                itertools::Itertools::multi_cartesian_product(
                    factor_powers.into_iter().map(|p_pows| p_pows.into_iter()),
                )
                .map(move |prime_power_factors| {
                    structure
                        .borrow()
                        .object_product(prime_power_factors.iter().collect())
                        .clone()
                }),
            )
        }
    }

    fn count_divisors(&self) -> Option<Natural> {
        let factors = self.factor_powers();
        let mut count = Natural::from(1u8);
        for (_p, k) in factors {
            count *= k + Natural::ONE;
        }
        Some(count)
    }

    fn gcd(a: Self, b: Self) -> Self {
        let structure = common_structure::<Self::Structure>(
            a.factored_structure().borrow(),
            b.factored_structure().borrow(),
        );
        let factor_powers = a
            .into_factor_powers()
            .into_iter()
            .filter_map(|(p, pk)| {
                for (q, qk) in b.factor_powers() {
                    if structure.prime_object_equivalent(&p, q) {
                        return Some((p, std::cmp::min(pk, qk.clone())));
                    }
                }
                None
            })
            .collect();
        Self::from_factor_powers(structure, factor_powers)
    }
}
