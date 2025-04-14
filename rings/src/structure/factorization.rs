use super::structure::*;
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;
use std::borrow::Borrow;
use std::fmt::Debug;
use std::fmt::Display;

pub trait FactoredAbstractStructure<F: FactoredAbstract>: Structure {
    /// A type used to hold the prime objects.
    type PrimeObject: Clone + Debug;
    /// A type used to hold any object.
    type Object: Clone + Debug;

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
}

pub trait FactoredAbstract: Debug + Clone + Sized {
    /// A structure for working with the above types.
    type Structure: FactoredAbstractStructure<Self>;

    fn factored_structure<'a>(&'a self) -> impl 'a + Borrow<Self::Structure>;

    fn from_factor_powers_impl(
        structure: Self::Structure,
        factor_powers: Vec<(
            <Self::Structure as FactoredAbstractStructure<Self>>::PrimeObject,
            Natural,
        )>,
    ) -> Self;

    fn from_factor_powers(
        structure: Self::Structure,
        factor_powers: Vec<(
            <Self::Structure as FactoredAbstractStructure<Self>>::PrimeObject,
            Natural,
        )>,
    ) -> Self {
        #[cfg(debug_assertions)]
        {
            let mut ps: Vec<&<Self::Structure as FactoredAbstractStructure<Self>>::PrimeObject> =
                vec![];
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
        prime: <Self::Structure as FactoredAbstractStructure<Self>>::PrimeObject,
    ) -> Self {
        debug_assert!(structure.object_is_prime(&prime));
        Self::from_factor_powers(structure, vec![(prime, Natural::ONE)])
    }

    fn factor_powers(
        &self,
    ) -> Vec<(
        &<Self::Structure as FactoredAbstractStructure<Self>>::PrimeObject,
        &Natural,
    )>;

    fn into_factor_powers(
        self,
    ) -> Vec<(
        <Self::Structure as FactoredAbstractStructure<Self>>::PrimeObject,
        Natural,
    )>;

    fn factor_list(
        &self,
    ) -> Vec<&<Self::Structure as FactoredAbstractStructure<Self>>::PrimeObject> {
        let mut factors = vec![];
        for (p, k) in self.factor_powers() {
            let k: usize = k.try_into().unwrap();
            for _ in 0..k {
                factors.push(p);
            }
        }
        factors
    }

    fn into_factor_list(
        self,
    ) -> Vec<<Self::Structure as FactoredAbstractStructure<Self>>::PrimeObject> {
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
    ) -> Vec<&<Self::Structure as FactoredAbstractStructure<Self>>::PrimeObject> {
        self.factor_powers().into_iter().map(|(p, _)| p).collect()
    }

    fn into_squarefree_factor_list(
        self,
    ) -> Vec<<Self::Structure as FactoredAbstractStructure<Self>>::PrimeObject> {
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

    fn expanded(&self) -> <Self::Structure as FactoredAbstractStructure<Self>>::Object;

    fn expanded_squarefree(&self) -> <Self::Structure as FactoredAbstractStructure<Self>>::Object {
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
    ) -> Box<dyn Iterator<Item = <Self::Structure as FactoredAbstractStructure<Self>>::Object> + 'a>
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

#[derive(Debug, Clone)]
pub struct Factored<RS: UniqueFactorizationStructure> {
    ring: RS,
    unit: RS::Set,
    // all prime factors should satisfy is_fav_assoc()
    factors: Vec<(RS::Set, Natural)>,
}

impl<RS: UniqueFactorizationStructure> Display for Factored<RS>
where
    RS::Set: Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.unit)?;
        for (factor, k) in &self.factors {
            write!(f, " * (")?;
            write!(f, "{}", factor)?;
            write!(f, ")")?;
            if k != &Natural::from(1u8) {
                write!(f, "^")?;
                write!(f, "{}", k)?;
            }
        }
        Ok(())
    }
}

impl<RS: FactorableStructure> Factored<RS> {
    pub fn check_invariants(&self) -> Result<(), &'static str> {
        if !self.ring.is_unit(&self.unit) {
            return Err("unit must be a unit");
        }
        for (p, k) in &self.factors {
            if k == &Natural::ZERO {
                return Err("prime powers must not be zero");
            }
            if !self.ring.is_fav_assoc(p) {
                return Err("prime factor must be their favoriate associate");
            }
            if !self.ring.is_irreducible(p) {
                return Err("prime factor must be irreducible");
            }

            let mut i = Natural::from(0u8);
            while &i < k {
                i += Natural::from(1u8);
            }
        }
        Ok(())
    }
}

impl<RS: UniqueFactorizationStructure> FactoredAbstractStructure<Factored<RS>> for RS {
    type PrimeObject = RS::Set;

    type Object = RS::Set;

    fn object_divides(&self, a: &Self::Object, b: &Self::Object) -> bool {
        self.divisible(a, b)
    }

    fn object_is_prime(&self, object: &Self::PrimeObject) -> bool {
        match self.try_is_irreducible(object) {
            Some(answer) => answer,
            None => {
                // unsure, so we are allowed to return true here
                true
            }
        }
    }

    fn prime_to_object(&self, prime: Self::PrimeObject) -> Self::Object {
        prime
    }

    fn object_product(&self, objects: Vec<&Self::Object>) -> Self::Object {
        self.product(objects)
    }
}

impl<RS: UniqueFactorizationStructure> Factored<RS> {
    pub(crate) fn equal(a: &Self, b: &Self) -> bool {
        let ring = common_structure::<RS>(&a.ring, &b.ring);
        if !ring.equal(&a.unit, &b.unit) {
            false
        } else {
            //every a_factor is a b_factor
            for (a_factor, _a_power) in &a.factors {
                debug_assert!(ring.is_fav_assoc(a_factor));
                if !b
                    .factors
                    .iter()
                    .any(|(b_factor, _b_power)| ring.equal(a_factor, b_factor))
                {
                    return false;
                }
            }
            //every b_factor is an a_factor
            for (b_factor, _b_power) in &b.factors {
                debug_assert!(ring.is_fav_assoc(b_factor));
                if !a
                    .factors
                    .iter()
                    .any(|(a_factor, _a_power)| ring.equal(a_factor, b_factor))
                {
                    return false;
                }
            }
            //the powers of the factors are equal
            for (a_factor, a_power) in &a.factors {
                for (b_factor, b_power) in &b.factors {
                    if ring.equal(a_factor, b_factor) {
                        if a_power != b_power {
                            return false;
                        }
                    }
                }
            }
            true
        }
    }

    pub fn from_unit(ring: RS, unit: RS::Set) -> Self {
        debug_assert!(ring.is_unit(&unit));
        Self {
            ring,
            unit,
            factors: vec![],
        }
    }

    pub fn from_unit_and_factor_powers_unchecked(
        ring: RS,
        unit: RS::Set,
        factor_powers: Vec<(RS::Set, Natural)>,
    ) -> Self {
        Self {
            ring,
            unit,
            factors: factor_powers,
        }
    }

    pub fn from_unit_and_factor_powers(
        ring: RS,
        mut unit: RS::Set,
        factor_powers: Vec<(RS::Set, Natural)>,
    ) -> Self {
        debug_assert!(ring.is_unit(&unit));
        let mut factor_powers = factor_powers
            .into_iter()
            .filter(|(_, k)| k > &Natural::ZERO)
            .collect::<Vec<_>>();
        for i in 0..factor_powers.len() {
            let (p, k) = &mut factor_powers[i];
            let (u, q) = ring.factor_fav_assoc(p);
            *p = q;
            unit = ring.mul(&unit, &ring.nat_pow(&u, k));
        }
        for (p, _) in &factor_powers {
            debug_assert!(ring.is_fav_assoc(p));
        }
        Self {
            ring,
            unit,
            factors: factor_powers,
        }
    }

    pub fn mul_mut(&mut self, other: Self) {
        *self = Self::mul(self.clone(), other)
    }

    pub fn ring(&self) -> &RS {
        &self.ring
    }

    pub fn into_unit_and_factor_powers(self) -> (RS::Set, Vec<(RS::Set, Natural)>) {
        (self.unit, self.factors)
    }

    fn mul_by_unchecked(&mut self, p: RS::Set, k: Natural) {
        for (q, t) in &mut self.factors {
            match (self.ring.div(&p, q), self.ring.div(q, &p)) {
                (Ok(u), Ok(v)) => {
                    if self.ring.is_unit(&u) && self.ring.is_unit(&v) {
                        //q = v*p so q^k = v^kp^k and this is what we are multiplying by
                        self.ring
                            .mul_mut(&mut self.unit, &self.ring.nat_pow(&v, &k));
                        *t += k;
                        return;
                    }
                }
                _ => {}
            }
        }
        self.factors.push((p, k));
    }
}

impl<RS: UniqueFactorizationStructure> FactoredAbstract for Factored<RS> {
    type Structure = RS;

    fn factored_structure<'a>(&'a self) -> impl 'a + Borrow<Self::Structure> {
        &self.ring
    }

    fn from_factor_powers_impl(
        ring: Self::Structure,
        factor_powers: Vec<(
            <Self::Structure as FactoredAbstractStructure<Self>>::PrimeObject,
            Natural,
        )>,
    ) -> Self {
        let one = ring.one();
        Self::from_unit_and_factor_powers(ring, one, factor_powers)
    }

    fn factor_powers(
        &self,
    ) -> Vec<(
        &<Self::Structure as FactoredAbstractStructure<Self>>::PrimeObject,
        &Natural,
    )> {
        self.factors.iter().map(|(p, k)| (p, k)).collect()
    }

    fn into_factor_powers(
        self,
    ) -> Vec<(
        <Self::Structure as FactoredAbstractStructure<Self>>::PrimeObject,
        Natural,
    )> {
        self.factors
    }

    fn expanded(&self) -> <Self::Structure as FactoredAbstractStructure<Self>>::Object {
        let mut ans = self.unit.clone();
        for (p, k) in &self.factors {
            self.ring.mul_mut(&mut ans, &self.ring.nat_pow(p, k));
        }
        ans
    }

    fn mul(mut a: Self, b: Self) -> Self {
        let ring = common_structure::<RS>(&a.ring, &b.ring);
        ring.mul_mut(&mut a.unit, &b.unit);
        for (p, k) in b.factors {
            a.mul_by_unchecked(p, k)
        }
        a
    }
}

/*

impl<RS: UniqueFactorizationStructure> Factored<RS> {
    pub fn ring(&self) -> &RS {
        &self.ring
    }

    //change to a new ring structure type
    //the new ring structure type must represent the same ring structure
    pub(crate) fn change_ring_unsafe<NewRS: UniqueFactorizationStructure<Set = RS::Set>>(
        self,
        ring: NewRS,
    ) -> Factored<NewRS> {
        Factored {
            ring,
            unit: self.unit,
            factors: self.factors,
        }
    }

    pub fn make_squarefree(self) -> Self {
        let one = self.ring.one();
        Self {
            ring: self.ring,
            unit: one,
            factors: self
                .factors
                .into_iter()
                .map(|(f, _k)| (f, Natural::ONE))
                .collect(),
        }
    }

    pub fn expand(&self) -> RS::Set {
        let mut ans = self.unit.clone();
        for (p, k) in &self.factors {
            self.ring.mul_mut(&mut ans, &self.ring.nat_pow(p, k));
        }
        ans
    }

    pub fn expand_squarefree(&self) -> RS::Set {
        self.clone().make_squarefree().expand()
    }

    pub fn new_unchecked(ring: RS, unit: RS::Set, factors: Vec<(RS::Set, Natural)>) -> Self {
        Self {
            ring,
            unit,
            factors,
        }
    }

    pub fn equal(a: &Self, b: &Self) -> bool {
        let ring = common_structure::<RS>(&a.ring, &b.ring);
        if !ring.equal(a.unit(), b.unit()) {
            false
        } else {
            //every a_factor is a b_factor
            for (a_factor, _a_power) in a.factors() {
                debug_assert!(ring.is_fav_assoc(a_factor));
                if !b
                    .factors()
                    .into_iter()
                    .any(|(b_factor, _b_power)| ring.equal(a_factor, b_factor))
                {
                    return false;
                }
            }
            //every b_factor is an a_factor
            for (b_factor, _b_power) in b.factors() {
                debug_assert!(ring.is_fav_assoc(b_factor));
                if !a
                    .factors()
                    .into_iter()
                    .any(|(a_factor, _a_power)| ring.equal(a_factor, b_factor))
                {
                    return false;
                }
            }
            //the powers of the factors are equal
            for (a_factor, a_power) in a.factors() {
                for (b_factor, b_power) in b.factors() {
                    if ring.equal(a_factor, b_factor) {
                        if a_power != b_power {
                            return false;
                        }
                    }
                }
            }
            true
        }
    }

    pub fn unit(&self) -> &RS::Set {
        &self.unit
    }

    pub fn factors(&self) -> &Vec<(RS::Set, Natural)> {
        &self.factors
    }

    pub fn into_factors(self) -> Vec<(RS::Set, Natural)> {
        self.factors
    }

    pub fn factors_list(&self) -> Vec<RS::Set> {
        let mut factors = vec![];
        for (factor, power) in &self.factors {
            let mut i = Natural::ZERO;
            while i < *power {
                factors.push(factor.clone());
                i += Natural::ONE;
            }
        }
        factors
    }

    pub fn unit_and_factors(self) -> (RS::Set, Vec<(RS::Set, Natural)>) {
        (self.unit, self.factors)
    }

    pub fn is_irreducible(&self) -> bool {
        if self.factors.len() == 1 {
            let (_p, k) = self.factors.iter().next().unwrap();
            if k == &Natural::ONE {
                return true;
            }
        }
        false
    }

    fn mul_by_unchecked(&mut self, p: RS::Set, k: Natural) {
        for (q, t) in &mut self.factors {
            match (self.ring.div(&p, q), self.ring.div(q, &p)) {
                (Ok(u), Ok(v)) => {
                    if self.ring.is_unit(&u) && self.ring.is_unit(&v) {
                        //q = v*p so q^k = v^kp^k and this is what we are multiplying by
                        self.ring
                            .mul_mut(&mut self.unit, &self.ring.nat_pow(&v, &k));
                        *t += k;
                        return;
                    }
                }
                _ => {}
            }
        }
        self.factors.push((p, k));
    }

    pub fn mul_mut(&mut self, other: Self) {
        let ring = common_structure::<RS>(&self.ring, &other.ring);
        ring.mul_mut(&mut self.unit, &other.unit);
        for (p, k) in other.factors {
            self.mul_by_unchecked(p, k);
        }
    }

    pub fn mul(mut a: Self, b: Self) -> Self {
        let ring = common_structure::<RS>(&a.ring, &b.ring);
        ring.mul_mut(&mut a.unit, &b.unit);
        for (p, k) in b.factors {
            a.mul_by_unchecked(p, k)
        }
        a
    }

    pub fn pow(mut self, k: &Natural) -> Self {
        for (_factor, power) in self.factors.iter_mut() {
            *power *= k;
        }
        self
    }

    pub fn factored_one(ring: RS) -> Self {
        Factored {
            unit: ring.one(),
            ring,
            factors: vec![],
        }
    }

    pub fn factored_irreducible_unchecked(ring: RS, elem: RS::Set) -> Self {
        let (unit, assoc) = ring.factor_fav_assoc(&elem);
        Factored {
            ring,
            unit,
            factors: vec![(assoc, Natural::from(1u8))],
        }
    }

    pub fn factored_irreducible_power_unchecked(ring: RS, elem: RS::Set, k: Natural) -> Self {
        debug_assert!(k >= Natural::ONE);
        let (unit, assoc) = ring.factor_fav_assoc(&elem);
        Factored {
            ring,
            unit,
            factors: vec![(assoc, k)],
        }
    }

    pub fn factored_unit_unchecked(ring: RS, unit: RS::Set) -> Self {
        Factored {
            ring,
            unit,
            factors: vec![],
        }
    }

    pub fn divisors<'a>(&'a self) -> Box<dyn Iterator<Item = RS::Set> + 'a> {
        if self.factors.len() == 0 {
            Box::new(vec![self.ring.one()].into_iter())
        } else {
            let mut factor_powers = vec![];
            for (p, k) in &self.factors {
                let j = factor_powers.len();
                factor_powers.push(vec![]);
                let mut p_pow = self.ring.one();
                let mut i = Natural::from(0u8);
                while &i <= k {
                    factor_powers[j].push(p_pow.clone());
                    p_pow = self.ring.mul(&p_pow, &p);
                    i += Natural::from(1u8);
                }
            }

            Box::new(
                itertools::Itertools::multi_cartesian_product(
                    factor_powers.into_iter().map(|p_pows| p_pows.into_iter()),
                )
                .map(|prime_power_factors| {
                    self.ring
                        .product(prime_power_factors.iter().collect())
                        .clone()
                }),
            )
        }
    }

    pub fn count_divisors(&self) -> Option<Natural> {
        let mut count = Natural::from(1u8);
        for (_p, k) in &self.factors {
            count *= k + Natural::from(1u8);
        }
        Some(count)
    }

    pub fn gcd(mut a: Self, b: Self) -> RS::Set {
        let ring = common_structure::<RS>(&a.ring, &b.ring);
        a.factors = a
            .factors
            .into_iter()
            .filter_map(|(p, pk)| {
                for (q, qk) in &b.factors {
                    if ring.equal(&p, q) {
                        return Some((p, std::cmp::min(pk, qk.clone())));
                    }
                }
                None
            })
            .collect();
        a.unit = ring.one();
        a.expand()
    }
}
*/

#[derive(Debug)]
pub enum FindFactorResult<RS: RingStructure> {
    Irreducible,
    Composite(RS::Set, RS::Set),
}

pub fn factorize_by_find_factor<RS: UniqueFactorizationStructure>(
    ring: &RS,
    elem: RS::Set,
    partial_factor: &impl Fn(RS::Set) -> FindFactorResult<RS>,
) -> Factored<RS> {
    debug_assert!(!ring.is_zero(&elem));
    if ring.is_unit(&elem) {
        Factored::from_unit(ring.clone().into(), elem)
    } else {
        debug_assert!(!ring.is_unit(&elem));
        match partial_factor(elem.clone()) {
            FindFactorResult::Composite(g, h) => Factored::mul(
                factorize_by_find_factor(ring, g, partial_factor),
                factorize_by_find_factor(ring, h, partial_factor),
            ),
            FindFactorResult::Irreducible => {
                //f is irreducible
                Factored::from_prime(ring.clone().into(), elem)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn factorization_invariants() {
        let f = Factored::from_unit_and_factor_powers_unchecked(
            Integer::structure(),
            Integer::from(-1),
            vec![
                (Integer::from(2), Natural::from(2u8)),
                (Integer::from(3), Natural::from(1u8)),
            ],
        );
        f.check_invariants().unwrap();

        let f =
            Factored::from_unit_and_factor_powers(Integer::structure(), Integer::from(1), vec![]);
        f.check_invariants().unwrap();

        let f = Factored::from_unit_and_factor_powers_unchecked(
            Integer::structure(),
            Integer::from(-1),
            vec![
                (Integer::from(2), Natural::from(2u8)),
                (Integer::from(3), Natural::from(1u8)),
                (Integer::from(5), Natural::from(0u8)),
            ],
        );
        assert_eq!(
            f.check_invariants().is_ok(),
            false,
            "can't have a power of zero"
        );

        let f = Factored::from_unit_and_factor_powers_unchecked(
            Integer::structure(),
            Integer::from(3),
            vec![(Integer::from(2), Natural::from(2u8))],
        );
        assert_eq!(f.check_invariants().is_ok(), false, "unit should be a unit");

        let f = Factored::from_unit_and_factor_powers_unchecked(
            Integer::structure(),
            Integer::from(1),
            vec![
                (Integer::from(0), Natural::from(1u8)),
                (Integer::from(3), Natural::from(1u8)),
            ],
        );
        assert_eq!(
            f.check_invariants().is_ok(),
            false,
            "prime factors must not be zero"
        );

        let f = Factored::from_unit_and_factor_powers_unchecked(
            Integer::structure(),
            Integer::from(-1),
            vec![
                (Integer::from(4), Natural::from(1u8)),
                (Integer::from(3), Natural::from(1u8)),
            ],
        );
        assert_eq!(
            f.check_invariants().is_ok(),
            false,
            "prime factors must be prime"
        );

        let f = Factored::from_unit_and_factor_powers_unchecked(
            Integer::structure(),
            Integer::from(-1),
            vec![
                (Integer::from(-2), Natural::from(2u8)),
                (Integer::from(3), Natural::from(1u8)),
            ],
        );
        assert_eq!(
            f.check_invariants().is_ok(),
            false,
            "prime factors must be favoriate associate"
        );
    }

    #[test]
    fn test_count_divisors() {
        for a in 1..25 {
            println!("a = {}", a);
            let b = Integer::from(a);
            println!("b = {}", b);
            let fs = Integer::structure().factor(&b).unwrap();
            println!("fs = {}", fs);
            assert_eq!(
                fs.count_divisors().unwrap(),
                Natural::from(fs.divisors().collect::<Vec<Integer>>().len())
            );
        }
    }
}
