use super::structure::*;
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;
use std::fmt::Debug;
use std::fmt::Display;

pub trait FactoredAbstractStructure<F: FactoredAbstract>: Structure {
    /// A type used to hold the prime objects.
    type PrimeObject: Clone + Debug;
    /// A type used to hold any object.
    type Object: Clone + Debug;

    fn object_divides(&self, a: &Self::Object, b: &Self::Object) -> bool;

    fn object_equal(&self, a: &Self::Object, b: &Self::Object) -> bool {
        self.object_divides(a, b) && self.object_divides(b, a)
    }

    fn prime_object_equal(&self, a: &Self::PrimeObject, b: &Self::PrimeObject) -> bool {
        self.object_equal(
            &self.prime_to_object(a.clone()),
            &self.prime_to_object(b.clone()),
        )
    }

    fn object_is_prime(&self, object: &Self::PrimeObject) -> bool;

    fn prime_to_object(&self, prime: Self::PrimeObject) -> Self::Object;

    fn object_mul(&self, a: &Self::Object, b: &Self::Object) -> Self::Object;
}

pub trait FactoredAbstract: Sized {
    /// A structure for working with the above types.
    type Structure: FactoredAbstractStructure<Self>;

    fn from_factor_powers_unchecked(
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
                assert!(!ps.iter().any(|q| structure.prime_object_equal(p, q)));
                ps.push(p);
            }
        }
        Self::from_factor_powers_unchecked(
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

    fn expanded(&self) -> <Self::Structure as FactoredAbstractStructure<Self>>::Object {
        compile_error!("hi");
        todo!()
    }
}

#[derive(Debug, Clone)]
pub struct Factored<RS: UniqueFactorizationStructure> {
    ring: RS,
    unit: RS::Set,
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
        Factored::factored_unit_unchecked(ring.clone().into(), elem)
    } else {
        debug_assert!(!ring.is_unit(&elem));
        match partial_factor(elem.clone()) {
            FindFactorResult::Composite(g, h) => Factored::mul(
                factorize_by_find_factor(ring, g, partial_factor),
                factorize_by_find_factor(ring, h, partial_factor),
            ),
            FindFactorResult::Irreducible => {
                //f is irreducible
                Factored::factored_irreducible_unchecked(ring.clone().into(), elem)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn factorization_invariants() {
        let f = Factored::new_unchecked(
            Integer::structure(),
            Integer::from(-1),
            vec![
                (Integer::from(2), Natural::from(2u8)),
                (Integer::from(3), Natural::from(1u8)),
            ],
        );
        f.check_invariants().unwrap();

        let f = Factored::new_unchecked(Integer::structure(), Integer::from(1), vec![]);
        f.check_invariants().unwrap();

        let f = Factored::new_unchecked(
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

        let f = Factored::new_unchecked(
            Integer::structure(),
            Integer::from(3),
            vec![(Integer::from(2), Natural::from(2u8))],
        );
        assert_eq!(f.check_invariants().is_ok(), false, "unit should be a unit");

        let f = Factored::new_unchecked(
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

        let f = Factored::new_unchecked(
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

        let f = Factored::new_unchecked(
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
