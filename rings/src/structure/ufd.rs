use super::*;
use algebraeon_nzq::Natural;
use algebraeon_sets::structure::*;
use std::{borrow::Borrow, fmt::Display};

pub trait UniqueFactorizationSignature: FavoriteAssociateSignature {
    /// Try to determine if a is irreducible. May fail to produce an answer.
    fn try_is_irreducible(&self, a: &Self::Set) -> Option<bool>;
}

pub trait FactorableSignature: UniqueFactorizationSignature {
    //a UFD with an explicit algorithm to compute unique factorizations
    fn factor(&self, a: &Self::Set) -> Option<FactoredElement<Self>>;

    fn is_irreducible(&self, a: &Self::Set) -> bool {
        match self.factor(a) {
            None => false, //zero is not irreducible
            Some(factored) => factored.is_prime(),
        }
    }

    fn gcd_by_factor(&self, a: &Self::Set, b: &Self::Set) -> Self::Set
// where
    //     Self: FactoredAbstractStructure<Factored<Self>, Object = Self::Set, PrimeObject = Self::Set>,
        // Factored<Self>: FactoredAbstract<Structure = Self>,
    {
        match (self.factor(a), self.factor(b)) {
            (Some(factored_a), Some(factored_b)) => {
                FactoredElement::gcd(factored_a, factored_b).expanded()
            }
            (None, Some(_)) => b.clone(),
            (Some(_), None) => a.clone(),
            (None, None) => self.zero(),
        }
    }
}

pub trait MetaFactorableSignature: MetaFavoriteAssociate
where
    Self::Signature: FactorableSignature,
{
    fn factor(&self) -> Option<FactoredElement<Self::Signature>> {
        Self::structure().factor(self)
    }

    fn is_irreducible(&self) -> bool {
        Self::structure().is_irreducible(self)
    }

    fn gcd_by_factor(a: &Self, b: &Self) -> Self {
        Self::structure().gcd_by_factor(a, b)
    }
}
impl<R: MetaRing> MetaFactorableSignature for R where Self::Signature: FactorableSignature<Set = R> {}

impl<FS: FieldSignature> UniqueFactorizationSignature for FS {
    fn try_is_irreducible(&self, a: &Self::Set) -> Option<bool> {
        Some(self.is_irreducible(a))
    }
}

impl<FS: FieldSignature> FactorableSignature for FS {
    fn factor(&self, a: &Self::Set) -> Option<FactoredElement<Self>> {
        if self.is_zero(a) {
            None
        } else {
            Some(FactoredElement::from_unit_and_factor_powers(
                self.clone(),
                a.clone(),
                vec![],
            ))
        }
    }
}

#[derive(Debug, Clone)]
pub struct FactoredElement<RS: UniqueFactorizationSignature> {
    ring: RS,
    unit: RS::Set,
    // all prime factors should satisfy is_fav_assoc()
    factors: Vec<(RS::Set, Natural)>,
}

impl<RS: UniqueFactorizationSignature> Display for FactoredElement<RS>
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

impl<RS: FactorableSignature> FactoredElement<RS> {
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

impl<RS: UniqueFactorizationSignature> FactoredSignature<FactoredElement<RS>> for RS {
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

impl<RS: UniqueFactorizationSignature> FactoredElement<RS> {
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

impl<RS: UniqueFactorizationSignature> Factored for FactoredElement<RS> {
    type Structure = RS;

    fn factored_structure<'a>(&'a self) -> impl 'a + Borrow<Self::Structure> {
        &self.ring
    }

    fn from_factor_powers_impl(
        ring: Self::Structure,
        factor_powers: Vec<(
            <Self::Structure as FactoredSignature<Self>>::PrimeObject,
            Natural,
        )>,
    ) -> Self {
        let one = ring.one();
        Self::from_unit_and_factor_powers(ring, one, factor_powers)
    }

    fn factor_powers(
        &self,
    ) -> Vec<(
        &<Self::Structure as FactoredSignature<Self>>::PrimeObject,
        &Natural,
    )> {
        self.factors.iter().map(|(p, k)| (p, k)).collect()
    }

    fn into_factor_powers(
        self,
    ) -> Vec<(
        <Self::Structure as FactoredSignature<Self>>::PrimeObject,
        Natural,
    )> {
        self.factors
    }

    fn expanded(&self) -> <Self::Structure as FactoredSignature<Self>>::Object {
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

#[derive(Debug)]
pub enum FindFactorResult<RS: RingSignature> {
    Irreducible,
    Composite(RS::Set, RS::Set),
}

pub fn factorize_by_find_factor<RS: UniqueFactorizationSignature>(
    ring: &RS,
    elem: RS::Set,
    partial_factor: &impl Fn(RS::Set) -> FindFactorResult<RS>,
) -> FactoredElement<RS> {
    debug_assert!(!ring.is_zero(&elem));
    if ring.is_unit(&elem) {
        FactoredElement::from_unit(ring.clone().into(), elem)
    } else {
        debug_assert!(!ring.is_unit(&elem));
        match partial_factor(elem.clone()) {
            FindFactorResult::Composite(g, h) => FactoredElement::mul(
                factorize_by_find_factor(ring, g, partial_factor),
                factorize_by_find_factor(ring, h, partial_factor),
            ),
            FindFactorResult::Irreducible => {
                //f is irreducible
                FactoredElement::from_prime(ring.clone().into(), elem)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_nzq::Integer;

    #[test]
    fn factorization_invariants() {
        let f = FactoredElement::from_unit_and_factor_powers_unchecked(
            Integer::structure(),
            Integer::from(-1),
            vec![
                (Integer::from(2), Natural::from(2u8)),
                (Integer::from(3), Natural::from(1u8)),
            ],
        );
        f.check_invariants().unwrap();

        let f = FactoredElement::from_unit_and_factor_powers(
            Integer::structure(),
            Integer::from(1),
            vec![],
        );
        f.check_invariants().unwrap();

        let f = FactoredElement::from_unit_and_factor_powers_unchecked(
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

        let f = FactoredElement::from_unit_and_factor_powers_unchecked(
            Integer::structure(),
            Integer::from(3),
            vec![(Integer::from(2), Natural::from(2u8))],
        );
        assert_eq!(f.check_invariants().is_ok(), false, "unit should be a unit");

        let f = FactoredElement::from_unit_and_factor_powers_unchecked(
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

        let f = FactoredElement::from_unit_and_factor_powers_unchecked(
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

        let f = FactoredElement::from_unit_and_factor_powers_unchecked(
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
