use super::*;
use algebraeon_nzq::Natural;
use algebraeon_sets::structure::*;
use std::{fmt::Display, marker::PhantomData};

pub trait UniqueFactorizationSignature: FavoriteAssociateSignature {
    /// Try to determine if a is irreducible. May fail to produce an answer.
    fn try_is_irreducible(&self, a: &Self::Set) -> Option<bool>;

    fn factorizations<'a>(&'a self) -> FactoredRingElementStructure<Self, &'a Self> {
        FactoredRingElementStructure::new(self)
    }

    fn into_factorizations(self) -> FactoredRingElementStructure<Self, Self> {
        FactoredRingElementStructure::new(self)
    }
}
pub trait MetaUniqueFactorizationSignature: MetaFavoriteAssociate
where
    Self::Signature: UniqueFactorizationSignature,
{
    fn try_is_irreducible(&self) -> Option<bool> {
        Self::structure().try_is_irreducible(self)
    }

    fn factorizations() -> FactoredRingElementStructure<Self::Signature, Self::Signature> {
        Self::structure().into_factorizations()
    }
}
impl<R: MetaRing> MetaUniqueFactorizationSignature for R where
    Self::Signature: UniqueFactorizationSignature<Set = R>
{
}

pub trait FactorableSignature: UniqueFactorizationSignature {
    //a UFD with an explicit algorithm to compute unique factorizations
    fn factor(&self, a: &Self::Set) -> Option<FactoredRingElement<Self::Set>>;

    fn is_irreducible(&self, a: &Self::Set) -> bool {
        match self.factor(a) {
            None => false, //zero is not irreducible
            Some(factored) => self.factorizations().is_prime_unchecked(&factored),
        }
    }

    fn gcd_by_factor(&self, a: &Self::Set, b: &Self::Set) -> Self::Set
// where
    //     Self: FactoredAbstractStructure<Factored<Self>, Object = Self::Set, PrimeObject = Self::Set>,
        // Factored<Self>: FactoredAbstract<Structure = Self>,
    {
        match (self.factor(a), self.factor(b)) {
            (Some(factored_a), Some(factored_b)) => self
                .factorizations()
                .expanded(&self.factorizations().gcd(&factored_a, &factored_b)),
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
    fn factor(&self) -> Option<FactoredRingElement<Self>> {
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
    fn factor(&self, a: &Self::Set) -> Option<FactoredRingElement<Self::Set>> {
        if self.is_zero(a) {
            None
        } else {
            Some(
                self.factorizations()
                    .from_unit_and_factor_powers(a.clone(), vec![]),
            )
        }
    }
}

#[derive(Debug, Clone)]
pub struct FactoredRingElement<Element> {
    unit: Element,
    // all prime factors should satisfy is_fav_assoc()
    factors: Vec<(Element, Natural)>,
}

impl<Element> FactoredRingElement<Element> {
    pub fn into_unit_and_factor_powers(self) -> (Element, Vec<(Element, Natural)>) {
        (self.unit, self.factors)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FactoredRingElementStructure<
    RS: UniqueFactorizationSignature,
    RSB: BorrowedStructure<RS>,
> {
    _ring: PhantomData<RS>,
    ring: RSB,
}

impl<RS: UniqueFactorizationSignature, RSB: BorrowedStructure<RS>>
    FactoredRingElementStructure<RS, RSB>
{
    pub fn new(ring: RSB) -> Self {
        Self {
            _ring: PhantomData::default(),
            ring,
        }
    }

    pub fn ring(&self) -> &RS {
        self.ring.borrow()
    }
}

impl<RS: UniqueFactorizationSignature, RSB: BorrowedStructure<RS>> Signature
    for FactoredRingElementStructure<RS, RSB>
{
}

impl<RS: UniqueFactorizationSignature, RSB: BorrowedStructure<RS>> SetSignature
    for FactoredRingElementStructure<RS, RSB>
{
    type Set = FactoredRingElement<RS::Set>;

    fn is_element(&self, x: &Self::Set) -> bool {
        if !self.ring().is_unit(&x.unit) {
            // unit must be a unit
            return false;
        }
        for (p, k) in &x.factors {
            if k == &Natural::ZERO {
                // prime powers must not be zero
                return false;
            }
            if !self.ring().is_fav_assoc(p) {
                // prime factor must be their favoriate associate
                return false;
            }
            if self.ring().try_is_irreducible(p) == Some(false) {
                // prime factor must be irreducible
                return false;
            }

            let mut i = Natural::from(0u8);
            while &i < k {
                i += Natural::from(1u8);
            }
        }
        return true;
    }
}

impl<RS: UniqueFactorizationSignature, RSB: BorrowedStructure<RS>> ToStringSignature
    for FactoredRingElementStructure<RS, RSB>
where
    Self::Set: std::fmt::Display,
{
    fn to_string(&self, elem: &Self::Set) -> String {
        format!("{}", elem)
    }
}

impl<RS: UniqueFactorizationSignature, RSB: BorrowedStructure<RS>> EqSignature
    for FactoredRingElementStructure<RS, RSB>
{
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        let ring = self.ring();
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
}

impl<Element> Display for FactoredRingElement<Element>
where
    Element: Display,
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

impl<RS: UniqueFactorizationSignature, RSB: BorrowedStructure<RS>> FactoredSignature
    for FactoredRingElementStructure<RS, RSB>
{
    type PrimeObject = RS::Set;
    type Object = RS::Set;

    fn object_divides(&self, a: &Self::Object, b: &Self::Object) -> bool {
        self.ring().divisible(a, b)
    }

    fn try_object_is_prime(&self, object: &Self::PrimeObject) -> Option<bool> {
        self.ring().try_is_irreducible(object)
    }

    fn prime_into_object(&self, prime: Self::PrimeObject) -> Self::Object {
        prime
    }

    fn object_product(&self, objects: Vec<&Self::Object>) -> Self::Object {
        self.ring().product(objects)
    }

    fn new_powers_unchecked(&self, factor_powers: Vec<(Self::PrimeObject, Natural)>) -> Self::Set {
        self.from_unit_and_factor_powers(self.ring().one(), factor_powers)
    }

    fn to_powers_unchecked<'a>(
        &self,
        a: &'a Self::Set,
    ) -> Vec<(&'a Self::PrimeObject, &'a Natural)> {
        a.factors.iter().map(|(p, k)| (p, k)).collect()
    }

    fn into_powers_unchecked(&self, a: Self::Set) -> Vec<(Self::PrimeObject, Natural)> {
        a.factors
    }

    fn expanded(&self, a: &Self::Set) -> Self::Object {
        debug_assert!(self.is_element(a));
        let mut ans = a.unit.clone();
        for (p, k) in &a.factors {
            self.ring().mul_mut(&mut ans, &self.ring().nat_pow(p, k));
        }
        ans
    }

    fn mul(&self, mut a: Self::Set, b: Self::Set) -> Self::Set {
        debug_assert!(self.is_element(&a));
        debug_assert!(self.is_element(&b));
        self.ring().mul_mut(&mut a.unit, &b.unit);
        for (p, k) in b.factors {
            self.mul_by_unchecked(&mut a, p, k);
        }
        a
    }
}

impl<RS: UniqueFactorizationSignature, RSB: BorrowedStructure<RS>>
    FactoredRingElementStructure<RS, RSB>
{
    pub fn from_unit(&self, unit: RS::Set) -> FactoredRingElement<RS::Set> {
        debug_assert!(self.ring().is_unit(&unit));
        FactoredRingElement {
            unit,
            factors: vec![],
        }
    }

    pub(crate) fn from_unit_and_factor_powers_unchecked(
        &self,
        unit: RS::Set,
        factor_powers: Vec<(RS::Set, Natural)>,
    ) -> FactoredRingElement<RS::Set> {
        FactoredRingElement {
            unit,
            factors: factor_powers,
        }
    }

    pub fn from_unit_and_factor_powers(
        &self,
        mut unit: RS::Set,
        factor_powers: Vec<(RS::Set, Natural)>,
    ) -> FactoredRingElement<RS::Set> {
        debug_assert!(self.ring().is_unit(&unit));
        let mut factor_powers = factor_powers
            .into_iter()
            .filter(|(_, k)| k > &Natural::ZERO)
            .collect::<Vec<_>>();
        for i in 0..factor_powers.len() {
            let (p, k) = &mut factor_powers[i];
            let (u, q) = self.ring().factor_fav_assoc(p);
            *p = q;
            unit = self.ring().mul(&unit, &self.ring().nat_pow(&u, k));
        }
        for (p, _) in &factor_powers {
            debug_assert!(self.ring().is_fav_assoc(p));
        }
        FactoredRingElement {
            unit,
            factors: factor_powers,
        }
    }

    pub fn mul_mut(&self, a: &mut FactoredRingElement<RS::Set>, b: FactoredRingElement<RS::Set>) {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(&b));
        *a = self.mul(a.clone(), b)
    }

    fn mul_by_unchecked(&self, a: &mut FactoredRingElement<RS::Set>, p: RS::Set, k: Natural) {
        for (q, t) in &mut a.factors {
            match (self.ring().div(&p, q), self.ring().div(q, &p)) {
                (Ok(u), Ok(v)) => {
                    if self.ring().is_unit(&u) && self.ring().is_unit(&v) {
                        //q = v*p so q^k = v^kp^k and this is what we are multiplying by
                        self.ring()
                            .mul_mut(&mut a.unit, &self.ring().nat_pow(&v, &k));
                        *t += k;
                        return;
                    }
                }
                _ => {}
            }
        }
        a.factors.push((p, k));
    }
}

#[derive(Debug)]
pub enum FindFactorResult<Element> {
    Irreducible,
    Composite(Element, Element),
}

pub fn factorize_by_find_factor<RS: UniqueFactorizationSignature>(
    ring: &RS,
    elem: RS::Set,
    partial_factor: &impl Fn(RS::Set) -> FindFactorResult<RS::Set>,
) -> FactoredRingElement<RS::Set> {
    debug_assert!(!ring.is_zero(&elem));
    if ring.is_unit(&elem) {
        ring.factorizations().from_unit(elem)
    } else {
        debug_assert!(!ring.is_unit(&elem));
        match partial_factor(elem.clone()) {
            FindFactorResult::Composite(g, h) => ring.factorizations().mul(
                factorize_by_find_factor(ring, g, partial_factor),
                factorize_by_find_factor(ring, h, partial_factor),
            ),
            FindFactorResult::Irreducible => {
                //f is irreducible
                ring.factorizations().new_prime(elem)
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
        let f = Integer::factorizations().from_unit_and_factor_powers_unchecked(
            Integer::from(-1),
            vec![
                (Integer::from(2), Natural::from(2u8)),
                (Integer::from(3), Natural::from(1u8)),
            ],
        );
        assert!(Integer::factorizations().is_element(&f));

        let f = Integer::factorizations().from_unit_and_factor_powers(Integer::from(1), vec![]);
        assert!(Integer::factorizations().is_element(&f));

        let f = Integer::factorizations().from_unit_and_factor_powers_unchecked(
            Integer::from(-1),
            vec![
                (Integer::from(2), Natural::from(2u8)),
                (Integer::from(3), Natural::from(1u8)),
                (Integer::from(5), Natural::from(0u8)),
            ],
        );
        assert!(
            !Integer::factorizations().is_element(&f),
            "can't have a power of zero"
        );

        let f = Integer::factorizations().from_unit_and_factor_powers_unchecked(
            Integer::from(3),
            vec![(Integer::from(2), Natural::from(2u8))],
        );
        assert!(
            !Integer::factorizations().is_element(&f),
            "unit should be a unit"
        );

        let f = Integer::factorizations().from_unit_and_factor_powers_unchecked(
            Integer::from(1),
            vec![
                (Integer::from(0), Natural::from(1u8)),
                (Integer::from(3), Natural::from(1u8)),
            ],
        );
        assert!(
            !Integer::factorizations().is_element(&f),
            "prime factors must not be zero"
        );

        let f = Integer::factorizations().from_unit_and_factor_powers_unchecked(
            Integer::from(-1),
            vec![
                (Integer::from(4), Natural::from(1u8)),
                (Integer::from(3), Natural::from(1u8)),
            ],
        );
        assert!(
            !Integer::factorizations().is_element(&f),
            "prime factors must be prime"
        );

        let f = Integer::factorizations().from_unit_and_factor_powers_unchecked(
            Integer::from(-1),
            vec![
                (Integer::from(-2), Natural::from(2u8)),
                (Integer::from(3), Natural::from(1u8)),
            ],
        );
        assert!(
            !Integer::factorizations().is_element(&f),
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
                Integer::factorizations().count_divisors(&fs),
                Natural::from(
                    Integer::factorizations()
                        .divisors(&fs)
                        .collect::<Vec<Integer>>()
                        .len()
                )
            );
        }
    }
}
