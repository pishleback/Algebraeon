use super::*;
use algebraeon_nzq::Natural;
use algebraeon_sets::structure::*;
use std::borrow::Borrow;
use std::{
    fmt::{Debug, Display},
    marker::PhantomData,
};

#[derive(Debug, Clone)]
pub struct FactoredRingElement<Element> {
    unit: Element,
    // all prime factors should satisfy is_fav_assoc()
    powers: Vec<(Element, Natural)>,
}

impl<Element> FactoredRingElement<Element> {
    pub fn new_unit(unit: Element) -> Self {
        Self {
            unit,
            powers: vec![],
        }
    }

    pub fn unit(&self) -> &Element {
        &self.unit
    }

    pub fn powers(&self) -> &Vec<(Element, Natural)> {
        &self.powers
    }

    pub fn from_unit_and_powers(unit: Element, powers: Vec<(Element, Natural)>) -> Self {
        Self { unit, powers }
    }

    pub fn into_unit_and_powers(self) -> (Element, Vec<(Element, Natural)>) {
        (self.unit, self.powers)
    }
}

pub trait RingFactorizationsSignature<
    Ring: UniqueFactorizationDomainSignature,
    RingB: BorrowedStructure<Ring>,
>:
    SetSignature<Set = FactoredRingElement<Ring::Set>>
    + FactoredSignature<Object = Ring::Set, PrimeObject = Ring::Set>
{
    fn ring(&self) -> &Ring;

    fn from_unit(&self, unit: Ring::Set) -> FactoredRingElement<Ring::Set>;

    fn from_unit_and_factor_powers_unchecked(
        &self,
        unit: Ring::Set,
        factor_powers: Vec<(Ring::Set, Natural)>,
    ) -> FactoredRingElement<Ring::Set>;

    fn from_unit_and_factor_powers(
        &self,
        unit: Ring::Set,
        factor_powers: Vec<(Ring::Set, Natural)>,
    ) -> FactoredRingElement<Ring::Set>;

    fn mul_mut(&self, a: &mut FactoredRingElement<Ring::Set>, b: FactoredRingElement<Ring::Set>);

    fn mul_by_unchecked(&self, a: &mut FactoredRingElement<Ring::Set>, p: Ring::Set, k: Natural);
}

pub trait UniqueFactorizationDomainSignature: FavoriteAssociateSignature {
    type FactorOrdering: OrdSignature<Set = Self::Set>;
    type Factorizations<SelfB: BorrowedStructure<Self>>: RingFactorizationsSignature<Self, SelfB>;

    fn factorizations<'a>(&'a self) -> Self::Factorizations<&'a Self>;

    fn into_factorizations(self) -> Self::Factorizations<Self>;

    fn factor_ordering(&self) -> impl Borrow<Self::FactorOrdering>;

    fn debug_try_is_irreducible(&self, a: &Self::Set) -> Option<bool>;
}

pub trait MetaUniqueFactorizationSignature: MetaType
where
    Self::Signature: UniqueFactorizationDomainSignature,
{
    fn debug_try_is_irreducible(&self) -> Option<bool> {
        Self::structure().debug_try_is_irreducible(self)
    }
}
impl<T: MetaType> MetaUniqueFactorizationSignature for T where
    T::Signature: UniqueFactorizationDomainSignature
{
}

pub trait FactorableSignature: UniqueFactorizationDomainSignature {
    fn factor(&self, element: &Self::Set) -> Option<FactoredRingElement<Self::Set>>;

    fn is_irreducible(&self, element: &Self::Set) -> bool {
        if let Some(factored) = self.factor(element) {
            self.factorizations().is_prime(&factored)
        } else {
            false
        }
    }

    fn gcd_by_factor(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
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

pub trait MetaFactorableSignature: MetaType
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
impl<T: MetaType> MetaFactorableSignature for T where T::Signature: FactorableSignature {}

impl<FS: FieldSignature> UniqueFactorizationDomainSignature for FS {
    type FactorOrdering = EmptySetStructure<Self::Set>;
    type Factorizations<SelfB: BorrowedStructure<Self>> =
        FieldElementFactorizationsStructure<Self, SelfB>;

    fn factorizations<'a>(&'a self) -> Self::Factorizations<&'a Self> {
        FieldElementFactorizationsStructure::new(self)
    }

    fn into_factorizations(self) -> Self::Factorizations<Self> {
        FieldElementFactorizationsStructure::new(self)
    }

    fn factor_ordering(&self) -> impl Borrow<Self::FactorOrdering> {
        EmptySetStructure::new()
    }

    fn debug_try_is_irreducible(&self, _: &Self::Set) -> Option<bool> {
        Some(false)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FieldElementFactorizationsStructure<FS: FieldSignature, FSB: BorrowedStructure<FS>> {
    _field: PhantomData<FS>,
    field: FSB,
}

impl<FS: FieldSignature, FSB: BorrowedStructure<FS>> FieldElementFactorizationsStructure<FS, FSB> {
    pub fn new(field: FSB) -> Self {
        Self {
            _field: PhantomData::default(),
            field,
        }
    }

    pub fn field(&self) -> &FS {
        self.field.borrow()
    }
}

impl<FS: FieldSignature, FSB: BorrowedStructure<FS>> Signature
    for FieldElementFactorizationsStructure<FS, FSB>
{
}

impl<FS: FieldSignature, FSB: BorrowedStructure<FS>> SetSignature
    for FieldElementFactorizationsStructure<FS, FSB>
{
    type Set = FactoredRingElement<FS::Set>;

    fn is_element(&self, x: &Self::Set) -> bool {
        self.field().is_unit(x.unit()) && x.powers().is_empty()
    }
}

impl<FS: FieldSignature, FSB: BorrowedStructure<FS>> FactoredSignature
    for FieldElementFactorizationsStructure<FS, FSB>
{
    type Object = FS::Set;
    type PrimeObject = FS::Set;

    fn new_powers_unchecked(&self, powers: Vec<(Self::PrimeObject, Natural)>) -> Self::Set {
        FactoredRingElement::from_unit_and_powers(self.field().one(), powers)
    }

    fn to_powers_unchecked<'a>(
        &self,
        a: &'a Self::Set,
    ) -> Vec<(&'a Self::PrimeObject, &'a Natural)> {
        a.powers().iter().map(|(p, k)| (p, k)).collect()
    }

    fn into_powers_unchecked(&self, powers: Self::Set) -> Vec<(Self::PrimeObject, Natural)> {
        powers.into_unit_and_powers().1
    }

    fn object_product(&self, objects: Vec<&Self::Object>) -> Self::Object {
        self.field().product(objects)
    }

    fn object_divides(&self, a: &Self::Object, b: &Self::Object) -> bool {
        self.field().is_zero(a) || self.field().is_unit(b)
    }

    fn try_object_is_prime(&self, object: &Self::PrimeObject) -> Option<bool> {
        self.field().debug_try_is_irreducible(object)
    }

    fn prime_into_object(&self, prime: Self::PrimeObject) -> Self::Object {
        prime
    }

    fn mul(&self, a: Self::Set, b: Self::Set) -> Self::Set {
        debug_assert!(self.is_element(&a));
        debug_assert!(self.is_element(&b));
        FactoredRingElement::new_unit(self.field().mul(a.unit(), b.unit()))
    }

    fn expanded(&self, a: &Self::Set) -> Self::Object {
        debug_assert!(self.is_element(a));
        a.unit().clone()
    }
}

impl<FS: FieldSignature, FSB: BorrowedStructure<FS>> RingFactorizationsSignature<FS, FSB>
    for FieldElementFactorizationsStructure<FS, FSB>
{
    fn ring(&self) -> &FS {
        self.field()
    }

    fn from_unit(&self, unit: FS::Set) -> FactoredRingElement<FS::Set> {
        FactoredRingElement::new_unit(unit)
    }

    fn from_unit_and_factor_powers_unchecked(
        &self,
        unit: FS::Set,
        powers: Vec<(FS::Set, Natural)>,
    ) -> FactoredRingElement<FS::Set> {
        FactoredRingElement::from_unit_and_powers(unit, powers)
    }

    fn from_unit_and_factor_powers(
        &self,
        unit: FS::Set,
        powers: Vec<(FS::Set, Natural)>,
    ) -> FactoredRingElement<FS::Set> {
        let f = FactoredRingElement::from_unit_and_powers(unit, powers);
        debug_assert!(self.is_element(&f));
        f
    }

    fn mul_mut(&self, a: &mut FactoredRingElement<FS::Set>, b: FactoredRingElement<FS::Set>) {
        *a = self.mul(a.clone(), b);
    }

    fn mul_by_unchecked(&self, _a: &mut FactoredRingElement<FS::Set>, _p: FS::Set, _k: Natural) {
        panic!();
    }
}

impl<FS: FieldSignature> FactorableSignature for FS {
    fn factor(&self, element: &FS::Set) -> Option<FactoredRingElement<FS::Set>> {
        if self.is_zero(element) {
            None
        } else {
            Some(FactoredRingElement::from_unit_and_powers(
                element.clone(),
                vec![],
            ))
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FactoredRingElementStructure<
    RS: UniqueFactorizationDomainSignature,
    RSB: BorrowedStructure<RS>,
> {
    _ring: PhantomData<RS>,
    ring: RSB,
}

impl<RS: UniqueFactorizationDomainSignature, RSB: BorrowedStructure<RS>>
    FactoredRingElementStructure<RS, RSB>
{
    pub fn new(ring: RSB) -> Self {
        Self {
            _ring: PhantomData::default(),
            ring,
        }
    }
}

impl<RS: UniqueFactorizationDomainSignature, RSB: BorrowedStructure<RS>> Signature
    for FactoredRingElementStructure<RS, RSB>
{
}

impl<RS: UniqueFactorizationDomainSignature, RSB: BorrowedStructure<RS>> SetSignature
    for FactoredRingElementStructure<RS, RSB>
{
    type Set = FactoredRingElement<RS::Set>;

    fn is_element(&self, x: &Self::Set) -> bool {
        if !self.ring().is_unit(&x.unit) {
            // unit must be a unit
            return false;
        }
        for (p, k) in &x.powers {
            if k == &Natural::ZERO {
                // prime powers must not be zero
                return false;
            }
            if !self.ring().is_fav_assoc(p) {
                // prime factor must be their favoriate associate
                return false;
            }
            if self.ring().debug_try_is_irreducible(p) == Some(false) {
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

impl<RS: UniqueFactorizationDomainSignature, RSB: BorrowedStructure<RS>> ToStringSignature
    for FactoredRingElementStructure<RS, RSB>
where
    RS: ToStringSignature,
{
    fn to_string(&self, elem: &Self::Set) -> String {
        format!("{:?}", elem)
    }
}

impl<RS: UniqueFactorizationDomainSignature, RSB: BorrowedStructure<RS>> EqSignature
    for FactoredRingElementStructure<RS, RSB>
{
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        let ring = self.ring();
        if !ring.equal(&a.unit, &b.unit) {
            false
        } else {
            //every a_factor is a b_factor
            for (a_factor, _a_power) in &a.powers {
                debug_assert!(ring.is_fav_assoc(a_factor));
                if !b
                    .powers
                    .iter()
                    .any(|(b_factor, _b_power)| ring.equal(a_factor, b_factor))
                {
                    return false;
                }
            }
            //every b_factor is an a_factor
            for (b_factor, _b_power) in &b.powers {
                debug_assert!(ring.is_fav_assoc(b_factor));
                if !a
                    .powers
                    .iter()
                    .any(|(a_factor, _a_power)| ring.equal(a_factor, b_factor))
                {
                    return false;
                }
            }
            //the powers of the factors are equal
            for (a_factor, a_power) in &a.powers {
                for (b_factor, b_power) in &b.powers {
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
        for (factor, k) in &self.powers {
            write!(f, " * (")?;
            write!(f, "{}", factor)?;
            write!(f, ")")?;
            if k != &Natural::ONE {
                write!(f, "^")?;
                write!(f, "{}", k)?;
            }
        }
        Ok(())
    }
}

impl<RS: UniqueFactorizationDomainSignature, RSB: BorrowedStructure<RS>> FactoredSignature
    for FactoredRingElementStructure<RS, RSB>
{
    type PrimeObject = RS::Set;
    type Object = RS::Set;

    fn object_divides(&self, a: &Self::Object, b: &Self::Object) -> bool {
        self.ring().divisible(a, b)
    }

    fn try_object_is_prime(&self, object: &Self::PrimeObject) -> Option<bool> {
        self.ring().debug_try_is_irreducible(object)
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
        a.powers.iter().map(|(p, k)| (p, k)).collect()
    }

    fn into_powers_unchecked(&self, a: Self::Set) -> Vec<(Self::PrimeObject, Natural)> {
        a.powers
    }

    fn expanded(&self, a: &Self::Set) -> Self::Object {
        debug_assert!(self.is_element(a));
        let mut ans = a.unit.clone();
        for (p, k) in &a.powers {
            self.ring().mul_mut(&mut ans, &self.ring().nat_pow(p, k));
        }
        ans
    }

    fn mul(&self, mut a: Self::Set, b: Self::Set) -> Self::Set {
        debug_assert!(self.is_element(&a));
        debug_assert!(self.is_element(&b));
        self.ring().mul_mut(&mut a.unit, &b.unit);
        for (p, k) in b.powers {
            self.mul_by_unchecked(&mut a, p, k);
        }
        a
    }
}

impl<RS: UniqueFactorizationDomainSignature, RSB: BorrowedStructure<RS>>
    RingFactorizationsSignature<RS, RSB> for FactoredRingElementStructure<RS, RSB>
{
    fn ring(&self) -> &RS {
        self.ring.borrow()
    }

    fn from_unit(&self, unit: RS::Set) -> FactoredRingElement<RS::Set> {
        debug_assert!(self.ring().is_unit(&unit));
        FactoredRingElement {
            unit,
            powers: vec![],
        }
    }

    fn from_unit_and_factor_powers_unchecked(
        &self,
        unit: RS::Set,
        powers: Vec<(RS::Set, Natural)>,
    ) -> FactoredRingElement<RS::Set> {
        FactoredRingElement {
            unit,
            powers: powers,
        }
    }

    fn from_unit_and_factor_powers(
        &self,
        mut unit: RS::Set,
        factor_powers: Vec<(RS::Set, Natural)>,
    ) -> FactoredRingElement<RS::Set> {
        debug_assert!(self.ring().is_unit(&unit));
        let mut powers = factor_powers
            .into_iter()
            .filter(|(_, k)| k > &Natural::ZERO)
            .collect::<Vec<_>>();
        for i in 0..powers.len() {
            let (p, k) = &mut powers[i];
            let (u, q) = self.ring().factor_fav_assoc(p);
            *p = q;
            unit = self.ring().mul(&unit, &self.ring().nat_pow(&u, k));
        }
        for (p, _) in &powers {
            debug_assert!(self.ring().is_fav_assoc(p));
        }
        FactoredRingElement { unit, powers }
    }

    fn mul_mut(&self, a: &mut FactoredRingElement<RS::Set>, b: FactoredRingElement<RS::Set>) {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(&b));
        *a = self.mul(a.clone(), b)
    }

    fn mul_by_unchecked(&self, a: &mut FactoredRingElement<RS::Set>, p: RS::Set, k: Natural) {
        for (q, t) in &mut a.powers {
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
        a.powers.push((p, k));
    }
}

#[derive(Debug)]
pub enum FindFactorResult<Element> {
    Irreducible,
    Composite(Element, Element),
}

pub fn factorize_by_find_factor<RS: UniqueFactorizationDomainSignature>(
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
        let f = Integer::structure()
            .factorizations()
            .from_unit_and_factor_powers_unchecked(
                Integer::from(-1),
                vec![
                    (Integer::from(2), Natural::from(2u8)),
                    (Integer::from(3), Natural::from(1u8)),
                ],
            );
        assert!(Integer::structure().factorizations().is_element(&f));

        let f = Integer::structure()
            .factorizations()
            .from_unit_and_factor_powers(Integer::from(1), vec![]);
        assert!(Integer::structure().factorizations().is_element(&f));

        let f = Integer::structure()
            .factorizations()
            .from_unit_and_factor_powers_unchecked(
                Integer::from(-1),
                vec![
                    (Integer::from(2), Natural::from(2u8)),
                    (Integer::from(3), Natural::from(1u8)),
                    (Integer::from(5), Natural::from(0u8)),
                ],
            );
        assert!(
            !Integer::structure().factorizations().is_element(&f),
            "can't have a power of zero"
        );

        let f = Integer::structure()
            .factorizations()
            .from_unit_and_factor_powers_unchecked(
                Integer::from(3),
                vec![(Integer::from(2), Natural::from(2u8))],
            );
        assert!(
            !Integer::structure().factorizations().is_element(&f),
            "unit should be a unit"
        );

        let f = Integer::structure()
            .factorizations()
            .from_unit_and_factor_powers_unchecked(
                Integer::from(1),
                vec![
                    (Integer::from(0), Natural::from(1u8)),
                    (Integer::from(3), Natural::from(1u8)),
                ],
            );
        assert!(
            !Integer::structure().factorizations().is_element(&f),
            "prime factors must not be zero"
        );

        let f = Integer::structure()
            .factorizations()
            .from_unit_and_factor_powers_unchecked(
                Integer::from(-1),
                vec![
                    (Integer::from(4), Natural::from(1u8)),
                    (Integer::from(3), Natural::from(1u8)),
                ],
            );
        assert!(
            !Integer::structure().factorizations().is_element(&f),
            "prime factors must be prime"
        );

        let f = Integer::structure()
            .factorizations()
            .from_unit_and_factor_powers_unchecked(
                Integer::from(-1),
                vec![
                    (Integer::from(-2), Natural::from(2u8)),
                    (Integer::from(3), Natural::from(1u8)),
                ],
            );
        assert!(
            !Integer::structure().factorizations().is_element(&f),
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
            println!("fs = {:?}", fs);
            assert_eq!(
                Integer::structure().factorizations().count_divisors(&fs),
                Natural::from(
                    Integer::structure()
                        .factorizations()
                        .divisors(&fs)
                        .collect::<Vec<Integer>>()
                        .len()
                )
            );
        }
    }
}
