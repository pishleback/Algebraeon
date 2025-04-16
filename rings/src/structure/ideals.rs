use super::*;
use algebraeon_nzq::Natural;
use algebraeon_sets::structure::common_structure;
use std::fmt::Debug;

pub trait IdealStructure: IntegralDomainStructure {
    type Ideal: Debug + Clone;
}
pub trait MetaIdealStructure: MetaIntegralDomain
where
    Self::Structure: IdealStructure,
{
}
impl<R: MetaRing> MetaIdealStructure for R where Self::Structure: IdealStructure<Set = R> {}

pub trait IdealArithmeticStructure: IdealStructure {
    // The zero ideal i.e. contains only 0
    fn zero_ideal(&self) -> Self::Ideal {
        self.principal_ideal(&self.zero())
    }
    // The unit ideal i.e. contains everything
    fn unit_ideal(&self) -> Self::Ideal {
        self.principal_ideal(&self.one())
    }
    /// The ideal generated by a list of elements
    fn generated_ideal(&self, elems: Vec<impl Into<Self::Set>>) -> Self::Ideal {
        self.ideal_sum(
            elems
                .into_iter()
                .map(|elem| self.principal_ideal(&elem.into()))
                .collect(),
        )
    }
    /// The principal ideal generated by a
    fn principal_ideal(&self, a: &Self::Set) -> Self::Ideal;
    /// Are the ideals equal?
    fn ideal_equal(&self, a: &Self::Ideal, b: &Self::Ideal) -> bool {
        self.ideal_contains(a, b) && self.ideal_contains(b, a)
    }
    /// Is the ideal the zero ideal
    fn ideal_is_zero(&self, a: &Self::Ideal) -> bool {
        self.ideal_equal(a, &self.zero_ideal())
    }
    /// Does a contain b i.e. does a divide b
    fn ideal_contains(&self, a: &Self::Ideal, b: &Self::Ideal) -> bool;
    /// Intersection of ideals
    fn ideal_intersect(&self, a: &Self::Ideal, b: &Self::Ideal) -> Self::Ideal;
    // Sum of two ideals
    fn ideal_add(&self, a: &Self::Ideal, b: &Self::Ideal) -> Self::Ideal;
    // Sum of many ideals
    fn ideal_sum(&self, ideals: Vec<impl Into<Self::Ideal>>) -> Self::Ideal {
        let mut total = self.zero_ideal();
        for i in ideals {
            total = self.ideal_add(&total, &i.into());
        }
        total
    }
    // Product of two ideals
    fn ideal_mul(&self, a: &Self::Ideal, b: &Self::Ideal) -> Self::Ideal;
    // Sum of many ideals
    fn ideal_product(&self, ideals: Vec<impl Into<Self::Ideal>>) -> Self::Ideal {
        let mut total = self.unit_ideal();
        for i in ideals {
            total = self.ideal_mul(&total, &i.into());
        }
        total
    }
}
pub trait MetaIdealArithmeticStructure: MetaIdealStructure
where
    Self::Structure: IdealArithmeticStructure,
{
    // The zero ideal i.e. contains only 0
    fn zero_ideal(&self) -> <Self::Structure as IdealStructure>::Ideal {
        Self::structure().zero_ideal()
    }
    // The unit ideal i.e. contains everything
    fn unit_ideal(&self) -> <Self::Structure as IdealStructure>::Ideal {
        Self::structure().unit_ideal()
    }
    /// The ideal generated by a list of elements
    fn generated_ideal(elems: Vec<impl Into<Self>>) -> <Self::Structure as IdealStructure>::Ideal {
        Self::structure().generated_ideal(elems)
    }
    /// The principal ideal generated by a
    fn principal_ideal(&self) -> <Self::Structure as IdealStructure>::Ideal {
        Self::structure().principal_ideal(self)
    }
    /// Are the ideals equal?
    fn ideal_equal(
        a: &<Self::Structure as IdealStructure>::Ideal,
        b: &<Self::Structure as IdealStructure>::Ideal,
    ) -> bool {
        Self::structure().ideal_equal(a, b)
    }
    /// Is the ideal the zero ideal
    fn ideal_is_zero(a: &<Self::Structure as IdealStructure>::Ideal) -> bool {
        Self::structure().ideal_is_zero(a)
    }
    /// Does a contain b i.e. does a divide b
    fn ideal_contains(
        a: &<Self::Structure as IdealStructure>::Ideal,
        b: &<Self::Structure as IdealStructure>::Ideal,
    ) -> bool {
        Self::structure().ideal_contains(a, b)
    }
    /// Intersection of ideals
    fn ideal_intersect(
        a: &<Self::Structure as IdealStructure>::Ideal,
        b: &<Self::Structure as IdealStructure>::Ideal,
    ) -> <Self::Structure as IdealStructure>::Ideal {
        Self::structure().ideal_intersect(a, b)
    }
    // Sum of two ideals
    fn ideal_add(
        a: &<Self::Structure as IdealStructure>::Ideal,
        b: &<Self::Structure as IdealStructure>::Ideal,
    ) -> <Self::Structure as IdealStructure>::Ideal {
        Self::structure().ideal_add(a, b)
    }
    // Sum of many ideals
    fn ideal_sum(
        ideals: Vec<impl Into<<Self::Structure as IdealStructure>::Ideal>>,
    ) -> <Self::Structure as IdealStructure>::Ideal {
        Self::structure().ideal_sum(ideals)
    }
    // Product of two ideals
    fn ideal_mul(
        a: &<Self::Structure as IdealStructure>::Ideal,
        b: &<Self::Structure as IdealStructure>::Ideal,
    ) -> <Self::Structure as IdealStructure>::Ideal {
        Self::structure().ideal_mul(a, b)
    }
    // Sum of many ideals
    fn ideal_product(
        ideals: Vec<impl Into<<Self::Structure as IdealStructure>::Ideal>>,
    ) -> <Self::Structure as IdealStructure>::Ideal {
        Self::structure().ideal_product(ideals)
    }
}
impl<R: MetaRing> MetaIdealArithmeticStructure for R where
    Self::Structure: IdealArithmeticStructure<Set = R>
{
}

pub trait PrincipalIdealDomainStructure: IdealStructure {
    fn ideal_generator(&self, ideal: &Self::Ideal) -> Self::Set;
}
pub trait MetaPrincipalIdealDomainStructure: MetaIdealStructure
where
    Self::Structure: PrincipalIdealDomainStructure,
{
    fn ideal_generator(ideal: &<Self::Structure as IdealStructure>::Ideal) -> Self {
        Self::structure().ideal_generator(ideal)
    }
}
impl<R: MetaRing> MetaPrincipalIdealDomainStructure for R where
    Self::Structure: PrincipalIdealDomainStructure<Set = R>
{
}

#[derive(Debug, Clone)]
pub struct DedekindDomainPrimeIdeal<RS: DedekindDomainStructure> {
    prime_ideal: RS::Ideal,
}

impl<RS: DedekindDomainStructure> DedekindDomainPrimeIdeal<RS> {
    pub fn from_ideal_unchecked(ideal: RS::Ideal) -> Self {
        Self { prime_ideal: ideal }
    }

    pub fn into_ideal(self) -> RS::Ideal {
        self.prime_ideal
    }

    pub fn ideal(&self) -> &RS::Ideal {
        &self.prime_ideal
    }
}

impl<RS: DedekindDomainStructure> FactoredAbstractStructure<DedekindDomainIdealFactorization<RS>>
    for RS
{
    type PrimeObject = DedekindDomainPrimeIdeal<RS>;

    type Object = RS::Ideal;

    fn object_divides(&self, a: &Self::Object, b: &Self::Object) -> bool {
        self.ideal_contains(a, b)
    }

    fn object_is_prime(&self, _object: &Self::PrimeObject) -> bool {
        true
    }

    fn prime_to_object(&self, prime: Self::PrimeObject) -> Self::Object {
        prime.into_ideal()
    }

    fn object_product(&self, objects: Vec<&Self::Object>) -> Self::Object {
        self.ideal_product(objects.into_iter().cloned().collect())
    }
}

#[derive(Debug, Clone)]
pub struct DedekindDomainIdealFactorization<RS: DedekindDomainStructure> {
    ring: RS,
    // The prime ideals should be distinct
    // All powers should be non-zero
    factors: Vec<(DedekindDomainPrimeIdeal<RS>, Natural)>,
}

impl<RS: DedekindDomainStructure> FactoredAbstract for DedekindDomainIdealFactorization<RS> {
    type Structure = RS;

    fn factored_structure<'a>(&'a self) -> impl 'a + std::borrow::Borrow<Self::Structure> {
        &self.ring
    }

    fn from_factor_powers_impl(
        structure: Self::Structure,
        factor_powers: Vec<(
            <Self::Structure as FactoredAbstractStructure<Self>>::PrimeObject,
            Natural,
        )>,
    ) -> Self {
        Self {
            ring: structure,
            factors: factor_powers,
        }
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
        self.ring.ideal_product(
            self.factor_list()
                .into_iter()
                .map(|p| p.ideal().clone())
                .collect(),
        )
    }

    fn mul(mut a: Self, b: Self) -> Self {
        let ring = common_structure::<RS>(a.factored_structure(), b.factored_structure());
        for (q, l) in b.into_factor_powers() {
            'SEARCH_A: {
                for (p, k) in &mut a.factors {
                    if ring.object_equivalent(p.ideal(), q.ideal()) {
                        *k += l;
                        break 'SEARCH_A;
                    }
                }
                a.factors.push((q, l));
            }
        }
        a
    }
}

/// A ring in which all ideals uniquely factor as a product of powers of prime ideals
pub trait DedekindDomainStructure: IdealArithmeticStructure {
    /// Return the largest power of prime_ideal which divides ideal
    fn largest_prime_ideal_factor_power(
        &self,
        prime_ideal: &DedekindDomainPrimeIdeal<Self>,
        ideal: &Self::Ideal,
    ) -> Natural {
        debug_assert!(!self.ideal_equal(prime_ideal.ideal(), &self.unit_ideal()));
        let mut k = Natural::ZERO;
        let mut prime_ideal_to_the_k_plus_one = prime_ideal.ideal().clone();
        while self.ideal_contains(&prime_ideal_to_the_k_plus_one, ideal) {
            k += Natural::ONE;
            prime_ideal_to_the_k_plus_one =
                self.ideal_mul(&prime_ideal_to_the_k_plus_one, prime_ideal.ideal())
        }
        k
    }
}

pub trait FactorableIdealsStructure: DedekindDomainStructure {
    fn factor_ideal(&self, ideal: &Self::Ideal) -> Option<DedekindDomainIdealFactorization<Self>>;
    fn is_prime_ideal(&self, ideal: &Self::Ideal) -> bool {
        if let Some(f) = self.factor_ideal(ideal) {
            f.is_prime()
        } else {
            false
        }
    }
}
pub trait MetaDedekindDomainStructure: MetaIdealArithmeticStructure
where
    Self::Structure: FactorableIdealsStructure,
{
    fn factor_ideal(
        ideal: &<Self::Structure as IdealStructure>::Ideal,
    ) -> Option<DedekindDomainIdealFactorization<Self::Structure>> {
        Self::structure().factor_ideal(ideal)
    }
    fn is_prime_ideal(ideal: &<Self::Structure as IdealStructure>::Ideal) -> bool {
        Self::structure().is_prime_ideal(ideal)
    }
}
impl<R: MetaRing> MetaDedekindDomainStructure for R where
    Self::Structure: FactorableIdealsStructure<Set = R>
{
}
