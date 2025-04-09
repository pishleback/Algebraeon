use algebraeon_nzq::Natural;

use super::*;

pub trait IdealStructure: IntegralDomainStructure {
    type Ideal;
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
    fn principal_ideal(a: &Self) -> <Self::Structure as IdealStructure>::Ideal {
        Self::structure().principal_ideal(a)
    }
    /// Are the ideals equal?
    fn ideal_equal(
        a: &<Self::Structure as IdealStructure>::Ideal,
        b: &<Self::Structure as IdealStructure>::Ideal,
    ) -> bool {
        Self::structure().ideal_equal(a, b)
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

pub struct DedekindDomainPrimeIdeal<RS: DedekindDomainStructure> {
    ideal: RS::Ideal,
}

pub struct DedekindDomainIdealFactorization<RS: DedekindDomainStructure> {
    factors: Vec<(DedekindDomainPrimeIdeal<RS>, Natural)>,
}

pub trait DedekindDomainStructure: IdealArithmeticStructure {
    fn factor_ideal(&self) -> DedekindDomainIdealFactorization<Self>;
    fn is_prime_ideal(&self, ideal: Self::Ideal) -> bool;
}
pub trait MetaDedekindDomainStructure: MetaIdealArithmeticStructure
where
    Self::Structure: DedekindDomainStructure,
{
}
impl<R: MetaRing> MetaDedekindDomainStructure for R where
    Self::Structure: DedekindDomainStructure<Set = R>
{
}
