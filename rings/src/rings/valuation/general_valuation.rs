use std::ops::{Add, Mul, Sub};

use crate::{
    rings::{
        integer::ideal::IntegerIdealsStructure,
        quotient::QuotientStructure,
        quotient_morphism::EuclideanDomainQuotienting,
        valuation::{Valuation, padic_rat_valuation},
    },
    structure::{
        AdditiveGroupSignature, FieldSignature, IdealsSignature, MetaRingEq,
        PrincipalSubringInclusion, RingHomomorphism, RingSignature,
    },
};
use algebraeon_nzq::{
    Integer, IntegerCanonicalStructure, Natural, RationalCanonicalStructure,
    traits::{Abs, Fraction},
};
use algebraeon_sets::structure::{InjectiveFunction, MetaType, SetSignature};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AdditiveValueGroup<T>
where
    T: Add<T, Output = T> + Ord,
{
    Infinity,
    Finite(T),
}

#[allow(clippy::non_canonical_partial_ord_impl)]
impl<T> PartialOrd for AdditiveValueGroup<T>
where
    T: Add<T, Output = T> + Ord,
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some({
            match (self, other) {
                (Self::Infinity, Self::Infinity) => std::cmp::Ordering::Equal,
                (Self::Infinity, Self::Finite(_)) => std::cmp::Ordering::Greater,
                (Self::Finite(_), Self::Infinity) => std::cmp::Ordering::Less,
                (Self::Finite(finite_self), Self::Finite(finite_other)) => {
                    finite_self.cmp(finite_other)
                }
            }
        })
    }
}
impl<T> Ord for AdditiveValueGroup<T>
where
    T: Add<T, Output = T> + Ord,
{
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other).unwrap()
    }
}

impl<T> Add<AdditiveValueGroup<T>> for AdditiveValueGroup<T>
where
    T: Add<T, Output = T> + Ord,
{
    type Output = AdditiveValueGroup<T>;
    fn add(self, other: Self) -> Self::Output {
        #[allow(clippy::unnested_or_patterns)]
        match (self, other) {
            (Self::Infinity, Self::Infinity)
            | (Self::Infinity, Self::Finite(_))
            | (Self::Finite(_), Self::Infinity) => Self::Infinity,
            (Self::Finite(a), Self::Finite(b)) => Self::Finite(a + b),
        }
    }
}

/*
Infinite - Infinite is undefined since it corresponds in valuation land to zero/zero
Infinite - Finite is Infinite since it corresponds in valuation land to zero/finite
Finite - Infinite is undefined since it corresponds in valuation land to finite/zero
*/
impl<T> Sub<AdditiveValueGroup<T>> for AdditiveValueGroup<T>
where
    T: Ord + Sub<T, Output = T> + Add<T, Output = T>,
{
    type Output = AdditiveValueGroup<T>;
    fn sub(self, other: Self) -> Self::Output {
        match (self, other) {
            (_, Self::Infinity) => {
                panic!("Can't subtract valuation Infinity")
            }
            (Self::Infinity, Self::Finite(_)) => Self::Infinity,
            (Self::Finite(a), Self::Finite(b)) => Self::Finite(a - b),
        }
    }
}

impl<T> From<Valuation> for AdditiveValueGroup<T>
where
    T: From<Integer> + Ord + Add<T, Output = T>,
{
    fn from(v: Valuation) -> Self {
        match v {
            Valuation::Infinity => Self::Infinity,
            Valuation::Finite(integer) => Self::Finite(integer.into()),
        }
    }
}

pub trait GeneralPlace {
    type DomainFieldSignature: FieldSignature;
    type ValuationGamma: Ord + Default + Add<Self::ValuationGamma, Output = Self::ValuationGamma>;

    fn nu(
        &self,
        r: &<Self::DomainFieldSignature as SetSignature>::Set,
    ) -> AdditiveValueGroup<Self::ValuationGamma>;
}

pub trait PlaceSignature : GeneralPlace {
    
    type ValuationRing: RingSignature;
    type MaximalIdealSignature: IdealsSignature<Self::ValuationRing, Self::ValuationRing>;
    type ResidueField: FieldSignature;

    fn try_to_valuation_ring(
        &self,
        r: &<Self::DomainFieldSignature as SetSignature>::Set,
    ) -> Result<<Self::ValuationRing as SetSignature>::Set, AdditiveValueGroup<Self::ValuationGamma>>
    {
        let zero_val = AdditiveValueGroup::Finite(Self::ValuationGamma::default());
        let r_val = self.nu(r);
        if r_val >= zero_val {
            Ok(self
                .valuation_inclusion()
                .try_preimage(r)
                .expect("The valuation was nonnegative so it should be in the valuation ring"))
        } else {
            Err(r_val)
        }
    }

    fn residue_field(&self) -> Self::ResidueField;

    fn valuation_inclusion(
        &self,
    ) -> impl RingHomomorphism<Self::ValuationRing, Self::DomainFieldSignature>
    + InjectiveFunction<Self::ValuationRing, Self::DomainFieldSignature>;

    fn valuation_quotient(&self) -> impl RingHomomorphism<Self::ValuationRing, Self::ResidueField>;
}

pub struct PAdicValuation(Natural);

impl GeneralPlace for PAdicValuation {
    type DomainFieldSignature = RationalCanonicalStructure;

    type ValuationGamma = Integer;

    fn nu(
        &self,
        r: &<Self::DomainFieldSignature as SetSignature>::Set,
    ) -> AdditiveValueGroup<Self::ValuationGamma>
    where
        <Self::ValuationGamma as MetaType>::Signature: AdditiveGroupSignature,
    {
        padic_rat_valuation(&self.0, r.clone()).into()
    }
}

impl PlaceSignature for PAdicValuation {

    // TODO: localize, this is only a smaller subset
    // don't see that implementor of RingSignature right now
    // so put something temporarily
    type ValuationRing = IntegerCanonicalStructure;

    // TODO: correspondingly change this ideal to be in the correct ambient ring
    type MaximalIdealSignature = IntegerIdealsStructure<Self::ValuationRing>;

    type ResidueField = QuotientStructure<IntegerCanonicalStructure, true>;
    fn residue_field(&self) -> Self::ResidueField {
        QuotientStructure::<IntegerCanonicalStructure, true>::new_field(
            IntegerCanonicalStructure {},
            self.0.clone().into(),
        )
    }

    // TODO: extend this inclusion to also take in sfuff without p in denominator
    fn valuation_inclusion(
        &self,
    ) -> impl RingHomomorphism<Self::ValuationRing, Self::DomainFieldSignature>
    + InjectiveFunction<Self::ValuationRing, Self::DomainFieldSignature> {
        PrincipalSubringInclusion::<RationalCanonicalStructure>::new(RationalCanonicalStructure {})
    }

    // TODO: modify to take into account localization
    fn valuation_quotient(&self) -> impl RingHomomorphism<Self::ValuationRing, Self::ResidueField> {
        EuclideanDomainQuotienting::new(IntegerCanonicalStructure {}, self.residue_field())
    }
}

#[derive(PartialEq, Eq, Debug, Clone)]
pub struct LogRatio<T>
where
    T: Ord + Mul<T, Output = T>,
{
    numerator: T,
    denominator: T,
}

impl<T> PartialOrd for LogRatio<T>
where
    T: Ord + Mul<T, Output = T> + Clone,
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        let cross_mult_0 = self.numerator.clone() * other.denominator.clone();
        let cross_mult_1 = other.numerator.clone() * self.denominator.clone();
        cross_mult_1.partial_cmp(&cross_mult_0)
    }
}

impl<T> Ord for LogRatio<T>
where
    T: Ord + Mul<T, Output = T> + Clone,
{
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other).unwrap()
    }
}

impl Default for LogRatio<Natural> {
    fn default() -> Self {
        Self {
            numerator: 1u8.into(),
            denominator: 1u8.into(),
        }
    }
}

impl<T: Ord + Mul<T, Output = T>> Add<Self> for LogRatio<T> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            numerator: self.numerator * rhs.numerator,
            denominator: self.denominator * rhs.denominator,
        }
    }
}
pub struct ArchimedeanPlace;

impl GeneralPlace for ArchimedeanPlace {
    type DomainFieldSignature = RationalCanonicalStructure;

    type ValuationGamma = LogRatio<Natural>;

    fn nu(
        &self,
        r: &<Self::DomainFieldSignature as SetSignature>::Set,
    ) -> AdditiveValueGroup<Self::ValuationGamma> {
        if r.is_zero() {
            return AdditiveValueGroup::Infinity;
        }
        let (numerator, denominator) = r.abs().numerator_and_denominator();
        let numerator = numerator.abs();
        AdditiveValueGroup::Finite(LogRatio {
            numerator,
            denominator,
        })
    }
}
