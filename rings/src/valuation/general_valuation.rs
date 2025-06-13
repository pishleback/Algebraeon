use crate::{
    integer::ideal::IntegerIdealsStructure,
    localization::{LocalizationInclusionFoF, LocalizationPrime, LocalizationResidueField},
    structure::{
        AdditiveGroupSignature, EuclideanDomainQuotienting, FieldSignature,
        IdealsArithmeticSignature, MetaRingEq, QuotientStructure, RingHomomorphism, RingSignature,
        RingToIdealsSignature, RingToQuotientFieldSignature,
    },
    valuation::{Valuation, padic_rat_valuation},
};
use algebraeon_nzq::{
    Integer, IntegerCanonicalStructure, Natural, Rational, RationalCanonicalStructure,
    traits::{Abs, Fraction},
};
use algebraeon_sets::structure::{InjectiveFunction, MetaType, SetSignature};
use std::ops::{Add, Mul, Sub};

/// `AdditiveGroup<T>` is an extension of
/// a totally ordered abelian monoid `T` with an
/// `Infinity` which is greater than everything else and
/// absorbs the sum as `Infintity + x = x + Infinity = Infinity`
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

/// A type implementing `PreAdditiveValuation`
/// has the `nu` map to `AdditiveValueGroup<T>`
/// but is not necessarily obeying the ultra-metric property
/// For example, -`log` of the ordinary Archimedean absolute value
/// obeys
/// `-log(|ab|) = -log(|a|) + -log(|b|)`
/// but not
/// `-log(|a+b|)` being bound by the minimum of `-log(|a|)` and `-log(|b|)`
pub trait PreAdditiveValuation {
    type DomainFieldSignature: FieldSignature;
    type ValuationGamma: Ord + Default + Add<Self::ValuationGamma, Output = Self::ValuationGamma>;

    fn nu(
        &self,
        r: &<Self::DomainFieldSignature as SetSignature>::Set,
    ) -> AdditiveValueGroup<Self::ValuationGamma>;
}

/// A type implementing `AdditiveValuation`
/// does obeying the ultra-metric property
/// For example, the `padic_rat_valuation`
/// In this situation, we have a valuation ring
/// with unique maximal ideal by taking the inverse
/// image of `>=0` and `>0` of the `nu` map.
/// `0` is the `default` for the `T` of `AdditiveValueGroup<T>`
pub trait AdditiveValuation: PreAdditiveValuation {
    type ValuationRing: RingSignature;
    type ResidueField: FieldSignature;

    /// Evaluate the `nu` for this `r`
    /// If it is `>=0`, then it is in the valuation ring
    /// so we get an `Ok(Self::ValuationRing::Set)`
    /// If it is not, then `Err(AdditiveValueGroup<Self::ValuationGamma>)`
    /// showing that it was `<0`
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

    /// Give the signature for the residue field which is the quotient of the valuation ring
    /// by the maximal ideal determined by the preimage of `>0` under `nu`
    fn residue_field(&self) -> Self::ResidueField;

    /// Give the homomorphism for the valuation ring including into the field
    fn valuation_inclusion(
        &self,
    ) -> impl RingHomomorphism<Self::ValuationRing, Self::DomainFieldSignature>
    + InjectiveFunction<Self::ValuationRing, Self::DomainFieldSignature>;

    /// Give the homomorphism for the quotient map from the valuation ring to `self.residue_field()`
    fn valuation_quotient(&self) -> impl RingHomomorphism<Self::ValuationRing, Self::ResidueField>;
}

pub struct PAdicValuation(Natural);

impl PreAdditiveValuation for PAdicValuation {
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

impl AdditiveValuation for PAdicValuation {
    type ValuationRing = LocalizationPrime<
        IntegerCanonicalStructure,
        IntegerCanonicalStructure,
        IntegerIdealsStructure<IntegerCanonicalStructure>,
    >;

    type ResidueField =
        QuotientStructure<IntegerCanonicalStructure, IntegerCanonicalStructure, true>;

    fn residue_field(&self) -> Self::ResidueField {
        Integer::structure().into_quotient_field_unchecked(self.0.clone().into())
    }

    fn valuation_inclusion(
        &self,
    ) -> impl RingHomomorphism<Self::ValuationRing, Self::DomainFieldSignature>
    + InjectiveFunction<Self::ValuationRing, Self::DomainFieldSignature> {
        let z_to_q = Rational::structure().into_principal_subring_inclusion();
        let z = Integer::structure();
        let all_is = z.into_ideals();
        let prime_ideal = all_is.principal_ideal(&self.0.clone().into());
        let source =
            LocalizationPrime::new_unchecked(IntegerCanonicalStructure {}, all_is, prime_ideal);
        LocalizationInclusionFoF::new(source, z_to_q)
    }

    fn valuation_quotient(&self) -> impl RingHomomorphism<Self::ValuationRing, Self::ResidueField> {
        let z = IntegerCanonicalStructure {};
        let all_is = z.into_ideals();
        let prime_ideal = all_is.principal_ideal(&self.0.clone().into());
        let source =
            LocalizationPrime::new_unchecked(IntegerCanonicalStructure {}, all_is, prime_ideal);
        let target = self.residue_field();
        let map_restricted = EuclideanDomainQuotienting::new(IntegerCanonicalStructure {}, target);
        LocalizationResidueField::<
            IntegerCanonicalStructure,
            IntegerCanonicalStructure,
            IntegerIdealsStructure<IntegerCanonicalStructure>,
            Self::ResidueField,
            EuclideanDomainQuotienting<IntegerCanonicalStructure, true>,
        >::new(source, map_restricted)
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

impl PreAdditiveValuation for ArchimedeanPlace {
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
