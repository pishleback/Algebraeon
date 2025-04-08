use crate::polynomial::Polynomial;

use super::*;
use algebraeon_nzq::Integer;
use algebraeon_sets::structure::*;

pub trait RingHomomorphismStructure<Domain: RingStructure, Range: RingStructure>:
    FunctionStructure<Domain, Range>
{
}

pub struct PrincipalSubringInclusion<Ring: RingStructure> {
    inclusion: Morphism<CannonicalStructure<Integer>, Ring>,
}

impl<Ring: RingStructure> PrincipalSubringInclusion<Ring> {
    pub fn new(ring: Ring) -> Self {
        Self {
            inclusion: Morphism::new(Integer::structure().clone(), ring),
        }
    }
}

impl<Ring: RingStructure> MorphismStructure<CannonicalStructure<Integer>, Ring>
    for PrincipalSubringInclusion<Ring>
{
    fn domain(&self) -> &CannonicalStructure<Integer> {
        self.inclusion.domain()
    }

    fn range(&self) -> &Ring {
        self.inclusion.range()
    }
}

impl<Ring: RingStructure> FunctionStructure<CannonicalStructure<Integer>, Ring>
    for PrincipalSubringInclusion<Ring>
{
    fn image(&self, x: &Integer) -> <Ring as SetStructure>::Set {
        self.range().from_int(x)
    }
}

impl<Ring: CharZeroStructure> InjectiveFunctionStructure<CannonicalStructure<Integer>, Ring>
    for PrincipalSubringInclusion<Ring>
{
    fn try_preimage(
        &self,
        x: &<Ring as SetStructure>::Set,
    ) -> Option<<CannonicalStructure<Integer> as SetStructure>::Set> {
        self.range().try_to_int(x)
    }
}

impl<Ring: RingStructure> RingHomomorphismStructure<CannonicalStructure<Integer>, Ring>
    for PrincipalSubringInclusion<Ring>
{
}

pub trait FieldOfFractionsInclusionStructure<Ring: RingStructure, Field: FieldStructure>:
    RingHomomorphismStructure<Ring, Field> + InjectiveFunctionStructure<Ring, Field>
{
    fn numerator_and_denominator(&self, a: &Field::Set) -> (Ring::Set, Ring::Set);
    fn numerator(&self, a: &Field::Set) -> Ring::Set {
        self.numerator_and_denominator(a).0
    }
    fn denominator(&self, a: &Field::Set) -> Ring::Set {
        self.numerator_and_denominator(a).1
    }
}

pub trait FiniteDimensionalFieldExtension<F: FieldStructure, K: FieldStructure>:
    RingHomomorphismStructure<F, K> + InjectiveFunctionStructure<F, K>
{
    fn degree(&self) -> usize;
    fn norm(&self, a: &K::Set) -> F::Set;
    fn trace(&self, a: &K::Set) -> F::Set;
    fn min_poly(&self, a: &K::Set) -> Polynomial<F::Set>;
}

/// Given a commuting square of injective ring homomorphisms
///
/// Q → K
/// ↑   ↑
/// Z → R
///
/// such that
///  - Q is the field of fractions of Z
///  - Q → K is a finite dimensional field extension
///
/// This trait expresses that R is the integral closure of Z in K
pub trait IntegralClosureSquare<
    Z: IntegralDomainStructure,
    R: IntegralDomainStructure,
    Q: FieldStructure,
    K: FieldStructure,
    ZR: RingHomomorphismStructure<Z, R> + InjectiveFunctionStructure<Z, R>,
    QK: FiniteDimensionalFieldExtension<Q, K>,
    ZQ: FieldOfFractionsInclusionStructure<Z, Q>,
    RK: RingHomomorphismStructure<R, K> + InjectiveFunctionStructure<R, K>,
>
{
    fn z_ring(&self) -> &Z;
    fn r_ring(&self) -> &R;
    fn q_field(&self) -> &Q;
    fn k_field(&self) -> &K;

    fn z_to_r(&self) -> &ZR;
    fn q_to_k(&self) -> &QK;
    fn z_to_q(&self) -> &ZQ;
    fn r_to_k(&self) -> &RK;

    /// Every element of K is a fraction of elements of R
    fn r_to_k_field_of_fractions(&self) -> impl FieldOfFractionsInclusionStructure<R, K>;

    /// Compose two squares
    ///
    /// Q → K → L
    /// ↑   ↑   ↑
    /// Z → R → S
    ///
    /// into
    ///
    /// Q → L
    /// ↑   ↑
    /// Z → S
    ///
    /// This is always possible because
    ///  - K is the field of fractions of R.
    ///  - The integral closure of R in L is equal to the integral closure of Z in L. Both are given by S.
    fn compose(square1: Self, square2: Self) -> Self;
}
