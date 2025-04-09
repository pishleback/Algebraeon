use crate::polynomial::Polynomial;

use super::*;
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;

pub trait RingHomomorphism<Domain: RingStructure, Range: RingStructure>:
    Function<Domain, Range>
{
}

impl<
    A: RingStructure,
    B: RingStructure,
    C: RingStructure,
    AB: RingHomomorphism<A, B>,
    BC: RingHomomorphism<B, C>,
> RingHomomorphism<A, C> for CompositionMorphism<A, B, C, AB, BC>
{
}

/// The unique ring homomorphism Z -> R of the integers into any ring R
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PrincipalSubringInclusion<Ring: RingStructure> {
    integer_structure: IntegerCanonicalStructure,
    ring: Ring,
}

impl<Ring: RingStructure> PrincipalSubringInclusion<Ring> {
    pub fn new(ring: Ring) -> Self {
        Self {
            integer_structure: Integer::structure(),
            ring,
        }
    }
}

impl<Ring: RingStructure> Morphism<IntegerCanonicalStructure, Ring>
    for PrincipalSubringInclusion<Ring>
{
    fn domain(&self) -> &IntegerCanonicalStructure {
        &self.integer_structure
    }

    fn range(&self) -> &Ring {
        &self.ring
    }
}

impl<Ring: RingStructure> Function<IntegerCanonicalStructure, Ring>
    for PrincipalSubringInclusion<Ring>
{
    fn image(&self, x: &Integer) -> <Ring as SetStructure>::Set {
        self.range().from_int(x)
    }
}

impl<Ring: CharZeroStructure> InjectiveFunction<IntegerCanonicalStructure, Ring>
    for PrincipalSubringInclusion<Ring>
{
    fn try_preimage(
        &self,
        x: &<Ring as SetStructure>::Set,
    ) -> Option<<IntegerCanonicalStructure as SetStructure>::Set> {
        self.range().try_to_int(x)
    }
}

impl<Ring: RingStructure> RingHomomorphism<IntegerCanonicalStructure, Ring>
    for PrincipalSubringInclusion<Ring>
{
}

/// The inclusion of an integral domain into its field of fractions
pub trait FieldOfFractionsInclusion<Ring: RingStructure, Field: FieldStructure>:
    RingHomomorphism<Ring, Field> + InjectiveFunction<Ring, Field>
{
    fn numerator_and_denominator(&self, a: &Field::Set) -> (Ring::Set, Ring::Set);
    fn numerator(&self, a: &Field::Set) -> Ring::Set {
        self.numerator_and_denominator(a).0
    }
    fn denominator(&self, a: &Field::Set) -> Ring::Set {
        self.numerator_and_denominator(a).1
    }
}

/// A finite dimensional field extension F -> K
pub trait FiniteDimensionalFieldExtension<F: FieldStructure, K: FieldStructure>:
    RingHomomorphism<F, K> + InjectiveFunction<F, K>
{
    fn degree(&self) -> usize;
    fn norm(&self, a: &K::Set) -> F::Set;
    fn trace(&self, a: &K::Set) -> F::Set;
    /// The monic minimal polynomial of a
    fn min_poly(&self, a: &K::Set) -> Polynomial<F::Set>;
}

/// A seperable finite dimensional field extension F -> K
pub trait SeperableFiniteDimensionalFieldExtension<F: FieldStructure, K: FieldStructure>:
    FiniteDimensionalFieldExtension<F, K>
{
}
