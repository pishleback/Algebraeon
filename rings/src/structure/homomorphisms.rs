use std::collections::{HashMap, HashSet};

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

impl<Ring: CharZeroRingStructure> InjectiveFunction<IntegerCanonicalStructure, Ring>
    for PrincipalSubringInclusion<Ring>
{
    fn try_preimage(&self, x: &<Ring as SetStructure>::Set) -> Option<Integer> {
        self.range().try_to_int(x)
    }
}

impl<Ring: RingStructure> RingHomomorphism<IntegerCanonicalStructure, Ring>
    for PrincipalSubringInclusion<Ring>
{
}

/// The unique field embedding Q -> K of the rationals into any field of characteristic zero
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PrincipalRationalSubfieldInclusion<Field: CharZeroFieldStructure> {
    rational_structure: RationalCanonicalStructure,
    field: Field,
}

impl<Field: CharZeroFieldStructure> PrincipalRationalSubfieldInclusion<Field> {
    pub fn new(field: Field) -> Self {
        Self {
            rational_structure: Rational::structure(),
            field,
        }
    }
}

impl<Field: CharZeroFieldStructure> Morphism<RationalCanonicalStructure, Field>
    for PrincipalRationalSubfieldInclusion<Field>
{
    fn domain(&self) -> &RationalCanonicalStructure {
        &self.rational_structure
    }

    fn range(&self) -> &Field {
        &self.field
    }
}

impl<Field: CharZeroFieldStructure> Function<RationalCanonicalStructure, Field>
    for PrincipalRationalSubfieldInclusion<Field>
{
    fn image(&self, x: &Rational) -> <Field as SetStructure>::Set {
        self.range().from_rat(x).unwrap()
    }
}

impl<Field: CharZeroFieldStructure> InjectiveFunction<RationalCanonicalStructure, Field>
    for PrincipalRationalSubfieldInclusion<Field>
{
    fn try_preimage(&self, x: &<Field as SetStructure>::Set) -> Option<Rational> {
        self.range().try_to_rat(x)
    }
}

impl<Field: CharZeroFieldStructure> RingHomomorphism<RationalCanonicalStructure, Field>
    for PrincipalRationalSubfieldInclusion<Field>
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

/// Represent all ring homomorphisms from `domain` to `range`
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RingHomomorphisms<Domain: RingStructure, Range: RingStructure> {
    domain: Domain,
    range: Range,
}

impl<Domain: RingStructure, Range: RingStructure> RingHomomorphisms<Domain, Range> {
    pub fn new(domain: Domain, range: Range) -> Self {
        Self { domain, range }
    }
}

impl<Domain: RingStructure, Range: RingStructure> Structure for RingHomomorphisms<Domain, Range> {}

impl<Domain: FreeRingStructure, Range: RingStructure> SetStructure
    for RingHomomorphisms<Domain, Range>
{
    type Set = HashMap<Domain::Generator, Range::Set>;

    fn is_element(&self, x: &Self::Set) -> bool {
        self.domain.free_generators() == x.keys().cloned().collect::<HashSet<_>>()
            && x.iter().all(|(_, v)| self.range.is_element(v))
    }
}
