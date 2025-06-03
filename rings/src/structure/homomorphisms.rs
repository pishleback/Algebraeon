use super::*;
use crate::polynomial::Polynomial;
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;
use std::borrow::{Borrow, Cow};
use std::collections::{HashMap, HashSet};
use std::fmt::Debug;
use std::marker::PhantomData;

#[derive(Debug, Clone)]
pub struct RingHomomorphismRangeModuleStructure<
    'h,
    Domain: RingSignature,
    Range: RingSignature,
    Hom: RingHomomorphism<Domain, Range>,
> {
    _domain: PhantomData<Domain>,
    _range: PhantomData<Range>,
    hom: &'h Hom,
}

impl<'h, Domain: RingSignature, Range: RingSignature, Hom: RingHomomorphism<Domain, Range>>
    RingHomomorphismRangeModuleStructure<'h, Domain, Range, Hom>
{
    fn new(hom: &'h Hom) -> Self {
        Self {
            _domain: PhantomData,
            _range: PhantomData,
            hom,
        }
    }
}

impl<'h, Domain: RingSignature, Range: RingSignature, Hom: RingHomomorphism<Domain, Range>>
    PartialEq for RingHomomorphismRangeModuleStructure<'h, Domain, Range, Hom>
{
    fn eq(&self, other: &Self) -> bool {
        std::ptr::eq(self.hom, other.hom)
    }
}

impl<'h, Domain: RingSignature, Range: RingSignature, Hom: RingHomomorphism<Domain, Range>> Eq
    for RingHomomorphismRangeModuleStructure<'h, Domain, Range, Hom>
{
}

impl<'h, Domain: RingSignature, Range: RingSignature, Hom: RingHomomorphism<Domain, Range>>
    Signature for RingHomomorphismRangeModuleStructure<'h, Domain, Range, Hom>
{
}

impl<'h, Domain: RingSignature, Range: RingSignature, Hom: RingHomomorphism<Domain, Range>>
    SetSignature for RingHomomorphismRangeModuleStructure<'h, Domain, Range, Hom>
{
    type Set = Range::Set;

    fn is_element(&self, x: &Self::Set) -> Result<(), String> {
        self.hom.range().is_element(x)
    }
}

impl<'h, Domain: RingSignature, Range: RingSignature, Hom: RingHomomorphism<Domain, Range>>
    EqSignature for RingHomomorphismRangeModuleStructure<'h, Domain, Range, Hom>
{
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        self.hom.range().equal(a, b)
    }
}

impl<'h, Domain: RingSignature, Range: RingSignature, Hom: RingHomomorphism<Domain, Range>>
    AdditiveMonoidSignature for RingHomomorphismRangeModuleStructure<'h, Domain, Range, Hom>
{
    fn zero(&self) -> Self::Set {
        self.hom.range().zero()
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.hom.range().add(a, b)
    }
}

impl<'h, Domain: RingSignature, Range: RingSignature, Hom: RingHomomorphism<Domain, Range>>
    AdditiveGroupSignature for RingHomomorphismRangeModuleStructure<'h, Domain, Range, Hom>
{
    fn neg(&self, a: &Self::Set) -> Self::Set {
        self.hom.range().neg(a)
    }
}

impl<'h, Domain: RingSignature, Range: RingSignature, Hom: RingHomomorphism<Domain, Range>>
    SemiModuleSignature<Domain> for RingHomomorphismRangeModuleStructure<'h, Domain, Range, Hom>
{
    fn ring(&self) -> &Domain {
        self.hom.domain()
    }

    fn scalar_mul(&self, x: &Domain::Set, a: &Self::Set) -> Self::Set {
        self.hom.range().mul(&self.hom.image(x), a)
    }
}

impl<
    'h,
    Domain: RingSignature,
    Range: RingSignature,
    Hom: FiniteRankFreeRingExtension<Domain, Range>,
> FreeModuleSignature<Domain> for RingHomomorphismRangeModuleStructure<'h, Domain, Range, Hom>
{
    type Basis = Hom::Basis;

    fn basis_set(&self) -> impl Borrow<Self::Basis> {
        self.hom.basis_set()
    }

    fn to_component<'a>(
        &self,
        b: &<Self::Basis as SetSignature>::Set,
        v: &'a Self::Set,
    ) -> Cow<'a, Domain::Set> {
        self.hom.to_component(b, v)
    }

    fn from_component(&self, b: &<Self::Basis as SetSignature>::Set, r: &Domain::Set) -> Self::Set {
        self.hom.from_component(b, r)
    }
}

pub trait RingHomomorphism<Domain: RingSignature, Range: RingSignature>:
    Function<Domain, Range>
{
    fn range_module_structure<'h>(
        &'h self,
    ) -> RingHomomorphismRangeModuleStructure<'h, Domain, Range, Self> {
        RingHomomorphismRangeModuleStructure::new(self)
    }
}

impl<
    A: RingSignature,
    B: RingSignature,
    C: RingSignature,
    AB: RingHomomorphism<A, B>,
    BC: RingHomomorphism<B, C>,
> RingHomomorphism<A, C> for CompositionMorphism<A, B, C, AB, BC>
{
}

/// The unique ring homomorphism Z -> R of the integers into any ring R
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PrincipalSubringInclusion<Ring: RingSignature> {
    integer_structure: IntegerCanonicalStructure,
    ring: Ring,
}

impl<Ring: RingSignature> PrincipalSubringInclusion<Ring> {
    pub fn new(ring: Ring) -> Self {
        Self {
            integer_structure: Integer::structure(),
            ring,
        }
    }
}

impl<Ring: RingSignature> Morphism<IntegerCanonicalStructure, Ring>
    for PrincipalSubringInclusion<Ring>
{
    fn domain(&self) -> &IntegerCanonicalStructure {
        &self.integer_structure
    }

    fn range(&self) -> &Ring {
        &self.ring
    }
}

impl<Ring: RingSignature> Function<IntegerCanonicalStructure, Ring>
    for PrincipalSubringInclusion<Ring>
{
    fn image(&self, x: &Integer) -> <Ring as SetSignature>::Set {
        self.range().from_int(x)
    }
}

impl<Ring: CharZeroRingSignature> InjectiveFunction<IntegerCanonicalStructure, Ring>
    for PrincipalSubringInclusion<Ring>
{
    fn try_preimage(&self, x: &<Ring as SetSignature>::Set) -> Option<Integer> {
        self.range().try_to_int(x)
    }
}

impl<Ring: RingSignature> RingHomomorphism<IntegerCanonicalStructure, Ring>
    for PrincipalSubringInclusion<Ring>
{
}

/// The unique field embedding Q -> K of the rationals into any field of characteristic zero
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PrincipalRationalSubfieldInclusion<
    Field: CharZeroFieldSignature,
    FieldB: BorrowedStructure<Field>,
> {
    rational_structure: RationalCanonicalStructure,
    _field: PhantomData<Field>,
    field: FieldB,
}

impl<Field: CharZeroFieldSignature, FieldB: BorrowedStructure<Field>>
    PrincipalRationalSubfieldInclusion<Field, FieldB>
{
    pub fn new(field: FieldB) -> Self {
        Self {
            rational_structure: Rational::structure(),
            _field: PhantomData,
            field,
        }
    }
}

impl<Field: CharZeroFieldSignature, FieldB: BorrowedStructure<Field>>
    Morphism<RationalCanonicalStructure, Field>
    for PrincipalRationalSubfieldInclusion<Field, FieldB>
{
    fn domain(&self) -> &RationalCanonicalStructure {
        &self.rational_structure
    }

    fn range(&self) -> &Field {
        self.field.borrow()
    }
}

impl<Field: CharZeroFieldSignature, FieldB: BorrowedStructure<Field>>
    Function<RationalCanonicalStructure, Field>
    for PrincipalRationalSubfieldInclusion<Field, FieldB>
{
    fn image(&self, x: &Rational) -> <Field as SetSignature>::Set {
        self.range().from_rat(x).unwrap()
    }
}

impl<Field: CharZeroFieldSignature, FieldB: BorrowedStructure<Field>>
    InjectiveFunction<RationalCanonicalStructure, Field>
    for PrincipalRationalSubfieldInclusion<Field, FieldB>
{
    fn try_preimage(&self, x: &<Field as SetSignature>::Set) -> Option<Rational> {
        self.range().try_to_rat(x)
    }
}

impl<Field: CharZeroFieldSignature, FieldB: BorrowedStructure<Field>>
    RingHomomorphism<RationalCanonicalStructure, Field>
    for PrincipalRationalSubfieldInclusion<Field, FieldB>
{
}

/// The inclusion of an integral domain into its field of fractions
pub trait FieldOfFractionsInclusion<Ring: RingSignature, Field: FieldSignature>:
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

/// An injective homomorphism of rings A -> B such that B is a finite rank free A module
pub trait FiniteRankFreeRingExtension<A: RingSignature, B: RingSignature>:
    RingHomomorphism<A, B> + InjectiveFunction<A, B>
{
    type Basis: FiniteSetSignature;
    fn basis_set(&self) -> impl Borrow<Self::Basis>;
    fn to_component<'a>(
        &self,
        b: &<Self::Basis as SetSignature>::Set,
        v: &'a B::Set,
    ) -> Cow<'a, A::Set>;
    fn from_component(&self, b: &<Self::Basis as SetSignature>::Set, r: &A::Set) -> B::Set;
}

/// A finite dimensional field extension F -> K
pub trait FiniteDimensionalFieldExtension<F: FieldSignature, K: FieldSignature>:
    FiniteRankFreeRingExtension<F, K>
{
    fn degree(&self) -> usize {
        self.basis_set().borrow().size()
    }
    fn norm(&self, a: &K::Set) -> F::Set;
    fn trace(&self, a: &K::Set) -> F::Set;
    /// The monic minimal polynomial of a
    fn min_poly(&self, a: &K::Set) -> Polynomial<F::Set>;
}

/// A separable finite dimensional field extension F -> K
pub trait SeparableFiniteDimensionalFieldExtension<F: FieldSignature, K: FieldSignature>:
    FiniteDimensionalFieldExtension<F, K>
{
}

/// Represent all ring homomorphisms from `domain` to `range`
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RingHomomorphisms<Domain: RingSignature, Range: RingSignature> {
    domain: Domain,
    range: Range,
}

impl<Domain: RingSignature, Range: RingSignature> RingHomomorphisms<Domain, Range> {
    pub fn new(domain: Domain, range: Range) -> Self {
        Self { domain, range }
    }
}

impl<Domain: RingSignature, Range: RingSignature> Signature for RingHomomorphisms<Domain, Range> {}

impl<Domain: FreeRingSignature, Range: RingSignature> SetSignature
    for RingHomomorphisms<Domain, Range>
{
    type Set = HashMap<Domain::Generator, Range::Set>;

    fn is_element(&self, x: &Self::Set) -> Result<(), String> {
        if self.domain.free_generators() != x.keys().cloned().collect::<HashSet<_>>() {
            return Err("missing key".to_string());
        }

        for v in x.values() {
            self.range.is_element(v)?;
        }

        Ok(())
    }
}
