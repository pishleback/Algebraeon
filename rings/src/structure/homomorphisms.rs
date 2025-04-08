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
            inclusion: Morphism::new(Integer::structure().as_ref().clone(), ring),
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

impl<Ring: RingStructure> InjectiveFunctionStructure<CannonicalStructure<Integer>, Ring>
    for PrincipalSubringInclusion<Ring>
{
    fn try_preimage(
        &self,
        x: &<Ring as SetStructure>::Set,
    ) -> Option<<CannonicalStructure<Integer> as SetStructure>::Set> {
        todo!()
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

// pub trait FieldOfFractionsStructure<RS: IntegralDomainStructure>: FieldStructure {
//     fn base_ring_structure(&self) -> Rc<RS>;
//     fn from_base_ring(&self, elem: <RS as SetStructure>::Set) -> Self::Set;
//     fn numerator(&self, elem: &Self::Set) -> <RS as SetStructure>::Set;
//     fn denominator(&self, elem: &Self::Set) -> <RS as SetStructure>::Set;
//     fn as_base_ring(&self, elem: Self::Set) -> Option<<RS as SetStructure>::Set> {
//         let base_ring = self.base_ring_structure();
//         if base_ring.equal(&self.denominator(&elem), &base_ring.one()) {
//             Some(self.numerator(&elem))
//         } else {
//             None
//         }
//     }
// }
// pub trait MetaFieldOfFractions<RS: IntegralDomainStructure>: MetaRing
// where
//     Self::Structure: FieldOfFractionsStructure<RS>,
// {
//     fn base_ring_structure() -> Rc<RS> {
//         Self::structure().base_ring_structure()
//     }
//     fn from_base_ring(a: <RS as SetStructure>::Set) -> Self {
//         Self::structure().from_base_ring(a)
//     }
//     fn numerator(&self) -> <RS as SetStructure>::Set {
//         Self::structure().numerator(self)
//     }
//     fn denominator(&self) -> <RS as SetStructure>::Set {
//         Self::structure().denominator(self)
//     }
//     fn as_base_ring(self) -> Option<<RS as SetStructure>::Set> {
//         Self::structure().as_base_ring(self)
//     }
// }
// impl<R: MetaRing, RS: IntegralDomainStructure> MetaFieldOfFractions<RS> for R where
//     Self::Structure: FieldOfFractionsStructure<RS, Set = R>
// {
// }
