use crate::{matrix::Matrix, structure::*};
use algebraeon_sets::sets::FunctionsStructure;
use algebraeon_structures::*;
use std::{borrow::Cow, cmp::Ordering};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FinitelyFreeModuleStructure<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> {
    functions: FunctionsStructure<Set, SetB, Ring, RingB>,
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
    pub fn new(set: SetB, ring: RingB) -> Self {
        Self {
            functions: FunctionsStructure::new(set, ring),
        }
    }
}
pub trait RingToFinitelyFreeModuleSignature: RingSignature {
    fn free_module<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>>(
        &self,
        set: SetB,
    ) -> FinitelyFreeModuleStructure<Set, SetB, Self, &Self> {
        FinitelyFreeModuleStructure::new(set, self)
    }

    fn into_free_module<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>>(
        self,
        set: SetB,
    ) -> FinitelyFreeModuleStructure<Set, SetB, Self, Self> {
        FinitelyFreeModuleStructure::new(set, self)
    }
}
impl<Ring: RingSignature> RingToFinitelyFreeModuleSignature for Ring {}

pub trait SetToFinitelyFreeModuleSignature: EnumeratedOrdFiniteSetSignature {
    fn free_module<Ring: RingSignature, RingB: BorrowedStructure<Ring>>(
        &self,
        ring: RingB,
    ) -> FinitelyFreeModuleStructure<Self, &Self, Ring, RingB> {
        FinitelyFreeModuleStructure::new(self, ring)
    }

    fn into_free_module<Ring: RingSignature, RingB: BorrowedStructure<Ring>>(
        self,
        ring: RingB,
    ) -> FinitelyFreeModuleStructure<Self, Self, Ring, RingB> {
        FinitelyFreeModuleStructure::new(self, ring)
    }
}
impl<Set: EnumeratedOrdFiniteSetSignature> SetToFinitelyFreeModuleSignature for Set {}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
    pub fn set(&self) -> &Set {
        self.functions.domain()
    }

    pub fn ring(&self) -> &Ring {
        self.functions.range()
    }

    pub fn functions_restructure(&self) -> &FunctionsStructure<Set, SetB, Ring, RingB> {
        &self.functions
    }

    pub fn into_functions_restructure(self) -> FunctionsStructure<Set, SetB, Ring, RingB> {
        self.functions
    }

    pub fn to_col(&self, v: &<Self as SetSignature>::Elem) -> Matrix<Ring::Elem> {
        debug_assert!(self.validate_element(v).is_ok());
        Matrix::construct(self.rank(), 1, |r, _| v[r].clone())
    }

    pub fn to_row(&self, v: &<Self as SetSignature>::Elem) -> Matrix<Ring::Elem> {
        debug_assert!(self.validate_element(v).is_ok());
        Matrix::construct(1, self.rank(), |_, c| v[c].clone())
    }

    pub fn from_row(&self, m: &Matrix<Ring::Elem>) -> <Self as SetSignature>::Elem {
        debug_assert_eq!(m.rows(), 1);
        debug_assert_eq!(m.cols(), self.rank());
        (0..self.rank())
            .map(|i| m.at(0, i).unwrap().clone())
            .collect()
    }

    pub fn from_col(&self, m: &Matrix<Ring::Elem>) -> <Self as SetSignature>::Elem {
        debug_assert_eq!(m.cols(), 1);
        debug_assert_eq!(m.rows(), self.rank());
        (0..self.rank())
            .map(|i| m.at(i, 0).unwrap().clone())
            .collect()
    }

    pub fn basis_element(&self, i: usize) -> <Self as SetSignature>::Elem {
        debug_assert!(i < self.rank());
        (0..self.rank())
            .map(|j| {
                if i == j {
                    self.ring().one()
                } else {
                    self.ring().zero()
                }
            })
            .collect()
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> Signature for FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> SetSignature for FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
    type Elem = Vec<Ring::Elem>;

    fn validate_element(&self, v: &Self::Elem) -> Result<(), String> {
        self.functions_restructure().validate_element(v)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature + EqSignature,
    RingB: BorrowedStructure<Ring>,
> EqSignature for FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
    fn equal(&self, v: &Self::Elem, w: &Self::Elem) -> bool {
        self.functions_restructure().equal(v, w)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature + OrdSignature,
    RingB: BorrowedStructure<Ring>,
> PartialOrdSignature for FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
    fn partial_cmp(&self, v: &Self::Elem, w: &Self::Elem) -> Option<Ordering> {
        self.functions_restructure().partial_cmp(v, w)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature + OrdSignature,
    RingB: BorrowedStructure<Ring>,
> OrdSignature for FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
    fn cmp(&self, v: &Self::Elem, w: &Self::Elem) -> Ordering {
        self.functions_restructure().cmp(v, w)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature + EnumeratedOrdFiniteSetSignature,
    RingB: BorrowedStructure<Ring>,
> CountableSetSignature for FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
    fn into_generate_all_elements(self) -> impl Iterator<Item = Self::Elem> {
        self.into_functions_restructure()
            .into_generate_all_elements()
    }

    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> {
        self.functions_restructure().generate_all_elements()
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature + EnumeratedOrdFiniteSetSignature,
    RingB: BorrowedStructure<Ring>,
> FiniteSetSignature for FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
    fn list_all_elements(&self) -> Vec<Self::Elem> {
        self.functions_restructure().list_all_elements()
    }

    fn size(&self) -> Natural {
        self.functions_restructure().size()
    }

    fn generate_random_elements(&self, seed: u64) -> impl Iterator<Item = Self::Elem> {
        self.functions_restructure().generate_random_elements(seed)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature + EnumeratedOrdFiniteSetSignature,
    RingB: BorrowedStructure<Ring>,
> EnumeratedOrdFiniteSetSignature for FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
    fn list_all_elements_ordered(&self) -> Vec<Self::Elem> {
        self.functions_restructure().list_all_elements_ordered()
    }

    fn element_to_enumeration(&self, elem: &Self::Elem) -> Natural {
        self.functions_restructure().element_to_enumeration(elem)
    }

    fn enumeration_to_element(&self, num: &Natural) -> Option<Self::Elem> {
        self.functions_restructure().enumeration_to_element(num)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> RinglikeSpecializationSignature for FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> ZeroSignature for FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
    fn zero(&self) -> Self::Elem {
        (0..self.rank()).map(|_| self.ring().zero()).collect()
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> AdditionSignature for FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
    fn add(&self, v: &Self::Elem, w: &Self::Elem) -> Self::Elem {
        debug_assert!(self.validate_element(v).is_ok());
        debug_assert!(self.validate_element(w).is_ok());
        (0..self.rank())
            .map(|i| self.ring().add(&v[i], &w[i]))
            .collect()
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> CancellativeAdditionSignature for FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
    fn try_sub(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.sub(a, b))
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> TryNegateSignature for FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
    fn try_neg(&self, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.neg(a))
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> AdditiveMonoidSignature for FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> AdditiveGroupSignature for FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
    fn neg(&self, v: &Self::Elem) -> Self::Elem {
        debug_assert!(self.validate_element(v).is_ok());
        v.iter().map(|r| self.ring().neg(r)).collect()
    }

    fn sub(&self, v: &Self::Elem, w: &Self::Elem) -> Self::Elem {
        debug_assert!(self.validate_element(v).is_ok());
        debug_assert!(self.validate_element(w).is_ok());
        (0..self.rank())
            .map(|i| self.ring().sub(&v[i], &w[i]))
            .collect()
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> SemiModuleSignature<Ring> for FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
    fn ring(&self) -> &Ring {
        self.functions.range()
    }

    fn scalar_mul(&self, v: &Self::Elem, r: &Ring::Elem) -> Self::Elem {
        debug_assert!(self.validate_element(v).is_ok());
        v.iter().map(|s| self.ring().mul(r, s)).collect()
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> FreeModuleSignature<Set, Ring> for FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
    fn basis_set(&self) -> impl std::borrow::Borrow<Set> {
        self.set()
    }

    fn to_component<'a>(&self, b: &Set::Elem, v: &'a Self::Elem) -> Cow<'a, Ring::Elem> {
        let b: usize = self.set().element_to_enumeration(b).try_into().unwrap();
        debug_assert!(b < self.rank());
        Cow::Borrowed(&v[b])
    }

    fn from_component(&self, b: &Set::Elem, r: &<Ring>::Elem) -> Self::Elem {
        let b: usize = self.set().element_to_enumeration(b).try_into().unwrap();
        debug_assert!(b < self.rank());
        let mut element = self.zero();
        element[b] = r.clone();
        element
    }
}

#[cfg(test)]
mod tests {
    use crate::finite_fields::quaternary_field::QuaternaryField;

    use super::*;
    use algebraeon_sets::sets::EnumeratedFiniteSetStructure;

    #[test]
    fn enumeration() {
        algebraeon_structures::assert_enumerated_ord_finite_set!(
            FinitelyFreeModuleStructure::new(
                EnumeratedFiniteSetStructure::new(3),
                QuaternaryField::structure(),
            ),
            64
        );
    }

    #[test]
    fn test_finite_rank_modules() {
        let m = FinitelyFreeModuleStructure::new(
            EnumeratedFiniteSetStructure::new(3),
            Integer::structure(),
        );

        let a = m.basis_element(0);
        let b = m.basis_element(1);
        let c = m.basis_element(2);

        assert_eq!(
            m.add(&m.neg(&b), &m.add(&a, &b)),
            vec![Integer::from(1), Integer::from(0), Integer::from(0)]
        );

        assert_eq!(
            m.add(&m.add(&a, &b), &m.add(&b, &c)),
            vec![Integer::from(1), Integer::from(2), Integer::from(1)]
        );

        assert_eq!(
            m.scalar_mul(&a, &5.into()),
            vec![Integer::from(5), Integer::from(0), Integer::from(0)]
        );

        assert_eq!(m.basis_vecs(), vec![a, b, c]);
    }
}
