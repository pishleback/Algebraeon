use crate::{
    linear::finitely_free_module::FinitelyFreeModuleStructure, matrix::Matrix, structure::*,
};
use algebraeon_sets::sets::ConstSizeFunctionsStructure;
use algebraeon_structures::*;
use std::{borrow::Cow, cmp::Ordering};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConstFinitelyFreeModuleStructure<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> {
    functions: ConstSizeFunctionsStructure<N, Set, SetB, Ring, RingB>,
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> ConstFinitelyFreeModuleStructure<N, Set, SetB, Ring, RingB>
{
    pub fn new(set: SetB, ring: RingB) -> Self {
        Self {
            functions: ConstSizeFunctionsStructure::new(set, ring),
        }
    }

    pub fn forget_const(&self) -> FinitelyFreeModuleStructure<Set, &Set, Ring, &Ring> {
        FinitelyFreeModuleStructure::new(self.functions.domain(), self.functions.range())
    }

    pub fn into_forget_const(self) -> FinitelyFreeModuleStructure<Set, SetB, Ring, RingB> {
        let (set, ring) = self.functions.into_forget_const().into_domain_and_range();
        FinitelyFreeModuleStructure::new(set, ring)
    }
}
pub trait RingToConstFinitelyFreeModuleSignature: RingSignature {
    fn free_module<
        const N: usize,
        Set: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
        SetB: BorrowedStructure<Set>,
    >(
        &self,
        set: SetB,
    ) -> ConstFinitelyFreeModuleStructure<N, Set, SetB, Self, &Self> {
        ConstFinitelyFreeModuleStructure::new(set, self)
    }

    fn into_free_module<
        const N: usize,
        Set: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
        SetB: BorrowedStructure<Set>,
    >(
        self,
        set: SetB,
    ) -> ConstFinitelyFreeModuleStructure<N, Set, SetB, Self, Self> {
        ConstFinitelyFreeModuleStructure::new(set, self)
    }
}
impl<Ring: RingSignature> RingToConstFinitelyFreeModuleSignature for Ring {}

pub trait SetToConstFinitelyFreeModuleSignature<const N: usize>:
    ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature
{
    fn free_module<Ring: RingSignature, RingB: BorrowedStructure<Ring>>(
        &self,
        ring: RingB,
    ) -> ConstFinitelyFreeModuleStructure<N, Self, &Self, Ring, RingB> {
        ConstFinitelyFreeModuleStructure::new(self, ring)
    }

    fn into_free_module<Ring: RingSignature, RingB: BorrowedStructure<Ring>>(
        self,
        ring: RingB,
    ) -> ConstFinitelyFreeModuleStructure<N, Self, Self, Ring, RingB> {
        ConstFinitelyFreeModuleStructure::new(self, ring)
    }
}
impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature>
    SetToConstFinitelyFreeModuleSignature<N> for Set
{
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> ConstFinitelyFreeModuleStructure<N, Set, SetB, Ring, RingB>
{
    pub fn set(&self) -> &Set {
        self.functions.domain()
    }

    pub fn ring(&self) -> &Ring {
        self.functions.range()
    }

    pub fn functions_restructure(&self) -> &ConstSizeFunctionsStructure<N, Set, SetB, Ring, RingB> {
        &self.functions
    }

    pub fn into_functions_restructure(
        self,
    ) -> ConstSizeFunctionsStructure<N, Set, SetB, Ring, RingB> {
        self.functions
    }

    pub fn to_col(&self, v: &<Self as SetSignature>::Elem) -> Matrix<Ring::Elem> {
        self.forget_const().to_col(&v.to_vec())
    }

    pub fn to_row(&self, v: &<Self as SetSignature>::Elem) -> Matrix<Ring::Elem> {
        self.forget_const().to_row(&v.to_vec())
    }

    pub fn from_row(&self, m: &Matrix<Ring::Elem>) -> <Self as SetSignature>::Elem {
        self.forget_const().from_row(m).try_into().unwrap()
    }

    pub fn from_col(&self, m: &Matrix<Ring::Elem>) -> <Self as SetSignature>::Elem {
        self.forget_const().from_col(m).try_into().unwrap()
    }

    pub fn basis_element(&self, i: usize) -> <Self as SetSignature>::Elem {
        self.forget_const().basis_element(i).try_into().unwrap()
    }
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> Signature for ConstFinitelyFreeModuleStructure<N, Set, SetB, Ring, RingB>
{
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> SetSignature for ConstFinitelyFreeModuleStructure<N, Set, SetB, Ring, RingB>
{
    type Elem = [Ring::Elem; N];

    fn validate_element(&self, v: &Self::Elem) -> Result<(), String> {
        self.functions_restructure().validate_element(v)
    }
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature + EqSignature,
    RingB: BorrowedStructure<Ring>,
> EqSignature for ConstFinitelyFreeModuleStructure<N, Set, SetB, Ring, RingB>
{
    fn equal(&self, v: &Self::Elem, w: &Self::Elem) -> bool {
        self.functions_restructure().equal(v, w)
    }
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature + OrdSignature,
    RingB: BorrowedStructure<Ring>,
> PartialOrdSignature for ConstFinitelyFreeModuleStructure<N, Set, SetB, Ring, RingB>
{
    fn partial_cmp(&self, v: &Self::Elem, w: &Self::Elem) -> Option<Ordering> {
        self.functions_restructure().partial_cmp(v, w)
    }
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature + OrdSignature,
    RingB: BorrowedStructure<Ring>,
> OrdSignature for ConstFinitelyFreeModuleStructure<N, Set, SetB, Ring, RingB>
{
    fn cmp(&self, v: &Self::Elem, w: &Self::Elem) -> Ordering {
        self.functions_restructure().cmp(v, w)
    }
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature + EnumeratedOrdFiniteSetSignature,
    RingB: BorrowedStructure<Ring>,
> CountableSetSignature for ConstFinitelyFreeModuleStructure<N, Set, SetB, Ring, RingB>
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
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature + EnumeratedOrdFiniteSetSignature,
    RingB: BorrowedStructure<Ring>,
> FiniteSetSignature for ConstFinitelyFreeModuleStructure<N, Set, SetB, Ring, RingB>
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
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature + EnumeratedOrdFiniteSetSignature,
    RingB: BorrowedStructure<Ring>,
> EnumeratedOrdFiniteSetSignature for ConstFinitelyFreeModuleStructure<N, Set, SetB, Ring, RingB>
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
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> RinglikeSpecializationSignature for ConstFinitelyFreeModuleStructure<N, Set, SetB, Ring, RingB>
{
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> ZeroSignature for ConstFinitelyFreeModuleStructure<N, Set, SetB, Ring, RingB>
{
    fn zero(&self) -> Self::Elem {
        std::array::from_fn(|_| self.ring().zero())
    }
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> AdditionSignature for ConstFinitelyFreeModuleStructure<N, Set, SetB, Ring, RingB>
{
    fn add(&self, v: &Self::Elem, w: &Self::Elem) -> Self::Elem {
        debug_assert!(self.validate_element(v).is_ok());
        debug_assert!(self.validate_element(w).is_ok());
        std::array::from_fn(|i| self.ring().add(&v[i], &w[i]))
    }
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> CancellativeAdditionSignature for ConstFinitelyFreeModuleStructure<N, Set, SetB, Ring, RingB>
{
    fn try_sub(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.sub(a, b))
    }
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> TryNegateSignature for ConstFinitelyFreeModuleStructure<N, Set, SetB, Ring, RingB>
{
    fn try_neg(&self, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.neg(a))
    }
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> AdditiveMonoidSignature for ConstFinitelyFreeModuleStructure<N, Set, SetB, Ring, RingB>
{
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> AdditiveGroupSignature for ConstFinitelyFreeModuleStructure<N, Set, SetB, Ring, RingB>
{
    fn neg(&self, v: &Self::Elem) -> Self::Elem {
        debug_assert!(self.validate_element(v).is_ok());
        std::array::from_fn(|i| self.ring().neg(&v[i]))
    }

    fn sub(&self, v: &Self::Elem, w: &Self::Elem) -> Self::Elem {
        debug_assert!(self.validate_element(v).is_ok());
        debug_assert!(self.validate_element(w).is_ok());
        std::array::from_fn(|i| self.ring().sub(&v[i], &w[i]))
    }
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> SemiModuleSignature<Ring> for ConstFinitelyFreeModuleStructure<N, Set, SetB, Ring, RingB>
{
    fn ring(&self) -> &Ring {
        self.functions.range()
    }

    fn scalar_mul(&self, v: &Self::Elem, r: &Ring::Elem) -> Self::Elem {
        debug_assert!(self.validate_element(v).is_ok());
        std::array::from_fn(|i| self.ring().mul(r, &v[i]))
    }
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> FreeModuleSignature<Set, Ring> for ConstFinitelyFreeModuleStructure<N, Set, SetB, Ring, RingB>
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
    use algebraeon_sets::sets::ConstSizeEnumeratedFiniteSetStructure;

    #[test]
    fn enumeration() {
        algebraeon_structures::assert_enumerated_ord_finite_set!(
            ConstFinitelyFreeModuleStructure::new(
                ConstSizeEnumeratedFiniteSetStructure::<3>::new(),
                QuaternaryField::structure(),
            ),
            64
        );
    }

    #[test]
    fn test_finite_rank_modules() {
        let m = ConstFinitelyFreeModuleStructure::new(
            ConstSizeEnumeratedFiniteSetStructure::<3>::new(),
            Integer::structure(),
        );

        let a = m.basis_element(0);
        let b = m.basis_element(1);
        let c = m.basis_element(2);

        assert_eq!(
            m.add(&m.neg(&b), &m.add(&a, &b)),
            [Integer::from(1), Integer::from(0), Integer::from(0)]
        );

        assert_eq!(
            m.add(&m.add(&a, &b), &m.add(&b, &c)),
            [Integer::from(1), Integer::from(2), Integer::from(1)]
        );

        assert_eq!(
            m.scalar_mul(&a, &5.into()),
            [Integer::from(5), Integer::from(0), Integer::from(0)]
        );

        assert_eq!(m.basis_vecs(), vec![a, b, c]);
    }
}
