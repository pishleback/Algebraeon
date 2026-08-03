use crate::{
    linear::finitely_free_module::FinitelyFreeModuleStructure, matrix::Matrix, structure::*,
};
use algebraeon_sets::sets::{ConstSizeFunctionsStructure, Function};
use algebraeon_structures::*;
use std::{borrow::Cow, cmp::Ordering, rc::Rc};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConstFinitelyFreeModuleStructure<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature,
> {
    functions: Rc<ConstSizeFunctionsStructure<N, Set, Ring>>,
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature,
> ConstFinitelyFreeModuleStructure<N, Set, Ring>
{
    pub fn new(set: Rc<Set>, ring: Rc<Ring>) -> Rc<Self> {
        Self {
            functions: ConstSizeFunctionsStructure::new(set, ring),
        }
        .into()
    }

    pub fn forget_const(&self) -> Rc<FinitelyFreeModuleStructure<Set, Ring>> {
        FinitelyFreeModuleStructure::new(
            self.functions.domain().clone(),
            self.functions.range().clone(),
        )
    }
}
pub trait RingToConstFinitelyFreeModuleSignature: RingSignature {
    fn free_module<
        const N: usize,
        Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    >(
        self: &Rc<Self>,
        set: Rc<Set>,
    ) -> Rc<ConstFinitelyFreeModuleStructure<N, Set, Self>> {
        ConstFinitelyFreeModuleStructure::new(set, self.clone())
    }
}
impl<Ring: RingSignature> RingToConstFinitelyFreeModuleSignature for Ring {}

pub trait SetToConstFinitelyFreeModuleSignature<const N: usize>:
    ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature
{
    fn free_module<Ring: RingSignature>(
        self: &Rc<Self>,
        ring: Rc<Ring>,
    ) -> Rc<ConstFinitelyFreeModuleStructure<N, Self, Ring>> {
        ConstFinitelyFreeModuleStructure::new(self.clone(), ring)
    }
}
impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature>
    SetToConstFinitelyFreeModuleSignature<N> for Set
{
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature,
> ConstFinitelyFreeModuleStructure<N, Set, Ring>
{
    pub fn set(&self) -> &Rc<Set> {
        self.functions.domain()
    }

    pub fn ring(&self) -> &Rc<Ring> {
        self.functions.range()
    }

    pub fn functions_restructure(&self) -> Rc<ConstSizeFunctionsStructure<N, Set, Ring>> {
        self.functions.clone()
    }

    pub fn to_col(&self, v: &<Self as SetSignature>::Elem) -> Matrix<Ring::Elem> {
        self.forget_const().to_col(&v.into())
    }

    pub fn to_row(&self, v: &<Self as SetSignature>::Elem) -> Matrix<Ring::Elem> {
        self.forget_const().to_row(&v.into())
    }

    pub fn from_row(self: &Rc<Self>, m: &Matrix<Ring::Elem>) -> <Self as SetSignature>::Elem {
        self.forget_const().from_row(m).try_into().unwrap()
    }

    pub fn from_col(self: &Rc<Self>, m: &Matrix<Ring::Elem>) -> <Self as SetSignature>::Elem {
        self.forget_const().from_col(m).try_into().unwrap()
    }

    pub fn basis_element(&self, i: usize) -> <Self as SetSignature>::Elem {
        self.forget_const().basis_element(i).try_into().unwrap()
    }
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature,
> Signature for ConstFinitelyFreeModuleStructure<N, Set, Ring>
{
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature,
> SetSignature for ConstFinitelyFreeModuleStructure<N, Set, Ring>
{
    type Elem = Function<N, Set::Elem, Ring::Elem>;

    fn validate_element(self: &Rc<Self>, v: &Self::Elem) -> Result<(), String> {
        self.functions_restructure().validate_element(v)
    }
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + EqSignature,
> EqSignature for ConstFinitelyFreeModuleStructure<N, Set, Ring>
{
    fn equal(self: &Rc<Self>, v: &Self::Elem, w: &Self::Elem) -> bool {
        self.functions_restructure().equal(v, w)
    }
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + OrdSignature,
> PartialOrdSignature for ConstFinitelyFreeModuleStructure<N, Set, Ring>
{
    fn partial_cmp(self: &Rc<Self>, v: &Self::Elem, w: &Self::Elem) -> Option<Ordering> {
        self.functions_restructure().partial_cmp(v, w)
    }
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + OrdSignature,
> OrdSignature for ConstFinitelyFreeModuleStructure<N, Set, Ring>
{
    fn cmp(self: &Rc<Self>, v: &Self::Elem, w: &Self::Elem) -> Ordering {
        self.functions_restructure().cmp(v, w)
    }
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + OrderedFiniteSetSignature,
> CountableSetSignature for ConstFinitelyFreeModuleStructure<N, Set, Ring>
{
    fn generate_all_elements(self: Rc<Self>) -> impl Iterator<Item = Self::Elem> {
        self.functions_restructure().generate_all_elements()
    }
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + OrderedFiniteSetSignature,
> FiniteSetSignature for ConstFinitelyFreeModuleStructure<N, Set, Ring>
{
    fn list_all_elements(self: &Rc<Self>) -> Vec<Self::Elem> {
        self.functions_restructure().list_all_elements()
    }

    fn size(self: &Rc<Self>) -> Natural {
        self.functions_restructure().size()
    }

    fn generate_random_elements(self: Rc<Self>, seed: u64) -> impl Iterator<Item = Self::Elem> {
        self.functions_restructure()
            .clone()
            .generate_random_elements(seed)
    }
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + OrderedFiniteSetSignature,
> OrderedFiniteSetSignature for ConstFinitelyFreeModuleStructure<N, Set, Ring>
{
    fn list_all_elements_ordered(self: &Rc<Self>) -> Vec<Self::Elem> {
        self.functions_restructure().list_all_elements_ordered()
    }

    fn element_to_enumeration(self: &Rc<Self>, elem: &Self::Elem) -> Natural {
        self.functions_restructure().element_to_enumeration(elem)
    }

    fn enumeration_to_element(self: &Rc<Self>, num: &Natural) -> Option<Self::Elem> {
        self.functions_restructure().enumeration_to_element(num)
    }
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature,
> RinglikeSpecializationSignature for ConstFinitelyFreeModuleStructure<N, Set, Ring>
{
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature,
> ZeroSignature for ConstFinitelyFreeModuleStructure<N, Set, Ring>
{
    fn zero(self: &Rc<Self>) -> Self::Elem {
        std::array::from_fn(|_| self.ring().zero()).into()
    }
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature,
> AdditionSignature for ConstFinitelyFreeModuleStructure<N, Set, Ring>
{
    fn add(self: &Rc<Self>, v: &Self::Elem, w: &Self::Elem) -> Self::Elem {
        debug_assert!(self.validate_element(v).is_ok());
        debug_assert!(self.validate_element(w).is_ok());
        std::array::from_fn(|i| self.ring().add(&v[i], &w[i])).into()
    }
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature,
> CancellativeAdditionSignature for ConstFinitelyFreeModuleStructure<N, Set, Ring>
{
    fn try_sub(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.sub(a, b))
    }
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature,
> TryNegateSignature for ConstFinitelyFreeModuleStructure<N, Set, Ring>
{
    fn try_neg(self: &Rc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.neg(a))
    }
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature,
> AdditiveMonoidSignature for ConstFinitelyFreeModuleStructure<N, Set, Ring>
{
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature,
> AdditiveGroupSignature for ConstFinitelyFreeModuleStructure<N, Set, Ring>
{
    fn neg(self: &Rc<Self>, v: &Self::Elem) -> Self::Elem {
        debug_assert!(self.validate_element(v).is_ok());
        std::array::from_fn(|i| self.ring().neg(&v[i])).into()
    }

    fn sub(self: &Rc<Self>, v: &Self::Elem, w: &Self::Elem) -> Self::Elem {
        debug_assert!(self.validate_element(v).is_ok());
        debug_assert!(self.validate_element(w).is_ok());
        std::array::from_fn(|i| self.ring().sub(&v[i], &w[i])).into()
    }
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature,
> SemiModuleSignature<Ring> for ConstFinitelyFreeModuleStructure<N, Set, Ring>
{
    fn ring(self: &Rc<Self>) -> Rc<Ring> {
        self.functions.range().clone()
    }

    fn scalar_mul(self: &Rc<Self>, v: &Self::Elem, r: &Ring::Elem) -> Self::Elem {
        debug_assert!(self.validate_element(v).is_ok());
        std::array::from_fn(|i| self.ring().mul(r, &v[i])).into()
    }
}

impl<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature,
> FreeModuleSignature<Set, Ring> for ConstFinitelyFreeModuleStructure<N, Set, Ring>
{
    fn basis_set(self: &Rc<Self>) -> Rc<Set> {
        self.set().clone()
    }

    fn to_component<'a>(self: &Rc<Self>, b: &Set::Elem, v: &'a Self::Elem) -> Cow<'a, Ring::Elem> {
        let b: usize = self.set().element_to_enumeration(b).try_into().unwrap();
        debug_assert!(b < self.rank());
        Cow::Borrowed(&v[b])
    }

    fn from_component(self: &Rc<Self>, b: &Set::Elem, r: &<Ring>::Elem) -> Self::Elem {
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
            [Integer::from(1), Integer::from(0), Integer::from(0)].into()
        );

        assert_eq!(
            m.add(&m.add(&a, &b), &m.add(&b, &c)),
            [Integer::from(1), Integer::from(2), Integer::from(1)].into()
        );

        assert_eq!(
            m.scalar_mul(&a, &5.into()),
            [Integer::from(5), Integer::from(0), Integer::from(0)].into()
        );

        assert_eq!(m.basis_vecs(), vec![a, b, c]);
    }
}
