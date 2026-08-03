use crate::{matrix::Matrix, structure::*};
use algebraeon_sets::sets::FunctionsStructure;
use algebraeon_structures::*;
use std::{borrow::Cow, cmp::Ordering, rc::Rc};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FinitelyFreeModuleStructure<Set: OrderedFiniteSetSignature, Ring: RingSignature> {
    functions: Rc<FunctionsStructure<Set, Ring>>,
}

impl<Set: OrderedFiniteSetSignature, Ring: RingSignature> FinitelyFreeModuleStructure<Set, Ring> {
    pub fn new(set: Rc<Set>, ring: Rc<Ring>) -> Rc<Self> {
        Self {
            functions: FunctionsStructure::new(set, ring),
        }
        .into()
    }
}
pub trait RingToFinitelyFreeModuleSignature: RingSignature {
    fn free_module<Set: OrderedFiniteSetSignature>(
        self: &Rc<Self>,
        set: Rc<Set>,
    ) -> Rc<FinitelyFreeModuleStructure<Set, Self>> {
        FinitelyFreeModuleStructure::new(set, self.clone())
    }
}
impl<Ring: RingSignature> RingToFinitelyFreeModuleSignature for Ring {}

pub trait SetToFinitelyFreeModuleSignature: OrderedFiniteSetSignature {
    fn free_module<Ring: RingSignature>(
        self: &Rc<Self>,
        ring: Rc<Ring>,
    ) -> Rc<FinitelyFreeModuleStructure<Self, Ring>> {
        FinitelyFreeModuleStructure::new(self.clone(), ring)
    }
}
impl<Set: OrderedFiniteSetSignature> SetToFinitelyFreeModuleSignature for Set {}

impl<Set: OrderedFiniteSetSignature, Ring: RingSignature> FinitelyFreeModuleStructure<Set, Ring> {
    pub fn set(&self) -> Rc<Set> {
        self.functions.domain().clone()
    }

    pub fn ring(&self) -> Rc<Ring> {
        self.functions.range().clone()
    }

    pub fn functions_restructure(&self) -> &Rc<FunctionsStructure<Set, Ring>> {
        &self.functions
    }

    pub fn to_col(self: &Rc<Self>, v: &<Self as SetSignature>::Elem) -> Matrix<Ring::Elem> {
        debug_assert!(self.validate_element(v).is_ok());
        Matrix::construct(self.rank(), 1, |r, _| v[r].clone())
    }

    pub fn to_row(self: &Rc<Self>, v: &<Self as SetSignature>::Elem) -> Matrix<Ring::Elem> {
        debug_assert!(self.validate_element(v).is_ok());
        Matrix::construct(1, self.rank(), |_, c| v[c].clone())
    }

    pub fn from_row(self: &Rc<Self>, m: &Matrix<Ring::Elem>) -> <Self as SetSignature>::Elem {
        debug_assert_eq!(m.rows(), 1);
        debug_assert_eq!(m.cols(), self.rank());
        (0..self.rank())
            .map(|i| m.at(0, i).unwrap().clone())
            .collect()
    }

    pub fn from_col(self: &Rc<Self>, m: &Matrix<Ring::Elem>) -> <Self as SetSignature>::Elem {
        debug_assert_eq!(m.cols(), 1);
        debug_assert_eq!(m.rows(), self.rank());
        (0..self.rank())
            .map(|i| m.at(i, 0).unwrap().clone())
            .collect()
    }

    pub fn basis_element(self: &Rc<Self>, i: usize) -> <Self as SetSignature>::Elem {
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

impl<Set: OrderedFiniteSetSignature, Ring: RingSignature> Signature
    for FinitelyFreeModuleStructure<Set, Ring>
{
}

impl<Set: OrderedFiniteSetSignature, Ring: RingSignature> SetSignature
    for FinitelyFreeModuleStructure<Set, Ring>
{
    type Elem = Vec<Ring::Elem>;

    fn validate_element(self: &Rc<Self>, v: &Self::Elem) -> Result<(), String> {
        self.functions_restructure().validate_element(v)
    }
}

impl<Set: OrderedFiniteSetSignature, Ring: RingSignature + EqSignature> EqSignature
    for FinitelyFreeModuleStructure<Set, Ring>
{
    fn equal(self: &Rc<Self>, v: &Self::Elem, w: &Self::Elem) -> bool {
        self.functions_restructure().equal(v, w)
    }
}

impl<Set: OrderedFiniteSetSignature, Ring: RingSignature + OrdSignature> PartialOrdSignature
    for FinitelyFreeModuleStructure<Set, Ring>
{
    fn partial_cmp(self: &Rc<Self>, v: &Self::Elem, w: &Self::Elem) -> Option<Ordering> {
        self.functions_restructure().partial_cmp(v, w)
    }
}

impl<Set: OrderedFiniteSetSignature, Ring: RingSignature + OrdSignature> OrdSignature
    for FinitelyFreeModuleStructure<Set, Ring>
{
    fn cmp(self: &Rc<Self>, v: &Self::Elem, w: &Self::Elem) -> Ordering {
        self.functions_restructure().cmp(v, w)
    }
}

impl<Set: OrderedFiniteSetSignature, Ring: RingSignature + OrderedFiniteSetSignature>
    CountableSetSignature for FinitelyFreeModuleStructure<Set, Ring>
{
    fn generate_all_elements(self: Rc<Self>) -> impl Iterator<Item = Self::Elem> {
        self.functions_restructure().clone().generate_all_elements()
    }
}

impl<Set: OrderedFiniteSetSignature, Ring: RingSignature + OrderedFiniteSetSignature>
    FiniteSetSignature for FinitelyFreeModuleStructure<Set, Ring>
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

impl<Set: OrderedFiniteSetSignature, Ring: RingSignature + OrderedFiniteSetSignature>
    OrderedFiniteSetSignature for FinitelyFreeModuleStructure<Set, Ring>
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

impl<Set: OrderedFiniteSetSignature, Ring: RingSignature> RinglikeSpecializationSignature
    for FinitelyFreeModuleStructure<Set, Ring>
{
}

impl<Set: OrderedFiniteSetSignature, Ring: RingSignature> ZeroSignature
    for FinitelyFreeModuleStructure<Set, Ring>
{
    fn zero(self: &Rc<Self>) -> Self::Elem {
        (0..self.rank()).map(|_| self.ring().zero()).collect()
    }
}

impl<Set: OrderedFiniteSetSignature, Ring: RingSignature> AdditionSignature
    for FinitelyFreeModuleStructure<Set, Ring>
{
    fn add(self: &Rc<Self>, v: &Self::Elem, w: &Self::Elem) -> Self::Elem {
        debug_assert!(self.validate_element(v).is_ok());
        debug_assert!(self.validate_element(w).is_ok());
        (0..self.rank())
            .map(|i| self.ring().add(&v[i], &w[i]))
            .collect()
    }
}

impl<Set: OrderedFiniteSetSignature, Ring: RingSignature> CancellativeAdditionSignature
    for FinitelyFreeModuleStructure<Set, Ring>
{
    fn try_sub(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.sub(a, b))
    }
}

impl<Set: OrderedFiniteSetSignature, Ring: RingSignature> TryNegateSignature
    for FinitelyFreeModuleStructure<Set, Ring>
{
    fn try_neg(self: &Rc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.neg(a))
    }
}

impl<Set: OrderedFiniteSetSignature, Ring: RingSignature> AdditiveMonoidSignature
    for FinitelyFreeModuleStructure<Set, Ring>
{
}

impl<Set: OrderedFiniteSetSignature, Ring: RingSignature> AdditiveGroupSignature
    for FinitelyFreeModuleStructure<Set, Ring>
{
    fn neg(self: &Rc<Self>, v: &Self::Elem) -> Self::Elem {
        debug_assert!(self.validate_element(v).is_ok());
        v.iter().map(|r| self.ring().neg(r)).collect()
    }

    fn sub(self: &Rc<Self>, v: &Self::Elem, w: &Self::Elem) -> Self::Elem {
        debug_assert!(self.validate_element(v).is_ok());
        debug_assert!(self.validate_element(w).is_ok());
        (0..self.rank())
            .map(|i| self.ring().sub(&v[i], &w[i]))
            .collect()
    }
}

impl<Set: OrderedFiniteSetSignature, Ring: RingSignature> SemiModuleSignature<Ring>
    for FinitelyFreeModuleStructure<Set, Ring>
{
    fn ring(self: &Rc<Self>) -> Rc<Ring> {
        self.functions.range().clone()
    }

    fn scalar_mul(self: &Rc<Self>, v: &Self::Elem, r: &Ring::Elem) -> Self::Elem {
        debug_assert!(self.validate_element(v).is_ok());
        v.iter().map(|s| self.ring().mul(r, s)).collect()
    }
}

impl<Set: OrderedFiniteSetSignature, Ring: RingSignature> FreeModuleSignature<Set, Ring>
    for FinitelyFreeModuleStructure<Set, Ring>
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
