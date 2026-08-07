//! Given two groups N and H and a left action ϕ : H -> Aut(N)
//! we get a semidirect product group G = N ⋊ H with group operation
//! `(n1, h1) * (n2, h2) = (n1 * ϕ_h1(n2), h1 * h2)`.
//! Conceptually `ϕ_h(n) * h = h * n` i.e. `ϕ_h` is the conjugation action `n -> h * n * h^-1`
//! N forms a normal subgroup of G and the quotient `G/N` is isomorphic to `H`.

use crate::*;
use std::{marker::PhantomData, sync::Arc};

#[derive(Debug, Clone)]
pub struct SemidirectProductElem<ElemN, ElemH> {
    // The composition order is important here. This is n*h not h*n.
    pub n: ElemN,
    pub h: ElemH,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SemidirectProductStructure<
    GroupN: GroupSignature,
    GroupH: GroupSignature,
    Phi: LeftGroupActionSignature<GroupH, GroupN>,
> {
    _n: PhantomData<GroupN>,
    _h: PhantomData<GroupH>,
    action: Arc<Phi>,
}

impl<GroupN: GroupSignature, GroupH: GroupSignature, Phi: LeftGroupActionSignature<GroupH, GroupN>>
    SemidirectProductStructure<GroupN, GroupH, Phi>
{
    pub fn new(action: Arc<Phi>) -> Arc<Self> {
        Self {
            _n: PhantomData,
            _h: PhantomData,
            action,
        }
        .into()
    }
}

impl<GroupN: GroupSignature, GroupH: GroupSignature, Phi: LeftGroupActionSignature<GroupH, GroupN>>
    SemidirectProductStructure<GroupN, GroupH, Phi>
{
    pub fn group_n(&self) -> Arc<GroupN> {
        self.action.set()
    }

    pub fn group_h(&self) -> Arc<GroupH> {
        self.action.group()
    }
}

impl<GroupN: GroupSignature, GroupH: GroupSignature, Phi: LeftGroupActionSignature<GroupH, GroupN>>
    Signature for SemidirectProductStructure<GroupN, GroupH, Phi>
{
}

impl<GroupN: GroupSignature, GroupH: GroupSignature, Phi: LeftGroupActionSignature<GroupH, GroupN>>
    SetSignature for SemidirectProductStructure<GroupN, GroupH, Phi>
{
    type Elem = SemidirectProductElem<GroupN::Elem, GroupH::Elem>;

    fn validate_element(self: &Arc<Self>, x: &Self::Elem) -> Result<(), String> {
        self.group_n().validate_element(&x.n)?;
        self.group_h().validate_element(&x.h)?;
        Ok(())
    }
}

impl<
    GroupN: GroupSignature + EqSignature,
    GroupH: GroupSignature + EqSignature,
    Phi: LeftGroupActionSignature<GroupH, GroupN>,
> EqSignature for SemidirectProductStructure<GroupN, GroupH, Phi>
{
    fn equal(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> bool {
        self.group_n().equal(&a.n, &b.n) && self.group_h().equal(&a.h, &b.h)
    }
}

impl<GroupN: GroupSignature, GroupH: GroupSignature, Phi: LeftGroupActionSignature<GroupH, GroupN>>
    IdentitySignature for SemidirectProductStructure<GroupN, GroupH, Phi>
{
    fn identity(self: &Arc<Self>) -> Self::Elem {
        SemidirectProductElem {
            n: self.group_n().identity(),
            h: self.group_h().identity(),
        }
    }
}

impl<GroupN: GroupSignature, GroupH: GroupSignature, Phi: LeftGroupActionSignature<GroupH, GroupN>>
    CompositionSignature for SemidirectProductStructure<GroupN, GroupH, Phi>
{
    fn compose(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        SemidirectProductElem {
            n: self.group_n().compose(&a.n, &self.action.apply(&a.h, &b.n)),
            h: self.group_h().compose(&a.h, &b.h),
        }
    }
}

impl<GroupN: GroupSignature, GroupH: GroupSignature, Phi: LeftGroupActionSignature<GroupH, GroupN>>
    AssociativeCompositionSignature for SemidirectProductStructure<GroupN, GroupH, Phi>
{
}

impl<GroupN: GroupSignature, GroupH: GroupSignature, Phi: LeftGroupActionSignature<GroupH, GroupN>>
    MonoidSignature for SemidirectProductStructure<GroupN, GroupH, Phi>
{
}

impl<GroupN: GroupSignature, GroupH: GroupSignature, Phi: LeftGroupActionSignature<GroupH, GroupN>>
    RightCancellativeCompositionSignature for SemidirectProductStructure<GroupN, GroupH, Phi>
{
    fn try_right_difference(
        self: &Arc<Self>,
        a: &Self::Elem,
        b: &Self::Elem,
    ) -> Option<Self::Elem> {
        Some(self.compose(a, &self.inverse(b)))
    }
}

impl<GroupN: GroupSignature, GroupH: GroupSignature, Phi: LeftGroupActionSignature<GroupH, GroupN>>
    LeftCancellativeCompositionSignature for SemidirectProductStructure<GroupN, GroupH, Phi>
{
    fn try_left_difference(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.compose(&self.inverse(b), a))
    }
}

impl<GroupN: GroupSignature, GroupH: GroupSignature, Phi: LeftGroupActionSignature<GroupH, GroupN>>
    TryRightInverseSignature for SemidirectProductStructure<GroupN, GroupH, Phi>
{
    fn try_right_inverse(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.inverse(a))
    }
}

impl<GroupN: GroupSignature, GroupH: GroupSignature, Phi: LeftGroupActionSignature<GroupH, GroupN>>
    TryLeftInverseSignature for SemidirectProductStructure<GroupN, GroupH, Phi>
{
    fn try_left_inverse(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.inverse(a))
    }
}

impl<GroupN: GroupSignature, GroupH: GroupSignature, Phi: LeftGroupActionSignature<GroupH, GroupN>>
    TryInverseSignature for SemidirectProductStructure<GroupN, GroupH, Phi>
{
    fn try_inverse(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.inverse(a))
    }
}

impl<GroupN: GroupSignature, GroupH: GroupSignature, Phi: LeftGroupActionSignature<GroupH, GroupN>>
    GroupSignature for SemidirectProductStructure<GroupN, GroupH, Phi>
{
    fn inverse(self: &Arc<Self>, a: &Self::Elem) -> Self::Elem {
        let h_inv = self.group_h().inverse(&a.h);
        SemidirectProductElem {
            n: self.action.apply(&h_inv, &self.group_n().inverse(&a.n)),
            h: h_inv,
        }
    }
}

impl<GroupN: GroupSignature, GroupH: GroupSignature, Phi: LeftGroupActionSignature<GroupH, GroupN>>
    SemidirectProductStructure<GroupN, GroupH, Phi>
{
    pub fn new_n(
        self: &Arc<Self>,
        n: &GroupN::Elem,
    ) -> SemidirectProductElem<GroupN::Elem, GroupH::Elem> {
        debug_assert!(self.group_n().is_element(n));
        SemidirectProductElem {
            n: n.clone(),
            h: self.group_h().identity(),
        }
    }

    pub fn new_h(
        self: &Arc<Self>,
        h: &GroupH::Elem,
    ) -> SemidirectProductElem<GroupN::Elem, GroupH::Elem> {
        debug_assert!(self.group_h().is_element(h));
        SemidirectProductElem {
            n: self.group_n().identity(),
            h: h.clone(),
        }
    }

    pub fn new_n_compose_h(
        self: &Arc<Self>,
        n: &GroupN::Elem,
        h: &GroupH::Elem,
    ) -> SemidirectProductElem<GroupN::Elem, GroupH::Elem> {
        debug_assert!(self.group_h().is_element(h));
        debug_assert!(self.group_n().is_element(n));
        SemidirectProductElem {
            n: n.clone(),
            h: h.clone(),
        }
    }

    pub fn new_h_compose_n(
        self: &Arc<Self>,
        h: &GroupH::Elem,
        n: &GroupN::Elem,
    ) -> SemidirectProductElem<GroupN::Elem, GroupH::Elem> {
        debug_assert!(self.group_h().is_element(h));
        debug_assert!(self.group_n().is_element(n));
        self.compose(&self.new_h(h), &self.new_n(n))
    }

    pub fn h_quotient_project(
        self: &Arc<Self>,
        g: &SemidirectProductElem<GroupN::Elem, GroupH::Elem>,
    ) -> GroupH::Elem {
        debug_assert!(self.is_element(g));
        g.h.clone()
    }

    pub fn n_compose_h(
        self: &Arc<Self>,
        g: &SemidirectProductElem<GroupN::Elem, GroupH::Elem>,
    ) -> (GroupN::Elem, GroupH::Elem) {
        debug_assert!(self.is_element(g));
        (g.n.clone(), g.h.clone())
    }

    pub fn h_compose_n(
        self: &Arc<Self>,
        g: &SemidirectProductElem<GroupN::Elem, GroupH::Elem>,
    ) -> (GroupH::Elem, GroupN::Elem) {
        debug_assert!(self.is_element(g));
        (g.h.clone(), self.action.apply_inverse(&g.h, &g.n))
    }
}
