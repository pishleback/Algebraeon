use crate::*;
use std::marker::PhantomData;

/// The grouplike structure obtained from another grouplike structure where the order of composition is reversed
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OppositeMagmaStructure<G: CompositionSignature, GB: BorrowedStructure<G>> {
    _magma: PhantomData<G>,
    magma: GB,
}

impl<G: CompositionSignature, GB: BorrowedStructure<G>> OppositeMagmaStructure<G, GB> {
    pub fn new(magma: GB) -> Self {
        Self {
            _magma: PhantomData,
            magma,
        }
    }
}

impl<G: CompositionSignature, GB: BorrowedStructure<G>> OppositeMagmaStructure<G, GB> {
    pub fn magma(&self) -> &G {
        self.magma.borrow()
    }
}

impl<G: GroupSignature, GB: BorrowedStructure<G>> OppositeMagmaStructure<G, GB> {
    pub fn group(&self) -> &G {
        self.magma()
    }
}

pub trait MagmaToOppositeSignature: CompositionSignature {
    fn opposite(&self) -> OppositeMagmaStructure<Self, &Self> {
        OppositeMagmaStructure::new(self)
    }

    fn into_opposite(self) -> OppositeMagmaStructure<Self, Self> {
        OppositeMagmaStructure::new(self)
    }
}
impl<G: CompositionSignature> MagmaToOppositeSignature for G {}

impl<G: CompositionSignature, GB: BorrowedStructure<G>> Signature
    for OppositeMagmaStructure<G, GB>
{
}

impl<G: CompositionSignature, GB: BorrowedStructure<G>> SetSignature
    for OppositeMagmaStructure<G, GB>
{
    type Elem = G::Elem;

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        self.magma().validate_element(x)
    }
}

impl<G: CompositionSignature + IdentitySignature, GB: BorrowedStructure<G>> IdentitySignature
    for OppositeMagmaStructure<G, GB>
{
    fn identity(&self) -> Self::Elem {
        self.magma().identity()
    }
}

impl<G: CompositionSignature, GB: BorrowedStructure<G>> CompositionSignature
    for OppositeMagmaStructure<G, GB>
{
    fn compose(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.magma().compose(b, a)
    }
}

impl<G: AssociativeCompositionSignature, GB: BorrowedStructure<G>> AssociativeCompositionSignature
    for OppositeMagmaStructure<G, GB>
{
}

impl<G: MonoidSignature, GB: BorrowedStructure<G>> MonoidSignature
    for OppositeMagmaStructure<G, GB>
{
}

impl<G: LeftCancellativeCompositionSignature, GB: BorrowedStructure<G>>
    RightCancellativeCompositionSignature for OppositeMagmaStructure<G, GB>
{
    fn try_right_difference(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        self.magma().try_left_difference(a, b)
    }
}

impl<G: RightCancellativeCompositionSignature, GB: BorrowedStructure<G>>
    LeftCancellativeCompositionSignature for OppositeMagmaStructure<G, GB>
{
    fn try_left_difference(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        self.magma().try_right_difference(a, b)
    }
}

impl<G: TryLeftInverseSignature, GB: BorrowedStructure<G>> TryRightInverseSignature
    for OppositeMagmaStructure<G, GB>
{
    fn try_right_inverse(&self, a: &Self::Elem) -> Option<Self::Elem> {
        self.magma().try_left_inverse(a)
    }
}

impl<G: TryRightInverseSignature, GB: BorrowedStructure<G>> TryLeftInverseSignature
    for OppositeMagmaStructure<G, GB>
{
    fn try_left_inverse(&self, a: &Self::Elem) -> Option<Self::Elem> {
        self.magma().try_right_inverse(a)
    }
}

impl<G: TryInverseSignature, GB: BorrowedStructure<G>> TryInverseSignature
    for OppositeMagmaStructure<G, GB>
{
    fn try_inverse(&self, a: &Self::Elem) -> Option<Self::Elem> {
        self.magma().try_inverse(a)
    }
}

impl<G: GroupSignature, GB: BorrowedStructure<G>> GroupSignature for OppositeMagmaStructure<G, GB> {
    fn inverse(&self, a: &Self::Elem) -> Self::Elem {
        self.group().inverse(a)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OppositeGroupAction<Action: Signature, ActionB: BorrowedStructure<Action>> {
    _action: PhantomData<Action>,
    action: ActionB,
}

impl<Action: Signature, ActionB: BorrowedStructure<Action>> OppositeGroupAction<Action, ActionB> {
    fn new(action: ActionB) -> Self {
        Self {
            _action: PhantomData,
            action,
        }
    }
}

impl<Action: Signature, ActionB: BorrowedStructure<Action>> Signature
    for OppositeGroupAction<Action, ActionB>
{
}

pub trait LeftActionToOppositeRightActionSignature<Group: GroupSignature, Set: SetSignature>:
    Signature
{
    fn opposite(&self) -> OppositeGroupAction<Self, &Self> {
        OppositeGroupAction::new(self)
    }

    fn into_opposite(self) -> OppositeGroupAction<Self, Self> {
        OppositeGroupAction::new(self)
    }
}
impl<Group: GroupSignature, Set: SetSignature, Action: LeftGroupActionSignature<Group, Set>>
    LeftActionToOppositeRightActionSignature<Group, Set> for Action
{
}

pub trait RightActionToOppositeLeftActionSignature<Set: SetSignature, Group: GroupSignature>:
    Signature
{
    fn opposite(&self) -> OppositeGroupAction<Self, &Self> {
        OppositeGroupAction::new(self)
    }

    fn into_opposite(self) -> OppositeGroupAction<Self, Self> {
        OppositeGroupAction::new(self)
    }
}
impl<Group: GroupSignature, Set: SetSignature, Action: RightGroupActionSignature<Set, Group>>
    RightActionToOppositeLeftActionSignature<Set, Group> for Action
{
}

impl<
    Group: GroupSignature,
    Set: SetSignature,
    Action: LeftGroupActionSignature<Group, Set>,
    ActionB: BorrowedStructure<Action>,
> RightGroupActionSignature<Set, Group> for OppositeGroupAction<Action, ActionB>
{
    fn group(&self) -> &Group {
        self.action.borrow().group()
    }

    fn set(&self) -> &Set {
        self.action.borrow().set()
    }

    fn right_apply(
        &self,
        g: &<Group>::Elem,
        x: &<Set as SetSignature>::Elem,
    ) -> <Set as SetSignature>::Elem {
        self.action.borrow().left_apply(g, x)
    }
}

impl<
    Group: GroupSignature,
    Set: SetSignature,
    Action: RightGroupActionSignature<Set, Group>,
    ActionB: BorrowedStructure<Action>,
> LeftGroupActionSignature<Group, Set> for OppositeGroupAction<Action, ActionB>
{
    fn group(&self) -> &Group {
        self.action.borrow().group()
    }

    fn set(&self) -> &Set {
        self.action.borrow().set()
    }

    fn left_apply(
        &self,
        g: &<Group>::Elem,
        x: &<Set as SetSignature>::Elem,
    ) -> <Set as SetSignature>::Elem {
        self.action.borrow().right_apply(g, x)
    }
}
