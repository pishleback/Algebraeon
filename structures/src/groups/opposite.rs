use crate::*;
use std::rc::Rc;

/// The grouplike structure obtained from another grouplike structure where the order of composition is reversed
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OppositeMagmaStructure<G: CompositionSignature> {
    magma: Rc<G>,
}

impl<G: CompositionSignature> OppositeMagmaStructure<G> {
    pub fn new(magma: Rc<G>) -> Self {
        Self { magma }
    }
}

impl<G: CompositionSignature> OppositeMagmaStructure<G> {
    pub fn magma(self: &Rc<Self>) -> &Rc<G> {
        &self.magma
    }
}

impl<G: GroupSignature> OppositeMagmaStructure<G> {
    pub fn group(self: &Rc<Self>) -> &Rc<G> {
        self.magma()
    }
}

pub trait MagmaToOppositeSignature: CompositionSignature {
    fn opposite(self: &Rc<Self>) -> OppositeMagmaStructure<Self> {
        OppositeMagmaStructure::new(self.clone())
    }
}
impl<G: CompositionSignature> MagmaToOppositeSignature for G {}

impl<G: CompositionSignature> Signature for OppositeMagmaStructure<G> {}

impl<G: CompositionSignature> SetSignature for OppositeMagmaStructure<G> {
    type Elem = G::Elem;

    fn validate_element(self: &Rc<Self>, x: &Self::Elem) -> Result<(), String> {
        self.magma().validate_element(x)
    }
}

impl<G: CompositionSignature + IdentitySignature> IdentitySignature for OppositeMagmaStructure<G> {
    fn identity(self: &Rc<Self>) -> Self::Elem {
        self.magma().identity()
    }
}

impl<G: CompositionSignature> CompositionSignature for OppositeMagmaStructure<G> {
    fn compose(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.magma().compose(b, a)
    }
}

impl<G: AssociativeCompositionSignature> AssociativeCompositionSignature
    for OppositeMagmaStructure<G>
{
}

impl<G: MonoidSignature> MonoidSignature for OppositeMagmaStructure<G> {}

impl<G: LeftCancellativeCompositionSignature> RightCancellativeCompositionSignature
    for OppositeMagmaStructure<G>
{
    fn try_right_difference(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        self.magma().try_left_difference(a, b)
    }
}

impl<G: RightCancellativeCompositionSignature> LeftCancellativeCompositionSignature
    for OppositeMagmaStructure<G>
{
    fn try_left_difference(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        self.magma().try_right_difference(a, b)
    }
}

impl<G: TryLeftInverseSignature> TryRightInverseSignature for OppositeMagmaStructure<G> {
    fn try_right_inverse(self: &Rc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        self.magma().try_left_inverse(a)
    }
}

impl<G: TryRightInverseSignature> TryLeftInverseSignature for OppositeMagmaStructure<G> {
    fn try_left_inverse(self: &Rc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        self.magma().try_right_inverse(a)
    }
}

impl<G: TryInverseSignature> TryInverseSignature for OppositeMagmaStructure<G> {
    fn try_inverse(self: &Rc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        self.magma().try_inverse(a)
    }
}

impl<G: GroupSignature> GroupSignature for OppositeMagmaStructure<G> {
    fn inverse(self: &Rc<Self>, a: &Self::Elem) -> Self::Elem {
        self.group().inverse(a)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OppositeGroupAction<Action: Signature> {
    action: Rc<Action>,
}

impl<Action: Signature> OppositeGroupAction<Action> {
    pub fn new(action: Rc<Action>) -> Rc<Self> {
        Self { action }.into()
    }

    pub fn action(self: &Rc<Self>) -> &Rc<Action> {
        &self.action
    }
}

impl<Action: Signature> Signature for OppositeGroupAction<Action> {}

pub trait LeftActionToOppositeRightActionSignature<Group: GroupSignature, Set: SetSignature>:
    Signature
{
    fn opposite(self: &Rc<Self>) -> Rc<OppositeGroupAction<Self>> {
        OppositeGroupAction::new(self.clone())
    }
}
impl<Group: GroupSignature, Set: SetSignature, Action: LeftGroupActionSignature<Group, Set>>
    LeftActionToOppositeRightActionSignature<Group, Set> for Action
{
}

pub trait RightActionToOppositeLeftActionSignature<Set: SetSignature, Group: GroupSignature>:
    Signature
{
    fn opposite(self: &Rc<Self>) -> Rc<OppositeGroupAction<Self>> {
        OppositeGroupAction::new(self.clone())
    }
}
impl<Group: GroupSignature, Set: SetSignature, Action: RightGroupActionSignature<Set, Group>>
    RightActionToOppositeLeftActionSignature<Set, Group> for Action
{
}

impl<Group: GroupSignature, Set: SetSignature, Action: LeftGroupActionSignature<Group, Set>>
    RightGroupActionSignature<Set, Group> for OppositeGroupAction<Action>
{
    fn group(self: &Rc<Self>) -> &Rc<Group> {
        self.action().group()
    }

    fn set(self: &Rc<Self>) -> &Rc<Set> {
        self.action().set()
    }

    fn apply(
        self: &Rc<Self>,
        g: &<Group>::Elem,
        x: &<Set as SetSignature>::Elem,
    ) -> <Set as SetSignature>::Elem {
        self.action().apply(&self.group().inverse(g), x)
    }
}

impl<Group: GroupSignature, Set: SetSignature, Action: RightGroupActionSignature<Set, Group>>
    LeftGroupActionSignature<Group, Set> for OppositeGroupAction<Action>
{
    fn group(self: &Rc<Self>) -> &Rc<Group> {
        self.action().group()
    }

    fn set(self: &Rc<Self>) -> &Rc<Set> {
        self.action().set()
    }

    fn apply(
        self: &Rc<Self>,
        g: &<Group>::Elem,
        x: &<Set as SetSignature>::Elem,
    ) -> <Set as SetSignature>::Elem {
        self.action().apply(&self.group().inverse(g), x)
    }
}
