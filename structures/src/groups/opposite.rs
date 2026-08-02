use crate::*;
use std::sync::Arc;

/// The grouplike structure obtained from another grouplike structure where the order of composition is reversed
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OppositeMagmaStructure<G: CompositionSignature> {
    magma: Arc<G>,
}

impl<G: CompositionSignature> OppositeMagmaStructure<G> {
    pub fn new(magma: Arc<G>) -> Self {
        Self { magma }
    }
}

impl<G: CompositionSignature> OppositeMagmaStructure<G> {
    pub fn magma(self: &Arc<Self>) -> &Arc<G> {
        &self.magma
    }
}

impl<G: GroupSignature> OppositeMagmaStructure<G> {
    pub fn group(self: &Arc<Self>) -> &Arc<G> {
        self.magma()
    }
}

pub trait MagmaToOppositeSignature: CompositionSignature {
    fn opposite(self: &Arc<Self>) -> OppositeMagmaStructure<Self> {
        OppositeMagmaStructure::new(self.clone())
    }
}
impl<G: CompositionSignature> MagmaToOppositeSignature for G {}

impl<G: CompositionSignature> Signature for OppositeMagmaStructure<G> {}

impl<G: CompositionSignature> SetSignature for OppositeMagmaStructure<G> {
    type Elem = G::Elem;

    fn validate_element(self: &Arc<Self>, x: &Self::Elem) -> Result<(), String> {
        self.magma().validate_element(x)
    }
}

impl<G: CompositionSignature + IdentitySignature> IdentitySignature for OppositeMagmaStructure<G> {
    fn identity(self: &Arc<Self>) -> Self::Elem {
        self.magma().identity()
    }
}

impl<G: CompositionSignature> CompositionSignature for OppositeMagmaStructure<G> {
    fn compose(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
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
    fn try_right_difference(
        self: &Arc<Self>,
        a: &Self::Elem,
        b: &Self::Elem,
    ) -> Option<Self::Elem> {
        self.magma().try_left_difference(a, b)
    }
}

impl<G: RightCancellativeCompositionSignature> LeftCancellativeCompositionSignature
    for OppositeMagmaStructure<G>
{
    fn try_left_difference(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        self.magma().try_right_difference(a, b)
    }
}

impl<G: TryLeftInverseSignature> TryRightInverseSignature for OppositeMagmaStructure<G> {
    fn try_right_inverse(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        self.magma().try_left_inverse(a)
    }
}

impl<G: TryRightInverseSignature> TryLeftInverseSignature for OppositeMagmaStructure<G> {
    fn try_left_inverse(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        self.magma().try_right_inverse(a)
    }
}

impl<G: TryInverseSignature> TryInverseSignature for OppositeMagmaStructure<G> {
    fn try_inverse(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        self.magma().try_inverse(a)
    }
}

impl<G: GroupSignature> GroupSignature for OppositeMagmaStructure<G> {
    fn inverse(self: &Arc<Self>, a: &Self::Elem) -> Self::Elem {
        self.group().inverse(a)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OppositeGroupAction<Action: Signature> {
    action: Arc<Action>,
}

impl<Action: Signature> OppositeGroupAction<Action> {
    pub fn new(action: Arc<Action>) -> Arc<Self> {
        Self { action }.into()
    }

    pub fn action(self: &Arc<Self>) -> &Arc<Action> {
        &self.action
    }
}

impl<Action: Signature> Signature for OppositeGroupAction<Action> {}

pub trait LeftActionToOppositeRightActionSignature<Group: GroupSignature, Set: SetSignature>:
    Signature
{
    fn opposite(self: &Arc<Self>) -> Arc<OppositeGroupAction<Self>> {
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
    fn opposite(self: &Arc<Self>) -> Arc<OppositeGroupAction<Self>> {
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
    fn group(self: &Arc<Self>) -> &Arc<Group> {
        self.action().group()
    }

    fn set(self: &Arc<Self>) -> &Arc<Set> {
        self.action().set()
    }

    fn apply(
        self: &Arc<Self>,
        g: &<Group>::Elem,
        x: &<Set as SetSignature>::Elem,
    ) -> <Set as SetSignature>::Elem {
        self.action().apply(&self.group().inverse(g), x)
    }
}

impl<Group: GroupSignature, Set: SetSignature, Action: RightGroupActionSignature<Set, Group>>
    LeftGroupActionSignature<Group, Set> for OppositeGroupAction<Action>
{
    fn group(self: &Arc<Self>) -> &Arc<Group> {
        self.action().group()
    }

    fn set(self: &Arc<Self>) -> &Arc<Set> {
        self.action().set()
    }

    fn apply(
        self: &Arc<Self>,
        g: &<Group>::Elem,
        x: &<Set as SetSignature>::Elem,
    ) -> <Set as SetSignature>::Elem {
        self.action().apply(&self.group().inverse(g), x)
    }
}
