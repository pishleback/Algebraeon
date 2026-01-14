use algebraeon_macros::CanonicalStructure;
use algebraeon_sets::structure::{EqSignature, MetaType, SetSignature, Signature};

use crate::structure::{
    AssociativeCompositionSignature, CompositionSignature, GroupSignature, IdentitySignature,
    LeftCancellativeCompositionSignature, MonoidSignature, RightCancellativeCompositionSignature,
    TryInverseSignature, TryLeftInverseSignature, TryRightInverseSignature,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, CanonicalStructure)]
#[canonical_structure(eq)]
pub enum C2 {
    Identity,
    Flip,
}

impl CompositionSignature for C2CanonicalStructure {
    fn compose(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        match (a, b) {
            (C2::Identity, C2::Identity) => C2::Identity,
            (C2::Identity, C2::Flip) => C2::Flip,
            (C2::Flip, C2::Identity) => C2::Flip,
            (C2::Flip, C2::Flip) => C2::Identity,
        }
    }
}

impl AssociativeCompositionSignature for C2CanonicalStructure {}

impl LeftCancellativeCompositionSignature for C2CanonicalStructure {
    fn try_left_difference(&self, a: &Self::Set, b: &Self::Set) -> Option<Self::Set> {
        Some(self.compose(&self.inverse(b), a))
    }
}

impl RightCancellativeCompositionSignature for C2CanonicalStructure {
    fn try_right_difference(&self, a: &Self::Set, b: &Self::Set) -> Option<Self::Set> {
        Some(self.compose(a, &self.inverse(b)))
    }
}

impl IdentitySignature for C2CanonicalStructure {
    fn identity(&self) -> Self::Set {
        C2::Identity
    }
}

impl MonoidSignature for C2CanonicalStructure {}

impl TryLeftInverseSignature for C2CanonicalStructure {
    fn try_left_inverse(&self, a: &Self::Set) -> Option<Self::Set> {
        Some(self.inverse(a))
    }
}

impl TryRightInverseSignature for C2CanonicalStructure {
    fn try_right_inverse(&self, a: &Self::Set) -> Option<Self::Set> {
        Some(self.inverse(a))
    }
}

impl TryInverseSignature for C2CanonicalStructure {
    fn try_inverse(&self, a: &Self::Set) -> Option<Self::Set> {
        Some(self.inverse(a))
    }
}

impl GroupSignature for C2CanonicalStructure {
    fn inverse(&self, a: &Self::Set) -> Self::Set {
        *a
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::{MetaCompositionSignature, MetaGroupSignature, MetaIdentitySignature};

    #[test]
    fn test_c2() {
        debug_assert_eq!(C2::identity(), C2::Identity);
        debug_assert_eq!(C2::Identity.inverse(), C2::Identity);
        debug_assert_eq!(C2::Flip.inverse(), C2::Flip);
        debug_assert_eq!(C2::compose(&C2::Identity, &C2::Identity), C2::Identity);
        debug_assert_eq!(C2::compose(&C2::Flip, &C2::Identity), C2::Flip);
        debug_assert_eq!(C2::compose(&C2::Identity, &C2::Flip), C2::Flip);
        debug_assert_eq!(C2::compose(&C2::Flip, &C2::Flip), C2::Identity);
    }
}
