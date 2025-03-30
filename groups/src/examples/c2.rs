use crate::structure::Group;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum C2 {
    Identity,
    Flip,
}

impl Group for C2 {
    fn identity() -> Self {
        Self::Identity
    }

    fn inverse(self) -> Self {
        self
    }

    fn compose_mut(&mut self, other: &Self) {
        match other {
            C2::Identity => {}
            C2::Flip => match self {
                C2::Identity => {
                    *self = Self::Flip;
                }
                C2::Flip => *self = Self::Identity,
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c2() {
        debug_assert_eq!(C2::identity(), C2::Identity);
        debug_assert_eq!(C2::Identity.inverse(), C2::Identity);
        debug_assert_eq!(C2::Flip.inverse(), C2::Flip);
        debug_assert_eq!(C2::compose(C2::Identity, C2::Identity), C2::Identity);
        debug_assert_eq!(C2::compose(C2::Flip, C2::Identity), C2::Flip);
        debug_assert_eq!(C2::compose(C2::Identity, C2::Flip), C2::Flip);
        debug_assert_eq!(C2::compose(C2::Flip, C2::Flip), C2::Identity);
    }
}
