use super::group::Group;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum C2 {
    Ident,
    Flip,
}

impl Group for C2 {
    fn identity() -> Self {
        Self::Ident
    }

    fn inverse(self) -> Self {
        self
    }

    fn compose_mut(&mut self, other: &Self) {
        match (other) {
            C2::Ident => {}
            C2::Flip => match self {
                C2::Ident => {
                    *self = Self::Flip;
                }
                C2::Flip => *self = Self::Ident,
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c2() {
        debug_assert_eq!(C2::identity(), C2::Ident);
        debug_assert_eq!(C2::Ident.inverse(), C2::Ident);
        debug_assert_eq!(C2::Flip.inverse(), C2::Flip);
        debug_assert_eq!(C2::compose(C2::Ident, C2::Ident), C2::Ident);
        debug_assert_eq!(C2::compose(C2::Flip, C2::Ident), C2::Flip);
        debug_assert_eq!(C2::compose(C2::Ident, C2::Flip), C2::Flip);
        debug_assert_eq!(C2::compose(C2::Flip, C2::Flip), C2::Ident);
    }
}
