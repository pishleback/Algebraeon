use crate::*;
use algebraeon_macros::signature_meta_trait;
use std::borrow::Borrow;
use std::collections::HashSet;

/// A set with a binary operation of composition.
#[signature_meta_trait]
pub trait CompositionSignature: SetSignature {
    fn compose(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem;
    fn compose_mut(&self, a: &mut Self::Elem, b: &Self::Elem) {
        *a = self.compose(a, b);
    }
}

/// When composition is associative.
#[signature_meta_trait]
pub trait AssociativeCompositionSignature: CompositionSignature {
    /// Returns `None` if the list is empty.
    fn compose_nonempty_list(&self, mut elems: Vec<impl Borrow<Self::Elem>>) -> Option<Self::Elem> {
        let mut total = elems.pop()?.borrow().clone();
        for elem in elems {
            total = self.compose(&total, elem.borrow());
        }
        Some(total)
    }
}

/// When composition is commutative.
#[signature_meta_trait]
pub trait CommutativeCompositionSignature: CompositionSignature {}

/// When `compose(a, x)` = `compose(a, y)` implies `x` = `y` for all `a`, `x`, `y`.
#[signature_meta_trait]
pub trait LeftCancellativeCompositionSignature: CompositionSignature {
    /// Try to find `x` such that `a` = `compose(b, x)`.
    fn try_left_difference(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem>;
}

/// When `compose(x, a)` = `compose(y, a)` implies `x` = `y` for all `a`, `x`, `y`.
#[signature_meta_trait]
pub trait RightCancellativeCompositionSignature: CompositionSignature {
    /// Try to find `x` such that `a` = `compose(x, b)`.
    fn try_right_difference(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem>;
}

#[signature_meta_trait]
pub trait CancellativeCompositionSignature:
    CompositionSignature + LeftCancellativeCompositionSignature + RightCancellativeCompositionSignature
{
    fn try_difference(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem>;
}
impl<S: CancellativeCompositionSignature> LeftCancellativeCompositionSignature for S {
    fn try_left_difference(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        self.try_difference(a, b)
    }
}
impl<S: CancellativeCompositionSignature> RightCancellativeCompositionSignature for S {
    fn try_right_difference(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        self.try_difference(a, b)
    }
}

/// A set with a special element `e` called the identity element.
#[signature_meta_trait]
pub trait IdentitySignature: SetSignature {
    /// Returns the identity element `e`.
    fn identity(&self) -> Self::Elem;
}

/// When the solution to `compose(x, a)` = `e` for `x` given `a` is unique whenever it exists.
#[signature_meta_trait]
pub trait TryLeftInverseSignature: IdentitySignature + CompositionSignature {
    /// Return `x` such that `compose(x, a)` = `e` or `None` if no such `x` exists.
    fn try_left_inverse(&self, a: &Self::Elem) -> Option<Self::Elem>;
}

/// When the solution to `compose(a, x)` = `e` for `x` given `a` is unique whenever it exists.
#[signature_meta_trait]
pub trait TryRightInverseSignature: IdentitySignature + CompositionSignature {
    /// Return `x` such that `compose(a, x)` = `e` or `None` if no such `x` exists.
    fn try_right_inverse(&self, a: &Self::Elem) -> Option<Self::Elem>;
}

/// When the solution to `compose(x, a)` = `compose(a, x)` = `e` for `x` given `a` is unique whenever it exists.
#[signature_meta_trait]
pub trait TryInverseSignature: IdentitySignature + CompositionSignature {
    /// Return `x` such that `compose(x, a)` = `compose(a, x)` = `e` or `None` if no such `x` exists.
    ///
    /// Note, whenever `try_inverse` returns `Some`, `try_left_inverse` and `try_right_inverse` must also return the same value.
    /// Also, whenever `try_left_inverse` and `try_right_inverse` both return a value, it must be the same value an `try_inverse` must also return that same value.
    fn try_inverse(&self, a: &Self::Elem) -> Option<Self::Elem>;
}

impl<S: TryInverseSignature + CommutativeCompositionSignature> TryLeftInverseSignature for S {
    fn try_left_inverse(&self, a: &Self::Elem) -> Option<Self::Elem> {
        self.try_inverse(a)
    }
}

impl<S: TryInverseSignature + CommutativeCompositionSignature> TryRightInverseSignature for S {
    fn try_right_inverse(&self, a: &Self::Elem) -> Option<Self::Elem> {
        self.try_inverse(a)
    }
}

/// When `compose(x, e)` = `compose(e, x)` = `x` for all `x`.
#[signature_meta_trait]
pub trait MonoidSignature: IdentitySignature + AssociativeCompositionSignature {
    fn compose_list(&self, elems: Vec<impl Borrow<Self::Elem>>) -> Self::Elem {
        if elems.is_empty() {
            self.identity()
        } else {
            self.compose_nonempty_list(elems).unwrap()
        }
    }

    fn nat_pow(&self, a: &Self::Elem, n: &Natural) -> Self::Elem {
        if *n == Natural::ZERO {
            self.identity()
        } else if *n == Natural::ONE {
            a.clone()
        } else {
            debug_assert!(*n >= Natural::TWO);
            let bits: Vec<_> = n.bits().collect();
            let mut pows = vec![a.clone()];
            while pows.len() < bits.len() {
                pows.push(self.compose(pows.last().unwrap(), pows.last().unwrap()));
            }
            let count = bits.len();
            debug_assert_eq!(count, pows.len());
            let mut ans = self.identity();
            for i in 0..count {
                if bits[i] {
                    ans = self.compose(&ans, &pows[i]);
                }
            }
            ans
        }
    }
}

/// When inverses always exist.
#[signature_meta_trait]
pub trait GroupSignature:
    MonoidSignature
    + TryInverseSignature
    + TryLeftInverseSignature
    + TryRightInverseSignature
    + LeftCancellativeCompositionSignature
    + RightCancellativeCompositionSignature
{
    fn inverse(&self, a: &Self::Elem) -> Self::Elem;

    fn int_pow(&self, a: &Self::Elem, n: &Integer) -> Self::Elem {
        #[allow(clippy::comparison_chain)]
        if *n == Integer::ZERO {
            self.identity()
        } else if *n > Integer::ZERO {
            self.nat_pow(a, &n.abs())
        } else {
            self.nat_pow(&self.inverse(a), &n.abs())
        }
    }

    fn generated_finite_subgroup(&self, gens: Vec<Self::Elem>) -> FiniteSubgroup<Self::Elem>
    where
        Self::Elem: std::hash::Hash + Eq,
    {
        //generate subgroup by adding all generated elements
        let mut sg = HashSet::new();
        sg.insert(self.identity());

        let mut boundary = vec![self.identity()];
        let mut next_boundary = vec![];
        let mut y;
        while !boundary.is_empty() {
            println!("{}", sg.len());
            for x in &boundary {
                for g in &gens {
                    y = self.compose(x, g);
                    if !sg.contains(&y) {
                        sg.insert(y.clone());
                        next_boundary.push(y);
                    }
                }
            }
            boundary = next_boundary.clone();
            next_boundary = vec![];
        }

        FiniteSubgroup {
            elems: sg.into_iter().collect(),
        }
    }
}

#[signature_meta_trait]
pub trait AbelianGroupSignature:
    GroupSignature + CommutativeCompositionSignature + CancellativeCompositionSignature
{
}
impl<S: GroupSignature + CommutativeCompositionSignature + CancellativeCompositionSignature>
    AbelianGroupSignature for S
{
}

#[derive(Debug, Clone)]
pub struct FiniteSubgroup<Set> {
    elems: Vec<Set>,
}

impl<Set> FiniteSubgroup<Set> {
    pub fn size(&self) -> usize {
        self.elems.len()
    }

    pub fn elements(&self) -> impl Iterator<Item = &Set> {
        self.elems.iter()
    }
}

/// A left group action on a set
pub trait LeftGroupActionSignature<Group: GroupSignature, Set: SetSignature>: Signature {
    fn group(&self) -> &Group;
    fn set(&self) -> &Set;
    fn left_apply(&self, g: &Group::Elem, x: &Set::Elem) -> Set::Elem;
}

/// A right group action on a set
pub trait RightGroupActionSignature<Set: SetSignature, Group: GroupSignature>: Signature {
    fn group(&self) -> &Group;
    fn set(&self) -> &Set;
    fn right_apply(&self, g: &Group::Elem, x: &Set::Elem) -> Set::Elem;
}
