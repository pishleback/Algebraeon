use algebraeon_macros::{signature_meta_trait, skip_meta};
use algebraeon_nzq::traits::Abs;
use algebraeon_nzq::{Integer, Natural};
use algebraeon_sets::structure::{MetaType, SetSignature};
use std::borrow::Borrow;
use std::collections::{HashMap, HashSet};

/// A set with a binary operation of composition.
#[signature_meta_trait]
pub trait CompositionSignature: SetSignature {
    fn compose(&self, a: &Self::Set, b: &Self::Set) -> Self::Set;
    fn compose_mut(&self, a: &mut Self::Set, b: &Self::Set) {
        *a = self.compose(a, b);
    }
}

/// When composition is associative.
#[signature_meta_trait]
pub trait AssociativeCompositionSignature: CompositionSignature {
    /// Returns `None` if the list is empty.
    fn compose_nonempty_list(&self, mut elems: Vec<impl Borrow<Self::Set>>) -> Option<Self::Set> {
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
    /// Try to find `x` such that `a` = `compose(b, x).`
    fn try_left_difference(&self, a: &Self::Set, b: &Self::Set) -> Option<Self::Set>;
}

/// When `compose(x, a)` = `compose(y, a)` implies `x` = `y` for all `a`, `x`, `y`.
#[signature_meta_trait]
pub trait RightCancellativeCompositionSignature: CompositionSignature {
    /// Try to find `x` such that `a` = `compose(x, b).`
    fn try_right_difference(&self, a: &Self::Set, b: &Self::Set) -> Option<Self::Set>;
}

#[signature_meta_trait]
pub trait CancellativeCompositionSignature: CompositionSignature {}
impl<S: LeftCancellativeCompositionSignature + RightCancellativeCompositionSignature>
    CancellativeCompositionSignature for S
{
}

/// A set with a special element `e` called the identity element.
#[signature_meta_trait]
pub trait IdentitySignature: SetSignature {
    /// Returns the identity element `e`.
    fn identity(&self) -> Self::Set;
}

/// When `compose(x, e)` = `compose(e, x)` = `x` for all `x`.
#[signature_meta_trait]
pub trait MonoidSignature: IdentitySignature + AssociativeCompositionSignature {
    fn compose_list(&self, elems: Vec<impl Borrow<Self::Set>>) -> Self::Set {
        if elems.is_empty() {
            self.identity()
        } else {
            self.compose_nonempty_list(elems).unwrap()
        }
    }

    fn nat_pow(&self, a: &Self::Set, n: &Natural) -> Self::Set {
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

/// When the solution to `compose(x, a)` = `e` for `x` given `a` is unique whenever it exists.
#[signature_meta_trait]
pub trait TryLeftInverseSignature: MonoidSignature {
    /// Return `x` such that `compose(x, a)` = `e` or `None` if no such `x` exists.
    fn try_left_inverse(&self, a: &Self::Set) -> Option<Self::Set>;
}

/// When the solution to `compose(a, x)` = `e` for `x` given `a` is unique whenever it exists.
#[signature_meta_trait]
pub trait TryRightInverseSignature: MonoidSignature {
    /// Return `x` such that `compose(a, x)` = `e` or `None` if no such `x` exists.
    fn try_right_inverse(&self, a: &Self::Set) -> Option<Self::Set>;
}

/// When the solution to `compose(x, a)` = `compose(a, x)` = `e` for `x` given `a` is unique whenever it exists.
#[signature_meta_trait]
pub trait TryInverseSignature: MonoidSignature {
    /// Return `x` such that `compose(x, a)` = `compose(a, x)` = `e` or `None` if no such `x` exists.
    ///
    /// Note, whenever `try_inverse` returns `Some`, `try_left_inverse` and `try_right_inverse` must also return the same value.
    /// Also, whenever `try_left_inverse` and `try_right_inverse` both return a value, it must be the same value an `try_inverse` must also return that same value.
    fn try_inverse(&self, a: &Self::Set) -> Option<Self::Set>;
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
    fn inverse(&self, a: &Self::Set) -> Self::Set;

    fn int_pow(&self, a: &Self::Set, n: &Integer) -> Self::Set {
        #[allow(clippy::comparison_chain)]
        if *n == Integer::ZERO {
            self.identity()
        } else if *n > Integer::ZERO {
            self.nat_pow(a, &n.abs())
        } else {
            self.nat_pow(&self.inverse(a), &n.abs())
        }
    }

    fn generated_finite_subgroup_table(
        &self,
        generators: Vec<Self::Set>,
    ) -> (
        crate::composition_table::group::FiniteGroupMultiplicationTable,
        Vec<Self::Set>,
        HashMap<Self::Set, usize>,
    )
    where
        Self::Set: std::hash::Hash + Eq,
    {
        let mut n = 0;
        let mut idx_to_elem: Vec<Self::Set> = vec![];
        let mut elem_to_idx: HashMap<Self::Set, usize> = HashMap::new();
        let mut mul: Vec<Vec<Option<usize>>> = vec![];
        let mut to_mul: Vec<(usize, usize)> = vec![];

        macro_rules! add_elem {
            ($elem : expr) => {{
                debug_assert_eq!(idx_to_elem.len(), n);
                debug_assert_eq!(elem_to_idx.len(), n);
                debug_assert_eq!(mul.len(), n);
                for m in &mul {
                    debug_assert_eq!(m.len(), n);
                }
                if !elem_to_idx.contains_key(&$elem) {
                    n += 1;
                    let k = elem_to_idx.len();
                    idx_to_elem.push($elem.clone());
                    elem_to_idx.insert($elem, k);
                    for i in (0..k) {
                        mul[i].push(None);
                        to_mul.push((i, k));
                        to_mul.push((k, i));
                    }
                    mul.push(vec![None; k + 1]);
                    to_mul.push((k, k));
                    k
                } else {
                    *elem_to_idx.get(&$elem).unwrap()
                }
            }};
        }

        add_elem!(self.identity());
        for g in generators {
            add_elem!(g);
        }
        #[allow(clippy::manual_while_let_some)]
        while !to_mul.is_empty() {
            let (i, j) = to_mul.pop().unwrap();
            let k = add_elem!(self.compose(&idx_to_elem[i], &idx_to_elem[j]));
            debug_assert!(mul[i][j].is_none());
            mul[i][j] = Some(k);
        }
        drop(to_mul);
        let mul = mul
            .into_iter()
            .map(|m| m.into_iter().map(|x| x.unwrap()).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        let inv = idx_to_elem
            .iter()
            .map(|elem| *elem_to_idx.get(&self.inverse(elem)).unwrap())
            .collect::<Vec<_>>();

        let grp = crate::composition_table::group::FiniteGroupMultiplicationTable::new_unchecked(
            n, 0, inv, mul, None, None,
        );

        #[cfg(debug_assertions)]
        grp.check_state().unwrap();

        (grp, idx_to_elem, elem_to_idx)
    }

    fn generated_finite_subgroup(&self, gens: Vec<Self::Set>) -> FiniteSubgroup<Self::Set>
    where
        Self::Set: std::hash::Hash + Eq,
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
pub trait AbelianGroupSignature: GroupSignature + CommutativeCompositionSignature {}
impl<S: GroupSignature + CommutativeCompositionSignature> AbelianGroupSignature for S {}

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
