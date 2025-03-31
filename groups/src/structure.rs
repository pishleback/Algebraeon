use algebraeon_nzq::traits::Abs;
use algebraeon_nzq::*;
use itertools::Itertools;
use std::{
    borrow::Borrow,
    collections::{HashMap, HashSet},
    fmt::Debug,
    hash::Hash,
};

pub trait Group: Debug + Clone + PartialEq + Eq {
    fn identity() -> Self;

    fn inverse(self) -> Self;
    fn inverse_ref(&self) -> Self {
        self.clone().inverse()
    }

    fn compose_mut(&mut self, other: &Self);
    fn compose(mut a: Self, b: Self) -> Self {
        Self::compose_mut(&mut a, &b);
        a
    }
    fn compose_lref(a: &Self, b: Self) -> Self {
        Self::compose(a.clone(), b)
    }
    fn compose_rref(a: Self, b: &Self) -> Self {
        Self::compose(a, b.clone())
    }
    fn compose_refs(a: &Self, b: &Self) -> Self {
        Self::compose(a.clone(), b.clone())
    }

    fn compose_list(elems: Vec<impl Borrow<Self>>) -> Self {
        let mut ans = Self::identity();
        for elem in elems {
            ans.compose_mut(elem.borrow());
        }
        ans
    }

    fn nat_pow(&self, n: &Natural) -> Self {
        if *n == Natural::ZERO {
            Self::identity()
        } else if *n == Natural::ONE {
            self.clone()
        } else {
            debug_assert!(*n >= Natural::TWO);
            let bits: Vec<_> = n.bits().collect();
            let mut pows = vec![self.clone()];
            while pows.len() < bits.len() {
                pows.push(Self::compose_refs(
                    &pows.last().unwrap(),
                    &pows.last().unwrap(),
                ));
            }
            let count = bits.len();
            debug_assert_eq!(count, pows.len());
            let mut ans = Self::identity();
            for i in 0..count {
                if bits[i] {
                    ans.compose_mut(&pows[i]);
                }
            }
            ans
        }
    }

    fn int_pow(&self, n: &Integer) -> Self {
        if *n == Integer::ZERO {
            Self::identity()
        } else if *n > Integer::ZERO {
            self.nat_pow(&n.abs())
        } else {
            self.nat_pow(&n.abs()).inverse()
        }
    }

    fn generated_finite_subgroup_table(
        generators: Vec<Self>,
    ) -> (
        crate::composition_table::group::FiniteGroupMultiplicationTable,
        Vec<Self>,
        HashMap<Self, usize>,
    )
    where
        Self: std::hash::Hash,
    {
        let mut n = 0;
        let mut idx_to_elem: Vec<Self> = vec![];
        let mut elem_to_idx: HashMap<Self, usize> = HashMap::new();
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

        add_elem!(Self::identity());
        for g in generators {
            add_elem!(g);
        }
        while !to_mul.is_empty() {
            let (i, j) = to_mul.pop().unwrap().clone();
            let k = add_elem!(Self::compose_refs(&idx_to_elem[i], &idx_to_elem[j]));
            debug_assert!(mul[i][j].is_none());
            mul[i][j] = Some(k);
        }
        drop(to_mul);
        let mul = mul
            .into_iter()
            .map(|m| m.into_iter().map(|x| x.unwrap()).collect_vec())
            .collect_vec();
        let inv = idx_to_elem
            .iter()
            .map(|elem| *elem_to_idx.get(&Self::inverse_ref(elem)).unwrap())
            .collect_vec();

        let grp =
            crate::composition_table::group::FiniteGroupMultiplicationTable::new_unchecked(n, 0, inv, mul, None, None);

        #[cfg(debug_assertions)]
        grp.check_state().unwrap();

        (grp, idx_to_elem, elem_to_idx)
    }

    fn generated_finite_subgroup(gens: Vec<Self>) -> FiniteSubgroup<Self>
    where
        Self: Hash,
    {
        //generate subgroup by adding all generated elements
        let mut sg = HashSet::new();
        sg.insert(Self::identity());

        let mut boundary = vec![Self::identity()];
        let mut next_boundary = vec![];
        let mut y;
        while boundary.len() > 0 {
            println!("{}", sg.len());
            for x in &boundary {
                for g in &gens {
                    y = Self::compose_refs(x, g);
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

pub trait Pow<ExpT> {
    fn pow(&self, exp: ExpT) -> Self;
}
impl<G: Group> Pow<Natural> for G {
    fn pow(&self, exp: Natural) -> Self {
        self.nat_pow(&exp)
    }
}
impl<G: Group> Pow<u8> for G {
    fn pow(&self, exp: u8) -> Self {
        self.nat_pow(&Natural::from(exp))
    }
}
impl<G: Group> Pow<u16> for G {
    fn pow(&self, exp: u16) -> Self {
        self.nat_pow(&Natural::from(exp))
    }
}
impl<G: Group> Pow<u32> for G {
    fn pow(&self, exp: u32) -> Self {
        self.nat_pow(&Natural::from(exp))
    }
}
impl<G: Group> Pow<u64> for G {
    fn pow(&self, exp: u64) -> Self {
        self.nat_pow(&Natural::from(exp))
    }
}
impl<G: Group> Pow<Integer> for G {
    fn pow(&self, exp: Integer) -> Self {
        self.int_pow(&exp)
    }
}
impl<G: Group> Pow<i8> for G {
    fn pow(&self, exp: i8) -> Self {
        self.int_pow(&Integer::from(exp))
    }
}
impl<G: Group> Pow<i16> for G {
    fn pow(&self, exp: i16) -> Self {
        self.int_pow(&Integer::from(exp))
    }
}
impl<G: Group> Pow<i32> for G {
    fn pow(&self, exp: i32) -> Self {
        self.int_pow(&Integer::from(exp))
    }
}
impl<G: Group> Pow<i64> for G {
    fn pow(&self, exp: i64) -> Self {
        self.int_pow(&Integer::from(exp))
    }
}

#[derive(Debug, Clone)]
pub struct FiniteSubgroup<G: Group> {
    elems: Vec<G>,
}

impl<G: Group> FiniteSubgroup<G> {
    pub fn size(&self) -> usize {
        self.elems.len()
    }

    pub fn elements(&self) -> impl Iterator<Item = &G> {
        self.elems.iter()
    }
}
