use std::{borrow::Borrow, rc::Rc};

use super::function::*;

pub trait SetT: Borrow<Set> + Clone {}
impl<T: Borrow<Set> + Clone> SetT for T {}

#[derive(Debug)]
pub struct Set {
    n: usize,
}

impl Set {
    pub fn new(n: usize) -> Self {
        Set { n }
    }

    pub fn size(&self) -> usize {
        self.n
    }

    pub fn elems(&self) -> std::ops::Range<usize> {
        0..self.n
    }

    pub fn contains(&self, x: usize) -> bool {
        x < self.n
    }
}

struct DisjointUnion<FirstT: SetT, SecondT: SetT> {
    first: FirstT,
    second: SecondT,
}

impl<FirstT: SetT, SecondT: SetT> DisjointUnion<FirstT, SecondT> {
    pub fn to_set(
        &self,
    ) -> (
        Rc<Set>,
        Function<FirstT, Rc<Set>>,
        Function<SecondT, Rc<Set>>,
    ) {
        let a = self.first.borrow().size();
        let b = self.second.borrow().size();
        let n = a + b;
        let set = Rc::new(Set { n });
        (
            set.clone(),
            Function::new_unchecked(
                self.first.clone(),
                set.clone(),
                (0..a).map(|x| 0 + x).collect(),
            ),
            Function::new_unchecked(
                self.second.clone(),
                set.clone(),
                (0..b).map(|x| a + x).collect(),
            ),
        )
    }

    pub fn size(&self) -> usize {
        self.first.borrow().size() + self.second.borrow().size()
    }
}

struct CartesianProduct<FirstT: SetT, SecondT: SetT> {
    first: FirstT,
    second: SecondT,
}

impl<FirstT: SetT, SecondT: SetT> CartesianProduct<FirstT, SecondT> {
    pub fn to_set(
        &self,
    ) -> (
        Rc<Set>,
        Function<Rc<Set>, FirstT>,
        Function<Rc<Set>, SecondT>,
    ) {
        let a = self.first.borrow().size();
        let b = self.second.borrow().size();
        let n = a * b;
        let set = Rc::new(Set { n });
        (
            set.clone(),
            Function::new_unchecked(
                set.clone(),
                self.first.clone(),
                (0..n).map(|x| x % a).collect(),
            ),
            Function::new_unchecked(
                set.clone(),
                self.second.clone(),
                (0..n).map(|x| x / a).collect(),
            ),
        )
    }

    pub fn size(&self) -> usize {
        self.first.borrow().size() * self.second.borrow().size()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    pub fn test_disjoint_union() {
        let a_set = Set { n: 3 };
        let b_set = Set { n: 5 };
        let union = DisjointUnion {
            first: &a_set,
            second: &b_set,
        };
        assert_eq!(union.size(), a_set.size() + b_set.size());
        let (union_set, inc_a, inc_b) = union.to_set();
        assert_eq!(union_set.size(), union.size());
        assert!(Rc::ptr_eq(&union_set, &inc_a.range()));
        assert!(Rc::ptr_eq(&union_set, &inc_b.range()));
        inc_a.check_state().unwrap();
        inc_b.check_state().unwrap();
        assert!(inc_a.is_injective());
        assert!(inc_b.is_injective());
    }

    #[test]
    pub fn test_cartesian_product() {
        let a_set = Set { n: 3 };
        let b_set = Set { n: 5 };
        let prod = CartesianProduct {
            first: &a_set,
            second: &b_set,
        };
        assert_eq!(prod.size(), a_set.size() * b_set.size());
        let (prod_set, proj_a, proj_b) = prod.to_set();
        assert_eq!(prod_set.size(), prod.size());
        assert!(Rc::ptr_eq(&prod_set, &proj_a.domain()));
        assert!(Rc::ptr_eq(&prod_set, &proj_b.domain()));
        proj_a.check_state().unwrap();
        proj_b.check_state().unwrap();
        assert!(proj_a.is_surjective());
        assert!(proj_b.is_surjective());
    }
}
