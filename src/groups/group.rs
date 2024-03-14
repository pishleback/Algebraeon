use std::{fmt::Debug, ops::Mul};

use malachite_base::num::{arithmetic::traits::UnsignedAbs, logic::traits::BitIterable};
use malachite_nz::{integer::Integer, natural::Natural};

pub trait Group: Clone + PartialEq + Eq {
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

    fn compose_list(elems: Vec<&Self>) -> Self {
        let mut ans = Self::identity();
        for elem in elems {
            ans.compose_mut(elem);
        }
        ans
    }

    fn nat_pow(&self, n: &Natural) -> Self {
        if *n == 0 {
            Self::identity()
        } else if *n == 1 {
            self.clone()
        } else {
            debug_assert!(*n >= 2);
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
        if *n == 0 {
            Self::identity()
        } else if *n > 0 {
            self.nat_pow(&n.unsigned_abs())
        } else {
            self.nat_pow(&(-n).unsigned_abs())
        }
    }

    fn finite_generated_subgroup(generators: Vec<Self>) -> Vec<Self> {
        todo!()
    }

    fn finite_generated_subgroup_table(
        generators: Vec<Self>,
    ) -> super::super::finite_group_tables::group::Group {
        todo!()
    }
}

// pub struct MultiplicativeNotation<G: Group> {
//     elem: G,
// }

// impl<G: Group> Mul<MultiplicativeNotation<G>> for MultiplicativeNotation<G> {
//     type Output = MultiplicativeNotation<G>;

//     fn mul(self, rhs: MultiplicativeNotation<G>) -> Self::Output {
//         MultiplicativeNotation {
//             elem: G::compose(self.elem, rhs.elem),
//         }
//     }
// }

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
