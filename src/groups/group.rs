use std::fmt::Debug;

pub trait Group: Clone + PartialEq + Eq {
    fn identity() -> Self;

    fn inverse(self) -> Self;
    fn inverse_ref(&self) -> Self {
        self.clone().inverse()
    }

    fn compose(a: Self, b: Self) -> Self;
    fn compose_lref(a: &Self, b: Self) -> Self {
        Self::compose(a.clone(), b)
    }
    fn compose_rref(a: Self, b: &Self) -> Self {
        Self::compose(a, b.clone())
    }
    fn compose_refs(a: &Self, b: &Self) -> Self {
        Self::compose(a.clone(), b.clone())
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
