use std::collections::HashMap;

use itertools::Itertools;
use malachite_base::num::{arithmetic::traits::DivRem, basic::traits::Zero};
use malachite_nz::integer::Integer;

use crate::group::Group;

use super::c2::C2;

#[derive(Debug, Clone)]
pub struct Cycle<const N: usize> {
    cyc: Vec<usize>,
}

impl<const N: usize> Cycle<N> {
    pub fn new(cyc: Vec<usize>) -> Result<Self, &'static str> {
        let mut present = [false; N];
        for i in &cyc {
            if present[*i] {
                return Err("Duplicate element in cycle");
            }
            present[*i] = true;
        }
        Ok(Self { cyc })
    }

    pub fn len(&self) -> usize {
        self.cyc.len()
    }
}

impl<const N: usize> std::convert::From<Cycle<N>> for Permutation<N> {
    fn from(cyc: Cycle<N>) -> Self {
        let mut perm = [0; N];
        for i in 0..N {
            perm[i] = i;
        }
        let n = cyc.cyc.len();
        for i in 0..n {
            perm[cyc.cyc[i]] = cyc.cyc[(i + 1) % n];
        }
        Self { perm }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Permutation<const N: usize> {
    perm: [usize; N],
}

impl<const N: usize> TryFrom<super::super::permutation::Permutation> for Permutation<N> {
    type Error = ();

    fn try_from(value: super::super::permutation::Permutation) -> Result<Self, Self::Error> {
        let n = value.n();
        if N < n {
            Err(())
        } else {
            let mut perm = [0; N];
            for i in 0..N {
                perm[i] = value.call(i);
            }
            Ok(Self { perm })
        }
    }
}

impl<const N: usize> Permutation<N> {
    pub fn new(perm: [usize; N]) -> Result<Self, &'static str> {
        //check that the numbers in forward are 0, 1, ..., n-1 in some order
        let mut present = [false; N];
        for i in &perm {
            if !(*i < N) {
                return Err("Permutation value out of range");
            }
            present[*i] = true;
        }
        for is_present in present {
            if !is_present {
                return Err("Not a valid permutation");
            }
        }

        Ok(Self { perm })
    }

    pub fn new_from_cycles(cycles: Vec<Vec<usize>>) -> Result<Self, &'static str> {
        let mut parsed_cycles = vec![];
        for c in cycles {
            parsed_cycles.push(Cycle::<N>::new(c)?);
        }
        Ok(Self::compose_list(
            parsed_cycles
                .into_iter()
                .map(|c| Into::<Permutation<N>>::into(c))
                .collect(),
        ))
    }

    pub fn call(&self, x: usize) -> Result<usize, &'static str> {
        if !(x < self.perm.len()) {
            return Err("argument too large");
        }
        Ok(self.perm[x])
    }

    pub fn sign(&self) -> C2 {
        let mut s = C2::identity();
        for i in 0..N {
            for j in 0..i {
                if (i < j) != (self.call(i) < self.call(j)) {
                    s.compose_mut(&C2::Flip);
                }
            }
        }
        s
    }

    pub fn disjoint_cycles(&self) -> Vec<Cycle<N>> {
        let mut missing: std::collections::HashSet<usize> = (0..N).collect();
        let mut cycles = vec![];
        while missing.len() > 0 {
            let mut cycle = vec![];
            let x = *missing.iter().min().unwrap();
            let mut i = x;
            loop {
                cycle.push(i);
                missing.remove(&i);
                i = self.perm[i];
                if i == x {
                    break;
                }
            }
            cycles.push(Cycle { cyc: cycle });
        }
        cycles
    }

    pub fn cycle_shape(&self) -> Vec<usize> {
        let mut shape = self.disjoint_cycles().iter().map(|c| c.len()).collect_vec();
        shape.sort();
        shape
    }

    pub fn all_permutations() -> impl Iterator<Item = Self> {
        (0..N).permutations(N).map(|perm| Self {
            perm: perm.try_into().unwrap(),
        })
    }

    pub fn symmetric_composition_table() -> (
        crate::composition_table::group::Group,
        Vec<Self>,
        HashMap<Self, usize>,
    ) {
        Self::generated_finite_subgroup_table(Self::all_permutations().collect())
    }

    pub fn alternating_composition_table() -> (
        crate::composition_table::group::Group,
        Vec<Self>,
        HashMap<Self, usize>,
    ) {
        Self::generated_finite_subgroup_table(
            Self::all_permutations()
                .filter(|p| p.sign() == C2::Ident)
                .collect(),
        )
    }
}

impl<const N: usize> std::fmt::Display for Permutation<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut cycles = self.disjoint_cycles();
        cycles.retain(|cycle| cycle.len() != 1);

        if cycles.len() == 0 {
            f.write_str("()");
        }

        let string = cycles
            .iter()
            .map(|cycle| {
                "(".to_owned()
                    + &cycle
                        .cyc
                        .iter()
                        .map(|x| x.to_string())
                        .collect::<Vec<String>>()
                        .join(" ")
                    + ")"
            })
            .collect::<Vec<String>>()
            .join(" ");

        f.write_str(string.as_str())
    }
}

impl<const N: usize> Group for Permutation<N> {
    fn identity() -> Self {
        let mut perm = [0; N];
        for i in 0..N {
            perm[i] = i;
        }
        Self { perm }
    }

    fn inverse(self) -> Self {
        let mut inv_perm = [0; N];
        for (i, j) in self.perm.into_iter().enumerate() {
            inv_perm[j] = i;
        }
        Self { perm: inv_perm }
    }

    fn compose_refs(a: &Self, b: &Self) -> Self {
        let mut comp_perm = [0; N];
        for i in 0..N {
            comp_perm[i] = a.perm[b.perm[i]];
        }
        Self { perm: comp_perm }
    }

    fn compose_mut(&mut self, other: &Self) {
        *self = Self::compose_refs(self, other);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_composition() {
        let a = Permutation::new([0, 2, 1]).unwrap();
        let b = Permutation::new([1, 2, 0]).unwrap();
        let c = Permutation::new([2, 1, 0]).unwrap();

        println!("a = {}", a);
        println!("b = {}", b);
        println!("ab = {}", Permutation::compose_refs(&a, &b));
        println!("c = {}", c);

        assert_eq!(Permutation::compose(a, b), c);
    }

    #[test]
    fn test_sign() {
        let a = Permutation::new([0, 2, 1]).unwrap();
        let b = Permutation::new([1, 2, 0]).unwrap();
        let c = Permutation::new([2, 1, 0]).unwrap();

        println!("a = {}", a);
        println!("b = {}", b);
        println!("c = {}", c);

        assert_eq!(a.sign(), C2::Flip);
        assert_eq!(b.sign(), C2::Ident);
        assert_eq!(c.sign(), C2::Flip);
    }

    #[test]
    fn test_all_permutations() {
        assert_eq!(Permutation::<0>::all_permutations().collect_vec().len(), 1);
        assert_eq!(Permutation::<1>::all_permutations().collect_vec().len(), 1);
        assert_eq!(Permutation::<2>::all_permutations().collect_vec().len(), 2);
        assert_eq!(Permutation::<3>::all_permutations().collect_vec().len(), 6);
        assert_eq!(Permutation::<4>::all_permutations().collect_vec().len(), 24);
        assert_eq!(
            Permutation::<5>::all_permutations().collect_vec().len(),
            120
        );
    }
}
