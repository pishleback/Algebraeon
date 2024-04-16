use std::collections::{HashMap, HashSet};
use std::hash::Hash;

use itertools::Itertools;

use super::{examples::c2::C2, group::Group};

#[derive(Debug, Clone)]
pub struct Cycle {
    cyc: Vec<usize>,
}

impl Cycle {
    pub fn new(cyc: Vec<usize>) -> Result<Self, &'static str> {
        let mut present = HashSet::new();
        for i in &cyc {
            if present.contains(i) {
                return Err("Duplicate element in cycle");
            }
            present.insert(i);
        }
        Ok(Self { cyc })
    }

    fn new_unchecked(cyc: Vec<usize>) -> Self {
        Self { cyc }
    }

    pub fn len(&self) -> usize {
        self.cyc.len()
    }
}

impl std::convert::From<Cycle> for Permutation {
    fn from(cyc: Cycle) -> Self {
        let n = *cyc.cyc.iter().max().unwrap_or(&0);
        let mut perm = vec![0; n];
        for i in 0..n {
            perm[i] = i;
        }
        for i in 0..cyc.cyc.len() {
            perm[cyc.cyc[i]] = cyc.cyc[(i + 1) % n];
        }
        Self::new_unchecked(perm)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Permutation {
    perm: Vec<usize>,
}

impl<const N: usize> From<super::examples::symmetric::Permutation<N>> for Permutation {
    fn from(value: super::examples::symmetric::Permutation<N>) -> Self {
        Self::new_unchecked((0..N).map(|i| value.call(i).unwrap()).collect())
    }
}

fn reduced(mut perm: Vec<usize>) -> Vec<usize> {
    while !perm.is_empty() {
        if *perm.last().unwrap() + 1 == perm.len() {
            perm.pop().unwrap();
        } else {
            break;
        }
    }
    perm
}

impl Permutation {
    pub fn new(perm: Vec<usize>) -> Result<Self, &'static str> {
        let n = perm.len();
        //check that the numbers in forward are 0, 1, ..., n-1 in some order
        let mut present = (0..n).map(|_| false).collect_vec();
        for i in &perm {
            if !(*i < n) {
                return Err("Permutation value out of range");
            }
            present[*i] = true;
        }
        for is_present in present {
            if !is_present {
                return Err("Not a valid permutation");
            }
        }

        Ok(Self::new_unchecked(perm))
    }

    pub fn new_unchecked(perm: Vec<usize>) -> Self {
        Self {
            perm: reduced(perm),
        }
    }

    pub fn n(&self) -> usize {
        self.perm.len()
    }

    pub fn call(&self, x: usize) -> usize {
        match self.perm.get(x) {
            Some(y) => *y,
            None => x,
        }
    }

    pub fn sign(&self) -> C2 {
        let mut s = C2::identity();
        for i in 0..self.perm.len() {
            for j in 0..i {
                if (i < j) != (self.call(i) < self.call(j)) {
                    s.compose_mut(&C2::Flip);
                }
            }
        }
        s
    }

    pub fn disjoint_cycles(&self) -> Vec<Cycle> {
        let n = self.perm.len();
        let mut missing: std::collections::HashSet<usize> = (0..n).collect();
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
            if cycle.len() >= 2 {
                cycles.push(Cycle { cyc: cycle });
            }
        }
        cycles
    }

    pub fn cycle_shape(&self) -> Vec<usize> {
        let mut shape = self.disjoint_cycles().iter().map(|c| c.len()).collect_vec();
        shape.sort();
        shape
    }

    pub fn all_permutations(n: usize) -> impl Iterator<Item = Self> {
        (0..n).permutations(n).map(|perm| Self::new_unchecked(perm))
    }

    pub fn symmetric_composition_table(
        n: usize,
    ) -> (
        crate::groups::composition_table::group::Group,
        Vec<Self>,
        HashMap<Self, usize>,
    ) {
        Self::generated_finite_subgroup_table(
            Self::all_permutations(n)
                .filter(|p| p.cycle_shape() == vec![2])
                .collect(),
        )
    }

    pub fn alternating_composition_table(
        n: usize,
    ) -> (
        crate::groups::composition_table::group::Group,
        Vec<Self>,
        HashMap<Self, usize>,
    ) {
        Self::generated_finite_subgroup_table(
            Self::all_permutations(n)
                .filter(|p| p.sign() == C2::Ident)
                .collect(),
        )
    }
}

// impl PartialEq for Permutation {
//     fn eq(&self, other: &Self) -> bool {
//         let n = std::cmp::max(self.perm.len(), other.perm.len());
//         (0..n).all(|i| self.call(i) == other.call(i))
//     }
// }

// impl Eq for Permutation {}

// impl Hash for Permutation {
//     fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
//         let mut self_reduced = self.clone();
//         self_reduced.reduce();
//         self_reduced.perm.hash(state);
//     }
// }

impl Group for Permutation {
    fn identity() -> Self {
        Self { perm: vec![] }
    }

    fn inverse(self) -> Self {
        let mut inv_perm = vec![0; self.perm.len()];
        for (i, j) in self.perm.into_iter().enumerate() {
            inv_perm[j] = i;
        }
        Self::new_unchecked(inv_perm)
    }

    fn compose_refs(a: &Self, b: &Self) -> Self {
        let n = std::cmp::max(a.perm.len(), b.perm.len());
        Self::new_unchecked((0..n).map(|i| a.call(b.call(i))).collect())
    }

    fn compose_mut(&mut self, other: &Self) {
        *self = Self::compose_refs(self, other);
    }
}

impl std::fmt::Display for Permutation {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_composition() {
        let a = Permutation::new(vec![0, 2, 1, 3]).unwrap();
        let b = Permutation::new(vec![1, 2, 0, 3]).unwrap();
        let c = Permutation::new(vec![2, 1, 0]).unwrap();

        println!("a = {}", a);
        println!("b = {}", b);
        println!("ab = {}", Permutation::compose_refs(&a, &b));
        println!("c = {}", c);

        assert_eq!(Permutation::compose(a, b), c);
    }

    // #[test]
    // fn test_dcn() {
    //     let a = Permutation::new(vec![0, 2, 1, 5, 3, 4]).unwrap();
    //     assert_eq!(
    //         a.disjoint_cycles(),
    //         vec![Cycle::new(vec![1, 2]), Cycle::new(vec![3, 4, 5])]
    //     );
    // }

    #[test]
    fn test_sign() {
        let a = Permutation::new(vec![0, 2, 1]).unwrap();
        let b = Permutation::new(vec![1, 2, 0]).unwrap();
        let c = Permutation::new(vec![2, 1, 0]).unwrap();

        println!("a = {}", a);
        println!("b = {}", b);
        println!("c = {}", c);

        assert_eq!(a.sign(), C2::Flip);
        assert_eq!(b.sign(), C2::Ident);
        assert_eq!(c.sign(), C2::Flip);
    }

    #[test]
    fn test_all_permutations() {
        assert_eq!(Permutation::all_permutations(0).collect_vec().len(), 1);
        assert_eq!(Permutation::all_permutations(1).collect_vec().len(), 1);
        assert_eq!(Permutation::all_permutations(2).collect_vec().len(), 2);
        assert_eq!(Permutation::all_permutations(3).collect_vec().len(), 6);
        assert_eq!(Permutation::all_permutations(4).collect_vec().len(), 24);
        assert_eq!(Permutation::all_permutations(5).collect_vec().len(), 120);
    }
}
