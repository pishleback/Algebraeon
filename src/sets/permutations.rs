use super::super::finite_group_tables::group::Group;
use std::collections::HashMap;

#[derive(PartialEq, Eq, Hash, Debug, Clone)]
pub struct Permutation(Vec<usize>);

impl Permutation {
    pub fn new(p: Vec<usize>) -> Result<Permutation, ()> {
        let n = p.len();
        let mut present = vec![false; n];
        for i in &p {
            if !(*i < n) {
                return Err(());
            }
            present[*i] = true;
        }
        for is_present in present {
            if !is_present {
                return Err(());
            }
        }
        Ok(Permutation(p))
    }

    pub fn new_ident(n: usize) -> Permutation {
        Permutation((0..n).collect())
    }

    pub fn compose(&self, other: &Permutation) -> Result<Permutation, &'static str> {
        let n = self.0.len();
        if !(n == other.0.len()) {
            return Err("size mismatch");
        };
        return Ok(Self(other.0.iter().map(|i| self.0[*i]).collect()));
    }

    pub fn invert(&self) -> Self {
        Self(
            (0..self.0.len())
                .map(|i| {
                    for (a, b) in self.0.iter().enumerate() {
                        if i == *b {
                            return a;
                        }
                    }
                    panic!();
                })
                .collect(),
        )
    }

    pub fn call(&self, x: usize) -> Result<usize, &'static str> {
        if !(x < self.0.len()) {
            return Err("argument too large");
        }
        Ok(self.0[x])
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    //true = even sign
    //false = odd sign
    pub fn sign(&self) -> bool {
        let mut s = true;
        for i in 0..self.0.len() {
            for j in 0..i {
                if (i < j) != (self.call(i) < self.call(j)) {
                    s = !s;
                }
            }
        }
        s
    }

    pub fn disjoint_cycle_notation(&self) -> String {
        let mut missing = std::collections::HashSet::new();
        for i in 0..self.len() {
            missing.insert(i);
        }
        let mut cycles = vec![];
        while missing.len() > 0 {
            let mut cycle = vec![];
            let x = *missing.iter().min().unwrap();
            let mut i = x;
            loop {
                cycle.push(i);
                missing.remove(&i);
                i = self.call(i).unwrap();
                if i == x {
                    break;
                }
            }
            cycles.push(cycle);
        }

        cycles.retain(|cycle| cycle.len() != 1);

        if cycles.len() == 0 {
            return "ident".to_string();
        }

        return cycles
            .iter()
            .map(|cycle| {
                "(".to_owned()
                    + &cycle
                        .iter()
                        .map(|x| x.to_string())
                        .collect::<Vec<String>>()
                        .join(" ")
                    + ")"
            })
            .collect::<Vec<String>>()
            .join(" ");
    }
}

pub fn all_perms(k: usize) -> Vec<Permutation> {
    if k == 0 {
        return vec![Permutation(vec![])];
    } else {
        let prev_perms = all_perms(k - 1);
        let mut next_perms: Vec<Permutation> = vec![];
        for i in 0..k {
            for pp in prev_perms.iter() {
                let mut next_perm = vec![];
                for j in 0..i {
                    next_perm.push(pp.call(j).unwrap() + 1)
                }
                next_perm.push(0);
                for j in i..k - 1 {
                    next_perm.push(pp.call(j).unwrap() + 1)
                }
                next_perms.push(Permutation(next_perm));
            }
        }
        return next_perms;
    }
}

pub fn alternating_group_structure(
    n: usize,
) -> (Group, Vec<Permutation>, HashMap<Permutation, usize>) {
    use super::permutations::*;
    let alt_perms: Vec<Permutation> = all_perms(n).into_iter().filter(|p| p.sign()).collect();
    (
        Group::from_raw_model_unchecked(
            alt_perms.clone(),
            || Permutation::new_ident(n),
            |x: Permutation| x.invert(),
            |x: Permutation, y: Permutation| x.compose(&y).unwrap(),
            Some(n <= 3),           //A0, A1, A2, A3 are abelian, the rest are not
            Some(n == 3 || 5 <= n), //A3, A5, A6, A7, ... are simple, the rest are not
        ),
        alt_perms.clone(),
        alt_perms
            .clone()
            .into_iter()
            .enumerate()
            .map(|(idx, perm)| (perm, idx))
            .collect(),
    )
}

pub fn symmetric_group_structure(
    n: usize,
) -> (Group, Vec<Permutation>, HashMap<Permutation, usize>) {
    let perms: Vec<Permutation> = all_perms(n);
    (
        Group::from_raw_model_unchecked(
            perms.clone(),
            || Permutation::new_ident(n),
            |x: Permutation| x.invert(),
            |x: Permutation, y: Permutation| x.compose(&y).unwrap(),
            Some(n <= 2), //S0, S1, S2 are abelian, the rest are not
            Some(n == 2), //S2 is simple, the rest are not
        ),
        perms.clone(),
        perms
            .clone()
            .into_iter()
            .enumerate()
            .map(|(idx, perm)| (perm, idx))
            .collect(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn alternating() {
        let mut size = 1;
        for k in 0..6 {
            let (mut grp, perms, _elems) = alternating_group_structure(k);
            grp.cache_conjugacy_classes();
            assert_eq!(grp.elems().len(), perms.len());
            //|A0|=1 |A1|=1 |A2|=1 |A3|=3 |A4|=12 |A5|=30 ...
            if k <= 2 {
                assert_eq!(grp.elems().len(), 1);
            } else {
                assert_eq!(grp.elems().len(), size / 2);
            }
            grp.check_state().unwrap();
            size *= k + 1;
        }
    }

    #[test]
    fn symmetric() {
        let mut size = 1;
        for k in 0..6 {
            let (mut grp, perms, _elems) = symmetric_group_structure(k);
            grp.cache_conjugacy_classes();
            assert_eq!(grp.elems().len(), perms.len());
            assert_eq!(grp.elems().len(), size);
            grp.check_state().unwrap();
            size *= k + 1;
        }
    }
}
