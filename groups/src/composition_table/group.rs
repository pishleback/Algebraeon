use super::normal_subgroup::NormalSubgroup;
use super::partition::GroupPartition;
use super::subgroup::Subgroup;
use super::subset::Subset;
use algebraeon_sets::combinatorics::Partition;
use rayon::prelude::*;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fmt::Debug;
use std::hash::Hash;

#[derive(Debug)]
pub struct FiniteGroupMultiplicationTable {
    n: usize,
    ident: usize,
    inv: Vec<usize>,
    mul: Vec<Vec<usize>>,
    conjugacy_classes: Option<Partition>,
    is_abelian: Option<bool>,
    is_simple: Option<bool>,
}

impl FiniteGroupMultiplicationTable {
    pub fn check_state(&self) -> Result<(), &'static str> {
        //check ident
        if self.ident >= self.n {
            return Err("bad ident elem");
        }
        //check inv
        if self.inv.len() != self.n {
            return Err("bad inv len");
        }
        for x in &self.inv {
            if *x >= self.n {
                return Err("bad inv elem");
            }
        }
        //check mul
        if self.mul.len() != self.n {
            return Err("bad mul left len");
        }
        for m in &self.mul {
            if m.len() != self.n {
                return Err("bad mul right len");
            }
            for x in m {
                if *x >= self.n {
                    return Err("bad mul elem");
                }
            }
        }
        //identity axiom
        for x in 0..self.n {
            if !(x == self.mul[x][self.ident] && x == self.mul[self.ident][x]) {
                return Err("identity axiom failed");
            }
        }
        //inv axiom
        for x in 0..self.n {
            if !(self.ident == self.mul[self.inv[x]][x] && self.ident == self.mul[x][self.inv[x]]) {
                return Err("inverse axiom failed");
            }
        }
        //assoc axiom
        #[allow(clippy::nonminimal_bool)]
        for x in 0..self.n {
            for y in 0..self.n {
                for z in 0..self.n {
                    if !(self.mul[x][self.mul[y][z]] == self.mul[self.mul[x][y]][z]) {
                        return Err("assoc axiom failed");
                    }
                }
            }
        }

        //cached is_abelian if present
        if let Some(claimed_is_abelian) = self.is_abelian {
            if claimed_is_abelian != self.compute_is_abelian() {
                return Err("incorrect is_abelian flag");
            }
        }

        //check is_simple
        if let Some(claimed_is_simple) = self.is_simple {
            if (self.normal_subgroups().len() == 2) != claimed_is_simple {
                return Err("is_simple flag is incorrect");
            }
        }

        Ok(())
    }

    pub fn new(
        n: usize,
        ident: usize,
        inv: Vec<usize>,
        mul: Vec<Vec<usize>>,
    ) -> Result<Self, &'static str> {
        let grp = FiniteGroupMultiplicationTable {
            n,
            ident,
            inv,
            mul,
            is_abelian: None,
            conjugacy_classes: None,
            is_simple: None,
        };
        match grp.check_state() {
            Ok(()) => Ok(grp),
            Err(msg) => Err(msg),
        }
    }

    pub fn new_unchecked(
        n: usize,
        ident: usize,
        inv: Vec<usize>,
        mul: Vec<Vec<usize>>,
        is_abelian: Option<bool>,
        is_simple: Option<bool>,
    ) -> Self {
        Self {
            n,
            ident,
            inv,
            mul,
            is_abelian,
            conjugacy_classes: None,
            is_simple,
        }
    }

    pub fn mul(&self, x: usize, y: usize) -> usize {
        self.mul[x][y]
    }

    pub fn inv(&self, x: usize) -> usize {
        self.inv[x]
    }

    pub fn ident(&self) -> usize {
        self.ident
    }

    pub fn from_raw_model<T: PartialEq + Eq + Hash + Clone + Debug>(
        elems: Vec<T>,
        ident: impl Fn() -> T,
        inv: impl Fn(T) -> T,
        mul: impl Fn(T, T) -> T,
    ) -> Result<Self, &'static str> {
        let grp = Self::from_raw_model_unchecked(elems, ident, inv, mul, None, None);
        match grp.check_state() {
            Ok(()) => Ok(grp),
            Err(msg) => Err(msg),
        }
    }

    pub fn from_raw_model_unchecked<T: PartialEq + Eq + Hash + Clone + Debug>(
        elems: Vec<T>,
        ident: impl Fn() -> T,
        inv: impl Fn(T) -> T,
        mul: impl Fn(T, T) -> T,
        is_abelian: Option<bool>,
        is_simple: Option<bool>,
    ) -> Self {
        let n = elems.len();
        let mut idx_map: HashMap<T, usize> = HashMap::new();
        for (i, x) in elems.iter().enumerate() {
            idx_map.insert(x.clone(), i);
        }

        Self::new_unchecked(
            n,
            idx_map[&ident()],
            {
                let mut inv_lookup = vec![];
                for x in &elems {
                    inv_lookup.push(idx_map[&inv(x.clone())]);
                }
                inv_lookup
            },
            {
                let mut mul_lookup = vec![];
                for (i, x) in elems.iter().enumerate() {
                    mul_lookup.push(vec![]);
                    for y in &elems {
                        mul_lookup[i].push(idx_map[&mul(x.clone(), y.clone())]);
                    }
                }
                mul_lookup
            },
            is_abelian,
            is_simple,
        )
    }

    pub fn size(&self) -> usize {
        self.n
    }

    pub fn elems(&self) -> std::ops::Range<usize> {
        0..self.n
    }

    #[allow(clippy::ptr_arg)]
    pub fn mul_many(&self, elems: &Vec<usize>) -> usize {
        if elems.is_empty() {
            return self.ident;
        }
        let mut prod = elems[0];
        for i in 1..elems.len() {
            prod = self.mul[prod][elems[i]];
        }
        prod
    }

    pub fn order(&self, x: usize) -> Result<usize, ()> {
        if x >= self.n {
            return Err(());
        }
        let mut y = x;
        let mut ord = 1;
        while y != self.ident {
            y = self.mul[x][y];
            ord += 1;
            debug_assert!(ord <= self.n);
        }
        Ok(ord)
    }

    fn compute_is_abelian(&self) -> bool {
        for x in 0..self.n {
            for y in 0..x {
                if self.mul[x][y] != self.mul[y][x] {
                    return false;
                }
            }
        }
        true
    }

    pub fn is_abelian(&self) -> bool {
        match &self.is_abelian {
            Some(flag) => *flag,
            None => self.compute_is_abelian(),
        }
    }

    fn compute_conjugacy_classes(&self) -> Partition {
        let mut unclassified_elems = self.elems().collect::<HashSet<_>>();
        let mut classes = vec![];
        let mut lookup = vec![0; self.n];
        while !unclassified_elems.is_empty() {
            let x = unclassified_elems.iter().next().unwrap();
            let mut class = HashSet::new();
            for g in self.elems() {
                class.insert(self.mul[self.mul[g][*x]][self.inv[g]]);
            }
            for y in &class {
                unclassified_elems.remove(y);
                lookup[*y] = classes.len();
            }
            classes.push(class);
        }
        Partition::new_unchecked(classes, lookup)
    }

    pub fn cache_conjugacy_classes(&mut self) {
        self.conjugacy_classes = Some(self.compute_conjugacy_classes());
    }

    pub fn conjugacy_class(&mut self, x: usize) -> Result<Subset, ()> {
        if x >= self.n {
            return Err(());
        }
        self.cache_conjugacy_classes();
        if let Some(conj_state) = &self.conjugacy_classes {
            Ok(Subset::new_unchecked(
                self,
                conj_state.class_containing(x).clone(),
            ))
        } else {
            panic!()
        }
    }

    pub fn conjugacy_classes(&self) -> GroupPartition {
        GroupPartition {
            group: self,
            partition: match &self.conjugacy_classes {
                Some(state) => state.clone(),
                None => self.compute_conjugacy_classes(),
            },
        }
    }

    //returns a hashmap of subgroups and a minimal generating set
    fn subgroups_impl(&self, only_normal: bool) -> Vec<(Subgroup, Subset)> {
        //choose one generator of each cyclic subgroup
        //The partition of elements where x~y iff <x>=<y> is such that an equivelence class is either a subset of a subgroup or doesnot intersect the subgroup
        let mut distinguished_gens = vec![];
        let mut subgroups: HashMap<Subgroup, Subset> = HashMap::new();
        for x in self.elems() {
            let singleton_x_subset = Subset::new_unchecked(self, HashSet::from_iter(vec![x]));
            let cyclic_sg = if only_normal {
                singleton_x_subset.normal_closure().unwrap().to_subgroup()
            } else {
                singleton_x_subset.generated_subgroup().unwrap()
            };
            #[allow(clippy::map_entry)]
            if !subgroups.contains_key(&cyclic_sg) {
                subgroups.insert(cyclic_sg, singleton_x_subset.clone());
                distinguished_gens.push(x);
            }
        }

        //generate all subgroups by itteratively adding generators
        let mut boundary = HashMap::new();
        for (sg, gens) in subgroups.clone() {
            boundary.insert(sg, gens);
        }
        let mut next_boundary = HashMap::new();
        while !boundary.is_empty() {
            #[allow(clippy::for_kv_map)]
            for (_sg, sg_gens) in &boundary {
                //compute and sort the next subgroups in parallel threads
                for (new_sg, new_gens) in distinguished_gens
                    .par_iter()
                    .map(|dgen| {
                        let mut new_gens = sg_gens.clone();
                        new_gens.add_elem(*dgen).unwrap();
                        let new_sg = if only_normal {
                            new_gens.normal_closure().unwrap().to_subgroup()
                        } else {
                            new_gens.generated_subgroup().unwrap()
                        };
                        (new_sg, new_gens)
                    })
                    .collect::<Vec<(Subgroup, Subset)>>()
                {
                    //collect the computed subgroups in the main thread
                    #[allow(clippy::map_entry)]
                    if !subgroups.contains_key(&new_sg) {
                        next_boundary.insert(new_sg.clone(), new_gens.clone());
                        subgroups.insert(new_sg, new_gens);
                    }
                }
            }
            boundary = next_boundary;
            next_boundary = HashMap::new();
        }
        #[allow(clippy::map_identity)]
        return subgroups
            .into_iter()
            .map(|(elems, gens)| (elems, gens))
            .collect();
    }

    pub fn subgroups(&self) -> Vec<(Subgroup, Subset)> {
        self.subgroups_impl(false)
    }

    pub fn normal_subgroups(&self) -> Vec<(NormalSubgroup, Subset)> {
        self.subgroups_impl(true)
            .into_iter()
            .map(|(subgroup, gens)| (NormalSubgroup::new_unchecked(subgroup), gens))
            .collect()
    }
}

pub fn direct_product_structure(
    group_one: &FiniteGroupMultiplicationTable,
    group_two: &FiniteGroupMultiplicationTable,
) -> FiniteGroupMultiplicationTable {
    let m = group_one.n;
    let n = group_two.n;

    let single_to_pair = |i: usize| -> (usize, usize) { (i % m, i / m) };
    let pair_to_single = |i: usize, j: usize| -> usize { i + j * m };

    FiniteGroupMultiplicationTable {
        n: m * n,
        ident: pair_to_single(group_one.ident, group_two.ident),
        inv: (0..m * n)
            .map(|x| {
                let (i, j) = single_to_pair(x);
                pair_to_single(group_one.inv[i], group_two.inv[j])
            })
            .collect(),
        mul: (0..m * n)
            .map(|x| {
                (0..m * n)
                    .map(|y| {
                        let (ix, jx) = single_to_pair(x);
                        let (iy, jy) = single_to_pair(y);
                        pair_to_single(group_one.mul[ix][iy], group_two.mul[jx][jy])
                    })
                    .collect()
            })
            .collect(),
        conjugacy_classes: None,
        is_abelian: None,
        is_simple: None,
    }
}

pub mod examples {
    use crate::free_group::todd_coxeter::FinitelyGeneratedGroupPresentation;

    use super::{FiniteGroupMultiplicationTable, direct_product_structure};

    pub fn trivial_group_structure() -> FiniteGroupMultiplicationTable {
        cyclic_group_structure(1)
    }

    pub fn cyclic_group_structure(n: usize) -> FiniteGroupMultiplicationTable {
        FiniteGroupMultiplicationTable::from_raw_model_unchecked(
            (0..n).collect(),
            || 0,
            |x: usize| (n - x) % n,
            |x: usize, y: usize| (x + y) % n,
            Some(true),
            None,
        )
    }

    pub fn klein_four_structure() -> FiniteGroupMultiplicationTable {
        direct_product_structure(&cyclic_group_structure(2), &cyclic_group_structure(2))
    }

    pub fn dihedral_group_structure(n: usize) -> FiniteGroupMultiplicationTable {
        // dihedral group using the presentation
        // <a b : a^2 = b^2 = (ab)^n = e>
        assert!(1 <= n);

        let mut grp = FinitelyGeneratedGroupPresentation::new();
        let a = grp.add_generator();
        let b = grp.add_generator();
        grp.add_relation(a.pow(2));
        grp.add_relation(b.pow(2));
        grp.add_relation((&a * &b).pow(n as isize));
        let mut grp = grp.into_finite_group();
        grp.is_abelian = Some(n <= 2);
        grp.is_simple = Some(n <= 1);
        grp
    }

    pub fn quaternion_group_structure() -> FiniteGroupMultiplicationTable {
        // quaternion group using the presentation
        // <-1 i j k : (-1)^2 = 1  i^2 = j^2 = k^2 = ijk = -1>
        let mut grp = FinitelyGeneratedGroupPresentation::new();
        let a = grp.add_generator();
        let i = grp.add_generator();
        let j = grp.add_generator();
        let k = grp.add_generator();
        grp.add_relation(a.pow(2));
        grp.add_two_sided_relation(i.pow(2), a.clone());
        grp.add_two_sided_relation(j.pow(2), a.clone());
        grp.add_two_sided_relation(k.pow(2), a.clone());
        grp.add_two_sided_relation(&i * &j * &k, a.clone());
        let mut grp = grp.into_finite_group();
        grp.is_abelian = Some(false);
        grp.is_simple = Some(false);
        grp
    }

    pub fn symmetric_group_structure(n: usize) -> FiniteGroupMultiplicationTable {
        super::super::super::permutation::Permutation::symmetric_composition_table(n).0
    }

    pub fn alternating_group_structure(n: usize) -> FiniteGroupMultiplicationTable {
        super::super::super::permutation::Permutation::alternating_composition_table(n).0
    }
}

#[cfg(test)]
mod group_tests {
    use super::*;

    #[test]
    fn test_cyclic() {
        for k in [1, 2, 3, 81, 91, 97, 100, 128] {
            let mut grp = examples::cyclic_group_structure(k);
            grp.cache_conjugacy_classes();
            grp.check_state().unwrap();
            assert_eq!(grp.elems().len(), k);
        }
    }

    #[test]
    fn test_dihedral() {
        for k in [1, 2, 3, 16, 17, 18, 50] {
            let mut grp = examples::dihedral_group_structure(k);
            grp.cache_conjugacy_classes();
            grp.check_state().unwrap();
            assert_eq!(grp.elems().len(), 2 * k);
        }
    }

    #[test]
    fn test_quaternion() {
        let mut grp = examples::quaternion_group_structure();
        assert_eq!(grp.size(), 8);
        grp.cache_conjugacy_classes();
        grp.check_state().unwrap();
    }

    #[test]
    fn test_direct_product() {
        let grp1 = examples::dihedral_group_structure(5);
        let grp2 = examples::dihedral_group_structure(4);
        let grp3 = direct_product_structure(&grp1, &grp2);
        grp3.check_state().unwrap();
    }

    #[test]
    fn test_conjugacy_class_count() {
        for (grp, num_ccls) in vec![
            (examples::cyclic_group_structure(12), 12),
            (examples::cyclic_group_structure(13), 13),
            (examples::dihedral_group_structure(1), 2),
            (examples::dihedral_group_structure(12), 9),
            (examples::dihedral_group_structure(13), 8),
            (examples::symmetric_group_structure(1), 1), //only the identity
            (examples::symmetric_group_structure(2), 2), //ident, 2-cycle
            (examples::symmetric_group_structure(3), 3), //ident, 2-cycle, 3-cycle
            (examples::symmetric_group_structure(4), 5), //ident, 2-cycle, 3-cycle, 4-cycle, (2, 2)-cycle
            (examples::symmetric_group_structure(5), 7), //ident, 2-cycle, 3-cycle, 4-cycle, 5-cycle, (2, 2)-cycle, (2, 3)-cycle
                                                         // (examples::symmetric_group_structure(6), 11), //ident, 2, 3, 4, 5, 6, (2, 2), (2, 3), (2, 4), (3, 3), (2, 2, 2)
        ] {
            assert_eq!(grp.conjugacy_classes().size(), num_ccls);
        }
    }
}
