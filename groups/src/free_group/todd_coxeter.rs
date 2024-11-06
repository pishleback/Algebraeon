use core::num;

use crate::{group::Group, permutation::*};

#[derive(Clone, Copy)]
enum Neighbor {
    None(),
    Coset(usize),
}

struct SchreierGraph {
    num_gens: usize,
    idents: Vec<usize>,
    neighbors: Vec<Vec<Neighbor>>,
}

impl SchreierGraph {
    pub fn new(num_gens: usize) -> Self {
        Self {
            num_gens,
            idents: vec![],
            neighbors: vec![],
        }
    }

    fn find_coset(&self, mut c: usize) -> usize {
        while self.idents[c] != c {
            c = self.idents[c];
        }
        return c;
    }

    fn new_coset(&mut self) -> usize {
        let c = self.idents.len();
        self.idents.push(c);
        let mut new_nbs = vec![];
        for _ in 0..2 * self.num_gens {
            new_nbs.push(Neighbor::None());
        }
        self.neighbors.push(new_nbs);
        return c;
    }

    fn unify(&mut self, c1: usize, c2: usize) {
        let mut c1 = self.find_coset(c1);
        let mut c2 = self.find_coset(c2);
        if c1 == c2 {
            return;
        }
        if c2 < c1 {
            (c1, c2) = (c2, c1);
        }
        self.idents[c2] = c1;
        for d in 0..2 * self.num_gens {
            let n1 = self.neighbors[c1][d];
            let n2 = self.neighbors[c2][d];
            match n1 {
                Neighbor::None() => {
                    self.neighbors[c1][d] = n2;
                }
                Neighbor::Coset(n1c) => match n2 {
                    Neighbor::None() => {}
                    Neighbor::Coset(n2c) => {
                        self.unify(n1c, n2c);
                    }
                },
            }
        }
    }

    fn follow(&mut self, mut c: usize, d: usize) -> usize {
        c = self.find_coset(c);
        match self.neighbors[c][d] {
            Neighbor::None() => {
                let nc = self.new_coset();
                self.neighbors[c][d] = Neighbor::Coset(nc);
                return nc;
            }
            Neighbor::Coset(nc) => {
                return self.find_coset(nc);
            }
        }
    }

    fn follow_path(&mut self, mut c: usize, ds: &Vec<usize>) -> usize {
        c = self.find_coset(c);
        for d in ds.iter().rev() {
            c = self.follow(c, *d);
        }
        return c;
    }
}

/*
n = num_gens
now a generator is represented by a signed integer.
The generators are 0, 2, 4, ..., 2n-2 with inverses 1, 3, 5, ..., 2n-1
 */
pub fn enumerate_cosets(
    num_gens: usize,
    rels: Vec<Vec<usize>>,
    subgens: Vec<Vec<usize>>,
) -> (usize, Vec<Permutation>) {
    //impose inverse relations
    let mut full_rels = rels.clone();
    for i in 0..num_gens {
        full_rels.push(vec![2 * i + 1, 2 * i]);
    }
    let rels = full_rels;

    //check rels
    for word in &rels {
        for a in word {
            assert!(*a < 2 * num_gens);
        }
    }
    //check subgens
    for word in &subgens {
        for a in word {
            assert!(*a < 2 * num_gens);
        }
    }

    let mut scg = SchreierGraph::new(num_gens);
    let start = scg.new_coset();

    for subgen in subgens {
        let end = scg.follow_path(start, &subgen);
        scg.unify(end, start);
    }

    let mut to_visit = 0;
    while to_visit < scg.idents.len() {
        let c = scg.find_coset(to_visit);
        if c == to_visit {
            for rel in &rels {
                let b = scg.follow_path(c, rel);
                scg.unify(b, c);
            }
        }
        to_visit += 1
    }

    enum CosetIndexEntry {
        None,
        Index(usize),
    }

    let mut coset_index_lookup = vec![];
    let mut cosets = vec![];
    for (i, c) in scg.idents.iter().enumerate() {
        if i == *c {
            coset_index_lookup.push(CosetIndexEntry::Index(cosets.len()));
            cosets.push(*c);
        } else {
            coset_index_lookup.push(CosetIndexEntry::None);
        }
    }

    let mut perms = vec![];
    for g in 0..num_gens {
        perms.push(
            Permutation::new(
                cosets
                    .iter()
                    .map(|c| match coset_index_lookup[scg.follow(*c, 2 * g)] {
                        CosetIndexEntry::None => {
                            panic!()
                        }
                        CosetIndexEntry::Index(i) => {
                            debug_assert_eq!(cosets[i], scg.follow(*c, 2 * g));
                            i
                        }
                    })
                    .collect(),
            )
            .unwrap(),
        );
    }

    return (cosets.len(), perms);
}

pub fn enumerate_elements(num_gens: usize, rels: Vec<Vec<usize>>) -> (usize, Vec<Permutation>) {
    enumerate_cosets(num_gens, rels, vec![])
}

pub fn enumerate_group(
    num_gens: usize,
    rels: Vec<Vec<usize>>,
) -> super::super::composition_table::group::Group {
    let (n, gen_perms) = enumerate_cosets(num_gens, rels, vec![]);
    let inv_gen_perms = gen_perms
        .iter()
        .map(|perm| perm.inverse_ref())
        .collect::<Vec<Permutation>>();

    //write each element as a word in the n generators
    let mut paths: Vec<(bool, Vec<usize>)> = vec![];
    for i in 0..n {
        paths.push((i == 0, vec![]));
    }

    let mut boundary = vec![0];
    let mut new_boundary = vec![];
    while boundary.len() > 0 {
        for b_idx in boundary {
            let (b_done, b_path) = paths[b_idx].clone();
            debug_assert!(b_done);
            for g in 0..num_gens {
                let c = gen_perms[g].call(b_idx);
                if !paths[c].0 {
                    let mut c_path = b_path.clone();
                    c_path.push(g);
                    paths[c] = (true, c_path);
                    new_boundary.push(c);
                }
            }
        }
        boundary = new_boundary;
        new_boundary = vec![];
    }

    //remove the done flags from paths and check that the values make sense
    for p in &paths {
        assert!(p.0);
    }
    let paths = paths
        .into_iter()
        .map(|(_done, path)| path)
        .collect::<Vec<Vec<usize>>>();
    for p in &paths {
        for g in p {
            debug_assert!(*g < num_gens);
        }
    }

    super::super::composition_table::group::Group::new_unchecked(
        n,
        0,
        (0..n)
            .map(|x| {
                let mut y = 0;
                for g in paths[x].iter().rev() {
                    y = inv_gen_perms[*g].call(y);
                }
                y
            })
            .collect(),
        (0..n)
            .map(|x| {
                (0..n)
                    .map(|y| {
                        let mut z = 0;
                        for g in paths[x].iter() {
                            z = gen_perms[*g].call(z);
                        }
                        for g in paths[y].iter() {
                            z = gen_perms[*g].call(z);
                        }
                        z
                    })
                    .collect()
            })
            .collect(),
        None,
        None,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trivial() {
        let (n, perms) = enumerate_elements(0, vec![]);
        assert_eq!(n, 1);
        assert_eq!(perms.len(), 0);
    }

    #[test]
    fn alternating_5() {
        let (n, _perms) = enumerate_elements(
            2,
            vec![vec![0, 0, 0, 0, 0], vec![2, 2], vec![0, 2, 0, 2, 0, 2]],
        );
        assert_eq!(n, 60);

        let (n, _perms) = enumerate_cosets(
            2,
            vec![vec![0, 0, 0, 0, 0], vec![2, 2], vec![0, 2, 0, 2, 0, 2]],
            vec![vec![0, 2]],
        );
        assert_eq!(n, 20);
    }

    #[test]
    fn small_coxeter() {
        //icosahedral symmetry (3, 5)
        let (n, _perms) = enumerate_elements(
            3,
            vec![
                vec![0, 0],
                vec![2, 2],
                vec![4, 4],
                vec![0, 2, 0, 2, 0, 2, 0, 2, 0, 2],
                vec![2, 4, 2, 4, 2, 4],
                vec![0, 4, 0, 4],
            ],
        );
        assert_eq!(n, 120);
    }

    #[test]
    fn large_coxeter() {
        //600-cell symmetry (3, 3, 5)
        let (n, _perms) = enumerate_elements(
            4,
            vec![
                vec![0, 0],
                vec![2, 2],
                vec![4, 4],
                vec![6, 6],
                vec![0, 2, 0, 2, 0, 2, 0, 2, 0, 2],
                vec![2, 4, 2, 4, 2, 4],
                vec![4, 6, 4, 6, 4, 6],
                vec![0, 4, 0, 4],
                vec![0, 6, 0, 6],
                vec![2, 6, 2, 6],
            ],
        );
        assert_eq!(n, 14400);
    }

    #[test]
    fn unexpected_trivial() {
        // <a, b | bab^-1=a^2, aba=bab> is the trivial group
        let (n, _perms) = enumerate_elements(2, vec![vec![2, 0, 3, 1, 1], vec![0, 2, 0, 3, 1, 3]]);
        assert_eq!(n, 1);
    }
}
