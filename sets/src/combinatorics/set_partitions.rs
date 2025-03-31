use std::collections::{HashMap, HashSet};

use indexmap::IndexMap;

#[derive(Clone, Debug)]
pub struct Partition {
    partition: Vec<HashSet<usize>>, // the partition
    lookup: Vec<usize>,             // for each element, the index of its part
}

impl Partition {
    fn check_state(&self) -> Result<(), &'static str> {
        let mut present = HashMap::new();
        let n = self.lookup.len();
        for (idx, part) in self.partition.iter().enumerate() {
            if part.len() == 0 {
                return Err("Partition contains an empty part");
            }
            for &x in part {
                if n <= x {
                    return Err("Partition contains element which is too big");
                }
                if present.contains_key(&x) {
                    return Err("Duplicate element in partition");
                }
                present.insert(x, idx);
            }
        }
        for x in 0..n {
            if !present.contains_key(&x) {
                return Err("Missing element from partition");
            }
            if present.get(&x).unwrap() != &self.lookup[x] {
                return Err("Incorrect entry in lookup");
            }
        }
        Ok(())
    }

    pub fn new_unchecked(partition: Vec<HashSet<usize>>, lookup: Vec<usize>) -> Self {
        let partition = Self { partition, lookup };
        #[cfg(debug_assertions)]
        partition.check_state().unwrap();
        partition
    }

    pub fn new_from_function<T: Clone + Eq + std::hash::Hash>(
        n: usize,
        f: impl Fn(usize) -> T,
    ) -> (Self, Vec<T>) {
        let mut t_lookup = vec![];
        for x in 0..n {
            t_lookup.push(f(x));
        }
        let mut t_partition = IndexMap::new();
        for x in 0..n {
            let t = &t_lookup[x];
            if !t_partition.contains_key(&t) {
                t_partition.insert(t, vec![x]);
            } else {
                t_partition.get_mut(&t).unwrap().push(x)
            }
        }

        let lookup = (0..n)
            .map(|x| t_partition.get_index_of(&t_lookup[x]).unwrap())
            .collect();
        let partition = t_partition
            .iter()
            .map(|(_t, part)| part.iter().cloned().collect())
            .collect();

        let partition = Partition::new_unchecked(partition, lookup);
        #[cfg(debug_assertions)]
        partition.check_state().unwrap();
        (
            partition,
            t_partition
                .into_iter()
                .map(|(t, _part)| t.clone())
                .collect(),
        )
    }

    pub fn project(&self, x: usize) -> usize {
        self.lookup[x]
    }

    pub fn class_containing(&self, x: usize) -> &HashSet<usize> {
        self.get_class(self.project(x))
    }

    pub fn get_class(&self, i: usize) -> &HashSet<usize> {
        &self.partition[i]
    }

    pub fn num_elements(&self) -> usize {
        self.lookup.len()
    }

    pub fn num_classes(&self) -> usize {
        self.partition.len()
    }

    pub fn size(&self) -> usize {
        self.partition.len()
    }
}

#[cfg(test)]
mod partition_tests {
    use super::*;

    #[test]
    fn partition_check_bad_state() {
        //not a covering set
        let p = Partition {
            partition: vec![
                vec![0, 2].into_iter().collect(),
                vec![3, 5].into_iter().collect(),
            ],
            lookup: vec![0, 0, 0, 1, 1, 1],
        };
        match p.check_state() {
            Ok(()) => assert!(false),
            Err(_) => {}
        }

        //not disjoint
        let p = Partition {
            partition: vec![
                vec![0, 1, 2, 3].into_iter().collect(),
                vec![2, 3, 4, 5].into_iter().collect(),
            ],
            lookup: vec![0, 0, 0, 0, 1, 1],
        };
        match p.check_state() {
            Ok(()) => assert!(false),
            Err(_) => {}
        }

        //lookup values too big
        let p = Partition {
            partition: vec![
                vec![0, 1, 2].into_iter().collect(),
                vec![3, 4, 5].into_iter().collect(),
            ],
            lookup: vec![0, 0, 0, 1, 1, 2],
        };
        match p.check_state() {
            Ok(()) => assert!(false),
            Err(_) => {}
        }

        //incorrect lookup values
        let p = Partition {
            partition: vec![
                vec![0, 1, 2].into_iter().collect(),
                vec![3, 4, 5].into_iter().collect(),
            ],
            lookup: vec![0, 0, 1, 1, 1, 1],
        };
        match p.check_state() {
            Ok(()) => assert!(false),
            Err(_) => {}
        }
    }

    #[test]
    fn from_function() {
        let (p, _ts) = Partition::new_from_function(6, |x| x % 2);
        println!("p = {:?}", p);
        assert_eq!(p.num_elements(), 6);
        assert_eq!(p.num_classes(), 2);
    }
}
