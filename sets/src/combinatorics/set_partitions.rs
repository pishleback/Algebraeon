use indexmap::IndexMap;
use std::collections::HashSet;

#[derive(Clone, Debug)]
pub struct Partition {
    partition: Vec<HashSet<usize>>, // the partition
    lookup: Vec<usize>,             // for each element, the index of its part
}

impl Partition {
    #[cfg(debug_assertions)]
    fn check_state(&self) -> Result<(), &'static str> {
        use std::collections::HashMap;
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

#[derive(Debug, Clone)]
pub struct Element {
    x: usize,
    cum_x: usize,
    pivot: bool,
}

#[derive(Debug, Clone)]
pub struct LexographicPartitionsNumPartsInRange {
    // how many elements in the set
    n: usize,
    // min and max number of parts in the partition
    min_x: usize,
    max_x: usize,
    elements: Vec<Element>,
    finished: bool,
}

impl LexographicPartitionsNumPartsInRange {
    #[cfg(debug_assertions)]
    fn check(&self) -> Result<(), ()> {
        // check invariants
        if !self.finished {
            assert_eq!(self.elements.len(), self.n);
            assert_eq!(self.elements[0].x, 0);
            assert_eq!(self.elements[0].cum_x, 0);
            assert_eq!(self.elements[0].pivot, true);
            let mut cum_max = 0;
            for i in 1..self.n {
                if self.elements[i].x <= cum_max {
                    assert_eq!(self.elements[i].cum_x, cum_max);
                    assert_eq!(self.elements[i].pivot, false);
                } else if self.elements[i].x == cum_max + 1 {
                    cum_max += 1;
                    assert_eq!(self.elements[i].cum_x, cum_max);
                    assert_eq!(self.elements[i].pivot, true);
                } else {
                    panic!();
                }
            }
            cum_max += 1;
            assert!(self.min_x <= cum_max);
            assert!(cum_max <= self.max_x);
        }
        Ok(())
    }

    pub fn new(n: usize, min_x: usize, max_x: usize) -> Self {
        let mut elements = vec![];
        for i in 0..n {
            elements.push(Element {
                x: 0,
                cum_x: 0,
                pivot: i == 0,
            })
        }
        let mut s = Self {
            n,
            min_x,
            max_x,
            elements,
            finished: false,
        };
        if (n == 0 && min_x > 0) || (n > 0 && max_x == 0) || (n < min_x) || (min_x > max_x) {
            s.finished = true;
        }
        if n > 0 {
            s.reset_tail(0);
        }
        s
    }

    fn reset_tail(&mut self, j: usize) {
        let cum_max_j = self.elements[j].cum_x;
        // if !self.elements[j].pivot {
        //     self.elements[j].part = 0;
        // }
        for i in (j + 1)..self.n {
            let rev_i = self.n - i;
            let x = if rev_i <= self.min_x {
                let x = self.min_x - rev_i;
                if x > cum_max_j { x } else { 0 }
            } else {
                0
            };
            self.elements[i] = Element {
                x,
                cum_x: if x == 0 { cum_max_j } else { x },
                pivot: x != 0,
            };
        }
        #[cfg(debug_assertions)]
        self.check().unwrap();
    }
}

impl Iterator for LexographicPartitionsNumPartsInRange {
    type Item = Vec<usize>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            None
        } else {
            let next = (0..self.n).map(|i| self.elements[i].x).collect();
            'SEARCH: {
                for i in (0..self.n).rev() {
                    if !self.elements[i].pivot {
                        let max = self.elements[i].cum_x;
                        let x = &mut self.elements[i].x;
                        if *x + 1 < self.max_x {
                            if *x < max {
                                *x += 1;
                                self.reset_tail(i);
                                break 'SEARCH;
                            } else if *x == max {
                                *x += 1;
                                self.elements[i].cum_x += 1;
                                self.elements[i].pivot = true;
                                self.reset_tail(i);
                                break 'SEARCH;
                            }
                        }
                    }
                }
                self.finished = true;
            }
            Some(next)
        }
    }
}

pub fn set_partitions_eq(n: usize, x: usize) -> impl Iterator<Item = Vec<usize>> {
    LexographicPartitionsNumPartsInRange::new(n, x, x)
}

pub fn set_partitions_le(n: usize, x: usize) -> impl Iterator<Item = Vec<usize>> {
    LexographicPartitionsNumPartsInRange::new(n, 0, x)
}

pub fn set_partitions_ge(n: usize, x: usize) -> impl Iterator<Item = Vec<usize>> {
    LexographicPartitionsNumPartsInRange::new(n, x, n)
}

pub fn set_partitions_range(
    n: usize,
    min_x: usize,
    max_x: usize,
) -> impl Iterator<Item = Vec<usize>> {
    LexographicPartitionsNumPartsInRange::new(n, min_x, max_x)
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

    #[test]
    fn generate_set_partitions() {
        assert_eq!(
            LexographicPartitionsNumPartsInRange::new(0, 0, 0)
                .collect::<Vec<_>>()
                .len(),
            1
        );
        assert_eq!(
            LexographicPartitionsNumPartsInRange::new(0, 1, 1)
                .collect::<Vec<_>>()
                .len(),
            0
        );
        assert_eq!(
            LexographicPartitionsNumPartsInRange::new(0, 2, 2)
                .collect::<Vec<_>>()
                .len(),
            0
        );
        assert_eq!(
            LexographicPartitionsNumPartsInRange::new(0, 3, 3)
                .collect::<Vec<_>>()
                .len(),
            0
        );

        assert_eq!(
            LexographicPartitionsNumPartsInRange::new(1, 0, 0)
                .collect::<Vec<_>>()
                .len(),
            0
        );
        assert_eq!(
            LexographicPartitionsNumPartsInRange::new(1, 1, 1)
                .collect::<Vec<_>>()
                .len(),
            1
        );
        assert_eq!(
            LexographicPartitionsNumPartsInRange::new(1, 2, 2)
                .collect::<Vec<_>>()
                .len(),
            0
        );
        assert_eq!(
            LexographicPartitionsNumPartsInRange::new(1, 3, 3)
                .collect::<Vec<_>>()
                .len(),
            0
        );

        assert_eq!(
            LexographicPartitionsNumPartsInRange::new(2, 0, 0)
                .collect::<Vec<_>>()
                .len(),
            0
        );
        assert_eq!(
            LexographicPartitionsNumPartsInRange::new(2, 1, 1)
                .collect::<Vec<_>>()
                .len(),
            1
        );
        assert_eq!(
            LexographicPartitionsNumPartsInRange::new(2, 2, 2)
                .collect::<Vec<_>>()
                .len(),
            1
        );
        assert_eq!(
            LexographicPartitionsNumPartsInRange::new(2, 3, 3)
                .collect::<Vec<_>>()
                .len(),
            0
        );

        assert_eq!(
            LexographicPartitionsNumPartsInRange::new(3, 0, 0)
                .collect::<Vec<_>>()
                .len(),
            0
        );
        assert_eq!(
            LexographicPartitionsNumPartsInRange::new(3, 1, 1)
                .collect::<Vec<_>>()
                .len(),
            1
        );
        assert_eq!(
            LexographicPartitionsNumPartsInRange::new(3, 2, 2)
                .collect::<Vec<_>>()
                .len(),
            3
        );
        assert_eq!(
            LexographicPartitionsNumPartsInRange::new(3, 3, 3)
                .collect::<Vec<_>>()
                .len(),
            1
        );

        assert_eq!(
            LexographicPartitionsNumPartsInRange::new(4, 5, 3)
                .collect::<Vec<_>>()
                .len(),
            0
        );
    }
}
