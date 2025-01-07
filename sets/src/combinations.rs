/// Provides an iterator over the k-element subsets of an n-element set in lexographic order such that elements of the set can be excluded part way through.
#[derive(Debug)]
pub struct LexicographicCombinationsWithRemovals {
    n: usize,
    all_items_idx: Vec<Option<usize>>,
    remaining_items: Vec<usize>,
    subset: Vec<usize>,
    finished: bool,
}


impl LexicographicCombinationsWithRemovals {
    pub fn new(n: usize, k: usize) -> Self {
        Self {
            n,
            all_items_idx: (0..n).map(|x| Some(x)).collect(),
            remaining_items: (0..n).collect(),
            subset: (0..k).collect(),
            finished: k >= n,
        }
    }

    /// Exclude a from any further iterations.
    pub fn exclude(&mut self, a: usize) {
        assert!(a < self.all_items_idx.len());
        if self.all_items_idx[a].is_none() {
            // Already excluded
            return;
        }

        while self.subset.iter().any(|x| self.remaining_items[*x] == a) {
            match self.next() {
                Some(_) => {}
                None => {
                    break;
                }
            }
        }

        if !self.finished {
            let a_idx = self.all_items_idx[a].unwrap();
            self.remaining_items.remove(a_idx);
            self.all_items_idx[a] = None;
            for i in (a + 1)..self.all_items_idx.len() {
                match self.all_items_idx[i].as_mut() {
                    Some(idx) => {
                        *idx -= 1;
                    }
                    None => {}
                }
            }
            for x in self.subset.iter_mut() {
                match (*x).cmp(&a_idx) {
                    std::cmp::Ordering::Less => {}
                    std::cmp::Ordering::Equal => {
                        unreachable!()
                    }
                    std::cmp::Ordering::Greater => {
                        *x -= 1;
                    }
                }
            }
            self.n -= 1;
        }
    }
}

impl Iterator for LexicographicCombinationsWithRemovals {
    type Item = Vec<usize>;

    fn next(&mut self) -> Option<Self::Item> {
        let next = match self.finished {
            true => None,
            false => Some(
                self.subset
                    .iter()
                    .map(|i| self.remaining_items[*i])
                    .collect(),
            ),
        };

        let k = self.subset.len();
        let i = 'FIND_I: {
            for (i, s) in self.subset.iter().rev().enumerate() {
                if *s < self.n - i - 1 {
                    break 'FIND_I Some(i);
                }
            }
            None
        };
        match i {
            Some(i) => {
                self.subset[k - 1 - i] += 1;
                for j in 0..i {
                    self.subset[k - i + j] = self.subset[k - i - 1] + 1 + j;
                }
            }
            None => {
                self.finished = true;
            }
        }

        next
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    pub fn run() {
        let mut c = LexicographicCombinationsWithRemovals::new(7, 3);
        for _ in 0..19 {
            let x = c.next().unwrap();
            println!("{:?}", x);
        }

        println!("rm 4");
        c.exclude(4);
        println!("rm 5");
        c.exclude(5);

        for _ in 0..2 {
            let x = c.next().unwrap();
            println!("{:?}", x);
        }

        assert_eq!(c.next(), None);
    }

    #[test]
    pub fn run2() {
        let mut c = LexicographicCombinationsWithRemovals::new(7, 3);
        c.exclude(0);
        c.exclude(1);
        c.exclude(2);
        c.exclude(3);
        c.exclude(4);
        debug_assert_eq!(c.next(), None);
    }
}
