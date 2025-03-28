/// Provides an iterator over the k-element subsets of {0, 1, ..., n-1} in lexographic order such that some elements of {0, 1, ..., n-1} can be excluded from future subsets at any point during the iteration.
#[derive(Debug)]
pub struct LexicographicCombinationsWithRemovals {
    n: usize,
    all_items_idx: Vec<Option<usize>>,
    remaining_items: Vec<usize>,
    subset: Vec<usize>,
    finished: bool,
}

impl LexicographicCombinationsWithRemovals {
    /// Constructor
    pub fn new(n: usize, k: usize) -> Self {
        Self {
            n,
            all_items_idx: (0..n).map(|x| Some(x)).collect(),
            remaining_items: (0..n).collect(),
            subset: (0..k).collect(),
            finished: k > n,
        }
    }

    /// Exclude the provided element from appearing in any further subsets.
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
            true => {
                return None;
            }
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

/// Returns all size k subsets of {0, 1, ..., n-1}.
/// ```
/// use algebraeon_sets::combinatorics::combinations;
/// assert_eq!(combinations(5, 3).collect::<Vec<_>>(), vec![
///     vec![0, 1, 2],
///     vec![0, 1, 3],
///     vec![0, 1, 4],
///     vec![0, 2, 3],
///     vec![0, 2, 4],
///     vec![0, 3, 4],
///     vec![1, 2, 3],
///     vec![1, 2, 4],
///     vec![1, 3, 4],
///     vec![2, 3, 4],
/// ]);
/// ```
pub fn combinations(n: usize, k: usize) -> impl Iterator<Item = Vec<usize>> {
    println!("{} {}", n, k);
    LexicographicCombinationsWithRemovals::new(n, k)
}

/// Returns all size k subsets of items.
/// /// ```
/// use algebraeon_sets::combinatorics::subsets;
/// assert_eq!(subsets(vec!["a", "b", "c"], 2).collect::<Vec<_>>(), vec![
///     vec!["a", "b"],
///     vec!["a", "c"],
///     vec!["b", "c"],
/// ]);
/// ```
pub fn subsets<'a, T: 'a + Clone>(items: Vec<T>, k: usize) -> impl 'a + Iterator<Item = Vec<T>> {
    combinations(items.len(), k)
        .map(move |subset| subset.into_iter().map(|idx| items[idx].clone()).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    pub fn lexographic_subsets_with_removals_test_edge_cases() {
        let x = combinations(0, 1).collect::<Vec<_>>();
        println!("{:?}", x);
        assert_eq!(x.len(), 0);

        let x = combinations(3, 5).collect::<Vec<_>>();
        println!("{:?}", x);
        assert_eq!(x.len(), 0);

        let x = combinations(3, 3).collect::<Vec<_>>();
        println!("{:?}", x);
        assert_eq!(x.len(), 1);
    }

    #[test]
    pub fn lexographic_subsets_with_removals_test_1() {
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
    pub fn lexographic_subsets_with_removals_test_2() {
        let mut c = LexicographicCombinationsWithRemovals::new(7, 3);
        c.exclude(0);
        c.exclude(1);
        c.exclude(2);
        c.exclude(3);
        c.exclude(4);
        debug_assert_eq!(c.next(), None);
    }

    #[test]
    pub fn run() {
        println!("{:?}", combinations(5, 3).collect::<Vec<_>>());
        assert_eq!(combinations(5, 3).collect::<Vec<_>>().len(), 10)
    }
}
