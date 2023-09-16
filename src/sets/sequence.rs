use super::function::*;
use super::set::*;

pub fn count_sequences(n: usize, x: usize) -> usize {
    let mut ans = 1;
    for _i in 0..n {
        ans *= x;
    }
    ans
}

struct SequenceIterator<N: SetT, X: SetT> {
    set_n: N,
    set_x: X,
    counters: Vec<usize>,
    done: bool,
}

impl<N: SetT, X: SetT> SequenceIterator<N, X> {
    pub fn new(set_n: N, set_x: X) -> Self {
        let size_n = set_n.borrow().size();
        SequenceIterator {
            set_n,
            set_x,
            counters: vec![0; size_n],
            done: false,
        }
    }
}

impl<N: SetT, X: SetT> Iterator for SequenceIterator<N, X> {
    type Item = Function<N, X>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        let next = Some(Function::new_unchecked(
            self.set_n.clone(),
            self.set_x.clone(),
            self.counters.clone(),
        ));

        'exit_for: {
            for i in 0..self.counters.len() {
                self.counters[i] += 1;
                if self.counters[i] == self.set_x.borrow().size() {
                    self.counters[i] = 0;
                } else {
                    break 'exit_for;
                }
            }
            self.done = true;
        }

        return next;
    }
}

pub fn count_distinct_sequences(n: usize, x: usize) -> usize {
    let mut ans = 1;
    for i in 0..n {
        ans *= x - i;
        if ans == 0 {
            return 0;
        }
    }
    ans
}

struct DistinctSequenceIterator<N: SetT, X: SetT> {
    set_n: N,
    set_x: X,
    counters: Vec<usize>,
    done: bool,
}

impl<N: SetT, X: SetT> DistinctSequenceIterator<N, X> {
    pub fn new(set_n: N, set_x: X) -> Self {
        let size_n = set_n.borrow().size();
        DistinctSequenceIterator {
            set_n,
            set_x,
            counters: vec![0; size_n],
            done: false,
        }
    }
}

impl<N: SetT, X: SetT> Iterator for DistinctSequenceIterator<N, X> {
    type Item = Function<N, X>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.set_n.borrow().size() > self.set_x.borrow().size() {
            return None;
        }

        if self.done {
            return None;
        }

        fn multiindex_inj(n: usize, multi_idx: &Vec<usize>) -> Vec<usize> {
            let mut possible_img: Vec<usize> = (0..n).collect();
            let mut func = vec![];
            for idx in multi_idx {
                func.push(possible_img.remove(*idx));
            }
            return func;
        }

        let next = Some(Function::new_unchecked(
            self.set_n.clone(),
            self.set_x.clone(),
            multiindex_inj(self.set_x.borrow().size(), &self.counters),
        ));

        'exit_for: {
            for i in 0..self.counters.len() {
                self.counters[i] += 1;
                if self.counters[i] == self.set_x.borrow().size() - i {
                    self.counters[i] = 0;
                } else {
                    break 'exit_for;
                }
            }
            self.done = true;
        }

        return next;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    pub fn test_sequences() {
        let set_n = Set::new(3);
        let set_x = Set::new(5);
        let seqs: Vec<Function<&Set, &Set>> = SequenceIterator::new(&set_n, &set_x).collect();
        assert_eq!(
            seqs.len(),
            count_sequences(set_n.size(), set_x.size())
        );
        for f in seqs {
            f.check_state().unwrap();
        }
    }

    #[test]
    pub fn test_distinct_sequences() {
        let set_n = Set::new(3);
        let set_x = Set::new(5);
        let seqs: Vec<Function<&Set, &Set>> =
            DistinctSequenceIterator::new(&set_n, &set_x).collect();
        assert_eq!(
            seqs.len(),
            count_distinct_sequences(set_n.size(), set_x.size())
        );
        for f in seqs {
            f.check_state().unwrap();
        }

        let set_n = Set::new(5);
        let set_x = Set::new(3);
        let seqs: Vec<Function<&Set, &Set>> =
            DistinctSequenceIterator::new(&set_n, &set_x).collect();
        assert_eq!(
            seqs.len(),
            count_distinct_sequences(set_n.size(), set_x.size())
        );
        for f in seqs {
            f.check_state().unwrap();
        }
    }
}
