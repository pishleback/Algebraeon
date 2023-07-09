#![allow(dead_code)]

use std::rc::Rc;

use super::function::*;
use super::set::*;

struct Sequences<N: SetT, X: SetT> {
    set_n: N,
    set_x: X,
}

impl<N: SetT, X: SetT> Sequences<N, X> {
    pub fn to_set(&self) -> (Rc<Set>, impl Iterator<Item = Function<N, X>>) {
        let size = self.size();
        let set = Rc::new(Set::new(size));

        struct SequenceIterator<N: SetT, X: SetT> {
            set_n: N,
            set_x: X,
            counters: Vec<usize>,
            done: bool,
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

        (
            set.clone(),
            SequenceIterator::<N, X> {
                set_n: self.set_n.clone(),
                set_x: self.set_x.clone(),
                done: false,
                counters: vec![0; self.set_n.borrow().size()],
            },
        )
    }

    pub fn size(&self) -> usize {
        let n = self.set_n.borrow().size();
        let x = self.set_x.borrow().size();
        let mut ans = 1;
        for _i in 0..n {
            ans *= x;
        }
        ans
    }
}

pub fn twelvefold<N: SetT, X: SetT>(
    set_n: N,
    set_x: X,
    injective: bool,
    surjective: bool,
    unordered_n: bool,
    unordered_x: bool,
) -> impl Iterator<Item = Function<N, X>> {
    match (injective, surjective, unordered_n, unordered_x) {
        (false, false, false, false) => Sequences { set_n, set_x }.to_set().1,
        (true, false, false, false) => Sequences { set_n, set_x }.to_set().1,
        _ => Sequences { set_n, set_x }.to_set().1,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    pub fn test_sequences() {
        let set_n = Set::new(3);
        let set_x = Set::new(5);
        let seq = Sequences {
            set_n: &set_n,
            set_x: &set_x,
        };
        assert_eq!(seq.size(), 125);
        let (seq_set, seqs_iter) = seq.to_set();
        let seqs: Vec<Function<&Set, &Set>> = seqs_iter.collect();
        assert_eq!(seq.size(), seq_set.size());
        assert_eq!(seq.size(), seqs.len());
        for f in seqs {
            f.check_state().unwrap();
        }
    }

    #[test]
    pub fn test_twelvefold() {
        {
            let set_n = Set::new(3);
            let set_x = Set::new(5);
            let funcs: Vec<Function<&Set, &Set>> =
                twelvefold(&set_n, &set_x, false, false, false, false).collect();
            assert_eq!(5 * 5 * 5, funcs.len());
        }

        {
            let set_n = Set::new(3);
            let set_x = Set::new(5);
            let funcs: Vec<Function<&Set, &Set>> =
                twelvefold(&set_n, &set_x, true, false, false, false).collect();
            assert_eq!(5 * 4 * 3, funcs.len());
        }
    }
}
