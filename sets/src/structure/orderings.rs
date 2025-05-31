use super::SetSignature;
use std::{borrow::Borrow, cmp::Ordering};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MergedSource {
    First,
    Second,
}

struct VecMerger<'s, X, O: OrdSignature, K: Fn(&X) -> &O::Set> {
    ordering: &'s O,
    i: usize,
    a: Vec<Option<X>>,
    j: usize,
    b: Vec<Option<X>>,
    key: K,
}

impl<'s, X, O: OrdSignature + 's, K: Fn(&X) -> &O::Set> VecMerger<'s, X, O, K> {
    fn new(ordering: &'s O, a: Vec<X>, b: Vec<X>, key: K) -> Self {
        Self {
            ordering,
            i: 0,
            a: a.into_iter().map(|x| Some(x)).collect(),
            j: 0,
            b: b.into_iter().map(|x| Some(x)).collect(),
            key,
        }
    }
}

impl<'s, X, O: OrdSignature + 's, K: Fn(&X) -> &O::Set> Iterator for VecMerger<'s, X, O, K> {
    type Item = (MergedSource, X);

    fn next(&mut self) -> Option<Self::Item> {
        match (self.i < self.a.len(), self.j < self.b.len()) {
            (true, true) => {
                match self.ordering.cmp(
                    (self.key)(self.a[self.i].as_ref().unwrap()).borrow(),
                    (self.key)(self.b[self.j].as_ref().unwrap()).borrow(),
                ) {
                    Ordering::Less | Ordering::Equal => {
                        let a_item = self.a[self.i].take().unwrap();
                        self.i += 1;
                        Some((MergedSource::First, a_item))
                    }
                    Ordering::Greater => {
                        let b_item = self.b[self.j].take().unwrap();
                        self.j += 1;
                        Some((MergedSource::Second, b_item))
                    }
                }
            }
            (true, false) => {
                let a_item = self.a[self.i].take().unwrap();
                self.i += 1;
                Some((MergedSource::First, a_item))
            }
            (false, true) => {
                let b_item = self.b[self.j].take().unwrap();
                self.j += 1;
                Some((MergedSource::Second, b_item))
            }
            (false, false) => None,
        }
    }
}

pub trait OrdSignature: SetSignature {
    fn cmp(&self, a: &Self::Set, b: &Self::Set) -> Ordering;

    fn is_sorted(&self, a: Vec<impl Borrow<Self::Set>>) -> bool {
        for i in 1..a.len() {
            if self.cmp(&a[i - 1].borrow(), &a[i].borrow()) == Ordering::Greater {
                return false;
            }
        }
        true
    }

    fn is_sorted_by<X>(&self, a: Vec<X>, key: impl Fn(&X) -> &Self::Set) -> bool {
        self.is_sorted(a.iter().map(key).collect())
    }

    fn merge_sorted<'s, S: Borrow<Self::Set> + 's>(
        &'s self,
        a: Vec<S>,
        b: Vec<S>,
    ) -> impl Iterator<Item = (MergedSource, S)> {
        self.merge_sorted_by(a, b, |x| (*x).borrow())
    }

    fn merge_sorted_by<'s, X>(
        &'s self,
        a: Vec<X>,
        b: Vec<X>,
        key: impl Fn(&X) -> &Self::Set,
    ) -> impl Iterator<Item = (MergedSource, X)> {
        debug_assert!(self.is_sorted_by(a.iter().collect(), |x| key(x)));
        debug_assert!(self.is_sorted_by(b.iter().collect(), |x| key(x)));
        VecMerger::new(self, a, b, key)
    }

    fn sort<S: Borrow<Self::Set>>(&self, mut a: Vec<S>) -> Vec<S> {
        match a.len() {
            0 | 1 => a,
            2 => match self.cmp(a[0].borrow(), a[1].borrow()) {
                Ordering::Less | Ordering::Equal => a,
                Ordering::Greater => {
                    a.swap(0, 1);
                    a
                }
            },
            n => {
                // merge sort
                let k = n / 2;

                let a1 = a.split_off(k);
                let a0 = a;

                let a0_sorted = self.sort(a0);
                let a1_sorted = self.sort(a1);

                self.merge_sorted(a0_sorted, a1_sorted)
                    .into_iter()
                    .map(|(_, s)| s)
                    .collect()
            }
        }
    }

    fn sort_by<X>(&self, mut a: Vec<X>, key: impl Fn(&X) -> &Self::Set) -> Vec<X> {
        match a.len() {
            0 | 1 => a,
            2 => match self.cmp(key(&a[0]).borrow(), key(&a[1]).borrow()) {
                Ordering::Less | Ordering::Equal => a,
                Ordering::Greater => {
                    a.swap(0, 1);
                    a
                }
            },
            n => {
                // merge sort
                let k = n / 2;

                let a1 = a.split_off(k);
                let a0 = a;

                let a0_sorted = self.sort_by(a0, &key);
                let a1_sorted = self.sort_by(a1, &key);

                self.merge_sorted_by(a0_sorted, a1_sorted, key)
                    .into_iter()
                    .map(|(_, s)| s)
                    .collect()
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::structure::*;

    use super::*;

    #[derive(Debug, Clone, PartialEq, Eq)]
    struct UsizeStructure {}

    impl Signature for UsizeStructure {}

    impl SetSignature for UsizeStructure {
        type Set = usize;

        fn is_element(&self, _x: &Self::Set) -> bool {
            true
        }
    }

    impl EqSignature for UsizeStructure {
        fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
            a == b
        }
    }

    impl OrdSignature for UsizeStructure {
        fn cmp(&self, a: &Self::Set, b: &Self::Set) -> Ordering {
            usize::cmp(a, b)
        }
    }

    #[test]
    fn ordering_structure() {
        let s = UsizeStructure {};

        assert_eq!(s.cmp(&2, &2), Ordering::Equal);
        assert_eq!(s.cmp(&1, &2), Ordering::Less);
        assert_eq!(s.cmp(&2, &1), Ordering::Greater);

        assert!(s.is_sorted(Vec::<usize>::new()));
        assert!(s.is_sorted(vec![7]));
        assert!(s.is_sorted(vec![1, 2]));
        assert!(!s.is_sorted(vec![2, 1]));
        assert!(s.is_sorted(vec![1, 2, 3]));
        assert!(s.is_sorted(vec![1, 1, 1]));
        assert!(!s.is_sorted(vec![1, 2, 1]));

        assert_eq!(
            s.merge_sorted(Vec::<usize>::new(), Vec::<usize>::new())
                .collect::<Vec<_>>(),
            vec![]
        );
        {
            let a = vec![2, 4, 5, 7];
            let b = vec![1, 3, 6, 7, 8];
            let c = s.merge_sorted(a, b).collect::<Vec<_>>();
            assert_eq!(
                c,
                vec![
                    (MergedSource::Second, 1),
                    (MergedSource::First, 2),
                    (MergedSource::Second, 3),
                    (MergedSource::First, 4),
                    (MergedSource::First, 5),
                    (MergedSource::Second, 6),
                    (MergedSource::First, 7),
                    (MergedSource::Second, 7),
                    (MergedSource::Second, 8)
                ]
            );
        }

        assert_eq!(s.sort(Vec::<usize>::new()), Vec::<usize>::new());
        assert_eq!(s.sort(vec![7]), vec![7]);
        assert_eq!(s.sort(vec![7, 7]), vec![7, 7]);
        assert_eq!(s.sort(vec![1, 2]), vec![1, 2]);
        assert_eq!(s.sort(vec![2, 1]), vec![1, 2]);
        assert_eq!(s.sort(vec![7, 6, 5, 4, 3, 2, 1]), vec![1, 2, 3, 4, 5, 6, 7]);
        assert_eq!(
            s.sort(vec![3, 3, 2, 2, 2, 1, 1, 1, 1]),
            vec![1, 1, 1, 1, 2, 2, 2, 3, 3]
        );
    }
}
