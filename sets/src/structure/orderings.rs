use super::EqSignature;
use algebraeon_macros::{signature_meta_trait, skip_meta};
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

use super::MetaType;

#[signature_meta_trait]
pub trait PartialOrdSignature: EqSignature {
    fn partial_cmp(&self, a: &Self::Set, b: &Self::Set) -> Option<Ordering>;
}

#[signature_meta_trait]
pub trait OrdSignature: PartialOrdSignature {
    fn cmp(&self, a: &Self::Set, b: &Self::Set) -> Ordering;

    fn max(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        let c = self.cmp(a, b);
        match c {
            Ordering::Less | Ordering::Equal => b.clone(),
            Ordering::Greater => a.clone(),
        }
    }

    fn min(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        let c = self.cmp(a, b);
        match c {
            Ordering::Less | Ordering::Equal => a.clone(),
            Ordering::Greater => b.clone(),
        }
    }

    fn is_sorted(&self, a: &[impl Borrow<Self::Set>]) -> bool {
        for i in 1..a.len() {
            if self.cmp(a[i - 1].borrow(), a[i].borrow()) == Ordering::Greater {
                return false;
            }
        }
        true
    }

    fn is_sorted_by_key<X>(&self, a: &[X], key: impl Fn(&X) -> &Self::Set) -> bool {
        self.is_sorted(&a.iter().map(key).collect::<Vec<_>>())
    }

    fn is_sorted_and_unique(&self, a: &[impl Borrow<Self::Set>]) -> bool {
        for i in 1..a.len() {
            if self.cmp(a[i - 1].borrow(), a[i].borrow()) != Ordering::Less {
                return false;
            }
        }
        true
    }

    fn is_sorted_and_unique_by_key<X>(&self, a: &[X], key: impl Fn(&X) -> &Self::Set) -> bool {
        self.is_sorted_and_unique(&a.iter().map(key).collect::<Vec<_>>())
    }

    fn binary_search(&self, v: &[impl Borrow<Self::Set>], target: &Self::Set) -> bool {
        self.binary_search_by_key(v, target, |x| x.borrow())
            .is_some()
    }

    fn binary_search_by_key<'x, X>(
        &self,
        v: &'x [X],
        target: &Self::Set,
        key: impl Fn(&X) -> &Self::Set,
    ) -> Option<&'x X> {
        debug_assert!(self.is_sorted_by_key(v, &key));
        if v.is_empty() {
            return None;
        }
        let mut a = 0;
        let mut b = v.len() - 1;

        if self.equal(target, key(&v[a])) {
            return Some(&v[a]);
        }

        if self.equal(target, key(&v[b])) {
            return Some(&v[b]);
        }

        while b - a >= 2 {
            println!("{:?} {:?}", a, b);

            let m = (a + b) / 2;
            let m_key = key(&v[m]);
            match self.cmp(target, m_key) {
                Ordering::Less => b = m,
                Ordering::Equal => {
                    return Some(&v[m]);
                }
                Ordering::Greater => a = m,
            }
        }

        None
    }

    #[skip_meta]
    fn merge_sorted<'s, S: Borrow<Self::Set> + 's>(
        &'s self,
        a: Vec<S>,
        b: Vec<S>,
    ) -> impl Iterator<Item = (MergedSource, S)> {
        self.merge_sorted_by_key(a, b, |x| (*x).borrow())
    }

    #[skip_meta]
    fn merge_sorted_by_key<X>(
        &self,
        a: Vec<X>,
        b: Vec<X>,
        key: impl Fn(&X) -> &Self::Set,
    ) -> impl Iterator<Item = (MergedSource, X)> {
        debug_assert!(self.is_sorted_by_key(&a.iter().collect::<Vec<_>>(), |x| key(x)));
        debug_assert!(self.is_sorted_by_key(&b.iter().collect::<Vec<_>>(), |x| key(x)));
        VecMerger::new(self, a, b, key)
    }

    fn sort<S: Borrow<Self::Set>>(&self, a: Vec<S>) -> Vec<S> {
        self.sort_by_key(a, &|x| x.borrow())
    }

    #[skip_meta]
    fn sort_by_key<X>(&self, mut a: Vec<X>, key: &impl Fn(&X) -> &Self::Set) -> Vec<X> {
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

                let a0_sorted = self.sort_by_key(a0, key);
                let a1_sorted = self.sort_by_key(a1, key);

                self.merge_sorted_by_key(a0_sorted, a1_sorted, key)
                    .map(|(_, s)| s)
                    .collect()
            }
        }
    }

    fn sort_by_cached_by<X: 'static>(&self, a: Vec<X>, key: impl Fn(&X) -> Self::Set) -> Vec<X>
    where
        Self::Set: 'static,
    {
        let mut a = a.into_iter().map(|x| (x, None)).collect::<Vec<_>>();
        for i in 0..a.len() {
            let (x, k) = a.get_mut(i).unwrap();
            *k = Some(key(x));
        }
        self.sort_by_key(a, &|(_, k)| k.as_ref().unwrap())
            .into_iter()
            .map(|(x, _)| x)
            .collect()
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

        fn is_element(&self, _x: &Self::Set) -> Result<(), String> {
            Ok(())
        }
    }

    impl EqSignature for UsizeStructure {
        fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
            a == b
        }
    }

    impl PartialOrdSignature for UsizeStructure {
        fn partial_cmp(&self, a: &Self::Set, b: &Self::Set) -> Option<Ordering> {
            Some(self.cmp(a, b))
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

        assert!(s.is_sorted(&Vec::<usize>::new()));
        assert!(s.is_sorted(&[7]));
        assert!(s.is_sorted(&[1, 2]));
        assert!(!s.is_sorted(&[2, 1]));
        assert!(s.is_sorted(&[1, 2, 3]));
        assert!(s.is_sorted(&[1, 1, 1]));
        assert!(!s.is_sorted(&[1, 2, 1]));

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
        assert_eq!(
            s.sort_by_cached_by(vec![3, 3, 2, 2, 2, 1, 1, 1, 1], |x| 5 - x),
            vec![3, 3, 2, 2, 2, 1, 1, 1, 1]
        );

        let v = vec![1, 2, 3, 3, 3, 4, 4, 5, 7, 8, 9, 10, 11];
        assert!(!s.binary_search(&v, &0));
        assert!(s.binary_search(&v, &1));
        assert!(s.binary_search(&v, &3));
        assert!(s.binary_search(&v, &4));
        assert!(s.binary_search(&v, &5));
        assert!(!s.binary_search(&v, &6));
        assert!(s.binary_search(&v, &7));
        assert!(s.binary_search(&v, &11));
        assert!(!s.binary_search(&v, &12));
    }
}
