#![feature(prelude_import)]
#![allow(
    clippy::uninlined_format_args,
    clippy::to_string_trait_impl,
    clippy::missing_panics_doc,
    clippy::must_use_candidate
)]
#[macro_use]
extern crate std;
#[prelude_import]
use std::prelude::rust_2024::*;
pub mod approximations {
    use crate::structure::{OrdSignature, SetSignature};
    /// A set of subsets of some set.
    pub trait SubsetsSignature: SetSignature {}
    /// A set of points in some space, defined in terms of a sequence of subsets whose intersection is the point being approximated.
    pub trait ApproximatePointsSignature: SetSignature {
        type Precision: OrdSignature;
        type OpenSubsetsStructure: SubsetsSignature;
        /// An open subset containing the point
        fn open_neighbourhood(
            &self,
            approx_point: &Self::Set,
        ) -> <Self::OpenSubsetsStructure as SetSignature>::Set;
        /// A value which decreases with refinements of the subset containing the point.
        fn precision(
            &self,
            approx_point: &Self::Set,
        ) -> <Self::Precision as SetSignature>::Set;
        /// Refine approx_point so that the value of precision() decreases
        fn refine(&self, approx_point: &mut Self::Set);
        /// Refine approx_point until its precision is at most the given value.
        fn refine_to(
            &self,
            approx_point: &mut Self::Set,
            precision: &<Self::Precision as SetSignature>::Set,
        );
    }
}
pub mod combinatorics {
    //! Contains combinatorial counting and enumeration algorithms.
    mod number_compositions {
        pub fn compositions_sized(
            n: usize,
            x: usize,
        ) -> impl Iterator<Item = Vec<usize>> {
            let c: Box<dyn Iterator<Item = Vec<usize>>> = if n == 0 && x == 0 {
                Box::new(
                    <[_]>::into_vec(::alloc::boxed::box_new([::alloc::vec::Vec::new()]))
                        .into_iter(),
                )
            } else if n == 0 || x == 0 {
                Box::new(::alloc::vec::Vec::new().into_iter())
            } else {
                Box::new(
                    super::subsets(n - 1, x - 1)
                        .map(move |mut s| {
                            s.push(n - 1);
                            let mut c = ::alloc::vec::Vec::new();
                            c.push(s[0] + 1);
                            for i in 1..x {
                                c.push(s[i] - s[i - 1]);
                            }
                            c
                        }),
                )
            };
            c
        }
        pub fn compositions_sized_zero(
            n: usize,
            x: usize,
        ) -> impl Iterator<Item = Vec<usize>> {
            compositions_sized(n + x, x).map(|c| c.into_iter().map(|i| i - 1).collect())
        }
        pub fn compositions(n: usize) -> impl Iterator<Item = Vec<usize>> {
            (0..=n).flat_map(move |x| compositions_sized(n, x))
        }
    }
    mod number_partitions {
        use itertools::Itertools;
        use std::collections::{BTreeMap, HashSet};
        struct NumPartitionIterator<P: Fn(usize) -> bool + Clone> {
            n: usize,
            x: usize,
            predicate: P,
            first: usize,
            min: usize,
            rest: Option<Box<NumPartitionIterator<P>>>,
        }
        impl<P: Fn(usize) -> bool + Clone> NumPartitionIterator<P> {
            pub fn new(n: usize, x: usize, predicate: P) -> Self {
                Self {
                    n,
                    x,
                    predicate,
                    first: 1,
                    min: 1,
                    rest: None,
                }
            }
        }
        impl<P: Fn(usize) -> bool + Clone> Iterator for NumPartitionIterator<P> {
            type Item = Vec<usize>;
            fn next(&mut self) -> Option<Self::Item> {
                if self.n < self.x {
                    None
                } else if self.n == 0 && self.x == 0 {
                    if self.first == 1 {
                        self.first += 1;
                        Some(::alloc::vec::Vec::new())
                    } else {
                        None
                    }
                } else if self.n == 0 || self.x == 0 {
                    None
                } else if self.x == 1 {
                    if self.first <= self.n {
                        self.first = self.n + 1;
                        if (self.predicate)(self.n) {
                            Some(<[_]>::into_vec(::alloc::boxed::box_new([self.n])))
                        } else {
                            None
                        }
                    } else {
                        None
                    }
                } else {
                    if self.first > self.n {
                        return None;
                    }
                    if self.rest.is_none() {
                        self.rest = Some(
                            Box::new(NumPartitionIterator {
                                n: self.n - self.first,
                                x: self.x - 1,
                                predicate: self.predicate.clone(),
                                min: self.min + 1,
                                first: self.first,
                                rest: None,
                            }),
                        );
                    }
                    if !(self.predicate)(self.first) {
                        self.first += 1;
                        self.rest = None;
                        return self.next();
                    }
                    if let Some(rest_part) = self
                        .rest
                        .as_mut()
                        .unwrap()
                        .as_mut()
                        .next()
                        .as_mut()
                    {
                        let mut part = <[_]>::into_vec(
                            ::alloc::boxed::box_new([self.first]),
                        );
                        part.append(rest_part);
                        Some(part)
                    } else {
                        self.first += 1;
                        self.rest = None;
                        self.next()
                    }
                }
            }
        }
        /// Returns all partitions of n into exactly x parts such that each part in the partition satisfies the predicate and all parts are non-zero.
        /// ```
        /// use algebraeon_sets::combinatorics::num_partitions_sized_predicated;
        /// assert_eq!(num_partitions_sized_predicated(6, 2, |k| k % 2 == 1).collect::<Vec<_>>(), vec![
        ///     vec![1, 5],
        ///     vec![3, 3],
        /// ]);
        /// assert_eq!(num_partitions_sized_predicated(6, 3, |k| k <= 3).collect::<Vec<_>>(), vec![
        ///     vec![1, 2, 3],
        ///     vec![2, 2, 2],
        /// ]);
        /// ```
        pub fn num_partitions_sized_predicated<P: Fn(usize) -> bool + Clone>(
            n: usize,
            x: usize,
            predicate: P,
        ) -> impl Iterator<Item = Vec<usize>> {
            NumPartitionIterator::new(n, x, predicate)
        }
        /// Returns all partitions of n into exactly x parts such that each part in the partition satisfies the predicate where parts may be zero.
        /// ```
        /// use algebraeon_sets::combinatorics::num_partitions_sized_zero_predicated;
        /// assert_eq!(num_partitions_sized_zero_predicated(6, 2, |k| k % 2 == 1).collect::<Vec<_>>(), vec![
        ///     vec![1, 5],
        ///     vec![3, 3],
        /// ]);
        /// assert_eq!(num_partitions_sized_zero_predicated(6, 3, |k| k <= 3).collect::<Vec<_>>(), vec![
        ///     vec![0, 3, 3],
        ///     vec![1, 2, 3],
        ///     vec![2, 2, 2],
        /// ]);
        /// ```
        pub fn num_partitions_sized_zero_predicated<P: Fn(usize) -> bool + Copy>(
            n: usize,
            x: usize,
            predicate: P,
        ) -> impl Iterator<Item = Vec<usize>> {
            NumPartitionIterator::new(n + x, x, move |k| predicate(k - 1))
                .map(|part| part.into_iter().map(|k| k - 1).collect::<Vec<usize>>())
        }
        /// Returns all partitions of n such that each part in the partition satisfies the predicate.
        /// ```
        /// use algebraeon_sets::combinatorics::num_partitions_predicated;
        /// assert_eq!(num_partitions_predicated(6, |k| k % 2 == 1).collect::<Vec<_>>(), vec![
        ///     vec![1, 5],
        ///     vec![3, 3],
        ///     vec![1, 1, 1, 3],
        ///     vec![1, 1, 1, 1, 1, 1],
        /// ]);
        /// assert_eq!(num_partitions_predicated(4, |k| k <= 2).collect::<Vec<_>>(), vec![
        ///     vec![2, 2],
        ///     vec![1, 1, 2],
        ///     vec![1, 1, 1, 1],
        /// ]);
        /// ```
        pub fn num_partitions_predicated<P: Fn(usize) -> bool + Clone>(
            n: usize,
            predicate: P,
        ) -> impl Iterator<Item = Vec<usize>> {
            (1..=n)
                .flat_map(move |x| num_partitions_sized_predicated::<
                    P,
                >(n, x, predicate.clone()))
        }
        /// Returns all partitions of n into exactly x parts where all parts are non-zero.
        /// ```
        /// use algebraeon_sets::combinatorics::num_partitions_sized;
        /// assert_eq!(num_partitions_sized(8, 3).collect::<Vec<_>>(), vec![
        ///     vec![1, 1, 6],
        ///     vec![1, 2, 5],
        ///     vec![1, 3, 4],
        ///     vec![2, 2, 4],
        ///     vec![2, 3, 3],
        /// ]);
        /// ```
        pub fn num_partitions_sized(
            n: usize,
            x: usize,
        ) -> impl Iterator<Item = Vec<usize>> {
            num_partitions_sized_predicated(n, x, |_| true)
        }
        /// Returns all partitions of n into exactly x parts where parts may be zero.
        /// ```
        /// use algebraeon_sets::combinatorics::num_partitions_sized_zero;
        /// assert_eq!(num_partitions_sized_zero(4, 3).collect::<Vec<_>>(), vec![
        ///     vec![0, 0, 4],
        ///     vec![0, 1, 3],
        ///     vec![0, 2, 2],
        ///     vec![1, 1, 2],
        /// ]);
        /// ```
        pub fn num_partitions_sized_zero(
            n: usize,
            x: usize,
        ) -> impl Iterator<Item = Vec<usize>> {
            num_partitions_sized_zero_predicated(n, x, |_| true)
        }
        /// Returns all partitions of n.
        /// ```
        /// use algebraeon_sets::combinatorics::num_partitions;
        /// assert_eq!(num_partitions(4).collect::<Vec<_>>(), vec![
        ///     vec![4],
        ///     vec![1, 3],
        ///     vec![2, 2],
        ///     vec![1, 1, 2],
        ///     vec![1, 1, 1, 1]
        /// ]);
        /// ```
        pub fn num_partitions(n: usize) -> impl Iterator<Item = Vec<usize>> {
            num_partitions_predicated(n, |_| true)
        }
        /// Return all partitions of n where parts come from a pool of allowed part sizes.
        /// Two different parts in the part size pool count as different parts.
        /// Parts may appear zero times or multiple times.
        /// ```
        /// use algebraeon_sets::combinatorics::num_partitions_part_pool;
        /// assert_eq!(num_partitions_part_pool(8, vec![2, 3, 3]).collect::<Vec<_>>(), vec![
        ///     vec![0, 1, 1],
        ///     vec![0, 1, 2],
        ///     vec![0, 2, 1],
        ///     vec![0, 2, 2],
        ///     vec![0, 0, 0, 0]
        /// ]);
        /// ```
        #[allow(clippy::needless_pass_by_value)]
        pub fn num_partitions_part_pool(
            n: usize,
            part_pool: Vec<usize>,
        ) -> impl Iterator<Item = Vec<usize>> {
            let mut allowed_sizes: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
            for (idx, size) in part_pool.iter().enumerate() {
                #[allow(clippy::unwrap_or_default)]
                allowed_sizes.entry(*size).or_insert(::alloc::vec::Vec::new()).push(idx);
            }
            let allowed_sizes_set = allowed_sizes
                .keys()
                .copied()
                .collect::<HashSet<_>>();
            num_partitions_predicated(n, move |k| allowed_sizes_set.contains(&k))
                .flat_map(move |partition| {
                    let mut counts: BTreeMap<usize, usize> = allowed_sizes
                        .keys()
                        .map(|x| (*x, 0))
                        .collect();
                    for p in partition {
                        *counts.get_mut(&p).unwrap() += 1;
                    }
                    counts
                        .into_iter()
                        .map(|(size, count)| {
                            (0..count)
                                .map(|_| allowed_sizes.get(&size).unwrap().iter().copied())
                                .multi_cartesian_product()
                                .collect::<Vec<_>>()
                        })
                        .multi_cartesian_product()
                        .map(|z| z.into_iter().flatten().collect::<Vec<_>>())
                        .collect::<Vec<_>>()
                })
        }
    }
    mod set_partitions {
        use indexmap::IndexMap;
        use itertools::Itertools;
        use std::collections::HashSet;
        pub struct Partition {
            partition: Vec<HashSet<usize>>,
            lookup: Vec<usize>,
        }
        #[automatically_derived]
        impl ::core::clone::Clone for Partition {
            #[inline]
            fn clone(&self) -> Partition {
                Partition {
                    partition: ::core::clone::Clone::clone(&self.partition),
                    lookup: ::core::clone::Clone::clone(&self.lookup),
                }
            }
        }
        #[automatically_derived]
        impl ::core::fmt::Debug for Partition {
            #[inline]
            fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
                ::core::fmt::Formatter::debug_struct_field2_finish(
                    f,
                    "Partition",
                    "partition",
                    &self.partition,
                    "lookup",
                    &&self.lookup,
                )
            }
        }
        impl Partition {
            fn check_state(&self) -> Result<(), &'static str> {
                use std::collections::HashMap;
                let mut present = HashMap::new();
                let n = self.lookup.len();
                for (idx, part) in self.partition.iter().enumerate() {
                    if part.is_empty() {
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
            pub fn new_unchecked(
                partition: Vec<HashSet<usize>>,
                lookup: Vec<usize>,
            ) -> Self {
                let partition = Self { partition, lookup };
                partition.check_state().unwrap();
                partition
            }
            pub fn new_from_function<T: Clone + Eq + std::hash::Hash>(
                n: usize,
                f: impl Fn(usize) -> T,
            ) -> (Self, Vec<T>) {
                let mut t_lookup = ::alloc::vec::Vec::new();
                for x in 0..n {
                    t_lookup.push(f(x));
                }
                let mut t_partition: IndexMap<_, Vec<usize>> = IndexMap::new();
                #[allow(clippy::needless_range_loop)]
                for x in 0..n {
                    let t = &t_lookup[x];
                    if t_partition.contains_key(&t) {
                        t_partition.get_mut(&t).unwrap().push(x);
                    } else {
                        t_partition
                            .insert(t, <[_]>::into_vec(::alloc::boxed::box_new([x])));
                    }
                }
                let lookup = (0..n)
                    .map(|x| t_partition.get_index_of(&t_lookup[x]).unwrap())
                    .collect();
                let partition = t_partition
                    .iter()
                    .map(|(_t, part)| part.iter().copied().collect())
                    .collect();
                let partition = Partition::new_unchecked(partition, lookup);
                partition.check_state().unwrap();
                (
                    partition,
                    t_partition.into_iter().map(|(t, _part)| t.clone()).collect(),
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
        pub struct Element {
            x: usize,
            cum_x: usize,
            pivot: bool,
        }
        #[automatically_derived]
        impl ::core::fmt::Debug for Element {
            #[inline]
            fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
                ::core::fmt::Formatter::debug_struct_field3_finish(
                    f,
                    "Element",
                    "x",
                    &self.x,
                    "cum_x",
                    &self.cum_x,
                    "pivot",
                    &&self.pivot,
                )
            }
        }
        #[automatically_derived]
        impl ::core::clone::Clone for Element {
            #[inline]
            fn clone(&self) -> Element {
                Element {
                    x: ::core::clone::Clone::clone(&self.x),
                    cum_x: ::core::clone::Clone::clone(&self.cum_x),
                    pivot: ::core::clone::Clone::clone(&self.pivot),
                }
            }
        }
        pub struct LexicographicPartitionsNumPartsInRange {
            n: usize,
            min_x: usize,
            max_x: usize,
            elements: Vec<Element>,
            finished: bool,
        }
        #[automatically_derived]
        impl ::core::fmt::Debug for LexicographicPartitionsNumPartsInRange {
            #[inline]
            fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
                ::core::fmt::Formatter::debug_struct_field5_finish(
                    f,
                    "LexicographicPartitionsNumPartsInRange",
                    "n",
                    &self.n,
                    "min_x",
                    &self.min_x,
                    "max_x",
                    &self.max_x,
                    "elements",
                    &self.elements,
                    "finished",
                    &&self.finished,
                )
            }
        }
        #[automatically_derived]
        impl ::core::clone::Clone for LexicographicPartitionsNumPartsInRange {
            #[inline]
            fn clone(&self) -> LexicographicPartitionsNumPartsInRange {
                LexicographicPartitionsNumPartsInRange {
                    n: ::core::clone::Clone::clone(&self.n),
                    min_x: ::core::clone::Clone::clone(&self.min_x),
                    max_x: ::core::clone::Clone::clone(&self.max_x),
                    elements: ::core::clone::Clone::clone(&self.elements),
                    finished: ::core::clone::Clone::clone(&self.finished),
                }
            }
        }
        impl LexicographicPartitionsNumPartsInRange {
            #[allow(clippy::unnecessary_wraps)]
            fn check(&self) -> Result<(), ()> {
                if !self.finished {
                    match (&self.elements.len(), &self.n) {
                        (left_val, right_val) => {
                            if !(*left_val == *right_val) {
                                let kind = ::core::panicking::AssertKind::Eq;
                                ::core::panicking::assert_failed(
                                    kind,
                                    &*left_val,
                                    &*right_val,
                                    ::core::option::Option::None,
                                );
                            }
                        }
                    };
                    match (&self.elements[0].x, &0) {
                        (left_val, right_val) => {
                            if !(*left_val == *right_val) {
                                let kind = ::core::panicking::AssertKind::Eq;
                                ::core::panicking::assert_failed(
                                    kind,
                                    &*left_val,
                                    &*right_val,
                                    ::core::option::Option::None,
                                );
                            }
                        }
                    };
                    match (&self.elements[0].cum_x, &0) {
                        (left_val, right_val) => {
                            if !(*left_val == *right_val) {
                                let kind = ::core::panicking::AssertKind::Eq;
                                ::core::panicking::assert_failed(
                                    kind,
                                    &*left_val,
                                    &*right_val,
                                    ::core::option::Option::None,
                                );
                            }
                        }
                    };
                    if !self.elements[0].pivot {
                        ::core::panicking::panic(
                            "assertion failed: self.elements[0].pivot",
                        )
                    }
                    let mut cum_max = 0;
                    for i in 1..self.n {
                        if self.elements[i].x <= cum_max {
                            match (&self.elements[i].cum_x, &cum_max) {
                                (left_val, right_val) => {
                                    if !(*left_val == *right_val) {
                                        let kind = ::core::panicking::AssertKind::Eq;
                                        ::core::panicking::assert_failed(
                                            kind,
                                            &*left_val,
                                            &*right_val,
                                            ::core::option::Option::None,
                                        );
                                    }
                                }
                            };
                            if !!self.elements[i].pivot {
                                ::core::panicking::panic(
                                    "assertion failed: !self.elements[i].pivot",
                                )
                            }
                        } else if self.elements[i].x == cum_max + 1 {
                            cum_max += 1;
                            match (&self.elements[i].cum_x, &cum_max) {
                                (left_val, right_val) => {
                                    if !(*left_val == *right_val) {
                                        let kind = ::core::panicking::AssertKind::Eq;
                                        ::core::panicking::assert_failed(
                                            kind,
                                            &*left_val,
                                            &*right_val,
                                            ::core::option::Option::None,
                                        );
                                    }
                                }
                            };
                            if !self.elements[i].pivot {
                                ::core::panicking::panic(
                                    "assertion failed: self.elements[i].pivot",
                                )
                            }
                        } else {
                            ::core::panicking::panic("explicit panic");
                        }
                    }
                    cum_max += 1;
                    if !(self.min_x <= cum_max) {
                        ::core::panicking::panic(
                            "assertion failed: self.min_x <= cum_max",
                        )
                    }
                    if !(cum_max <= self.max_x) {
                        ::core::panicking::panic(
                            "assertion failed: cum_max <= self.max_x",
                        )
                    }
                }
                Ok(())
            }
            pub fn new(n: usize, min_x: usize, max_x: usize) -> Self {
                let mut elements = ::alloc::vec::Vec::new();
                for i in 0..n {
                    elements
                        .push(Element {
                            x: 0,
                            cum_x: 0,
                            pivot: i == 0,
                        });
                }
                let mut s = Self {
                    n,
                    min_x,
                    max_x,
                    elements,
                    finished: false,
                };
                if (n == 0 && min_x > 0) || (n > 0 && max_x == 0) || (n < min_x)
                    || (min_x > max_x)
                {
                    s.finished = true;
                }
                if n > 0 {
                    s.reset_tail(0);
                }
                s
            }
            fn reset_tail(&mut self, j: usize) {
                let cum_max_j = self.elements[j].cum_x;
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
                self.check().unwrap();
            }
        }
        impl Iterator for LexicographicPartitionsNumPartsInRange {
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
                                    #[allow(clippy::comparison_chain)]
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
        pub fn set_partitions_eq(
            n: usize,
            x: usize,
        ) -> impl Iterator<Item = Vec<usize>> {
            LexicographicPartitionsNumPartsInRange::new(n, x, x)
        }
        pub fn set_partitions_le(
            n: usize,
            x: usize,
        ) -> impl Iterator<Item = Vec<usize>> {
            LexicographicPartitionsNumPartsInRange::new(n, 0, x)
        }
        pub fn set_partitions_ge(
            n: usize,
            x: usize,
        ) -> impl Iterator<Item = Vec<usize>> {
            LexicographicPartitionsNumPartsInRange::new(n, x, n)
        }
        pub fn set_partitions_range(
            n: usize,
            min_x: usize,
            max_x: usize,
        ) -> impl Iterator<Item = Vec<usize>> {
            LexicographicPartitionsNumPartsInRange::new(n, min_x, max_x)
        }
        pub fn set_compositions_eq(
            n: usize,
            x: usize,
        ) -> impl Iterator<Item = Vec<usize>> {
            (0..x)
                .permutations(x)
                .flat_map(move |perm| {
                    set_partitions_eq(n, x)
                        .map(move |partition| {
                            partition.into_iter().map(|i| perm[i]).collect()
                        })
                })
        }
    }
    mod subsets {
        /// Provides an iterator over the k-element subsets of {0, 1, ..., n-1} in lexicographic order such that some elements of {0, 1, ..., n-1} can be excluded from future subsets at any point during the iteration.
        pub struct LexicographicSubsetsWithRemovals {
            n: usize,
            all_items_idx: Vec<Option<usize>>,
            remaining_items: Vec<usize>,
            subset: Vec<usize>,
            finished: bool,
        }
        #[automatically_derived]
        impl ::core::fmt::Debug for LexicographicSubsetsWithRemovals {
            #[inline]
            fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
                ::core::fmt::Formatter::debug_struct_field5_finish(
                    f,
                    "LexicographicSubsetsWithRemovals",
                    "n",
                    &self.n,
                    "all_items_idx",
                    &self.all_items_idx,
                    "remaining_items",
                    &self.remaining_items,
                    "subset",
                    &self.subset,
                    "finished",
                    &&self.finished,
                )
            }
        }
        impl LexicographicSubsetsWithRemovals {
            /// Constructor
            pub fn new(n: usize, k: usize) -> Self {
                Self {
                    n,
                    all_items_idx: (0..n).map(Some).collect(),
                    remaining_items: (0..n).collect(),
                    subset: (0..k).collect(),
                    finished: k > n,
                }
            }
            /// Exclude the provided element from appearing in any further subsets.
            pub fn exclude(&mut self, a: usize) {
                if !(a < self.all_items_idx.len()) {
                    ::core::panicking::panic(
                        "assertion failed: a < self.all_items_idx.len()",
                    )
                }
                if self.all_items_idx[a].is_none() {
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
                        if let Some(idx) = self.all_items_idx[i].as_mut() {
                            *idx -= 1;
                        }
                    }
                    for x in &mut self.subset {
                        match (*x).cmp(&a_idx) {
                            std::cmp::Ordering::Less => {}
                            std::cmp::Ordering::Equal => {
                                ::core::panicking::panic(
                                    "internal error: entered unreachable code",
                                )
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
        impl Iterator for LexicographicSubsetsWithRemovals {
            type Item = Vec<usize>;
            fn next(&mut self) -> Option<Self::Item> {
                let next = if self.finished {
                    return None;
                } else {
                    Some(self.subset.iter().map(|i| self.remaining_items[*i]).collect())
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
        /// use algebraeon_sets::combinatorics::subsets;
        /// assert_eq!(subsets(5, 3).collect::<Vec<_>>(), vec![
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
        pub fn subsets(n: usize, k: usize) -> impl Iterator<Item = Vec<usize>> {
            LexicographicSubsetsWithRemovals::new(n, k)
        }
        /// Returns all subsets of {0, 1, ..., n-1}.
        /// ```
        /// use algebraeon_sets::combinatorics::all_subsets;
        /// assert_eq!(all_subsets(3).collect::<Vec<_>>(), vec![
        ///     vec![],
        ///     vec![0],
        ///     vec![1],
        ///     vec![0, 1],
        ///     vec![2],
        ///     vec![0, 2],
        ///     vec![1, 2],
        ///     vec![0, 1, 2],
        /// ]);
        /// ```
        pub fn all_subsets(n: usize) -> impl Iterator<Item = Vec<usize>> {
            (0usize..(1 << n))
                .map(move |i| (0..n).filter(|j| i & (1 << j) != 0).collect())
        }
        /// Returns all size k subsets of items.
        /// ```
        /// use algebraeon_sets::combinatorics::subsets_of_vec;
        /// assert_eq!(subsets_of_vec(vec!["a", "b", "c"], 2).collect::<Vec<_>>(), vec![
        ///     vec!["a", "b"],
        ///     vec!["a", "c"],
        ///     vec!["b", "c"],
        /// ]);
        /// ```
        pub fn subsets_of_vec<'a, T: 'a + Clone>(
            items: Vec<T>,
            k: usize,
        ) -> impl 'a + Iterator<Item = Vec<T>> {
            subsets(items.len(), k)
                .map(move |subset| {
                    subset.into_iter().map(|idx| items[idx].clone()).collect()
                })
        }
    }
    mod twelvefold_way {}
    pub use number_compositions::compositions;
    pub use number_compositions::compositions_sized;
    pub use number_compositions::compositions_sized_zero;
    pub use number_partitions::num_partitions;
    pub use number_partitions::num_partitions_part_pool;
    pub use number_partitions::num_partitions_predicated;
    pub use number_partitions::num_partitions_sized;
    pub use number_partitions::num_partitions_sized_predicated;
    pub use number_partitions::num_partitions_sized_zero;
    pub use number_partitions::num_partitions_sized_zero_predicated;
    pub use set_partitions::Partition;
    pub use set_partitions::set_compositions_eq;
    pub use set_partitions::set_partitions_eq;
    pub use set_partitions::set_partitions_ge;
    pub use set_partitions::set_partitions_le;
    pub use set_partitions::set_partitions_range;
    pub use subsets::LexicographicSubsetsWithRemovals;
    pub use subsets::all_subsets;
    pub use subsets::subsets;
    pub use subsets::subsets_of_vec;
}
pub mod structure {
    //! Abstractions over sets with certain structure.
    //!
    //! The structure framework used by `algebraeon_rings` is established here.
    mod empty_set {
        use crate::structure::orderings::PartialOrdSignature;
        use super::{EqSignature, OrdSignature, SetSignature, Signature};
        use std::fmt::Debug;
        use std::marker::PhantomData;
        pub struct EmptySetStructure<Set> {
            _set: PhantomData<Set>,
        }
        impl<Set> Clone for EmptySetStructure<Set> {
            fn clone(&self) -> Self {
                Self { _set: PhantomData }
            }
        }
        impl<Set> Debug for EmptySetStructure<Set> {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                f.debug_struct("EmptySetStructure").finish()
            }
        }
        impl<Set> PartialEq for EmptySetStructure<Set> {
            fn eq(&self, _: &Self) -> bool {
                true
            }
        }
        impl<Set> Eq for EmptySetStructure<Set> {}
        impl<Set> Default for EmptySetStructure<Set> {
            fn default() -> Self {
                Self { _set: PhantomData }
            }
        }
        impl<Set: Send + Sync> Signature for EmptySetStructure<Set> {}
        impl<Set: Debug + Clone + Send + Sync> SetSignature for EmptySetStructure<Set> {
            type Set = Set;
            fn is_element(&self, _: &Self::Set) -> Result<(), String> {
                Err("Empty set has no elements".to_string())
            }
        }
        impl<Set: Debug + Clone + Send + Sync> EqSignature for EmptySetStructure<Set> {
            fn equal(&self, _: &Self::Set, _: &Self::Set) -> bool {
                {
                    ::core::panicking::panic_fmt(
                        format_args!("Empty set had no elements to compare for equality"),
                    );
                }
            }
        }
        impl<Set: Debug + Clone + Send + Sync> PartialOrdSignature
        for EmptySetStructure<Set> {
            fn partial_cmp(
                &self,
                a: &Self::Set,
                b: &Self::Set,
            ) -> Option<std::cmp::Ordering> {
                Some(self.cmp(a, b))
            }
        }
        impl<Set: Debug + Clone + Send + Sync> OrdSignature for EmptySetStructure<Set> {
            fn cmp(&self, _: &Self::Set, _: &Self::Set) -> std::cmp::Ordering {
                {
                    ::core::panicking::panic_fmt(
                        format_args!("Empty set had no elements to compare for ordering"),
                    );
                }
            }
        }
    }
    mod finite_set {
        use super::{EqSignature, OrdSignature, SetSignature, Signature};
        use crate::structure::{
            CountableSetSignature, FiniteSetSignature, orderings::PartialOrdSignature,
        };
        use std::fmt::Debug;
        pub struct EnumeratedFiniteSetStructure {
            n: usize,
        }
        #[automatically_derived]
        impl ::core::fmt::Debug for EnumeratedFiniteSetStructure {
            #[inline]
            fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
                ::core::fmt::Formatter::debug_struct_field1_finish(
                    f,
                    "EnumeratedFiniteSetStructure",
                    "n",
                    &&self.n,
                )
            }
        }
        #[automatically_derived]
        impl ::core::clone::Clone for EnumeratedFiniteSetStructure {
            #[inline]
            fn clone(&self) -> EnumeratedFiniteSetStructure {
                EnumeratedFiniteSetStructure {
                    n: ::core::clone::Clone::clone(&self.n),
                }
            }
        }
        #[automatically_derived]
        impl ::core::marker::StructuralPartialEq for EnumeratedFiniteSetStructure {}
        #[automatically_derived]
        impl ::core::cmp::PartialEq for EnumeratedFiniteSetStructure {
            #[inline]
            fn eq(&self, other: &EnumeratedFiniteSetStructure) -> bool {
                self.n == other.n
            }
        }
        #[automatically_derived]
        impl ::core::cmp::Eq for EnumeratedFiniteSetStructure {
            #[inline]
            #[doc(hidden)]
            #[coverage(off)]
            fn assert_receiver_is_total_eq(&self) -> () {
                let _: ::core::cmp::AssertParamIsEq<usize>;
            }
        }
        impl EnumeratedFiniteSetStructure {
            pub fn new(n: usize) -> Self {
                Self { n }
            }
        }
        impl Signature for EnumeratedFiniteSetStructure {}
        impl SetSignature for EnumeratedFiniteSetStructure {
            type Set = usize;
            fn is_element(&self, x: &Self::Set) -> Result<(), String> {
                if x >= &self.n {
                    return Err("Too big to be an element".to_string());
                }
                Ok(())
            }
        }
        impl EqSignature for EnumeratedFiniteSetStructure {
            fn equal(&self, x: &Self::Set, y: &Self::Set) -> bool {
                x == y
            }
        }
        impl PartialOrdSignature for EnumeratedFiniteSetStructure {
            fn partial_cmp(
                &self,
                x: &Self::Set,
                y: &Self::Set,
            ) -> Option<std::cmp::Ordering> {
                Some(self.cmp(x, y))
            }
        }
        impl OrdSignature for EnumeratedFiniteSetStructure {
            fn cmp(&self, x: &Self::Set, y: &Self::Set) -> std::cmp::Ordering {
                x.cmp(y)
            }
        }
        impl CountableSetSignature for EnumeratedFiniteSetStructure {
            fn generate_all_elements(&self) -> impl Iterator<Item = Self::Set> + Clone {
                0..self.n
            }
        }
        impl FiniteSetSignature for EnumeratedFiniteSetStructure {
            fn size(&self) -> usize {
                self.n
            }
        }
    }
    mod morphism {
        use super::{
            CountableSetSignature, EqSignature, FiniteSetSignature, SetSignature,
            Signature,
        };
        use itertools::Itertools;
        use std::borrow::Borrow;
        use std::fmt::Debug;
        use std::marker::PhantomData;
        pub trait Morphism<
            Domain: Signature,
            Range: Signature,
        >: Debug + Clone + Send + Sync {
            fn domain(&self) -> &Domain;
            fn range(&self) -> &Range;
        }
        /// A morphism from an object to itself
        pub trait Endomorphism<X: Signature>: Morphism<X, X> {}
        impl<X: Signature, T: Morphism<X, X>> Endomorphism<X> for T {}
        pub trait Function<
            Domain: SetSignature,
            Range: SetSignature,
        >: Morphism<Domain, Range> {
            fn image(&self, x: &Domain::Set) -> Range::Set;
        }
        /// A function from a set into itself
        pub trait Endofunction<X: SetSignature + EqSignature>: Function<X, X> {
            /// check if an element is fixed
            fn is_fixed_point(&self, x: X::Set) -> bool {
                self.domain().equal(&self.image(&x), &x)
            }
        }
        pub trait InjectiveFunction<
            Domain: SetSignature,
            Range: SetSignature,
        >: Function<Domain, Range> {
            fn try_preimage(&self, y: &Range::Set) -> Option<Domain::Set>;
        }
        pub trait BijectiveFunction<
            Domain: SetSignature,
            Range: SetSignature,
        >: InjectiveFunction<Domain, Range> {
            fn preimage(&self, y: &Range::Set) -> Domain::Set {
                self.try_preimage(y).unwrap()
            }
        }
        /// A permutation is a bijective function from a set to itself
        pub trait Permutation<X: SetSignature>: BijectiveFunction<X, X> {}
        impl<X: SetSignature, T: BijectiveFunction<X, X>> Permutation<X> for T {}
        /// The identity morphism X -> X
        pub struct IdentityMorphism<X: Signature> {
            x: X,
        }
        #[automatically_derived]
        impl<X: ::core::fmt::Debug + Signature> ::core::fmt::Debug
        for IdentityMorphism<X> {
            #[inline]
            fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
                ::core::fmt::Formatter::debug_struct_field1_finish(
                    f,
                    "IdentityMorphism",
                    "x",
                    &&self.x,
                )
            }
        }
        #[automatically_derived]
        impl<X: ::core::clone::Clone + Signature> ::core::clone::Clone
        for IdentityMorphism<X> {
            #[inline]
            fn clone(&self) -> IdentityMorphism<X> {
                IdentityMorphism {
                    x: ::core::clone::Clone::clone(&self.x),
                }
            }
        }
        #[automatically_derived]
        impl<X: Signature> ::core::marker::StructuralPartialEq for IdentityMorphism<X> {}
        #[automatically_derived]
        impl<X: ::core::cmp::PartialEq + Signature> ::core::cmp::PartialEq
        for IdentityMorphism<X> {
            #[inline]
            fn eq(&self, other: &IdentityMorphism<X>) -> bool {
                self.x == other.x
            }
        }
        #[automatically_derived]
        impl<X: ::core::cmp::Eq + Signature> ::core::cmp::Eq for IdentityMorphism<X> {
            #[inline]
            #[doc(hidden)]
            #[coverage(off)]
            fn assert_receiver_is_total_eq(&self) -> () {
                let _: ::core::cmp::AssertParamIsEq<X>;
            }
        }
        impl<X: Signature> IdentityMorphism<X> {
            pub fn new(x: X) -> Self {
                Self { x }
            }
        }
        impl<X: Signature> Morphism<X, X> for IdentityMorphism<X> {
            fn domain(&self) -> &X {
                &self.x
            }
            fn range(&self) -> &X {
                &self.x
            }
        }
        impl<X: SetSignature> Function<X, X> for IdentityMorphism<X> {
            fn image(&self, x: &X::Set) -> X::Set {
                x.clone()
            }
        }
        impl<X: SetSignature> InjectiveFunction<X, X> for IdentityMorphism<X> {
            fn try_preimage(&self, x: &X::Set) -> Option<X::Set> {
                Some(x.clone())
            }
        }
        impl<X: SetSignature> BijectiveFunction<X, X> for IdentityMorphism<X> {}
        /// The composition A -> B -> C of two morphisms A -> B and B -> C
        pub struct CompositionMorphism<
            A: Signature,
            B: Signature,
            C: Signature,
            AB: Morphism<A, B>,
            BC: Morphism<B, C>,
        > {
            a: PhantomData<A>,
            b: PhantomData<B>,
            c: PhantomData<C>,
            a_to_b: AB,
            b_to_c: BC,
        }
        #[automatically_derived]
        impl<
            A: ::core::fmt::Debug + Signature,
            B: ::core::fmt::Debug + Signature,
            C: ::core::fmt::Debug + Signature,
            AB: ::core::fmt::Debug + Morphism<A, B>,
            BC: ::core::fmt::Debug + Morphism<B, C>,
        > ::core::fmt::Debug for CompositionMorphism<A, B, C, AB, BC> {
            #[inline]
            fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
                ::core::fmt::Formatter::debug_struct_field5_finish(
                    f,
                    "CompositionMorphism",
                    "a",
                    &self.a,
                    "b",
                    &self.b,
                    "c",
                    &self.c,
                    "a_to_b",
                    &self.a_to_b,
                    "b_to_c",
                    &&self.b_to_c,
                )
            }
        }
        #[automatically_derived]
        impl<
            A: ::core::clone::Clone + Signature,
            B: ::core::clone::Clone + Signature,
            C: ::core::clone::Clone + Signature,
            AB: ::core::clone::Clone + Morphism<A, B>,
            BC: ::core::clone::Clone + Morphism<B, C>,
        > ::core::clone::Clone for CompositionMorphism<A, B, C, AB, BC> {
            #[inline]
            fn clone(&self) -> CompositionMorphism<A, B, C, AB, BC> {
                CompositionMorphism {
                    a: ::core::clone::Clone::clone(&self.a),
                    b: ::core::clone::Clone::clone(&self.b),
                    c: ::core::clone::Clone::clone(&self.c),
                    a_to_b: ::core::clone::Clone::clone(&self.a_to_b),
                    b_to_c: ::core::clone::Clone::clone(&self.b_to_c),
                }
            }
        }
        #[automatically_derived]
        impl<
            A: Signature,
            B: Signature,
            C: Signature,
            AB: Morphism<A, B>,
            BC: Morphism<B, C>,
        > ::core::marker::StructuralPartialEq for CompositionMorphism<A, B, C, AB, BC> {}
        #[automatically_derived]
        impl<
            A: ::core::cmp::PartialEq + Signature,
            B: ::core::cmp::PartialEq + Signature,
            C: ::core::cmp::PartialEq + Signature,
            AB: ::core::cmp::PartialEq + Morphism<A, B>,
            BC: ::core::cmp::PartialEq + Morphism<B, C>,
        > ::core::cmp::PartialEq for CompositionMorphism<A, B, C, AB, BC> {
            #[inline]
            fn eq(&self, other: &CompositionMorphism<A, B, C, AB, BC>) -> bool {
                self.a == other.a && self.b == other.b && self.c == other.c
                    && self.a_to_b == other.a_to_b && self.b_to_c == other.b_to_c
            }
        }
        #[automatically_derived]
        impl<
            A: ::core::cmp::Eq + Signature,
            B: ::core::cmp::Eq + Signature,
            C: ::core::cmp::Eq + Signature,
            AB: ::core::cmp::Eq + Morphism<A, B>,
            BC: ::core::cmp::Eq + Morphism<B, C>,
        > ::core::cmp::Eq for CompositionMorphism<A, B, C, AB, BC> {
            #[inline]
            #[doc(hidden)]
            #[coverage(off)]
            fn assert_receiver_is_total_eq(&self) -> () {
                let _: ::core::cmp::AssertParamIsEq<PhantomData<A>>;
                let _: ::core::cmp::AssertParamIsEq<PhantomData<B>>;
                let _: ::core::cmp::AssertParamIsEq<PhantomData<C>>;
                let _: ::core::cmp::AssertParamIsEq<AB>;
                let _: ::core::cmp::AssertParamIsEq<BC>;
            }
        }
        impl<
            A: Signature,
            B: Signature,
            C: Signature,
            AB: Morphism<A, B>,
            BC: Morphism<B, C>,
        > CompositionMorphism<A, B, C, AB, BC> {
            pub fn new(a_to_b: AB, b_to_c: BC) -> Self {
                match (&a_to_b.range(), &b_to_c.domain()) {
                    (left_val, right_val) => {
                        if !(*left_val == *right_val) {
                            let kind = ::core::panicking::AssertKind::Eq;
                            ::core::panicking::assert_failed(
                                kind,
                                &*left_val,
                                &*right_val,
                                ::core::option::Option::None,
                            );
                        }
                    }
                };
                Self {
                    a: PhantomData,
                    b: PhantomData,
                    c: PhantomData,
                    a_to_b,
                    b_to_c,
                }
            }
            pub fn a(&self) -> &A {
                self.a_to_b.domain()
            }
            pub fn b(&self) -> &B {
                let b = self.a_to_b.range();
                if true {
                    match (&b, &self.b_to_c.domain()) {
                        (left_val, right_val) => {
                            if !(*left_val == *right_val) {
                                let kind = ::core::panicking::AssertKind::Eq;
                                ::core::panicking::assert_failed(
                                    kind,
                                    &*left_val,
                                    &*right_val,
                                    ::core::option::Option::None,
                                );
                            }
                        }
                    };
                }
                b
            }
            pub fn c(&self) -> &C {
                self.b_to_c.range()
            }
            pub fn a_to_b(&self) -> &AB {
                &self.a_to_b
            }
            pub fn b_to_c(&self) -> &BC {
                &self.b_to_c
            }
        }
        impl<
            A: Signature,
            B: Signature,
            C: Signature,
            AB: Morphism<A, B>,
            BC: Morphism<B, C>,
        > Morphism<A, C> for CompositionMorphism<A, B, C, AB, BC> {
            fn domain(&self) -> &A {
                self.a()
            }
            fn range(&self) -> &C {
                self.c()
            }
        }
        impl<
            A: SetSignature,
            B: SetSignature,
            C: SetSignature,
            AB: Function<A, B>,
            BC: Function<B, C>,
        > Function<A, C> for CompositionMorphism<A, B, C, AB, BC> {
            fn image(&self, x: &A::Set) -> C::Set {
                self.b_to_c.image(&self.a_to_b.image(x))
            }
        }
        impl<
            A: SetSignature,
            B: SetSignature,
            C: SetSignature,
            AB: InjectiveFunction<A, B>,
            BC: InjectiveFunction<B, C>,
        > InjectiveFunction<A, C> for CompositionMorphism<A, B, C, AB, BC> {
            fn try_preimage(&self, x: &C::Set) -> Option<A::Set> {
                self.a_to_b.try_preimage(&self.b_to_c.try_preimage(x)?)
            }
        }
        impl<
            A: SetSignature,
            B: SetSignature,
            C: SetSignature,
            AB: BijectiveFunction<A, B>,
            BC: BijectiveFunction<B, C>,
        > BijectiveFunction<A, C> for CompositionMorphism<A, B, C, AB, BC> {
            fn preimage(&self, x: &C::Set) -> A::Set {
                self.a_to_b.preimage(&self.b_to_c.preimage(x))
            }
        }
        /// Represent all functions from `domain` to `range`
        pub struct Functions<Domain: SetSignature, Range: SetSignature> {
            domain: Domain,
            range: Range,
        }
        #[automatically_derived]
        impl<
            Domain: ::core::fmt::Debug + SetSignature,
            Range: ::core::fmt::Debug + SetSignature,
        > ::core::fmt::Debug for Functions<Domain, Range> {
            #[inline]
            fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
                ::core::fmt::Formatter::debug_struct_field2_finish(
                    f,
                    "Functions",
                    "domain",
                    &self.domain,
                    "range",
                    &&self.range,
                )
            }
        }
        #[automatically_derived]
        impl<
            Domain: ::core::clone::Clone + SetSignature,
            Range: ::core::clone::Clone + SetSignature,
        > ::core::clone::Clone for Functions<Domain, Range> {
            #[inline]
            fn clone(&self) -> Functions<Domain, Range> {
                Functions {
                    domain: ::core::clone::Clone::clone(&self.domain),
                    range: ::core::clone::Clone::clone(&self.range),
                }
            }
        }
        #[automatically_derived]
        impl<
            Domain: SetSignature,
            Range: SetSignature,
        > ::core::marker::StructuralPartialEq for Functions<Domain, Range> {}
        #[automatically_derived]
        impl<
            Domain: ::core::cmp::PartialEq + SetSignature,
            Range: ::core::cmp::PartialEq + SetSignature,
        > ::core::cmp::PartialEq for Functions<Domain, Range> {
            #[inline]
            fn eq(&self, other: &Functions<Domain, Range>) -> bool {
                self.domain == other.domain && self.range == other.range
            }
        }
        #[automatically_derived]
        impl<
            Domain: ::core::cmp::Eq + SetSignature,
            Range: ::core::cmp::Eq + SetSignature,
        > ::core::cmp::Eq for Functions<Domain, Range> {
            #[inline]
            #[doc(hidden)]
            #[coverage(off)]
            fn assert_receiver_is_total_eq(&self) -> () {
                let _: ::core::cmp::AssertParamIsEq<Domain>;
                let _: ::core::cmp::AssertParamIsEq<Range>;
            }
        }
        impl<Domain: SetSignature, Range: SetSignature> Functions<Domain, Range> {
            pub fn new(domain: Domain, range: Range) -> Self {
                Self { domain, range }
            }
        }
        impl<Domain: SetSignature, Range: SetSignature> Signature
        for Functions<Domain, Range> {}
        impl<Domain: FiniteSetSignature, Range: EqSignature> SetSignature
        for Functions<Domain, Range> {
            type Set = Vec<Range::Set>;
            fn is_element(&self, x: &Self::Set) -> Result<(), String> {
                if x.len() != self.domain.size() {
                    return Err("Incorrect vector length".to_string());
                }
                for y in x {
                    self.range.is_element(y)?;
                }
                Ok(())
            }
        }
        impl<
            Domain: FiniteSetSignature,
            Range: EqSignature + FiniteSetSignature,
        > CountableSetSignature for Functions<Domain, Range> {
            fn generate_all_elements(&self) -> impl Iterator<Item = Self::Set> + Clone {
                (0..self.domain.size())
                    .map(|_| self.range.list_all_elements())
                    .multi_cartesian_product()
            }
        }
        impl<
            Domain: FiniteSetSignature,
            Range: EqSignature + FiniteSetSignature,
        > FiniteSetSignature for Functions<Domain, Range> {}
        /// The set of all endofunctions on a finite set X: functions X → X
        pub struct FiniteSetEndofunctions<X: FiniteSetSignature + EqSignature> {
            set: X,
        }
        #[automatically_derived]
        impl<X: ::core::fmt::Debug + FiniteSetSignature + EqSignature> ::core::fmt::Debug
        for FiniteSetEndofunctions<X> {
            #[inline]
            fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
                ::core::fmt::Formatter::debug_struct_field1_finish(
                    f,
                    "FiniteSetEndofunctions",
                    "set",
                    &&self.set,
                )
            }
        }
        #[automatically_derived]
        impl<
            X: ::core::clone::Clone + FiniteSetSignature + EqSignature,
        > ::core::clone::Clone for FiniteSetEndofunctions<X> {
            #[inline]
            fn clone(&self) -> FiniteSetEndofunctions<X> {
                FiniteSetEndofunctions {
                    set: ::core::clone::Clone::clone(&self.set),
                }
            }
        }
        #[automatically_derived]
        impl<X: FiniteSetSignature + EqSignature> ::core::marker::StructuralPartialEq
        for FiniteSetEndofunctions<X> {}
        #[automatically_derived]
        impl<
            X: ::core::cmp::PartialEq + FiniteSetSignature + EqSignature,
        > ::core::cmp::PartialEq for FiniteSetEndofunctions<X> {
            #[inline]
            fn eq(&self, other: &FiniteSetEndofunctions<X>) -> bool {
                self.set == other.set
            }
        }
        #[automatically_derived]
        impl<X: ::core::cmp::Eq + FiniteSetSignature + EqSignature> ::core::cmp::Eq
        for FiniteSetEndofunctions<X> {
            #[inline]
            #[doc(hidden)]
            #[coverage(off)]
            fn assert_receiver_is_total_eq(&self) -> () {
                let _: ::core::cmp::AssertParamIsEq<X>;
            }
        }
        impl<X: FiniteSetSignature + EqSignature> FiniteSetEndofunctions<X> {
            pub fn new(set: X) -> Self {
                Self { set }
            }
        }
        impl<X: FiniteSetSignature + EqSignature> Signature
        for FiniteSetEndofunctions<X> {}
        impl<X: FiniteSetSignature + EqSignature> SetSignature
        for FiniteSetEndofunctions<X> {
            type Set = Vec<X::Set>;
            fn is_element(&self, f: &Self::Set) -> Result<(), String> {
                if f.len() != self.set.size() {
                    return Err(
                        "Function must have one value per element in the domain."
                            .to_string(),
                    );
                }
                for y in f {
                    self.set.is_element(y)?;
                }
                Ok(())
            }
        }
        impl<X: FiniteSetSignature + EqSignature> CountableSetSignature
        for FiniteSetEndofunctions<X> {
            fn generate_all_elements(&self) -> impl Iterator<Item = Self::Set> + Clone {
                (0..self.set.size())
                    .map(|_| self.set.list_all_elements())
                    .multi_cartesian_product()
            }
        }
        impl<X: FiniteSetSignature + EqSignature> FiniteSetSignature
        for FiniteSetEndofunctions<X> {}
        pub trait BorrowedMorphism<
            Domain: Signature,
            Range: Signature,
            M: Morphism<Domain, Range>,
        >: Borrow<M> + Clone + Debug + Send + Sync {}
        impl<
            Domain: Signature,
            Range: Signature,
            M: Morphism<Domain, Range>,
            BM: Borrow<M> + Clone + Debug + Send + Sync,
        > BorrowedMorphism<Domain, Range, M> for BM {}
    }
    mod orderings {
        use algebraeon_macros::{meta_signature, skip};
        use super::EqSignature;
        use std::{borrow::Borrow, cmp::Ordering};
        pub enum MergedSource {
            First,
            Second,
        }
        #[automatically_derived]
        impl ::core::fmt::Debug for MergedSource {
            #[inline]
            fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
                ::core::fmt::Formatter::write_str(
                    f,
                    match self {
                        MergedSource::First => "First",
                        MergedSource::Second => "Second",
                    },
                )
            }
        }
        #[automatically_derived]
        impl ::core::clone::Clone for MergedSource {
            #[inline]
            fn clone(&self) -> MergedSource {
                *self
            }
        }
        #[automatically_derived]
        impl ::core::marker::Copy for MergedSource {}
        #[automatically_derived]
        impl ::core::marker::StructuralPartialEq for MergedSource {}
        #[automatically_derived]
        impl ::core::cmp::PartialEq for MergedSource {
            #[inline]
            fn eq(&self, other: &MergedSource) -> bool {
                let __self_discr = ::core::intrinsics::discriminant_value(self);
                let __arg1_discr = ::core::intrinsics::discriminant_value(other);
                __self_discr == __arg1_discr
            }
        }
        #[automatically_derived]
        impl ::core::cmp::Eq for MergedSource {
            #[inline]
            #[doc(hidden)]
            #[coverage(off)]
            fn assert_receiver_is_total_eq(&self) -> () {}
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
        impl<'s, X, O: OrdSignature + 's, K: Fn(&X) -> &O::Set> Iterator
        for VecMerger<'s, X, O, K> {
            type Item = (MergedSource, X);
            fn next(&mut self) -> Option<Self::Item> {
                match (self.i < self.a.len(), self.j < self.b.len()) {
                    (true, true) => {
                        match self
                            .ordering
                            .cmp(
                                (self.key)(self.a[self.i].as_ref().unwrap()).borrow(),
                                (self.key)(self.b[self.j].as_ref().unwrap()).borrow(),
                            )
                        {
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
        pub trait PartialOrdSignature: EqSignature {
            fn partial_cmp(&self, a: &Self::Set, b: &Self::Set) -> Option<Ordering>;
        }
        pub trait MetaPartialOrdSignature: MetaType
        where
            Self::Signature: PartialOrdSignature,
        {
            fn partial_cmp(a: &Self, b: &Self) -> Option<Ordering> {
                Self::structure().partial_cmp(a, b)
            }
        }
        impl<T> MetaPartialOrdSignature for T
        where
            T: MetaType,
            T::Signature: PartialOrdSignature,
        {}
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
            fn is_sorted_by_key<X>(
                &self,
                a: &[X],
                key: impl Fn(&X) -> &Self::Set,
            ) -> bool {
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
            fn is_sorted_and_unique_by_key<X>(
                &self,
                a: &[X],
                key: impl Fn(&X) -> &Self::Set,
            ) -> bool {
                self.is_sorted_and_unique(&a.iter().map(key).collect::<Vec<_>>())
            }
            fn binary_search(
                &self,
                v: &[impl Borrow<Self::Set>],
                target: &Self::Set,
            ) -> bool {
                self.binary_search_by_key(v, target, |x| x.borrow()).is_some()
            }
            fn binary_search_by_key<'x, X>(
                &self,
                v: &'x [X],
                target: &Self::Set,
                key: impl Fn(&X) -> &Self::Set,
            ) -> Option<&'x X> {
                if true {
                    if !self.is_sorted_by_key(v, &key) {
                        ::core::panicking::panic(
                            "assertion failed: self.is_sorted_by_key(v, &key)",
                        )
                    }
                }
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
                    {
                        ::std::io::_print(format_args!("{0:?} {1:?}\n", a, b));
                    };
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
            fn merge_sorted<'s, S: Borrow<Self::Set> + 's>(
                &'s self,
                a: Vec<S>,
                b: Vec<S>,
            ) -> impl Iterator<Item = (MergedSource, S)> {
                self.merge_sorted_by_key(a, b, |x| (*x).borrow())
            }
            fn merge_sorted_by_key<X>(
                &self,
                a: Vec<X>,
                b: Vec<X>,
                key: impl Fn(&X) -> &Self::Set,
            ) -> impl Iterator<Item = (MergedSource, X)> {
                if true {
                    if !self.is_sorted_by_key(&a.iter().collect::<Vec<_>>(), |x| key(x))
                    {
                        ::core::panicking::panic(
                            "assertion failed: self.is_sorted_by_key(&a.iter().collect::<Vec<_>>(), |x| key(x))",
                        )
                    }
                }
                if true {
                    if !self.is_sorted_by_key(&b.iter().collect::<Vec<_>>(), |x| key(x))
                    {
                        ::core::panicking::panic(
                            "assertion failed: self.is_sorted_by_key(&b.iter().collect::<Vec<_>>(), |x| key(x))",
                        )
                    }
                }
                VecMerger::new(self, a, b, key)
            }
            fn sort<S: Borrow<Self::Set>>(&self, a: Vec<S>) -> Vec<S> {
                self.sort_by_key(a, &|x| x.borrow())
            }
            fn sort_by_key<X>(
                &self,
                mut a: Vec<X>,
                key: &impl Fn(&X) -> &Self::Set,
            ) -> Vec<X> {
                match a.len() {
                    0 | 1 => a,
                    2 => {
                        match self.cmp(key(&a[0]).borrow(), key(&a[1]).borrow()) {
                            Ordering::Less | Ordering::Equal => a,
                            Ordering::Greater => {
                                a.swap(0, 1);
                                a
                            }
                        }
                    }
                    n => {
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
            fn sort_by_cached_by<X: 'static>(
                &self,
                a: Vec<X>,
                key: impl Fn(&X) -> Self::Set,
            ) -> Vec<X>
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
        pub trait MetaOrdSignature: MetaType
        where
            Self::Signature: OrdSignature,
        {
            fn cmp(a: &Self, b: &Self) -> Ordering {
                Self::structure().cmp(a, b)
            }
            fn max(a: &Self, b: &Self) -> Self {
                Self::structure().max(a, b)
            }
            fn min(a: &Self, b: &Self) -> Self {
                Self::structure().min(a, b)
            }
            fn is_sorted(a: &[impl Borrow<Self>]) -> bool {
                Self::structure().is_sorted(a)
            }
            fn is_sorted_by_key<X>(a: &[X], key: impl Fn(&X) -> &Self) -> bool {
                Self::structure().is_sorted_by_key(a, key)
            }
            fn is_sorted_and_unique(a: &[impl Borrow<Self>]) -> bool {
                Self::structure().is_sorted_and_unique(a)
            }
            fn is_sorted_and_unique_by_key<X>(
                a: &[X],
                key: impl Fn(&X) -> &Self,
            ) -> bool {
                Self::structure().is_sorted_and_unique_by_key(a, key)
            }
            fn binary_search(v: &[impl Borrow<Self>], target: &Self) -> bool {
                Self::structure().binary_search(v, target)
            }
            fn binary_search_by_key<'x, X>(
                v: &'x [X],
                target: &Self,
                key: impl Fn(&X) -> &Self,
            ) -> Option<&'x X> {
                Self::structure().binary_search_by_key(v, target, key)
            }
            fn merge_sorted_by_key<X>(
                a: Vec<X>,
                b: Vec<X>,
                key: impl Fn(&X) -> &Self,
            ) -> impl Iterator<Item = (MergedSource, X)> {
                Self::structure().merge_sorted_by_key(a, b, key)
            }
            fn sort<S: Borrow<Self>>(a: Vec<S>) -> Vec<S> {
                Self::structure().sort(a)
            }
            fn sort_by_cached_by<X: 'static>(
                a: Vec<X>,
                key: impl Fn(&X) -> Self,
            ) -> Vec<X>
            where
                Self: 'static,
            {
                Self::structure().sort_by_cached_by(a, key)
            }
        }
        impl<T> MetaOrdSignature for T
        where
            T: MetaType,
            T::Signature: OrdSignature,
        {}
    }
    mod pairs {
        use crate::structure::Signature;
        use super::{EqSignature, SetSignature};
        use std::fmt::Debug;
        /// The set of Pairs
        pub struct PairsStructure<S> {
            set: S,
        }
        #[automatically_derived]
        impl<S: ::core::fmt::Debug> ::core::fmt::Debug for PairsStructure<S> {
            #[inline]
            fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
                ::core::fmt::Formatter::debug_struct_field1_finish(
                    f,
                    "PairsStructure",
                    "set",
                    &&self.set,
                )
            }
        }
        #[automatically_derived]
        impl<S: ::core::clone::Clone> ::core::clone::Clone for PairsStructure<S> {
            #[inline]
            fn clone(&self) -> PairsStructure<S> {
                PairsStructure {
                    set: ::core::clone::Clone::clone(&self.set),
                }
            }
        }
        #[automatically_derived]
        impl<S> ::core::marker::StructuralPartialEq for PairsStructure<S> {}
        #[automatically_derived]
        impl<S: ::core::cmp::PartialEq> ::core::cmp::PartialEq for PairsStructure<S> {
            #[inline]
            fn eq(&self, other: &PairsStructure<S>) -> bool {
                self.set == other.set
            }
        }
        #[automatically_derived]
        impl<S: ::core::cmp::Eq> ::core::cmp::Eq for PairsStructure<S> {
            #[inline]
            #[doc(hidden)]
            #[coverage(off)]
            fn assert_receiver_is_total_eq(&self) -> () {
                let _: ::core::cmp::AssertParamIsEq<S>;
            }
        }
        impl<S: SetSignature> PairsStructure<S> {
            pub fn new(set: S) -> Self {
                Self { set }
            }
            /// Construct a new pair from two elements of a set
            pub fn new_pair(
                &self,
                a: S::Set,
                b: S::Set,
            ) -> Result<(S::Set, S::Set), String> {
                self.set.is_element(&a)?;
                self.set.is_element(&b)?;
                Ok((a, b))
            }
        }
        impl<S: SetSignature> Signature for PairsStructure<S> {}
        impl<S: SetSignature> SetSignature for PairsStructure<S> {
            type Set = (S::Set, S::Set);
            fn is_element(&self, x: &Self::Set) -> Result<(), String> {
                self.set.is_element(&x.0)?;
                self.set.is_element(&x.1)?;
                Ok(())
            }
        }
        impl<S: SetSignature + EqSignature> EqSignature for PairsStructure<S> {
            fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
                self.set.equal(&a.0, &b.0) && self.set.equal(&a.1, &b.1)
            }
        }
        /// The set of unordered Pairs of distinct elements
        pub struct UnorderedPairs<Set> {
            set: Set,
        }
        #[automatically_derived]
        impl<Set: ::core::fmt::Debug> ::core::fmt::Debug for UnorderedPairs<Set> {
            #[inline]
            fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
                ::core::fmt::Formatter::debug_struct_field1_finish(
                    f,
                    "UnorderedPairs",
                    "set",
                    &&self.set,
                )
            }
        }
        #[automatically_derived]
        impl<Set: ::core::clone::Clone> ::core::clone::Clone for UnorderedPairs<Set> {
            #[inline]
            fn clone(&self) -> UnorderedPairs<Set> {
                UnorderedPairs {
                    set: ::core::clone::Clone::clone(&self.set),
                }
            }
        }
        #[automatically_derived]
        impl<Set> ::core::marker::StructuralPartialEq for UnorderedPairs<Set> {}
        #[automatically_derived]
        impl<Set: ::core::cmp::PartialEq> ::core::cmp::PartialEq
        for UnorderedPairs<Set> {
            #[inline]
            fn eq(&self, other: &UnorderedPairs<Set>) -> bool {
                self.set == other.set
            }
        }
        #[automatically_derived]
        impl<Set: ::core::cmp::Eq> ::core::cmp::Eq for UnorderedPairs<Set> {
            #[inline]
            #[doc(hidden)]
            #[coverage(off)]
            fn assert_receiver_is_total_eq(&self) -> () {
                let _: ::core::cmp::AssertParamIsEq<Set>;
            }
        }
        impl<Set: SetSignature> UnorderedPairs<Set> {
            pub fn new(set: Set) -> Self {
                Self { set }
            }
        }
        pub struct UnorderedPair<T>(T, T);
        #[automatically_derived]
        impl<T: ::core::fmt::Debug> ::core::fmt::Debug for UnorderedPair<T> {
            #[inline]
            fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
                ::core::fmt::Formatter::debug_tuple_field2_finish(
                    f,
                    "UnorderedPair",
                    &self.0,
                    &&self.1,
                )
            }
        }
        #[automatically_derived]
        impl<T: ::core::clone::Clone> ::core::clone::Clone for UnorderedPair<T> {
            #[inline]
            fn clone(&self) -> UnorderedPair<T> {
                UnorderedPair(
                    ::core::clone::Clone::clone(&self.0),
                    ::core::clone::Clone::clone(&self.1),
                )
            }
        }
        #[automatically_derived]
        impl<T> ::core::marker::StructuralPartialEq for UnorderedPair<T> {}
        #[automatically_derived]
        impl<T: ::core::cmp::PartialEq> ::core::cmp::PartialEq for UnorderedPair<T> {
            #[inline]
            fn eq(&self, other: &UnorderedPair<T>) -> bool {
                self.0 == other.0 && self.1 == other.1
            }
        }
        #[automatically_derived]
        impl<T: ::core::cmp::Eq> ::core::cmp::Eq for UnorderedPair<T> {
            #[inline]
            #[doc(hidden)]
            #[coverage(off)]
            fn assert_receiver_is_total_eq(&self) -> () {
                let _: ::core::cmp::AssertParamIsEq<T>;
            }
        }
        #[automatically_derived]
        impl<T: ::core::hash::Hash> ::core::hash::Hash for UnorderedPair<T> {
            #[inline]
            fn hash<__H: ::core::hash::Hasher>(&self, state: &mut __H) -> () {
                ::core::hash::Hash::hash(&self.0, state);
                ::core::hash::Hash::hash(&self.1, state)
            }
        }
        impl<S: SetSignature + EqSignature> UnorderedPairs<S> {
            pub fn new_pair(
                &self,
                a: &S::Set,
                b: &S::Set,
            ) -> Result<UnorderedPair<S::Set>, String> {
                if self.set.equal(a, b) {
                    Err("UnorderedPair elements must be distinct".to_string())
                } else {
                    Ok(UnorderedPair(a.clone(), b.clone()))
                }
            }
        }
        impl<S: SetSignature> Signature for UnorderedPairs<S> {}
        impl<S: SetSignature + EqSignature> SetSignature for UnorderedPairs<S> {
            type Set = UnorderedPair<S::Set>;
            fn is_element(&self, x: &Self::Set) -> Result<(), String> {
                self.set.is_element(&x.0)?;
                self.set.is_element(&x.1)?;
                if self.set.equal(&x.0, &x.1) {
                    Err("UnorderedPair elements must be distinct".to_string())
                } else {
                    Ok(())
                }
            }
        }
        impl<S: SetSignature + EqSignature> EqSignature for UnorderedPairs<S> {
            fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
                self.set.equal(&a.0, &b.0) && self.set.equal(&a.1, &b.1)
                    || self.set.equal(&a.0, &b.1) && self.set.equal(&a.1, &b.0)
            }
        }
    }
    mod singleton_set {
        use super::{EqSignature, OrdSignature, SetSignature, Signature};
        use crate::structure::{
            CountableSetSignature, FiniteSetSignature, orderings::PartialOrdSignature,
        };
        use std::fmt::Debug;
        pub struct SingletonSetStructure {}
        #[automatically_derived]
        impl ::core::default::Default for SingletonSetStructure {
            #[inline]
            fn default() -> SingletonSetStructure {
                SingletonSetStructure {}
            }
        }
        #[automatically_derived]
        impl ::core::fmt::Debug for SingletonSetStructure {
            #[inline]
            fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
                ::core::fmt::Formatter::write_str(f, "SingletonSetStructure")
            }
        }
        #[automatically_derived]
        impl ::core::clone::Clone for SingletonSetStructure {
            #[inline]
            fn clone(&self) -> SingletonSetStructure {
                SingletonSetStructure {}
            }
        }
        #[automatically_derived]
        impl ::core::marker::StructuralPartialEq for SingletonSetStructure {}
        #[automatically_derived]
        impl ::core::cmp::PartialEq for SingletonSetStructure {
            #[inline]
            fn eq(&self, other: &SingletonSetStructure) -> bool {
                true
            }
        }
        #[automatically_derived]
        impl ::core::cmp::Eq for SingletonSetStructure {
            #[inline]
            #[doc(hidden)]
            #[coverage(off)]
            fn assert_receiver_is_total_eq(&self) -> () {}
        }
        impl Signature for SingletonSetStructure {}
        impl SetSignature for SingletonSetStructure {
            type Set = ();
            fn is_element(&self, _: &Self::Set) -> Result<(), String> {
                Ok(())
            }
        }
        impl EqSignature for SingletonSetStructure {
            fn equal(&self, x: &Self::Set, y: &Self::Set) -> bool {
                x == y
            }
        }
        impl PartialOrdSignature for SingletonSetStructure {
            fn partial_cmp(
                &self,
                a: &Self::Set,
                b: &Self::Set,
            ) -> Option<std::cmp::Ordering> {
                Some(self.cmp(a, b))
            }
        }
        impl OrdSignature for SingletonSetStructure {
            fn cmp(&self, x: &Self::Set, y: &Self::Set) -> std::cmp::Ordering {
                x.cmp(y)
            }
        }
        impl CountableSetSignature for SingletonSetStructure {
            fn generate_all_elements(&self) -> impl Iterator<Item = Self::Set> + Clone {
                [()].into_iter()
            }
        }
        impl FiniteSetSignature for SingletonSetStructure {
            fn size(&self) -> usize {
                1
            }
        }
    }
    #[allow(clippy::module_inception)]
    mod structure {
        use paste::paste;
        use rand::{Rng, SeedableRng, rngs::StdRng};
        use std::{borrow::Borrow, fmt::Debug};
        pub trait Signature: Clone + Debug + Eq + Send + Sync {}
        /// Instances of a type implementing this trait represent
        /// a set of elements of type `Self::Set` with some
        /// structure, for example, the structure of a ring.
        pub trait SetSignature: Signature {
            type Set: Clone + Debug + Send + Sync;
            /// Some instances of `Self::Set` may not be valid to represent elements of this set.
            /// Return `true` if `x` is a valid element and `false` if not.
            fn is_element(&self, x: &Self::Set) -> Result<(), String>;
        }
        pub trait MetaType: Clone + Debug {
            type Signature: SetSignature<Set = Self>;
            fn structure() -> Self::Signature;
        }
        pub trait ToStringSignature: SetSignature {
            fn to_string(&self, elem: &Self::Set) -> String;
        }
        pub trait EqSignature: SetSignature {
            fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool;
        }
        pub trait CountableSetSignature: SetSignature {
            /// Yield distinct elements of the set such that every element eventually appears.
            /// Always yields elements in the same order.
            fn generate_all_elements(&self) -> impl Iterator<Item = Self::Set> + Clone;
        }
        pub trait FiniteSetSignature: CountableSetSignature {
            /// A list of all elements in the set.
            /// Always returns elements in the same order.
            fn list_all_elements(&self) -> Vec<Self::Set> {
                self.generate_all_elements().collect()
            }
            fn size(&self) -> usize {
                self.list_all_elements().len()
            }
            fn generate_random_elements(
                &self,
                seed: u64,
            ) -> impl Iterator<Item = Self::Set> + Clone {
                let rng = StdRng::seed_from_u64(seed);
                FiniteSetRandomElementGenerator::<Self, StdRng> {
                    all_elements: self.list_all_elements(),
                    rng,
                }
            }
        }
        pub trait MaybeFiniteSetSignature: SetSignature {
            type FiniteSetStructure: FiniteSetSignature<Set = Self::Set>;
            #[allow(clippy::result_unit_err)]
            fn finite_set_structure(&self) -> Result<Self::FiniteSetStructure, ()>;
        }
        pub struct FiniteSetRandomElementGenerator<S: FiniteSetSignature, R: Rng> {
            all_elements: Vec<S::Set>,
            rng: R,
        }
        #[automatically_derived]
        impl<
            S: ::core::fmt::Debug + FiniteSetSignature,
            R: ::core::fmt::Debug + Rng,
        > ::core::fmt::Debug for FiniteSetRandomElementGenerator<S, R>
        where
            S::Set: ::core::fmt::Debug,
        {
            #[inline]
            fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
                ::core::fmt::Formatter::debug_struct_field2_finish(
                    f,
                    "FiniteSetRandomElementGenerator",
                    "all_elements",
                    &self.all_elements,
                    "rng",
                    &&self.rng,
                )
            }
        }
        #[automatically_derived]
        impl<
            S: ::core::clone::Clone + FiniteSetSignature,
            R: ::core::clone::Clone + Rng,
        > ::core::clone::Clone for FiniteSetRandomElementGenerator<S, R>
        where
            S::Set: ::core::clone::Clone,
        {
            #[inline]
            fn clone(&self) -> FiniteSetRandomElementGenerator<S, R> {
                FiniteSetRandomElementGenerator {
                    all_elements: ::core::clone::Clone::clone(&self.all_elements),
                    rng: ::core::clone::Clone::clone(&self.rng),
                }
            }
        }
        impl<S: FiniteSetSignature, R: Rng> Iterator
        for FiniteSetRandomElementGenerator<S, R> {
            type Item = S::Set;
            fn next(&mut self) -> Option<Self::Item> {
                if self.all_elements.is_empty() {
                    None
                } else {
                    let idx = self.rng.random_range(0..self.all_elements.len());
                    Some(self.all_elements[idx].clone())
                }
            }
        }
        pub trait BorrowedSet<S>: Borrow<S> + Clone + Debug + Send + Sync {}
        impl<S, BS: Borrow<S> + Clone + Debug + Send + Sync> BorrowedSet<S> for BS {}
        pub trait BorrowedStructure<
            S: Signature,
        >: Borrow<S> + Clone + Debug + Eq + Send + Sync {}
        impl<
            S: Signature,
            BS: Borrow<S> + Clone + Debug + Eq + Send + Sync,
        > BorrowedStructure<S> for BS {}
    }
    pub use algebraeon_macros::CanonicalStructure;
    pub use empty_set::EmptySetStructure;
    pub use finite_set::EnumeratedFiniteSetStructure;
    pub use morphism::{
        BijectiveFunction, BorrowedMorphism, CompositionMorphism, Endofunction,
        Endomorphism, FiniteSetEndofunctions, Function, Functions, IdentityMorphism,
        InjectiveFunction, Morphism, Permutation,
    };
    pub use orderings::{OrdSignature, PartialOrdSignature};
    pub use pairs::{PairsStructure, UnorderedPair, UnorderedPairs};
    pub use singleton_set::SingletonSetStructure;
    pub use structure::{
        BorrowedSet, BorrowedStructure, CountableSetSignature, EqSignature,
        FiniteSetSignature, MaybeFiniteSetSignature, MetaType, SetSignature, Signature,
        ToStringSignature,
    };
}
