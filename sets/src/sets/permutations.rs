use algebraeon_structures::*;
use std::fmt::Display;
use std::marker::PhantomData;
use std::ops::Mul;

#[derive(Debug, Clone)]
pub struct Cycle<T> {
    pub cycle: Vec<T>,
}

#[derive(Debug, Clone)]
pub struct FinitelySupportedPermutation<T> {
    pub perm: Vec<(T, T)>,
}

/*

impl<T: PartialEq + Eq + Clone> Cycle<T> {
    fn new_unchecked(cycle: Vec<T>) -> Self {
        debug_assert_eq!(
            cycle.len(),
            cycle.iter().cloned().collect::<HashSet<_>>().len()
        );
        Self { cycle }
    }

    pub fn len(&self) -> usize {
        self.cycle.len()
    }
}

impl<T: Display> Display for Cycle<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "(")?;
        for (i, x) in self.cycle.iter().enumerate() {
            if i > 0 {
                write!(f, " ")?;
            }
            write!(f, "{}", x)?;
        }
        write!(f, ")")?;
        Ok(())
    }
}

impl<T: PartialEq + Eq + Hash> PartialEq for FinitelySupportedPermutation<T> {
    fn eq(&self, other: &Self) -> bool {
        self.right == other.right
    }
}

impl<T: PartialEq + Eq + Hash> Eq for FinitelySupportedPermutation<T> {}

impl<T> FinitelySupportedPermutation<T> {
    pub fn support_size(&self) -> usize {
        self.perm.len()
    }
}

impl<T: PartialEq + Eq + Hash + Clone> From<Cycle<T>> for FinitelySupportedPermutation<T> {
    fn from(cycle: Cycle<T>) -> Self {
        let n = cycle.len();
        Self::new_unchecked(
            (0..n)
                .map(|i| (cycle.cycle[i].clone(), cycle.cycle[(i + 1) % n].clone()))
                .collect(),
        )
    }
}

impl<T: Display + PartialEq + Eq + Hash + Clone> Display for FinitelySupportedPermutation<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.support_size() == 0 {
            write!(f, "()")?;
        } else {
            for cycle in self.disjoint_cycles() {
                write!(f, "{}", cycle)?;
            }
        }
        Ok(())
    }
}

impl<T: PartialEq + Eq + Hash + Clone> FinitelySupportedPermutation<T> {
    pub fn identity() -> Self {
        Self {
            perm: vec![],
            right: HashMap::new(),
            left: HashMap::new(),
        }
    }

    fn new_unchecked(perm: Vec<(T, T)>) -> Self {
        let right = perm.iter().cloned().collect::<HashMap<T, T>>();
        let left = perm
            .iter()
            .cloned()
            .map(|(a, b)| (b, a))
            .collect::<HashMap<T, T>>();
        debug_assert_eq!(perm.len(), right.len());
        debug_assert_eq!(perm.len(), left.len());
        Self { perm, left, right }
    }

    pub fn new(perm: Vec<(T, T)>) -> Result<Self, ()> {
        let right = perm.iter().cloned().collect::<HashMap<T, T>>();
        let left = perm
            .iter()
            .cloned()
            .map(|(a, b)| (b, a))
            .collect::<HashMap<T, T>>();
        if perm.len() != right.len() || perm.len() != left.len() {
            return Err(());
        }
        Ok(Self { perm, left, right })
    }

    pub fn map_injective_unchecked<S: PartialEq + Eq + Hash + Clone>(
        self,
        f: impl Fn(T) -> S,
    ) -> FinitelySupportedPermutation<S> {
        FinitelySupportedPermutation::new_unchecked(
            self.perm.into_iter().map(|(a, b)| (f(a), f(b))).collect(),
        )
    }

    // pub fn from_fn(f: impl Fn(T) -> T) -> Self
    // where
    //     T: Enumerated,
    // {
    //     let points = T::points().collect::<Vec<_>>();
    //     let images = T::points().map(f).collect::<Vec<_>>();
    //     debug_assert_eq!(points.len(), T::N);
    //     debug_assert_eq!(images.len(), T::N);
    //     assert_eq!(T::N, images.iter().collect::<HashSet<_>>().len());
    //     Self::from_perm_unchecked(points.into_iter().zip(images).collect())
    // }

    pub fn new_swap(t1: &T, t2: &T) -> Self {
        if t1 == t2 {
            Self::identity()
        } else {
            Self::new_unchecked(vec![(t1.clone(), t2.clone()), (t2.clone(), t1.clone())])
        }
    }

    pub fn new_cycle(ts: Vec<&T>) -> Self {
        // Check the items are unique
        let n = ts.len();
        assert_eq!(ts.iter().collect::<HashSet<_>>().len(), n);
        Self::new_unchecked(
            (0..n)
                .map(|i| (ts[i].clone(), ts[(i + 1) % n].clone()))
                .collect(),
        )
    }

    pub fn inverse(self) -> Self {
        Self {
            perm: self.perm.into_iter().map(|(a, b)| (b, a)).collect(),
            right: self.left,
            left: self.right,
        }
    }

    pub fn apply<'a>(&'a self, t: &'a T) -> &'a T {
        self.right.get(t).unwrap_or(t)
    }

    pub fn apply_inverse<'a>(&'a self, t: &'a T) -> &'a T {
        self.left.get(t).unwrap_or(t)
    }

    pub fn disjoint_cycles(&self) -> Vec<Cycle<T>> {
        let mut cycles = vec![];
        let mut used = HashSet::new();
        for t in self.right.keys() {
            if !used.contains(t) {
                let mut cycle = vec![];
                let mut s = t;
                loop {
                    cycle.push(s.clone());
                    used.insert(s);
                    s = self.right.get(s).unwrap();
                    if s == t {
                        break;
                    }
                }
                if cycle.len() >= 2 {
                    cycles.push(Cycle::new_unchecked(cycle));
                }
            }
        }
        cycles
    }

    pub fn is_even(&self) -> bool {
        let mut is_even = true;
        for cycle in self.disjoint_cycles() {
            if cycle.len() % 2 == 0 {
                is_even = !is_even;
            }
        }
        is_even
    }

    pub fn is_odd(&self) -> bool {
        !self.is_even()
    }
}

impl<T: PartialEq + Eq + Hash> Mul<&FinitelySupportedPermutation<T>>
    for &FinitelySupportedPermutation<T>
where
    T: Clone,
{
    type Output = FinitelySupportedPermutation<T>;

    fn mul(self, other: &FinitelySupportedPermutation<T>) -> Self::Output {
        FinitelySupportedPermutation::new_unchecked(
            self.right
                .keys()
                .collect::<HashSet<_>>()
                .union(&other.right.keys().collect::<HashSet<_>>())
                .map(|t| ((*t).clone(), other.apply(self.apply(t)).clone()))
                .filter(|(a, b)| a != b)
                .collect(),
        )
    }
}
*/

impl<Elem: MetaType> MetaType for FinitelySupportedPermutation<Elem>
where
    Elem::Signature: OrdSignature,
{
    type Signature = FinitelySupportedPermutationsStructure<Elem::Signature, Elem::Signature>;

    fn structure() -> Self::Signature {
        FinitelySupportedPermutationsStructure::new(Elem::structure())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FinitelySupportedPermutationsStructure<Set: SetSignature, SetB: BorrowedStructure<Set>> {
    _set: PhantomData<Set>,
    set: SetB,
}

impl<Set: SetSignature, SetB: BorrowedStructure<Set>>
    FinitelySupportedPermutationsStructure<Set, SetB>
{
    pub fn new(set: SetB) -> Self {
        Self {
            _set: PhantomData,
            set,
        }
    }
}

impl<Set: SetSignature, SetB: BorrowedStructure<Set>> Signature
    for FinitelySupportedPermutationsStructure<Set, SetB>
{
}

impl<Set: OrdSignature, SetB: BorrowedStructure<Set>> SetSignature
    for FinitelySupportedPermutationsStructure<Set, SetB>
{
    type Elem = FinitelySupportedPermutation<Set::Elem>;

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        todo!()
    }
}

// impl<Set: SetSignature, SetB: BorrowedStructure<Set>> CompositionSignature
//     for FinitelySupportedPermutationsStructure<Set, SetB>
// {
//     fn compose(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
//         todo!()
//     }
// }

// impl<Set: SetSignature, SetB: BorrowedStructure<Set>> AssociativeCompositionSignature
//     for FinitelySupportedPermutationsStructure<Set, SetB>
// {
// }

// impl<Set: SetSignature, SetB: BorrowedStructure<Set>> LeftCancellativeCompositionSignature
//     for FinitelySupportedPermutationsStructure<Set, SetB>
// {
//     fn try_left_difference(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
//         Some(self.compose(&self.inverse(b), a))
//     }
// }

// impl<Set: SetSignature, SetB: BorrowedStructure<Set>> RightCancellativeCompositionSignature
//     for FinitelySupportedPermutationsStructure<Set, SetB>
// {
//     fn try_right_difference(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
//         Some(self.compose(a, &self.inverse(b)))
//     }
// }

// impl<Set: SetSignature, SetB: BorrowedStructure<Set>> IdentitySignature
//     for FinitelySupportedPermutationsStructure<Set, SetB>
// {
//     fn identity(&self) -> Self::Elem {
//         todo!()
//     }
// }

// impl<Set: SetSignature, SetB: BorrowedStructure<Set>> MonoidSignature
//     for FinitelySupportedPermutationsStructure<Set, SetB>
// {
// }

// impl<Set: SetSignature, SetB: BorrowedStructure<Set>> TryLeftInverseSignature
//     for FinitelySupportedPermutationsStructure<Set, SetB>
// {
//     fn try_left_inverse(&self, a: &Self::Elem) -> Option<Self::Elem> {
//         Some(self.inverse(a))
//     }
// }

// impl<Set: SetSignature, SetB: BorrowedStructure<Set>> TryRightInverseSignature
//     for FinitelySupportedPermutationsStructure<Set, SetB>
// {
//     fn try_right_inverse(&self, a: &Self::Elem) -> Option<Self::Elem> {
//         Some(self.inverse(a))
//     }
// }

// impl<Set: SetSignature, SetB: BorrowedStructure<Set>> TryInverseSignature
//     for FinitelySupportedPermutationsStructure<Set, SetB>
// {
//     fn try_inverse(&self, a: &Self::Elem) -> Option<Self::Elem> {
//         Some(self.inverse(a))
//     }
// }

// impl<Set: SetSignature, SetB: BorrowedStructure<Set>> GroupSignature
//     for FinitelySupportedPermutationsStructure<Set, SetB>
// {
//     fn inverse(&self, a: &Self::Elem) -> Self::Elem {
//         todo!()
//     }
// }

#[cfg(test)]
mod partition_tests {
    use super::*;

    #[test]
    fn test_valid() {
        let x = FinitelySupportedPermutation {
            perm: vec![(1, 2), (2, 3), (3, 1)],
        };

        println!("{:?}", x.validate_element());
    }

    // #[test]
    // fn test_permutations() {
    //     let x = FinitelySupportedPermutation::new_cycle(vec![&1, &2, &3]);
    //     let y = FinitelySupportedPermutation::new_swap(&4, &5);

    //     debug_assert_eq!(x.is_even(), true);
    //     debug_assert_eq!(y.is_even(), false);
    // }
}
