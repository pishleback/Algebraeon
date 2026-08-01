use crate::*;
use std::{cmp::Ordering, marker::PhantomData};

/// A sized finite set from an unsized finite set
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CartesianProductSetStructure<
    Set0: SetSignature,
    Set0B: BorrowedStructure<Set0>,
    Set1: SetSignature,
    Set1B: BorrowedStructure<Set1>,
> {
    _set_0: PhantomData<Set0>,
    set_0: Set0B,
    _set_1: PhantomData<Set1>,
    set_1: Set1B,
}

impl<
    Set0: SetSignature,
    Set0B: BorrowedStructure<Set0>,
    Set1: SetSignature,
    Set1B: BorrowedStructure<Set1>,
> CartesianProductSetStructure<Set0, Set0B, Set1, Set1B>
{
    pub fn new(set_0: Set0B, set_1: Set1B) -> Self {
        Self {
            _set_0: PhantomData,
            set_0,
            _set_1: PhantomData,
            set_1,
        }
    }

    pub fn set_0(&self) -> &Set0 {
        self.set_0.borrow()
    }

    pub fn set_1(&self) -> &Set1 {
        self.set_1.borrow()
    }
}

impl<
    Set0: SetSignature,
    Set0B: BorrowedStructure<Set0>,
    Set1: SetSignature,
    Set1B: BorrowedStructure<Set1>,
> Signature for CartesianProductSetStructure<Set0, Set0B, Set1, Set1B>
{
}

impl<
    Set0: SetSignature,
    Set0B: BorrowedStructure<Set0>,
    Set1: SetSignature,
    Set1B: BorrowedStructure<Set1>,
> SetSignature for CartesianProductSetStructure<Set0, Set0B, Set1, Set1B>
{
    type Elem = (Set0::Elem, Set1::Elem);

    fn validate_element(&self, elem: &Self::Elem) -> Result<(), String> {
        self.set_0().validate_element(&elem.0)?;
        self.set_1().validate_element(&elem.1)?;
        Ok(())
    }
}

impl<
    Set0: EqSignature,
    Set0B: BorrowedStructure<Set0>,
    Set1: EqSignature,
    Set1B: BorrowedStructure<Set1>,
> EqSignature for CartesianProductSetStructure<Set0, Set0B, Set1, Set1B>
{
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        self.set_0().equal(&a.0, &b.0) && self.set_1().equal(&a.1, &b.1)
    }
}

impl<
    Set0: OrdSignature,
    Set0B: BorrowedStructure<Set0>,
    Set1: OrdSignature,
    Set1B: BorrowedStructure<Set1>,
> PartialOrdSignature for CartesianProductSetStructure<Set0, Set0B, Set1, Set1B>
{
    fn partial_cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Ordering> {
        Some(self.cmp(a, b))
    }
}

/// Compare element 0 first, then element 1
impl<
    Set0: OrdSignature,
    Set0B: BorrowedStructure<Set0>,
    Set1: OrdSignature,
    Set1B: BorrowedStructure<Set1>,
> OrdSignature for CartesianProductSetStructure<Set0, Set0B, Set1, Set1B>
{
    fn cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        match self.set_0().cmp(&a.0, &b.0) {
            Ordering::Less => Ordering::Less,
            Ordering::Greater => Ordering::Greater,
            Ordering::Equal => self.set_1().cmp(&a.1, &b.1),
        }
    }
}

// TODO: allow this over a product of countable sets
impl<
    Set0: FiniteSetSignature,
    Set0B: BorrowedStructure<Set0>,
    Set1: FiniteSetSignature,
    Set1B: BorrowedStructure<Set1>,
> CountableSetSignature for CartesianProductSetStructure<Set0, Set0B, Set1, Set1B>
{
    fn into_generate_all_elements(self) -> impl Iterator<Item = Self::Elem> {
        self.list_all_elements().into_iter()
    }

    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> {
        self.list_all_elements().into_iter()
    }
}

impl<
    Set0: FiniteSetSignature,
    Set0B: BorrowedStructure<Set0>,
    Set1: FiniteSetSignature,
    Set1B: BorrowedStructure<Set1>,
> FiniteSetSignature for CartesianProductSetStructure<Set0, Set0B, Set1, Set1B>
{
    fn list_all_elements(&self) -> Vec<Self::Elem> {
        let mut all = vec![];
        for elem_0 in self.set_0().list_all_elements() {
            for elem_1 in self.set_1().list_all_elements() {
                all.push((elem_0.clone(), elem_1));
            }
        }
        all
    }

    fn size(&self) -> Natural {
        self.set_0().size() * self.set_1().size()
    }

    fn generate_random_elements(&self, seed: u64) -> impl Iterator<Item = Self::Elem> {
        let mut set0_rand = self.set_0().generate_random_elements(seed);
        let mut set1_rand = self.set_1().generate_random_elements(seed + 1);
        (0usize..).map(move |_| (set0_rand.next().unwrap(), set1_rand.next().unwrap()))
    }
}

impl<
    Set0: OrderedFiniteSetSignature,
    Set0B: BorrowedStructure<Set0>,
    Set1: OrderedFiniteSetSignature,
    Set1B: BorrowedStructure<Set1>,
> OrderedFiniteSetSignature for CartesianProductSetStructure<Set0, Set0B, Set1, Set1B>
{
    fn list_all_elements_ordered(&self) -> Vec<Self::Elem> {
        self.list_all_elements()
    }

    fn element_to_enumeration(&self, elem: &Self::Elem) -> Natural {
        self.set_0().element_to_enumeration(&elem.0) * self.set_1().size()
            + self.set_1().element_to_enumeration(&elem.1)
    }

    fn enumeration_to_element(&self, num: &Natural) -> Option<Self::Elem> {
        if *num < self.size() {
            let (q, r) = num.div_mod(self.set_1().size());
            Some((
                self.set_0().enumeration_to_element(&q).unwrap(),
                self.set_1().enumeration_to_element(&r).unwrap(),
            ))
        } else {
            None
        }
    }
}

// Todo: but needs https://github.com/rust-lang/rust/issues/76560
// For now you must call .try_into_sized::<N01>().unwrap() to get a sized cartesian product
// impl<
//     const N0: usize,
//     const N1: usize,
//     Set0: FiniteSetSizedSignature<N0>,
//     Set0B: BorrowedStructure<Set0>,
//     Set1: FiniteSetSizedSignature<N1>,
//     Set1B: BorrowedStructure<Set1>,
// > FiniteSetSizedSignature<{ N0 * N1 }> for CartesianProductSetStructure<Set0, Set0B, Set1, Set1B>
// {
// }

impl<Elem0: MetaType, Elem1: MetaType> MetaType for (Elem0, Elem1) {
    type Signature = CartesianProductSetStructure<
        Elem0::Signature,
        Elem0::Signature,
        Elem1::Signature,
        Elem1::Signature,
    >;

    fn structure() -> Self::Signature {
        CartesianProductSetStructure::new(Elem0::structure(), Elem1::structure())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test() {
        let set_0 = i32::structure().into_const_size_finite_subset([1, 2, 3, 4, 5]);
        let set_1 = i32::structure().into_const_size_finite_subset([6, 7, 8]);
        let set_01 = CartesianProductSetStructure::new(set_0, set_1);

        assert!(set_01.is_element(&(1, 6)));
        assert!(set_01.is_element(&(5, 8)));
        assert!(!set_01.is_element(&(1, 1)));
        assert!(!set_01.is_element(&(6, 6)));
    }

    #[test]
    fn enumeration() {
        let set_0 = i32::structure().into_const_size_finite_subset([1, 2, 3, 4, 5]);
        let set_1 = i32::structure().into_const_size_finite_subset([6, 7, 8]);
        let set_01 = CartesianProductSetStructure::new(set_0, set_1);
        assert_enumerated_ord_finite_set!(set_01, 15);
    }
}
