use algebraeon_macros::signature_meta_trait;
use algebraeon_sets::sets::{
    FiniteSetToFinitelySupportedPermutationsStructure, FiniteSubsetByOrd,
    FinitelySupportedPermutation, SetToFiniteSubsetsByOrdSignature,
    SetToFixedSizeFiniteSubsetsByOrdSignature,
};
use algebraeon_structures::*;
use std::{cmp::Ordering, marker::PhantomData};

/// A 2-element subset of a 6-element set
#[derive(Debug, Clone)]
pub struct Duad<Point> {
    // must have p1 < p2
    pub(crate) points: [Point; 2],
}

impl<Point: MetaType> MetaType for Duad<Point>
where
    Point::Signature: ConstSizeFiniteSetSignature<6> + EnumeratedOrdFiniteSetSignature,
{
    type Signature = DuadsStructure<Point::Signature, Point::Signature>;

    fn structure() -> Self::Signature {
        debug_assert_eq!(Point::structure().size(), Natural::from(6usize));
        DuadsStructure::new(Point::structure())
    }
}

impl<Elem> TryFrom<FiniteSubsetByOrd<Elem>> for Duad<Elem> {
    type Error = &'static str;

    fn try_from(subset: FiniteSubsetByOrd<Elem>) -> Result<Self, Self::Error> {
        let mut elems = subset.elems.into_iter();
        // the elems of subset should be ordered and distinct already
        let duad = Self {
            points: [
                elems.next().ok_or("subset not big enough")?,
                elems.next().ok_or("subset not big enough")?,
            ],
        };
        if elems.next().is_some() {
            return Err("subset too big");
        }
        Ok(duad)
    }
}

impl<Elem: Clone> From<Duad<Elem>> for FiniteSubsetByOrd<Elem> {
    fn from(duad: Duad<Elem>) -> Self {
        FiniteSubsetByOrd {
            elems: duad.points.to_vec(),
        }
    }
}

impl<Elem: Clone> From<&Duad<Elem>> for FiniteSubsetByOrd<Elem> {
    fn from(duad: &Duad<Elem>) -> Self {
        FiniteSubsetByOrd {
            elems: duad.points.to_vec(),
        }
    }
}

/// The 15-element set of duads on a 6-element set
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DuadsStructure<
    Set: ConstSizeFiniteSetSignature<6> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
> {
    _set: PhantomData<Set>,
    set: SetB,
}

impl<
    Set: ConstSizeFiniteSetSignature<6> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
> DuadsStructure<Set, SetB>
{
    pub fn new(set: SetB) -> Self {
        debug_assert_eq!(set.borrow().size(), Natural::from(6usize));
        Self {
            _set: PhantomData,
            set,
        }
    }

    pub fn set(&self) -> &Set {
        self.set.borrow()
    }
}

pub trait SetToDuadsSignature:
    ConstSizeFiniteSetSignature<6> + EnumeratedOrdFiniteSetSignature
{
    fn duads(&self) -> DuadsStructure<Self, &Self> {
        DuadsStructure::new(self)
    }

    fn into_duads(self) -> DuadsStructure<Self, Self> {
        DuadsStructure::new(self)
    }
}
impl<Set: ConstSizeFiniteSetSignature<6> + EnumeratedOrdFiniteSetSignature> SetToDuadsSignature
    for Set
{
}

impl<
    Set: ConstSizeFiniteSetSignature<6> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
> Signature for DuadsStructure<Set, SetB>
{
}

impl<
    Set: ConstSizeFiniteSetSignature<6> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
> SetSignature for DuadsStructure<Set, SetB>
{
    type Elem = Duad<Set::Elem>;

    fn validate_element(&self, duad: &Self::Elem) -> Result<(), String> {
        if !self.set().cmp(&duad.points[0], &duad.points[1]).is_lt() {
            return Err("invalid duad".to_string());
        }
        Ok(())
    }
}

impl<
    Set: ConstSizeFiniteSetSignature<6> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
> EqSignature for DuadsStructure<Set, SetB>
{
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        self.set().finite_subsets().equal(&a.into(), &b.into())
    }
}

impl<
    Set: ConstSizeFiniteSetSignature<6> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
> PartialOrdSignature for DuadsStructure<Set, SetB>
{
    fn partial_cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Ordering> {
        self.set()
            .finite_subsets()
            .partial_cmp(&a.into(), &b.into())
    }
}

impl<
    Set: ConstSizeFiniteSetSignature<6> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
> OrdSignature for DuadsStructure<Set, SetB>
{
    fn cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        self.set().finite_subsets().cmp(&a.into(), &b.into())
    }
}

impl<
    Set: ConstSizeFiniteSetSignature<6> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
> CountableSetSignature for DuadsStructure<Set, SetB>
{
    fn into_generate_all_elements(self) -> impl Iterator<Item = Self::Elem> {
        self.set()
            .clone()
            .into_fixed_size_finite_subsets(2)
            .into_generate_all_elements()
            .map(|subset| Duad::try_from(subset).unwrap())
    }

    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> {
        self.clone().into_generate_all_elements()
    }
}

impl<
    Set: ConstSizeFiniteSetSignature<6> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
> FiniteSetSignature for DuadsStructure<Set, SetB>
{
    fn size(&self) -> Natural {
        debug_assert_eq!(
            self.set().fixed_size_finite_subsets(2).size(),
            Natural::from(15usize)
        );
        Natural::from(15usize)
    }
}

impl<
    Set: ConstSizeFiniteSetSignature<6> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
> ConstSizeFiniteSetSignature<15> for DuadsStructure<Set, SetB>
{
}

impl<
    Set: ConstSizeFiniteSetSignature<6> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
> EnumeratedOrdFiniteSetSignature for DuadsStructure<Set, SetB>
{
    fn list_all_elements_ordered(&self) -> Vec<Self::Elem> {
        self.list_all_elements()
    }

    fn element_to_enumeration(&self, elem: &Self::Elem) -> Natural {
        self.set()
            .fixed_size_finite_subsets(2)
            .element_to_enumeration(&elem.into())
    }

    fn enumeration_to_element(&self, num: &Natural) -> Option<Self::Elem> {
        self.set()
            .fixed_size_finite_subsets(2)
            .enumeration_to_element(num)
            .map(|subset| Duad::try_from(subset).unwrap())
    }
}

pub enum DuadOverlapResult<Set: EnumeratedOrdFiniteSetSignature> {
    Equal,
    Disjoint,
    UniqueCommonPoint(Set::Elem),
}

impl<Set: EnumeratedOrdFiniteSetSignature> DuadOverlapResult<Set> {
    pub fn is_equal(&self) -> bool {
        match self {
            DuadOverlapResult::Equal => true,
            DuadOverlapResult::Disjoint | DuadOverlapResult::UniqueCommonPoint(_) => false,
        }
    }

    pub fn is_disjoint(&self) -> bool {
        match self {
            DuadOverlapResult::Disjoint => true,
            DuadOverlapResult::Equal | DuadOverlapResult::UniqueCommonPoint(_) => false,
        }
    }

    pub fn unwrap_common_point(&self) -> &Set::Elem {
        match self {
            DuadOverlapResult::Equal | DuadOverlapResult::Disjoint => panic!(),
            DuadOverlapResult::UniqueCommonPoint(point) => point,
        }
    }

    pub fn into_unwrap_common_point(self) -> Set::Elem {
        match self {
            DuadOverlapResult::Equal | DuadOverlapResult::Disjoint => panic!(),
            DuadOverlapResult::UniqueCommonPoint(point) => point,
        }
    }
}

impl<
    Set: ConstSizeFiniteSetSignature<6> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
> DuadsStructure<Set, SetB>
{
    pub fn duad(
        &self,
        point_1: Set::Elem,
        point_2: Set::Elem,
    ) -> Result<Duad<Set::Elem>, &'static str> {
        match self.set().cmp(&point_1, &point_2) {
            Ordering::Equal => Err("points are not distinct"),
            Ordering::Less => {
                let duad = Duad {
                    points: [point_1, point_2],
                };
                debug_assert!(self.is_element(&duad));
                Ok(duad)
            }
            Ordering::Greater => {
                let duad = Duad {
                    points: [point_2, point_1],
                };
                debug_assert!(self.is_element(&duad));
                Ok(duad)
            }
        }
    }

    pub fn overlap(&self, d1: &Duad<Set::Elem>, d2: &Duad<Set::Elem>) -> DuadOverlapResult<Set> {
        let mut common = vec![];
        for item in self.set().merge_sorted_and_unique(
            vec![&d1.points[0], &d1.points[1]],
            vec![&d2.points[0], &d2.points[1]],
        ) {
            match item {
                MergedUniqueSource::First(_) | MergedUniqueSource::Second(_) => {}
                MergedUniqueSource::Both(p1, p2) => {
                    debug_assert!(self.set().equal(p1, p2));
                    common.push(p1);
                }
            }
        }
        match common.len() {
            0 => DuadOverlapResult::Disjoint,
            1 => DuadOverlapResult::UniqueCommonPoint(common.into_iter().next().unwrap().clone()),
            2 => DuadOverlapResult::Equal,
            _ => {
                unreachable!()
            }
        }
    }

    /// Convert a duad into a 2-swap
    pub fn to_permutation(
        &self,
        duad: &Duad<Set::Elem>,
    ) -> FinitelySupportedPermutation<Set::Elem> {
        self.set()
            .permutations()
            .new_swap(duad.points[0].clone(), duad.points[1].clone())
            .unwrap()
    }

    /// Convert a 2-swap into a duad, or return None if the permutation is not a 2-swap
    pub fn try_from_permutation(
        &self,
        swap: &FinitelySupportedPermutation<Set::Elem>,
    ) -> Option<Duad<Set::Elem>> {
        let disjoint_cycles = self.set().permutations().disjoint_cycles(swap);
        if disjoint_cycles.len() != 1 {
            return None;
        }
        let unique_cycle = disjoint_cycles.into_iter().next().unwrap();
        if unique_cycle.len() != 2 {
            return None;
        }
        let mut unique_cycle = unique_cycle.into_iter();
        Some(
            self.duad(unique_cycle.next().unwrap(), unique_cycle.next().unwrap())
                .unwrap(),
        )
    }
}

#[signature_meta_trait]
pub trait SetPermutationAsDuadPermutation<
    Set: ConstSizeFiniteSetSignature<6> + EnumeratedOrdFiniteSetSignature,
>: PermutationsSignature<Set>
{
    fn duad_image(&self, set_perm: &Self::Elem, duad: &Duad<Set::Elem>) -> Duad<Set::Elem> {
        let set = self.set();
        debug_assert_eq!(set.size(), Natural::from(6usize));
        let duads = set.duads();
        duads
            .duad(
                self.image(set_perm, &duad.points[0]),
                self.image(set_perm, &duad.points[1]),
            )
            .unwrap()
    }

    fn duad_action(&self, set_perm: &Self::Elem) -> FinitelySupportedPermutation<Duad<Set::Elem>> {
        let set = self.set();
        debug_assert_eq!(set.size(), Natural::from(6usize));
        let duads = set.duads();
        let duad_perms = duads.permutations();
        duad_perms
            .new_perm(
                duads
                    .list_all_elements()
                    .into_iter()
                    .map(|from| {
                        let to = self.duad_image(set_perm, &from);
                        (from, to)
                    })
                    .collect(),
            )
            .unwrap()
    }
}
impl<
    Set: ConstSizeFiniteSetSignature<6> + EnumeratedOrdFiniteSetSignature,
    SetPerms: PermutationsSignature<Set>,
> SetPermutationAsDuadPermutation<Set> for SetPerms
{
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_sets::sets::{
        FiniteSetToFinitelySupportedPermutationsStructure, SetToConstSizeFiniteSubsetByOrdSignature,
    };
    use algebraeon_structures::MetaType;
    use std::collections::HashMap;

    #[test]
    fn test_enumeration() {
        algebraeon_structures::assert_enumerated_ord_finite_set!(
            i32::structure()
                .into_const_size_finite_subset([1, 2, 3, 4, 5, 6])
                .into_duads(),
            15
        );
    }

    #[test]
    fn test_permutation() {
        let set = i32::structure().into_const_size_finite_subset([1, 2, 3, 4, 5, 6]);
        let set_perms = set.permutations();
        let duads = set.duads();
        let duad_perms = duads.permutations();

        assert_eq!(
            duad_perms
                .cycle_shape(&set_perms.duad_action(&set_perms.new_cycle(vec![1, 2, 3]).unwrap())),
            HashMap::from([(3, 4)])
        );

        assert_eq!(
            set_perms.cycle_shape(&duads.to_permutation(&duads.duad(1, 2).unwrap())),
            HashMap::from([(2, 1)])
        );

        assert!(
            duads.equal(
                &duads
                    .try_from_permutation(&set_perms.new_swap(2, 3).unwrap())
                    .unwrap(),
                &duads.duad(2, 3).unwrap()
            )
        );

        assert!(
            &duads
                .try_from_permutation(&set_perms.new_cycle(vec![2, 3, 4]).unwrap())
                .is_none()
        );
    }
}
