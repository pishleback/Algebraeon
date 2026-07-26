use algebraeon_rings::{
    linear::{
        finitely_free_module::{FinitelyFreeModuleStructure, RingToFinitelyFreeModuleSignature},
        finitely_free_submodule::FinitelyFreeSubmoduleStructure,
    },
    num_theory::modulo::const_naive::{Modulo, ModuloCanonicalStructure},
    structure::{
        AdditionSignature, AdditiveGroupSignature, AdditiveMonoidSignature,
        CancellativeAdditionSignature, RinglikeSpecializationSignature, SemiModuleSignature,
        TryNegateSignature, ZeroSignature,
    },
};
use algebraeon_structures::*;
use std::cmp::Ordering;

/// The extended binary Golay code on a 24-element set {0, 1, ..., 22, 23 (inf)}
///
///   0    1      2    3      4    5
/// +----+----+ +----+----+ +----+----+
/// | ∞  | 0  | | 1  | 11 | | 2  | 22 |  0
/// +----+----+ +----+----+ +----+----+
/// | 19 | 3  | | 20 | 4  | | 10 | 18 |  1
/// +----+----+ +----+----+ +----+----+
/// | 15 | 6  | | 14 | 16 | | 17 | 8  |  a
/// +----+----+ +----+----+ +----+----+
/// | 5  | 9  | | 21 | 13 | | 7  | 12 |  b
/// +----+----+ +----+----+ +----+----+
///
/// The numbering fairly arbitrary, but is chosen such that the group PSL(2, F32) acting on the points is a subgroup of M24
/// The columns are ordered 0..6 and the hexacode on the columns is related to the structure of the extended binary Golay code
/// The rows are labelled by the elements of F4
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExtendedBinaryGolayCodeStructure<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<24>,
    SetB: BorrowedStructure<Set>,
> {
    ebgc: FinitelyFreeSubmoduleStructure<
        Set,
        SetB,
        ModuloCanonicalStructure<2>,
        ModuloCanonicalStructure<2>,
        FinitelyFreeModuleStructure<
            Set,
            SetB,
            ModuloCanonicalStructure<2>,
            ModuloCanonicalStructure<2>,
        >,
    >,
}

pub trait SetToExtendedBinaryGolayCodeSignature:
    EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<24>
{
    fn extended_binary_golay_code(&self) -> ExtendedBinaryGolayCodeStructure<Self, &Self> {
        ExtendedBinaryGolayCodeStructure::new(self)
    }

    fn into_extended_binary_golay_code(self) -> ExtendedBinaryGolayCodeStructure<Self, Self> {
        ExtendedBinaryGolayCodeStructure::new(self)
    }
}
impl<Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<24>>
    SetToExtendedBinaryGolayCodeSignature for Set
{
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<24>,
    SetB: BorrowedStructure<Set>,
> ExtendedBinaryGolayCodeStructure<Set, SetB>
{
    pub fn new(set: SetB) -> Self {
        let full_space = Modulo::structure().into_free_module(set);
        let ebgc_subspace = full_space.generated_submodule(vec![todo!()]);
        Self {
            ebgc: FinitelyFreeSubmoduleStructure::new(full_space, ebgc_subspace),
        }
    }

    pub fn set(&self) -> &Set {
        self.ebgc.set()
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<24>,
    SetB: BorrowedStructure<Set>,
> Signature for ExtendedBinaryGolayCodeStructure<Set, SetB>
{
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<24>,
    SetB: BorrowedStructure<Set>,
> SetSignature for ExtendedBinaryGolayCodeStructure<Set, SetB>
{
    type Elem = Vec<Modulo<2>>;

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        self.ebgc.validate_element(x)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<24>,
    SetB: BorrowedStructure<Set>,
> EqSignature for ExtendedBinaryGolayCodeStructure<Set, SetB>
{
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        self.ebgc.equal(a, b)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<24>,
    SetB: BorrowedStructure<Set>,
> PartialOrdSignature for ExtendedBinaryGolayCodeStructure<Set, SetB>
{
    fn partial_cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Ordering> {
        self.ebgc.partial_cmp(a, b)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<24>,
    SetB: BorrowedStructure<Set>,
> OrdSignature for ExtendedBinaryGolayCodeStructure<Set, SetB>
{
    fn cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        self.ebgc.cmp(a, b)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<24>,
    SetB: BorrowedStructure<Set>,
> CountableSetSignature for ExtendedBinaryGolayCodeStructure<Set, SetB>
{
    fn into_generate_all_elements(self) -> impl Iterator<Item = Self::Elem> {
        self.ebgc.into_generate_all_elements()
    }

    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> {
        self.ebgc.generate_all_elements()
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<24>,
    SetB: BorrowedStructure<Set>,
> FiniteSetSignature for ExtendedBinaryGolayCodeStructure<Set, SetB>
{
    fn list_all_elements(&self) -> Vec<Self::Elem> {
        self.ebgc.list_all_elements()
    }

    fn size(&self) -> Natural {
        self.ebgc.size()
    }

    fn generate_random_elements(&self, seed: u64) -> impl Iterator<Item = Self::Elem> {
        self.ebgc.generate_random_elements(seed)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<24>,
    SetB: BorrowedStructure<Set>,
> EnumeratedOrdFiniteSetSignature for ExtendedBinaryGolayCodeStructure<Set, SetB>
{
    fn list_all_elements_ordered(&self) -> Vec<Self::Elem> {
        self.ebgc.list_all_elements_ordered()
    }

    fn element_to_enumeration(&self, elem: &Self::Elem) -> Natural {
        self.ebgc.element_to_enumeration(elem)
    }

    fn enumeration_to_element(&self, num: &Natural) -> Option<Self::Elem> {
        self.ebgc.enumeration_to_element(num)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<24>,
    SetB: BorrowedStructure<Set>,
> RinglikeSpecializationSignature for ExtendedBinaryGolayCodeStructure<Set, SetB>
{
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<24>,
    SetB: BorrowedStructure<Set>,
> ZeroSignature for ExtendedBinaryGolayCodeStructure<Set, SetB>
{
    fn zero(&self) -> Self::Elem {
        self.ebgc.zero()
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<24>,
    SetB: BorrowedStructure<Set>,
> AdditionSignature for ExtendedBinaryGolayCodeStructure<Set, SetB>
{
    fn add(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.ebgc.add(a, b)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<24>,
    SetB: BorrowedStructure<Set>,
> CancellativeAdditionSignature for ExtendedBinaryGolayCodeStructure<Set, SetB>
{
    fn try_sub(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        self.ebgc.try_sub(a, b)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<24>,
    SetB: BorrowedStructure<Set>,
> TryNegateSignature for ExtendedBinaryGolayCodeStructure<Set, SetB>
{
    fn try_neg(&self, a: &Self::Elem) -> Option<Self::Elem> {
        self.ebgc.try_neg(a)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<24>,
    SetB: BorrowedStructure<Set>,
> AdditiveMonoidSignature for ExtendedBinaryGolayCodeStructure<Set, SetB>
{
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<24>,
    SetB: BorrowedStructure<Set>,
> AdditiveGroupSignature for ExtendedBinaryGolayCodeStructure<Set, SetB>
{
    fn neg(&self, a: &Self::Elem) -> Self::Elem {
        self.ebgc.neg(a)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<24>,
    SetB: BorrowedStructure<Set>,
> SemiModuleSignature<ModuloCanonicalStructure<2>> for ExtendedBinaryGolayCodeStructure<Set, SetB>
{
    fn ring(&self) -> &ModuloCanonicalStructure<2> {
        self.ebgc.ring()
    }

    fn scalar_mul(&self, a: &Self::Elem, x: &Modulo<2>) -> Self::Elem {
        self.ebgc.scalar_mul(a, x)
    }
}

#[cfg(test)]
mod tests {
    use std::println;

    use super::*;
    use algebraeon_rings::finite_fields::quaternary_field::QuaternaryField;
    use algebraeon_sets::sets::{
        CartesianProductSetStructure, FiniteSetToFiniteSetSizedSignature,
        SetToFiniteSubsetByOrdSizedSignature,
    };

    #[test]
    fn test() {
        let set = i32::structure()
            .into_finite_subset_sized(std::array::from_fn::<_, 24, _>(|i| i as i32));
        let ebgc = set.extended_binary_golay_code();

        println!("{:?}", ebgc);

        for x in ebgc.list_all_elements() {
            println!("{:?}", x);
        }
    }

    #[test]
    fn test_enumeration() {
        let set = i32::structure()
            .into_finite_subset_sized(std::array::from_fn::<_, 24, _>(|i| i as i32));
        let ebgc = set.extended_binary_golay_code();
        let codewords = ebgc.list_all_elements_ordered();
        assert_eq!(codewords.len(), 2usize.pow(12));
        assert_eq!(ebgc.size(), Natural::from(2usize.pow(12)));

        // codewords are all valid
        for v in &codewords {
            println!("{:?}", v);
            assert!(ebgc.validate_element(v).is_ok());
        }

        // enumeration is correct
        for (i, v) in codewords.iter().enumerate() {
            assert_eq!(Natural::from(i), ebgc.element_to_enumeration(v));
            assert!(ebgc.equal(&ebgc.enumeration_to_element(&Natural::from(i)).unwrap(), v));
        }
        assert!(
            ebgc.enumeration_to_element(&Natural::from(2usize.pow(12)))
                .is_none()
        );
    }

    #[test]
    fn test_6set_times_quaternary_field() {
        let set6 = i32::structure().into_finite_subset_sized([1, 2, 3, 4, 5, 6]);

        let set24 = CartesianProductSetStructure::new(QuaternaryField::structure(), set6);

        for x in set24.list_all_elements_ordered() {
            println!("{:?}", x);
        }

        let ebgc = set24
            .try_into_sized::<24>()
            .unwrap()
            .extended_binary_golay_code();
    }
}
