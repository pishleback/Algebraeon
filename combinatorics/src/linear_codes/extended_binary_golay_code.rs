//! The hexacode on the standard ordered syntheme

use crate::linear_codes::ordered_syntheme::{
    OrderedSynthemePoint, OrderedSynthemePointCanonicalStructure,
};
use algebraeon_rings::{
    finite_fields::quaternary_field::QuaternaryField,
    linear::{
        const_finitely_free_module::{
            ConstFinitelyFreeModuleStructure, RingToConstFinitelyFreeModuleSignature,
        },
        finitely_free_submodule::FinitelyFreeSubmoduleStructure,
        finitely_free_submodules::FinitelyFreeSubmodule,
    },
    num_theory::modulo::const_naive::Modulo,
    structure::FinitelyFreeModuleSignature,
};
use algebraeon_sets::sets::{CartesianProductSetStructure, Function};
use algebraeon_structures::*;
use ambassador::Delegate;
use std::sync::OnceLock;

type F2 = Modulo<2>;
type F2Structure = <F2 as MetaType>::Signature;
const ZERO: F2 = F2::new(0);
const ONE: F2 = F2::new(1);

type F4 = QuaternaryField;
type F4Structure = <F4 as MetaType>::Signature;

#[derive(Delegate, Debug, Clone, PartialEq, Eq)]
#[delegate(Signature)]
#[delegate(SetSignature)]
#[delegate(EqSignature)]
#[delegate(PartialOrdSignature)]
#[delegate(OrdSignature)]
#[delegate(CountableSetSignature)]
#[delegate(FiniteSetSignature)]
#[delegate(EnumeratedOrdFiniteSetSignature)]
pub struct MogPointStructure(
    CartesianProductSetStructure<
        F4Structure,
        F4Structure,
        OrderedSynthemePointCanonicalStructure,
        OrderedSynthemePointCanonicalStructure,
    >,
);
impl Default for MogPointStructure {
    fn default() -> Self {
        Self(CartesianProductSetStructure::new(
            F4::structure(),
            OrderedSynthemePoint::structure(),
        ))
    }
}
impl ConstSizeFiniteSetSignature<24> for MogPointStructure {}
type MogPoint = <MogPointStructure as SetSignature>::Elem;

pub type LabelledPoints<Elem> = Function<24, MogPoint, Elem>;
pub type Vector = LabelledPoints<F2>;

type AmbientSpace = ConstFinitelyFreeModuleStructure<
    24,
    MogPointStructure,
    MogPointStructure,
    F2Structure,
    F2Structure,
>;

struct ExtendedBinaryGolayCodeCache {
    subspace:
        FinitelyFreeSubmoduleStructure<MogPointStructure, F2Structure, AmbientSpace, AmbientSpace>,
}

static EXTENDED_BINARY_GOLAY_CODE_CACHE: OnceLock<ExtendedBinaryGolayCodeCache> = OnceLock::new();

fn cache() -> &'static ExtendedBinaryGolayCodeCache {
    EXTENDED_BINARY_GOLAY_CODE_CACHE.get_or_init(|| {
        let points = MogPointStructure::default();
        let space = F2::structure().into_free_module(points);
        const Z: F2 = ZERO;
        const O: F2 = ONE;
        #[rustfmt::skip]
        let subspace = space.generated_submodule(vec![
           &[
            O, O, Z, Z, Z, Z,
            O, O, Z, Z, Z, Z,
            O, O, Z, Z, Z, Z,
            O, O, Z, Z, Z, Z,
        ].into(), &[
            O, Z, O, Z, Z, Z,
            O, Z, O, Z, Z, Z,
            O, Z, O, Z, Z, Z,
            O, Z, O, Z, Z, Z,
        ].into(), &[
            O, Z, Z, O, Z, Z,
            O, Z, Z, O, Z, Z,
            O, Z, Z, O, Z, Z,
            O, Z, Z, O, Z, Z,
        ].into(), &[
            O, Z, Z, Z, O, Z,
            O, Z, Z, Z, O, Z,
            O, Z, Z, Z, O, Z,
            O, Z, Z, Z, O, Z,
        ].into(), &[
            O, Z, Z, Z, Z, O,
            O, Z, Z, Z, Z, O,
            O, Z, Z, Z, Z, O,
            O, Z, Z, Z, Z, O,
        ].into(), &[
            Z, O, Z, Z, Z, Z,
            O, Z, O, O, O, O,
            O, Z, Z, Z, Z, Z,
            O, Z, Z, Z, Z, Z,
        ].into(), &[
            Z, O, Z, Z, Z, Z,
            O, Z, Z, Z, Z, Z,
            O, Z, O, O, O, O,
            O, Z, Z, Z, Z, Z,
        ].into(), &[
            Z, O, Z, Z, Z, Z,
            O, Z, Z, Z, Z, Z,
            O, Z, Z, Z, Z, Z,
            O, Z, O, O, O, O,
        ].into(), &[
            Z, Z, O, Z, Z, Z,
            O, O, Z, O, Z, Z,
            O, Z, Z, Z, O, Z,
            O, Z, Z, Z, Z, O,
        ].into(), &[
            Z, Z, O, Z, Z, Z,
            O, Z, Z, Z, Z, O,
            O, O, Z, O, Z, Z,
            O, Z, Z, Z, O, Z,
        ].into(), &[
            Z, Z, O, Z, Z, Z,
            O, O, Z, O, Z, Z,
            Z, O, Z, Z, Z, O,
            Z, O, Z, Z, O, Z,
        ].into(), &[
            Z, Z, O, Z, Z, Z,
            Z, O, Z, Z, O, Z,
            O, O, Z, O, Z, Z,
            Z, O, Z, Z, Z, O,
        ].into()]);
        let subspace_structure = FinitelyFreeSubmoduleStructure::new(space, subspace);
        debug_assert_eq!(subspace_structure.rank(), 12);
        ExtendedBinaryGolayCodeCache {
            subspace: subspace_structure,
        }
    })
}

/// The 24 dimensional vector space structure over F2 with basis given by the points of the MOG
pub fn space_structure() -> &'static AmbientSpace {
    cache().subspace.module()
}

/// The 12 dimensional vector subspace given by the extended binary Golay code
pub fn extended_binary_golay_code_subspace() -> &'static FinitelyFreeSubmodule<F2> {
    cache().subspace.submodule()
}

/// The 12 dimensional vector subspace structure given by the extended binary Golay code
pub fn extended_binary_golay_code_subspace_structure() -> &'static FinitelyFreeSubmoduleStructure<
    MogPointStructure,
    F2Structure,
    AmbientSpace,
    AmbientSpace,
> {
    &cache().subspace
}

pub fn is_codeword(vector: &LabelledPoints<F2>) -> bool {
    extended_binary_golay_code_subspace_structure().is_element(vector)
}

pub fn weight(vector: &LabelledPoints<F2>) -> usize {
    vector.iter().map(|x| if *x == ZERO { 0 } else { 1 }).sum()
}

pub fn is_octad(vector: &LabelledPoints<F2>) -> bool {
    is_codeword(vector) && weight(vector) == 8
}

#[cfg(test)]
mod tests {
    use crate::linear_codes::{
        hexacode::HexacodeVector,
        ordered_syntheme::{OrderedSynthemePair, OrderedSynthemeSide},
    };

    use super::*;

    #[test]
    fn weight_distribution() {
        let mut wt0 = 0;
        let mut wt8 = 0;
        let mut wt12 = 0;
        let mut wt16 = 0;
        let mut wt24 = 0;
        for x in extended_binary_golay_code_subspace_structure().list_all_elements() {
            assert_eq!(is_octad(&x), weight(&x) == 8);
            match weight(&x) {
                0 => {
                    wt0 += 1;
                }
                8 => {
                    wt8 += 1;
                }
                12 => {
                    wt12 += 1;
                }
                16 => {
                    wt16 += 1;
                }
                24 => {
                    wt24 += 1;
                }
                _ => unreachable!(),
            }
        }
        assert_eq!(wt0, 1);
        assert_eq!(wt8, 759);
        assert_eq!(wt12, 2576);
        assert_eq!(wt16, 759);
        assert_eq!(wt24, 1);
        // total gives all 2^12 codewords
        assert_eq!(wt0 + wt8 + wt12 + wt16 + wt24, 1usize << 12);
    }

    #[test]
    fn test() {
        let v = HexacodeVector::new(|p| match (p.pair, p.side) {
            (OrderedSynthemePair::Left, OrderedSynthemeSide::Left) => F4::Zero,
            (OrderedSynthemePair::Left, OrderedSynthemeSide::Right) => F4::One,
            (OrderedSynthemePair::Middle, OrderedSynthemeSide::Left) => F4::Alpha,
            (OrderedSynthemePair::Middle, OrderedSynthemeSide::Right) => F4::Beta,
            (OrderedSynthemePair::Right, OrderedSynthemeSide::Left) => F4::One,
            (OrderedSynthemePair::Right, OrderedSynthemeSide::Right) => F4::Zero,
        });
        assert_eq!(
            v,
            [F4::Zero, F4::One, F4::Alpha, F4::Beta, F4::One, F4::Zero].into()
        );
    }
}
