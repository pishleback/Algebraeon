//! The hexacode on the standard ordered syntheme

use crate::linear_codes::ordered_syntheme::OrderedSynthemePoint;
use algebraeon_macros::CanonicalStructure;
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
    structure::{FinitelyFreeModuleSignature, MetaAdditionSignature, MetaZeroEqSignature},
};
use algebraeon_sets::sets::Function;
use algebraeon_structures::*;
use cantor::Finite;
use std::{
    ops::{Add, BitAnd, BitOr},
    sync::{Arc, OnceLock},
};

type F2 = Modulo<2>;
type F2Structure = <F2 as MetaType>::Signature;
const ZERO: F2 = F2::new(0);
const ONE: F2 = F2::new(1);

type F4 = QuaternaryField;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Finite, CanonicalStructure)]
#[canonical_structure(eq, partial_ord, ord, finite, ord_finite)]
pub struct Point {
    row: F4,
    col: OrderedSynthemePoint,
}

impl ConstSizeFiniteSetSignature<24> for PointCanonicalStructure {}

pub type LabelledPoints<Elem> = Function<24, Point, Elem>;

#[derive(Clone, PartialEq, Eq)]
pub struct Vector(LabelledPoints<F2>);

impl std::fmt::Debug for Vector {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use std::fmt::Write;
        let mut s = String::new();
        for (i, row) in F4::list_all_elements_ordered().into_iter().enumerate() {
            if i != 0 {
                write!(&mut s, " ")?;
            }
            for col in OrderedSynthemePoint::list_all_elements_ordered() {
                let p = Point { row, col };
                if self.at(&p).is_zero() {
                    write!(&mut s, "0")?;
                } else {
                    write!(&mut s, "1")?;
                }
            }
        }
        f.debug_tuple("Vector")
            .field(&format_args!("{}", s))
            .finish()
    }
}

impl From<LabelledPoints<F2>> for Vector {
    fn from(value: LabelledPoints<F2>) -> Self {
        Self(value)
    }
}

impl From<Vector> for LabelledPoints<F2> {
    fn from(value: Vector) -> Self {
        value.0
    }
}

impl<'a> From<&'a Vector> for &'a LabelledPoints<F2> {
    fn from(value: &'a Vector) -> Self {
        &value.0
    }
}

impl From<[F2; 24]> for Vector {
    fn from(value: [F2; 24]) -> Self {
        LabelledPoints::from(value).into()
    }
}

impl Vector {
    pub fn new(f: impl FnMut(Point) -> F2) -> Vector {
        Vector(LabelledPoints::new(f))
    }

    pub fn at(&self, point: &Point) -> &F2 {
        &self.0[point.element_to_enumeration().try_into().unwrap()]
    }

    pub fn at_mut(&mut self, point: &Point) -> &mut F2 {
        &mut self.0[point.element_to_enumeration().try_into().unwrap()]
    }
}

impl Add<&Vector> for &Vector {
    type Output = Vector;

    fn add(self, other: &Vector) -> Self::Output {
        LabelledPoints::<F2>::new(|p| F2::add(self.0.image(&p), other.0.image(&p))).into()
    }
}

impl BitAnd<&Vector> for &Vector {
    type Output = Vector;

    fn bitand(self, other: &Vector) -> Self::Output {
        LabelledPoints::<F2>::new(|p| {
            match (self.0.image(&p).is_zero(), other.0.image(&p).is_zero()) {
                (true, true) | (true, false) | (false, true) => ZERO,
                (false, false) => ONE,
            }
        })
        .into()
    }
}

impl BitOr<&Vector> for &Vector {
    type Output = Vector;

    fn bitor(self, other: &Vector) -> Self::Output {
        LabelledPoints::<F2>::new(|p| {
            match (self.0.image(&p).is_zero(), other.0.image(&p).is_zero()) {
                (true, true) => ZERO,
                (false, false) | (true, false) | (false, true) => ONE,
            }
        })
        .into()
    }
}

impl Vector {
    pub fn is_codeword(&self) -> bool {
        extended_binary_golay_code_subspace_structure().is_element(self.into())
    }

    pub fn weight(&self) -> usize {
        self.0.iter().map(|x| if *x == ZERO { 0 } else { 1 }).sum()
    }

    pub fn is_octad(&self) -> bool {
        self.is_codeword() && self.weight() == 8
    }

    pub fn is_foursome(&self) -> bool {
        self.weight() == 4
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OrderedSextet(LabelledPoints<OrderedSynthemePoint>);

impl From<LabelledPoints<OrderedSynthemePoint>> for OrderedSextet {
    fn from(value: LabelledPoints<OrderedSynthemePoint>) -> Self {
        Self(value)
    }
}

impl From<OrderedSextet> for LabelledPoints<OrderedSynthemePoint> {
    fn from(value: OrderedSextet) -> Self {
        value.0
    }
}

impl<'a> From<&'a OrderedSextet> for &'a LabelledPoints<OrderedSynthemePoint> {
    fn from(value: &'a OrderedSextet) -> Self {
        &value.0
    }
}

type AmbientSpace = ConstFinitelyFreeModuleStructure<24, PointCanonicalStructure, F2Structure>;

struct ExtendedBinaryGolayCodeCache {
    subspace:
        Arc<FinitelyFreeSubmoduleStructure<PointCanonicalStructure, F2Structure, AmbientSpace>>,
}

static EXTENDED_BINARY_GOLAY_CODE_CACHE: OnceLock<ExtendedBinaryGolayCodeCache> = OnceLock::new();

fn cache() -> &'static ExtendedBinaryGolayCodeCache {
    EXTENDED_BINARY_GOLAY_CODE_CACHE.get_or_init(|| {
        let points = Point::structure();
        let space = F2::structure().free_module(&points);
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
pub fn extended_binary_golay_code_subspace_structure()
-> Arc<FinitelyFreeSubmoduleStructure<PointCanonicalStructure, F2Structure, AmbientSpace>> {
    cache().subspace.clone()
}

pub fn all_codewords() -> Vec<Vector> {
    extended_binary_golay_code_subspace_structure()
        .list_all_elements()
        .into_iter()
        .map(Vector::from)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_sets::sets::SetToConstSizePermutationsStructure;

    #[test]
    fn weight_distribution() {
        let mut wt0 = 0;
        let mut wt8 = 0;
        let mut wt12 = 0;
        let mut wt16 = 0;
        let mut wt24 = 0;
        for x in all_codewords() {
            assert_eq!(x.is_octad(), x.weight() == 8);
            match x.weight() {
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
    fn vector_ops() {
        const Z: F2 = ZERO;
        const O: F2 = ONE;

        let a: Vector = [
            O, O, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z,
        ]
        .into();

        let b: Vector = [
            O, Z, O, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z,
        ]
        .into();

        assert_eq!(
            &a + &b,
            [
                Z, O, O, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z,
            ]
            .into()
        );

        assert_eq!(
            &a & &b,
            [
                O, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z,
            ]
            .into()
        );

        assert_eq!(
            &a | &b,
            [
                O, O, O, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z,
            ]
            .into()
        );
    }

    #[test]
    fn domain_permutation_action() {
        const Z: F2 = ZERO;
        const O: F2 = ONE;

        let a = LabelledPoints::<F2>::from([
            O, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z,
        ]);

        let b = LabelledPoints::<F2>::from([
            Z, O, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z,
        ]);

        let c = LabelledPoints::<F2>::from([
            Z, Z, O, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z,
        ]);

        let cycle012 = Point::structure()
            .const_size_permutations()
            .new_cycle(vec![
                Point::structure()
                    .enumeration_to_element(&Natural::from(0u32))
                    .unwrap(),
                Point::structure()
                    .enumeration_to_element(&Natural::from(1u32))
                    .unwrap(),
                Point::structure()
                    .enumeration_to_element(&Natural::from(2u32))
                    .unwrap(),
            ])
            .unwrap();

        // The right action given by precomposition with domain elements moves labels according to the inverse permutation
        assert_eq!(
            LabelledPoints::<F2>::structure()
                .domain_precomposition_const_size_permutation_action()
                .apply(&cycle012, &a),
            c
        );

        // The left action given by precomposition of the inverse with domain elements moves labels according to the permutation
        assert_eq!(
            LabelledPoints::<F2>::structure()
                .output_const_size_permutation_action()
                .apply(&cycle012, &a),
            b
        );
    }
}
