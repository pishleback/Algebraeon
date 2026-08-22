//! The extended binary golay code on the cartesian product of the standard ordered syntheme and the finite field of 4 elements

use crate::golay_codes::{
    hexacode::{self, HexacodeVector},
    ordered_syntheme::{OrderedSynthemePair, OrderedSynthemePoint, OrderedSynthemeSide},
};
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
    matrix::Matrix,
    num_theory::modulo::const_naive::Modulo,
    structure::{
        FinitelyFreeModuleSignature, MetaAdditionSignature, MetaOneSignature,
        MetaTryReciprocalSignature, MetaZeroEqSignature, MetaZeroSignature, ZeroSignature,
    },
};
use algebraeon_sets::sets::{
    ConstSizePermutation, Function, SetToConstSizeFunctionsToSignature,
    SetToConstSizePermutationAction, SetToConstSizePermutationsStructure,
};
use algebraeon_structures::*;
use cantor::Finite;
use std::{
    borrow::Borrow,
    collections::{BTreeMap, HashSet},
    ops::{Add, BitAnd, BitOr},
    sync::{Arc, OnceLock},
};

type F2 = Modulo<2>;
type F2Structure = <F2 as MetaType>::Signature;
const ZERO: F2 = F2::new(0);
const ONE: F2 = F2::new(1);

type F4 = QuaternaryField;

// This numbering is chosen such that the group PSL(2, F32) acting on the points is a subgroup of M24
//   0    1      2    3      4    5
// +----+----+ +----+----+ +----+----+
// | ∞  | 0  | | 1  | 11 | | 2  | 22 |  0
// +----+----+ +----+----+ +----+----+
// | 19 | 3  | | 20 | 4  | | 10 | 18 |  1
// +----+----+ +----+----+ +----+----+
// | 15 | 6  | | 14 | 16 | | 17 | 8  |  a
// +----+----+ +----+----+ +----+----+
// | 5  | 9  | | 21 | 13 | | 7  | 12 |  b
// +----+----+ +----+----+ +----+----+

/// This is the extended binary Golay code on a 24-element set.
/// With the 24-element set numbered as shown, the extended binary Golay code is the set of all vectors over F2 such that
///  - The parity of the top row is equal to the parity of every column
///  - The length 6 vector over F4, formed by summing for each column the elements of F4 for which the entry is 1, is a hexacodeword
///
///   0    1      2    3      4    5
/// +----+----+ +----+----+ +----+----+
/// | 0  | 1  | | 2  | 3  | | 4  | 5  |  0
/// +----+----+ +----+----+ +----+----+
/// | 6  | 7  | | 8  | 9  | | 10 | 11 |  1
/// +----+----+ +----+----+ +----+----+
/// | 12 | 13 | | 14 | 15 | | 16 | 17 |  a
/// +----+----+ +----+----+ +----+----+
/// | 18 | 19 | | 20 | 21 | | 22 | 23 |  b
/// +----+----+ +----+----+ +----+----+
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Finite, CanonicalStructure)]
#[canonical_structure(eq, partial_ord, ord, finite, ord_finite)]
pub struct Point {
    pub row: F4,
    pub col: OrderedSynthemePoint,
}

impl ConstSizeFiniteSetSignature<24> for PointCanonicalStructure {}

pub type LabelledPoints<Elem> = Function<24, Point, Elem>;

#[derive(CanonicalStructure, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Vector {
    pub coords: LabelledPoints<F2>,
}

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
    fn from(coords: LabelledPoints<F2>) -> Self {
        Self { coords }
    }
}

impl From<Vector> for LabelledPoints<F2> {
    fn from(value: Vector) -> Self {
        value.coords
    }
}

impl<'a> From<&'a Vector> for &'a LabelledPoints<F2> {
    fn from(value: &'a Vector) -> Self {
        &value.coords
    }
}

impl From<[F2; 24]> for Vector {
    fn from(value: [F2; 24]) -> Self {
        LabelledPoints::from(value).into()
    }
}

impl Vector {
    pub fn from_fn(f: impl FnMut(Point) -> F2) -> Vector {
        Vector {
            coords: LabelledPoints::new(f),
        }
    }

    pub fn at(&self, point: &Point) -> &F2 {
        &self.coords[point.element_to_enumeration().try_into().unwrap()]
    }

    pub fn at_mut(&mut self, point: &Point) -> &mut F2 {
        &mut self.coords[point.element_to_enumeration().try_into().unwrap()]
    }
}

impl Add<&Vector> for &Vector {
    type Output = Vector;

    fn add(self, other: &Vector) -> Self::Output {
        LabelledPoints::<F2>::new(|p| F2::add(self.coords.image(&p), other.coords.image(&p))).into()
    }
}

impl Add<Vector> for &Vector {
    type Output = Vector;

    fn add(self, other: Vector) -> Self::Output {
        self.add(&other)
    }
}

impl Add<&Vector> for Vector {
    type Output = Vector;

    fn add(self, other: &Vector) -> Self::Output {
        (&self).add(other)
    }
}

impl Add<Vector> for Vector {
    type Output = Vector;

    fn add(self, other: Vector) -> Self::Output {
        (&self).add(&other)
    }
}

impl BitAnd<&Vector> for &Vector {
    type Output = Vector;

    fn bitand(self, other: &Vector) -> Self::Output {
        LabelledPoints::<F2>::new(|p| {
            match (
                self.coords.image(&p).is_zero(),
                other.coords.image(&p).is_zero(),
            ) {
                (true, true) | (true, false) | (false, true) => ZERO,
                (false, false) => ONE,
            }
        })
        .into()
    }
}

impl BitAnd<Vector> for &Vector {
    type Output = Vector;

    fn bitand(self, other: Vector) -> Self::Output {
        self.bitand(&other)
    }
}

impl BitAnd<&Vector> for Vector {
    type Output = Vector;

    fn bitand(self, other: &Vector) -> Self::Output {
        (&self).bitand(other)
    }
}

impl BitAnd<Vector> for Vector {
    type Output = Vector;

    fn bitand(self, other: Vector) -> Self::Output {
        (&self).bitand(&other)
    }
}

impl BitOr<&Vector> for &Vector {
    type Output = Vector;

    fn bitor(self, other: &Vector) -> Self::Output {
        LabelledPoints::<F2>::new(|p| {
            match (
                self.coords.image(&p).is_zero(),
                other.coords.image(&p).is_zero(),
            ) {
                (true, true) => ZERO,
                (false, false) | (true, false) | (false, true) => ONE,
            }
        })
        .into()
    }
}

impl BitOr<Vector> for &Vector {
    type Output = Vector;

    fn bitor(self, other: Vector) -> Self::Output {
        self.bitor(&other)
    }
}

impl BitOr<&Vector> for Vector {
    type Output = Vector;

    fn bitor(self, other: &Vector) -> Self::Output {
        (&self).bitor(other)
    }
}

impl BitOr<Vector> for Vector {
    type Output = Vector;

    fn bitor(self, other: Vector) -> Self::Output {
        (&self).bitor(&other)
    }
}

impl Vector {
    pub fn zero() -> Self {
        Self {
            coords: ebgc_structure().zero(),
        }
    }

    pub fn is_codeword(&self) -> bool {
        ebgc_structure().is_element(self.into())
    }

    pub fn weight(&self) -> usize {
        self.coords
            .iter()
            .map(|(_, x)| if *x == ZERO { 0 } else { 1 })
            .sum()
    }

    pub fn is_octad(&self) -> bool {
        self.is_codeword() && self.weight() == 8
    }

    pub fn is_foursome(&self) -> bool {
        self.weight() == 4
    }

    pub fn contains_point(&self, p: &Point) -> bool {
        self.coords.image(p).equal(&F2::one())
    }

    pub fn points(&self) -> impl Iterator<Item = Point> {
        Point::generate_all_elements().filter(|p| self.contains_point(p))
    }

    pub fn from_points(points: &Vec<impl Borrow<Point>>) -> Self {
        let mut coords = LabelledPoints::new(|_| F2::zero());
        for p in points {
            *coords.image_mut(p.borrow()) = F2::one();
        }
        Self { coords }
    }

    pub fn to_row(&self) -> Matrix<F2> {
        Matrix::construct(1, 24, |r, c| {
            debug_assert_eq!(r, 0);
            if self.contains_point(&Point::enumeration_to_element(&c.into()).unwrap()) {
                F2::one()
            } else {
                F2::zero()
            }
        })
    }

    pub fn from_row(m: &Matrix<F2>) -> Option<Self> {
        if m.rows() == 1 && m.cols() == 24 {
            Some(Self::from_fn(|p| {
                let i: usize = p.element_to_enumeration().try_into().unwrap();
                m.at(0, i).unwrap().clone()
            }))
        } else {
            None
        }
    }

    pub fn to_col(&self) -> Matrix<F2> {
        self.to_row().transpose()
    }

    pub fn from_col(m: &Matrix<F2>) -> Option<Self> {
        Self::from_row(&m.transpose_ref())
    }

    /// As a 4x6 matrix representing an element of the MOG
    pub fn to_mat(&self) -> Matrix<F2> {
        Matrix::construct(4, 6, |r, c| {
            if self.contains_point(&Point::enumeration_to_element(&(c + 6 * r).into()).unwrap()) {
                ONE
            } else {
                ZERO
            }
        })
    }

    /// From a 4x6 matrix representing an element of the MOG
    pub fn from_mat(m: &Matrix<F2>) -> Option<Self> {
        if m.rows() == 4 && m.cols() == 6 {
            Some(Self::from_fn(|p| {
                let i: usize = p.element_to_enumeration().try_into().unwrap();
                let (r, c) = (i / 6, i % 6);
                m.at(r, c).unwrap().clone()
            }))
        } else {
            None
        }
    }

    pub fn permute(&self, perm: &ConstSizePermutation<24, Point>) -> Self {
        Self {
            coords: Point::structure()
                .const_size_functions_to(&F2::structure())
                .output_const_size_permutation_action()
                .apply(perm, &self.coords),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OrderedSextet {
    // 6 foursomes with any pair forming an octact
    inner: LabelledPoints<OrderedSynthemePoint>,
}

impl OrderedSextet {
    pub fn permute(&self, permutation: &ConstSizePermutation<6, OrderedSynthemePoint>) -> Self {
        Self {
            inner: Point::structure()
                .const_size_functions_to(&OrderedSynthemePoint::structure())
                .left_action_from_range_action(
                    OrderedSynthemePoint::structure().const_size_permutation_action(),
                )
                .apply(permutation, &self.inner),
        }
    }

    pub fn foursomes(&self) -> hexacode::LabelledPoints<Vector> {
        let mut foursomes = hexacode::LabelledPoints::new(|_| Vector::zero());
        for p in Point::generate_all_elements() {
            *foursomes.image_mut(self.inner.image(&p)).at_mut(&p) = F2::one();
        }
        foursomes
    }

    pub fn from_foursomes(foursomes: hexacode::LabelledPoints<Vector>) -> Self {
        for foursome in foursomes.images() {
            debug_assert_eq!(foursome.weight(), 4);
        }
        let mut labels = LabelledPoints::new(|_| OrderedSynthemePoint {
            side: OrderedSynthemeSide::Left,
            pair: OrderedSynthemePair::Left,
        });
        for (i, foursome) in foursomes.images().enumerate() {
            let q = OrderedSynthemePoint::enumeration_to_element(&i.into()).unwrap();
            for p in foursome.points() {
                *labels.image_mut(&p) = q;
            }
        }
        let s = OrderedSextet { inner: labels };
        #[cfg(debug_assertions)]
        s.validate().unwrap();
        s
    }

    fn validate(&self) -> Result<(), String> {
        let foursomes = self.foursomes();
        let syntheme_points = OrderedSynthemePoint::list_all_elements();
        for q in &syntheme_points {
            if !foursomes.image(q).clone().is_foursome() {
                return Err("Expected a foursome".into());
            }
        }
        let (root, rest) = syntheme_points.split_at(1);
        let root = &root[0];
        for q in rest {
            if !(foursomes.image(root) | foursomes.image(q)).is_octad() {
                return Err("Expected an octad".into());
            }
        }
        Ok(())
    }
}

impl TryFrom<LabelledPoints<OrderedSynthemePoint>> for OrderedSextet {
    type Error = String;

    fn try_from(inner: LabelledPoints<OrderedSynthemePoint>) -> Result<Self, Self::Error> {
        let s = Self { inner };
        s.validate()?;
        Ok(s)
    }
}

impl From<OrderedSextet> for LabelledPoints<OrderedSynthemePoint> {
    fn from(value: OrderedSextet) -> Self {
        value.inner
    }
}

impl<'a> From<&'a OrderedSextet> for &'a LabelledPoints<OrderedSynthemePoint> {
    fn from(value: &'a OrderedSextet) -> Self {
        &value.inner
    }
}

#[derive(Debug, Clone)]
pub struct Sextet {
    // 6 foursomes with any pair forming an octact
    inner: OrderedSextet,
}

impl Sextet {
    pub fn foursomes(&self) -> [Vector; 6] {
        let foursomes = self.inner.foursomes();
        std::array::from_fn(|i| {
            let q = OrderedSynthemePoint::enumeration_to_element(&i.into()).unwrap();
            foursomes.image(&q).clone()
        })
    }

    pub fn from_foursomes(foursomes: hexacode::LabelledPoints<Vector>) -> Self {
        Self {
            inner: OrderedSextet::from_foursomes(foursomes),
        }
    }

    #[allow(unused)]
    fn validate(&self) -> Result<(), String> {
        self.inner.validate()
    }

    pub fn orderings(&self) -> impl Iterator<Item = OrderedSextet> {
        let root = self.inner.clone();
        OrderedSynthemePoint::structure()
            .const_size_permutations()
            .generate_all_elements()
            .map(move |perm| root.permute(&perm))
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OrderedSextetLabelling {
    sextet: OrderedSextet,
    f4_labels: LabelledPoints<F4>,
}

impl OrderedSextetLabelling {
    pub fn standard_labelling() -> Self {
        let s = Self {
            sextet: OrderedSextet {
                inner: Function::new(|p: Point| p.col),
            },
            f4_labels: Function::new(|p: Point| p.row),
        };
        #[cfg(debug_assertions)]
        s.validate().unwrap();
        s
    }

    pub fn foursomes(&self) -> hexacode::LabelledPoints<Vector> {
        self.sextet.foursomes()
    }

    pub fn point_foursomes(&self) -> &LabelledPoints<OrderedSynthemePoint> {
        &self.sextet.inner
    }

    #[allow(unused)]
    fn validate(&self) -> Result<(), String> {
        self.sextet.validate()?;
        for foursome in self.foursomes() {
            debug_assert_eq!(foursome.weight(), 4);
            let mut counts = [0usize; 4];
            for p in foursome.points() {
                let i: usize = self
                    .f4_labels
                    .image(&p)
                    .element_to_enumeration()
                    .try_into()
                    .unwrap();
                counts[i] += 1;
            }
            for c in counts {
                if c != 1 {
                    return Err("Expected one of each element of F4".to_string());
                }
            }
        }
        Ok(())
    }

    pub fn mog_isomorphism(
        &self,
    ) -> (
        ConstSizePermutation<24, Point>,
        ConstSizePermutation<24, Point>,
    ) {
        let mut to_mog = [0; 24];
        let mut from_mog = [0; 24];
        for (c, foursome) in self.sextet.foursomes().into_iter().enumerate() {
            for p in foursome.points() {
                let r: usize = self
                    .f4_labels
                    .image(&p)
                    .element_to_enumeration()
                    .try_into()
                    .unwrap();
                let p: usize = p.element_to_enumeration().try_into().unwrap();
                let q = c + 6 * r;
                to_mog[p] = q;
                from_mog[q] = p;
            }
        }
        (
            ConstSizePermutation::new_fn(|p: &Point| {
                let i: usize = p.element_to_enumeration().try_into().unwrap();
                Point::enumeration_to_element(&to_mog[i].into()).unwrap()
            })
            .unwrap(),
            ConstSizePermutation::new_fn(|p: &Point| {
                let i: usize = p.element_to_enumeration().try_into().unwrap();
                Point::enumeration_to_element(&from_mog[i].into()).unwrap()
            })
            .unwrap(),
        )
    }

    pub fn permute_foursomes(self, perm: &ConstSizePermutation<6, OrderedSynthemePoint>) -> Self {
        Self {
            sextet: self.sextet.permute(perm),
            f4_labels: self.f4_labels,
        }
    }

    pub fn labels(&self) -> &LabelledPoints<F4> {
        &self.f4_labels
    }

    pub fn add_vector(self, vector: &HexacodeVector) -> Self {
        let point_foursomes = &self.sextet.inner;
        Self {
            f4_labels: LabelledPoints::new(|point: Point| {
                let foursome = point_foursomes.image(&point);
                *self.f4_labels.get(&point) + *vector.at(foursome)
            }),
            sextet: self.sextet,
        }
    }

    pub fn scalar_mul(self, lambda: F4) -> Self {
        assert_ne!(lambda, F4::Zero);
        Self {
            f4_labels: LabelledPoints::new(|point: Point| {
                // use lambda.reciprocal() here because we want to permute the points not the labels
                *self.f4_labels.get(&point) * lambda.try_reciprocal().unwrap()
            }),
            sextet: self.sextet,
        }
    }

    pub fn conjugate(self) -> Self {
        Self {
            f4_labels: LabelledPoints::new(|point: Point| self.f4_labels.get(&point).conjugate()),
            sextet: self.sextet,
        }
    }

    pub fn permutation_to_standard_labelling(&self) -> ConstSizePermutation<24, Point> {
        ConstSizePermutation::new_fn(|p| Point {
            row: *self.f4_labels.image(p),
            col: *self.sextet.inner.image(p),
        })
        .unwrap()
    }
}

type AmbientSpace = ConstFinitelyFreeModuleStructure<24, PointCanonicalStructure, F2Structure>;

struct ExtendedBinaryGolayCodeCache {
    subspace:
        Arc<FinitelyFreeSubmoduleStructure<PointCanonicalStructure, F2Structure, AmbientSpace>>,
    octad_from_5_points: BTreeMap<Vector, Vector>,
}

static EBGC_CACHE: OnceLock<ExtendedBinaryGolayCodeCache> = OnceLock::new();

fn ebgc_basis() -> [Vector; 12] {
    const Z: F2 = ZERO;
    const O: F2 = ONE;
    [
        [
            O, O, Z, Z, Z, Z, O, O, Z, Z, Z, Z, O, O, Z, Z, Z, Z, O, O, Z, Z, Z, Z,
        ],
        [
            O, Z, O, Z, Z, Z, O, Z, O, Z, Z, Z, O, Z, O, Z, Z, Z, O, Z, O, Z, Z, Z,
        ],
        [
            O, Z, Z, O, Z, Z, O, Z, Z, O, Z, Z, O, Z, Z, O, Z, Z, O, Z, Z, O, Z, Z,
        ],
        [
            O, Z, Z, Z, O, Z, O, Z, Z, Z, O, Z, O, Z, Z, Z, O, Z, O, Z, Z, Z, O, Z,
        ],
        [
            O, Z, Z, Z, Z, O, O, Z, Z, Z, Z, O, O, Z, Z, Z, Z, O, O, Z, Z, Z, Z, O,
        ],
        [
            Z, O, Z, Z, Z, Z, O, Z, O, O, O, O, O, Z, Z, Z, Z, Z, O, Z, Z, Z, Z, Z,
        ],
        [
            Z, O, Z, Z, Z, Z, O, Z, Z, Z, Z, Z, O, Z, O, O, O, O, O, Z, Z, Z, Z, Z,
        ],
        [
            Z, O, Z, Z, Z, Z, O, Z, Z, Z, Z, Z, O, Z, Z, Z, Z, Z, O, Z, O, O, O, O,
        ],
        [
            Z, Z, O, Z, Z, Z, O, O, Z, O, Z, Z, O, Z, Z, Z, O, Z, O, Z, Z, Z, Z, O,
        ],
        [
            Z, Z, O, Z, Z, Z, O, Z, Z, Z, Z, O, O, O, Z, O, Z, Z, O, Z, Z, Z, O, Z,
        ],
        [
            Z, Z, O, Z, Z, Z, O, O, Z, O, Z, Z, Z, O, Z, Z, Z, O, Z, O, Z, Z, O, Z,
        ],
        [
            Z, Z, O, Z, Z, Z, Z, O, Z, Z, O, Z, O, O, Z, O, Z, Z, Z, O, Z, Z, Z, O,
        ],
    ]
    .into_iter()
    .map(|b| Vector::from(Function::from(b)))
    .collect::<Vec<_>>()
    .try_into()
    .unwrap()
}

fn cache() -> &'static ExtendedBinaryGolayCodeCache {
    EBGC_CACHE.get_or_init(|| {
        let points = Point::structure();
        let space = F2::structure().free_module(&points);
        let subspace =
            space.generated_submodule(ebgc_basis().into_iter().map(Function::from).collect());
        let subspace_structure = FinitelyFreeSubmoduleStructure::new(space, subspace);
        debug_assert_eq!(subspace_structure.rank(), 12);

        let mut octad_from_5_points = BTreeMap::new();
        for vector in subspace_structure
            .clone()
            .generate_all_elements()
            .map(Vector::from)
        {
            if vector.weight() == 8 {
                let points = vector.points().collect::<Vec<_>>();
                debug_assert_eq!(points.len(), 8);
                for a in 0..8 {
                    for b in 0..a {
                        for c in 0..b {
                            for d in 0..c {
                                for e in 0..d {
                                    octad_from_5_points.insert(
                                        Vector::from_points(&vec![
                                            &points[a], &points[b], &points[c], &points[d],
                                            &points[e],
                                        ]),
                                        vector.clone(),
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }

        ExtendedBinaryGolayCodeCache {
            subspace: subspace_structure,
            octad_from_5_points,
        }
    })
}

/// The 24 dimensional vector space structure over F2 with basis given by the points of the MOG
pub fn ebgc_ambient_space_structure() -> &'static AmbientSpace {
    cache().subspace.module()
}

/// The 12 dimensional vector subspace given by the extended binary Golay code
pub fn ebgc() -> &'static FinitelyFreeSubmodule<F2> {
    cache().subspace.submodule()
}

/// The 12 dimensional vector subspace structure given by the extended binary Golay code
pub fn ebgc_structure()
-> Arc<FinitelyFreeSubmoduleStructure<PointCanonicalStructure, F2Structure, AmbientSpace>> {
    cache().subspace.clone()
}

pub fn all_ebgc_codewords() -> Vec<Vector> {
    ebgc_structure()
        .list_all_elements()
        .into_iter()
        .map(Vector::from)
        .collect()
}

pub fn complete_octad(five_pts: &Vector) -> Vector {
    debug_assert_eq!(five_pts.weight(), 5);
    cache().octad_from_5_points.get(five_pts).unwrap().clone()
}

pub fn complete_sextet(four_pts: Vector) -> Sextet {
    let mut other_pts = (0usize..24).collect::<HashSet<_>>();
    for pt in four_pts.points() {
        other_pts.remove(&pt.element_to_enumeration().try_into().unwrap());
    }
    let mut foursomes = vec![four_pts.clone()];
    for i in 0..5 {
        debug_assert_eq!(other_pts.len(), 20 - 4 * i);
        let five_pts = Vector::from_points(
            &other_pts
                .iter()
                .take(1)
                .map(|i| Point::enumeration_to_element(&Natural::from(*i)).unwrap())
                .chain(four_pts.points())
                .collect(),
        );
        let octad = complete_octad(&five_pts);
        foursomes.push(Vector::from_points(
            &octad
                .points()
                .filter(|pt| !four_pts.contains_point(pt))
                .inspect(|pt| {
                    let i = pt.element_to_enumeration().try_into().unwrap();
                    debug_assert!(other_pts.contains(&i));
                    other_pts.remove(&i);
                })
                .collect(),
        ));
    }
    debug_assert!(other_pts.is_empty());
    debug_assert_eq!(foursomes.len(), 6);
    Sextet::from_foursomes(foursomes.try_into().unwrap())
}

/// Complete a labelling of an ordered sextet
/// T1: [x, ?, ?, ?]
/// T2: [y, z, ?, ?]
/// T3: [w, ?, ?, ?]
/// T4: [?, ?, ?, ?]
/// T5: [?, ?, ?, ?]
/// T6: [?, ?, ?, ?]
/// where
///  - x is labelled 0
///  - y is labelled 0
///  - z is labelled 1
///  - w is labelled alpha
pub fn complete_sextet_labelling(
    sextet: &OrderedSextet,
    x: &Point,
    y: &Point,
    z: &Point,
    w: &Point,
    alpha: F4,
) -> OrderedSextetLabelling {
    let foursomes = sextet.foursomes();
    let osp = |i: usize| OrderedSynthemePoint::enumeration_to_element(&i.into()).unwrap();
    assert!(foursomes.image(&osp(0)).contains_point(x));
    assert!(foursomes.image(&osp(1)).contains_point(y));
    assert!(foursomes.image(&osp(1)).contains_point(z));
    assert!(foursomes.image(&osp(2)).contains_point(w));
    assert_ne!(y, z);
    let mut labels = LabelledPoints::new(|_| F4::Zero);
    debug_assert_eq!(labels.image(x), &F4::Zero);
    debug_assert_eq!(labels.image(y), &F4::Zero);
    *labels.image_mut(z) = F4::One;
    *labels.image_mut(w) = alpha;

    let t0 = foursomes.image(&osp(0));
    let t1 = foursomes.image(&osp(1));
    let t2 = foursomes.image(&osp(2));
    let t3 = foursomes.image(&osp(3));
    let t4 = foursomes.image(&osp(4));
    let t5 = foursomes.image(&osp(5));

    let _ = t1; //It's not used

    #[allow(clippy::items_after_statements)]
    fn take_unique_pt(v: Vector) -> Point {
        let mut pts = v.points();
        let pt = pts.next().unwrap();
        debug_assert_eq!(pts.next(), None);
        pt
    }

    // Complete the hexacodeword (0, 1, alpha, ?, ?, ?)
    let (beta, gamma, delta) = match alpha {
        F4::Zero => (F4::One, F4::Alpha, F4::Beta),
        F4::One => (F4::Zero, F4::Beta, F4::Alpha),
        F4::Alpha => (F4::Beta, F4::Zero, F4::One),
        F4::Beta => (F4::Alpha, F4::One, F4::Zero),
    };
    // Use the octad containing (T1 \ {x}) U {z, w} to label 1 point in each of T3, T4, T5, T6
    let octad = complete_octad(&Vector::from_points(
        &t0.points()
            .filter(|p| p != x)
            .chain(vec![z.clone(), w.clone()])
            .collect(),
    ));
    let (w2, w3, w4, w5) = (
        take_unique_pt(&octad & t2),
        take_unique_pt(&octad & t3),
        take_unique_pt(&octad & t4),
        take_unique_pt(&octad & t5),
    );
    debug_assert_eq!(&w2, w);
    *labels.image_mut(&w3) = beta;
    *labels.image_mut(&w4) = gamma;
    *labels.image_mut(&w5) = delta;

    // Use the sextet formed by completing (T1 \ {x}) U {y} to label the rest of T3 U T4 U T5 U T6
    let t2345 = t2 | t3 | t4 | t5;
    debug_assert_eq!(t2345.weight(), 16);
    for (_i, li, wi) in [
        (2, alpha, &w2),
        (3, beta, &w3),
        (4, gamma, &w4),
        (5, delta, &w5),
    ] {
        let octad = complete_octad(&Vector::from_points(
            &t0.points()
                .filter(|p| p != x)
                .chain(vec![y.clone(), wi.clone()])
                .collect(),
        ));
        for p in (octad & &t2345).points() {
            *labels.image_mut(&p) = li;
        }
    }
    debug_assert_eq!(labels.image(&w2), &alpha);
    debug_assert_eq!(labels.image(&w3), &beta);
    debug_assert_eq!(labels.image(&w4), &gamma);
    debug_assert_eq!(labels.image(&w5), &delta);

    //Find the point labelled 0 in T4 and the three points not labelled 0 in T5
    let mut final_four = vec![];
    for p in t5.points() {
        if *labels.image(&p) != F4::Zero {
            final_four.push(p);
        }
    }
    for p in t4.points() {
        if *labels.image(&p) == F4::Zero {
            final_four.push(p);
        }
    }
    //Complete these 4 to an octad using each point in T3. The hexacodewords have labels llll00 so we can complete the labelling in T1 U T2
    debug_assert_eq!(final_four.len(), 4);
    for q in t3.points() {
        let l = *labels.image(&q);
        let octad = complete_octad(&Vector::from_points(
            &final_four.iter().cloned().chain(vec![q]).collect(),
        ));
        for p in octad.points() {
            if !t4.contains_point(&p) && !t5.contains_point(&p) {
                *labels.image_mut(&p) = l;
            }
        }
    }
    debug_assert_eq!(labels.image(x), &F4::Zero);
    debug_assert_eq!(labels.image(y), &F4::Zero);
    debug_assert_eq!(labels.image(z), &F4::One);
    debug_assert_eq!(labels.image(&w2), &alpha);
    debug_assert_eq!(labels.image(&w3), &beta);
    debug_assert_eq!(labels.image(&w4), &gamma);
    debug_assert_eq!(labels.image(&w5), &delta);

    #[cfg(debug_assertions)]
    {
        for i in 0..6 {
            let t = foursomes.image(&osp(i));
            let mut zero_count: usize = 0;
            let mut one_count: usize = 0;
            let mut alpha_count: usize = 0;
            let mut beta_count: usize = 0;
            for p in t.points() {
                match labels.image(&p) {
                    F4::Zero => zero_count += 1,
                    F4::One => one_count += 1,
                    F4::Alpha => alpha_count += 1,
                    F4::Beta => beta_count += 1,
                }
            }
            assert_eq!(zero_count, 1);
            assert_eq!(one_count, 1);
            assert_eq!(alpha_count, 1);
            assert_eq!(beta_count, 1);
        }
    }
    OrderedSextetLabelling {
        sextet: sextet.clone(),
        f4_labels: labels,
    }
}

#[allow(clippy::large_enum_variant)]
#[derive(Debug, Clone)]
pub enum NearestCodewordsResult {
    Unique { codeword: Vector, distance: usize },
    Six { codewords: [Vector; 6] },
}

impl NearestCodewordsResult {
    pub fn distance(&self) -> usize {
        match self {
            NearestCodewordsResult::Unique { distance, .. } => *distance,
            NearestCodewordsResult::Six { .. } => 4,
        }
    }
}

pub fn nearest_ebgc_codeword(vector: &Vector) -> NearestCodewordsResult {
    let mut dist_4_codewords = vec![];
    for codeword in all_ebgc_codewords() {
        let diff = vector + &codeword;
        let distance = diff.weight();
        if distance <= 3 {
            debug_assert!(dist_4_codewords.is_empty());
            return NearestCodewordsResult::Unique {
                codeword: codeword.clone(),
                distance,
            };
        } else if distance == 4 {
            dist_4_codewords.push(codeword);
        }
    }
    debug_assert_eq!(dist_4_codewords.len(), 6);
    NearestCodewordsResult::Six {
        codewords: std::array::from_fn(|i| dist_4_codewords[i].clone()),
    }
}

pub trait EbgcPointPermutation {
    fn is_ebgc_automorphism(&self) -> bool;
}

impl EbgcPointPermutation for ConstSizePermutation<24, Point> {
    fn is_ebgc_automorphism(&self) -> bool {
        for b in ebgc_basis() {
            if !Vector::from(
                ebgc_ambient_space_structure()
                    .functions_restructure()
                    .domain_precomposition_const_size_permutation_action()
                    .apply(self, &b.into()),
            )
            .is_codeword()
            {
                return false;
            }
        }
        true
    }
}

/// Return an element of M24 sending (a, b, c, d, e) to the following points of the MOG
///
/// a e  - -  - -
/// b -  - -  - -
/// c -  - -  - -
/// d -  - -  - -
///
/// There may be more than one such permutation, there is no guarantee which one this returns.
///
/// # Errors
/// If a, b, c, d, e are not all mutually distinct
#[allow(clippy::result_unit_err)]
pub fn find_5_transitive_permutation_to_standard_points(
    a: &Point,
    b: &Point,
    c: &Point,
    d: &Point,
    e: &Point,
) -> Result<ConstSizePermutation<24, Point>, ()> {
    /*
    Strategy:
    1. Take an ordered sextet starting with the tetrad {a, b, c, d}
    2. Using a sextet labelling, move {a, b, c, d} to the first standard tetrad
    3. Use hexacode automorphisms to put a, b, c, d into their places:
      - Add hexacodewords to put a in its place
      - Multiply by a scalar from F4 to put b in its place
      - Optionally conjugate to put c and d in their place
    4. An alternating permutation of the last 5 standard tetrads puts e in the second standard tetrad
    5. Adding a hexacodeword of the form 0x---- puts e in its place
     */

    let given_points = vec![a, b, c, d, e];
    let given_sorted_points = Point::structure().sort(given_points);
    if !Point::structure().is_sorted_and_unique(&given_sorted_points) {
        return Err(());
    }

    let first_tetrad = OrderedSynthemePoint::enumeration_to_element(&Natural::ZERO).unwrap();
    let second_tetrad = OrderedSynthemePoint::enumeration_to_element(&Natural::ONE).unwrap();
    let third_tetrad = OrderedSynthemePoint::enumeration_to_element(&Natural::TWO).unwrap();
    let fourth_tetrad = OrderedSynthemePoint::enumeration_to_element(&Natural::from(3u8)).unwrap();

    // make a labelling with a, b, c, d in the first sextet and e in the second sextet
    let mut ordered_sextet_labelling = {
        let mut ordered_sextet = complete_sextet(Vector::from_points(&vec![a, b, c, d]))
            .orderings()
            .next()
            .unwrap();

        // make {a, b, c, d} the first tetrad
        let abcd_tetrad = ordered_sextet.inner.image(a);
        if *abcd_tetrad != first_tetrad {
            ordered_sextet = ordered_sextet
                .permute(&ConstSizePermutation::new_swap(first_tetrad, *abcd_tetrad).unwrap());
        }

        // put e in the second tetrad
        let e_tetrad = ordered_sextet.inner.image(e);
        if *e_tetrad != second_tetrad {
            ordered_sextet = ordered_sextet
                .permute(&ConstSizePermutation::new_swap(second_tetrad, *e_tetrad).unwrap());
        }

        // any point in the second tetrad different from e
        let f = {
            let ordered_sextet_foursomes = ordered_sextet.foursomes();
            let mut second_tetrad_points = ordered_sextet_foursomes.image(&second_tetrad).points();
            let f = second_tetrad_points.next().unwrap();
            if *e != f {
                f
            } else {
                second_tetrad_points.next().unwrap()
            }
        };

        // any point in the third tetrad
        let g = ordered_sextet
            .foursomes()
            .image(&third_tetrad)
            .points()
            .next()
            .unwrap();

        complete_sextet_labelling(&ordered_sextet, a, e, &f, &g, F4::Zero)
    };

    debug_assert!(
        ordered_sextet_labelling
            .permutation_to_standard_labelling()
            .is_ebgc_automorphism()
    );

    // add a vector to give `a` the label `F4::Zero`
    let hexacode_vector = HexacodeVector::new(|q| match q.pair {
        OrderedSynthemePair::Left => F4::One,
        OrderedSynthemePair::Middle => F4::One,
        OrderedSynthemePair::Right => F4::Zero,
    });
    let a_label = *ordered_sextet_labelling.f4_labels.image(a);
    ordered_sextet_labelling = ordered_sextet_labelling.add_vector(&(a_label * &hexacode_vector));
    debug_assert_eq!(ordered_sextet_labelling.f4_labels.image(a), &F4::Zero);
    debug_assert!(
        ordered_sextet_labelling
            .permutation_to_standard_labelling()
            .is_ebgc_automorphism()
    );

    // multiply by a scalar to give `b` the label `F4::ONE`
    let b_label = *ordered_sextet_labelling.f4_labels.image(b);
    ordered_sextet_labelling = ordered_sextet_labelling.scalar_mul(b_label);
    debug_assert_eq!(ordered_sextet_labelling.f4_labels.image(a), &F4::Zero);
    debug_assert_eq!(ordered_sextet_labelling.f4_labels.image(b), &F4::One);
    debug_assert!(
        ordered_sextet_labelling
            .permutation_to_standard_labelling()
            .is_ebgc_automorphism()
    );

    // conjugate if necessary to give `c` and `d` the labels `F4::Alpha` and `F4::Beta` respectively
    if *ordered_sextet_labelling.f4_labels.image(c) != F4::Alpha {
        ordered_sextet_labelling = ordered_sextet_labelling.conjugate();
        ordered_sextet_labelling = ordered_sextet_labelling.permute_foursomes(
            &ConstSizePermutation::new_swap(third_tetrad, fourth_tetrad).unwrap(),
        );
    }
    debug_assert_eq!(ordered_sextet_labelling.f4_labels.image(a), &F4::Zero);
    debug_assert_eq!(ordered_sextet_labelling.f4_labels.image(b), &F4::One);
    debug_assert_eq!(ordered_sextet_labelling.f4_labels.image(c), &F4::Alpha);
    debug_assert_eq!(ordered_sextet_labelling.f4_labels.image(d), &F4::Beta);
    debug_assert!(
        ordered_sextet_labelling
            .permutation_to_standard_labelling()
            .is_ebgc_automorphism()
    );

    // add a vector to give `e` the label `F4::Zero`
    let hexacode_vector = HexacodeVector::new(|q| match q.pair {
        OrderedSynthemePair::Left | OrderedSynthemePair::Middle => match q.side {
            OrderedSynthemeSide::Left => F4::Zero,
            OrderedSynthemeSide::Right => F4::One,
        },
        OrderedSynthemePair::Right => match q.side {
            OrderedSynthemeSide::Left => F4::Alpha,
            OrderedSynthemeSide::Right => F4::Beta,
        },
    });
    let e_label = *ordered_sextet_labelling.f4_labels.image(e);
    ordered_sextet_labelling = ordered_sextet_labelling.add_vector(&(e_label * &hexacode_vector));
    debug_assert_eq!(ordered_sextet_labelling.f4_labels.image(a), &F4::Zero);
    debug_assert_eq!(ordered_sextet_labelling.f4_labels.image(b), &F4::One);
    debug_assert_eq!(ordered_sextet_labelling.f4_labels.image(c), &F4::Alpha);
    debug_assert_eq!(ordered_sextet_labelling.f4_labels.image(d), &F4::Beta);
    debug_assert_eq!(ordered_sextet_labelling.f4_labels.image(e), &F4::Zero);

    let aut = ordered_sextet_labelling.permutation_to_standard_labelling();
    debug_assert!(aut.is_ebgc_automorphism());
    Ok(aut)
}

/// Return an element of M24 sending 5 distinct points to 5 other distinct points
///
/// # Errors
/// If any of the 5 inputs are not mutually distinct or if any of the 5 outputs are not mutually distinct
#[allow(clippy::result_unit_err)]
pub fn find_5_transitive_permutation(
    a: (&Point, &Point),
    b: (&Point, &Point),
    c: (&Point, &Point),
    d: (&Point, &Point),
    e: (&Point, &Point),
) -> Result<ConstSizePermutation<24, Point>, ()> {
    Ok(
        find_5_transitive_permutation_to_standard_points(a.1, b.1, c.1, d.1, e.1)?
            .inverse()
            .compose(&find_5_transitive_permutation_to_standard_points(
                a.0, b.0, c.0, d.0, e.0,
            )?),
    )
}

#[cfg(test)]
mod tests {
    use crate::golay_codes::hexacode::hexacode_ambient_space_structure;

    use super::*;
    use algebraeon_groups::examples::c2::C2;
    use algebraeon_rings::{
        linear::monomial_transformations::{
            GaloisMonomialTransformationsSignature, MonomialTransformationsSupersetSignature,
        },
        structure::FreeModuleSignature,
    };
    use algebraeon_sets::sets::SetToConstSizePermutationsStructure;

    #[test]
    fn weight_distribution() {
        let mut wt0 = 0;
        let mut wt8 = 0;
        let mut wt12 = 0;
        let mut wt16 = 0;
        let mut wt24 = 0;
        for x in all_ebgc_codewords() {
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
    fn test_complete_octad() {
        let p = |i: usize| Point::enumeration_to_element(&i.into()).unwrap();
        let v = Vector::from_points(&vec![p(0), p(1), p(2), p(3), p(4)]);
        let w = complete_octad(&v);
        assert!(w.is_octad());
    }

    #[test]
    fn test_complete_sextet() {
        const Z: F2 = ZERO;
        const O: F2 = ONE;

        let a: Vector = [
            O, O, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, O, O,
        ]
        .into();

        let s = complete_sextet(a);

        println!("{:#?}", s);
    }

    #[test]
    fn test_complete_sextet_labelling() {
        let sextet = OrderedSextet {
            inner: LabelledPoints::new(|p| p.col),
        };

        let p = |i: usize| Point::enumeration_to_element(&i.into()).unwrap();

        println!("{:#?}", sextet);

        let sextet_labelling =
            complete_sextet_labelling(&sextet, &p(0), &p(1), &p(7), &p(2), F4::Zero);

        println!("{:#?}", sextet_labelling);
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

    #[test]
    fn test_is_ebgc_automorphism() {
        let p = |i: usize| -> Point {
            Point::structure()
                .enumeration_to_element(&Natural::from(i))
                .unwrap()
        };

        let perm1 = Point::structure()
            .const_size_permutations()
            .new_cycles(vec![vec![p(0), p(1), p(2)]])
            .unwrap();
        assert!(!perm1.is_ebgc_automorphism());

        let perm2 = Point::structure()
            .const_size_permutations()
            .new_cycles(vec![
                vec![p(0), p(1)],
                vec![p(2), p(3)],
                vec![p(6), p(7)],
                vec![p(8), p(9)],
                vec![p(12), p(13)],
                vec![p(14), p(15)],
                vec![p(18), p(19)],
                vec![p(20), p(21)],
            ])
            .unwrap();
        assert!(perm2.is_ebgc_automorphism());
    }

    #[test]
    fn automorphism_from_hexacode_automorphism() {
        let hexacode_perms = hexacode_ambient_space_structure()
            .basis_set()
            .const_size_permutations();

        let hexacode_aut = hexacode_ambient_space_structure().galois_monomial_transformations();

        let p = |n: u8| -> OrderedSynthemePoint {
            OrderedSynthemePoint::enumeration_to_element(&n.into()).unwrap()
        };

        let aut1 = hexacode_aut.new_galois_automorphism(&C2::Flip);
        let aut2 = hexacode_aut.new_permutation(
            &hexacode_perms
                .new_cycles(vec![vec![p(0), p(1)], vec![p(2), p(3)], vec![p(4), p(5)]])
                .unwrap(),
        );
        let aut3 = hexacode_aut
            .new_scalars(&[F4::One, F4::One, F4::One, F4::One, F4::One, F4::One].into());

        println!("{:?} {:?} {:?}", aut1, aut2, aut3);

        todo!();
    }

    #[test]
    fn validate_ordered_sextet() {
        assert!(
            OrderedSextet::try_from(LabelledPoints::new(|p| {
                let i: usize = p.element_to_enumeration().try_into().unwrap();
                let q = |i: usize| -> OrderedSynthemePoint {
                    OrderedSynthemePoint::enumeration_to_element(&i.into()).unwrap()
                };
                match i {
                    0 | 1 | 6 | 7 => q(0),
                    2 | 3 | 8 | 9 => q(1),
                    4 | 5 | 10 | 11 => q(2),
                    12 | 13 | 18 | 19 => q(3),
                    14 | 15 | 20 | 21 => q(4),
                    16 | 17 | 22 | 23 => q(5),
                    _ => unreachable!(),
                }
            }))
            .is_ok()
        );

        assert!(
            OrderedSextet::try_from(LabelledPoints::new(|p| {
                let i: usize = p.element_to_enumeration().try_into().unwrap();
                let q = |i: usize| -> OrderedSynthemePoint {
                    OrderedSynthemePoint::enumeration_to_element(&i.into()).unwrap()
                };
                match i {
                    2 | 1 | 6 | 7 => q(0),
                    0 | 3 | 8 | 9 => q(1),
                    4 | 5 | 10 | 11 => q(2),
                    12 | 13 | 18 | 19 => q(3),
                    14 | 15 | 20 | 21 => q(4),
                    16 | 17 | 22 | 23 => q(5),
                    _ => unreachable!(),
                }
            }))
            .is_err()
        );
    }

    #[test]
    fn test_to_mat_from_mat_inv() {
        let x = Vector::from_fn(|p| {
            let i: usize = p.element_to_enumeration().try_into().unwrap();
            if i.is_multiple_of(5) {
                F2::one()
            } else {
                F2::zero()
            }
        });
        assert_eq!(x.to_mat(), Vector::from_mat(&x.to_mat()).unwrap().to_mat());
    }

    #[test]
    fn test_5_transitive() {
        let p = |n: usize| -> Point { Point::enumeration_to_element(&n.into()).unwrap() };

        let perm =
            find_5_transitive_permutation_to_standard_points(&p(7), &p(9), &p(21), &p(6), &p(5))
                .unwrap();

        assert!(perm.is_ebgc_automorphism());
        assert_eq!(perm.image(&p(7)), p(0));
        assert_eq!(perm.image(&p(9)), p(6));
        assert_eq!(perm.image(&p(21)), p(12));
        assert_eq!(perm.image(&p(6)), p(18));
        assert_eq!(perm.image(&p(5)), p(1));
    }
}
