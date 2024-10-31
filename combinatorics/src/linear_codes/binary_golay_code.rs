use std::{
    borrow::Borrow,
    collections::{HashMap, HashSet},
};

use derivative::Derivative;
use itertools::Itertools;
use malachite_base::num::logic::traits::BitIterable;

use orthoclase_groups::examples::symmetric::*;
use orthoclase_rings::{
    linear::matrix::Matrix,
    number::{modulo::Modulo, quaternary_field::QuaternaryField},
    ring_structure::cannonical::*,
};

pub type BinaryField = Modulo<2>;
pub const ZERO: BinaryField = BinaryField::new(0);
pub const ONE: BinaryField = BinaryField::new(1);

#[derive(Copy, Clone, PartialEq, Eq, Hash)]
pub struct Vector24 {
    /*
    Use the first 24 bits to represent an vector in the MOG as follows.
    The last 8 bits must be zero
     0  1  2  3  4  5
     6  7  8  8  10 11
     12 13 14 15 16 17
     18 19 20 21 22 23
    */
    points: u32,
}
impl core::fmt::Debug for Vector24 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut bits = self.points.bits();
        let mut pts = vec![];
        for i in 0..24usize {
            if match bits.next() {
                Some(b) => b,
                None => false,
            } {
                pts.push(i);
            }
        }
        f.debug_struct("Vector").field("points", &pts).finish()
    }
}
impl Vector24 {
    pub fn new(contains: impl Fn(usize) -> bool) -> Self {
        let mut points = 0;
        for i in 0..24usize {
            match contains(i) {
                true => points |= 1 << i,
                false => {}
            }
        }
        Self { points }
    }
    pub fn zero() -> Self {
        Self { points: 0 }
    }
    pub fn from_points(pts: impl IntoIterator<Item = usize>) -> Self {
        let mut data = 0;
        for i in pts {
            assert!(i < 24);
            data |= 1 << i
        }
        Self { points: data }
    }
    pub fn points<'a>(&'a self) -> impl Iterator<Item = usize> + 'a {
        (0..24).filter(|i| self.contains(*i))
    }
    pub fn contains(&self, pt: usize) -> bool {
        assert!(pt < 24);
        self.points & (1 << pt) != 0
    }
    pub fn into_row(&self) -> Matrix<BinaryField> {
        Matrix::construct(1, 24, |r, c| {
            debug_assert_eq!(r, 0);
            match { self.contains(c) } {
                true => ONE,
                false => ZERO,
            }
        })
    }
    pub fn from_row(m: &Matrix<BinaryField>) -> Option<Self> {
        if m.rows() == 1 && m.cols() == 24 {
            Some(Self::new(|i| m.at(0, i).unwrap() == &ONE))
        } else {
            None
        }
    }
    pub fn into_col(&self) -> Matrix<BinaryField> {
        self.into_row().transpose()
    }
    pub fn from_col(m: &Matrix<BinaryField>) -> Option<Self> {
        Self::from_row(&m.transpose_ref())
    }
    /// As a 4x6 matrix representing an element of the MOG
    pub fn into_mat(&self) -> Matrix<BinaryField> {
        Matrix::construct(4, 6, |r, c| match self.contains(c + 6 * r) {
            true => ONE,
            false => ZERO,
        })
    }
    /// From a 4x6 matrix representing an element of the MOG
    pub fn from_mat(m: &Matrix<BinaryField>) -> Option<Self> {
        if m.rows() == 4 && m.cols() == 6 {
            Some(Self::new(|i| {
                let (r, c) = (i / 6, i % 6);
                m.at(r, c).unwrap() == &ONE
            }))
        } else {
            None
        }
    }
}

impl std::ops::Add<Vector24> for Vector24 {
    type Output = Self;
    fn add(self, other: Vector24) -> Self::Output {
        Self {
            points: self.points ^ other.points,
        }
    }
}

impl std::ops::BitAnd<Vector24> for Vector24 {
    type Output = Self;
    fn bitand(self, other: Vector24) -> Self::Output {
        Self {
            points: self.points & other.points,
        }
    }
}

impl std::ops::BitOr<Vector24> for Vector24 {
    type Output = Self;
    fn bitor(self, other: Vector24) -> Self::Output {
        Self {
            points: self.points | other.points,
        }
    }
}

impl Vector24 {
    pub fn weight(&self) -> usize {
        let mut t = 0;
        for i in 0..24 {
            if self.contains(i) {
                t += 1;
            }
        }
        t
    }
}

macro_rules! define_special_vector_struct {
    ($name:ident) => {
        #[derive(Derivative, Clone, Copy)]
        #[derivative(Debug)]
        pub struct $name<'g, V: Borrow<Vector24>> {
            #[derivative(Debug = "ignore")]
            ebgc: &'g ExtendedBinaryGolayCode,
            vec: V,
        }
        impl<'g, V: Borrow<Vector24>> From<$name<'g, V>> for Vector24 {
            fn from(octad: $name<V>) -> Self {
                *octad.vec.borrow()
            }
        }
        impl<'g, V: Borrow<Vector24>> $name<'g, V> {
            pub fn points<'a>(&'a self) -> impl 'a + Iterator<Item = usize> {
                self.vec.borrow().points()
            }
            pub fn contains(&self, i: usize) -> bool {
                self.vec.borrow().contains(i)
            }
        }
    };
}
define_special_vector_struct!(Codeword);
define_special_vector_struct!(Octad);
define_special_vector_struct!(Foursome);

#[derive(Derivative, Clone)]
#[derivative(Debug)]
pub struct OrderedSextet<'g> {
    #[derivative(Debug = "ignore")]
    ebgc: &'g ExtendedBinaryGolayCode,
    parts: [Vector24; 6],
}
impl<'g> OrderedSextet<'g> {
    pub fn foursome(&self, i: usize) -> Foursome<Vector24> {
        assert!(i < 6);
        Foursome {
            ebgc: self.ebgc,
            vec: self.parts[i],
        }
    }
}

#[derive(Derivative, Clone)]
#[derivative(Debug)]
pub struct UnorderedSextet<'g> {
    #[derivative(Debug = "ignore")]
    ebgc: &'g ExtendedBinaryGolayCode,
    parts: [Vector24; 6],
}
impl<'g> UnorderedSextet<'g> {
    pub fn orderings<'a>(&'a self) -> impl Iterator<Item = OrderedSextet> + 'a {
        (0..6).permutations(6).map(|perm| {
            let mut parts = [Vector24::zero(); 6];
            for (i, j) in perm.into_iter().enumerate() {
                parts[i] = self.parts[j];
            }
            OrderedSextet {
                ebgc: self.ebgc,
                parts,
            }
        })
    }
}

#[derive(Derivative, Clone)]
#[derivative(Debug)]
pub struct SextetLabelling<'g> {
    #[derivative(Debug = "ignore")]
    ebgc: &'g ExtendedBinaryGolayCode,
    sextet: &'g OrderedSextet<'g>,
    labels: [QuaternaryField; 24],
}

impl<'g> SextetLabelling<'g> {
    /// Return a pair of permutations (sigma, tau)
    /// sigma permutes the coordinates mapping this binary golay code to the miracle octad generator
    /// tau permutes the coordinates mapping the miracle octad generator to this binary golay code
    pub fn mog_isomorphism(&self) -> (Permutation<24>, Permutation<24>) {
        let mut to_mog = [0; 24];
        let mut from_mog = [0; 24];
        for c in 0..6 {
            for p in self.sextet.foursome(c).points() {
                let l = self.labels[p];
                let r = match l {
                    QuaternaryField::Zero => 0,
                    QuaternaryField::One => 1,
                    QuaternaryField::Alpha => 2,
                    QuaternaryField::Beta => 3,
                };
                let q = c + 6 * r;
                to_mog[p] = q;
                from_mog[q] = p;
            }
        }
        (
            Permutation::new(to_mog).unwrap(),
            Permutation::new(from_mog).unwrap(),
        )
    }
}

/// Represents a length=24 dimension=12 linear code over F2 with minimum weight (at least) 8
#[derive(Debug)]
pub struct ExtendedBinaryGolayCode {
    blocks: HashMap<Vector24, Vector24>,
}

impl ExtendedBinaryGolayCode {
    pub fn from_row_basis_matrix(m: Matrix<BinaryField>) -> Result<Self, &'static str> {
        if m.rows() != 12 || m.cols() != 24 {
            return Err("Matrix has the wrong dimensions");
        }
        if m.rank() != 12 {
            return Err("Matrix does not have full rank");
        }
        let basis = (0..12)
            .map(|i| Vector24::from_row(&m.get_row(i)).unwrap())
            .collect::<Vec<_>>();

        let mut blocks = HashMap::new();

        // Iterate over the span of the basis
        for coeffs in (0..12).map(|_| vec![false, true]).multi_cartesian_product() {
            let mut v = Vector24::zero();
            for (i, c) in coeffs.into_iter().enumerate() {
                if c {
                    v = v + basis[i];
                }
            }
            let wt = v.weight();
            if 0 < wt && wt < 8 {
                return Err("Matrix span contains vector(s) of non-zero weight <8");
            } else if wt == 8 {
                for pts in v.points().combinations(5) {
                    blocks.insert(Vector24::from_points(pts), v);
                }
            }
        }
        debug_assert_eq!(blocks.len(), 42504); // 24 choose 5
        Ok(Self { blocks })
    }

    pub fn from_col_basis_matrix(m: Matrix<BinaryField>) -> Result<Self, &'static str> {
        Self::from_row_basis_matrix(m.transpose())
    }
}

impl ExtendedBinaryGolayCode {
    pub fn complete_octad(&self, five_pts: Vector24) -> Octad<&Vector24> {
        assert_eq!(five_pts.weight(), 5);
        let oct = self.blocks.get(&five_pts).unwrap();
        debug_assert_eq!(oct.weight(), 8);
        Octad {
            ebgc: self,
            vec: oct,
        }
    }

    pub fn complete_sextet(&self, four_pts: Vector24) -> UnorderedSextet {
        let mut other_pts = (0..24).collect::<HashSet<_>>();
        for pt in four_pts.points() {
            other_pts.remove(&pt);
        }
        let mut parts = [Vector24::zero(); 6];
        parts[0] = four_pts;
        for i in 0..5 {
            debug_assert_eq!(other_pts.len(), 20 - 4 * i);
            let five_pts =
                Vector24::from_points(other_pts.iter().take(1).cloned().chain(four_pts.points()));
            let octad = self.complete_octad(five_pts);
            parts[i + 1] = Vector24::from_points(
                Vector24::from(octad)
                    .points()
                    .filter(|pt| !four_pts.contains(*pt))
                    .map(|pt| {
                        debug_assert!(other_pts.contains(&pt));
                        other_pts.remove(&pt);
                        pt
                    }),
            );
        }
        debug_assert!(other_pts.is_empty());
        UnorderedSextet { ebgc: self, parts }
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
    pub fn complete_labelling<'g>(
        &'g self,
        sextet: &'g OrderedSextet<'g>,
        x: usize,
        y: usize,
        z: usize,
        w: usize,
        alpha: QuaternaryField,
    ) -> SextetLabelling<'g> {
        assert!(std::ptr::eq(sextet.ebgc, self));
        assert!(sextet.foursome(0).contains(x));
        assert!(sextet.foursome(1).contains(y));
        assert!(sextet.foursome(1).contains(z));
        assert!(sextet.foursome(2).contains(w));
        assert_ne!(y, z);
        let mut labels = [QuaternaryField::Zero; 24];
        debug_assert_eq!(labels[x], QuaternaryField::Zero);
        debug_assert_eq!(labels[y], QuaternaryField::Zero);
        labels[z] = QuaternaryField::One;
        labels[w] = alpha;

        let t0 = Vector24::from(sextet.foursome(0));
        let t1 = Vector24::from(sextet.foursome(1));
        let t2 = Vector24::from(sextet.foursome(2));
        let t3 = Vector24::from(sextet.foursome(3));
        let t4 = Vector24::from(sextet.foursome(4));
        let t5 = Vector24::from(sextet.foursome(5));

        let _ = t1; //It's not used

        fn take_unique_pt(v: Vector24) -> usize {
            let mut pts = v.points();
            let pt = pts.next().unwrap();
            debug_assert_eq!(pts.next(), None);
            pt
        }

        // Complete the hexacodeword (0, 1, alpha, ?, ?, ?)
        let (beta, gamma, delta) = match alpha {
            QuaternaryField::Zero => (
                QuaternaryField::One,
                QuaternaryField::Alpha,
                QuaternaryField::Beta,
            ),
            QuaternaryField::One => (
                QuaternaryField::Zero,
                QuaternaryField::Beta,
                QuaternaryField::Alpha,
            ),
            QuaternaryField::Alpha => (
                QuaternaryField::Beta,
                QuaternaryField::Zero,
                QuaternaryField::One,
            ),
            QuaternaryField::Beta => (
                QuaternaryField::Alpha,
                QuaternaryField::One,
                QuaternaryField::Zero,
            ),
        };
        // Use the octad containing (T1 \ {x}) U {z, w} to label 1 point in each of T3, T4, T5, T6
        let octad = self.complete_octad(Vector24::from_points(
            t0.points().filter(|p| *p != x).chain(vec![z, w]),
        ));
        let (w2, w3, w4, w5) = (
            take_unique_pt(Vector24::from(octad) & t2),
            take_unique_pt(Vector24::from(octad) & t3),
            take_unique_pt(Vector24::from(octad) & t4),
            take_unique_pt(Vector24::from(octad) & t5),
        );
        debug_assert_eq!(w2, w);
        labels[w3] = beta;
        labels[w4] = gamma;
        labels[w5] = delta;

        // Use the sextet formed by completing (T1 \ {x}) U {y} to label the rest of T3 U T4 U T5 U T6
        let t2345 = t2 | t3 | t4 | t5;
        debug_assert_eq!(t2345.weight(), 16);
        for (_i, li, wi) in vec![
            (2, alpha, w2),
            (3, beta, w3),
            (4, gamma, w4),
            (5, delta, w5),
        ] {
            let octad = Vector24::from(self.complete_octad(Vector24::from_points(
                t0.points().filter(|p| *p != x).chain(vec![y, wi]),
            )));
            for p in (octad & t2345).points() {
                labels[p] = li;
            }
        }
        debug_assert_eq!(labels[w2], alpha);
        debug_assert_eq!(labels[w3], beta);
        debug_assert_eq!(labels[w4], gamma);
        debug_assert_eq!(labels[w5], delta);

        //Find the point labelled 0 in T4 and the three points not labelled 0 in T5
        let mut final_four = vec![];
        for p in t5.points() {
            if labels[p] != QuaternaryField::Zero {
                final_four.push(p);
            }
        }
        for p in t4.points() {
            if labels[p] == QuaternaryField::Zero {
                final_four.push(p);
            }
        }
        //Complete these 4 to an octad using each point in T3. The hexacodewords have labels llll00 so we can complete the labelling in T1 U T2
        debug_assert_eq!(final_four.len(), 4);
        for q in t3.points() {
            let l = labels[q];
            let octad = self.complete_octad(Vector24::from_points(
                final_four.iter().cloned().chain(vec![q]),
            ));
            for p in octad.points() {
                if !t4.contains(p) && !t5.contains(p) {
                    labels[p] = l;
                }
            }
        }
        debug_assert_eq!(labels[x], QuaternaryField::Zero);
        debug_assert_eq!(labels[y], QuaternaryField::Zero);
        debug_assert_eq!(labels[z], QuaternaryField::One);
        debug_assert_eq!(labels[w2], alpha);
        debug_assert_eq!(labels[w3], beta);
        debug_assert_eq!(labels[w4], gamma);
        debug_assert_eq!(labels[w5], delta);

        #[cfg(debug_assertions)]
        {
            for i in 0..6 {
                let t = sextet.foursome(i);
                let mut zero_count: usize = 0;
                let mut one_count: usize = 0;
                let mut alpha_count: usize = 0;
                let mut beta_count: usize = 0;
                for p in t.points() {
                    match labels[p] {
                        QuaternaryField::Zero => zero_count += 1,
                        QuaternaryField::One => one_count += 1,
                        QuaternaryField::Alpha => alpha_count += 1,
                        QuaternaryField::Beta => beta_count += 1,
                    }
                }
                assert_eq!(zero_count, 1);
                assert_eq!(one_count, 1);
                assert_eq!(alpha_count, 1);
                assert_eq!(beta_count, 1);
            }
        }
        SextetLabelling {
            ebgc: self,
            sextet,
            labels,
        }
    }
}

fn ebgc_from_string(values: &'static str) -> ExtendedBinaryGolayCode {
    let rows = values.split("\n").collect::<Vec<_>>();
    for row in &rows {
        assert_eq!(row.len(), 24);
    }
    ExtendedBinaryGolayCode::from_row_basis_matrix(Matrix::construct(12, 24, |r, c| {
        match rows.iter().nth(r).unwrap().chars().nth(c).unwrap() {
            '0' => BinaryField::zero(),
            '1' => BinaryField::one(),
            _ => panic!(),
        }
    }))
    .unwrap()
}

pub fn ebgc_wiki() -> ExtendedBinaryGolayCode {
    ebgc_from_string(
        "\
100000000000100111110001
010000000000010011111010
001000000000001001111101
000100000000100100111110
000010000000110010011101
000001000000111001001110
000000100000111100100101
000000010000111110010010
000000001000011111001001
000000000100001111100110
000000000010010101010111
000000000001101010101011",
    )
}

pub fn ebgc_qr11() -> ExtendedBinaryGolayCode {
    ebgc_from_string(
        "\
100000000000101000111011
010000000000110100011101
001000000000011010001111
000100000000101101000111
000010000000110110100011
000001000000111011010001
000000100000011101101001
000000010000001110110101
000000001000000111011011
000000000100100011101101
000000000010010001110111
000000000001111111111110",
    )
}

pub fn ebgc_mog() -> ExtendedBinaryGolayCode {
    ebgc_from_string(
        "\
110000110000110000110000
101000101000101000101000
100100100100100100100100
100010100010100010100010
100001100001100001100001
010000101111100000100000
010000100000101111100000
010000100000100000101111
001000110100100010100001
001000100001110100100010
001000110100010001010010
001000010010110100010001",
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn making_ebgc() {
        ebgc_wiki();
        ebgc_qr11();
        ebgc_mog();
    }

    #[test]
    fn sextet_labelling() {
        let g24 = ebgc_mog();
        let unordered_sextet = g24.complete_sextet(Vector24::from_points(vec![3, 4, 7, 8]));
        let sextet = unordered_sextet.orderings().next().unwrap();
        let x = Vector24::from(sextet.foursome(0)).points().next().unwrap();
        let (y, z) = {
            let foursome = Vector24::from(sextet.foursome(1));
            let mut pts = foursome.points();
            let y = pts.next().unwrap();
            let z = pts.next().unwrap();
            (y, z)
        };
        let w = Vector24::from(sextet.foursome(2)).points().next().unwrap();

        g24.complete_labelling(&sextet, x, y, z, w, QuaternaryField::Zero);
        g24.complete_labelling(&sextet, x, y, z, w, QuaternaryField::One);
        g24.complete_labelling(&sextet, x, y, z, w, QuaternaryField::Alpha);
        g24.complete_labelling(&sextet, x, y, z, w, QuaternaryField::Beta);
    }

    #[test]
    fn test_into_mat_from_mat_inv() {
        let x = Vector24::new(|i| i % 5 == 0);
        assert_eq!(
            x.into_mat(),
            Vector24::from_mat(&x.into_mat()).unwrap().into_mat()
        );
    }
}
