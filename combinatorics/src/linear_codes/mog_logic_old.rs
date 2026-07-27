#![allow(dead_code)]

pub mod hexacode {
    use super::finite_field_4::Point as F4Point;
    use super::traits::{Enumerated, Labelled};

    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    pub enum Side {
        Left,
        Right,
    }

    impl Enumerated for Side {
        // 0 1

        const N: usize = 2;

        fn usize_to_point(i: usize) -> Result<Self, ()> {
            match i {
                0 => Ok(Self::Left),
                1 => Ok(Self::Right),
                _ => Err(()),
            }
        }

        fn point_to_usize(&self) -> usize {
            match self {
                Self::Left => 0,
                Self::Right => 1,
            }
        }
    }

    impl Side {
        pub fn flip(self) -> Self {
            match self {
                Self::Left => Self::Right,
                Self::Right => Self::Left,
            }
        }
    }

    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    pub enum Pair {
        Left,
        Middle,
        Right,
    }

    impl Enumerated for Pair {
        // 0 1 2

        const N: usize = 3;

        fn usize_to_point(i: usize) -> Result<Self, ()> {
            match i {
                0 => Ok(Self::Left),
                1 => Ok(Self::Middle),
                2 => Ok(Self::Right),
                _ => Err(()),
            }
        }

        fn point_to_usize(&self) -> usize {
            match self {
                Self::Left => 0,
                Self::Middle => 1,
                Self::Right => 2,
            }
        }
    }

    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    pub struct Point {
        pub side: Side,
        pub pair: Pair,
    }

    impl Enumerated for Point {
        // 0 1  2 3  4 5

        const N: usize = 6;

        fn usize_to_point(i: usize) -> Result<Self, ()> {
            if i < 6 {
                Ok(Self {
                    side: Side::usize_to_point(i % 2).unwrap(),
                    pair: Pair::usize_to_point(i / 2).unwrap(),
                })
            } else {
                Err(())
            }
        }

        fn point_to_usize(&self) -> usize {
            self.side.point_to_usize() + 2 * self.pair.point_to_usize()
        }
    }

    pub type Vector = Labelled<Point, F4Point>;

    impl Vector {
        fn component(&self, p: Point) -> F4Point {
            *self.get(p)
        }
    }
}

pub mod miracle_octad_generator {
    use super::finite_field_4::Point as F4Point;
    use super::{
        hexacode,
        permutation::Permutation,
        traits::{Enumerated, Labelled},
    };
    use std::{
        collections::HashSet,
        ops::{Add, BitAnd, BitOr},
        vec,
    };

    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
    pub struct Point {
        pub col: hexacode::Point,
        pub row: F4Point,
    }

    impl Enumerated for Point {
        /*
        MOG numbering:
        0  |  0  1   2  3   4  5
        1  |  6  7   8  9   10 11
        a  |  12 13  14 15  16 17
        b  |  18 19  20 21  22 23
        */

        const N: usize = 24;

        fn usize_to_point(i: usize) -> Result<Self, ()> {
            if i < 24 {
                Ok(Self {
                    col: hexacode::Point::usize_to_point(i % 6).unwrap(),
                    row: F4Point::usize_to_point(i / 6).unwrap(),
                })
            } else {
                Err(())
            }
        }

        fn point_to_usize(&self) -> usize {
            self.col.point_to_usize() + 6 * self.row.point_to_usize()
        }
    }

    pub type Vector = Labelled<Point, bool>;

    impl Add<&Vector> for &Vector {
        type Output = Vector;

        fn add(self, other: &Vector) -> Self::Output {
            Vector::from_fn(|p| self.contains_point(p) ^ other.contains_point(p))
        }
    }

    impl BitAnd<&Vector> for &Vector {
        type Output = Vector;

        fn bitand(self, other: &Vector) -> Self::Output {
            Vector::from_fn(|p| self.contains_point(p) && other.contains_point(p))
        }
    }

    impl BitOr<&Vector> for &Vector {
        type Output = Vector;

        fn bitor(self, other: &Vector) -> Self::Output {
            Vector::from_fn(|p| self.contains_point(p) || other.contains_point(p))
        }
    }

    impl PartialOrd for Vector {
        fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
            Some(self.cmp(other))
        }
    }

    impl Ord for Vector {
        fn cmp(&self, other: &Self) -> std::cmp::Ordering {
            let points = hexacode::Point::points()
                .flat_map(|col| F4Point::points().map(move |row| Point { col, row }))
                .collect::<Vec<_>>();
            points
                .iter()
                .map(|p| self.get(*p))
                .collect::<Vec<_>>()
                .cmp(&points.iter().map(|p| other.get(*p)).collect::<Vec<_>>())
        }
    }

    impl Vector {
        pub fn zero() -> Self {
            Self::new_constant(false)
        }

        pub fn contains(&self, other: &Self) -> bool {
            for p in Point::points() {
                if !self.contains_point(p) && other.contains_point(p) {
                    return false;
                }
            }
            true
        }

        pub fn contains_point(&self, p: Point) -> bool {
            *self.get(p)
        }

        pub fn weight(&self) -> usize {
            let mut w = 0;
            for p in Point::points() {
                if self.contains_point(p) {
                    w += 1;
                }
            }
            w
        }

        pub fn points(&self) -> impl Iterator<Item = Point> {
            Point::points().filter(|p| self.contains_point(*p))
        }

        pub fn from_points(points: impl Iterator<Item = Point>) -> Self {
            let mut vector = Self::zero();
            for point in points {
                vector.set(point, true);
            }
            vector
        }
    }

    #[derive(Debug, Clone)]
    pub struct OrderedSextet {
        foursomes: Labelled<hexacode::Point, Vector>,
    }

    impl OrderedSextet {
        pub fn point_foursomes(&self) -> Labelled<Point, hexacode::Point> {
            let mut labels: Labelled<Point, Option<hexacode::Point>> = Labelled::new_constant(None);
            for (h, foursome) in self.foursomes.iter() {
                for p in foursome.points() {
                    labels.set(p, Some(h));
                }
            }
            labels.apply_fn(|x| x.unwrap())
        }

        pub fn from_foursomes(foursomes: Labelled<hexacode::Point, Vector>) -> Self {
            let hexacode_points = hexacode::Point::points().collect::<Vec<_>>();
            for h in &hexacode_points {
                assert_eq!(foursomes.get(*h).weight(), 4);
            }
            for i in 0..hexacode_points.len() {
                for j in 0..i {
                    assert_eq!(
                        (foursomes.get(hexacode_points[i]) & foursomes.get(hexacode_points[j]))
                            .weight(),
                        0
                    );
                }
            }
            Self { foursomes }
        }

        pub fn foursome(&self, foursome: hexacode::Point) -> &Vector {
            self.foursomes.get(foursome)
        }

        pub fn permute(self, perm: &Permutation<hexacode::Point>) -> Self {
            Self {
                foursomes: self.foursomes.permute(perm),
            }
        }
    }

    // Each foursome of the sexet labelled with F4 defining an isomorphism with the MOG
    #[derive(Debug, Clone)]
    pub struct OrderedSextetLabelling {
        sextet: OrderedSextet,
        labels: Labelled<Point, F4Point>,
    }

    impl OrderedSextetLabelling {
        pub fn labels(&self) -> &Labelled<Point, F4Point> {
            &self.labels
        }

        pub fn foursomes(&self) -> Labelled<Point, hexacode::Point> {
            self.sextet.point_foursomes()
        }

        pub fn permute_foursomes(self, perm: &Permutation<hexacode::Point>) -> Self {
            Self {
                sextet: self.sextet.permute(perm),
                labels: self.labels,
            }
        }

        pub fn add_vector(self, vector: hexacode::Vector) -> Self {
            let point_foursomes = self.sextet.point_foursomes();
            Self {
                sextet: self.sextet,
                labels: Labelled::from_fn(|point: Point| {
                    let foursome = point_foursomes.get(point);
                    *self.labels.get(point) + *vector.get(*foursome)
                }),
            }
        }

        pub fn scalar_mul(self, lambda: F4Point) -> Self {
            assert_ne!(lambda, F4Point::Zero);
            Self {
                sextet: self.sextet,
                labels: self.labels.apply_fn(|value| {
                    // use lambda.inverse() here because we want to permute the points not the labels
                    *value * lambda.inverse().unwrap()
                }),
            }
        }

        pub fn conjugate(self) -> Self {
            Self {
                sextet: self.sextet,
                labels: self.labels.apply_fn(|value| value.conjugate()),
            }
        }
    }

    pub struct BinaryGolayCode {
        basis: Vec<Vector>,
        codewords: HashSet<Vector>,
    }

    impl Default for BinaryGolayCode {
        fn default() -> Self {
            // Construct a basis for the binary golay code on the 6x4 MOG grid
            let basis = {
                let mut basis = vec![];

                let col1 = (hexacode::Side::Left, hexacode::Pair::Left);
                for col2 in [
                    (hexacode::Side::Right, hexacode::Pair::Left),
                    (hexacode::Side::Left, hexacode::Pair::Middle),
                    (hexacode::Side::Right, hexacode::Pair::Middle),
                    (hexacode::Side::Left, hexacode::Pair::Right),
                    (hexacode::Side::Right, hexacode::Pair::Right),
                ] {
                    basis.push(Vector::from_fn(|p| {
                        (p.col.side, p.col.pair) == col1 || (p.col.side, p.col.pair) == col2
                    }));
                }

                for val in [F4Point::One, F4Point::Alpha, F4Point::Beta] {
                    basis.push(Vector::from_fn(|p| match p.col.pair {
                        hexacode::Pair::Left => match p.col.side {
                            hexacode::Side::Left => p.row != F4Point::Zero,
                            hexacode::Side::Right => p.row == F4Point::Zero,
                        },
                        hexacode::Pair::Middle | hexacode::Pair::Right => p.row == val,
                    }))
                }

                basis.push(Vector::from_fn(|p| match (p.col.side, p.col.pair) {
                    (hexacode::Side::Left, hexacode::Pair::Left) => p.row != F4Point::Zero,
                    (hexacode::Side::Right, hexacode::Pair::Left) => p.row == F4Point::One,
                    (hexacode::Side::Left, hexacode::Pair::Middle) => p.row == F4Point::Zero,
                    (hexacode::Side::Right, hexacode::Pair::Middle) => p.row == F4Point::One,
                    (hexacode::Side::Left, hexacode::Pair::Right) => p.row == F4Point::Alpha,
                    (hexacode::Side::Right, hexacode::Pair::Right) => p.row == F4Point::Beta,
                }));

                basis.push(Vector::from_fn(|p| match (p.col.side, p.col.pair) {
                    (hexacode::Side::Left, hexacode::Pair::Left) => p.row != F4Point::Zero,
                    (hexacode::Side::Right, hexacode::Pair::Left) => p.row == F4Point::Alpha,
                    (hexacode::Side::Left, hexacode::Pair::Middle) => p.row == F4Point::Zero,
                    (hexacode::Side::Right, hexacode::Pair::Middle) => p.row == F4Point::Alpha,
                    (hexacode::Side::Left, hexacode::Pair::Right) => p.row == F4Point::Beta,
                    (hexacode::Side::Right, hexacode::Pair::Right) => p.row == F4Point::One,
                }));

                basis.push(Vector::from_fn(|p| match (p.col.side, p.col.pair) {
                    (hexacode::Side::Left, hexacode::Pair::Left) => p.row == F4Point::One,
                    (hexacode::Side::Right, hexacode::Pair::Left) => p.row != F4Point::Zero,
                    (hexacode::Side::Left, hexacode::Pair::Middle) => p.row == F4Point::Zero,
                    (hexacode::Side::Right, hexacode::Pair::Middle) => p.row == F4Point::One,
                    (hexacode::Side::Left, hexacode::Pair::Right) => p.row == F4Point::Beta,
                    (hexacode::Side::Right, hexacode::Pair::Right) => p.row == F4Point::Alpha,
                }));

                basis.push(Vector::from_fn(|p| match (p.col.side, p.col.pair) {
                    (hexacode::Side::Left, hexacode::Pair::Left) => p.row == F4Point::Alpha,
                    (hexacode::Side::Right, hexacode::Pair::Left) => p.row != F4Point::Zero,
                    (hexacode::Side::Left, hexacode::Pair::Middle) => p.row == F4Point::Zero,
                    (hexacode::Side::Right, hexacode::Pair::Middle) => p.row == F4Point::Alpha,
                    (hexacode::Side::Left, hexacode::Pair::Right) => p.row == F4Point::One,
                    (hexacode::Side::Right, hexacode::Pair::Right) => p.row == F4Point::Beta,
                }));

                basis
            };

            // Compute the span of the basis to obtain the codewords in the binary golay code
            let mut codewords = HashSet::new();
            for b in 0usize..(1 << basis.len()) {
                let mut codeword = Vector::zero();
                #[allow(clippy::needless_range_loop)]
                for i in 0..basis.len() {
                    if b & (1usize << i) != 0 {
                        codeword = &codeword + &basis[i];
                    }
                }
                debug_assert!(!codewords.contains(&codeword));
                codewords.insert(codeword);
            }

            // Sanity checks
            for vector in &codewords {
                debug_assert!(
                    vector.weight() == 0
                        || vector.weight() == 8
                        || vector.weight() == 12
                        || vector.weight() == 16
                        || vector.weight() == 24
                );
            }
            debug_assert_eq!(basis.len(), 12);
            debug_assert_eq!(codewords.len(), 1usize << 12);

            Self { basis, codewords }
        }
    }

    impl BinaryGolayCode {
        pub fn is_codeword(&self, vector: &Vector) -> bool {
            self.codewords.contains(vector)
        }

        pub fn is_octad(&self, vector: &Vector) -> bool {
            vector.weight() == 8 && self.codewords.contains(vector)
        }

        // If the input vector has weight 5, return the unique octad containing it
        // Otherwise, return an Err
        pub fn complete_octad(&self, vector: &Vector) -> Result<Vector, ()> {
            if vector.weight() != 5 {
                return Err(());
            }
            for codeword in &self.codewords {
                if codeword.weight() == 8 && codeword.contains(vector) {
                    return Ok(codeword.clone());
                }
            }
            unreachable!()
        }

        pub fn complete_sextet(&self, vector: &Vector) -> Result<HashSet<Vector>, ()> {
            if vector.weight() != 4 {
                return Err(());
            }
            let mut sextet = HashSet::new();
            for codeword in &self.codewords {
                let offset = codeword + vector;
                if offset.weight() == 4 {
                    sextet.insert(offset);
                }
            }
            debug_assert_eq!(sextet.len(), 6);
            Ok(sextet)
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
        pub fn complete_labelling(
            &self,
            sextet: OrderedSextet,
            x: Point,
            y: Point,
            z: Point,
            w: Point,
            alpha: F4Point,
        ) -> OrderedSextetLabelling {
            // T1 T2  T3 T4  T5 T6
            let t0 = sextet.foursome(hexacode::Point::usize_to_point(0).unwrap());
            let t1 = sextet.foursome(hexacode::Point::usize_to_point(1).unwrap());
            let t2 = sextet.foursome(hexacode::Point::usize_to_point(2).unwrap());
            let t3 = sextet.foursome(hexacode::Point::usize_to_point(3).unwrap());
            let t4 = sextet.foursome(hexacode::Point::usize_to_point(4).unwrap());
            let t5 = sextet.foursome(hexacode::Point::usize_to_point(5).unwrap());

            // Validate inputs
            assert!(t0.contains_point(x));
            assert!(t1.contains_point(y));
            assert!(t1.contains_point(z));
            assert!(t2.contains_point(w));
            assert_ne!(y, z);
            let mut labels = Labelled::new_constant(F4Point::Zero);
            debug_assert_eq!(*labels.get(x), F4Point::Zero);
            debug_assert_eq!(*labels.get(y), F4Point::Zero);
            labels.set(z, F4Point::One);
            labels.set(w, alpha);

            let _ = t1; //It's not used

            #[allow(clippy::items_after_statements)]
            fn take_unique_pt(v: Vector) -> Point {
                let mut pts = v.points();
                let pt = pts.next().unwrap();
                debug_assert_eq!(pts.next(), None);
                pt
            }

            // Complete the hexacodeword (0, 1, alpha, ?, ?, ?) -> (0, 1, alpha, beta, gamma, delta)
            let (beta, gamma, delta) = match alpha {
                F4Point::Zero => (F4Point::One, F4Point::Alpha, F4Point::Beta),
                F4Point::One => (F4Point::Zero, F4Point::Beta, F4Point::Alpha),
                F4Point::Alpha => (F4Point::Beta, F4Point::Zero, F4Point::One),
                F4Point::Beta => (F4Point::Alpha, F4Point::One, F4Point::Zero),
            };
            // Use the octad containing (T1 \ {x}) U {z, w} to label 1 point in each of T3, T4, T5, T6
            let octad = self
                .complete_octad(&Vector::from_points(
                    t0.points().filter(|p| *p != x).chain(vec![z, w]),
                ))
                .unwrap();
            let (w2, w3, w4, w5) = (
                take_unique_pt(&octad & t2),
                take_unique_pt(&octad & t3),
                take_unique_pt(&octad & t4),
                take_unique_pt(&octad & t5),
            );
            debug_assert_eq!(w2, w);
            labels.set(w3, beta);
            labels.set(w4, gamma);
            labels.set(w5, delta);

            // Use the sextet formed by completing (T1 \ {x}) U {y} to label the rest of T3 U T4 U T5 U T6
            let t2345 = &(t2 | t3) | &(t4 | t5);
            debug_assert_eq!(t2345.weight(), 16);
            for (_i, li, wi) in [
                (2, alpha, w2),
                (3, beta, w3),
                (4, gamma, w4),
                (5, delta, w5),
            ] {
                let octad = self
                    .complete_octad(&Vector::from_points(
                        t0.points().filter(|p| *p != x).chain(vec![y, wi]),
                    ))
                    .unwrap();
                for p in (&octad & &t2345).points() {
                    labels.set(p, li);
                }
            }
            debug_assert_eq!(*labels.get(w2), alpha);
            debug_assert_eq!(*labels.get(w3), beta);
            debug_assert_eq!(*labels.get(w4), gamma);
            debug_assert_eq!(*labels.get(w5), delta);

            //Find the point labelled 0 in T4 and the three points not labelled 0 in T5
            let mut final_four = vec![];
            for p in t5.points() {
                if *labels.get(p) != F4Point::Zero {
                    final_four.push(p);
                }
            }
            for p in t4.points() {
                if *labels.get(p) == F4Point::Zero {
                    final_four.push(p);
                }
            }
            //Complete these 4 to an octad using each point in T3. The hexacodewords have labels llll00 so we can complete the labelling in T1 U T2
            debug_assert_eq!(final_four.len(), 4);
            for q in t3.points() {
                let l = *labels.get(q);
                let octad = self
                    .complete_octad(&Vector::from_points(
                        final_four.iter().copied().chain(vec![q]),
                    ))
                    .unwrap();
                for p in octad.points() {
                    if !t4.contains_point(p) && !t5.contains_point(p) {
                        labels.set(p, l);
                    }
                }
            }
            debug_assert_eq!(*labels.get(x), F4Point::Zero);
            debug_assert_eq!(*labels.get(y), F4Point::Zero);
            debug_assert_eq!(*labels.get(z), F4Point::One);
            debug_assert_eq!(*labels.get(w2), alpha);
            debug_assert_eq!(*labels.get(w3), beta);
            debug_assert_eq!(*labels.get(w4), gamma);
            debug_assert_eq!(*labels.get(w5), delta);

            #[cfg(debug_assertions)]
            {
                for t in [t0, t1, t2, t3, t4, t5] {
                    let mut zero_count: usize = 0;
                    let mut one_count: usize = 0;
                    let mut alpha_count: usize = 0;
                    let mut beta_count: usize = 0;
                    for p in t.points() {
                        match labels.get(p) {
                            F4Point::Zero => zero_count += 1,
                            F4Point::One => one_count += 1,
                            F4Point::Alpha => alpha_count += 1,
                            F4Point::Beta => beta_count += 1,
                        }
                    }
                    assert_eq!(zero_count, 1);
                    assert_eq!(one_count, 1);
                    assert_eq!(alpha_count, 1);
                    assert_eq!(beta_count, 1);
                }
            }

            OrderedSextetLabelling { sextet, labels }
        }
    }

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

    impl BinaryGolayCode {
        pub fn nearest_codeword(&self, vector: &Vector) -> NearestCodewordsResult {
            let mut dist_4_codewords = vec![];
            for codeword in &self.codewords {
                let diff = vector + codeword;
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
    }

    impl BinaryGolayCode {
        pub fn is_automorphism(&self, permutation: &Permutation<Point>) -> bool {
            self.basis
                .iter()
                .all(|b| self.codewords.contains(&b.permute(permutation)))
        }
    }
}
