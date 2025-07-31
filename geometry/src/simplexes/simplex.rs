use super::*;
use algebraeon_sets::structure::BorrowedStructure;
use itertools::Itertools;

#[derive(Clone)]
pub struct Simplex<
    FS: OrderedRingSignature + FieldSignature,
    FSB: BorrowedStructure<FS>,
    SP: Borrow<AffineSpace<FS, FSB>> + Clone,
> {
    ambient_space: SP,
    // points must be ordered w.r.t the ordering on vectors
    // points must be non-degenerate
    // points must belong to the ambient_space
    points: Vec<Vector<FS, FSB, SP>>,
}

#[allow(clippy::missing_fields_in_debug)]
impl<
    FS: OrderedRingSignature + FieldSignature,
    FSB: BorrowedStructure<FS>,
    SP: Borrow<AffineSpace<FS, FSB>> + Clone,
> std::fmt::Debug for Simplex<FS, FSB, SP>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Simplex")
            .field("points", &self.points)
            .finish()
    }
}

impl<
    FS: OrderedRingSignature + FieldSignature,
    FSB: BorrowedStructure<FS>,
    SP: Borrow<AffineSpace<FS, FSB>> + Clone,
> PartialEq for Simplex<FS, FSB, SP>
{
    fn eq(&self, other: &Self) -> bool {
        self.ambient_space.borrow() == other.ambient_space.borrow() && self.points == other.points
    }
}

impl<
    FS: OrderedRingSignature + FieldSignature,
    FSB: BorrowedStructure<FS>,
    SP: Borrow<AffineSpace<FS, FSB>> + Clone,
> Eq for Simplex<FS, FSB, SP>
{
}

impl<
    FS: OrderedRingSignature + FieldSignature,
    FSB: BorrowedStructure<FS>,
    SP: Borrow<AffineSpace<FS, FSB>> + Clone,
> Hash for Simplex<FS, FSB, SP>
where
    FS::Set: Hash,
{
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.points.hash(state);
    }
}

impl<
    FS: OrderedRingSignature + FieldSignature,
    FSB: BorrowedStructure<FS>,
    SP: Borrow<AffineSpace<FS, FSB>> + Clone,
> Simplex<FS, FSB, SP>
where
    SP: Clone,
{
    pub fn new(
        ambient_space: SP,
        mut points: Vec<Vector<FS, FSB, SP>>,
    ) -> Result<Self, &'static str> {
        for point in &points {
            assert_eq!(ambient_space.borrow(), point.ambient_space().borrow());
        }
        points.sort_unstable();
        if ambient_space
            .borrow()
            .are_points_affine_independent(points.iter().collect())
        {
            Ok(Self {
                ambient_space,
                points,
            })
        } else {
            Err("Can't make a simplex using degenerate points")
        }
    }

    pub fn ambient_space(&self) -> SP {
        self.ambient_space.clone()
    }

    pub fn n(&self) -> usize {
        self.points.len()
    }

    pub fn points(&self) -> &Vec<Vector<FS, FSB, SP>> {
        &self.points
    }

    pub fn into_points(self) -> Vec<Vector<FS, FSB, SP>> {
        self.points
    }

    pub fn point(&self, i: usize) -> &Vector<FS, FSB, SP> {
        &self.points[i]
    }

    pub fn skeleton(&self, skel_n: isize) -> Vec<Self> {
        if skel_n < 0 {
            vec![]
        } else {
            #[allow(clippy::cast_sign_loss)]
            let skel_n = skel_n as usize;
            let mut parts = vec![];
            for subset in (0..self.points.len()).combinations(skel_n) {
                let part = Self::new(
                    self.ambient_space.clone(),
                    subset.into_iter().map(|i| self.points[i].clone()).collect(),
                )
                .unwrap();
                parts.push(part);
            }
            parts
        }
    }

    pub fn vertices(&self) -> Vec<Self> {
        self.skeleton(1)
    }

    pub fn edges(&self) -> Vec<Self> {
        self.skeleton(2)
    }

    pub fn faces(&self) -> Vec<Self> {
        self.skeleton(3)
    }

    pub fn ridges(&self) -> Vec<Self> {
        #[allow(clippy::cast_possible_wrap)]
        self.skeleton(self.points.len() as isize - 2)
    }

    pub fn facets(&self) -> Vec<Self> {
        #[allow(clippy::cast_possible_wrap)]
        self.skeleton(self.points.len() as isize - 1)
    }

    pub fn facet(&self, k: usize) -> Self {
        assert!(k <= self.points.len());
        Self::new(self.ambient_space.clone(), {
            let mut facet_points = self.points.clone();
            facet_points.remove(k);
            facet_points
        })
        .unwrap()
    }

    pub fn sub_simplices(&self) -> Vec<Self> {
        self.points()
            .clone()
            .into_iter()
            .powerset()
            .map(|sub_points| Self::new(self.ambient_space.clone(), sub_points).unwrap())
            .collect()
    }

    pub fn sub_simplices_not_null(&self) -> Vec<Self> {
        self.sub_simplices()
            .into_iter()
            .filter(|spx| spx.n() != 0)
            .collect()
    }

    pub fn proper_sub_simplices_not_null(&self) -> Vec<Self> {
        self.sub_simplices()
            .into_iter()
            .filter(|spx| spx.n() != 0 && spx.n() != self.n())
            .collect()
    }

    #[allow(clippy::needless_pass_by_value)]
    pub fn sub_simplex(&self, pts: Vec<usize>) -> Self {
        Self::new(
            self.ambient_space(),
            pts.iter().map(|idx| self.points[*idx].clone()).collect(),
        )
        .unwrap()
    }

    pub fn oriented_facet(&self, k: usize) -> OrientedSimplex<FS, FSB, SP> {
        //return the oriented facet of self with negative side on the outside and positive side on the inside
        assert!(k <= self.points.len());
        let mut facet_points = self.points.clone();
        let other_pt = facet_points.remove(k);
        OrientedSimplex::new_with_positive_point(
            self.ambient_space.clone(),
            facet_points,
            &other_pt,
        )
        .unwrap()
    }

    pub fn oriented_facets(&self) -> Vec<OrientedSimplex<FS, FSB, SP>> {
        assert_eq!(self.ambient_space.borrow().affine_dimension(), self.n());
        (0..self.n()).map(|k| self.oriented_facet(k)).collect()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OrientationSide {
    Positive,
    Neutral,
    Negative,
}

/// a simplex spanning a hyperplane with a positive side and a negative side
/// in 3d space it is a triangle
/// in 2d space it is a line
/// in 1d space it is a point
/// in 0d space it is a null simplex. The orientation looses meaning but it is nice to still count this case.
#[derive(Clone)]
struct OrientedSimplexOrientation<
    FS: OrderedRingSignature + FieldSignature,
    FSB: BorrowedStructure<FS>,
    SP: Borrow<AffineSpace<FS, FSB>> + Clone,
> {
    flip: bool,                          // flip the orientation if necessary
    positive_point: Vector<FS, FSB, SP>, // a point on the positive side
}

impl<
    FS: OrderedRingSignature + FieldSignature,
    FSB: BorrowedStructure<FS>,
    SP: Borrow<AffineSpace<FS, FSB>> + Clone,
> std::fmt::Debug for OrientedSimplexOrientation<FS, FSB, SP>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("OrientedSimplexOrientation")
            .field("flip", &self.flip)
            .field("positive_point", &self.positive_point)
            .finish()
    }
}

#[derive(Clone)]
pub struct OrientedSimplex<
    FS: OrderedRingSignature + FieldSignature,
    FSB: BorrowedStructure<FS>,
    SP: Borrow<AffineSpace<FS, FSB>> + Clone,
> {
    simplex: Simplex<FS, FSB, SP>,
    orientation: Option<OrientedSimplexOrientation<FS, FSB, SP>>, //None iff simplex is null
}

impl<
    FS: OrderedRingSignature + FieldSignature,
    FSB: BorrowedStructure<FS>,
    SP: Borrow<AffineSpace<FS, FSB>> + Clone,
> std::fmt::Debug for OrientedSimplex<FS, FSB, SP>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("OrientedSimplex")
            .field("simplex", &self.simplex)
            .field("orientation", &self.orientation)
            .finish()
    }
}

impl<
    FS: OrderedRingSignature + FieldSignature,
    FSB: BorrowedStructure<FS>,
    SP: Borrow<AffineSpace<FS, FSB>> + Clone,
> OrientedSimplex<FS, FSB, SP>
{
    pub fn new_with_positive_point(
        ambient_space: SP,
        points: Vec<Vector<FS, FSB, SP>>,
        ref_point: &Vector<FS, FSB, SP>,
    ) -> Result<Self, &'static str> {
        assert_eq!(ref_point.ambient_space().borrow(), ambient_space.borrow());
        if points.len() != ambient_space.borrow().linear_dimension().unwrap() {
            return Err("Oriented simplex must have dimension one less than the ambient space");
        }
        let n = points.len();
        let simplex = Simplex::new(ambient_space, points)?;
        if n == 0 {
            Ok(Self {
                simplex,
                orientation: None,
            })
        } else {
            let mut guess = Self {
                simplex,
                orientation: Some(OrientedSimplexOrientation {
                    flip: false,
                    positive_point: ref_point.clone(),
                }),
            };
            match guess.classify_point(ref_point) {
                OrientationSide::Positive => Ok(guess),
                OrientationSide::Neutral => Err("The reference point lies on the hyperplane"),
                OrientationSide::Negative => {
                    let orientation = guess.orientation.as_mut().unwrap();
                    orientation.flip = !orientation.flip;
                    Ok(guess)
                }
            }
        }
    }

    pub fn new_with_negative_point(
        ambient_space: SP,
        points: Vec<Vector<FS, FSB, SP>>,
        ref_point: &Vector<FS, FSB, SP>,
    ) -> Result<Self, &'static str> {
        let mut ans = Self::new_with_positive_point(ambient_space, points, ref_point)?;
        ans.flip();
        Ok(ans)
    }

    pub fn positive_point(&self) -> Option<Vector<FS, FSB, SP>> {
        Some(self.orientation.as_ref()?.positive_point.clone())
    }

    pub fn negative_point(&self) -> Option<Vector<FS, FSB, SP>> {
        let positive_point = self.positive_point()?;
        Some({
            let pt: &Vector<FS, FSB, SP> = &self.simplex.points()[0];
            &(pt + pt) - &positive_point
        })
    }

    pub fn ambient_space(&self) -> SP {
        self.simplex().ambient_space()
    }

    pub fn simplex(&self) -> &Simplex<FS, FSB, SP> {
        &self.simplex
    }

    pub fn flip(&mut self) {
        let negative_point = self.negative_point();
        if let Some(OrientedSimplexOrientation {
            flip,
            positive_point,
        }) = &mut self.orientation
        {
            (*flip, *positive_point) = (!*flip, negative_point.unwrap());
        }
        if let Some(OrientedSimplexOrientation {
            flip: _flip,
            positive_point,
        }) = &self.orientation
        {
            debug_assert_eq!(
                self.classify_point(positive_point),
                OrientationSide::Positive
            );
        }
    }

    fn classify_point_quantitatively(&self, point: &Vector<FS, FSB, SP>) -> FS::Set {
        let space = self.ambient_space();
        let field = space.borrow().field();
        match &self.orientation {
            Some(OrientedSimplexOrientation {
                flip,
                positive_point: _,
            }) => match space.borrow().linear_dimension().unwrap() {
                0 => unreachable!(),
                d => {
                    let root = &self.simplex.points[0];
                    let mut vecs = (1..d)
                        .map(|i| &self.simplex.points[i] - root)
                        .collect::<Vec<_>>();
                    vecs.push(point - root);
                    let det = space.borrow().determinant(vecs.iter().collect());
                    match flip {
                        true => field.neg(&det),
                        false => det,
                    }
                }
            },
            None => field.zero(),
        }
    }

    pub fn classify_point(&self, point: &Vector<FS, FSB, SP>) -> OrientationSide {
        let space = self.ambient_space();
        let field = space.borrow().field();
        let value = self.classify_point_quantitatively(point);
        match field.ring_cmp(&value, &field.zero()) {
            std::cmp::Ordering::Less => OrientationSide::Negative,
            std::cmp::Ordering::Equal => OrientationSide::Neutral,
            std::cmp::Ordering::Greater => OrientationSide::Positive,
        }
    }

    pub fn into_oriented_hyperplane(self) -> OrientedHyperplane<FS, FSB, SP> {
        OrientedHyperplane {
            oriented_simplex: self,
        }
    }
}

#[derive(Clone)]
pub struct OrientedHyperplane<
    FS: OrderedRingSignature + FieldSignature,
    FSB: BorrowedStructure<FS>,
    SP: Borrow<AffineSpace<FS, FSB>> + Clone,
> {
    oriented_simplex: OrientedSimplex<FS, FSB, SP>,
}

impl<
    FS: OrderedRingSignature + FieldSignature,
    FSB: BorrowedStructure<FS>,
    SP: Borrow<AffineSpace<FS, FSB>> + Clone,
> std::fmt::Debug for OrientedHyperplane<FS, FSB, SP>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("OrientedHyperplane")
            .field("oriented_simplex", &self.oriented_simplex)
            .finish()
    }
}

// impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone>
//     From<OrientedSimplex<FS, SP>> for OrientedHyperplane<FS, SP>
// {
//     fn from(oriented_simplex: OrientedSimplex<FS, SP>) -> Self {
//         Self { oriented_simplex }
//     }
// }

pub enum OrientedHyperplaneIntersectLineSegmentResult<
    FS: OrderedRingSignature + FieldSignature,
    FSB: BorrowedStructure<FS>,
    SP: Borrow<AffineSpace<FS, FSB>> + Clone,
> {
    PositivePositive,
    NegativeNegative,
    PositiveNeutral,
    NegativeNeutral,
    NeutralPositive,
    NeutralNegative,
    NeutralNeutral,
    PositiveNegative {
        intersection_point: Vector<FS, FSB, SP>,
    },
    NegativePositive {
        intersection_point: Vector<FS, FSB, SP>,
    },
}

impl<
    FS: OrderedRingSignature + FieldSignature,
    FSB: BorrowedStructure<FS>,
    SP: Borrow<AffineSpace<FS, FSB>> + Clone,
> OrientedHyperplane<FS, FSB, SP>
{
    pub fn ambient_space(&self) -> SP {
        self.oriented_simplex.ambient_space()
    }

    pub fn classify_point(&self, point: &Vector<FS, FSB, SP>) -> OrientationSide {
        self.oriented_simplex.classify_point(point)
    }

    pub fn intersect_line(
        &self,
        a: &Vector<FS, FSB, SP>,
        b: &Vector<FS, FSB, SP>,
    ) -> OrientedHyperplaneIntersectLineSegmentResult<FS, FSB, SP> {
        let space = self.ambient_space();
        let field = space.borrow().field();

        let a_val = self.oriented_simplex.classify_point_quantitatively(a);
        let b_val = self.oriented_simplex.classify_point_quantitatively(b);

        match (
            field.ring_cmp(&a_val, &field.zero()),
            field.ring_cmp(&b_val, &field.zero()),
        ) {
            (std::cmp::Ordering::Less, std::cmp::Ordering::Greater) => {
                let t = field
                    .div(&a_val, &field.add(&a_val, &field.neg(&b_val)))
                    .unwrap();
                {
                    debug_assert_eq!(
                        self.oriented_simplex.classify_point(a),
                        OrientationSide::Negative
                    );
                    debug_assert_eq!(
                        self.oriented_simplex.classify_point(b),
                        OrientationSide::Positive
                    );
                    OrientedHyperplaneIntersectLineSegmentResult::NegativePositive {
                        intersection_point: a + &(b - a).scalar_mul(&t),
                    }
                }
            }
            (std::cmp::Ordering::Greater, std::cmp::Ordering::Less) => {
                let t = field
                    .div(&a_val, &field.add(&a_val, &field.neg(&b_val)))
                    .unwrap();
                {
                    debug_assert_eq!(
                        self.oriented_simplex.classify_point(a),
                        OrientationSide::Positive
                    );
                    debug_assert_eq!(
                        self.oriented_simplex.classify_point(b),
                        OrientationSide::Negative
                    );
                    OrientedHyperplaneIntersectLineSegmentResult::PositiveNegative {
                        intersection_point: a + &(b - a).scalar_mul(&t),
                    }
                }
            }
            (std::cmp::Ordering::Less, std::cmp::Ordering::Less) => {
                OrientedHyperplaneIntersectLineSegmentResult::NegativeNegative
            }
            (std::cmp::Ordering::Less, std::cmp::Ordering::Equal) => {
                OrientedHyperplaneIntersectLineSegmentResult::NegativeNeutral
            }
            (std::cmp::Ordering::Equal, std::cmp::Ordering::Less) => {
                OrientedHyperplaneIntersectLineSegmentResult::NeutralNegative
            }
            (std::cmp::Ordering::Equal, std::cmp::Ordering::Equal) => {
                OrientedHyperplaneIntersectLineSegmentResult::NeutralNeutral
            }
            (std::cmp::Ordering::Equal, std::cmp::Ordering::Greater) => {
                OrientedHyperplaneIntersectLineSegmentResult::NeutralPositive
            }
            (std::cmp::Ordering::Greater, std::cmp::Ordering::Equal) => {
                OrientedHyperplaneIntersectLineSegmentResult::PositiveNeutral
            }
            (std::cmp::Ordering::Greater, std::cmp::Ordering::Greater) => {
                OrientedHyperplaneIntersectLineSegmentResult::PositivePositive
            }
        }
    }

    pub fn flip(&mut self) {
        self.oriented_simplex.flip();
    }
}

/*
pub fn simplex_intersect_with_hyperplane<
    FS: OrderedRingStructure + FieldStructure,
    SP: Borrow<AffineSpace<FS>> + Clone,
>(
    simplex: &Simplex<FS, SP>,
    hyperplane: &OrientedSimplex<FS, SP>,
) -> Vec<Simplex<FS, SP>> {
    todo!()
}

pub fn simplex_intersect_positive_side_hyperplane<
    FS: OrderedRingStructure + FieldStructure,
    SP: Borrow<AffineSpace<FS>> + Clone,
>(
    simplex: &Simplex<FS, SP>,
    hyperplane: &OrientedSimplex<FS, SP>,
) -> Vec<Simplex<FS, SP>> {
    todo!()
}

pub fn simplex_intersect_negative_side_hyperplane<
    FS: OrderedRingStructure + FieldStructure,
    SP: Borrow<AffineSpace<FS>> + Clone,
>(
    simplex: &Simplex<FS, SP>,
    hyperplane: &OrientedSimplex<FS, SP>,
) -> Vec<Simplex<FS, SP>> {
    simplex_intersect_positive_side_hyperplane(simplex, &hyperplane.clone().flip())
}
*/

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_nzq::Rational;
    use algebraeon_sets::structure::*;

    #[test]
    fn make_simplex() {
        let space = AffineSpace::new_linear(Rational::structure(), 2);
        let v1 = Vector::new(&space, vec![Rational::from(1), Rational::from(1)]);
        let v2 = Vector::new(&space, vec![Rational::from(1), Rational::from(0)]);
        let v3 = Vector::new(&space, vec![Rational::from(0), Rational::from(1)]);
        let s = Simplex::new(&space, vec![v1, v2, v3]);
        assert!(s.is_ok());

        let space = AffineSpace::new_linear(Rational::structure(), 2);
        let v1 = Vector::new(&space, vec![Rational::from(0), Rational::from(0)]);
        let v2 = Vector::new(&space, vec![Rational::from(1), Rational::from(0)]);
        let v3 = Vector::new(&space, vec![Rational::from(2), Rational::from(0)]);
        let s = Simplex::new(&space, vec![v1, v2, v3]);
        assert!(s.is_err());
    }

    #[test]
    fn simplex_skeleton() {
        let space = AffineSpace::new_linear(Rational::structure(), 2);
        let v1 = Vector::new(&space, vec![Rational::from(1), Rational::from(1)]);
        let v2 = Vector::new(&space, vec![Rational::from(1), Rational::from(0)]);
        let v3 = Vector::new(&space, vec![Rational::from(0), Rational::from(1)]);
        let s = Simplex::new(&space, vec![v1, v2, v3]).unwrap();

        assert_eq!(s.skeleton(-2).len(), 0);
        assert_eq!(s.skeleton(-1).len(), 0);
        assert_eq!(s.skeleton(0).len(), 1);
        assert_eq!(s.vertices().len(), 3);
        assert_eq!(s.edges().len(), 3);
        assert_eq!(s.faces().len(), 1);
        assert_eq!(s.skeleton(4).len(), 0);
    }

    #[test]
    fn make_oriented_simplex() {
        let space = AffineSpace::new_linear(Rational::structure(), 2);
        let v1 = Vector::new(&space, vec![Rational::from(1), Rational::from(0)]);
        let v2 = Vector::new(&space, vec![Rational::from(0), Rational::from(1)]);
        let v3 = Vector::new(&space, vec![Rational::from(2), Rational::from(3)]);
        let s_pos =
            OrientedSimplex::new_with_positive_point(&space, vec![v1.clone(), v2.clone()], &v3)
                .unwrap();
        let s_neg = OrientedSimplex::new_with_negative_point(&space, vec![v1, v2], &v3).unwrap();

        assert_ne!(
            s_pos.orientation.unwrap().flip,
            s_neg.orientation.unwrap().flip
        );
    }
}
