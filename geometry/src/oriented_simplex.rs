use crate::{ambient_space::AffineSpace, coordinates::Vector, simplex::Simplex};
use algebraeon_rings::structure::{FieldSignature, OrderedRingSignature};

#[derive(Clone)]
pub struct OrientedSimplex<'f, FS: OrderedRingSignature + FieldSignature> {
    simplex: Simplex<'f, FS>,
    orientation: Option<OrientedSimplexOrientation<'f, FS>>, //None iff simplex is null
}

impl<'f, FS: OrderedRingSignature + FieldSignature> std::fmt::Debug for OrientedSimplex<'f, FS> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("OrientedSimplex")
            .field("simplex", &self.simplex)
            .field("orientation", &self.orientation)
            .finish()
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> OrientedSimplex<'f, FS> {
    pub fn new_with_positive_point(
        ambient_space: AffineSpace<'f, FS>,
        points: Vec<Vector<'f, FS>>,
        ref_point: &Vector<'f, FS>,
    ) -> Result<Self, &'static str> {
        assert_eq!(ref_point.ambient_space(), &ambient_space);
        if points.len() != ambient_space.linear_dimension().unwrap() {
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
        ambient_space: AffineSpace<'f, FS>,
        points: Vec<Vector<'f, FS>>,
        ref_point: &Vector<'f, FS>,
    ) -> Result<Self, &'static str> {
        let mut ans = Self::new_with_positive_point(ambient_space, points, ref_point)?;
        ans.flip();
        Ok(ans)
    }

    pub fn positive_point(&self) -> Option<Vector<'f, FS>> {
        Some(self.orientation.as_ref()?.positive_point.clone())
    }

    pub fn negative_point(&self) -> Option<Vector<'f, FS>> {
        let positive_point = self.positive_point()?;
        Some({
            let pt: &Vector<'f, FS> = &self.simplex.points()[0];
            &(pt + pt) - &positive_point
        })
    }

    pub fn ambient_space(&self) -> &AffineSpace<'f, FS> {
        self.simplex().ambient_space()
    }

    pub fn simplex(&self) -> &Simplex<'f, FS> {
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

    fn classify_point_quantitatively(&self, point: &Vector<'f, FS>) -> FS::Set {
        let space = self.ambient_space();
        let field = space.field();
        match &self.orientation {
            Some(OrientedSimplexOrientation {
                flip,
                positive_point: _,
            }) => match space.linear_dimension().unwrap() {
                0 => unreachable!(),
                d => {
                    let root = &self.simplex.points()[0];
                    let mut vecs = (1..d)
                        .map(|i| self.simplex.point(i) - root)
                        .collect::<Vec<_>>();
                    vecs.push(point - root);
                    let det = space.determinant(vecs.iter().collect());
                    match flip {
                        true => field.neg(&det),
                        false => det,
                    }
                }
            },
            None => field.zero(),
        }
    }

    pub fn classify_point(&self, point: &Vector<'f, FS>) -> OrientationSide {
        let space = self.ambient_space();
        let field = space.field();
        let value = self.classify_point_quantitatively(point);
        match field.ring_cmp(&value, &field.zero()) {
            std::cmp::Ordering::Less => OrientationSide::Negative,
            std::cmp::Ordering::Equal => OrientationSide::Neutral,
            std::cmp::Ordering::Greater => OrientationSide::Positive,
        }
    }

    pub fn into_oriented_hyperplane(self) -> OrientedHyperplane<'f, FS> {
        OrientedHyperplane {
            oriented_simplex: self,
        }
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
struct OrientedSimplexOrientation<'f, FS: OrderedRingSignature + FieldSignature> {
    flip: bool,                     // flip the orientation if necessary
    positive_point: Vector<'f, FS>, // a point on the positive side
}

impl<'f, FS: OrderedRingSignature + FieldSignature> std::fmt::Debug
    for OrientedSimplexOrientation<'f, FS>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("OrientedSimplexOrientation")
            .field("flip", &self.flip)
            .field("positive_point", &self.positive_point)
            .finish()
    }
}

#[derive(Clone)]
pub struct OrientedHyperplane<'f, FS: OrderedRingSignature + FieldSignature> {
    oriented_simplex: OrientedSimplex<'f, FS>,
}

impl<'f, FS: OrderedRingSignature + FieldSignature> std::fmt::Debug for OrientedHyperplane<'f, FS> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("OrientedHyperplane")
            .field("oriented_simplex", &self.oriented_simplex)
            .finish()
    }
}

// impl<FS: OrderedRingStructure + FieldStructure, AffineSpace<'f, FS>: Borrow<AffineSpace<FS>> + Clone>
//     From<OrientedSimplex<FS>> for OrientedHyperplane<FS>
// {
//     fn from(oriented_simplex: OrientedSimplex<FS>) -> Self {
//         Self { oriented_simplex }
//     }
// }

pub enum OrientedHyperplaneIntersectLineSegmentResult<'f, FS: OrderedRingSignature + FieldSignature>
{
    PositivePositive,
    NegativeNegative,
    PositiveNeutral,
    NegativeNeutral,
    NeutralPositive,
    NeutralNegative,
    NeutralNeutral,
    PositiveNegative { intersection_point: Vector<'f, FS> },
    NegativePositive { intersection_point: Vector<'f, FS> },
}

impl<'f, FS: OrderedRingSignature + FieldSignature> OrientedHyperplane<'f, FS> {
    pub fn ambient_space(&self) -> &AffineSpace<'f, FS> {
        self.oriented_simplex.ambient_space()
    }

    pub fn classify_point(&self, point: &Vector<'f, FS>) -> OrientationSide {
        self.oriented_simplex.classify_point(point)
    }

    pub fn intersect_line(
        &self,
        a: &Vector<'f, FS>,
        b: &Vector<'f, FS>,
    ) -> OrientedHyperplaneIntersectLineSegmentResult<'f, FS> {
        let space = self.ambient_space();
        let field = space.field();

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

#[cfg(test)]
mod tests {
    use algebraeon_nzq::Rational;

    use super::*;

    #[test]
    fn make_oriented_simplex() {
        let space = AffineSpace::new_linear(Rational::structure_ref(), 2);
        let v1 = Vector::new(space.clone(), vec![Rational::from(1), Rational::from(0)]);
        let v2 = Vector::new(space.clone(), vec![Rational::from(0), Rational::from(1)]);
        let v3 = Vector::new(space.clone(), vec![Rational::from(2), Rational::from(3)]);
        let s_pos = OrientedSimplex::new_with_positive_point(
            space.clone(),
            vec![v1.clone(), v2.clone()],
            &v3,
        )
        .unwrap();
        let s_neg =
            OrientedSimplex::new_with_negative_point(space.clone(), vec![v1, v2], &v3).unwrap();

        assert_ne!(
            s_pos.orientation.unwrap().flip,
            s_neg.orientation.unwrap().flip
        );
    }
}
