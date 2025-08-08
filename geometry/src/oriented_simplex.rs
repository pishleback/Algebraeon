use crate::{
    ambient_space::AffineSpace,
    simplex::Simplex,
    vector::{DotProduct, Vector},
};
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
        assert_eq!(ref_point.ambient_space(), ambient_space);
        if points.len() != ambient_space.linear_dimension().unwrap() {
            return Err("Oriented simplex must have dimension one less than the ambient space");
        }
        let n = points.len();
        if n == 0 {
            Ok(Self {
                simplex: ambient_space.simplex(points)?,
                orientation: None,
            })
        } else {
            let root = &points[0];
            let ref_vec = ref_point - root;
            let ref_normal = {
                let mut ref_normal = ref_vec;
                #[allow(clippy::needless_range_loop)]
                for i in 1..n {
                    let vec = &points[i] - root;
                    ref_normal = &ref_normal
                        - &vec.scalar_mul(
                            &ambient_space
                                .field()
                                .div(&vec.dot(&ref_normal), &vec.dot(&vec))
                                .unwrap(),
                        );
                }
                ref_normal
            };

            debug_assert!(!(0..n).all(|i| ambient_space.field().is_zero(ref_normal.coordinate(i))));
            #[allow(clippy::needless_range_loop)]
            for i in 1..n {
                debug_assert!(
                    ambient_space
                        .field()
                        .is_zero(&(&points[i] - root).dot(&ref_normal))
                );
            }

            let plane_point = root.clone();

            let mut guess = Self {
                simplex: ambient_space.simplex(points)?,
                orientation: Some(OrientedSimplexOrientation {
                    flip: false,
                    plane_point,
                    positive_normal: ref_normal,
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
        if let Some(orientation) = self.orientation.as_ref() {
            match orientation.flip {
                false => Some(&orientation.plane_point + &orientation.positive_normal),
                true => Some(&orientation.plane_point - &orientation.positive_normal),
            }
        } else {
            None
        }
    }

    pub fn negative_point(&self) -> Option<Vector<'f, FS>> {
        let positive_point = self.positive_point()?;
        Some({
            let pt: &Vector<'f, FS> = &self.simplex.points()[0];
            &(pt + pt) - &positive_point
        })
    }

    pub fn ambient_space(&self) -> AffineSpace<'f, FS> {
        self.simplex().ambient_space()
    }

    pub fn simplex(&self) -> &Simplex<'f, FS> {
        &self.simplex
    }

    pub fn into_simplex(self) -> Simplex<'f, FS> {
        self.simplex
    }

    pub fn into_oriented_hyperplane(self) -> OrientedHyperplane<'f, FS> {
        OrientedHyperplane {
            oriented_simplex: self,
        }
    }

    pub fn flip(&mut self) {
        if let Some(OrientedSimplexOrientation { flip, .. }) = &mut self.orientation {
            *flip = !*flip;
        }
        if let Some(OrientedSimplexOrientation {
            flip,
            plane_point,
            positive_normal,
        }) = &self.orientation
        {
            debug_assert_eq!(
                self.classify_point(&(plane_point + positive_normal)),
                match flip {
                    false => OrientationSide::Positive,
                    true => OrientationSide::Negative,
                }
            );
        }
    }

    fn classify_point_quantitatively(&self, point: &Vector<'f, FS>) -> FS::Set {
        match &self.orientation {
            Some(OrientedSimplexOrientation {
                flip,
                plane_point,
                positive_normal,
            }) => {
                let mut value = positive_normal.dot(&(point - plane_point));
                if *flip {
                    value = self.ambient_space().field().neg(&value);
                }
                value
            }
            None => self.ambient_space().field().zero(),
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
#[derive(Debug, Clone)]
struct OrientedSimplexOrientation<'f, FS: OrderedRingSignature + FieldSignature> {
    flip: bool,                      // flip the orientation if necessary
    plane_point: Vector<'f, FS>,     // A point on the hyperplane
    positive_normal: Vector<'f, FS>, // normal vector pointing out of the positive side
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
    pub fn ambient_space(&self) -> AffineSpace<'f, FS> {
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
        let v1 = space.vector([1, 0]);
        let v2 = space.vector([0, 1]);
        let v3 = space.vector([2, 3]);
        let s_pos =
            OrientedSimplex::new_with_positive_point(space, vec![v1.clone(), v2.clone()], &v3)
                .unwrap();
        let s_neg = OrientedSimplex::new_with_negative_point(space, vec![v1, v2], &v3).unwrap();

        assert_ne!(
            s_pos.orientation.unwrap().flip,
            s_neg.orientation.unwrap().flip
        );
    }
}
