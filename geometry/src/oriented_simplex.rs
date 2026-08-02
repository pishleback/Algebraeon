use crate::{
    ambient_space::AffineSpace,
    simplex::Simplex,
    vector::{DotProduct, Vector},
};
use algebraeon_rings::{
    matrix::{Matrix, MatrixStructure},
    structure::{FieldSignature, OrderedRingSignature},
};

#[derive(Clone)]
pub struct OrientedSimplex<FS: OrderedRingSignature + FieldSignature> {
    simplex: Simplex<FS>,
    orientation: Option<OrientedSimplexOrientation<FS>>, //None iff simplex is null
}

impl<FS: OrderedRingSignature + FieldSignature> std::fmt::Debug for OrientedSimplex<FS> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("OrientedSimplex")
            .field("simplex", &self.simplex)
            .field("orientation", &self.orientation)
            .finish()
    }
}

impl<FS: OrderedRingSignature + FieldSignature> OrientedSimplex<FS> {
    pub fn new_with_positive_point(
        ambient_space: &AffineSpace<FS>,
        points: Vec<Vector<FS>>,
        ref_point: &Vector<FS>,
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
            let positive_normal = {
                let mat = Matrix::construct(n - 1, n, |r, c| {
                    ambient_space
                        .field()
                        .sub(points[r + 1].coordinate(c), points[0].coordinate(c))
                });
                let kernel = MatrixStructure::<FS>::new(ambient_space.field().clone())
                    .col_kernel(mat)
                    .basis();
                if kernel.len() != 1 {
                    return Err("points are not affine independent");
                }
                ambient_space.vector(kernel.into_iter().next().unwrap())
            };
            let flip = match ambient_space.field().cmp(
                &(ref_point - root).dot(&positive_normal),
                &ambient_space.field().zero(),
            ) {
                std::cmp::Ordering::Less => true,
                std::cmp::Ordering::Equal => {
                    return Err("ref_point lines inside the hyperplane");
                }
                std::cmp::Ordering::Greater => false,
            };
            let plane_point = root.clone();

            Ok(Self {
                simplex: ambient_space.simplex(points)?,
                orientation: Some(OrientedSimplexOrientation {
                    flip,
                    plane_point,
                    positive_normal,
                }),
            })
        }
    }

    pub fn new_with_negative_point(
        ambient_space: &AffineSpace<FS>,
        points: Vec<Vector<FS>>,
        ref_point: &Vector<FS>,
    ) -> Result<Self, &'static str> {
        let mut ans = Self::new_with_positive_point(ambient_space, points, ref_point)?;
        ans.flip();
        Ok(ans)
    }

    pub fn positive_point(&self) -> Option<Vector<FS>> {
        if let Some(orientation) = self.orientation.as_ref() {
            match orientation.flip {
                false => Some(&orientation.plane_point + &orientation.positive_normal),
                true => Some(&orientation.plane_point - &orientation.positive_normal),
            }
        } else {
            None
        }
    }

    pub fn negative_point(&self) -> Option<Vector<FS>> {
        let positive_point = self.positive_point()?;
        Some({
            let pt: &Vector<FS> = &self.simplex.points()[0];
            &(pt + pt) - &positive_point
        })
    }

    pub fn ambient_space(&self) -> &AffineSpace<FS> {
        self.simplex().ambient_space()
    }

    pub fn simplex(&self) -> &Simplex<FS> {
        &self.simplex
    }

    pub fn into_simplex(self) -> Simplex<FS> {
        self.simplex
    }

    pub fn into_oriented_hyperplane(self) -> OrientedHyperplane<FS> {
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

    fn classify_point_quantitatively(&self, point: &Vector<FS>) -> FS::Elem {
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

    pub fn classify_point(&self, point: &Vector<FS>) -> OrientationSide {
        let space = self.ambient_space();
        let field = space.field();
        let value = self.classify_point_quantitatively(point);
        match field.cmp(&value, &field.zero()) {
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
struct OrientedSimplexOrientation<FS: OrderedRingSignature + FieldSignature> {
    flip: bool,                  // flip the orientation if necessary
    plane_point: Vector<FS>,     // A point on the hyperplane
    positive_normal: Vector<FS>, // normal vector pointing out of the positive side
}

#[derive(Clone)]
pub struct OrientedHyperplane<FS: OrderedRingSignature + FieldSignature> {
    oriented_simplex: OrientedSimplex<FS>,
}

impl<FS: OrderedRingSignature + FieldSignature> std::fmt::Debug for OrientedHyperplane<FS> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("OrientedHyperplane")
            .field("oriented_simplex", &self.oriented_simplex)
            .finish()
    }
}

// impl<FS: OrderedRingStructure + FieldStructure, AffineSpace< FS>: Borrow<AffineSpace<FS>> + Clone>
//     From<OrientedSimplex<FS>> for OrientedHyperplane<FS>
// {
//     fn from(oriented_simplex: OrientedSimplex<FS>) -> Self {
//         Self { oriented_simplex }
//     }
// }

pub enum OrientedHyperplaneIntersectLineSegmentResult<FS: OrderedRingSignature + FieldSignature> {
    PositivePositive,
    NegativeNegative,
    PositiveNeutral,
    NegativeNeutral,
    NeutralPositive,
    NeutralNegative,
    NeutralNeutral,
    PositiveNegative { intersection_point: Vector<FS> },
    NegativePositive { intersection_point: Vector<FS> },
}

impl<FS: OrderedRingSignature + FieldSignature> OrientedHyperplane<FS> {
    pub fn ambient_space(&self) -> &AffineSpace<FS> {
        self.oriented_simplex.ambient_space()
    }

    pub fn classify_point(&self, point: &Vector<FS>) -> OrientationSide {
        self.oriented_simplex.classify_point(point)
    }

    pub fn intersect_line(
        &self,
        a: &Vector<FS>,
        b: &Vector<FS>,
    ) -> OrientedHyperplaneIntersectLineSegmentResult<FS> {
        let space = self.ambient_space();
        let field = space.field();

        let a_val = self.oriented_simplex.classify_point_quantitatively(a);
        let b_val = self.oriented_simplex.classify_point_quantitatively(b);

        match (
            field.cmp(&a_val, &field.zero()),
            field.cmp(&b_val, &field.zero()),
        ) {
            (std::cmp::Ordering::Less, std::cmp::Ordering::Greater) => {
                let t = field
                    .try_divide(&a_val, &field.add(&a_val, &field.neg(&b_val)))
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
                    .try_divide(&a_val, &field.add(&a_val, &field.neg(&b_val)))
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
    use super::*;
    use algebraeon_structures::{MetaType, Rational};

    #[test]
    fn make_oriented_simplex() {
        let space = AffineSpace::new_linear(Rational::structure(), 2);
        let v1 = space.vector([1, 0]);
        let v2 = space.vector([0, 1]);
        let v3 = space.vector([2, 3]);
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
