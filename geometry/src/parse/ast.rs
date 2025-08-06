use crate::{
    ambient_space::AffineSpace,
    boolean_operations::{Difference, Intersect, Union},
    partial_simplicial_complex::PartialSimplicialComplex,
    simplex_collection::{InteriorOrBoundarySimplexCollection, LabelledSimplexCollection},
    vector::Vector,
};
use algebraeon_nzq::{Rational, RationalCanonicalStructure};
use std::str::FromStr;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Sign {
    Negative,
    Positive,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SignedValue {
    pub sign: Sign,
    pub value: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Point {
    pub coordinates: Vec<SignedValue>,
}

impl Point {
    fn dimension(&self) -> usize {
        self.coordinates.len()
    }

    fn to_vector(
        &self,
        space: AffineSpace<'static, RationalCanonicalStructure>,
    ) -> Vector<'static, RationalCanonicalStructure> {
        space.vector(self.coordinates.iter().map(|ast_value| {
            let mut value = Rational::from_str(&ast_value.value).unwrap();
            match ast_value.sign {
                Sign::Negative => value = -value,
                Sign::Positive => {}
            }
            value
        }))
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ShapeExpression {
    Union(Box<ShapeExpression>, Box<ShapeExpression>),
    Intersect(Box<ShapeExpression>, Box<ShapeExpression>),
    Difference(Box<ShapeExpression>, Box<ShapeExpression>),
    Point(Point),
    ConvexHull(Vec<Point>),
    ConvexHullInterior(Vec<Point>),
    ConvexHullBoundary(Vec<Point>),
    Lines(Vec<Point>),
    Loop(Vec<Point>),
    Polygon(Vec<Point>),
    PolygonInterior(Vec<Point>),
}

impl ShapeExpression {
    fn dimension(&self) -> Result<usize, ()> {
        let mut same_dim = None;
        for next_dim in match self {
            ShapeExpression::Union(left, right) => {
                vec![left.dimension(), right.dimension()]
            }
            ShapeExpression::Intersect(left, right) => {
                vec![left.dimension(), right.dimension()]
            }
            ShapeExpression::Difference(left, right) => {
                vec![left.dimension(), right.dimension()]
            }
            ShapeExpression::Point(point) => {
                vec![Ok(point.dimension())]
            }
            ShapeExpression::ConvexHull(points) => {
                points.iter().map(|pt| Ok(pt.dimension())).collect()
            }
            ShapeExpression::ConvexHullInterior(points) => {
                points.iter().map(|pt| Ok(pt.dimension())).collect()
            }
            ShapeExpression::ConvexHullBoundary(points) => {
                points.iter().map(|pt| Ok(pt.dimension())).collect()
            }
            ShapeExpression::Lines(points) => points.iter().map(|pt| Ok(pt.dimension())).collect(),
            ShapeExpression::Loop(points) => points.iter().map(|pt| Ok(pt.dimension())).collect(),
            ShapeExpression::Polygon(points) => {
                points.iter().map(|pt| Ok(pt.dimension())).collect()
            }
            ShapeExpression::PolygonInterior(points) => {
                points.iter().map(|pt| Ok(pt.dimension())).collect()
            }
        } {
            if let Ok(next_dim) = next_dim {
                if let Some(current_d) = same_dim {
                    if current_d != next_dim {
                        return Err(());
                    }
                } else {
                    same_dim = Some(next_dim)
                }
            } else {
                return Err(());
            }
        }
        if let Some(dim) = same_dim {
            Ok(dim)
        } else {
            Err(())
        }
    }

    fn to_partial_simplicial_complex(
        &self,
        space: AffineSpace<'static, RationalCanonicalStructure>,
    ) -> PartialSimplicialComplex<'static, RationalCanonicalStructure> {
        match self {
            ShapeExpression::Union(left, right) => left
                .to_partial_simplicial_complex(space)
                .union(&right.to_partial_simplicial_complex(space)),
            ShapeExpression::Intersect(left, right) => left
                .to_partial_simplicial_complex(space)
                .intersect(&right.to_partial_simplicial_complex(space)),
            ShapeExpression::Difference(left, right) => left
                .to_partial_simplicial_complex(space)
                .difference(&right.to_partial_simplicial_complex(space)),
            ShapeExpression::Point(point) => space
                .convex_hull(vec![point.to_vector(space)])
                .to_simplicial_complex()
                .into_forget_labels()
                .into_partial_simplicial_complex(),
            ShapeExpression::ConvexHull(points) => space
                .convex_hull(points.iter().map(|pt| pt.to_vector(space)).collect())
                .to_simplicial_complex()
                .into_forget_labels()
                .into_partial_simplicial_complex(),
            ShapeExpression::ConvexHullInterior(points) => space
                .convex_hull(points.iter().map(|pt| pt.to_vector(space)).collect())
                .to_simplicial_complex()
                .interior()
                .into_partial_simplicial_complex(),
            ShapeExpression::ConvexHullBoundary(points) => space
                .convex_hull(points.iter().map(|pt| pt.to_vector(space)).collect())
                .to_simplicial_complex()
                .boundary()
                .into_partial_simplicial_complex(),
            ShapeExpression::Lines(points) => {
                let mut shape = space
                    .convex_hull(vec![])
                    .to_simplicial_complex()
                    .into_forget_labels()
                    .into_partial_simplicial_complex();
                for i in 0..points.len() - 1 {
                    let p0 = points[i].to_vector(space);
                    let p1 = points[i + 1].to_vector(space);
                    shape = shape.union(
                        &space
                            .convex_hull(vec![p0, p1])
                            .to_simplicial_complex()
                            .into_forget_labels(),
                    );
                }
                shape
            }
            ShapeExpression::Loop(points) => {
                let mut shape = space
                    .convex_hull(vec![])
                    .to_simplicial_complex()
                    .into_forget_labels()
                    .into_partial_simplicial_complex();
                let n = points.len();
                for i in 0..n {
                    let p0 = points[i].to_vector(space);
                    let p1 = points[(i + 1) % n].to_vector(space);
                    shape = shape.union(
                        &space
                            .convex_hull(vec![p0, p1])
                            .to_simplicial_complex()
                            .into_forget_labels(),
                    );
                }
                shape
            }
            ShapeExpression::Polygon(points) => ShapeExpression::PolygonInterior(points.clone())
                .to_partial_simplicial_complex(space)
                .union(&ShapeExpression::Loop(points.clone()).to_partial_simplicial_complex(space)),
            ShapeExpression::PolygonInterior(points) => {
                assert_eq!(space.linear_dimension(), Some(2));

                let mut shape = space
                    .convex_hull(vec![])
                    .to_simplicial_complex()
                    .into_forget_labels()
                    .into_partial_simplicial_complex();

                if points.is_empty() {
                    shape
                } else {
                    let root = space.origin().unwrap();
                    let n = points.len();

                    for i in 0..n {
                        let p0 = points[i].to_vector(space);
                        let p1 = points[(i + 1) % n].to_vector(space);
                        let part = space
                            .convex_hull(vec![root.clone(), p0, p1])
                            .to_simplicial_complex()
                            .interior();
                        shape = shape.difference(&part).union(&part.difference(&shape));
                    }

                    PartialSimplicialComplex::new_unchecked(
                        space,
                        shape
                            .simplexes()
                            .into_iter()
                            .filter(|spx| spx.n() == 3)
                            .cloned()
                            .collect(),
                    )
                    .closure()
                    .interior()
                }
            }
        }
    }

    pub fn to_partial_simplicial_complex_root(
        &self,
    ) -> PartialSimplicialComplex<'static, RationalCanonicalStructure> {
        let n = self.dimension().unwrap();
        let space = AffineSpace::new_linear(Rational::structure_ref(), n);
        self.to_partial_simplicial_complex(space)
    }
}
