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
    ConvexHull(Vec<Point>),
    ConvexHullInterior(Vec<Point>),
    ConvexHullBoundary(Vec<Point>),
    Lines(Vec<Point>),
    Loop(Vec<Point>),
    Polygon(Vec<Point>),
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
            ShapeExpression::Lines(_points) => todo!(),
            ShapeExpression::Loop(_points) => todo!(),
            ShapeExpression::Polygon(_points) => todo!(),
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
