use std::collections::HashSet;

use malachite_q::Rational;

use crate::geometry::simplex::Simplex;

use super::{
    affine_coordinate_system::AffineSubspaceCoordinateSystem, oriented_simplex::OrientedSimplex,
    shape::Shape, simplicial_complex::SimplicialComplex, vector::Vector,
};

#[derive(Debug, Clone, Copy)]
pub enum ConvexSimplicialComplexPointResult {
    Inside,            //point lies inside the convex hull
    Boundary,          //point lies on the boundary of the convex hull
    Outside,           //point lines in the plane of the convex hull but is outside
    OutsideHyperplane, //point lies outside the plane of the convex hull
}

#[derive(Debug, Clone)]
enum ConvexSimplicialComplexCases {
    NonEmpty {
        subspace: AffineSubspaceCoordinateSystem, //an embedding of a vector space of dimension self.rank into the ambient space
        boundary: Vec<OrientedSimplex>, //oriented simplexes in the subspace of rank self.rank-1 with positive side on the outside
        interior: Vec<Simplex>,         //simplexes in the subspace which form the interior
                                        //furthermore, the convex simplicial complex should contain the embedding of the standard simplex in the AffineSubspace
    },
    Empty {
        dim: usize,
    },
}

#[derive(Debug, Clone)]
pub struct ConvexSimplicialComplex(ConvexSimplicialComplexCases);

impl ConvexSimplicialComplex {
    pub fn check(&self) -> Result<(), &'static str> {
        match self {
            ConvexSimplicialComplex(ConvexSimplicialComplexCases::NonEmpty {
                subspace,
                boundary,
                interior,
            }) => {
                SimplicialComplex::new(subspace.rank(), interior.clone()).check()?;

                // if interior.dim() != subspace.rank() {
                //     return Err("interior simplexes should live in the subspace coordinates");
                // }

                for simplex in boundary {
                    if simplex.dim() != subspace.rank() {
                        return Err("boundary simplexes should live in the subspace coordinates");
                    }
                }

                for simplex in interior {
                    if simplex.dim() != subspace.rank() {
                        return Err("interior simplexes should live in the subspace coordinates");
                    }

                    for b in boundary {
                        if !b.sign_point(&simplex.centroid()).is_le() {
                            return Err("centroid of interior simplexes should be on the negative side of every boundary simplex");
                        }

                        //not valid when subspace.rank() == 0
                        if subspace.rank() >= 1 {
                            if !b
                                .sign_point(&self.subspace_interior_point().unwrap())
                                .is_lt()
                            {
                                return Err("centroid of the standard simplex should be strictly inside each boundary when rank >= 1");
                            }
                        }
                    }
                }
            }
            ConvexSimplicialComplex(ConvexSimplicialComplexCases::Empty { dim }) => {}
        }

        Ok(())
    }

    pub fn as_simplicial_complex(&self) -> SimplicialComplex {
        match self {
            ConvexSimplicialComplex(ConvexSimplicialComplexCases::NonEmpty {
                subspace,
                boundary,
                interior,
            }) => subspace.simplicial_complex_image(&SimplicialComplex::new(
                subspace.rank(),
                interior.clone(),
            )),
            ConvexSimplicialComplex(ConvexSimplicialComplexCases::Empty { dim }) => {
                SimplicialComplex::empty(*dim)
            }
        }
    }

    #[deprecated(note = "for debug use only")]
    pub fn interior_as_shape(&self) -> Shape {
        self.as_simplicial_complex().as_shape()
    }

    #[deprecated(note = "for debug use only")]
    pub fn boundary_as_shape(&self) -> Shape {
        match self {
            ConvexSimplicialComplex(ConvexSimplicialComplexCases::NonEmpty {
                subspace,
                boundary,
                interior,
            }) => subspace.shape_image(&Shape::new(
                subspace.rank(),
                boundary.iter().map(|s| s.clone().as_simplex()).collect(),
            )),
            ConvexSimplicialComplex(ConvexSimplicialComplexCases::Empty { dim }) => {
                Shape::empty(*dim)
            }
        }
    }

    fn subspace_interior_point(&self) -> Option<Vector> {
        match self {
            ConvexSimplicialComplex(ConvexSimplicialComplexCases::NonEmpty {
                subspace,
                boundary,
                interior,
            }) => Some(Simplex::new_standard_simplex(subspace.rank()).centroid()),
            ConvexSimplicialComplex(ConvexSimplicialComplexCases::Empty { dim }) => None,
        }
    }

    fn interior_point(&self) -> Option<Vector> {
        match self {
            ConvexSimplicialComplex(ConvexSimplicialComplexCases::NonEmpty {
                subspace,
                boundary,
                interior,
            }) => Some(subspace.point_image(&self.subspace_interior_point().unwrap())),
            ConvexSimplicialComplex(ConvexSimplicialComplexCases::Empty { dim }) => None,
        }
    }

    pub fn new_empty(dim: usize) -> Self {
        Self(ConvexSimplicialComplexCases::Empty { dim })
    }

    pub fn new_from_simplex(simplex: &Simplex) -> Self {
        assert!(simplex.rank().is_some());
        let subspace = simplex.affine_subspace();
        let simplex_preimage = Simplex::new_standard_simplex(simplex.rank().unwrap());
        let boundary = simplex_preimage.oriented_facets();
        let interior = simplex_preimage.as_simplicial_complex().simplices();

        let ans = Self(ConvexSimplicialComplexCases::NonEmpty {
            subspace,
            boundary,
            interior,
        });

        debug_assert!(ans.check().is_ok());
        ans
    }

    pub fn dim(&self) -> usize {
        match self {
            ConvexSimplicialComplex(ConvexSimplicialComplexCases::NonEmpty {
                subspace,
                boundary,
                interior,
            }) => subspace.dim(),
            ConvexSimplicialComplex(ConvexSimplicialComplexCases::Empty { dim }) => *dim,
        }
    }

    pub fn rank(&self) -> Option<usize> {
        match self {
            ConvexSimplicialComplex(ConvexSimplicialComplexCases::NonEmpty {
                subspace,
                boundary,
                interior,
            }) => Some(subspace.rank()),
            ConvexSimplicialComplex(ConvexSimplicialComplexCases::Empty { dim }) => None,
        }
    }

    pub fn contains_point(&self, point: &Vector) -> ConvexSimplicialComplexPointResult {
        match self {
            ConvexSimplicialComplex(ConvexSimplicialComplexCases::NonEmpty {
                subspace,
                boundary,
                interior,
            }) => match subspace.point_preimage(point) {
                Some(point_preimage) => {
                    let mut sign = std::cmp::Ordering::Less;
                    for b in boundary {
                        sign = std::cmp::Ordering::max(sign, b.sign_point(&point_preimage));
                        if sign.is_gt() {
                            //an optimization - we can skip checking the other faces if the point is found to be outside
                            return ConvexSimplicialComplexPointResult::Outside;
                        }
                    }
                    match sign {
                        std::cmp::Ordering::Less => ConvexSimplicialComplexPointResult::Inside,
                        std::cmp::Ordering::Equal => ConvexSimplicialComplexPointResult::Boundary,
                        std::cmp::Ordering::Greater => ConvexSimplicialComplexPointResult::Outside,
                    }
                }
                None => ConvexSimplicialComplexPointResult::OutsideHyperplane,
            },
            ConvexSimplicialComplex(ConvexSimplicialComplexCases::Empty { dim }) => {
                ConvexSimplicialComplexPointResult::OutsideHyperplane
            }
        }
    }

    ///Return the convex hull of self and point by extending self.
    pub fn extend_by_point(&self, point: &Vector) -> Self {
        match self {
            ConvexSimplicialComplex(ConvexSimplicialComplexCases::NonEmpty {
                subspace,
                boundary,
                interior,
            }) => {
                match subspace.point_preimage(point) {
                    Some(point_preimage) => {
                        //find which boundary simplices are visible from the point
                        let mut visible_boundary = vec![];
                        let mut obstructed_boundary = vec![];
                        for b in boundary {
                            match b.sign_point(&point_preimage).is_gt() {
                                true => visible_boundary.push(b.clone()),
                                false => obstructed_boundary.push(b.clone()),
                            }
                        }

                        //find the facets of the visible_boundary which are at the horizon of what is visible
                        let mut horizon = HashSet::new();
                        for v in &visible_boundary {
                            for b in v.clone().as_simplex().facets() {
                                if horizon.contains(&b) {
                                    horizon.remove(&b);
                                } else {
                                    horizon.insert(b);
                                }
                            }
                        }

                        let mut extended_boundary = vec![];
                        let mut extended_interior = vec![];

                        //compute the extended boundary
                        extended_boundary.append(&mut obstructed_boundary);
                        for h in horizon {
                            extended_boundary.push(OrientedSimplex::from_simplex(
                                h.extend_by_point(point_preimage.clone()),
                                &self.subspace_interior_point().unwrap(),
                            ));
                        }

                        //the current interior
                        extended_interior.append(&mut interior.clone());

                        //cone of the visible boundary
                        let mut visible_boundary_completion = HashSet::new();
                        for b in visible_boundary.iter().map(|b| b.clone().as_simplex()) {
                            for bb in b.boundary_simplices() {
                                visible_boundary_completion.insert(bb);
                            }
                            visible_boundary_completion.insert(b);
                        }

                        for b in visible_boundary_completion {
                            extended_interior.push(b.extend_by_point(point_preimage.clone()));
                        }

                        //add the point to the new interior if it is not already inside
                        if !visible_boundary.is_empty() {
                            extended_interior
                                .push(Simplex::new(subspace.rank(), vec![point_preimage]));
                        }

                        let ans = Self(ConvexSimplicialComplexCases::NonEmpty {
                            subspace: subspace.clone(),
                            boundary: extended_boundary,
                            interior: extended_interior,
                        });

                        debug_assert!(ans.check().is_ok());

                        ans
                    }
                    None => {
                        // subspace C extended_subspace C ambient_space

                        let extended_interior_point =
                            Simplex::new_standard_simplex(subspace.rank() + 1).centroid();

                        let extended_subspace_into_ambient =
                            subspace.clone().extend_basis(point - subspace.origin());

                        let point_extended_preimage = Vector::new(
                            (0..subspace.rank() + 1)
                                .map(|i| match i == subspace.rank() {
                                    true => Rational::from(1),
                                    false => Rational::from(0),
                                })
                                .collect(),
                        );

                        let subspace_into_extended_subspace =
                            AffineSubspaceCoordinateSystem::cannonical_subspace(
                                subspace.rank(),
                                subspace.rank() + 1,
                            );

                        let mut extended_boundary = vec![];
                        let mut extended_interior = vec![];

                        for boundary_simplex in boundary {
                            // add its cone to the new boundary
                            extended_boundary.push(OrientedSimplex::from_simplex(
                                subspace_into_extended_subspace
                                    .simplex_image(&boundary_simplex.clone().as_simplex())
                                    .extend_by_point(point_extended_preimage.clone()),
                                &extended_interior_point,
                            ));
                        }

                        for interior_simplex in interior {
                            // add it to the new boundary if it is of maximal rank
                            if interior_simplex.rank().unwrap() == subspace.rank() {
                                extended_boundary.push(OrientedSimplex::from_simplex(
                                    subspace_into_extended_subspace.simplex_image(interior_simplex),
                                    &extended_interior_point,
                                ));
                            }

                            // add itself to the new interior
                            extended_interior.push(
                                subspace_into_extended_subspace.simplex_image(interior_simplex),
                            );

                            // add its cone to the new interior
                            extended_interior.push(
                                subspace_into_extended_subspace
                                    .simplex_image(interior_simplex)
                                    .extend_by_point(point_extended_preimage.clone()),
                            );
                        }

                        //add the point to the new interior
                        extended_interior.push(Simplex::new(
                            subspace.rank() + 1,
                            vec![point_extended_preimage],
                        ));

                        let ans = Self(ConvexSimplicialComplexCases::NonEmpty {
                            subspace: extended_subspace_into_ambient,
                            boundary: extended_boundary,
                            interior: extended_interior,
                        });

                        debug_assert!(ans.check().is_ok());
                        ans
                    }
                }
            }
            ConvexSimplicialComplex(ConvexSimplicialComplexCases::Empty { dim }) => {
                Self::new_from_simplex(&Simplex::new(*dim, vec![point.clone()]))
            }
        }
    }
}

pub fn convex_hull(dim: usize, points: Vec<Vector>) -> ConvexSimplicialComplex {
    for point in &points {
        assert_eq!(dim, point.dim());
    }
    let mut ch = ConvexSimplicialComplex::new_empty(dim);
    for point in points {
        // println!("{:?}", point);
        ch = ch.extend_by_point(&point);
    }
    ch
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use super::*;
    #[test]
    fn test_convex_simplicial_complex_point_extension() {
        let csc0 = ConvexSimplicialComplex::new_empty(2);
        assert!(csc0.check().is_ok());
        assert_eq!(csc0.rank(), None);

        let csc1 = csc0.extend_by_point(&Vector::new(vec![
            Rational::from_str("1").unwrap(),
            Rational::from_str("1").unwrap(),
        ]));
        assert!(csc1.check().is_ok());
        assert_eq!(csc1.rank(), Some(0));

        let csc2 = csc1.extend_by_point(&Vector::new(vec![
            Rational::from_str("1").unwrap(),
            Rational::from_str("1").unwrap(),
        ]));
        assert!(csc2.check().is_ok());
        assert_eq!(csc2.rank(), Some(0));

        let csc3 = csc2.extend_by_point(&Vector::new(vec![
            Rational::from_str("2").unwrap(),
            Rational::from_str("3").unwrap(),
        ]));
        assert!(csc3.check().is_ok());
        assert_eq!(csc3.rank(), Some(1));

        let csc4 = csc3.extend_by_point(&Vector::new(vec![
            Rational::from_str("4").unwrap(),
            Rational::from_str("7").unwrap(),
        ]));
        assert!(csc4.check().is_ok());
        assert_eq!(csc4.rank(), Some(1));

        let csc5 = csc4.extend_by_point(&Vector::new(vec![
            Rational::from_str("3").unwrap(),
            Rational::from_str("5").unwrap(),
        ]));
        assert!(csc5.check().is_ok());
        assert_eq!(csc5.rank(), Some(1));

        let csc6 = csc5.extend_by_point(&Vector::new(vec![
            Rational::from_str("0").unwrap(),
            Rational::from_str("6").unwrap(),
        ]));
        assert!(csc6.check().is_ok());
        assert_eq!(csc6.rank(), Some(2));

        let csc7 = csc6.extend_by_point(&Vector::new(vec![
            Rational::from_str("1").unwrap(),
            Rational::from_str("0").unwrap(),
        ]));
        assert!(csc7.check().is_ok());
        assert_eq!(csc7.rank(), Some(2));

        let csc8 = csc7.extend_by_point(&Vector::new(vec![
            Rational::from_str("4").unwrap(),
            Rational::from_str("1").unwrap(),
        ]));
        assert!(csc8.check().is_ok());
        assert_eq!(csc8.rank(), Some(2));
    }
}
