use crate::{
    geometry::vector::Vector,
    rings::linear::{
        lattice::{AffineLattice, LinearLattice},
        matrix::Matrix,
    },
};

use super::{shape::Shape, simplex::Simplex, simplicial_complex::SimplicialComplex};
use malachite_q::Rational;

#[derive(Debug, Clone)]
pub struct AffineSubspaceCoordinateSystem {
    dim: usize,         //the dimension of the ambient space
    origin: Vector,     //the origin of the affine subspace in the ambient space
    basis: Vec<Vector>, //a basis for the affine subspace relative to the origin
}

impl AffineSubspaceCoordinateSystem {
    pub fn check(&self) -> Result<(), &'static str> {
        if self.dim != self.origin.dim() {
            return Err("origin of affine subspace should live in the ambient space");
        }
        for v in &self.basis {
            if self.dim != self.origin.dim() {
                return Err("basis vector of affine subspace should live in the ambient space");
            }
        }

        if self.basis.len() != self.basis_matrix().rank() {
            return Err("affine subspace vectors should be linearly independent");
        }

        Ok(())
    }

    pub fn origin(&self) -> &Vector {
        &self.origin
    }

    //a coordinate system for the affine span of the given points
    pub fn span_of_points(dim: usize, mut points: Vec<&Vector>) -> Option<Self> {
        //put the points into the columns of a matrix
        let mat = Matrix::join_cols(
            dim,
            points.into_iter().map(|point| point.as_matrix()).collect(),
        );
        let afflat = mat.col_affine_span();
        Self::from_affine_lattice(afflat)
    }

    fn basis_matrix(&self) -> Matrix<Rational> {
        Matrix::construct(self.dim, self.basis.len(), |i, j| {
            self.basis[j].get_coord(i)
        })
    }

    pub fn cannonical_subspace(k: usize, n: usize) -> Self {
        // k <= n
        // The cannonical embedding of A^k into A^n by taking the first k out of n coordinates
        assert!(k <= n);
        Self {
            dim: n,
            origin: Vector::new((0..n).map(|_i| Rational::from(0)).collect()),
            basis: (0..k)
                .map(|i| {
                    Vector::new(
                        (0..n)
                            .map(|j| match i == j {
                                true => Rational::from(1),
                                false => Rational::from(0),
                            })
                            .collect(),
                    )
                })
                .collect(),
        }
    }

    pub fn from_affine_lattice(
        afflat: AffineLattice<Rational>,
    ) -> Option<Self> {
        let dim = afflat.rows();
        match afflat.to_offset_and_linear_lattice() {
            Some((offset, linlat)) => Some(Self {
                dim,
                origin: Vector::from_matrix(&offset),
                basis: linlat
                    .basis_matrices()
                    .into_iter()
                    .map(|v| Vector::from_matrix(&v))
                    .collect(),
            }),
            None => None,
        }
    }

    pub fn to_affine_lattice(&self) -> AffineLattice<Rational> {
        AffineLattice::from_offset_and_linear_lattice(
            self.dim,
            1,
            self.origin.as_matrix(),
            LinearLattice::from_basis(
                self.dim,
                1,
                (0..self.basis.len())
                    .map(|i| (&self.basis[i]).as_matrix())
                    .collect(),
            ),
        )
    }

    pub fn dim(&self) -> usize {
        self.dim
    }

    pub fn rank(&self) -> usize {
        self.basis.len()
    }

    pub fn point_image(&self, p: &Vector) -> Vector {
        assert_eq!(p.dim(), self.rank());
        &self.origin
            + &Vector::from_matrix(&Matrix::mul(&self.basis_matrix(), &p.as_matrix()).unwrap())
    }

    pub fn simplex_image(&self, s: &Simplex) -> Simplex {
        Simplex::new(
            self.dim(),
            s.points().iter().map(|p| self.point_image(p)).collect(),
        )
    }

    pub fn shape_image(&self, s: &Shape) -> Shape {
        Shape::new(
            self.dim(),
            s.simplices_ref()
                .into_iter()
                .map(|s| self.simplex_image(s))
                .collect(),
        )
    }

    pub fn simplicial_complex_image(&self, s: &SimplicialComplex) -> SimplicialComplex {
        SimplicialComplex::new(
            self.dim(),
            s.simplices_ref()
                .into_iter()
                .map(|s| self.simplex_image(s))
                .collect(),
        )
    }

    pub fn point_preimage(&self, p: &Vector) -> Option<Vector> {
        assert_eq!(p.dim(), self.dim());

        match self
            .basis_matrix()
            .col_solve(&(p - &self.origin).as_matrix())
        {
            Some(sol) => Some(Vector::from_matrix(&sol)),
            None => None,
        }
    }

    pub fn simplex_preimage(&self, s: &Simplex) -> Option<Simplex> {
        let mut new_points = vec![];
        for p in &s.points() {
            match self.point_preimage(p) {
                Some(q) => new_points.push(q),
                None => {
                    return None;
                }
            }
        }
        Some(Simplex::new(self.rank(), new_points))
    }

    pub fn shape_preimage(&self, s: &Shape) -> Option<Shape> {
        let mut new_simplices = vec![];
        for a in s.simplices_ref() {
            match self.simplex_preimage(a) {
                Some(b) => new_simplices.push(b),
                None => {
                    return None;
                }
            }
        }
        Some(Shape::new(self.rank(), new_simplices))
    }

    pub fn simplicial_complex_preimage(&self, s: &SimplicialComplex) -> Option<SimplicialComplex> {
        let mut new_simplices = vec![];
        for a in s.simplices_ref() {
            match self.simplex_preimage(a) {
                Some(b) => new_simplices.push(b),
                None => {
                    return None;
                }
            }
        }
        Some(SimplicialComplex::new(self.rank(), new_simplices))
    }

    pub fn extend_basis(mut self, new_basis_vector: Vector) -> Self {
        self.basis.push(new_basis_vector);
        debug_assert!(self.check().is_ok());
        self
    }
}

impl Simplex {
    pub fn affine_subspace(&self) -> AffineSubspaceCoordinateSystem {
        //the affine subspace coordinate system in which self is the standard simplex
        AffineSubspaceCoordinateSystem {
            dim: self.dim(),
            origin: self.points().iter().nth(0).unwrap().clone(),
            basis: (1..self.points().len())
                .map(|i| {
                    self.points().iter().nth(i).unwrap() - self.points().iter().nth(0).unwrap()
                })
                .collect(),
        }
    }
}

pub fn affine_span_of_simplices(
    dim: usize,
    simplices: Vec<&Simplex>,
) -> Option<AffineSubspaceCoordinateSystem> {
    let mut points = vec![];
    for simplex in &simplices {
        for point in simplex.points() {
            points.push(point);
        }
    }
    AffineSubspaceCoordinateSystem::span_of_points(dim, points.iter().collect())
}
