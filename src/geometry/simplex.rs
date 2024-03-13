use itertools::Itertools;
use malachite_q::Rational;

use crate::rings::lattice::{AffineLattice, LinearLattice};

use super::{
    affine_coordinate_system::AffineSubspaceCoordinateSystem,
    box_region::BoxRegion,
    oriented_simplex::OrientedSimplex,
    simplicial_complex::SimplicialComplex,
    vector::{are_points_nondegenerage, Vector},
};

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Simplex {
    dim: usize,
    vertices: Vec<Vector>, //ordered
}

impl Simplex {
    pub fn check(&self) -> Result<(), &'static str> {
        let mut sorted_vertices = self.vertices.clone();
        sorted_vertices.sort();
        if self.vertices != sorted_vertices {
            return Err("Simplex vertices are not sorted");
        }

        // if self.vertices.is_empty() {
        //     return Err("Simplex should have at least one vertex");
        // }

        for p in &self.vertices {
            if p.dim() != self.dim {
                return Err("Simplex should live in the same dimension as its vertices");
            }
        }

        if !are_points_nondegenerage(self.dim, self.vertices.iter().collect()) {
            return Err("Simplex is degenerate");
        }

        Ok(())
    }

    pub fn new(dim: usize, mut vertices: Vec<Vector>) -> Self {
        vertices.sort();
        let ans = Self { dim, vertices };
        ans.check().unwrap();
        ans
    }

    pub fn new_unchecked(dim: usize, mut vertices: Vec<Vector>) -> Self {
        vertices.sort();
        let ans = Self { dim, vertices };
        debug_assert!(ans.check().is_ok());
        ans
    }

    pub fn new_standard_simplex(rank: usize) -> Self {
        //the simplex whose verticies are 0 and the standard basis vectors
        let mut vertices = vec![];
        vertices.push(Vector::new((0..rank).map(|_j| Rational::from(0)).collect()));
        for i in 0..rank {
            vertices.push(Vector::new(
                (0..rank)
                    .map(|j| match i == j {
                        true => Rational::from(1),
                        false => Rational::from(0),
                    })
                    .collect(),
            ));
        }

        Simplex::new(rank, vertices)
    }

    pub fn try_new(dim: usize, mut vertices: Vec<Vector>) -> Option<Self> {
        if are_points_nondegenerage(dim, vertices.iter().collect()) {
            vertices.sort();
            Some(Self { dim, vertices })
        } else {
            None
        }
    }

    pub fn n(&self) -> usize {
        self.vertices.len()
    }

    pub fn rank(&self) -> Option<usize> {
        if self.vertices.is_empty() {
            None
        } else {
            Some(self.vertices.len() - 1)
        }
    }

    pub fn dim(&self) -> usize {
        self.dim
    }

    pub fn points(&self) -> Vec<Vector> {
        self.vertices.clone()
    }

    pub fn points_ref(&self) -> Vec<&Vector> {
        self.vertices.iter().collect()
    }

    pub fn bounding_box(&self) -> BoxRegion {
        BoxRegion::bounding_box(self.dim, &self.vertices)
    }

    #[deprecated]
    pub fn transform(self, new_dim: usize, f: &dyn Fn(Vector) -> Vector) -> Self {
        let new_vertices = self
            .vertices
            .into_iter()
            .map(|p| p.transform(f))
            .collect_vec();
        for pt in new_vertices.iter() {
            debug_assert_eq!(pt.dim(), new_dim);
        }
        Self::new(new_dim, new_vertices)
    }

    pub fn has_vertex(&self, pt: &Vector) -> bool {
        self.vertices.binary_search(pt).is_ok()
    }

    pub fn skeleton(&self, skel_n: usize) -> Vec<Simplex> {
        let mut parts = vec![];
        for skeleton_piece_index in (0..self.vertices.len()).combinations(skel_n) {
            let part = Simplex {
                dim: self.dim,
                vertices: skeleton_piece_index
                    .into_iter()
                    .map(|i| self.vertices[i].clone())
                    .collect(),
            };
            parts.push(part);
        }
        return parts;
    }

    pub fn vertices(&self) -> Vec<Simplex> {
        self.skeleton(0)
    }

    pub fn edges(&self) -> Vec<Simplex> {
        self.skeleton(1)
    }

    pub fn faces(&self) -> Vec<Simplex> {
        self.skeleton(2)
    }

    pub fn ridges(&self) -> Vec<Simplex> {
        self.skeleton(self.n() - 2)
    }

    pub fn facets(&self) -> Vec<Simplex> {
        self.skeleton(self.n() - 1)
    }

    pub fn facet(&self, k: usize) -> Simplex {
        assert!(k <= self.n());
        let facet = Simplex {
            dim: self.dim,
            vertices: (0..self.n())
                .filter(|i| i != &k)
                .map(|i| self.vertices[i].clone())
                .collect(),
        };
        debug_assert!(facet.check().is_ok());
        facet
    }

    pub fn oriented_facet(&self, k: usize) -> OrientedSimplex {
        //return the oriented facet of self with positive side on the outside and negative side on the inside
        assert_eq!(self.dim, self.rank().unwrap());
        assert!(k <= self.n());
        let oriented_facet = OrientedSimplex::new(self.facet(k));
        if self.dim() == 0 {
            //self.dim == self.rank == k == 0
            //in this case, we are looking for the oriented facet of a point which is the empty simplex
            //both orientations of the empty simplex are the same and no points are on either side. the unique point has signed distance 0 from it but does not lie on it
            oriented_facet
        } else {
            match oriented_facet.sign_point(&self.vertices[k]) {
                std::cmp::Ordering::Less => oriented_facet,
                std::cmp::Ordering::Equal => panic!(),
                std::cmp::Ordering::Greater => oriented_facet.flipped(),
            }
        }
    }

    pub fn oriented_facets(&self) -> Vec<OrientedSimplex> {
        assert_eq!(self.dim, self.rank().unwrap());
        (0..self.n()).map(|k| self.oriented_facet(k)).collect()
    }

    pub fn boundary_simplices(&self) -> Vec<Simplex> {
        let n = self.n();
        let total: u128 = 1 << n;
        let mut parts = vec![];
        for subset_mask in 1..total - 1 {
            let mut subset = Vec::new();
            for bit_pos in 0..n {
                if subset_mask & (1 << bit_pos) != 0 {
                    subset.push(bit_pos);
                }
            }
            let b = Simplex {
                dim: self.dim,
                vertices: subset
                    .into_iter()
                    .map(|i| self.vertices[i].clone())
                    .collect(),
            };
            debug_assert!(b.check().is_ok());
            parts.push(b);
        }
        parts
    }

    pub fn boundary(&self) -> SimplicialComplex {
        let ans = SimplicialComplex::new(self.dim, self.boundary_simplices());
        debug_assert!(ans.check().is_ok());
        ans
    }

    pub fn as_simplicial_complex(&self) -> SimplicialComplex {
        let mut parts = self.boundary_simplices();
        parts.push(self.clone());
        let ans = SimplicialComplex::new(self.dim, parts);
        debug_assert!(ans.check().is_ok());
        ans
    }

    #[deprecated(note = "use affine_subspace instead")]
    pub fn affine_span(&self) -> AffineLattice<Rational> {
        if self.vertices.len() == 0 {
            AffineLattice::empty(self.dim, 1)
        } else {
            AffineLattice::from_offset_and_linear_lattice(
                self.dim,
                1,
                self.vertices[0].as_matrix(),
                LinearLattice::from_basis(
                    self.dim,
                    1,
                    (1..self.vertices.len())
                        .map(|i| (&self.vertices[i] - &self.vertices[0]).as_matrix())
                        .collect(),
                ),
            )
        }
    }

    pub fn centroid(&self) -> Vector {
        let mut coords = (0..self.dim).map(|_i| Rational::from(0)).collect_vec();
        for pt in &self.vertices {
            for i in 0..self.dim {
                coords[i] += pt.get_coord(i);
            }
        }
        coords = coords
            .into_iter()
            .map(|c| c / Rational::from(self.n()))
            .collect();

        Vector::new(coords)
    }

    pub fn extend_by_point(mut self, point: Vector) -> Self {
        self.vertices.push(point);
        self.vertices.sort();
        debug_assert!(self.check().is_ok());
        self
    }

    // pub fn orthogonal_project(&self, point: &Point) -> Point {
    //     assert_ne!(self.n(), 0);
    //     let root = &self.vertices[0];
    //     let vecs = (1..self.n())
    //         .map(|i| &self.vertices[i] - root)
    //         .collect_vec();
    //     let pt_vec = point - root;

    //     let mut proj = root.clone();
    //     for vec in vecs {
    //         proj = &proj + &(&vec * &(pt_vec.dot(&vec) / vec.length_sq()));
    //     }

    //     proj
    // }

    // pub fn orthogonal_distance_sq(&self, point: &Point) -> Rational {
    //     let proj = self.orthogonal_project(point);
    //     (point - &proj).length_sq()
    // }
}
