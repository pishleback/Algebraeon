use core::panic;
use std::borrow::BorrowMut;
#[allow(dead_code)]
use std::collections::{HashMap, HashSet};

use itertools::Itertools;
use malachite_nz::natural::arithmetic::mul::poly_eval;
use malachite_q::arithmetic::simplest_rational_in_interval;
use malachite_q::Rational;
use rayon::iter::ParallelIterator;

use crate::rings::nzq::QQ;
use crate::rings::ring::Real;

use super::rings::lattice::*;
use super::rings::matrix::*;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Point {
    coords: Vec<Rational>,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Vector {
    coords: Vec<Rational>,
}

impl std::ops::Neg for Vector {
    type Output = Vector;

    fn neg(self) -> Self::Output {
        Vector {
            coords: self.coords.into_iter().map(|x| -x).collect(),
        }
    }
}
impl std::ops::Mul<&Vector> for &Rational {
    type Output = Vector;

    fn mul(self, other: &Vector) -> Self::Output {
        Vector {
            coords: other.coords.iter().map(|x| self * x).collect(),
        }
    }
}
impl std::ops::Mul<&Rational> for &Vector {
    type Output = Vector;

    fn mul(self, other: &Rational) -> Self::Output {
        Vector {
            coords: self.coords.iter().map(|x| x * other).collect(),
        }
    }
}
impl std::ops::Neg for &Vector {
    type Output = Vector;

    fn neg(self) -> Self::Output {
        -self.clone()
    }
}
impl std::ops::Add<&Vector> for &Vector {
    type Output = Vector;
    fn add(self, other: &Vector) -> Self::Output {
        let n = self.coords.len();
        if n == other.coords.len() {
            Vector {
                coords: (0..n).map(|i| &self.coords[i] + &other.coords[i]).collect(),
            }
        } else {
            panic!("Can only add vectors in the same dimension")
        }
    }
}
impl std::ops::Add<&Vector> for &Point {
    type Output = Point;
    fn add(self, other: &Vector) -> Self::Output {
        let n = self.coords.len();
        if n == other.coords.len() {
            Point {
                coords: (0..n).map(|i| &self.coords[i] + &other.coords[i]).collect(),
            }
        } else {
            panic!("Can only add a vector to a point in the same dimension")
        }
    }
}
impl std::ops::Sub<&Vector> for &Point {
    type Output = Point;
    fn sub(self, other: &Vector) -> Self::Output {
        let n = self.coords.len();
        if n == other.coords.len() {
            Point {
                coords: (0..n).map(|i| &self.coords[i] - &other.coords[i]).collect(),
            }
        } else {
            panic!("Can only subtract a vector from a point in the same dimension")
        }
    }
}
impl std::ops::Sub<&Vector> for &Vector {
    type Output = Vector;
    fn sub(self, other: &Vector) -> Self::Output {
        self + &(-other)
    }
}
impl std::ops::Sub<&Point> for &Point {
    type Output = Vector;
    fn sub(self, other: &Point) -> Self::Output {
        let n = self.coords.len();
        if n == other.coords.len() {
            Vector {
                coords: (0..n).map(|i| &self.coords[i] - &other.coords[i]).collect(),
            }
        } else {
            panic!("Can only add vectors in the same dimension")
        }
    }
}

impl Vector {
    pub fn new(coords: Vec<Rational>) -> Self {
        Self { coords }
    }

    pub fn dim(&self) -> usize {
        self.coords.len()
    }

    pub fn get_coord(&self, i: usize) -> Rational {
        if i < self.coords.len() {
            self.coords[i].clone()
        } else {
            Rational::from(0)
        }
    }

    pub fn as_point(self) -> Point {
        Point {
            coords: self.coords,
        }
    }

    pub fn as_matrix(&self) -> Matrix<Rational> {
        Matrix::construct(self.coords.len(), 1, |r, _c| self.coords[r].clone())
    }

    pub fn from_matrix(mat: &Matrix<Rational>) -> Self {
        assert_eq!(mat.cols(), 1);
        Self::new(
            (0..mat.rows())
                .map(|r| mat.at(r, 0).unwrap().clone())
                .collect(),
        )
    }

    pub fn dot(&self, other: &Self) -> Rational {
        let dim = self.dim();
        assert_eq!(dim, other.dim());
        let mut ans = Rational::from(0);
        for i in 0..dim {
            ans += &self.coords[i] * &other.coords[i];
        }
        ans
    }

    pub fn length_sq(&self) -> Rational {
        self.dot(self)
    }
}

impl Point {
    pub fn new(coords: Vec<Rational>) -> Self {
        Self { coords }
    }

    pub fn dim(&self) -> usize {
        self.coords.len()
    }

    pub fn get_coord(&self, i: usize) -> Rational {
        if i < self.coords.len() {
            self.coords[i].clone()
        } else {
            Rational::from(0)
        }
    }

    pub fn as_vector(self) -> Vector {
        Vector {
            coords: self.coords,
        }
    }

    pub fn as_matrix(&self) -> Matrix<Rational> {
        Matrix::construct(self.coords.len(), 1, |r, _c| self.coords[r].clone())
    }

    pub fn from_matrix(mat: &Matrix<Rational>) -> Self {
        assert_eq!(mat.cols(), 1);
        Self::new(
            (0..mat.rows())
                .map(|r| mat.at(r, 0).unwrap().clone())
                .collect(),
        )
    }

    #[deprecated]
    fn transform(self, f: &dyn Fn(Point) -> Point) -> Self {
        f(self)
    }
}

impl PartialOrd for Point {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        let dim = self.coords.len();
        if dim != other.coords.len() {
            return None;
        }
        for i in 0..dim {
            match self.coords[i].cmp(&other.coords[i]) {
                std::cmp::Ordering::Less => {
                    return Some(std::cmp::Ordering::Less);
                }
                std::cmp::Ordering::Equal => {}
                std::cmp::Ordering::Greater => {
                    return Some(std::cmp::Ordering::Greater);
                }
            }
        }
        Some(std::cmp::Ordering::Equal)
    }
}

impl Ord for Point {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match self.partial_cmp(other) {
            Some(ans) => ans,
            None => panic!("Cant compare points in different dimensions"),
        }
    }
}

fn are_points_nondegenerage(dim: usize, points: &Vec<Point>) -> bool {
    for point in points {
        debug_assert_eq!(dim, point.dim());
    }
    if points.len() >= 1 {
        let root = points[0].clone();
        let mut vecs = vec![];
        for i in 1..points.len() {
            vecs.push(&points[i] - &root);
        }
        let mat = Matrix::construct(dim, vecs.len(), |r, c| vecs[c].get_coord(r));
        if QQ_MAT.rank(mat) != vecs.len() {
            return false;
        }
    }
    true
}

#[derive(Debug, Clone)]
struct AffineSubspace {
    dim: usize,         //the dimension of the ambient space
    origin: Point,      //the origin of the affine subspace in the ambient space
    basis: Vec<Vector>, //a basis for the affine subspace relative to the origin
}

impl AffineSubspace {
    pub fn check(&self) -> Result<(), &'static str> {
        if self.dim != self.origin.dim() {
            return Err("origin of affine subspace should live in the ambient space");
        }
        for v in &self.basis {
            if self.dim != self.origin.dim() {
                return Err("basis vector of affine subspace should live in the ambient space");
            }
        }

        if self.basis.len() != QQ_MAT.rank(self.basis_matrix()) {
            return Err("affine subspace vectors should be linearly independent");
        }

        Ok(())
    }

    fn basis_matrix(&self) -> Matrix<Rational> {
        Matrix::construct(self.dim, self.basis.len(), |i, j| {
            self.basis[j].get_coord(i)
        })
    }

    fn cannonical_subspace(k: usize, n: usize) -> Self {
        // k <= n
        // The cannonical embedding of A^k into A^n by taking the first k out of n coordinates
        assert!(k <= n);
        Self {
            dim: n,
            origin: Point {
                coords: (0..n).map(|_i| Rational::from(0)).collect(),
            },
            basis: (0..k)
                .map(|i| Vector {
                    coords: (0..n)
                        .map(|j| match i == j {
                            true => Rational::from(1),
                            false => Rational::from(0),
                        })
                        .collect(),
                })
                .collect(),
        }
    }

    pub fn affine_span(&self) -> AffineLattice<Rational> {
        QQ_AFFLAT.from_offset_and_linear_lattice(
            self.dim,
            1,
            self.origin.as_matrix(),
            QQ_LINLAT.from_basis(
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

    pub fn point_image(&self, p: &Point) -> Point {
        assert_eq!(p.dim(), self.rank());
        &self.origin
            + &Vector::from_matrix(
                &QQ_MAT
                    .mul_refs(&self.basis_matrix(), &p.as_matrix())
                    .unwrap(),
            )
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

    pub fn point_preimage(&self, p: &Point) -> Option<Point> {
        assert_eq!(p.dim(), self.dim());

        match QQ_MAT.col_solve(&self.basis_matrix(), &(p - &self.origin).as_matrix()) {
            Some(sol) => Some(Point::from_matrix(&sol)),
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
        Some(Simplex::new(self.dim(), new_points))
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
        Some(Shape::new(self.dim(), new_simplices))
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
        Some(SimplicialComplex::new(self.dim(), new_simplices))
    }

    fn extend_basis(mut self, new_basis_vector: Vector) -> Self {
        self.basis.push(new_basis_vector);
        debug_assert!(self.check().is_ok());
        self
    }
}

struct Interval {
    start: Rational,
    end: Rational,
}

impl Interval {
    fn intersects(&self, other: &Self) -> bool {
        //do the closures intersect?
        self.start <= other.end && other.start <= self.end
    }
}

struct Box {
    intervals: Vec<Interval>,
}

impl Box {
    fn bounding_box(dim: usize, points: &Vec<Point>) -> Self {
        debug_assert_ne!(points.len(), 0);
        for point in points {
            debug_assert_eq!(dim, point.dim());
        }
        Self {
            intervals: (0..dim)
                .map(|i| Interval {
                    start: points.iter().map(|p| p.get_coord(i)).min().unwrap(),
                    end: points.iter().map(|p| p.get_coord(i)).max().unwrap(),
                })
                .collect(),
        }
    }

    fn intersects(&self, other: &Self) -> bool {
        let dim = self.intervals.len();
        assert_eq!(dim, other.intervals.len());
        (0..dim).any(|i| self.intervals[i].intersects(&other.intervals[i]))
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Simplex {
    dim: usize,
    vertices: Vec<Point>, //ordered
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

        if !are_points_nondegenerage(self.dim, &self.vertices) {
            return Err("Simplex is degenerate");
        }

        Ok(())
    }

    pub fn new(dim: usize, mut vertices: Vec<Point>) -> Self {
        vertices.sort();
        let ans = Self { dim, vertices };
        ans.check().unwrap();
        ans
    }

    pub fn new_standard_simplex(rank: usize) -> Self {
        //the simplex whose verticies are 0 and the standard basis vectors
        let mut vertices = vec![];
        vertices.push(Point::new((0..rank).map(|_j| Rational::from(0)).collect()));
        for i in 0..rank {
            vertices.push(Point::new(
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

    pub fn try_new(dim: usize, mut vertices: Vec<Point>) -> Option<Self> {
        if are_points_nondegenerage(dim, &vertices) {
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

    pub fn points(&self) -> Vec<Point> {
        self.vertices.clone()
    }

    fn bounding_box(&self) -> Box {
        Box::bounding_box(self.dim, &self.vertices)
    }

    #[deprecated]
    fn transform(self, new_dim: usize, f: &dyn Fn(Point) -> Point) -> Self {
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

    pub fn has_vertex(&self, pt: &Point) -> bool {
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

    fn oriented_facet(&self, k: usize) -> OrientedSimplex {
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

    fn oriented_facets(&self) -> Vec<OrientedSimplex> {
        assert_eq!(self.dim, self.rank().unwrap());
        (0..self.n()).map(|k| self.oriented_facet(k)).collect()
    }

    fn boundary_simplices(&self) -> Vec<Simplex> {
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
            QQ_AFFLAT.empty(self.dim, 1)
        } else {
            QQ_AFFLAT.from_offset_and_linear_lattice(
                self.dim,
                1,
                self.vertices[0].as_matrix(),
                QQ_LINLAT.from_basis(
                    self.dim,
                    1,
                    (1..self.vertices.len())
                        .map(|i| (&self.vertices[i] - &self.vertices[0]).as_matrix())
                        .collect(),
                ),
            )
        }
    }

    fn affine_subspace(&self) -> AffineSubspace {
        //the affine subspace coordinate system in which self is the standard simplex
        AffineSubspace {
            dim: self.dim,
            origin: self.vertices[0].clone(),
            basis: (1..self.vertices.len())
                .map(|i| &self.vertices[i] - &self.vertices[0])
                .collect(),
        }
    }

    pub fn centroid(&self) -> Point {
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

        Point { coords }
    }

    fn extend_by_point(mut self, point: Point) -> Self {
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

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct OrientedSimplex {
    simplex: Simplex, //every simplex has a natural orientation induced by the ordering of the points of the space
    flip: bool,       //whether the orientation is the natural one or its opposite
}

impl OrientedSimplex {
    pub fn check(&self) -> Result<(), &'static str> {
        self.simplex.check()?;

        if self.simplex.n() != self.simplex.dim {
            return Err(
                "OrientedSimplex should have dimension one less than the space it lives in",
            );
        }

        Ok(())
    }

    pub fn new(simplex: Simplex) -> Self {
        assert_eq!(simplex.dim, simplex.n());
        let ans = Self {
            simplex,
            flip: false,
        };
        ans.check().unwrap();
        ans
    }

    fn as_simplex(self) -> Simplex {
        self.simplex
    }

    fn from_simplex(simplex: Simplex, neg_pt: &Point) -> Self {
        assert_eq!(simplex.n(), simplex.dim());

        let ans = Self {
            simplex: simplex,
            flip: false,
        };

        debug_assert!(ans.check().is_ok());

        match ans.sign_point(&neg_pt) {
            std::cmp::Ordering::Less => ans,
            std::cmp::Ordering::Equal => panic!(),
            std::cmp::Ordering::Greater => ans.flipped(),
        }
    }

    fn from_points(dim: usize, points: Vec<Point>, neg_pt: &Point) -> Self {
        for point in &points {
            debug_assert_eq!(dim, point.dim());
        }
        debug_assert_eq!(dim, neg_pt.dim());
        debug_assert_eq!(dim, points.len());

        Self::from_simplex(Simplex::new(dim, points), neg_pt)
    }

    pub fn dim(&self) -> usize {
        self.simplex.dim
    }

    pub fn flipped(self) -> Self {
        Self {
            simplex: self.simplex,
            flip: !self.flip,
        }
    }

    pub fn det_point(&self, point: &Point) -> Rational {
        debug_assert_eq!(point.dim(), self.dim());
        if self.dim() == 0 {
            Rational::from(0)
        } else {
            self.det_vector(&(point - &self.simplex.vertices[0]))
        }
    }

    pub fn det_vector(&self, vec: &Vector) -> Rational {
        //vector relative to self.root
        let mat = Matrix::construct(self.simplex.dim, self.simplex.dim, |r, c| {
            if c == self.simplex.dim - 1 {
                vec.get_coord(r)
            } else {
                (&self.simplex.vertices[c + 1] - &self.simplex.vertices[0]).get_coord(r)
            }
        });
        let d = QQ_MAT.det(mat).unwrap();
        match self.flip {
            false => d,
            true => -d,
        }
    }

    pub fn sign_point(&self, point: &Point) -> std::cmp::Ordering {
        self.det_point(point).cmp(&Rational::from(0))
    }

    // pub fn rank(&self) -> usize {
    //     self.vecs.len()
    // }
}

//disjoint union of simplices
#[derive(Debug, Clone)]
pub struct Shape {
    dim: usize,
    simplices: Vec<Simplex>,
}

pub fn shape_disjoint_union(dim: usize, shapes: Vec<Shape>) -> Shape {
    let mut simplices = vec![];
    for mut shape in shapes {
        assert_eq!(shape.dim, dim);
        simplices.append(&mut shape.simplices);
    }
    let ans = Shape { dim, simplices };
    debug_assert!(ans.check().is_ok());
    ans
}

impl Shape {
    pub fn check(&self) -> Result<(), &'static str> {
        for simplex in &self.simplices {
            if !simplex.rank().is_some() {
                return Err("simplex in a shape musn't be the empty simplex");
            }

            if simplex.dim != self.dim {
                return Err("Simplex dim does not match shape dim");
            }
            simplex.check()?;
        }
        //TODO: check that the simplices are disjoint
        Ok(())
    }

    pub fn new(dim: usize, simplices: Vec<Simplex>) -> Self {
        let ans = Self { dim, simplices };
        debug_assert!(ans.check().is_ok());
        ans
    }

    pub fn dim(&self) -> usize {
        self.dim
    }

    pub fn empty(dim: usize) -> Self {
        Self {
            dim,
            simplices: vec![],
        }
    }

    pub fn is_empty(&self) -> bool {
        self.simplices.is_empty()
    }

    pub fn simplex(simplex: Simplex) -> Self {
        Self {
            dim: simplex.dim,
            simplices: vec![simplex],
        }
    }

    pub fn simplices(self) -> Vec<Simplex> {
        self.simplices
    }

    pub fn simplices_ref(&self) -> Vec<&Simplex> {
        self.simplices.iter().collect()
    }

    #[deprecated]
    fn transform(self, new_dim: usize, f: &dyn Fn(Point) -> Point) -> Self {
        Self::new(
            new_dim,
            self.simplices
                .into_iter()
                .map(|s| s.transform(new_dim, f))
                .collect(),
        )
    }

    fn points(&self) -> Vec<Point> {
        let mut all_points = vec![];
        for s in &self.simplices {
            for p in &s.vertices {
                all_points.push(p.clone());
            }
        }
        all_points
    }

    pub fn is_complete(&self) -> bool {
        let all_simplices: HashSet<_> = self.simplices.iter().collect();
        for s in &self.simplices {
            if s.n() >= 2 {
                for f in s.facets() {
                    if !all_simplices.contains(&f) {
                        return false;
                    }
                }
            }
        }
        true
    }

    pub fn completion(self) -> Shape {
        //does not always produce a valid shape
        let mut all_simplices = HashSet::new();
        for s in self.simplices {
            for b in s.boundary().simplices() {
                all_simplices.insert(b);
            }
            all_simplices.insert(s);
        }
        let ans = Shape {
            dim: self.dim,
            simplices: all_simplices.into_iter().collect(),
        };
        debug_assert!(ans.check().is_ok());
        ans
    }

    pub fn boundary(&self) -> Shape {
        let mut facets = HashSet::new();
        for simplex in &self.simplices {
            for facet in simplex.facets() {
                if facets.contains(&facet) {
                    facets.remove(&facet);
                } else {
                    facets.insert(facet);
                }
            }
        }
        let shape = Shape {
            dim: self.dim,
            simplices: facets.into_iter().collect(),
        };
        debug_assert!(shape.check().is_ok());
        shape
    }

    pub fn equal(&self, other: &Self) -> bool {
        let dim = self.dim;
        assert_eq!(dim, other.dim);
        self.symmetric_difference_nosimp(other).is_empty()
    }

    fn simplify(&self) -> Self {
        self.clone() //TODO
    }

    fn intersect_nosimp(&self, other: &Self) -> Self {
        let dim = self.dim;
        assert_eq!(dim, other.dim);
        let mut parts = vec![];
        for s in &self.simplices {
            for t in &other.simplices {
                let (inner, _outer) = cut_simplex_by_simplex(&s, &t);
                parts.push(inner);
            }
        }
        shape_disjoint_union(dim, parts)
    }

    fn union_nosimp(&self, other: &Self) -> Self {
        let dim = self.dim;
        assert_eq!(dim, other.dim);
        shape_disjoint_union(
            dim,
            vec![
                self.intersect_nosimp(other),
                self.symmetric_difference_nosimp(other),
            ],
        )
        // shape_union(dim, vec![self.subtract_nosimp(other), other.clone()])
    }

    fn subtract_nosimp(&self, other: &Self) -> Self {
        let dim = self.dim;
        assert_eq!(dim, other.dim);
        let mut ans = self.clone();
        for s in &other.simplices {
            let (_inner, outer) = cut_shape_by_simplex(&s, &ans);
            ans = outer;
        }
        ans
    }

    fn symmetric_difference_nosimp(&self, other: &Self) -> Self {
        let dim = self.dim;
        assert_eq!(dim, other.dim);
        shape_disjoint_union(
            dim,
            vec![self.subtract_nosimp(other), other.subtract_nosimp(self)],
        )
    }

    pub fn intersect(&self, other: &Self) -> Self {
        self.intersect_nosimp(other).simplify()
    }
    pub fn union(&self, other: &Self) -> Self {
        self.union_nosimp(other).simplify()
    }
    pub fn subtract(&self, other: &Self) -> Self {
        self.subtract_nosimp(other).simplify()
    }
    pub fn symmetric_difference(&self, other: &Self) -> Self {
        self.symmetric_difference_nosimp(other).simplify()
    }
}

fn interior_of_convex_shell(shape: &Shape) -> Shape {
    debug_assert!(shape.is_complete());

    //shape should be the full boundary of a convex polytope

    if shape.simplices.is_empty() {
        return Shape::empty(shape.dim);
    }

    let n = shape.simplices.iter().map(|s| s.n()).max().unwrap();
    //n=1: shape is the shell of a line i.e. two points
    //n=2: shape is the shell of a polygon
    //n=3: shape is the shell of a solid
    //...

    if n == 1 {
        for s in &shape.simplices {
            debug_assert_eq!(s.n(), 1);
        }

        //every 1-dimensional convex shape must either be empty, a point, or a line
        match shape.simplices.len() {
            0 => Shape::empty(shape.dim), //no points has empty interior
            1 => Shape::empty(shape.dim), //single point has empty interior
            2 => Shape::simplex(Simplex::new(
                shape.dim,
                vec![
                    shape.simplices[0].vertices[0].clone(),
                    shape.simplices[1].vertices[0].clone(),
                ],
            )), //pair of points has interior of a line
            _ => {
                panic!()
            }
        }
    } else {
        debug_assert!(n >= 2);

        //n = 2: interior of a polygon
        //n = 3: interior of a solid
        //n = 4: interior of a 4d solid
        //n = 5: ...

        //the idea is to pick a vertex and add simplexes in a fan pattern from that vertex to disjoint bits of the boundary
        //e.g. for a pentagon: add 3 triangles, and 2 edges
        //e.g. for an icosahedron: add 15 tetrahedra, 25 triangles, and 6 edges

        if shape.simplices.len() == 0 {
            Shape::empty(shape.dim)
        } else {
            // println!("CUT SIMPLEX");

            //1) choose a base point
            let root = &shape.simplices[0].vertices[0];
            // println!("root = {:?}", root);
            // for s in &shape.simplices {
            //     println!("shell: {:?}", s);
            // }

            //2) find all cells adjacent to root
            let adj_cell_spans: Vec<_> = shape
                .simplices
                .iter()
                .filter(|s| s.has_vertex(root))
                .map(|s| s.affine_span())
                .collect();

            //3) find which adjacent faces are vertex belongs to
            let mut point_degen = HashMap::new();
            for s in &shape.simplices {
                if s.n() == 1 {
                    let pt = &s.vertices[0];
                    debug_assert!(!point_degen.contains_key(pt));
                    let adj_cells: HashSet<usize> = adj_cell_spans
                        .iter()
                        .enumerate()
                        .filter(|(_idx, adj_cell_span)| {
                            QQ_AFFLAT.contains_point(adj_cell_span, pt.as_matrix())
                        })
                        .map(|(idx, _adj_cell_span)| idx)
                        .collect();
                    point_degen.insert(pt, adj_cells);
                }
            }
            //check that every vertex of every part of the shell is present - it should be since the shell is complete
            for s in &shape.simplices {
                for p in &s.vertices {
                    debug_assert!(point_degen.contains_key(p));
                }
            }
            // println!("{:?}", point_degen);

            //4) for each shell simplex, if not all vertices lie in some adjacent face span, fill it in
            let mut interior = Shape::empty(shape.dim);
            for s in &shape.simplices {
                debug_assert!(s.vertices.len() >= 1);
                let mut common = point_degen[&s.vertices[0]].clone();
                common = common
                    .into_iter()
                    .filter(|idx| {
                        (1..s.vertices.len())
                            .map(|i| &s.vertices[i])
                            .all(|pt| point_degen[pt].contains(idx))
                    })
                    .collect();
                if common.len() == 0 {
                    let filler = Simplex::new(shape.dim, {
                        let mut filler_pts = vec![root.clone()];
                        filler_pts.append(&mut s.vertices.clone());
                        // println!("filler_pts = {:?}", filler_pts);
                        filler_pts
                    });
                    // println!("filler simplex: {:?}", filler);
                    interior =
                        shape_disjoint_union(shape.dim, vec![interior, Shape::simplex(filler)]);
                }
            }
            interior
        }
    }
}

fn cut_simplex_by_plane(cut_plane: &OrientedSimplex, simplex: &Simplex) -> (Shape, Shape, Shape) {
    let dim = cut_plane.dim();
    assert_eq!(dim, simplex.dim);
    debug_assert!(simplex.n() <= dim + 1);

    match simplex.n() {
        0 => (Shape::empty(dim), Shape::empty(dim), Shape::empty(dim)),
        1 => match cut_plane.sign_point(&simplex.vertices[0]) {
            std::cmp::Ordering::Greater => (
                Shape::simplex(simplex.clone()),
                Shape::empty(dim),
                Shape::empty(dim),
            ),
            std::cmp::Ordering::Equal => (
                Shape::empty(dim),
                Shape::simplex(simplex.clone()),
                Shape::empty(dim),
            ),
            std::cmp::Ordering::Less => (
                Shape::empty(dim),
                Shape::empty(dim),
                Shape::simplex(simplex.clone()),
            ),
        },
        2 => {
            let (p, q) = (&simplex.vertices[0], &simplex.vertices[1]);
            match (cut_plane.sign_point(p), cut_plane.sign_point(q)) {
                (std::cmp::Ordering::Greater, std::cmp::Ordering::Greater)
                | (std::cmp::Ordering::Equal, std::cmp::Ordering::Greater)
                | (std::cmp::Ordering::Greater, std::cmp::Ordering::Equal) => (
                    Shape::simplex(simplex.clone()),
                    Shape::empty(dim),
                    Shape::empty(dim),
                ),
                (std::cmp::Ordering::Equal, std::cmp::Ordering::Equal) => (
                    Shape::empty(dim),
                    Shape::simplex(simplex.clone()),
                    Shape::empty(dim),
                ),
                (std::cmp::Ordering::Less, std::cmp::Ordering::Less)
                | (std::cmp::Ordering::Equal, std::cmp::Ordering::Less)
                | (std::cmp::Ordering::Less, std::cmp::Ordering::Equal) => (
                    Shape::empty(dim),
                    Shape::empty(dim),
                    Shape::simplex(simplex.clone()),
                ),
                (p_sign, q_sign) => {
                    let t = -cut_plane.det_point(&p) / cut_plane.det_vector(&(q - p));
                    let r = p + &(&t * &(q - p));
                    debug_assert!(0 < t && t < 1);
                    debug_assert_eq!(cut_plane.det_point(&r), Rational::from(0));
                    match (p_sign, q_sign) {
                        (std::cmp::Ordering::Greater, std::cmp::Ordering::Less) => (
                            Shape::simplex(Simplex::new(dim, vec![p.clone(), r.clone()])),
                            Shape::simplex(Simplex::new(dim, vec![r.clone()])),
                            Shape::simplex(Simplex::new(dim, vec![r, q.clone()])),
                        ),
                        (std::cmp::Ordering::Less, std::cmp::Ordering::Greater) => (
                            Shape::simplex(Simplex::new(dim, vec![q.clone(), r.clone()])),
                            Shape::simplex(Simplex::new(dim, vec![r.clone()])),
                            Shape::simplex(Simplex::new(dim, vec![r, p.clone()])),
                        ),
                        _ => panic!(),
                    }
                }
            }
        }
        n => {
            debug_assert!(n >= 3);
            let (open_left, hollow_middle, open_right) =
                cut_shape_by_plane(&cut_plane, &simplex.boundary().as_shape());

            let interior_middle = interior_of_convex_shell(&hollow_middle);
            let hollow_left = shape_disjoint_union(
                dim,
                vec![
                    open_left.clone(),
                    hollow_middle.clone(),
                    interior_middle.clone(),
                ],
            );
            let hollow_right = shape_disjoint_union(
                dim,
                vec![
                    open_right.clone(),
                    hollow_middle.clone(),
                    interior_middle.clone(),
                ],
            );
            let interior_left = interior_of_convex_shell(&hollow_left);
            let interior_right = interior_of_convex_shell(&hollow_right);
            (interior_left, interior_middle, interior_right)
        }
    }
}

fn cut_shape_by_plane(cut_plane: &OrientedSimplex, shape: &Shape) -> (Shape, Shape, Shape) {
    let dim = cut_plane.dim();
    assert_eq!(dim, shape.dim);
    let (mut left, mut middle, mut right) =
        (Shape::empty(dim), Shape::empty(dim), Shape::empty(dim));
    for simplex in &shape.simplices {
        let (s_left, s_middle, s_right) = cut_simplex_by_plane(&cut_plane, &simplex);
        left = shape_disjoint_union(dim, vec![left, s_left]);
        middle = shape_disjoint_union(dim, vec![middle, s_middle]);
        right = shape_disjoint_union(dim, vec![right, s_right]);
    }
    (left, middle, right)
}

//return (part of simplex inside cut_simplex, part of simplex outside cut_simplex)
pub fn cut_simplex_by_simplex(cut_simplex: &Simplex, simplex: &Simplex) -> (Shape, Shape) {
    let dim = cut_simplex.dim;
    assert!(dim >= 1);
    assert_eq!(dim, simplex.dim);

    if !cut_simplex
        .bounding_box()
        .intersects(&simplex.bounding_box())
    {
        return (Shape::empty(dim), Shape::simplex(simplex.clone()));
    }

    //cut in some hyperplanes to restrict to the affine subspace of cut_simplex
    match cut_simplex.affine_span().elems() {
        AffineLatticeElements::Empty() => panic!(),
        AffineLatticeElements::NonEmpty { offset, linlat } => {
            let flat_dim = QQ_LINLAT.rank(linlat);
            let offset_pt = Point::from_matrix(offset);
            let offset_vec = Vector::from_matrix(offset);

            let linhyperplanes = QQ_LINLAT.as_hyperplane_intersection(linlat);

            let mut inside = Shape::simplex(simplex.clone());
            let mut outside = Shape::empty(dim);

            for linhyperplane in linhyperplanes {
                //construct an oriented simplex which spans linhyperplane
                let mut pts = vec![];
                pts.push(offset_pt.clone());
                for c in 0..dim - 1 {
                    pts.push(
                        &offset_pt
                            + &Vector::from_matrix(&QQ_LINLAT.basis_matrix(&linhyperplane, c)),
                    );
                }

                let cut_plane = OrientedSimplex {
                    simplex: Simplex::new(dim, pts),
                    flip: false,
                };
                debug_assert!(cut_plane.check().is_ok());

                let (inside_left, inside_mid, inside_right) =
                    cut_shape_by_plane(&cut_plane, &inside);
                // let (outside_left, outside_mid, outside_right) =
                //     cut_shape_by_plane(&cut_plane, &outside);

                inside = inside_mid;
                outside = shape_disjoint_union(dim, vec![inside_left, inside_right, outside]);
            }

            //restrict to the span of cut_simplex and intersect in there using the facets as oriented hyperplanes in the subspace
            let mat = Matrix::join_cols(
                dim,
                (0..flat_dim)
                    .map(|i| QQ_LINLAT.basis_matrix(linlat, i))
                    .collect(),
            );

            if flat_dim != 0 {
                //apply so that cut simplex has full rank in the affine subspace
                // mat(pt - offset)
                //apply to put back into this space
                // (mat(!Ǝ)=x) + offset

                let flat_transform = &|pt: Point| {
                    Point::from_matrix(
                        &QQ_MAT
                            .col_solve(&mat, &(&pt - &offset_vec).as_matrix())
                            .unwrap(),
                    )
                };
                let flat_cut_simplex = cut_simplex.clone().transform(flat_dim, flat_transform);
                let mut flat_inside = inside.transform(flat_dim, flat_transform);
                let mut flat_outside = Shape::empty(flat_dim);
                for k in 0..flat_cut_simplex.n() {
                    let (flat_inside_out, flat_inside_bdry, flat_inside_in) =
                        cut_shape_by_plane(&flat_cut_simplex.oriented_facet(k), &flat_inside);
                    // let (flat_outside_out, flat_outside_bdry, flat_outside_in) =
                    //     cut_shape_by_plane(&simplex.oriented_facet(k), &flat_outside);

                    flat_inside = flat_inside_in;
                    flat_outside = shape_disjoint_union(
                        flat_dim,
                        vec![flat_inside_out, flat_inside_bdry, flat_outside],
                    );
                }

                let unflat_transform = &|pt: Point| {
                    &Point::from_matrix(&QQ_MAT.mul_refs(&mat, &pt.as_matrix()).unwrap())
                        + &offset_vec
                };
                let unflat_inside = flat_inside.transform(dim, unflat_transform);
                let unflat_outside = flat_outside.transform(dim, unflat_transform);

                (inside, outside) = (
                    unflat_inside,
                    shape_disjoint_union(dim, vec![unflat_outside, outside]),
                );
            }

            if inside.is_empty() {
                (inside, Shape::simplex(simplex.clone()))
            } else if outside.is_empty() {
                (Shape::simplex(simplex.clone()), outside)
            } else {
                (inside, outside)
            }
        }
    }
}

//return (part of shape inside cut_simplex, part of shape outside cut_simplex)
pub fn cut_shape_by_simplex(cut_simplex: &Simplex, shape: &Shape) -> (Shape, Shape) {
    let dim = cut_simplex.dim;
    assert_eq!(dim, shape.dim);
    let (mut inside, mut outside) = (Shape::empty(dim), Shape::empty(dim));
    for simplex in &shape.simplices {
        let (s_inside, s_outside) = cut_simplex_by_simplex(&cut_simplex, &simplex);
        inside = shape_disjoint_union(dim, vec![inside, s_inside]);
        outside = shape_disjoint_union(dim, vec![outside, s_outside]);
    }
    (inside, outside)
}

/*
pub fn convexhull_boundary(dim: usize, points: Vec<Point>) -> (Shape, usize) {
    //return None if the convex hull is flat
    for point in &points {
        assert_eq!(point.dim(), dim);
    }

    if points.len() == 0 {
        panic!("points should be non-empty");
    } else if points.len() == 1 {
        (Shape::simplex(Simplex::new(dim, points)), 0)
    } else {
        let (first, rest) = points.split_first().unwrap();

        let mat = Matrix::construct(dim, rest.len(), |r, c| {
            rest[c].get_coord(r) - first.get_coord(r)
        });
        let (_h, _u, _u_det, pivs) = QQ_MAT.row_hermite_algorithm(mat);
        let rank = pivs.len();

        let starting_simplex = Simplex::new(dim, {
            let mut starting_simplex_points = vec![first.clone()];
            for i in pivs {
                starting_simplex_points.push(rest[i].clone());
            }
            starting_simplex_points
        });
        debug_assert!(starting_simplex.check().is_ok());

        debug_assert!(rank <= dim);
        if rank == dim {
            let middle_point = starting_simplex.centroid();
            let mut facets: HashSet<_> = starting_simplex.oriented_facets().into_iter().collect();
            for pt in points.iter() {
                let visible = facets
                    .iter()
                    .filter(|f| f.sign_point(pt).is_ge())
                    .map(|f| f.clone())
                    .collect_vec();

                let mut horizon = HashSet::new();
                for v in visible {
                    for b in v.simplex.facets() {
                        if horizon.contains(&b) {
                            horizon.remove(&b);
                        } else {
                            horizon.insert(b);
                        }
                    }
                    facets.remove(&v);
                }

                for h in horizon {
                    let mut new_facet_points = h.vertices;
                    new_facet_points.push(pt.clone());
                    let new_facet =
                        OrientedSimplex::from_points(dim, new_facet_points, &middle_point);
                    new_facet.check().unwrap();
                    facets.insert(new_facet);
                }
            }

            (
                Shape {
                    dim,
                    simplices: facets.into_iter().map(|facet| facet.simplex).collect(),
                }
                .completion(),
                rank,
            )
        } else {
            match starting_simplex.affine_span().elems() {
                AffineLatticeElements::Empty() => panic!(),
                AffineLatticeElements::NonEmpty { offset, linlat } => {
                    let flat_dim = QQ_LINLAT.rank(linlat);
                    let offset_pt = Point::from_matrix(offset);
                    let offset_vec = Vector::from_matrix(offset);

                    //restrict to the span of starting_simplex and take convex hull there
                    let mat = Matrix::join_cols(
                        dim,
                        (0..flat_dim)
                            .map(|i| QQ_LINLAT.basis_matrix(linlat, i))
                            .collect(),
                    );

                    if flat_dim == 0 {
                        panic!("convex hull of a single point should already be taken care of");
                    } else {
                        //apply so that cut simplex has full rank in the affine subspace
                        // mat(pt - offset)
                        //apply to put back into this space
                        // (mat(!Ǝ)=x) + offset

                        let flat_transform = &|pt: Point| {
                            Point::from_matrix(
                                &QQ_MAT
                                    .col_solve(&mat, &(&pt - &offset_vec).as_matrix())
                                    .unwrap(),
                            )
                        };
                        let flat_points = points
                            .iter()
                            .map(|pt| pt.clone().transform(flat_transform))
                            .collect();
                        let (flat_hull_boundary, rank) = convexhull_boundary(flat_dim, flat_points);
                        debug_assert_eq!(rank, flat_dim);

                        let unflat_transform = &|pt: Point| {
                            &Point::from_matrix(&QQ_MAT.mul_refs(&mat, &pt.as_matrix()).unwrap())
                                + &offset_vec
                        };
                        let unflat_hull_boundary =
                            flat_hull_boundary.transform(dim, unflat_transform);

                        (unflat_hull_boundary, flat_dim)
                    }
                }
            }
        }
    }
}

pub fn convexhull_interior(dim: usize, points: Vec<Point>) -> (Shape, usize) {
    for point in &points {
        assert_eq!(point.dim(), dim);
    }
    let (shell, rank) = convexhull_boundary(dim, points);
    (interior_of_convex_shell(&shell), rank)
}

pub fn convexhull(dim: usize, points: Vec<Point>) -> Shape {
    for point in &points {
        assert_eq!(point.dim(), dim);
    }
    let (shell, _rank) = convexhull_boundary(dim, points);
    let interior = interior_of_convex_shell(&shell);
    shape_disjoint_union(dim, vec![shell, interior])
}
*/

#[derive(Debug, Clone)]
enum ConvexSimplicialComplexCases {
    NonEmpty {
        subspace: AffineSubspace, //an embedding of a vector space of dimension self.rank into the ambient space
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
            }) => subspace.shape_image(&Shape {
                dim: subspace.rank(),
                simplices: boundary.iter().map(|s| s.clone().as_simplex()).collect(),
            }),
            ConvexSimplicialComplex(ConvexSimplicialComplexCases::Empty { dim }) => {
                Shape::empty(*dim)
            }
        }
    }

    fn subspace_interior_point(&self) -> Option<Point> {
        match self {
            ConvexSimplicialComplex(ConvexSimplicialComplexCases::NonEmpty {
                subspace,
                boundary,
                interior,
            }) => Some(Simplex::new_standard_simplex(subspace.rank()).centroid()),
            ConvexSimplicialComplex(ConvexSimplicialComplexCases::Empty { dim }) => None,
        }
    }

    fn interior_point(&self) -> Option<Point> {
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

    ///Return the convex hull of self and point by extending self.
    pub fn extend_by_point(&self, point: &Point) -> Self {
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
                            for b in v.simplex.facets() {
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
                            extended_interior.push(Simplex {
                                dim: subspace.rank(),
                                vertices: vec![point_preimage],
                            });
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
                            subspace.clone().extend_basis(point - &subspace.origin);

                        let point_extended_preimage = Point {
                            coords: (0..subspace.rank() + 1)
                                .map(|i| match i == subspace.rank() {
                                    true => Rational::from(1),
                                    false => Rational::from(0),
                                })
                                .collect(),
                        };

                        let subspace_into_extended_subspace = AffineSubspace::cannonical_subspace(
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
                        extended_interior.push(Simplex {
                            dim: subspace.rank() + 1,
                            vertices: vec![point_extended_preimage],
                        });

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

pub fn convexhull(dim: usize, points: Vec<Point>) -> ConvexSimplicialComplex {
    for point in &points {
        assert_eq!(dim, point.dim());
    }
    let mut ch = ConvexSimplicialComplex::new_empty(dim);
    for point in points {
        println!("{:?}", point);
        ch = ch.extend_by_point(&point);
    }
    ch
}

// #[derive(Debug, Clone)]
// pub enum ConvexSimplicialComplex {
//     NonEmpty(NonemptyConvexSimplicialComplex),
//     Empty { dim: usize },
// }

// impl ConvexSimplicialComplex {
//     pub fn new_empty(dim: usize) -> Self {
//         ConvexSimplicialComplex::Empty{dim}
//     }

//     pub fn dim(&self) -> usize {
//         match self {
//             ConvexSimplicialComplex::NonEmpty(necsc) => necsc.dim(),
//             ConvexSimplicialComplex::Empty{dim} => *dim,
//         }
//     }

//     pub fn extend_by_point(&self, point: &Point) -> Self {
//         assert_eq!(self.dim(), point.dim());
//         match self {
//             ConvexSimplicialComplex::NonEmpty(necsc) => {
//                 ConvexSimplicialComplex::NonEmpty(necsc.extend_by_point(point))
//             }
//             ConvexSimplicialComplex::Empty(dim) => ConvexSimplicialComplex::NonEmpty(
//                 NonemptyConvexSimplicialComplex::new_from_simplex(&Simplex {
//                     dim: *dim,
//                     vertices: vec![point.clone()],
//                 }),
//             ),
//         }
//     }

//     #[deprecated(note = "for debug use only")]
//     pub fn interior_as_shape(&self) -> Shape {
//         match self {
//             ConvexSimplicialComplex::NonEmpty(nescs) => nescs.interior_as_shape(),
//             ConvexSimplicialComplex::Empty(dim) => Shape::empty(*dim),
//         }
//     }

//     #[deprecated(note = "for debug use only")]
//     pub fn boundary_as_shape(&self) -> Shape {
//         match self {
//             ConvexSimplicialComplex::NonEmpty(nescs) => nescs.boundary_as_shape(),
//             ConvexSimplicialComplex::Empty(dim) => Shape::empty(*dim),
//         }
//     }
// }

// #[derive(Debug, Clone)]
// pub struct PartialSimplicialComplex {
//     dim: usize,
//     simplices: Vec<Simplex>,
// }

// impl PartialSimplicialComplex {
//     pub fn as_shape(self) -> SimplexDisjointUnion {
//         SimplexDisjointUnion {
//             dim: self.dim,
//             simplices: self.simplices,
//         }
//     }

//     pub fn check(&self) -> Result<(), &'static str> {
//         self.clone().as_shape().check()?;
//         //TODO: check that every face of a simplex is a simplex
//         Ok(())
//     }

//     pub fn empty(dim: usize) -> Self {
//         Self {
//             dim,
//             simplices: vec![],
//         }
//     }

//     pub fn simplex(simplex: Simplex) -> Self {
//         Self {
//             dim: simplex.dim,
//             simplices: vec![simplex],
//         }
//     }

//     pub fn new(dim: usize, simplices: Vec<Simplex>) -> Self {
//         let ans = Self { dim, simplices };
//         debug_assert!(ans.check().is_ok());
//         ans
//     }

//     pub fn dim(&self) -> usize {
//         self.dim
//     }

//     pub fn is_empty(&self) -> bool {
//         self.simplices.is_empty()
//     }

//     pub fn simplices(self) -> Vec<Simplex> {
//         self.simplices
//     }
// }

#[derive(Debug, Clone)]
pub struct SimplicialComplex {
    dim: usize,
    simplices: Vec<Simplex>,
    point_simplex_lookup: HashMap<Point, Vec<Simplex>>,
}

impl SimplicialComplex {
    pub fn as_shape(self) -> Shape {
        Shape {
            dim: self.dim,
            simplices: self.simplices,
        }
    }

    pub fn check(&self) -> Result<(), &'static str> {
        self.clone().as_shape().check()?; //simplices are disjoint

        let simplices = self.simplices.iter().collect::<HashSet<_>>();
        for simplex in &self.simplices {
            for b in simplex.boundary().simplices_ref() {
                if !simplices.contains(b) {
                    return Err("simplicial complex is missing the boundary of some simplex");
                }
            }
        }

        //check self.point_simplex_lookup
        for (point, point_simplices) in &self.point_simplex_lookup {
            if point.dim() != self.dim() {
                return Err("point dim does not match self dim");
            }

            for simplex in point_simplices {
                if !simplices.contains(simplex) {
                    return Err(
                        "point simplex lookup contains a simplex which is not part of self",
                    );
                }

                if !simplex
                    .points()
                    .iter()
                    .collect::<HashSet<_>>()
                    .contains(point)
                {
                    return Err(
                        "point simplex lookup result is a simplex not containing the point",
                    );
                }
            }
        }

        for simplex in &self.simplices {
            for point in simplex.points() {
                if !self
                    .point_simplex_lookup
                    .get(&point)
                    .unwrap()
                    .iter()
                    .collect::<HashSet<_>>()
                    .contains(simplex)
                {
                    return Err("point simplex lookup is missing a point of a simplex");
                }
            }
        }

        Ok(())
    }

    pub fn empty(dim: usize) -> Self {
        Self {
            dim,
            simplices: vec![],
            point_simplex_lookup: HashMap::new(),
        }
    }

    pub fn new(dim: usize, simplices: Vec<Simplex>) -> Self {
        let mut point_simplex_lookup: HashMap<_, Vec<_>> = HashMap::new();

        for simplex in &simplices {
            for point in simplex.points() {
                if point_simplex_lookup.contains_key(&point) {
                    point_simplex_lookup
                        .get_mut(&point)
                        .unwrap()
                        .push(simplex.clone());
                } else {
                    point_simplex_lookup.insert(point, vec![simplex.clone()]);
                }
            }
        }

        let ans = Self {
            dim,
            simplices,
            point_simplex_lookup,
        };
        debug_assert!(ans.check().is_ok());
        ans
    }

    pub fn dim(&self) -> usize {
        self.dim
    }

    pub fn is_empty(&self) -> bool {
        self.simplices.is_empty()
    }

    pub fn simplices(self) -> Vec<Simplex> {
        self.simplices
    }

    pub fn simplices_ref(&self) -> Vec<&Simplex> {
        self.simplices.iter().collect()
    }

    // fn subset(&self) -> SubSimplicialComplex {
    //     SubSimplicialComplex {
    //         simplicial_complex: self,
    //         subset: vec![&self.simplices[0]],
    //     }
    // }

    pub fn interior_and_boundary(&self) -> (SubSimplicialComplex, SubSimplicialComplex) {
        /*
        let n be the dimension of the space self is living in
         - every simplex of rank n is part of the interior
         - a simplex of rank n-1 is the facet of at most 2 simplices of rank n, and is part of the interior if and only if it is the facet of exactly 2 simplices of rank n
         - a simplex of rank less or equal to n-2 is part of the interior iff it is in the boundary of some strictly higher rank simplex AND every strictly higher rank simplex containing it as part of the boundary is part of the interior
        */

        //each simplex x points at all simplicies y such that x is part of the boundary of y
        let mut inverse_boundary_lookup: HashMap<&Simplex, HashSet<&Simplex>> = self
            .simplices_ref()
            .into_iter()
            .map(|s| (s, HashSet::new()))
            .collect();
        for s in self.simplices_ref() {
            for b in s.boundary_simplices() {
                inverse_boundary_lookup.get_mut(&b).unwrap().insert(s);
            }
        }

        let n = self.dim();

        let mut interior = HashSet::new();
        let mut boundary = HashSet::new();

        let mut all = self.simplices.iter().collect_vec();
        all.sort_by_key(|s| std::cmp::Reverse(s.rank().unwrap())); //so that we process largest rank first
        for simplex in all {
            let r = simplex.rank().unwrap();
            if r == n {
                interior.insert(simplex);
            } else if r == n - 1 {
                match inverse_boundary_lookup.get(simplex).unwrap().len() {
                    0 | 1 => {
                        boundary.insert(simplex);
                    }
                    2 => {
                        interior.insert(simplex);
                    }
                    _ => panic!(
                        "rank n-1 simplex should be in the boundary of at most 2 rank n simplices"
                    ),
                }
            } else {
                debug_assert!(r < n - 1);
                let bdry = inverse_boundary_lookup.get(simplex).unwrap();
                if bdry.is_empty() {
                    boundary.insert(simplex);
                } else {
                    if bdry.iter().all(|b| interior.contains(*b)) {
                        interior.insert(simplex);
                    } else {
                        boundary.insert(simplex);
                    }
                }
            }
        }

        (
            SubSimplicialComplex {
                simplicial_complex: self,
                subset: interior,
            },
            SubSimplicialComplex {
                simplicial_complex: self,
                subset: boundary,
            },
        )
    }

    pub fn interior(&self) -> SubSimplicialComplex {
        self.interior_and_boundary().0
    }

    pub fn boundary(&self) -> SubSimplicialComplex {
        self.interior_and_boundary().1
    }

    // pub fn simplify(mut self) -> Self {
    //     println!("self = {:?}", self);

    //     for simplex in &self.simplices {
    //         println!("simplex = {:?}", simplex);
    //     }

    //     let mut all_points = HashSet::new();
    //     for simplex in &self.simplices {
    //         for point in simplex.points() {
    //             all_points.insert(point);
    //         }
    //     }

    //     for center_point in all_points {
    //         println!("center_point = {:?}", center_point);

    //         //compute all points adjacent to center_point
    //         let mut adj_points = HashSet::new();
    //         for simplex in &self.simplices {
    //             if simplex.points().iter().any(|pt| pt == &center_point) {
    //                 for adj_point in simplex.points() {
    //                     adj_points.insert(adj_point);
    //                 }
    //             }
    //         }

    //         println!("adj_points = {:?}", adj_points);
    //         let adj_points_convex_hull = convexhull(self.dim(), adj_points.into_iter().collect());

    //         println!("adj_points_convex_hull = {:?}", adj_points_convex_hull);
    //     }

    //     todo!();

    //     self
    // }
}

pub struct SubSimplicialComplex<'a> {
    simplicial_complex: &'a SimplicialComplex,
    subset: HashSet<&'a Simplex>,
}

impl<'a> SubSimplicialComplex<'a> {
    pub fn check(&self) -> Result<(), &'static str> {
        self.simplicial_complex.check()?;

        for simplex in &self.subset {
            if !self
                .simplicial_complex
                .simplices_ref()
                .into_iter()
                .any(|sc_simplex| std::ptr::eq(*simplex, sc_simplex))
            {
                return Err("simplices in the subset should reference simplices in the simplicial complex only");
            }
        }

        Ok(())
    }

    pub fn dim(&self) -> usize {
        self.simplicial_complex.dim()
    }

    pub fn as_shape(&self) -> Shape {
        Shape::new(
            self.dim(),
            self.subset.iter().map(|s| (*s).clone()).collect(),
        )
    }
}

#[cfg(test)]
mod geometry_tests {
    use std::str::FromStr;

    use super::*;

    #[test]
    fn test_cut_simplex() {
        let s = Simplex::new(
            3,
            vec![
                Point::new(vec![
                    Rational::from_str("0").unwrap(),
                    Rational::from_str("0").unwrap(),
                    Rational::from_str("1").unwrap(),
                ]),
                Point::new(vec![
                    Rational::from_str("0").unwrap(),
                    Rational::from_str("0").unwrap(),
                    Rational::from_str("-1").unwrap(),
                ]),
                Point::new(vec![
                    Rational::from_str("1").unwrap(),
                    Rational::from_str("2").unwrap(),
                    Rational::from_str("0").unwrap(),
                ]),
                Point::new(vec![
                    Rational::from_str("2").unwrap(),
                    Rational::from_str("1").unwrap(),
                    Rational::from_str("0").unwrap(),
                ]),
            ],
        );

        let h = OrientedSimplex::new(Simplex::new(
            3,
            vec![
                Point::new(vec![
                    Rational::from_str("1").unwrap(),
                    Rational::from_str("0").unwrap(),
                    Rational::from_str("0").unwrap(),
                ]),
                Point::new(vec![
                    Rational::from_str("1").unwrap(),
                    Rational::from_str("1").unwrap(),
                    Rational::from_str("0").unwrap(),
                ]),
                Point::new(vec![
                    Rational::from_str("1").unwrap(),
                    Rational::from_str("0").unwrap(),
                    Rational::from_str("1").unwrap(),
                ]),
            ],
        ));

        s.check().unwrap();
        h.check().unwrap();
        let (a, b, c) = cut_simplex_by_plane(&h, &s);
        a.check().unwrap();
        b.check().unwrap();
        c.check().unwrap();
        println!(
            "{:?}",
            a.clone()
                .simplices()
                .iter()
                .map(|s| s.n())
                .collect::<Vec<_>>()
        );
        println!(
            "{:?}",
            b.clone()
                .simplices()
                .iter()
                .map(|s| s.n())
                .collect::<Vec<_>>()
        );
        println!(
            "{:?}",
            c.clone()
                .simplices()
                .iter()
                .map(|s| s.n())
                .collect::<Vec<_>>()
        );
        let d = shape_disjoint_union(3, vec![a, b, c]);
        d.check().unwrap();
    }

    #[test]
    fn shape_equal() {
        let shape1 = Shape {
            dim: 2,
            simplices: vec![
                Simplex::new(
                    2,
                    vec![
                        Point::new(vec![
                            Rational::from_str("0").unwrap(),
                            Rational::from_str("0").unwrap(),
                        ]),
                        Point::new(vec![
                            Rational::from_str("1").unwrap(),
                            Rational::from_str("0").unwrap(),
                        ]),
                        Point::new(vec![
                            Rational::from_str("0").unwrap(),
                            Rational::from_str("1").unwrap(),
                        ]),
                    ],
                ),
                Simplex::new(
                    2,
                    vec![
                        Point::new(vec![
                            Rational::from_str("1").unwrap(),
                            Rational::from_str("1").unwrap(),
                        ]),
                        Point::new(vec![
                            Rational::from_str("1").unwrap(),
                            Rational::from_str("0").unwrap(),
                        ]),
                        Point::new(vec![
                            Rational::from_str("0").unwrap(),
                            Rational::from_str("1").unwrap(),
                        ]),
                    ],
                ),
                Simplex::new(
                    2,
                    vec![
                        Point::new(vec![
                            Rational::from_str("1").unwrap(),
                            Rational::from_str("0").unwrap(),
                        ]),
                        Point::new(vec![
                            Rational::from_str("0").unwrap(),
                            Rational::from_str("1").unwrap(),
                        ]),
                    ],
                ),
            ],
        };

        let shape2 = Shape {
            dim: 2,
            simplices: vec![
                Simplex::new(
                    2,
                    vec![
                        Point::new(vec![
                            Rational::from_str("1").unwrap(),
                            Rational::from_str("0").unwrap(),
                        ]),
                        Point::new(vec![
                            Rational::from_str("0").unwrap(),
                            Rational::from_str("0").unwrap(),
                        ]),
                        Point::new(vec![
                            Rational::from_str("1").unwrap(),
                            Rational::from_str("1").unwrap(),
                        ]),
                    ],
                ),
                Simplex::new(
                    2,
                    vec![
                        Point::new(vec![
                            Rational::from_str("0").unwrap(),
                            Rational::from_str("1").unwrap(),
                        ]),
                        Point::new(vec![
                            Rational::from_str("0").unwrap(),
                            Rational::from_str("0").unwrap(),
                        ]),
                        Point::new(vec![
                            Rational::from_str("1").unwrap(),
                            Rational::from_str("1").unwrap(),
                        ]),
                    ],
                ),
                Simplex::new(
                    2,
                    vec![
                        Point::new(vec![
                            Rational::from_str("0").unwrap(),
                            Rational::from_str("0").unwrap(),
                        ]),
                        Point::new(vec![
                            Rational::from_str("1").unwrap(),
                            Rational::from_str("1").unwrap(),
                        ]),
                    ],
                ),
            ],
        };

        assert!(shape1.equal(&shape2));
    }

    // #[test]
    // fn test_simplex_boundary_as_mesh() {
    //     let a = Simplex::new(
    //         2,
    //         vec![
    //             Point::new(vec![
    //                 Rational::from_str("0").unwrap(),
    //                 Rational::from_str("1").unwrap(),
    //             ]),
    //             Point::new(vec![
    //                 Rational::from_str("2").unwrap(),
    //                 Rational::from_str("1").unwrap(),
    //             ]),
    //             Point::new(vec![
    //                 Rational::from_str("1").unwrap(),
    //                 Rational::from_str("-1").unwrap(),
    //             ]),
    //         ],
    //     );

    //     debug_assert!(a.boundary_as_mesh().unwrap().check().is_ok());
    // }

    #[test]
    fn test_convex_simplicial_complex_point_extension() {
        let csc0 = ConvexSimplicialComplex::new_empty(2);
        assert!(csc0.check().is_ok());
        assert_eq!(csc0.rank(), None);

        let csc1 = csc0.extend_by_point(&Point::new(vec![
            Rational::from_str("1").unwrap(),
            Rational::from_str("1").unwrap(),
        ]));
        assert!(csc1.check().is_ok());
        assert_eq!(csc1.rank(), Some(0));

        let csc2 = csc1.extend_by_point(&Point::new(vec![
            Rational::from_str("1").unwrap(),
            Rational::from_str("1").unwrap(),
        ]));
        assert!(csc2.check().is_ok());
        assert_eq!(csc2.rank(), Some(0));

        let csc3 = csc2.extend_by_point(&Point::new(vec![
            Rational::from_str("2").unwrap(),
            Rational::from_str("3").unwrap(),
        ]));
        assert!(csc3.check().is_ok());
        assert_eq!(csc3.rank(), Some(1));

        let csc4 = csc3.extend_by_point(&Point::new(vec![
            Rational::from_str("4").unwrap(),
            Rational::from_str("7").unwrap(),
        ]));
        assert!(csc4.check().is_ok());
        assert_eq!(csc4.rank(), Some(1));

        let csc5 = csc4.extend_by_point(&Point::new(vec![
            Rational::from_str("3").unwrap(),
            Rational::from_str("5").unwrap(),
        ]));
        assert!(csc5.check().is_ok());
        assert_eq!(csc5.rank(), Some(1));

        let csc6 = csc5.extend_by_point(&Point::new(vec![
            Rational::from_str("0").unwrap(),
            Rational::from_str("6").unwrap(),
        ]));
        assert!(csc6.check().is_ok());
        assert_eq!(csc6.rank(), Some(2));

        let csc7 = csc6.extend_by_point(&Point::new(vec![
            Rational::from_str("1").unwrap(),
            Rational::from_str("0").unwrap(),
        ]));
        assert!(csc7.check().is_ok());
        assert_eq!(csc7.rank(), Some(2));

        let csc8 = csc7.extend_by_point(&Point::new(vec![
            Rational::from_str("4").unwrap(),
            Rational::from_str("1").unwrap(),
        ]));
        assert!(csc8.check().is_ok());
        assert_eq!(csc8.rank(), Some(2));
    }
}
