use core::panic;
#[allow(dead_code)]
use std::collections::{HashMap, HashSet};

use itertools::Itertools;
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

    pub fn skeleton(&self, skel_rank: usize) -> Vec<Simplex> {
        let mut parts = vec![];
        for skeleton_piece_index in (0..self.vertices.len()).combinations(skel_rank + 1) {
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
        self.skeleton(self.n() - 3)
    }

    pub fn facets(&self) -> Vec<Simplex> {
        self.skeleton(self.n() - 2)
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
        match oriented_facet.sign_point(&self.vertices[k]) {
            std::cmp::Ordering::Less => oriented_facet,
            std::cmp::Ordering::Equal => panic!(),
            std::cmp::Ordering::Greater => oriented_facet.flipped(),
        }
    }

    fn oriented_facets(&self) -> Vec<OrientedSimplex> {
        (0..self.n()).map(|k| self.oriented_facet(k)).collect()
    }

    pub fn boundary(&self) -> Shape {
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
        Shape {
            dim: self.dim,
            simplices: parts,
        }
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
        debug_assert_ne!(simplex.dim, 0);
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
        let ans = Self {
            simplex: simplex,
            flip: false,
        };

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
        self.det_vector(&(point - &self.simplex.vertices[0]))
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
                cut_shape_by_plane(&cut_plane, &simplex.boundary());

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

#[derive(Debug, Clone)]
pub struct ConvexSimplicialComplex {
    subspace: AffineSubspace, //an embedding of a vector space of dimension self.rank into the ambient space
    boundary: Vec<OrientedSimplex>, //oriented simplexes in the subspace of rank self.rank-1 with positive side on the outside
    interior: Vec<Simplex>,         //simplexes in the subspace which form the interior
                                    //furthermore, the convex simplicial complex should contain the embedding of the standard simplex in the AffineSubspace
}

impl ConvexSimplicialComplex {
    pub fn check(&self) -> Result<(), &'static str> {
        for simplex in &self.boundary {
            if simplex.dim() != self.subspace.rank() {
                return Err("boundary simplexes should live in the subspace coordinates");
            }
        }

        for simplex in &self.interior {
            if simplex.dim() != self.subspace.rank() {
                return Err("interior simplexes should live in the subspace coordinates");
            }

            for b in &self.boundary {
                if !b.sign_point(&simplex.centroid()).is_le() {
                    return Err("centroid of interior simplexes should be on the negative side of every boundary simplex");
                }
                if !b.sign_point(&self.subspace_interior_point()).is_le() {
                    return Err("centroid should contain the embedding of the standard simplex in the subspace, in particular it should contain its centroid");
                }
            }
        }

        Ok(())
    }

    pub fn as_shape(&self) -> Shape {
        let mut simplices = vec![];
        for simplex in &self.boundary {
            simplices.push(self.subspace.simplex_image(&simplex.clone().as_simplex()));
        }
        for simplex in &self.interior {
            simplices.push(self.subspace.simplex_image(simplex));
        }
        Shape::new(self.dim(), simplices)
    }

    fn subspace_interior_point(&self) -> Point {
        Simplex::new_standard_simplex(self.subspace.rank()).centroid()
    }

    fn interior_point(&self) -> Point {
        self.subspace.point_image(&self.subspace_interior_point())
    }

    pub fn new_from_simplex(simplex: &Simplex) -> Self {
        assert!(simplex.rank().is_some());
        let subspace = simplex.affine_subspace();
        let simplex_preimage = Simplex::new_standard_simplex(simplex.rank().unwrap());
        let boundary = simplex_preimage.oriented_facets();
        let interior = vec![simplex_preimage];

        let ans = Self {
            subspace,
            boundary,
            interior,
        };

        ans.check().unwrap();
        debug_assert!(ans.check().is_ok());
        ans
    }

    pub fn dim(&self) -> usize {
        self.subspace.dim()
    }

    pub fn rank(&self) -> usize {
        self.subspace.rank()
    }

    pub fn expand_by_point(&self, point: &Point) -> Self {
        println!("self = {:?}", self);

        match self.subspace.point_preimage(point) {
            Some(point_preimage) => {
                println!("flat case");
                println!("{:?}", point_preimage);

                todo!()
            }
            None => {
                // subspace C extended_subspace C ambient_space

                let extended_interior_point =
                    Simplex::new_standard_simplex(self.subspace.rank() + 1).centroid();

                let extended_subspace_into_ambient = self
                    .subspace
                    .clone()
                    .extend_basis(point - &self.subspace.origin);

                let point_extended_preimage = Point {
                    coords: (0..self.subspace.rank() + 1)
                        .map(|i| match i == self.subspace.rank() {
                            true => Rational::from(1),
                            false => Rational::from(0),
                        })
                        .collect(),
                };

                let subspace_into_extended_subspace = AffineSubspace::cannonical_subspace(
                    self.subspace.rank(),
                    self.subspace.rank() + 1,
                );

                println!("{:?}", extended_subspace_into_ambient);
                println!("{:?}", subspace_into_extended_subspace);

                println!("{:?}", point);
                println!("{:?}", point_extended_preimage);

                let mut extended_boundary = vec![];
                let mut extended_interior = vec![];

                for boundary_simplex in &self.boundary {
                    // add its cone to the new boundary
                    extended_boundary.push(OrientedSimplex::from_simplex(
                        subspace_into_extended_subspace
                            .simplex_image(&boundary_simplex.clone().as_simplex())
                            .extend_by_point(point_extended_preimage.clone()),
                        &extended_interior_point,
                    ));
                }

                for interior_simplex in &self.interior {
                    // add it to the new boundary
                    extended_boundary.push(OrientedSimplex::from_simplex(
                        subspace_into_extended_subspace.simplex_image(interior_simplex),
                        &extended_interior_point,
                    ));

                    // add its cone to the new interior
                    extended_interior.push(
                        subspace_into_extended_subspace
                            .simplex_image(interior_simplex)
                            .extend_by_point(point_extended_preimage.clone()),
                    );
                }

                let ans = Self {
                    subspace: extended_subspace_into_ambient,
                    boundary: extended_boundary,
                    interior: extended_interior,
                };

                debug_assert!(ans.check().is_ok());
                ans
            }
        }
    }
}

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

// #[derive(Debug, Clone)]
// pub struct SimplicialComplex {
//     dim: usize,
//     simplices: Vec<Simplex>,
// }

// impl SimplicialComplex {
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
}
