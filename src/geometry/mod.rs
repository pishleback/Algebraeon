use std::collections::{HashMap, HashSet};

use itertools::Itertools;
use malachite_q::Rational;
use rayon::prelude::IndexedParallelIterator;

use super::rings::lattice::*;
use super::rings::matrix::*;
use crate::rings::nzq::QQ;

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

    pub fn as_matrix(&self) -> Matrix<Rational> {
        Matrix::construct(self.coords.len(), 1, |r, c| self.coords[r].clone())
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

    pub fn as_matrix(&self) -> Matrix<Rational> {
        Matrix::construct(self.coords.len(), 1, |r, c| self.coords[r].clone())
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

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Simplex {
    dim: usize,
    points: Vec<Point>, //ordered
}

impl Simplex {
    pub fn check(&self) -> Result<(), &'static str> {
        let mut sorted_points = self.points.clone();
        sorted_points.sort();
        if self.points != sorted_points {
            return Err("Simplex points are not sorted");
        }

        for p in &self.points {
            if p.dim() != self.dim {
                return Err("Simplex should live in the same dimension as its points");
            }
        }

        if self.points.len() >= 1 {
            let root = self.points[0].clone();
            let mut vecs = vec![];
            for i in 1..self.points.len() {
                vecs.push(&self.points[i] - &root);
            }
            let mat = Matrix::construct(self.dim, vecs.len(), |r, c| vecs[c].get_coord(r));
            if QQ_MAT.rank(mat) != vecs.len() {
                return Err("Simplex is degenerate");
            }
        }

        Ok(())
    }

    pub fn new(dim: usize, mut points: Vec<Point>) -> Self {
        points.sort();
        let ans = Self { dim, points };
        ans.check().unwrap();
        ans
    }

    pub fn n(&self) -> usize {
        self.points.len()
    }

    pub fn has_vertex(&self, pt: &Point) -> bool {
        self.points.binary_search(pt).is_ok()
    }

    pub fn facets(&self) -> Vec<Simplex> {
        let mut facets = vec![];
        for i in 0..self.n() {
            let facet = Simplex {
                dim: self.dim,
                points: (0..self.n())
                    .filter(|j| &i != j)
                    .map(|j| self.points[j].clone())
                    .collect(),
            };
            debug_assert!(facet.check().is_ok());
            facets.push(facet);
        }
        facets
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
                points: subset.into_iter().map(|i| self.points[i].clone()).collect(),
            };
            debug_assert!(b.check().is_ok());
            parts.push(b);
        }
        Shape {
            dim: self.dim,
            simplices: parts,
        }
    }

    pub fn span(&self) -> AffineLattice<Rational> {
        if self.points.len() == 0 {
            QQ_AFFLAT.empty(self.dim, 1)
        } else {
            QQ_AFFLAT.from_offset_and_linear_lattice(
                self.dim,
                1,
                self.points[0].as_matrix(),
                QQ_LINLAT.from_basis(
                    self.dim,
                    1,
                    (1..self.points.len())
                        .map(|i| (&self.points[i] - &self.points[0]).as_matrix())
                        .collect(),
                ),
            )
        }
    }

    // pub fn oriented_hyperplanes(&self) -> Vec<OrientedHyperplane> {
    //     //one OrientedHyperplane for each facet such that
    //     //the intersection of the positive sides is the interior of the simplex
    //     todo!()
    // }
}

#[derive(Debug, Clone)]
pub struct OrientedHyperplane {
    dim: usize,
    root: Point,
    vecs: Vec<Point>,
}

impl OrientedHyperplane {
    pub fn check(&self) -> Result<(), &'static str> {
        if self.dim == 0 {
            return Err("Can't have a hyperplane in zero dimensional space");
        }

        if self.root.dim() != self.dim {
            return Err("Hyperplane should live in the same dimension as its root point");
        }

        for p in &self.vecs {
            if p.dim() != self.dim {
                return Err("Hyperplane should live in the same dimension as its vecs");
            }
        }

        if self.vecs.len() != self.dim - 1 {
            return Err("Hyperplane should have dimension one less than the space it lives in");
        }

        let dim = self.vecs.iter().map(|v| v.dim()).max().unwrap();
        let mat = Matrix::construct(dim, self.vecs.len(), |r, c| self.vecs[c].get_coord(r));
        if QQ_MAT.rank(mat) != self.vecs.len() {
            return Err("Hyperplane is degenerate");
        }

        Ok(())
    }

    pub fn new(dim: usize, root: Point, vecs: Vec<Point>) -> Self {
        debug_assert_ne!(dim, 0);
        let ans = Self { dim, root, vecs };
        ans.check().unwrap();
        ans
    }

    fn det_point(&self, point: &Point) -> Rational {
        self.det_vector(&(point - &self.root))
    }

    fn det_vector(&self, vec: &Vector) -> Rational {
        //vector relative to self.root
        assert_eq!(self.vecs.len() + 1, self.dim);

        let mat = Matrix::construct(self.dim, self.dim, |r, c| {
            if c == self.dim - 1 {
                vec.get_coord(r)
            } else {
                self.vecs[c].get_coord(r)
            }
        });
        QQ_MAT.det(mat).unwrap()
    }

    fn sign_point(&self, point: &Point) -> std::cmp::Ordering {
        self.det_point(point).cmp(&Rational::from(0))
    }

    pub fn rank(&self) -> usize {
        self.vecs.len()
    }
}

//disjoint union of simplices
#[derive(Debug, Clone)]
pub struct Shape {
    dim: usize,
    simplices: Vec<Simplex>,
}

fn shape_union(dim: usize, shapes: Vec<Shape>) -> Shape {
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
            if simplex.dim != self.dim {
                return Err("Simplex dim does not match shape dim");
            }
            if simplex.n() == 0 {
                return Err("Shape contains empty simplex");
            }
            simplex.check()?;
        }
        //TODO: check that the simplices are disjoint
        Ok(())
    }

    pub fn empty(dim: usize) -> Self {
        Self {
            dim,
            simplices: vec![],
        }
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

    pub fn skeleton(&self, n: usize) -> Self {
        Shape {
            dim: self.dim,
            simplices: self
                .simplices
                .iter()
                .filter(|s| s.n() == n)
                .map(|s| s.clone())
                .collect(),
        }
    }
}

fn interior_of_convex_hollow_shell(n: usize, shape: &Shape) -> Shape {
    debug_assert!(shape.is_complete());
    //n=1: shape is the shell of a line
    //n=2: shape is the shell of a polygon
    //n=3: shape is the shell of a solid
    //...

    //shape should be the full boundary of a convex polytope
    debug_assert!(n >= 1);
    for s in &shape.simplices {
        debug_assert!(s.n() <= n);
    }

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
                    shape.simplices[0].points[0].clone(),
                    shape.simplices[1].points[0].clone(),
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
            let root = &shape.simplices[0].points[0];
            // println!("root = {:?}", root);
            // for s in &shape.simplices {
            //     println!("shell: {:?}", s);
            // }

            //2) find all cells adjacent to root
            let adj_cell_spans: Vec<_> = shape
                .simplices
                .iter()
                .filter(|s| s.has_vertex(root))
                .map(|s| s.span())
                .collect();

            //3) find which adjacent faces are vertex belongs to
            let mut point_degen = HashMap::new();
            for s in &shape.simplices {
                if s.n() == 1 {
                    let pt = &s.points[0];
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
                for p in &s.points {
                    debug_assert!(point_degen.contains_key(p));
                }
            }
            // println!("{:?}", point_degen);

            //4) for each shell simplex, if not all vertices lie in some adjacent face span, fill it in
            let mut interior = Shape::empty(shape.dim);
            for s in &shape.simplices {
                debug_assert!(s.points.len() >= 1);
                let mut common = point_degen[&s.points[0]].clone();
                common = common
                    .into_iter()
                    .filter(|idx| {
                        (1..s.points.len())
                            .map(|i| &s.points[i])
                            .all(|pt| point_degen[pt].contains(idx))
                    })
                    .collect();
                if common.len() == 0 {
                    let filler = Simplex::new(shape.dim, {
                        let mut filler_pts = vec![root.clone()];
                        filler_pts.append(&mut s.points.clone());
                        // println!("filler_pts = {:?}", filler_pts);
                        filler_pts
                    });
                    // println!("filler simplex: {:?}", filler);
                    interior = shape_union(shape.dim, vec![interior, Shape::simplex(filler)]);
                }
            }
            interior
        }
    }
}

pub fn cut_simplex(plane: &OrientedHyperplane, simplex: &Simplex) -> (Shape, Shape, Shape) {
    let dim = plane.dim;
    assert_eq!(dim, simplex.dim);
    assert_eq!(dim, plane.rank() + 1);
    debug_assert!(simplex.n() <= dim + 1);

    match simplex.n() {
        0 => (Shape::empty(dim), Shape::empty(dim), Shape::empty(dim)),
        1 => match plane.sign_point(&simplex.points[0]) {
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
            let (p, q) = (&simplex.points[0], &simplex.points[1]);
            match (plane.sign_point(p), plane.sign_point(q)) {
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
                    let t = -plane.det_point(&p) / plane.det_vector(&(q - p));
                    let r = p + &(&t * &(q - p));
                    debug_assert!(0 < t && t < 1);
                    debug_assert_eq!(plane.det_point(&r), Rational::from(0));
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
            let (open_left, hollow_middle, open_right) = cut_shape(&plane, &simplex.boundary());

            let interior_middle = interior_of_convex_hollow_shell(n - 2, &hollow_middle);
            let hollow_left = shape_union(
                dim,
                vec![
                    open_left.clone(),
                    hollow_middle.clone(),
                    interior_middle.clone(),
                ],
            );
            let hollow_right = shape_union(
                dim,
                vec![
                    open_right.clone(),
                    hollow_middle.clone(),
                    interior_middle.clone(),
                ],
            );
            let interior_left = interior_of_convex_hollow_shell(n - 1, &hollow_left);
            let interior_right = interior_of_convex_hollow_shell(n - 1, &hollow_right);
            (interior_left, interior_middle, interior_right)
        }
    }
}

pub fn cut_shape(plane: &OrientedHyperplane, shape: &Shape) -> (Shape, Shape, Shape) {
    let dim = plane.dim;
    assert_eq!(dim, shape.dim);
    assert_eq!(dim, plane.rank() + 1);

    let (mut left, mut middle, mut right) =
        (Shape::empty(dim), Shape::empty(dim), Shape::empty(dim));
    for simplex in &shape.simplices {
        let (s_left, s_middle, s_right) = cut_simplex(&plane, &simplex);
        left = shape_union(dim, vec![left, s_left]);
        middle = shape_union(dim, vec![middle, s_middle]);
        right = shape_union(dim, vec![right, s_right]);
    }
    (left, middle, right)
}

pub fn quickhull_boundary(dim: usize, points: Vec<Point>) -> Option<Shape> {
    //implement https://en.wikipedia.org/wiki/Quickhull
    //return None if the convex hull is flat
    for point in &points {
        assert_eq!(point.dim(), dim);
    }
    todo!()
}

pub fn quickhull_interior(dim: usize, points: Vec<Point>) -> Shape {
    for point in &points {
        assert_eq!(point.dim(), dim);
    }
    //use quickhull_boundary and find its interior with interior_of_convex_hollow_shell
    //if quickhull_boundary is None, then the interior is empty
    todo!()
}

pub fn quickhull_complete(dim: usize, points: Vec<Point>) -> Shape {
    for point in &points {
        assert_eq!(point.dim(), dim);
    }
    //use quickhull_boundary and union with its interior found using interior_of_convex_hollow_shell
    //if quickhull_boundary is None i.e. the convex hull is flat, then try again but in a smaller dimension untill it works
    todo!()
}

#[cfg(test)]
mod geometry_tests {
    use std::str::FromStr;

    use super::*;

    #[test]
    fn test_something() {}
}
