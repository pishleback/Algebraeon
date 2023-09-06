use std::collections::{HashMap, HashSet};

use itertools::Itertools;
use malachite_q::Rational;
use rayon::prelude::IndexedParallelIterator;

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
}

#[derive(Debug, Clone)]
pub struct Simplex {
    dim: usize,
    points: Vec<Point>,
}

impl Simplex {
    pub fn check(&self) -> Result<(), &'static str> {
        use super::rings::matrix::*;

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

    pub fn new(dim: usize, points: Vec<Point>) -> Self {
        let ans = Self { dim, points };
        ans.check().unwrap();
        ans
    }

    pub fn n(&self) -> usize {
        self.points.len()
    }

    pub fn faces(&self) -> Vec<Simplex> {
        let mut faces = vec![];
        for i in 0..self.n() {
            faces.push(Simplex {
                dim: self.dim,
                points: (0..self.n())
                    .filter(|j| &i != j)
                    .map(|j| self.points[j].clone())
                    .collect(),
            });
        }
        faces
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
            parts.push(Simplex {
                dim: self.dim,
                points: subset.into_iter().map(|i| self.points[i].clone()).collect(),
            });
        }
        Shape {
            dim: self.dim,
            simplices: parts,
        }
    }
}

#[derive(Debug, Clone)]
pub struct OrientedHyperplane {
    dim: usize,
    root: Point,
    vecs: Vec<Point>,
}

impl OrientedHyperplane {
    pub fn check(&self) -> Result<(), &'static str> {
        use super::rings::matrix::*;

        if self.root.dim() != self.dim {
            return Err("Hyperplane should live in the same dimension as its root point");
        }
        for v in &self.vecs {
            if v.dim() != self.dim {
                return Err("Hyperplane should live in the same dimension as its vecs");
            }
        }

        if self.vecs.len() >= 1 {
            let dim = self.vecs.iter().map(|v| v.dim()).max().unwrap();
            let mat = Matrix::construct(dim, self.vecs.len(), |r, c| self.vecs[c].get_coord(r));
            if QQ_MAT.rank(mat) != self.vecs.len() {
                return Err("Hyperplane is degenerate");
            }
        }

        Ok(())
    }

    pub fn new(dim: usize, root: Point, vecs: Vec<Point>) -> Self {
        let ans = Self { dim, root, vecs };
        ans.check().unwrap();
        ans
    }

    fn det_point(&self, point: &Point) -> Rational {
        self.det_vector(&(point - &self.root))
    }

    fn det_vector(&self, vec: &Vector) -> Rational {
        //vector relative to self.root
        use super::rings::matrix::*;
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
                            Shape::simplex(Simplex {
                                dim,
                                points: vec![p.clone(), r.clone()],
                            }),
                            Shape::simplex(Simplex {
                                dim,
                                points: vec![r.clone()],
                            }),
                            Shape::simplex(Simplex {
                                dim,
                                points: vec![r, q.clone()],
                            }),
                        ),
                        (std::cmp::Ordering::Less, std::cmp::Ordering::Greater) => (
                            Shape::simplex(Simplex {
                                dim,
                                points: vec![q.clone(), r.clone()],
                            }),
                            Shape::simplex(Simplex {
                                dim,
                                points: vec![r.clone()],
                            }),
                            Shape::simplex(Simplex {
                                dim,
                                points: vec![r, p.clone()],
                            }),
                        ),
                        _ => panic!(),
                    }
                }
            }
        }
        n => {
            fn interior_of_hollow_shell(n: usize, shape: &Shape) -> Shape {
                //shape should be the full boundary of a convex polytope
                debug_assert!(n >= 1);
                println!("interior dim={:?} shape={:?}", n, shape);

                // {
                //     //sanity check Euler characteristic
                //     let mut euler: isize = 0;
                //     for s in &shape.simplices {
                //         println!("{:?} {:?}", n, s);
                //         if s.n() % 2 == 0 {
                //             euler -= 1;
                //         } else {
                //             euler += 1;
                //         }
                //     }
                //     if n % 2 == 0 {
                //         debug_assert_eq!(euler, 0);
                //     } else {
                //         debug_assert_eq!(euler, 2);
                //     }
                // }

                if n == 1 {
                    for s in &shape.simplices {
                        debug_assert_eq!(s.n(), 1);
                    }

                    //The only way this happens is when cutting a triangle
                    //Thus there are only 3 possibilites: shape is 0, 1, or 2 points
                    match shape.simplices.len() {
                        0 => Shape::empty(shape.dim), //no points has empty interior
                        1 => Shape::empty(shape.dim), //single point has empty interior
                        2 => Shape::simplex(Simplex {
                            dim: shape.dim,
                            points: vec![
                                shape.simplices[0].points[0].clone(),
                                shape.simplices[1].points[0].clone(),
                            ],
                        }), //pair of points has interior of a line
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

                    for s in &shape.simplices {
                        println!("{:?}", s);
                    }

                    println!("todo here");

                    todo!()
                }
            }

            debug_assert!(n >= 3);
            let (open_left, hollow_middle, open_right) = cut_shape(&plane, &simplex.boundary());

            let interior_middle = interior_of_hollow_shell(n - 2, &hollow_middle);
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
                    open_left.clone(),
                    hollow_middle.clone(),
                    interior_middle.clone(),
                ],
            );
            let interior_left = interior_of_hollow_shell(n - 1, &hollow_left);
            let interior_right = interior_of_hollow_shell(n - 1, &hollow_right);
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

#[cfg(test)]
mod geometry_tests {
    use std::str::FromStr;

    use super::*;

    #[test]
    fn test_something() {}
}
