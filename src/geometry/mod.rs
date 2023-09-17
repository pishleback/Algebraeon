#[allow(dead_code)]
use std::collections::{HashMap, HashSet};

use itertools::Itertools;
use malachite_q::Rational;

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

        if !are_points_nondegenerage(self.dim, &self.points) {
            return Err("Simplex is degenerate");
        }

        Ok(())
    }

    pub fn new(dim: usize, mut points: Vec<Point>) -> Self {
        points.sort();
        let ans = Self { dim, points };
        ans.check().unwrap();
        ans
    }

    pub fn try_new(dim: usize, mut points: Vec<Point>) -> Option<Self> {
        if are_points_nondegenerage(dim, &points) {
            points.sort();
            Some(Self { dim, points })
        } else {
            None
        }
    }

    pub fn n(&self) -> usize {
        self.points.len()
    }

    pub fn points(&self) -> Vec<Point> {
        self.points.clone()
    }

    fn bounding_box(&self) -> Box {
        Box::bounding_box(self.dim, &self.points)
    }

    fn transform(self, new_dim: usize, f: &dyn Fn(Point) -> Point) -> Self {
        let new_points = self
            .points
            .into_iter()
            .map(|p| p.transform(f))
            .collect_vec();
        for pt in new_points.iter() {
            debug_assert_eq!(pt.dim(), new_dim);
        }
        Self::new(new_dim, new_points)
    }

    pub fn has_vertex(&self, pt: &Point) -> bool {
        self.points.binary_search(pt).is_ok()
    }

    pub fn facet(&self, k: usize) -> Simplex {
        assert!(k <= self.n());
        let facet = Simplex {
            dim: self.dim,
            points: (0..self.n())
                .filter(|i| i != &k)
                .map(|i| self.points[i].clone())
                .collect(),
        };
        debug_assert!(facet.check().is_ok());
        facet
    }

    fn oriented_facet(&self, k: usize) -> OrientedSimplex {
        //oriented facets have positive side on the outside
        assert!(k <= self.n());
        let oriented_facet = OrientedSimplex::new(self.facet(k));
        match oriented_facet.sign_point(&self.points[k]) {
            std::cmp::Ordering::Less => oriented_facet,
            std::cmp::Ordering::Equal => panic!(),
            std::cmp::Ordering::Greater => oriented_facet.flipped(),
        }
    }

    pub fn facets(&self) -> Vec<Simplex> {
        (0..self.n()).map(|k| self.facet(k)).collect()
    }

    fn oriented_facets(&self) -> Vec<OrientedSimplex> {
        (0..self.n()).map(|k| self.oriented_facet(k)).collect()
    }

    pub fn ridges(&self) -> Vec<Simplex> {
        let mut ridges = vec![];
        for i in 0..self.n() {
            for j in i + 1..self.n() {
                let ridge = Simplex {
                    dim: self.dim,
                    points: (0..self.n())
                        .filter(|k| (k != &i) && (k != &j))
                        .map(|k| self.points[k].clone())
                        .collect(),
                };
                debug_assert!(ridge.check().is_ok());
                ridges.push(ridge);
            }
        }
        ridges
    }

    pub fn edges(&self) -> Vec<Simplex> {
        let mut edges = vec![];
        for i in 0..self.n() {
            for j in i + 1..self.n() {
                let edge = Simplex {
                    dim: self.dim,
                    points: vec![self.points[i].clone(), self.points[j].clone()],
                };
                debug_assert!(edge.check().is_ok());
                edges.push(edge);
            }
        }
        edges
    }

    // pub fn boundary_as_mesh(&self) -> Option<Mesh> {
    //     if self.n() == 0 {
    //         None
    //     } else {
    //         Some(Mesh {
    //             dim: self.dim,
    //             rank: self.n() - 1,
    //             simplices: HashMap::from([(
    //                 self.clone(),
    //                 (0..self.n())
    //                     .map(|_facet_idx| MeshAdjacency::Boundary)
    //                     .collect(),
    //             )]),
    //         })
    //     }
    // }

    // pub fn as_convex_hull(&self) -> ConvexHull {
    //     ConvexHull {
    //         dim: self.dim,
    //         points: self.points.clone(),
    //         facets: self.oriented_facets().into_iter().collect(),
    //         ridges: self.ridges().into_iter().collect(),
    //     }
    // }

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

    pub fn affine_span(&self) -> AffineLattice<Rational> {
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

    pub fn middle(&self) -> Point {
        let mut coords = (0..self.dim).map(|_i| Rational::from(0)).collect_vec();

        for pt in &self.points {
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

    // pub fn measure(&self) -> Rational {
    //not gooood need to check that its full rank
    //     if self.n() == 0 {
    //         Rational::from(0)
    //     } else {
    //         let mat = Matrix::construct(self.n() - 1, self.n() - 1, |r, c| {
    //             (&self.points[c + 1] - &self.points[0]).get_coord(r)
    //         });
    //         let mut vol = QQ_MAT.det(mat).unwrap();
    //         vol.abs_assign();
    //         vol = vol.div(Rational::from(factorial(Natural::from(self.n() - 1))));
    //         vol
    //     }
    // }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct OrientedSimplex {
    simplex: Simplex,
    flip: bool,
}

impl OrientedSimplex {
    pub fn check(&self) -> Result<(), &'static str> {
        self.simplex.check()?;

        if self.simplex.dim == 0 {
            return Err("Can't have a hyperplane in zero dimensional space");
        }

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

    fn from_points(dim: usize, points: Vec<Point>, neg_pt: &Point) -> Self {
        for point in &points {
            debug_assert_eq!(dim, point.dim());
        }
        debug_assert_eq!(dim, neg_pt.dim());
        debug_assert_eq!(dim, points.len());

        let ans = Self {
            simplex: Simplex::new(dim, points),
            flip: false,
        };

        match ans.sign_point(&neg_pt) {
            std::cmp::Ordering::Less => ans,
            std::cmp::Ordering::Equal => panic!(),
            std::cmp::Ordering::Greater => ans.flipped(),
        }
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
        self.det_vector(&(point - &self.simplex.points[0]))
    }

    pub fn det_vector(&self, vec: &Vector) -> Rational {
        //vector relative to self.root
        let mat = Matrix::construct(self.simplex.dim, self.simplex.dim, |r, c| {
            if c == self.simplex.dim - 1 {
                vec.get_coord(r)
            } else {
                (&self.simplex.points[c + 1] - &self.simplex.points[0]).get_coord(r)
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

pub fn shape_union(dim: usize, shapes: Vec<Shape>) -> Shape {
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

    fn transform(self, new_dim: usize, f: &dyn Fn(Point) -> Point) -> Self {
        Self::new(
            new_dim,
            self.simplices
                .into_iter()
                .map(|s| s.transform(new_dim, f))
                .collect(),
        )
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
        return self.clone();

        //connected decomposition (TODO!)
        //convex decomposition for each connected component

        let dim = self.dim;
        let mut points = HashSet::new();
        for simplex in &self.simplices {
            for point in &simplex.points {
                points.insert(point);
            }
        }
        let hull = convexhull(dim, points.clone().into_iter().map(|p| p.clone()).collect());
        let remaining = hull.subtract_nosimp(&self);
        if remaining.is_empty() {
            todo!("we are convex");
        } else {
            let mut canvas = super::drawing::canvas2d::Shape2dCanvas::new();
            canvas.draw_shape(self, (1.0, 1.0, 0.0));
            canvas.draw_shape(&hull.subtract_nosimp(&self), (0.0, 1.0, 0.0));
            crate::drawing::Canvas::run(canvas);
        }

        return self.clone();

        let dim = self.dim;
        let mut ans = self.clone();

        'main_loop: loop {
            for simplex in &ans.simplices {
                println!("pp");

                let mut local_simplices = vec![simplex.clone()];
                let mut other_simplices: HashSet<&Simplex> = ans.simplices.iter().collect();
                other_simplices.remove(simplex);

                let mut done = false;
                while !done {
                    done = true;
                    let local_convexhull = {
                        let mut points = HashSet::new();
                        for s in &local_simplices {
                            for p in s.points() {
                                points.insert(p.clone());
                            }
                        }
                        convexhull(dim, points.into_iter().collect())
                    };

                    for s in other_simplices.clone() {
                        let (inner, outer) = cut_shape_by_simplex(s, &local_convexhull);
                        if !inner.is_empty() {
                            other_simplices.remove(s);
                            local_simplices.push(s.clone());
                            done = false;
                        }
                    }

                    //now local shape is a subshape of ans
                    //and every convex subshape appears as some local shape in this loop
                    //note that a particular local shape may not be convex
                    let local_shape =
                        Shape::new(dim, local_simplices.iter().map(|s| (*s).clone()).collect());

                    if local_convexhull.simplices.len() < local_shape.simplices.len() {
                        if local_shape.equal(&local_convexhull) {
                            ans = shape_union(
                                dim,
                                vec![
                                    Shape::new(
                                        dim,
                                        other_simplices.iter().map(|s| (*s).clone()).collect(),
                                    ),
                                    local_convexhull,
                                ],
                            );
                            println!("simp {}", ans.simplices.len());
                            continue 'main_loop;
                        }
                    }
                }

                // let mut canvas = super::drawing::canvas2d::Shape2dCanvas::new();
                // canvas.draw_shape(&ans, (1.0, 1.0, 0.0));
                // canvas.draw_shape(&local_shape, (0.0, 1.0, 0.0));
                // canvas.draw_shape(
                //     &local_convexhull.symmetric_difference_nosimp(&local_shape),
                //     (0.0, 0.0, 1.0),
                // );
                // crate::drawing::Canvas::run(canvas);
            }
            break;
        }
        ans
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
        shape_union(dim, parts)
    }

    fn union_nosimp(&self, other: &Self) -> Self {
        let dim = self.dim;
        assert_eq!(dim, other.dim);
        shape_union(
            dim,
            vec![
                self.intersect_nosimp(other),
                self.symmetric_difference_nosimp(other),
            ],
        )
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
        shape_union(
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
    //n=1: shape is the shell of a line
    //n=2: shape is the shell of a polygon
    //n=3: shape is the shell of a solid
    //...

    //shape should be the full boundary of a convex polytope

    if shape.simplices.is_empty() {
        return Shape::empty(shape.dim);
    }

    let n = shape.simplices.iter().map(|s| s.n()).max().unwrap();
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
                .map(|s| s.affine_span())
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

fn cut_simplex_by_plane(cut_plane: &OrientedSimplex, simplex: &Simplex) -> (Shape, Shape, Shape) {
    let dim = cut_plane.dim();
    assert_eq!(dim, simplex.dim);
    debug_assert!(simplex.n() <= dim + 1);

    match simplex.n() {
        0 => (Shape::empty(dim), Shape::empty(dim), Shape::empty(dim)),
        1 => match cut_plane.sign_point(&simplex.points[0]) {
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
        left = shape_union(dim, vec![left, s_left]);
        middle = shape_union(dim, vec![middle, s_middle]);
        right = shape_union(dim, vec![right, s_right]);
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
                outside = shape_union(dim, vec![inside_left, inside_right, outside]);
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
                    flat_outside = shape_union(
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
                    shape_union(dim, vec![unflat_outside, outside]),
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
        inside = shape_union(dim, vec![inside, s_inside]);
        outside = shape_union(dim, vec![outside, s_outside]);
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
            let middle_point = starting_simplex.middle();
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
                    let mut new_facet_points = h.points;
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
    shape_union(dim, vec![shell, interior])
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
        let d = shape_union(3, vec![a, b, c]);
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
