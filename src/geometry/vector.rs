use malachite_q::Rational;

use crate::rings::matrix::{Matrix, QQ_MAT};

//represent a point in space
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Point {
    coords: Vec<Rational>,
}

//represent a vector in space
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
    pub fn transform(self, f: &dyn Fn(Point) -> Point) -> Self {
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

pub fn are_points_nondegenerage(dim: usize, points: Vec<&Point>) -> bool {
    for point in &points {
        debug_assert_eq!(dim, point.dim());
    }
    if points.len() >= 1 {
        let root = points[0].clone();
        let mut vecs = vec![];
        for i in 1..points.len() {
            vecs.push(points[i] - &root);
        }
        let mat = Matrix::construct(dim, vecs.len(), |r, c| vecs[c].get_coord(r));
        if QQ_MAT.rank(mat) != vecs.len() {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use super::*;

    #[test]
    fn test_arithmetic() {
        let a = Point::new(vec![
            Rational::from(0),
            Rational::from(2),
            Rational::from(3),
        ]);
        let b = Point::new(vec![
            Rational::from(2),
            Rational::from(-4),
            Rational::from(0),
        ]);
        let c = Point::new(vec![
            Rational::from(2),
            Rational::from(-2),
            Rational::from(3),
        ]);

        let u = Vector::new(vec![
            Rational::from(0),
            Rational::from(2),
            Rational::from(3),
        ]);
        let v = Vector::new(vec![
            Rational::from(2),
            Rational::from(-4),
            Rational::from(0),
        ]);
        let w = Vector::new(vec![
            Rational::from(2),
            Rational::from(-2),
            Rational::from(3),
        ]);

        assert_eq!(&u + &v, w);
        assert_eq!(&a + &v, c);
    }

    fn test_degenerate() {
        let a = Point::new(vec![
            Rational::from(0),
            Rational::from(2),
            Rational::from(3),
        ]);
        let b = Point::new(vec![
            Rational::from(2),
            Rational::from(-4),
            Rational::from(0),
        ]);
        let c = Point::new(vec![
            Rational::from(2),
            Rational::from(-2),
            Rational::from(3),
        ]);

        assert!(are_points_nondegenerage(3, vec![&a, &b, &c]));

        let a = Point::new(vec![
            Rational::from(0),
            Rational::from(2),
            Rational::from(3),
        ]);
        let b = Point::new(vec![
            Rational::from(0),
            Rational::from(-2),
            Rational::from(-3),
        ]);
        let c = Point::new(vec![
            Rational::from(0),
            Rational::from(4),
            Rational::from(6),
        ]);

        assert!(!are_points_nondegenerage(3, vec![&a, &b, &c]));
    }
}
