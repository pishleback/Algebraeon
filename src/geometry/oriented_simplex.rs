use malachite_q::Rational;

use crate::rings::matrix::Matrix;

use super::{
    simplex::Simplex,
    vector::Vector,
};

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct OrientedSimplex {
    simplex: Simplex, //every simplex has a natural orientation induced by the ordering of the points of the space
    flip: bool,       //whether the orientation is the natural one or its opposite
}

impl OrientedSimplex {
    pub fn check(&self) -> Result<(), &'static str> {
        self.simplex.check()?;

        if self.simplex.n() != self.simplex.dim() {
            return Err(
                "OrientedSimplex should have dimension one less than the space it lives in",
            );
        }

        Ok(())
    }

    pub fn new(simplex: Simplex) -> Self {
        assert_eq!(simplex.dim(), simplex.n());
        let ans = Self {
            simplex,
            flip: false,
        };
        ans.check().unwrap();
        ans
    }

    pub fn as_simplex(self) -> Simplex {
        self.simplex
    }

    pub fn from_simplex(simplex: Simplex, neg_pt: &Vector) -> Self {
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

    fn from_points(dim: usize, points: Vec<Vector>, neg_pt: &Vector) -> Self {
        for point in &points {
            debug_assert_eq!(dim, point.dim());
        }
        debug_assert_eq!(dim, neg_pt.dim());
        debug_assert_eq!(dim, points.len());

        Self::from_simplex(Simplex::new(dim, points), neg_pt)
    }

    pub fn dim(&self) -> usize {
        self.simplex.dim()
    }

    pub fn flipped(self) -> Self {
        Self {
            simplex: self.simplex,
            flip: !self.flip,
        }
    }

    pub fn det_point(&self, point: &Vector) -> Rational {
        debug_assert_eq!(point.dim(), self.dim());
        if self.dim() == 0 {
            Rational::from(0)
        } else {
            self.det_vector(&(point - self.simplex.points().iter().next().unwrap()))
        }
    }

    pub fn det_vector(&self, vec: &Vector) -> Rational {
        //vector relative to self.root
        let mat = Matrix::construct(self.dim(), self.dim(), |r, c| {
            if c == self.dim() - 1 {
                vec.get_coord(r)
            } else {
                (self.simplex.points().iter().nth(c + 1).unwrap()
                    - self.simplex.points().iter().next().unwrap())
                .get_coord(r)
            }
        });
        let d = mat.det().unwrap();
        match self.flip {
            false => d,
            true => -d,
        }
    }

    pub fn sign_point(&self, point: &Vector) -> std::cmp::Ordering {
        self.det_point(point).cmp(&Rational::from(0))
    }

    // pub fn rank(&self) -> usize {
    //     self.vecs.len()
    // }
}
