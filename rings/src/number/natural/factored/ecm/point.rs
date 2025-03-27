use algebraeon_nzq::{integer::Integer, natural::Natural};
use algebraeon_sets::structure::*;

use crate::structure::{quotient::QuotientStructure, structure::*};

/// Montgomery form of Points in an elliptic curve.
///
/// In this form, the addition and doubling of points
/// does not need any y-coordinate information thus
/// decreasing the number of operations.
/// Using Montgomery form we try to perform point addition
/// and doubling in least amount of multiplications.
///
/// The elliptic curve used here is of the form
/// `(E : b*y**2*z = x**3 + a*x**2*z + x*z**2)`.
/// The `a_24` parameter is equal to `(a + 2)/4`.
///
/// References
/// ----------
/// - http://www.hyperelliptic.org/tanja/SHARCS/talks06/Gaj.pdf
#[derive(Debug, Clone)]
pub struct Point {
    /// X coordinate of the Point
    pub x_cord: Integer,
    /// Z coordinate of the Point
    pub z_cord: Integer,
    /// Parameter of the elliptic curve in Montgomery form
    pub a_24: Integer,
    /// modulus
    pub modulus: QuotientStructure<CannonicalStructure<Integer>, false>,
}

impl Default for Point {
    fn default() -> Self {
        Self {
            x_cord: Default::default(),
            z_cord: Default::default(),
            a_24: Default::default(),
            modulus: QuotientStructure::new_ring(Integer::structure(), Integer::one()),
        }
    }
}

impl Point {
    /// Initial parameters for the Point struct.
    ///
    /// # Parameters
    ///
    /// - `x_cord`: X coordinate of the Point
    /// - `z_cord`: Z coordinate of the Point
    /// - `a_24`: Parameter of the elliptic curve in Montgomery form
    /// - `mod`: modulus
    pub fn new(
        x_cord: Integer,
        z_cord: Integer,
        a_24: Integer,
        modulus: QuotientStructure<CannonicalStructure<Integer>, false>,
    ) -> Point {
        Point {
            x_cord,
            z_cord,
            a_24,
            modulus,
        }
    }

    /// Adds two points `self` and `Q` where `diff = self - Q`.
    ///
    /// This algorithm requires 6 multiplications. The assumption is that `self.x_cord * Q.x_cord * (self.x_cord - Q.x_cord) != 0`.
    /// Using this algorithm speeds up the addition by reducing the number of multiplications required.
    ///
    /// The `mont_ladder` algorithm is constructed in a way that the difference between intermediate points is always equal to the initial point.
    /// So, we always know what the difference between the point is.
    ///
    /// # Parameters
    ///
    /// - `Q`: Point on the curve in Montgomery form.
    /// - `diff`: `self - Q`
    pub fn add(&self, q: &Point, diff: &Point) -> Point {
        let u = Integer::from(&self.x_cord - &self.z_cord) * Integer::from(&q.x_cord + &q.z_cord);
        let v = Integer::from(&self.x_cord + &self.z_cord) * Integer::from(&q.x_cord - &q.z_cord);
        let add = Integer::from(&u + &v);
        let subt = u - v;
        let x_cord = self.modulus.reduce(&diff.z_cord * &add * &add);
        let z_cord = self.modulus.reduce(&diff.x_cord * &subt * &subt);

        Point::new(x_cord, z_cord, self.a_24.clone(), self.modulus.clone())
    }

    /// Doubles a point in an elliptic curve in Montgomery form.
    pub fn double(&self) -> Point {
        let u =
            Integer::from(&self.x_cord + &self.z_cord) * Integer::from(&self.x_cord + &self.z_cord);
        let v =
            Integer::from(&self.x_cord - &self.z_cord) * Integer::from(&self.x_cord - &self.z_cord);
        let diff = Integer::from(&u - &v);
        let x_cord = self.modulus.reduce(u * &v);
        let z_cord = self.modulus.reduce((v + &self.a_24 * &diff) * diff);

        Point::new(x_cord, z_cord, self.a_24.clone(), self.modulus.clone())
    }

    /// Scalar multiplication of a point in Montgomery form
    /// using Montgomery Ladder Algorithm.
    /// A total of 11 multiplications are required in each step of this
    /// algorithm.
    ///
    /// # Parameters
    ///
    /// - `k`: The positive integer multiplier
    pub fn mont_ladder(&self, k: &Natural) -> Point {
        let mut q = self.clone();
        let mut r = self.double();
        for i in k.bits().rev().skip(1) {
            if i {
                q = r.add(&q, self);
                r = r.double();
            } else {
                r = q.add(&r, self);
                q = q.double();
            }
        }
        q
    }
}

impl PartialEq for Point {
    /// Two points are equal if X/Z of both points are equal.
    fn eq(&self, other: &Self) -> bool {
        if self.a_24 != other.a_24 {
            return false;
        }
        let modulus = &self.modulus;
        if modulus != &other.modulus {
            return false;
        }
        modulus.equal(
            &(modulus.inv(&self.z_cord).unwrap() * &self.x_cord),
            &(modulus.inv(&other.z_cord).unwrap() * &other.x_cord),
        )
    }
}

#[cfg(test)]
mod tests {
    use algebraeon_sets::structure::MetaType;

    use super::*;

    #[test]
    fn test_point_add() {
        let mod_29 = QuotientStructure::new_ring(Integer::structure(), Integer::from(29));

        let p1 = Point::new(11.into(), 16.into(), 7.into(), mod_29.clone());
        let p2 = Point::new(13.into(), 10.into(), 7.into(), mod_29.clone());
        let p3 = p2.add(&p1, &p1);

        assert_eq!(p3.x_cord, Integer::from(23));
        assert_eq!(p3.z_cord, Integer::from(17));
    }

    #[test]
    fn test_point_double() {
        let mod_29 = QuotientStructure::new_ring(Integer::structure(), Integer::from(29));

        let p1 = Point::new(11.into(), 16.into(), 7.into(), mod_29.clone());
        let p2 = p1.double();

        assert_eq!(p2.x_cord, Integer::from(13));
        assert_eq!(p2.z_cord, Integer::from(10));
    }

    #[test]
    fn test_point_mont_ladder() {
        let mod_29 = QuotientStructure::new_ring(Integer::structure(), Integer::from(29));

        let p1 = Point::new(11.into(), 16.into(), 7.into(), mod_29.clone());
        let p3 = p1.mont_ladder(&3u32.into());

        assert_eq!(p3.x_cord, Integer::from(23));
        assert_eq!(p3.z_cord, Integer::from(17));
    }

    #[test]
    fn test_point() {
        let mod_101 = QuotientStructure::new_ring(Integer::structure(), Integer::from(101));

        let a: Integer = 10.into();
        let a_24: Integer = (a + Integer::from(2)) * mod_101.inv(&Integer::from(4)).unwrap();

        let p1 = Point::new(10.into(), 17.into(), a_24.clone(), mod_101.clone());
        let p2 = p1.double();
        assert_eq!(
            p2,
            Point::new(68.into(), 56.into(), a_24.clone(), mod_101.clone())
        );
        let p4 = p2.double();
        assert_eq!(
            p4,
            Point::new(22.into(), 64.into(), a_24.clone(), mod_101.clone())
        );
        let p8 = p4.double();
        assert_eq!(
            p8,
            Point::new(71.into(), 95.into(), a_24.clone(), mod_101.clone())
        );
        let p16 = p8.double();
        assert_eq!(
            p16,
            Point::new(5.into(), 16.into(), a_24.clone(), mod_101.clone())
        );
        let p32 = p16.double();
        assert_eq!(
            p32,
            Point::new(33.into(), 96.into(), a_24.clone(), mod_101.clone())
        );

        // p3 = p2 + p1
        let p3 = p2.add(&p1, &p1);
        assert_eq!(
            p3,
            Point::new(1.into(), 61.into(), a_24.clone(), mod_101.clone())
        );
        // p5 = p3 + p2 or p4 + p1
        let p5 = p3.add(&p2, &p1);
        assert_eq!(
            p5,
            Point::new(49.into(), 90.into(), a_24.clone(), mod_101.clone())
        );
        assert_eq!(p5, p4.add(&p1, &p3));
        // # p6 = 2*p3
        let p6 = p3.double();
        assert_eq!(
            p6,
            Point::new(87.into(), 43.into(), a_24.clone(), mod_101.clone())
        );
        assert_eq!(p6, p4.add(&p2, &p2));
        // # p7 = p5 + p2
        let p7 = p5.add(&p2, &p3);
        assert_eq!(
            p7,
            Point::new(69.into(), 23.into(), a_24.clone(), mod_101.clone())
        );
        assert_eq!(p7, p4.add(&p3, &p1));
        assert_eq!(p7, p6.add(&p1, &p5));
        // # p9 = p5 + p4
        let p9 = p5.add(&p4, &p1);
        assert_eq!(p9, Point::new(56.into(), 99.into(), a_24, mod_101));
        assert_eq!(p9, p6.add(&p3, &p3));
        assert_eq!(p9, p7.add(&p2, &p5));
        assert_eq!(p9, p8.add(&p1, &p7));

        assert_eq!(p5, p1.mont_ladder(&5u32.into()));
        assert_eq!(p9, p1.mont_ladder(&9u32.into()));
        assert_eq!(p16, p1.mont_ladder(&16u32.into()));
        assert_eq!(p9, p3.mont_ladder(&3u32.into()));
    }
}
