use algebraeon_nzq::natural::Natural;

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
#[derive(Default, Debug, Clone)]
pub struct Point {
    /// X coordinate of the Point
    pub x_cord: Natural,
    /// Z coordinate of the Point
    pub z_cord: Natural,
    /// Parameter of the elliptic curve in Montgomery form
    pub a_24: Natural,
    /// modulus
    pub modulus: Natural,
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
    pub fn new(x_cord: Natural, z_cord: Natural, a_24: Natural, modulus: Natural) -> Point {
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
        let u = (&self.modulus + &self.x_cord - &self.z_cord) * (&q.x_cord + &q.z_cord);
        let v = (&self.x_cord + &self.z_cord) * (&self.modulus + &q.x_cord - &q.z_cord);
        let add = &u + &v;
        let subt = u + (-v) % &self.modulus;
        let x_cord = (&diff.z_cord * &add * &add) % &self.modulus;
        let z_cord = (&diff.x_cord * &subt * &subt) % &self.modulus;

        Point::new(x_cord, z_cord, self.a_24.clone(), self.modulus.clone())
    }

    /// Doubles a point in an elliptic curve in Montgomery form.
    pub fn double(&self) -> Point {
        let u = (&self.x_cord + &self.z_cord) * (&self.x_cord + &self.z_cord);
        let v = (&self.modulus + &self.x_cord - &self.z_cord)
            * (&self.modulus + &self.x_cord - &self.z_cord);
        let diff = &u + (-&v) % &self.modulus;
        let x_cord = (u * &v) % &self.modulus;
        let z_cord = ((v + &self.a_24 * &diff) * diff) % &self.modulus;

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

        (self.z_cord.mod_inv_ref(modulus).unwrap() * &self.x_cord) % modulus
            == (other.z_cord.mod_inv_ref(modulus).unwrap() * &other.x_cord) % modulus
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_point_add() {
        let p1 = Point::new(11u32.into(), 16u32.into(), 7u32.into(), 29u32.into());
        let p2 = Point::new(13u32.into(), 10u32.into(), 7u32.into(), 29u32.into());
        let p3 = p2.add(&p1, &p1);

        assert_eq!(p3.x_cord, Natural::from(23u32));
        assert_eq!(p3.z_cord, Natural::from(17u32));
    }

    #[test]
    fn test_point_double() {
        let p1 = Point::new(11u32.into(), 16u32.into(), 7u32.into(), 29u32.into());
        let p2 = p1.double();

        assert_eq!(p2.x_cord, Natural::from(13u32));
        assert_eq!(p2.z_cord, Natural::from(10u32));
    }

    #[test]
    fn test_point_mont_ladder() {
        let p1 = Point::new(11u32.into(), 16u32.into(), 7u32.into(), 29u32.into());
        let p3 = p1.mont_ladder(&3u32.into());

        assert_eq!(p3.x_cord, Natural::from(23u32));
        assert_eq!(p3.z_cord, Natural::from(17u32));
    }

    #[test]
    fn test_point() {
        let a: Natural = 10u32.into();
        let a_24: Natural =
            (a + Natural::from(2u32)) * Natural::from(4u32).mod_inv(&101u32.into()).unwrap();

        let p1 = Point::new(10u32.into(), 17u32.into(), a_24.clone(), 101u32.into());
        let p2 = p1.double();
        assert_eq!(
            p2,
            Point::new(68u32.into(), 56u32.into(), a_24.clone(), 101u32.into())
        );
        let p4 = p2.double();
        assert_eq!(
            p4,
            Point::new(22u32.into(), 64u32.into(), a_24.clone(), 101u32.into())
        );
        let p8 = p4.double();
        assert_eq!(
            p8,
            Point::new(71u32.into(), 95u32.into(), a_24.clone(), 101u32.into())
        );
        let p16 = p8.double();
        assert_eq!(
            p16,
            Point::new(5u32.into(), 16u32.into(), a_24.clone(), 101u32.into())
        );
        let p32 = p16.double();
        assert_eq!(
            p32,
            Point::new(33u32.into(), 96u32.into(), a_24.clone(), 101u32.into())
        );

        // p3 = p2 + p1
        let p3 = p2.add(&p1, &p1);
        assert_eq!(
            p3,
            Point::new(1u32.into(), 61u32.into(), a_24.clone(), 101u32.into())
        );
        // p5 = p3 + p2 or p4 + p1
        let p5 = p3.add(&p2, &p1);
        assert_eq!(
            p5,
            Point::new(49u32.into(), 90u32.into(), a_24.clone(), 101u32.into())
        );
        assert_eq!(p5, p4.add(&p1, &p3));
        // # p6 = 2*p3
        let p6 = p3.double();
        assert_eq!(
            p6,
            Point::new(87u32.into(), 43u32.into(), a_24.clone(), 101u32.into())
        );
        assert_eq!(p6, p4.add(&p2, &p2));
        // # p7 = p5 + p2
        let p7 = p5.add(&p2, &p3);
        assert_eq!(
            p7,
            Point::new(69u32.into(), 23u32.into(), a_24.clone(), 101u32.into())
        );
        assert_eq!(p7, p4.add(&p3, &p1));
        assert_eq!(p7, p6.add(&p1, &p5));
        // # p9 = p5 + p4
        let p9 = p5.add(&p4, &p1);
        assert_eq!(
            p9,
            Point::new(56u32.into(), 99u32.into(), a_24, 101u32.into())
        );
        assert_eq!(p9, p6.add(&p3, &p3));
        assert_eq!(p9, p7.add(&p2, &p5));
        assert_eq!(p9, p8.add(&p1, &p7));

        assert_eq!(p5, p1.mont_ladder(&5u32.into()));
        assert_eq!(p9, p1.mont_ladder(&9u32.into()));
        assert_eq!(p16, p1.mont_ladder(&16u32.into()));
        assert_eq!(p9, p3.mont_ladder(&3u32.into()));
    }
}
