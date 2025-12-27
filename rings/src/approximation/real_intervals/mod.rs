mod point;
mod subset;

pub use point::Point;
pub use point::PointInterface;
pub use point::continued_fraction::SimpleContinuedFractionPoint;
pub use point::rational::RationalPoint;
pub use subset::Subset;
pub use subset::SubsetsStructure;
pub use subset::subsets;

pub fn e() -> Point {
    Point::new(SimpleContinuedFractionPoint::from(
        crate::continued_fraction::eulers_constant(),
    ))
}

pub fn pi() -> Point {
    Point::new(point::pi::Pi::new())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::*;
    use algebraeon_nzq::{Integer, Natural, Rational};
    use std::str::FromStr;

    #[test]
    fn test_eulers_constant() {
        let e = e();

        // 2 < e < 3
        if let Subset::Interval(interval) = e.lock().rational_interval_neighbourhood() {
            assert_eq!(interval.a(), &Integer::from(2));
            assert_eq!(interval.b(), &Integer::from(3));
        } else {
            panic!()
        }

        e.lock().refine();

        // 8/3 < e < 3
        if let Subset::Interval(interval) = e.lock().rational_interval_neighbourhood() {
            assert_eq!(interval.a(), &Rational::from_str("8/3").unwrap());
            assert_eq!(interval.b(), &Integer::from(3));
        } else {
            panic!()
        }

        e.lock().refine();

        // 8/3 < e < 11/4
        if let Subset::Interval(interval) = e.lock().rational_interval_neighbourhood() {
            assert_eq!(interval.a(), &Rational::from_str("8/3").unwrap());
            assert_eq!(interval.b(), &Rational::from_str("11/4").unwrap());
        } else {
            panic!()
        }

        assert_eq!(e.as_f32(), std::f32::consts::E);
        assert_eq!(e.as_f64(), std::f64::consts::E);
    }

    #[test]
    fn test_pi() {
        let pi = pi();
        assert_eq!(pi.as_f32(), std::f32::consts::PI);
        assert_eq!(pi.as_f64(), std::f64::consts::PI);
    }

    #[test]
    fn test_rounding() {
        let pi = pi();
        assert_eq!(pi.floor(), Integer::from(3));
        assert_eq!(pi.round(), Integer::from(3));
        assert_eq!(pi.ceil(), Integer::from(4));

        let e = e();
        assert_eq!(e.floor(), Integer::from(2));
        assert_eq!(e.round(), Integer::from(3));
        assert_eq!(e.ceil(), Integer::from(3));
    }

    #[test]
    fn test_ring_opps() {
        let e = e();
        // compute e+(-e)^3
        let v = Point::add(&e, &Point::nat_pow(&Point::neg(&e), &Natural::from(3u32)));
        let f = v.as_f64();
        debug_assert_eq!(f, -17.367_255_094_728_623);
    }
}
