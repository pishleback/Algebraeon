mod point;
mod subset;

pub use subset::Subset;
pub use subset::SubsetsStructure;
pub use subset::subsets;

pub use point::Point;
pub use point::PointsStructure;
pub use point::continued_fraction::SimpleContinuedFractionPoint;
pub use point::pi::pi;
pub use point::point;
pub use point::points;
pub use point::rational::RationalPoint;

#[cfg(test)]
mod tests {
    use crate::approximation::real_intervals::point::pi::pi;
    use crate::structure::*;
    use crate::{continued_fraction::eulers_constant, structure::RealSubsetSignature};
    use algebraeon_nzq::{Integer, Natural, Rational};
    use std::str::FromStr;

    use super::*;

    #[test]
    fn test_eulers_constant() {
        let mut e = SimpleContinuedFractionPoint::from(eulers_constant());

        // 2 < e < 3
        if let Subset::Interval(interval) = e.rational_interval_neighbourhood() {
            assert_eq!(interval.a(), &Integer::from(2));
            assert_eq!(interval.b(), &Integer::from(3));
        } else {
            panic!()
        }

        e.refine();

        // 8/3 < e < 3
        if let Subset::Interval(interval) = e.rational_interval_neighbourhood() {
            assert_eq!(interval.a(), &Rational::from_str("8/3").unwrap());
            assert_eq!(interval.b(), &Integer::from(3));
        } else {
            panic!()
        }

        e.refine();

        // 8/3 < e < 11/4
        if let Subset::Interval(interval) = e.rational_interval_neighbourhood() {
            assert_eq!(interval.a(), &Rational::from_str("8/3").unwrap());
            assert_eq!(interval.b(), &Rational::from_str("11/4").unwrap());
        } else {
            panic!()
        }

        assert_eq!(points().as_f32(&point(e.clone())), std::f32::consts::E);
        assert_eq!(points().as_f64(&point(e)), std::f64::consts::E);
    }

    #[test]
    fn test_pi() {
        let pi = point(pi());
        assert_eq!(points().as_f32(&pi), std::f32::consts::PI);
        assert_eq!(points().as_f64(&pi), std::f64::consts::PI);
    }

    #[test]
    fn test_rounding() {
        let pi = point(pi());
        assert_eq!(points().floor(&pi), Integer::from(3));
        assert_eq!(points().round(&pi), Integer::from(3));
        assert_eq!(points().ceil(&pi), Integer::from(4));

        let e = point(SimpleContinuedFractionPoint::from(eulers_constant()));
        assert_eq!(points().floor(&e), Integer::from(2));
        assert_eq!(points().round(&e), Integer::from(3));
        assert_eq!(points().ceil(&e), Integer::from(3));
    }

    #[test]
    fn test_ring_opps() {
        let e = point(SimpleContinuedFractionPoint::from(eulers_constant()));
        // compute e+(-e)^3
        let v = points().add(
            &e,
            &points().nat_pow(&points().neg(&e), &Natural::from(3u32)),
        );
        let f = points().as_f64(&v);
        debug_assert_eq!(f, -17.367_255_094_728_623);
    }
}
