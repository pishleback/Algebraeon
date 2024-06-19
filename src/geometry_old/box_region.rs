use malachite_q::Rational;

use super::vector::Vector;


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

pub struct BoxRegion {
    intervals: Vec<Interval>,
}

impl BoxRegion {
    pub fn bounding_box(dim: usize, points: &Vec<Vector>) -> Self {
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

    pub fn intersects(&self, other: &Self) -> bool {
        let dim = self.intervals.len();
        assert_eq!(dim, other.intervals.len());
        (0..dim).any(|i| self.intervals[i].intersects(&other.intervals[i]))
    }
}