//TODO: make structs for intervals and boxes instead of using tuples of rationals

use crate::number::integer::*;
use crate::number::natural::*;
use crate::number::rational::*;

pub fn add_intervals(
    first: (&Rational, &Rational),
    second: (&Rational, &Rational),
) -> (Rational, Rational) {
    (first.0 + second.0, first.1 + second.1)
}

pub fn add_interval_rat(interval: (&Rational, &Rational), rat: &Rational) -> (Rational, Rational) {
    (interval.0 + rat, interval.1 + rat)
}

pub fn mul_intervals(
    first: (&Rational, &Rational),
    second: (&Rational, &Rational),
) -> (Rational, Rational) {
    let mut pts = vec![];
    for x in [first.0, first.1] {
        for y in [second.0, second.1] {
            pts.push(x * y);
        }
    }
    (
        pts.iter().min().unwrap().clone(),
        pts.iter().max().unwrap().clone(),
    )
}

pub fn mul_interval_rat(interval: (&Rational, &Rational), rat: &Rational) -> (Rational, Rational) {
    match rat.cmp(&Rational::ZERO) {
        std::cmp::Ordering::Less => (rat * interval.1, rat * interval.0),
        std::cmp::Ordering::Equal => panic!(),
        std::cmp::Ordering::Greater => (rat * interval.0, rat * interval.1),
    }
}
