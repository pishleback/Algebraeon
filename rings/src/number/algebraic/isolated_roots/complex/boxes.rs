//TODO: make structs for intervals and boxes instead of using tuples of rationals

use malachite_base::num::basic::traits::Zero;
use malachite_q::Rational;

pub fn add_boxes(
    first: (&Rational, &Rational, &Rational, &Rational),
    second: (&Rational, &Rational, &Rational, &Rational),
) -> (Rational, Rational, Rational, Rational) {
    (
        first.0 + second.0,
        first.1 + second.1,
        first.2 + second.2,
        first.3 + second.3,
    )
}

pub fn add_box_rat(
    box_region: (&Rational, &Rational, &Rational, &Rational),
    rat: &Rational,
) -> (Rational, Rational, Rational, Rational) {
    (
        box_region.0 + rat,
        box_region.1 + rat,
        box_region.2.clone(),
        box_region.3.clone(),
    )
}

pub fn mul_boxes(
    first: (&Rational, &Rational, &Rational, &Rational),
    second: (&Rational, &Rational, &Rational, &Rational),
) -> (Rational, Rational, Rational, Rational) {
    let mut pts_re = vec![];
    let mut pts_im = vec![];
    for (re1, im1) in [
        (first.0, first.2),
        (first.0, first.3),
        (first.1, first.2),
        (first.1, first.3),
    ] {
        for (re2, im2) in [
            (second.0, second.2),
            (second.0, second.3),
            (second.1, second.2),
            (second.1, second.3),
        ] {
            pts_re.push(re1 * re2 - im1 * im2);
            pts_im.push(re1 * im2 + im1 * re2);
        }
    }

    (
        pts_re.iter().min().unwrap().clone(),
        pts_re.into_iter().max().unwrap(),
        pts_im.iter().min().unwrap().clone(),
        pts_im.into_iter().max().unwrap(),
    )
}

pub fn mul_box_rat(
    box_region: (&Rational, &Rational, &Rational, &Rational),
    rat: &Rational,
) -> (Rational, Rational, Rational, Rational) {
    match rat.cmp(&Rational::ZERO) {
        std::cmp::Ordering::Less => (
            rat * box_region.1,
            rat * box_region.0,
            rat * box_region.3,
            rat * box_region.2,
        ),
        std::cmp::Ordering::Equal => panic!(),
        std::cmp::Ordering::Greater => (
            rat * box_region.0,
            rat * box_region.1,
            rat * box_region.2,
            rat * box_region.3,
        ),
    }
}
