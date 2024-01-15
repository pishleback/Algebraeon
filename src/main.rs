#![allow(dead_code, warnings, unused)]

#[macro_use]
extern crate glium;

use std::str::FromStr;

use drawing::canvas2d::*;
use drawing::Canvas;
use geometry::*;
use itertools::Itertools;
use malachite_base::num::conversion::traits::IntegerMantissaAndExponent;
use malachite_nz::integer::Integer;
use malachite_nz::natural::Natural;
use malachite_q::Rational;
use rand::Rng;
use rings::algebraic::*;
use rings::ergonomic::*;
use rings::multipoly::*;
use rings::nzq::*;
use rings::poly::*;
use rings::ring::*;

use crate::geometry::convex_simplicial_complex::*;
use crate::geometry::vector::*;

pub mod drawing;
pub mod geometry;
pub mod groups;
pub mod numbers;
pub mod rings;
pub mod sets;

fn main() {
    fn random_point(dim: usize, rad: f64) -> Vector {
        let mut rng = rand::thread_rng();
        Vector::new(
            (0..dim)
                .map(|_i| QQ.from_f64_approx(rng.gen_range(-rad..rad)))
                .collect(),
        )
    }

    let shape = convex_hull(
        2,
        (0..24)
            .map(|i| random_point(2, f64::sqrt((i + 1) as f64)))
            .collect(),
    )
    .as_simplicial_complex();

    // let shape = convex_hull(
    //     2,
    //     vec![
    //         Point::new(vec![Rational::from(1), Rational::from(1)]),
    //         Point::new(vec![Rational::from(0), Rational::from(0)]),
    //         Point::new(vec![Rational::from(2), Rational::from(2)]),
    //     ],
    // )
    // .as_simplicial_complex();

    let a = shape;
    let b = a.clone().simplify();
    // let shape = shape.as_shape();

    // let (a, b) = shape.interior_and_boundary();
    // a.check().unwrap();
    // b.check().unwrap();

    let mut canvas = drawing::canvas2d::Shape2dCanvas::new();
    canvas.draw_shape(&a.as_shape(), (1.0, 0.0, 0.0));
    canvas.draw_shape(&b.as_shape(), (0.0, 1.0, 0.0));
    canvas.run();
}

fn todo() {
    // let s = NaturalPrimeGenerator::new();
    // for p in s {
    //     println!("{}", p);
    // }

    // let x = Integer::from(10000);
    // let f = ZZ.factor(&x);
    // println!("{:?}", f);

    let f = ZZ_POLY.from_coeffs(vec![
        Integer::from(1),
        Integer::from(0),
        Integer::from(0),
        Integer::from(0),
        Integer::from(0),
        Integer::from(1),
    ]);
    let roots = ZZ_POLY.all_complex_roots(&f);

    let a = QQ_BAR.sum(roots.iter().collect());

    println!("{:?}", a);
}

fn main2() {
    // (x+11)*(x-222)*(x-3333)*(x^2+x+1)*(x+67)

    let f = ZZ_POLY.product(vec![
        &ZZ_POLY.from_coeffs(vec![Integer::from(11), Integer::from(1)]),
        &ZZ_POLY.from_coeffs(vec![Integer::from(-222), Integer::from(1)]),
        &ZZ_POLY.from_coeffs(vec![Integer::from(-3333), Integer::from(1)]),
        &ZZ_POLY.from_coeffs(vec![Integer::from(1), Integer::from(1), Integer::from(1)]),
        &ZZ_POLY.from_coeffs(vec![
            Integer::from(-1),
            Integer::from(-1),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(1),
        ]),
        &ZZ_POLY.from_coeffs(vec![Integer::from(67), Integer::from(1)]),
        &ZZ_POLY.from_coeffs(vec![
            Integer::from(1),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(1),
        ]),
    ]);
    println!("{}", ZZ_POLY.to_string(&f));

    let fs = ZZ_POLY.factorize_by_hensel_lifting(f);
    println!("{:?}", fs);

    // let modn = EuclideanQuotient::new_field_unchecked(ZZ, Integer::from(5));

    // println!("{:?}", modn.all_units());

    // let poly_modn = PolynomialRing::new(&modn);

    // // let f = poly_modn.from_coeffs(vec![
    // //     Integer::from(-1),
    // //     Integer::from(1),
    // //     Integer::from(0),
    // //     Integer::from(0),
    // //     Integer::from(0),
    // //     Integer::from(0),
    // //     Integer::from(0),
    // //     Integer::from(0),
    // //     Integer::from(1),
    // // ]);

    // // let fs = ZZ_POLY.factor(&f);
    // // println!("{:?}", fs);

    // let fs = poly_modn.factorize_by_trying_all_factors(f);
    // println!("{:?}", fs);
}
