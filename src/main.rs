#![allow(dead_code)]

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

mod drawing;
mod geometry;
mod groups;
mod numbers;
mod rings;
mod sets;

fn main() {
    // let canvas = drawing::canvas2d::Canvas::new();
    // canvas.run();

    

    fn make_shape() -> Shape {
        let mut rng = rand::thread_rng();

        let points = (0..6)
            .map(|i| {
                Point::new(vec![
                    Rational::from_sci_string_simplest(
                        rng.gen_range(0.0..10.0).to_string().as_str(),
                    )
                    .unwrap(),
                    Rational::from_sci_string_simplest(
                        rng.gen_range(0.0..10.0).to_string().as_str(),
                    )
                    .unwrap(),
                ])
            })
            .collect_vec();

        let shape = convexhull(2, points.clone());
        let pts_shape = Shape::new(
            2,
            points
                .iter()
                .map(|p| Simplex::new(2, vec![p.clone()]))
                .collect(),
        );

        shape
    }

    let a = make_shape();
    let b = make_shape();
    let c = a.union(&b);

    let mut canvas = drawing::canvas2d::Shape2dCanvas::new();
    canvas.draw_shape(&a, (1.0, 0.0, 0.0));
    canvas.draw_shape(&b, (0.0, 1.0, 0.0));
    canvas.draw_shape(&c, (1.0, 1.0, 0.0));
    canvas.run();

    return;

    /*

    let a = Simplex::new(
        2,
        vec![
            Point::new(vec![
                Rational::from_str("0").unwrap(),
                Rational::from_str("0").unwrap(),
            ]),
            Point::new(vec![
                Rational::from_str("2").unwrap(),
                Rational::from_str("0").unwrap(),
            ]),
            Point::new(vec![
                Rational::from_str("1").unwrap(),
                Rational::from_str("2").unwrap(),
            ]),
        ],
    );

    let b = Simplex::new(
        2,
        vec![
            Point::new(vec![
                Rational::from_str("0").unwrap(),
                Rational::from_str("1").unwrap(),
            ]),
            Point::new(vec![
                Rational::from_str("2").unwrap(),
                Rational::from_str("1").unwrap(),
            ]),
            Point::new(vec![
                Rational::from_str("1").unwrap(),
                Rational::from_str("-1").unwrap(),
            ]),
        ],
    );

    let a = Shape::simplex(a);

    let c = intersect_shape_full_simplex(&b, &a);

    println!("{:?}", c);

    for s in c.clone().simplices() {
        println!("{:?} {:?}", s.n(), s);
    }

    c.view2d();

    */

    /*
    let points = vec![
        Point::new(vec![
            Rational::from_str("0").unwrap(),
            Rational::from_str("0").unwrap(),
        ]),
        Point::new(vec![
            Rational::from_str("1").unwrap(),
            Rational::from_str("0").unwrap(),
        ]),
        Point::new(vec![
            Rational::from_str("0").unwrap(),
            Rational::from_str("1").unwrap(),
        ]),
    ];

    let b = convexhull(2, points);
    b.view2d();

    */

    /*
    let a = Simplex::new(
        3,
        vec![
            Point::new(vec![
                Rational::from_str("0").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("-1").unwrap(),
            ]),
            Point::new(vec![
                Rational::from_str("0").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("1").unwrap(),
            ]),
            Point::new(vec![
                Rational::from_str("2").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("0").unwrap(),
            ]),
        ],
    );

    let b = Simplex::new(
        3,
        vec![
            Point::new(vec![
                Rational::from_str("1").unwrap(),
                Rational::from_str("-1").unwrap(),
                Rational::from_str("0").unwrap(),
            ]),
            Point::new(vec![
                Rational::from_str("1").unwrap(),
                Rational::from_str("1").unwrap(),
                Rational::from_str("0").unwrap(),
            ]),
            Point::new(vec![
                Rational::from_str("3").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("0").unwrap(),
            ]),
        ],
    );

    let (c, d) = cut_simplex_by_simplex(&a, &b);

    println!("{:?}", c);
    println!("{:?}", d);

    */
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
