#![allow(dead_code, warnings, unused)]

#[macro_use]
extern crate glium;

use std::marker::PhantomData;
use std::str::FromStr;

// use drawing::canvas2d::*;
// use drawing::Canvas;
// use geometry::*;
use glium::glutin::event::MouseButton;
use groups::group::*;
use groups::permutation::*;
use itertools::Itertools;
use malachite_base::num::arithmetic::traits::Mod;
use malachite_base::num::conversion::traits::IntegerMantissaAndExponent;
use malachite_nz::integer::Integer;
use malachite_nz::natural::Natural;
use malachite_q::Rational;
use rand::Rng;
// use rings_old::ergonomic::*;
// use rings_old::numbers::algebraic::*;
// use rings_old::numbers::nzq::*;
// use rings_old::numbers::small_modulo::*;
// use rings_old::polynomial::multipoly::*;
// use rings_old::polynomial::poly::*;
// use rings_old::ring::*;
use rings::*;

// use sets::permutations::Permutation;

use crate::number::modulo::Modulo;
// use crate::geometry::convex_simplicial_complex::*;
// use crate::geometry::vector::*;
use crate::ring_structure::cannonical::*;
use crate::structure::*;

// pub mod drawing;
pub mod finite_group_tables;
// pub mod geometry;
pub mod groups;
pub mod numbers;
// pub mod rings_old;
pub mod rings;
pub mod sets;

/*
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

    drawing::canvas2d::Shape2dCanvas::run(|canvas| {
        canvas.draw_shape(&a.as_shape(), (1.0, 0.0, 0.0));
        canvas.draw_shape(&b.as_shape(), (0.0, 1.0, 0.0));
    });


    //where I left with the simplicial stuff:
    //simplification of simplicial complex is done
    //combining simplicial complexes is not done
    //the plan is to recursively compare pairs of simplexes and refine them until every pair intersects trivially
}
*/

// fn todo() {
//     // let s = NaturalPrimeGenerator::new();
//     // for p in s {
//     //     println!("{}", p);
//     // }

//     // let x = Integer::from(10000);
//     // let f = ZZ.factor(&x);
//     // println!("{:?}", f);

//     let f = Polynomial::from_coeffs(vec![
//         Integer::from(1),
//         Integer::from(0),
//         Integer::from(0),
//         Integer::from(0),
//         Integer::from(0),
//         Integer::from(1),
//     ]);
//     let roots = f.all_complex_roots();

//     for root in roots.iter() {
//         println!("{:?}", root);
//     }

//     let a = ComplexAlgebraic::sum(roots.iter().collect());

//     println!("{:?}", a);
// }

fn main() {
    // let x = Ergonomic::new(Polynomial::<Integer>::var());
    // let f = ((2 * x.pow(3) + 6 * x.pow(2) - 4) * (3 * x.pow(5) + 7 * x.pow(4) - 4)).elem();
    // println!("{}", f);
    // println!("{}", f.clone().factorize_by_kroneckers_method().unwrap());
    // println!("{}", f.clone().factorize_by_zassenhaus_algorithm().unwrap());

    let a = &Modulo::<5>::from(4).as_elem();
    let b = &Modulo::<5>::from(-3).as_elem();

    // let a = Integer::from(720);
    println!("{} + {} = {}", a, b, a + b);
}
