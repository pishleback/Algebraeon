#![allow(dead_code, warnings, unused)]

#[macro_use]
extern crate glium;

use std::marker::PhantomData;
use std::str::FromStr;
use std::task::Poll;

use malachite_nz::integer::Integer;
use malachite_nz::natural::Natural;
use malachite_q::Rational;
use orthoclase_all::drawing::canvas2d::*;
use orthoclase_all::drawing::Canvas;
use orthoclase_all::geometry_old::convex_simplicial_complex::*;
use orthoclase_all::geometry_old::vector::*;
use orthoclase_all::geometry_old::*;
use orthoclase_all::groups::group::*;
use orthoclase_all::groups::permutation::*;
use orthoclase_all::rings::linear::matrix::*;
use orthoclase_all::rings::number::algebraic::isolated_roots::*;
use orthoclase_all::rings::number::algebraic::number_field::*;
use orthoclase_all::rings::number::modulo::*;
use orthoclase_all::rings::polynomial::multipoly::*;
use orthoclase_all::rings::polynomial::polynomial::*;
use orthoclase_all::rings::ring_structure::cannonical::*;
use orthoclase_all::rings::ring_structure::cannonical::*;
use orthoclase_all::rings::ring_structure::quotient::*;
use orthoclase_all::rings::ring_structure::structure::*;
use orthoclase_all::rings::structure::*;

use orthoclase_all::rings::polynomial::polynomial::*;
use rand::Rng;

fn main() {
    fn random_point(dim: usize, rad: f64) -> Vector {
        let mut rng = rand::thread_rng();
        Vector::new(
            (0..dim)
                .map(|_i| Rational::from_f64_approx(rng.gen_range(-rad..rad)))
                .collect(),
        )
    }

    let shape = convex_hull(
        2,
        (0..5000)
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

    orthoclase_all::drawing::canvas2d::Shape2dCanvas::run(|canvas| {
        canvas.draw_shape(&a.as_shape(), (1.0, 0.0, 0.0));
        canvas.draw_shape(&b.as_shape(), (0.0, 1.0, 0.0));
    });

    //where I left with the simplicial stuff:
    //simplification of simplicial complex is done
    //combining simplicial complexes is not done
    //the plan is to recursively compare pairs of simplexes and refine them until every pair intersects trivially
}
