#![allow(dead_code, warnings, unused)]

use std::rc::Rc;

use malachite_nz::integer::Integer;
use malachite_nz::natural::Natural;
use malachite_q::arithmetic::traits::Approximate;
use malachite_q::Rational;
use orthoclase_all::drawing::canvas2d::*;
use orthoclase_all::drawing::Canvas;
use orthoclase_all::geometry::*;
use orthoclase_all::rings::ring_structure::cannonical::*;
use orthoclase_all::rings::structure::CannonicalStructure;
use orthoclase_all::rings::structure::StructuredType;
use rand::Rng;
use simplexes::ConvexHull;
use simplexes::OrientationSide;
use simplexes::OrientedSimplex;
use simplexes::Simplex;
use simplexes::SimplicialDisjointUnion;
use simplexes::VennResult;

fn main() {
    // let space = AffineSpace::new_linear(Rational::structure(), 2);
    // let p1 = Vector::new(&space, vec![Rational::from(0), Rational::from(0)]);
    // let p2 = Vector::new(&space, vec![Rational::from(1), Rational::from(0)]);
    // let p3 = Vector::new(&space, vec![Rational::from(0), Rational::from(1)]);

    // let s1 = Simplex::new(&space, vec![p1.clone()]).unwrap();
    // let s2 = Simplex::new(&space, vec![p1.clone(), p2.clone()]).unwrap();
    // let s3 = Simplex::new(&space, vec![p1.clone(), p2.clone(), p3.clone()]).unwrap();

    let space = AffineSpace::new_linear(Rational::structure(), 2);

    let a = ConvexHull::new(
        &space,
        vec![
            Vector::new(&space, vec![Rational::from(-1), Rational::from(-1)]),
            Vector::new(&space, vec![Rational::from(1), Rational::from(-1)]),
            Vector::new(&space, vec![Rational::from(-1), Rational::from(1)]),
            Vector::new(&space, vec![Rational::from(1), Rational::from(1)]),
        ],
    )
    .as_simplicial_complex()
    .entire;

    let b = ConvexHull::new(
        &space,
        vec![
            Vector::new(&space, vec![Rational::from(0), Rational::from(1)]),
            Vector::new(&space, vec![Rational::from(1), Rational::from(0)]),
            Vector::new(&space, vec![Rational::from(0), Rational::from(0)]),
        ],
    )
    .as_simplicial_complex()
    .entire;
    let a = SimplicialDisjointUnion::subtract(&a.clone().into(), &b.clone().into())
        .closure_as_simplicial_complex();

    let b = ConvexHull::new(
        &space,
        vec![
            Vector::new(&space, vec![Rational::from(-1), Rational::from(0)]),
            Vector::new(&space, vec![Rational::from(0), Rational::from(-1)]),
            Vector::new(&space, vec![Rational::from(-1), Rational::from(-1)]),
        ],
    )
    .as_simplicial_complex()
    .entire;
    let a = SimplicialDisjointUnion::subtract(&a.clone().into(), &b.clone().into())
        .closure_as_simplicial_complex();

    let b = ConvexHull::new(
        &space,
        vec![
            Vector::new(&space, vec![Rational::from(-1), Rational::from(0)]),
            Vector::new(&space, vec![Rational::from(0), Rational::from(1)]),
            Vector::new(&space, vec![Rational::from(-1), Rational::from(1)]),
        ],
    )
    .as_simplicial_complex()
    .entire;
    let a = SimplicialDisjointUnion::subtract(&a.clone().into(), &b.clone().into())
        .closure_as_simplicial_complex();

    let a = a.simplify();

    orthoclase_all::drawing::canvas2d::Diagram2dCanvas::run(|canvas| {
        // canvas.draw(&a, (1.0, 0.0, 0.0));
        // canvas.draw(&b, (0.0, 1.0, 0.0));
        canvas.draw(&a, (1.0, 1.0, 0.0));
    });
}
