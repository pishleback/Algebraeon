#![allow(dead_code, warnings, unused)]

use std::rc::Rc;

use geometry::AffineSpace;
use geometry::Vector;
use malachite_nz::integer::Integer;
use malachite_nz::natural::Natural;
use malachite_q::arithmetic::traits::Approximate;
use malachite_q::Rational;
use orthoclase_drawing::canvas::canvas2d::*;
use orthoclase_drawing::canvas::Canvas;
use orthoclase_geometry::*;
use orthoclase_rings::ring_structure::cannonical::*;
use orthoclase_rings::structure::CannonicalStructure;
use orthoclase_rings::structure::StructuredType;
use rand::Rng;
use orthoclase_geometry::geometry::simplexes::ConvexHull;
use orthoclase_geometry::geometry::simplexes::OrientationSide;
use orthoclase_geometry::geometry::simplexes::OrientedSimplex;
use orthoclase_geometry::geometry::simplexes::Simplex;
use orthoclase_geometry::geometry::simplexes::SimplicialDisjointUnion;
use orthoclase_geometry::geometry::simplexes::VennResult;

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

    //This is where I left off.
    //As you can see, it is not fully simplified
    //The solution will be to simplify points on flat parts of the boundary as if they are interior points of the boundary in a lower dimension

    orthoclase_drawing::canvas::canvas2d::Diagram2dCanvas::run(|canvas| {
        // canvas.draw(&a, (1.0, 0.0, 0.0));
        // canvas.draw(&b, (0.0, 1.0, 0.0));
        canvas.draw(&a, (1.0, 1.0, 0.0));
    });
}
