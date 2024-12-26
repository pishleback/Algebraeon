#![allow(dead_code, warnings, unused)]

use std::rc::Rc;

use malachite_nz::integer::Integer;
use malachite_nz::natural::Natural;
use malachite_q::arithmetic::traits::Approximate;
use malachite_q::Rational;
use orthoclase_drawing::canvas::canvas2d::*;
use orthoclase_drawing::canvas::Canvas;
use orthoclase_geometry::simplexes::ConvexHull;
use orthoclase_geometry::simplexes::OrientationSide;
use orthoclase_geometry::simplexes::OrientedSimplex;
use orthoclase_geometry::simplexes::Simplex;
use orthoclase_geometry::simplexes::SimplicialDisjointUnion;
use orthoclase_geometry::simplexes::VennResult;
use orthoclase_geometry::*;
use orthoclase_rings::ring_structure::cannonical::*;
use orthoclase_rings::structure::CannonicalStructure;
use orthoclase_rings::structure::StructuredType;
use rand::Rng;

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
            Vector::new(&space, vec![Rational::from(0), Rational::from(3)]),
            Vector::new(&space, vec![Rational::from(3), Rational::from(0)]),
            Vector::new(&space, vec![Rational::from(0), Rational::from(-3)]),
            Vector::new(&space, vec![Rational::from(-3), Rational::from(0)]),
        ],
    )
    .as_simplicial_complex()
    .entire;
    let x = a;

    let b = ConvexHull::new(
        &space,
        vec![
            Vector::new(&space, vec![Rational::from(-2), Rational::from(-2)]),
            Vector::new(&space, vec![Rational::from(2), Rational::from(-2)]),
            Vector::new(&space, vec![Rational::from(-2), Rational::from(2)]),
            Vector::new(&space, vec![Rational::from(2), Rational::from(2)]),
        ],
    )
    .as_simplicial_complex()
    .entire;
    let x = SimplicialDisjointUnion::union(&x.clone().into(), &b.clone().into())
        .closure_as_simplicial_complex();

    let c = ConvexHull::new(
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
    let x = SimplicialDisjointUnion::subtract(&x.clone().into(), &c.clone().into())
        .closure_as_simplicial_complex();

    let y = x.clone().simplify();

    // let mut y = ConvexHull::new_empty(&space);
    // y.extend_by_point(Vector::new(
    //     &space,
    //     vec![Rational::from(0), Rational::from(0)],
    // ));
    // y.extend_by_point(Vector::new(
    //     &space,
    //     vec![Rational::from(1), Rational::from(0)],
    // ));
    // y.extend_by_point(Vector::new(
    //     &space,
    //     vec![Rational::from(-1), Rational::from(0)],
    // ));
    // let y = y.as_simplicial_complex().entire;
    // let y = y.simplify();

    //This is where I left off.
    //As you can see, it is not fully simplified
    //The solution will be to simplify points on flat parts of the boundary as if they are interior points of the boundary in a lower dimension

    orthoclase_drawing::canvas::canvas2d::Diagram2dCanvas::run(|canvas| {
        // canvas.draw(&a, (1.0, 0.0, 0.0));
        // canvas.draw(&b, (0.0, 1.0, 0.0));
        // canvas.draw(&x, (1.0, 0.0, 0.0));
        canvas.draw(&y, (0.0, 1.0, 0.0));
    });
}
