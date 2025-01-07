#![allow(dead_code, warnings, unused)]

use std::rc::Rc;

use malachite_nz::integer::Integer;
use malachite_nz::natural::Natural;
use malachite_q::arithmetic::traits::Approximate;
use malachite_q::Rational;
use algebraeon_drawing::canvas::canvas2d::*;
use algebraeon_drawing::canvas::Canvas;
use algebraeon_geometry::simplexes::ConvexHull;
use algebraeon_geometry::simplexes::LabelledSimplicialDisjointUnion;
use algebraeon_geometry::simplexes::OrientationSide;
use algebraeon_geometry::simplexes::OrientedSimplex;
use algebraeon_geometry::simplexes::Simplex;
use algebraeon_geometry::*;
use algebraeon_sets::structure::*;
use rand::Rng;
use simplexes::LabelledSimplexCollection;

fn main() {
    // let space = AffineSpace::new_linear(Rational::structure(), 2);
    // let p1 = Vector::new(&space, vec![Rational::from(0), Rational::from(0)]);
    // let p2 = Vector::new(&space, vec![Rational::from(1), Rational::from(0)]);
    // let p3 = Vector::new(&space, vec![Rational::from(0), Rational::from(1)]);

    // let s1 = Simplex::new(&space, vec![p1.clone()]).unwrap();
    // let s2 = Simplex::new(&space, vec![p1.clone(), p2.clone()]).unwrap();
    // let s3 = Simplex::new(&space, vec![p1.clone(), p2.clone(), p3.clone()]).unwrap();

    let field = Rational::structure();

    let space = AffineSpace::new_linear(field, 2);

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
    .into_forget_labels();
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
    .into_forget_labels();
    let x = LabelledSimplicialDisjointUnion::union_raw(&(&x).into(), &(&b).into())
        .refine_to_partial_simplicial_complex();

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
    .into_forget_labels();
    let x = LabelledSimplicialDisjointUnion::subtract_raw(&(&x).into(), &(&c).into());

    let y = x.clone().refine_to_partial_simplicial_complex().simplify();

    // let y = x.clone().refine_to_partial_simplicial_complex().simplify();

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

    algebraeon_drawing::canvas::canvas2d::Diagram2dCanvas::run(|canvas| {
        // canvas.draw(&a, (1.0, 0.0, 0.0));
        // canvas.draw(&b, (0.0, 1.0, 0.0));
        // canvas.draw(&x, (1.0, 0.0, 0.0));
        canvas.draw(&y, (0.0, 1.0, 0.0));
    });
}
