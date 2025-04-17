#![allow(dead_code, warnings, unused)]

use std::rc::Rc;

use algebraeon_drawing::canvas::Canvas;
use algebraeon_drawing::canvas::canvas2d::*;
use algebraeon_geometry::simplexes::ConvexHull;
use algebraeon_geometry::simplexes::LabelledSimplicialDisjointUnion;
use algebraeon_geometry::simplexes::OrientationSide;
use algebraeon_geometry::simplexes::OrientedSimplex;
use algebraeon_geometry::simplexes::Simplex;
use algebraeon_geometry::*;
use algebraeon_nzq::*;
use algebraeon_rings::rings::isolated_algebraic::real::RealAlgebraic;
use algebraeon_rings::structure::PositiveRealNthRootSignature;
use algebraeon_rings::structure::RingSignature;
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

    let field = RealAlgebraic::structure();
    let space = AffineSpace::new_linear(field.clone(), 2);

    let sqrt2 = field
        .nth_root(&field.from_int(Integer::from(2)), 2)
        .unwrap();
    let sqrt3 = field
        .nth_root(&field.from_int(Integer::from(3)), 2)
        .unwrap();

    let a = ConvexHull::new(
        &space,
        vec![
            Vector::new(
                &space,
                vec![
                    field.from_int(Integer::from(0)),
                    field.from_int(Integer::from(0)),
                ],
            ),
            Vector::new(
                &space,
                vec![
                    field.from_int(Integer::from(0)),
                    field.from_int(Integer::from(1)),
                ],
            ),
            Vector::new(
                &space,
                vec![sqrt3.clone(), field.from_int(Integer::from(0))],
            ),
        ],
    )
    .as_simplicial_complex()
    .forget_labels();

    let b = ConvexHull::new(
        &space,
        vec![
            Vector::new(
                &space,
                vec![
                    field.from_int(Integer::from(0)),
                    field.from_int(Integer::from(0)),
                ],
            ),
            Vector::new(
                &space,
                vec![sqrt2.clone(), field.from_int(Integer::from(1))],
            ),
            Vector::new(
                &space,
                vec![sqrt2.clone(), field.from_int(Integer::from(0))],
            ),
        ],
    )
    .as_simplicial_complex()
    .into_forget_labels();

    let x = LabelledSimplicialDisjointUnion::union_raw(&(&a).into(), &(&b).into())
        .refine_to_partial_simplicial_complex()
        .closure();

    for spx in x.simplexes() {
        for pt in spx.points() {
            println!("{} {}", pt.coordinate(0), pt.coordinate(1));
        }
    }

    algebraeon_drawing::canvas::canvas2d::Diagram2dCanvas::run(|canvas| {
        // canvas.draw(&a, (1.0, 0.0, 0.0));
        // canvas.draw(&b, (0.0, 1.0, 0.0));
        // canvas.draw(&x, (1.0, 0.0, 0.0));
        canvas.draw(&x, (0.0, 1.0, 0.0));
    });
}
