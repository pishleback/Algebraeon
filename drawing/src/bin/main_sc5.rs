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
use orthoclase_geometry::*;
use orthoclase_rings::ring_structure::cannonical::*;
use orthoclase_rings::structure::CannonicalStructure;
use orthoclase_rings::structure::StructuredType;
use rand::Rng;
use simplexes::LabelledSimplicialComplex;

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

    let a = SimplicialDisjointUnion::from(
        &ConvexHull::new(
            &space,
            vec![
                Vector::new(&space, vec![Rational::from(0), Rational::from(3)]),
                Vector::new(&space, vec![Rational::from(3), Rational::from(0)]),
                Vector::new(&space, vec![Rational::from(0), Rational::from(-3)]),
                Vector::new(&space, vec![Rational::from(-3), Rational::from(0)]),
            ],
        )
        .as_simplicial_complex()
        .interior,
    );

    let b = SimplicialDisjointUnion::from(
        &ConvexHull::new(
            &space,
            vec![
                Vector::new(&space, vec![Rational::from(-2), Rational::from(-2)]),
                Vector::new(&space, vec![Rational::from(2), Rational::from(-2)]),
                Vector::new(&space, vec![Rational::from(-2), Rational::from(2)]),
                Vector::new(&space, vec![Rational::from(2), Rational::from(2)]),
            ],
        )
        .as_simplicial_complex()
        .entire,
    );

    let x = a.union_raw(&b);

    let c = SimplicialDisjointUnion::from(
        &ConvexHull::new(
            &space,
            vec![
                Vector::new(&space, vec![Rational::from(-5), Rational::from(0)]),
                Vector::new(&space, vec![Rational::from(5), Rational::from(1)]),
                Vector::new(&space, vec![Rational::from(0), Rational::from(2)]),
            ],
        )
        .as_simplicial_complex()
        .entire,
    );

    let x = x
        .union_raw(&c)
        .refine_to_partial_simplicial_complex()
        .closure_as_simplicial_complex();
    // let y = x.clone().simplify();

    orthoclase_drawing::canvas::canvas2d::Diagram2dCanvas::run(|canvas| {
        canvas.draw(&x, (1.0, 0.0, 0.0));
        // canvas.draw(&y, (0.0, 1.0, 0.0));
    });
}
