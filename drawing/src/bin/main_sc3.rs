#![allow(dead_code, warnings, unused)]

use std::rc::Rc;

use algebraeon_drawing::canvas::Canvas;
use algebraeon_drawing::canvas2d::Canvas2D;
use algebraeon_drawing::canvas2d::MouseWheelZoomCamera;
use algebraeon_drawing::canvas2d::shapes::Shape;
use algebraeon_drawing::canvas2d::shapes::simplicial_complex_shapes;
use algebraeon_drawing::colour::Colour;
use algebraeon_geometry::simplexes::ConvexHull;
use algebraeon_geometry::simplexes::LabelledSimplicialDisjointUnion;
use algebraeon_geometry::simplexes::OrientationSide;
use algebraeon_geometry::simplexes::OrientedSimplex;
use algebraeon_geometry::simplexes::Simplex;
use algebraeon_geometry::*;
use algebraeon_nzq::*;
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

    let space = AffineSpace::new_linear(Rational::structure_ref(), 2);

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

    let x = x.clone();
    let x = x.refine_to_partial_simplicial_complex();
    let x = x.simplify();

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

    let mut canvas = Canvas2D::new(Box::new(MouseWheelZoomCamera::new()));
    canvas.plot_shapes(
        [Shape::SetThickness(0.3)]
            .into_iter()
            .chain(simplicial_complex_shapes(
                &Colour::green(),
                &Colour::green().darken(),
                0.5,
                &x,
            )),
    );
    canvas.run();
}
