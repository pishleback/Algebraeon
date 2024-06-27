#![allow(dead_code, warnings, unused)]

use malachite_nz::integer::Integer;
use malachite_nz::natural::Natural;
use malachite_q::Rational;
use orthoclase_all::drawing::canvas2d::*;
use orthoclase_all::drawing::Canvas;
use orthoclase_all::geometry::*;
use orthoclase_all::rings::structure::StructuredType;
use simplexes::ConvexHull;
use simplexes::OrientationSide;
use simplexes::OrientedSimplex;
use simplexes::Simplex;

fn main() {
    // let space = AffineSpace::new_linear(Rational::structure(), 2);
    // let p1 = Vector::new(&space, vec![Rational::from(0), Rational::from(0)]);
    // let p2 = Vector::new(&space, vec![Rational::from(1), Rational::from(0)]);
    // let p3 = Vector::new(&space, vec![Rational::from(0), Rational::from(1)]);

    // let s1 = Simplex::new(&space, vec![p1.clone()]).unwrap();
    // let s2 = Simplex::new(&space, vec![p1.clone(), p2.clone()]).unwrap();
    // let s3 = Simplex::new(&space, vec![p1.clone(), p2.clone(), p3.clone()]).unwrap();

    let space = AffineSpace::new_linear(Rational::structure(), 2);
    let mut ch1 = ConvexHull::new(
        &space,
        vec![
            Vector::new(&space, vec![Rational::from(2), Rational::from(0)]),
            Vector::new(&space, vec![Rational::from(-2), Rational::from(-1)]),
            Vector::new(&space, vec![Rational::from(-1), Rational::from(1)]),
        ],
    );
    let ch2 = ConvexHull::new(
        &space,
        vec![
            Vector::new(&space, vec![Rational::from(0), Rational::from(-1)]),
            Vector::new(&space, vec![Rational::from(1), Rational::from(1)]),
            Vector::new(&space, vec![Rational::from(-1), Rational::from(0)]),
        ],
    );
    let ch3 = ch1.intersect(&ch2);

    orthoclase_all::drawing::canvas2d::Diagram2dCanvas::run(|canvas| {
        // canvas.draw_point((0.0, 0.0), (1.0, 0.0, 0.0));
        // canvas.draw_point((-1.0, -1.0), (1.0, 0.0, 0.0));
        // canvas.draw_point((1.0, -1.0), (1.0, 0.0, 0.0));
        // canvas.draw_point((-1.0, 1.0), (1.0, 0.0, 0.0));
        // canvas.draw_point((1.0, 1.0), (1.0, 0.0, 0.0));

        // canvas.draw_line((0.0, 0.0), (-1.0, -1.0), (1.0, 0.0, 0.0));
        // canvas.draw_line((0.0, 0.0), (1.0, -1.0), (1.0, 0.0, 0.0));
        // canvas.draw_line((0.0, 0.0), (-1.0, 1.0), (1.0, 0.0, 0.0));
        // canvas.draw_line((0.0, 0.0), (1.0, 1.0), (1.0, 0.0, 0.0));

        // canvas.draw_triangle((0.0, 0.0), (1.0, 1.0), (-1.0, 1.0), (1.0, 0.0, 0.0));
        // canvas.draw_triangle((0.0, 0.0), (1.0, 0.0), (1.0, -1.0), (1.0, 0.0, 0.0));
        // canvas.draw_triangle((0.0, 0.0), (-1.0, 0.0), (-1.0, -1.0), (1.0, 0.0, 0.0));

        // canvas.draw(ch.as_simplicial_complex().entire.as_ref(), (1.0, 1.0, 1.0));
        canvas.draw(ch1.as_simplicial_complex().entire.as_ref(), (1.0, 0.0, 0.0));
        canvas.draw(ch2.as_simplicial_complex().entire.as_ref(), (0.0, 1.0, 0.0));
        canvas.draw(ch3.as_simplicial_complex().entire.as_ref(), (1.0, 1.0, 0.0));
    });
}
