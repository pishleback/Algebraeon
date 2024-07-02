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
    let mut ch = ConvexHull::new(
        &space,
        vec![
            Vector::new(&space, vec![Rational::from(1), Rational::from(1)]),
            Vector::new(&space, vec![Rational::from(1), Rational::from(1)]),
            Vector::new(
                &space,
                vec![Rational::from(0), Rational::from(1) / Rational::from(2)],
            ),
            Vector::new(&space, vec![Rational::from(-1), Rational::from(0)]),
            Vector::new(&space, vec![Rational::from(-1), Rational::from(0)]),
            Vector::new(
                &space,
                vec![Rational::from(0), Rational::from(1) / Rational::from(2)],
            ),
            Vector::new(
                &space,
                vec![
                    Rational::from(1) / Rational::from(2),
                    Rational::from(3) / Rational::from(4),
                ],
            ),
            Vector::new(&space, vec![Rational::from(2), Rational::from(-1)]),
            Vector::new(&space, vec![Rational::from(0), Rational::from(-2)]),
            Vector::new(&space, vec![Rational::from(2), Rational::from(0)]),
            Vector::new(&space, vec![Rational::from(2), Rational::from(2)]),
        ],
    );

    let ospx = OrientedSimplex::new_with_positive_point(
        &space,
        vec![
            Vector::new(&space, vec![Rational::from(0), Rational::from(4)]),
            Vector::new(&space, vec![Rational::from(1), Rational::from(-4)]),
        ],
        &Vector::new(&space, vec![Rational::from(10), Rational::from(0)]),
    )
    .unwrap();

    let ohsp = ospx.clone().into_oriented_hyperplane();

    let smaller_ch_neutral = ch.intersect_with_oriented_hyperplane(&ohsp, OrientationSide::Neutral);
    let smaller_ch_pos = ch.intersect_with_oriented_hyperplane(&ohsp, OrientationSide::Positive);
    let smaller_ch_neg = ch.intersect_with_oriented_hyperplane(&ohsp, OrientationSide::Negative);

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
        canvas.draw(ospx.simplex(), (1.0, 0.0, 0.0));
        canvas.draw(
            &smaller_ch_neutral.as_simplicial_complex().entire,
            (0.0, 1.0, 0.0),
        );
        canvas.draw(
            &smaller_ch_pos.as_simplicial_complex().interior,
            (0.0, 1.0, 1.0),
        );
        canvas.draw(
            &smaller_ch_neg.as_simplicial_complex().interior,
            (0.0, 0.0, 1.0),
        );
    });
}
