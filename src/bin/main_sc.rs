#![allow(dead_code, warnings, unused)]

use malachite_nz::integer::Integer;
use malachite_nz::natural::Natural;
use malachite_q::Rational;
use orthoclase_all::drawing::canvas2d::*;
use orthoclase_all::drawing::Canvas;
use orthoclase_all::geometry::*;
use orthoclase_all::rings::structure::StructuredType;
use simplices::Simplex;

fn main() {
    let space = AffineSpace::new_linear(Rational::structure(), 2);
    let p1 = Vector::new(&space, vec![Rational::from(0), Rational::from(0)]);
    let p2 = Vector::new(&space, vec![Rational::from(1), Rational::from(0)]);
    let p3 = Vector::new(&space, vec![Rational::from(0), Rational::from(1)]);

    let s1 = Simplex::new(&space, vec![p1.clone()]).unwrap();
    let s2 = Simplex::new(&space, vec![p1.clone(), p2.clone()]).unwrap();
    let s3 = Simplex::new(&space, vec![p1.clone(), p2.clone(), p3.clone()]).unwrap();

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

        canvas.draw(s1, (1.0, 0.0, 0.0));
        canvas.draw(s2, (1.0, 0.0, 0.0));
        canvas.draw(s3, (1.0, 0.0, 0.0));
    });
}
