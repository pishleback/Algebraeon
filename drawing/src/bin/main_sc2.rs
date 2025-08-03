#![allow(dead_code, warnings, unused)]

use algebraeon_drawing::canvas::Canvas;
use algebraeon_drawing::canvas2d::Canvas2D;
use algebraeon_drawing::canvas2d::MouseWheelZoomCamera;
use algebraeon_drawing::canvas2d::shapes::Shape;
use algebraeon_drawing::canvas2d::shapes::simplicial_complex_shapes;
use algebraeon_drawing::colour::Colour;
use algebraeon_geometry::ambient_space::AffineSpace;
use algebraeon_geometry::convex_hull::ConvexHull;
use algebraeon_geometry::oriented_simplex::OrientationSide;
use algebraeon_geometry::oriented_simplex::OrientedSimplex;
use algebraeon_geometry::simplex_collection::LabelledSimplexCollection;
use algebraeon_geometry::simplicial_disjoint_union::SimplicialDisjointUnion;
use algebraeon_geometry::vector::Vector;
use algebraeon_geometry::*;
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;
use std::collections::HashSet;
use std::str::FromStr;

fn main() {
    // let space = AffineSpace::new_linear(Rational::structure(), 2);
    // let p1 = Vector::new( space, vec![Rational::from(0), Rational::from(0)]);
    // let p2 = Vector::new( space, vec![Rational::from(1), Rational::from(0)]);
    // let p3 = Vector::new( space, vec![Rational::from(0), Rational::from(1)]);

    // let s1 = Simplex::new( space, vec![p1.clone()]).unwrap();
    // let s2 = Simplex::new( space, vec![p1.clone(), p2.clone()]).unwrap();
    // let s3 = Simplex::new( space, vec![p1.clone(), p2.clone(), p3.clone()]).unwrap();

    let space = AffineSpace::new_linear(Rational::structure_ref(), 2);
    let mut ch = ConvexHull::new(
        space,
        vec![
            space.vector([1, 1]),
            space.vector([1, 1]),
            space.vector([Rational::ZERO, Rational::from_str("1/2").unwrap()]),
            space.vector([-1, 0]),
            space.vector([-1, 0]),
            space.vector([Rational::ZERO, Rational::from_str("1/2").unwrap()]),
            space.vector([
                Rational::from_str("1/2").unwrap(),
                Rational::from_str("3/4").unwrap(),
            ]),
            space.vector([2, -1]),
            space.vector([0, -2]),
            space.vector([2, 0]),
            space.vector([2, 2]),
        ],
    );

    let ospx = OrientedSimplex::new_with_positive_point(
        space,
        vec![space.vector([0, 4]), space.vector([1, -4])],
        &space.vector([Rational::from(10), Rational::from(0)]),
    )
    .unwrap();

    let ohsp = ospx.clone().into_oriented_hyperplane();

    let smaller_ch_neutral = ch.intersect_with_oriented_hyperplane(&ohsp, OrientationSide::Neutral);
    let smaller_ch_pos = ch.intersect_with_oriented_hyperplane(&ohsp, OrientationSide::Positive);
    let smaller_ch_neg = ch.intersect_with_oriented_hyperplane(&ohsp, OrientationSide::Negative);

    let mut canvas = Canvas2D::new(Box::new(MouseWheelZoomCamera::new()));
    canvas.plot_shapes(
        [Shape::SetThickness(0.3)]
            .into_iter()
            .chain(simplicial_complex_shapes(
                &Colour::green(),
                &Colour::green().darken(),
                0.5,
                &smaller_ch_neutral.as_simplicial_complex(),
            ))
            .chain(simplicial_complex_shapes(
                &Colour::cyan(),
                &Colour::cyan().darken(),
                0.5,
                &smaller_ch_pos.as_simplicial_complex(),
            ))
            .chain(simplicial_complex_shapes(
                &Colour::blue(),
                &Colour::blue().darken(),
                0.5,
                &smaller_ch_neg.as_simplicial_complex(),
            ))
            .chain(simplicial_complex_shapes(
                &Colour::red(),
                &Colour::red().darken(),
                0.5,
                &SimplicialDisjointUnion::new_unchecked(space, [ospx.simplex().clone()].into()),
            )),
    );
    canvas.run();
}
