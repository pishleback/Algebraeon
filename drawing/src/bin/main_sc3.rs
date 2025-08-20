#![allow(dead_code, warnings, unused)]

use algebraeon_drawing::canvas::Canvas;
use algebraeon_drawing::canvas2d::Canvas2D;
use algebraeon_drawing::canvas2d::MouseWheelZoomCamera;
use algebraeon_drawing::canvas2d::shapes::Shape;
use algebraeon_drawing::canvas2d::shapes::simplicial_complex_shapes;
use algebraeon_drawing::colour::Colour;
use algebraeon_geometry::ambient_space::AffineSpace;
use algebraeon_geometry::boolean_operations::Difference;
use algebraeon_geometry::boolean_operations::Union;
use algebraeon_geometry::convex_hull::ConvexHull;
use algebraeon_geometry::simplex_collection::InteriorOrBoundarySimplexCollection;
use algebraeon_geometry::simplex_collection::LabelledSimplexCollection;
use algebraeon_geometry::simplicial_disjoint_union::LabelledSimplicialDisjointUnion;
use algebraeon_geometry::vector::Vector;
use algebraeon_geometry::*;
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;
use rand::Rng;
use std::rc::Rc;

fn main() {
    unsafe {
        std::env::set_var("RUST_BACKTRACE", "1");
    }

    let space = AffineSpace::new_linear(Rational::structure_ref(), 2);

    let a = space
        .convex_hull(vec![
            space.vector([0, 3]),
            space.vector([3, 0]),
            space.vector([0, -3]),
            space.vector([-3, 0]),
        ])
        .to_simplicial_complex()
        .interior();
    let x = a;

    let b = space
        .convex_hull(vec![
            space.vector([-2, -2]),
            space.vector([2, -2]),
            space.vector([-2, 2]),
            space.vector([2, 2]),
        ])
        .to_simplicial_complex()
        .into_forget_labels();
    let x = x.union(&b);

    let c = space
        .convex_hull(vec![
            space.vector([-1, -1]),
            space.vector([1, -1]),
            space.vector([-1, 1]),
            space.vector([1, 1]),
        ])
        .to_simplicial_complex()
        .into_forget_labels();
    let x = x.difference(&c);

    let y = x.closure().difference(&x);

    let mut canvas = Canvas2D::new(Box::new(MouseWheelZoomCamera::new()));
    canvas.plot_shapes(
        [Shape::SetThickness(0.3)]
            .into_iter()
            .chain(simplicial_complex_shapes(
                &Colour::green(),
                &Colour::green(),
                0.5,
                &x,
            )),
    );
    canvas.plot_shapes(
        [Shape::SetThickness(0.3)]
            .into_iter()
            .chain(simplicial_complex_shapes(
                &Colour::red(),
                &Colour::red(),
                0.5,
                &y,
            )),
    );
    canvas.run();
}
