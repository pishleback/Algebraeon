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
use algebraeon_geometry::parse::parse_shape;
use algebraeon_geometry::simplex_collection::InteriorOrBoundarySimplexCollection;
use algebraeon_geometry::simplex_collection::LabelledSimplexCollection;
use algebraeon_geometry::simplicial_disjoint_union::LabelledSimplicialDisjointUnion;
use algebraeon_geometry::vector::Vector;
use algebraeon_geometry::*;
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;
use rand::Rng;

fn plot_shape(shape: impl LabelledSimplexCollection<'static, RationalCanonicalStructure, ()>) {
    let x = shape.into_partial_simplicial_complex();
    let y = x.closure().difference(&x);
    let z = y.closure().difference(&y);

    let mut canvas = Canvas2D::new(Box::new(MouseWheelZoomCamera::new()));
    canvas.plot_shapes(
        [Shape::SetThickness(0.3)]
            .into_iter()
            .chain(simplicial_complex_shapes(
                &Colour::green(),
                &Colour::green().darken(),
                1.0,
                &x,
            )),
    );
    canvas.plot_shapes(
        [Shape::SetThickness(0.3)]
            .into_iter()
            .chain(simplicial_complex_shapes(
                &Colour::red(),
                &Colour::red().darken(),
                1.0,
                &y,
            )),
    );
    canvas.plot_shapes(
        [Shape::SetThickness(0.3)]
            .into_iter()
            .chain(simplicial_complex_shapes(
                &Colour::green(),
                &Colour::green().darken(),
                1.0,
                &z,
            )),
    );
    canvas.run();
}

fn main() {
    let x = parse_shape(
        "ConvexHull((0, 0), (1, 0), (0, 1), (1, 1)) \\ ((1/2, 1/2) | (1/3, 1/3) | (1/4, 1/4) | (1/5, 1/5) | (1/6, 1/6) | (1/3, 4/5))",
    );

    let x = parse_shape(
        "ConvexHull((0, 0), (3, 0), (0, 3), (3, 3)) \\ ConvexHull((1, 1), (2, 1), (1, 2), (2, 2))",
    );

    let x = parse_shape("Polygon((0, 0), (0, 3), (1/2, -1), (3, 3), (3, 0), (5/2, 4))");

    let x = parse_shape(
        "(Polygon((0, 0), (4, 0), (4, 1), (3, 1), (3, -1), (2, -1), (2, 2), (1, 2)) | Polygon((5/2, 0/2), (4/2, 1/2), (3/2, 0/2), (4/2, -1/2))) \\ (Loop((0, 0), (1, 1), (11/4, -1/4)) \\ (1/3, 1/3))",
    );

    plot_shape(x);
}
