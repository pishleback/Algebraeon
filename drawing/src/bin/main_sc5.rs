#![allow(dead_code, warnings, unused)]

use algebraeon_drawing::canvas::Canvas;
use algebraeon_drawing::canvas2d::Canvas2D;
use algebraeon_drawing::canvas2d::MouseWheelZoomCamera;
use algebraeon_drawing::canvas2d::shapes::Shape;
use algebraeon_drawing::canvas2d::shapes::simplicial_complex_shapes;
use algebraeon_drawing::colour::Colour;
use algebraeon_geometry::ambient_space::AffineSpace;
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

fn main() {
    let space = AffineSpace::new_linear(Rational::structure_ref(), 2);

    let a = space
        .convex_hull(vec![
            space.vector([0, 3]),
            space.vector([3, 0]),
            space.vector([0, -3]),
            space.vector([-3, 0]),
        ])
        .to_simplicial_complex()
        .interior()
        .into_simplicial_disjoint_union();

    let b = space
        .convex_hull(vec![
            space.vector([-2, -2]),
            space.vector([2, -2]),
            space.vector([-2, 2]),
            space.vector([2, 2]),
        ])
        .to_simplicial_complex()
        .into_forget_labels()
        .into_simplicial_disjoint_union();

    let x = a.union(&b);

    let c = space
        .convex_hull(vec![
            space.vector([-5, 0]),
            space.vector([5, 1]),
            space.vector([0, 2]),
        ])
        .to_simplicial_complex()
        .into_forget_labels()
        .into_simplicial_disjoint_union();

    let x = x.union(&c);
    let y = x.clone().simplify();

    let mut canvas = Canvas2D::new(Box::new(MouseWheelZoomCamera::new()));
    canvas.plot_shapes(
        [Shape::SetThickness(0.3)]
            .into_iter()
            .chain(simplicial_complex_shapes(
                &Colour::green(),
                &Colour::green().darken(),
                0.5,
                &y,
            )),
    );
    canvas.run();
}
