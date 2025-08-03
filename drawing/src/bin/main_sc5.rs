#![allow(dead_code, warnings, unused)]

use algebraeon_drawing::canvas::Canvas;
use algebraeon_drawing::canvas2d::Canvas2D;
use algebraeon_drawing::canvas2d::MouseWheelZoomCamera;
use algebraeon_drawing::canvas2d::shapes::Shape;
use algebraeon_drawing::canvas2d::shapes::simplicial_complex_shapes;
use algebraeon_drawing::colour::Colour;
use algebraeon_geometry::ambient_space::AffineSpace;
use algebraeon_geometry::convex_hull::ConvexHull;
use algebraeon_geometry::coordinates::Vector;
use algebraeon_geometry::simplex_collection::LabelledSimplexCollection;
use algebraeon_geometry::simplicial_complex::InteriorOrBoundary;
use algebraeon_geometry::simplicial_disjoint_union::LabelledSimplicialDisjointUnion;
use algebraeon_geometry::*;
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;
use rand::Rng;

fn main() {
    // let space = AffineSpace::new_linear(Rational::structure(), 2);
    // let p1 = Vector::new(space.clone(), vec![Rational::from(0), Rational::from(0)]);
    // let p2 = Vector::new(space.clone(), vec![Rational::from(1), Rational::from(0)]);
    // let p3 = Vector::new(space.clone(), vec![Rational::from(0), Rational::from(1)]);

    // let s1 = Simplex::new(space.clone(), vec![p1.clone()]).unwrap();
    // let s2 = Simplex::new(space.clone(), vec![p1.clone(), p2.clone()]).unwrap();
    // let s3 = Simplex::new(space.clone(), vec![p1.clone(), p2.clone(), p3.clone()]).unwrap();

    let space = AffineSpace::new_linear(Rational::structure_ref(), 2);

    let a = LabelledSimplicialDisjointUnion::from(
        &ConvexHull::new(
            space,
            vec![
                space.vector([0, 3]),
                space.vector([3, 0]),
                space.vector([0, -3]),
                space.vector([-3, 0]),
            ],
        )
        .as_simplicial_complex()
        .subset_by_label(&InteriorOrBoundary::Interior),
    );

    let b = LabelledSimplicialDisjointUnion::from(
        &ConvexHull::new(
            space,
            vec![
                space.vector([-2, -2]),
                space.vector([2, -2]),
                space.vector([-2, 2]),
                space.vector([2, 2]),
            ],
        )
        .as_simplicial_complex()
        .into_forget_labels(),
    );

    let x = a.union_raw(&b);

    let c = LabelledSimplicialDisjointUnion::from(
        &ConvexHull::new(
            space,
            vec![
                space.vector([-5, 0]),
                space.vector([5, 1]),
                space.vector([0, 2]),
            ],
        )
        .as_simplicial_complex()
        .into_forget_labels(),
    );

    let x = x.union_raw(&c).refine_to_partial_simplicial_complex();
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
