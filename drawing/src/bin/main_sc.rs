#![allow(dead_code, warnings, unused)]

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
use algebraeon_rings::structure::*;
use algebraeon_sets::structure::*;
use rand::Rng;
use simplexes::LabelledSimplexCollection;
use std::rc::Rc;

fn main() {
    // let space = AffineSpace::new_linear(Rational::structure(), 2);
    // let p1 = Vector::new(&space, vec![Rational::from(0), Rational::from(0)]);
    // let p2 = Vector::new(&space, vec![Rational::from(1), Rational::from(0)]);
    // let p3 = Vector::new(&space, vec![Rational::from(0), Rational::from(1)]);

    // let s1 = Simplex::new(&space, vec![p1.clone()]).unwrap();
    // let s2 = Simplex::new(&space, vec![p1.clone(), p2.clone()]).unwrap();
    // let s3 = Simplex::new(&space, vec![p1.clone(), p2.clone(), p3.clone()]).unwrap();

    let space = AffineSpace::new_linear(Rational::structure_ref(), 2);

    fn random_point(
        space: &AffineSpace<'static, RationalCanonicalStructure>,
        rad: f64,
    ) -> Vector<'static, RationalCanonicalStructure> {
        let mut rng = rand::thread_rng();
        Vector::construct(space.clone(), |i| {
            Rational::from_f64_approx(rng.gen_range(-rad..rad)).approximate(&Natural::from(3u64))
        })
    }

    // let pt1 = Vector::new(space.clone(), vec![Rational::from(0), Rational::from(0)]);
    // let pt2 = Vector::new(space.clone(), vec![Rational::from(0), Rational::from(-1)]);
    // let pt3 = Vector::new(space.clone(), vec![Rational::from(0), Rational::from(1)]);
    // let pt4 = Vector::new(space.clone(), vec![Rational::from(1), Rational::from(0)]);

    // let spx1 = Simplex::new(space.clone(), vec![pt1]).unwrap();
    // let spx2 = Simplex::new(space.clone(), vec![pt2, pt3, pt4]).unwrap();

    // let VennResult {
    //     left: a,
    //     middle: b,
    //     right: c,
    // } = spx1.venn(&spx2);

    let ch1 = ConvexHull::new(
        space.clone(),
        (0..6)
            .map(|i| random_point(&space, (i + 1) as f64))
            .collect(),
    );
    let ch2 = ConvexHull::new(
        space.clone(),
        (0..6)
            .map(|i| random_point(&space, (i + 1) as f64))
            .collect(),
    );
    // let ch3 = ch1.intersect(&ch2);

    let sc1 = ch1.as_simplicial_complex().into_forget_labels();
    let sc2 = ch2.as_simplicial_complex().into_forget_labels();
    // let sc3 = ch3.as_simplicial_complex().entire;

    // let VennResult {
    //     left: a,
    //     middle: b,
    //     right: c,0

    // }

    let sc4 = LabelledSimplicialDisjointUnion::union_raw(&(&sc1).into(), &(&sc2).into());
    println!("done union");
    let sc5 = sc4.clone().refine_to_partial_simplicial_complex().closure();
    println!("done to sc");
    let sc6 = sc5.clone().simplify();
    println!("done simplify");

    let mut canvas = Canvas2D::new(Box::new(MouseWheelZoomCamera::new()));
    canvas.plot_shapes(
        [Shape::SetThickness(0.3)]
            .into_iter()
            .chain(simplicial_complex_shapes(
                &Colour::magenta(),
                &Colour::magenta().darken(),
                0.6,
                &sc1,
            ))
            .chain(simplicial_complex_shapes(
                &Colour::cyan(),
                &Colour::cyan().darken(),
                0.6,
                &sc2,
            ))
            .chain(simplicial_complex_shapes(
                &Colour::green(),
                &Colour::green().darken(),
                0.3,
                &sc6,
            )),
    );
    canvas.run();
}
