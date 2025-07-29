use algebraeon_drawing::{
    canvas::Canvas,
    canvas2d::{points::Point, rectangles::Rectangle, *},
};
use algebraeon_nzq::Integer;
use algebraeon_rings::{
    polynomial::{Polynomial, PolynomialFromStr},
    structure::MetaRealToFloat,
};

fn main() {
    simplelog::CombinedLogger::init(vec![simplelog::TermLogger::new(
        simplelog::LevelFilter::Trace,
        simplelog::Config::default(),
        simplelog::TerminalMode::Mixed,
        simplelog::ColorChoice::Auto,
    )])
    .unwrap();
    log::set_max_level(log::LevelFilter::Warn);

    let p = Polynomial::<Integer>::from_str("x^5 - x + 1", "x").unwrap();

    let mut canvas = Canvas2D::new(Box::new(MouseWheelZoomCamera::new()));

    canvas.plot_complex_polynomial(p.clone());

    let roots = p.all_complex_roots();

    canvas.plot_points(roots.iter().filter_map(|root| match root.isolate() {
        algebraeon_rings::isolated_algebraic::ComplexIsolatingRegion::Rational(r) => Some(Point {
            x: r.as_f32(),
            y: 0.0,
        }),
        algebraeon_rings::isolated_algebraic::ComplexIsolatingRegion::RealInterval(..) => None,
        algebraeon_rings::isolated_algebraic::ComplexIsolatingRegion::Box(..) => None,
    }));

    canvas.plot_rectangles(roots.iter().filter_map(|root| match root.isolate() {
        algebraeon_rings::isolated_algebraic::ComplexIsolatingRegion::Rational(..) => None,
        algebraeon_rings::isolated_algebraic::ComplexIsolatingRegion::RealInterval(a, b) => {
            Some(Rectangle {
                a: a.as_f32(),
                b: b.as_f32(),
                c: -0.02,
                d: 0.02,
            })
        }
        algebraeon_rings::isolated_algebraic::ComplexIsolatingRegion::Box(a, b, c, d) => {
            Some(Rectangle {
                a: a.as_f32(),
                b: b.as_f32(),
                c: c.as_f32(),
                d: d.as_f32(),
            })
        }
    }));
    canvas.run();
}
