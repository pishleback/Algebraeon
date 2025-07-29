use algebraeon_drawing::{
    canvas::Canvas,
    canvas2d::{points::Point, rectangles::Rectangle, *},
};
use algebraeon_nzq::Integer;
use algebraeon_rings::{
    polynomial::{Polynomial, PolynomialFromStr},
    structure::{MetaComplexSubset, MetaRealToFloat},
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

    canvas.plot_points(p.all_complex_roots().into_iter().map(|root| {
        let (x, y) = root.as_f32_real_and_imaginary_parts();
        Point { pos: [x, y] }
    }));

    canvas.plot_rectangles(p.all_complex_roots().into_iter().filter_map(
        |root| match root.isolate() {
            algebraeon_rings::isolated_algebraic::ComplexIsolatingRegion::Rational(r) => None,
            algebraeon_rings::isolated_algebraic::ComplexIsolatingRegion::RealInterval(a, b) => {
                None
            }
            algebraeon_rings::isolated_algebraic::ComplexIsolatingRegion::Box(a, b, c, d) => {
                Some(Rectangle {
                    a: a.as_f32(),
                    b: b.as_f32(),
                    c: c.as_f32(),
                    d: d.as_f32(),
                })
            }
        },
    ));
    canvas.run();
}
