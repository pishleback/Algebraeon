use algebraeon_drawing::{canvas::Canvas, canvas2d::*};
use algebraeon_nzq::Integer;
use algebraeon_rings::{
    polynomial::{Polynomial, PolynomialFromStr}, structure::MetaComplexSubset,
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
    // canvas.plot_test_pentagon();
    let mut points = canvas.points();
    for root in p.all_complex_roots() {
        let (x, y) = root.as_f32_real_and_imaginary_parts();
        println!("{x} {y}");
        points = points.add(x, y);
    }
    points.plot();

    canvas.run();
}
