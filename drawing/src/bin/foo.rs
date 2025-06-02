use algebraeon_drawing::{
    canvas::Canvas,
    canvas2d::{complex_polynomial::PolynomialPlot, *},
};
use algebraeon_nzq::Integer;
use algebraeon_rings::{polynomial::Polynomial, structure::IntoErgonomic};

fn main() {
    simplelog::CombinedLogger::init(vec![simplelog::TermLogger::new(
        simplelog::LevelFilter::Trace,
        simplelog::Config::default(),
        simplelog::TerminalMode::Mixed,
        simplelog::ColorChoice::Auto,
    )])
    .unwrap();
    log::set_max_level(log::LevelFilter::Warn);

    let x = Polynomial::<Integer>::var().into_ergonomic();
    let p = (x.pow(5) - x + 1).into_verbose();

    let mut canvas = Canvas2D::new(Box::new(MouseWheelZoomCamera::new()));

    canvas.add_item(PolynomialPlot::new(p));
    // canvas.add_item(Pentagon::new());

    canvas.run();
}
