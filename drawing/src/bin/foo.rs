use algebraeon_drawing::{
    canvas::Canvas,
    canvas2d::{complex_polynomial::PolynomialPlot, *},
};

fn main() {
    simplelog::CombinedLogger::init(vec![simplelog::TermLogger::new(
        simplelog::LevelFilter::Trace,
        simplelog::Config::default(),
        simplelog::TerminalMode::Mixed,
        simplelog::ColorChoice::Auto,
    )])
    .unwrap();
    log::set_max_level(log::LevelFilter::Debug);

    let mut canvas = Canvas2D::new(Box::new(MouseWheelZoomCamera::new()));

    canvas.add_item(PolynomialPlot::new());
    // canvas.add_item(Pentagon::new());

    canvas.run();
}
