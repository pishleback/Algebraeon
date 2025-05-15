use algebraeon_drawing::{canvas::Canvas, canvas2d::*};

fn main() {
    simplelog::CombinedLogger::init(vec![simplelog::TermLogger::new(
        simplelog::LevelFilter::Trace,
        simplelog::Config::default(),
        simplelog::TerminalMode::Mixed,
        simplelog::ColorChoice::Auto,
    )])
    .unwrap();
    log::set_max_level(log::LevelFilter::Debug);

    let canvas = Canvas2D::new(Box::new(MouseWheelZoomCamera::new()));

    // canvas.add_element(Box::new(TestElement::new()));

    canvas.run();
}
