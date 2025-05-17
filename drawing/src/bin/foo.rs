use algebraeon_drawing::{canvas::Canvas, canvas2d::{pentagon::{Pentagon, PentagonWgpu}, *}};

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

    canvas.add_item(Pentagon::new());

    // canvas.add_element(Box::new(TestElement::new()));

    canvas.run();
}
