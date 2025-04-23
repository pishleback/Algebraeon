use algebraeon_drawing::canvas2d::*;

struct TestElement {}

impl TestElement {
    pub fn new() -> Self {
        Self {}
    }
}

impl EventHandler for TestElement {
    fn tick(&mut self, dt: f64) {}
    fn event(&mut self, ev: &glium::glutin::event::Event<'_, ()>) {}
}

impl DrawElement for TestElement {
    fn draw(&mut self, display: &glium::Display) {}
}

fn main() {
    let mut canvas = Canvas2D::new(Box::new(MouseWheelZoomCamera::new()));

    canvas.add_element(Box::new(TestElement::new()));

    canvas.run(1000, 600);
}
