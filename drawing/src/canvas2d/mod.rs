use glium::{Display, glutin::event::Event};
use std::time::Instant;

pub trait EventHandler {
    fn tick(&mut self, dt: f64);
    fn event(&mut self, ev: &Event<'_, ()>);
}

pub trait Camera: EventHandler {}

pub trait DrawElement: EventHandler {
    fn draw(&mut self, display: &Display);
}

pub struct MouseWheelZoomCamera {
    // the x coordinate at the centre of the screen
    mid_x: f64,
    // the y coordinate at the centre of the screen
    mid_y: f64,
    // the square root of the visible area
    sqrt_area: f64,
}

impl MouseWheelZoomCamera {
    pub fn new() -> Self {
        Self {
            mid_x: 0.0,
            mid_y: 0.0,
            sqrt_area: 2.0,
        }
    }
}

impl EventHandler for MouseWheelZoomCamera {
    fn tick(&mut self, dt: f64) {}

    fn event(&mut self, ev: &Event<'_, ()>) {}
}

impl Camera for MouseWheelZoomCamera {}

pub struct Canvas2D {
    camera: Box<dyn Camera>,
    elements: Vec<Box<dyn DrawElement>>,
}

impl Canvas2D {
    pub fn new(camera: Box<dyn Camera>) -> Self {
        Self {
            camera,
            elements: vec![],
        }
    }

    pub fn add_element(&mut self, element: Box<dyn DrawElement>) {
        self.elements.push(element);
    }

    fn tick(&mut self, dt: f64) {
        for element in &mut self.elements {
            element.tick(dt);
        }
    }

    fn draw(&mut self, display: &Display) {
        for element in &mut self.elements {
            element.draw(display);
        }
    }

    fn event(&mut self, ev: &Event<'_, ()>) {
        for element in &mut self.elements {
            element.event(ev);
        }
    }

    pub fn run(mut self, display_width: i32, display_height: i32) {
        let event_loop = glium::glutin::event_loop::EventLoopBuilder::new().build();

        let wb = glium::glutin::window::WindowBuilder::new()
            .with_inner_size(glium::glutin::dpi::LogicalSize::new(
                display_width as f64,
                display_height as f64,
            ))
            .with_title("Hello world");

        let cb = glium::glutin::ContextBuilder::new();

        let display = glium::Display::new(wb, cb, &event_loop).unwrap();

        let mut prev_time = Instant::now();

        event_loop.run(move |ev, _, control_flow| {
            let time = Instant::now();
            let dt = (time - prev_time).as_secs_f64();
            prev_time = time;

            let mut stop = false;

            match &ev {
                glium::glutin::event::Event::WindowEvent { event, .. } => match event {
                    glium::glutin::event::WindowEvent::CloseRequested => {
                        stop = true;
                    }
                    // glium::glutin::event::WindowEvent::CursorMoved { position, .. } => {
                    //     state.mouse_pos = (position.x, position.y);
                    // }
                    // glium::glutin::event::WindowEvent::Resized(size) => {
                    //     state.display_size = (size.width, size.height);
                    // }
                    _ => {}
                },
                _ => {}
            }
            self.event(&ev);
            self.tick(dt);
            self.draw(&display);

            match stop {
                true => {
                    *control_flow = glium::glutin::event_loop::ControlFlow::Exit;
                }
                false => {
                    *control_flow = glium::glutin::event_loop::ControlFlow::Poll;
                }
            }
        });
    }
}
