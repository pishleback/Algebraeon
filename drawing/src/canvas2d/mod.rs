use glium::{
    Display, Frame, Surface,
    glutin::surface::WindowSurface,
    winit::{
        event::{Event, WindowEvent},
        event_loop::{EventLoop, EventLoopBuilder},
    },
};
use std::time::{Duration, Instant};

pub trait EventHandler {
    fn tick(&mut self, dt: &Duration);
    fn event(&mut self, ev: &Event<()>);
}

pub trait Camera: EventHandler {
    fn view_matrix_and_shift(&self, display_size: (u32, u32)) -> ([[f64; 2]; 2], [f64; 2]);
}

pub trait DrawElement: EventHandler {
    fn draw(
        &mut self,
        display: &Display<WindowSurface>,
        target: &mut Frame,
        display_size: (u32, u32),
        camera: &Box<dyn Camera>,
    );
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
            mid_x: 1.0,
            mid_y: 1.0,
            sqrt_area: 6.0,
        }
    }
}

impl EventHandler for MouseWheelZoomCamera {
    fn tick(&mut self, dt: &Duration) {}
    fn event(&mut self, ev: &Event<()>) {}
}

impl Camera for MouseWheelZoomCamera {
    fn view_matrix_and_shift(&self, display_size: (u32, u32)) -> ([[f64; 2]; 2], [f64; 2]) {
        let display_size = (display_size.0 as f64, display_size.1 as f64);
        let avg_side = (display_size.0 * display_size.1).sqrt();
        let x_mult = 2.0 * avg_side / (display_size.0 * self.sqrt_area);
        let y_mult = 2.0 * avg_side / (display_size.1 * self.sqrt_area);
        (
            [[x_mult, 0.0], [0.0, y_mult]],
            [-self.mid_x * x_mult, -self.mid_y * y_mult],
        )
    }
}

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

    fn tick(&mut self, dt: &Duration) {
        for element in &mut self.elements {
            element.tick(dt);
        }
    }

    fn draw(&mut self, display: &Display<WindowSurface>, display_size: (u32, u32)) {
        let mut target = display.draw();
        target.clear_color(0.0, 0.0, 0.0, 1.0);
        for element in &mut self.elements {
            element.draw(display, &mut target, display_size, &self.camera);
        }
        target.finish().unwrap();
    }

    fn event(&mut self, ev: &Event<()>) {
        for element in &mut self.elements {
            element.event(ev);
        }
    }

    pub fn run(mut self, display_width: u32, display_height: u32) {
        let event_loop = glium::winit::event_loop::EventLoop::builder()
            .build()
            .expect("event loop building");
        let (_window, display) = glium::backend::glutin::SimpleWindowBuilder::new()
            .with_inner_size(display_width, display_height)
            .with_title("Hello World")
            .build(&event_loop);

        let mut display_size = (display_width, display_height);

        let mut since_tick = Duration::from_secs(0);
        let mut prev_time = Instant::now();
        let tick_interval = Duration::from_secs(1 / 30);

        #[allow(deprecated)]
        event_loop
            .run(move |ev, window_target| {
                let time = Instant::now();
                let dt = time - prev_time;
                prev_time = time;
                since_tick += dt;

                match &ev {
                    glium::winit::event::Event::WindowEvent { event, .. } => match event {
                        glium::winit::event::WindowEvent::CloseRequested => {
                            window_target.exit();
                        }
                        // glium::glutin::event::WindowEvent::CursorMoved { position, .. } => {
                        //     state.mouse_pos = (position.x, position.y);
                        // }
                        glium::winit::event::WindowEvent::Resized(size) => {
                            display_size = (size.width, size.height);
                        }
                        _ => {}
                    },
                    glium::winit::event::Event::AboutToWait => {
                        if since_tick >= tick_interval {
                            self.tick(&dt);
                            self.draw(&display, display_size);
                            since_tick = Duration::from_secs(0);
                            window_target.set_control_flow(
                                glium::winit::event_loop::ControlFlow::WaitUntil(
                                    time + tick_interval,
                                ),
                            );
                        } else {
                            window_target.set_control_flow(
                                glium::winit::event_loop::ControlFlow::WaitUntil(
                                    time + tick_interval - since_tick,
                                ),
                            );
                        }
                    }
                    _ => {}
                }
                self.event(&ev);
            })
            .unwrap();
    }
}
