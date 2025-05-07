use std::{
    cell::OnceCell,
    sync::Arc,
    time::{Duration, Instant},
};
use wgpu::util::DeviceExt;
use winit::{
    application::ApplicationHandler,
    event::{Event, WindowEvent},
    event_loop::{ActiveEventLoop, ControlFlow, EventLoop},
    window::{Window, WindowId},
};

pub trait EventHandler {
    fn tick(&mut self, dt: &Duration);
    fn event(&mut self, ev: &Event<()>);
}

pub trait Camera: EventHandler {
    fn view_matrix_and_shift(&self, display_size: (u32, u32)) -> ([[f64; 2]; 2], [f64; 2]);
}

// pub trait DrawElement: EventHandler {
//     fn draw(
//         &mut self,
//         display: &Display<WindowSurface>,
//         target: &mut Frame,
//         display_size: (u32, u32),
//         camera: &Box<dyn Camera>,
//     );
// }

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
    // elements: Vec<Box<dyn DrawElement>>,
}

impl Canvas2D {
    pub fn new(camera: Box<dyn Camera>) -> Self {
        Self {
            camera,
            // elements: vec![],
        }
    }

    // pub fn add_element(&mut self, element: Box<dyn DrawElement>) {
    //     self.elements.push(element);
    // }

    fn tick(&mut self, dt: &Duration) {
        // for element in &mut self.elements {
        //     element.tick(dt);
        // }
    }

    // fn draw(&mut self, display: &Display<WindowSurface>, display_size: (u32, u32)) {
    //     let mut target = display.draw();
    //     target.clear_color(0.0, 0.0, 0.0, 1.0);
    //     for element in &mut self.elements {
    //         element.draw(display, &mut target, display_size, &self.camera);
    //     }
    //     target.finish().unwrap();
    // }

    fn event(&mut self, ev: &Event<()>) {
        // for element in &mut self.elements {
        //     element.event(ev);
        // }
    }

    pub fn run(mut self, display_width: u32, display_height: u32) {
        let event_loop = EventLoop::new().unwrap();
        event_loop.set_control_flow(ControlFlow::Poll);
        let mut app = Canvas2DApp {
            canvas: self,
            window: OnceCell::default(),
            window_state: None,
        };
        event_loop.run_app(&mut app).unwrap();
    }
}
