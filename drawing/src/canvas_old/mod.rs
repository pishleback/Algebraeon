use glium::{Display, backend::Facade, glutin::surface::WindowSurface, winit::event::Event};
use std::time::Instant;

pub mod canvas2d;
pub mod simplexes;

#[derive(Debug)]
pub struct State {
    mouse_pos: (f64, f64),
    display_size: (u32, u32),
}

pub trait Canvas {
    fn new(facade: &impl Facade) -> Self;
    fn tick(&mut self, state: &State, dt: f64);
    fn draw(&mut self, state: &State, display: &Display<WindowSurface>);
    fn event(&mut self, state: &State, ev: &Event<()>);
    fn run(init: impl FnOnce(&mut Self)) -> ()
    where
        Self: Sized + 'static,
    {
        let display_size = (1 * 1024, 1 * 768);

        let event_loop = glium::winit::event_loop::EventLoop::builder()
            .build()
            .expect("event loop building");
        let (_window, display) = glium::backend::glutin::SimpleWindowBuilder::new()
            .with_inner_size(1000, 600)
            .with_title("Hello World")
            .build(&event_loop);

        let mut canvas = Self::new(&display);
        init(&mut canvas);

        let mut state = State {
            mouse_pos: (0.0, 0.0),
            display_size,
        };

        let mut prev_time = Instant::now();

        #[allow(deprecated)]
        event_loop.run(move |ev, window_target| {
            let time = Instant::now();
            let dt = (time - prev_time).as_secs_f64();
            prev_time = time;

            let mut stop = false;

            //events
            match &ev {
                glium::winit::event::Event::WindowEvent { event, .. } => match event {
                    glium::winit::event::WindowEvent::CloseRequested => {
                        window_target.exit();
                    }
                    glium::winit::event::WindowEvent::CursorMoved { position, .. } => {
                        state.mouse_pos = (position.x, position.y);
                    }
                    glium::winit::event::WindowEvent::Resized(size) => {
                        state.display_size = (size.width, size.height);
                    }
                    _ => {}
                },
                _ => {}
            }
            canvas.event(&state, &ev);

            //update
            canvas.tick(&state, dt);

            //draw
            canvas.draw(&state, &display);
        }).unwrap();
    }
}
