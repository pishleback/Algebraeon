use glium::{Display, backend::Facade, glutin::event::Event};
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
    fn draw(&mut self, state: &State, display: &Display);
    fn event(&mut self, state: &State, ev: &Event<'_, ()>);
    fn run(init: impl FnOnce(&mut Self)) -> !
    where
        Self: Sized + 'static,
    {
        let display_size = (1 * 1024, 1 * 768);

        // 1. We start by creating the EventLoop, this can only be done once per process.
        // This also needs to happen on the main thread to make the program portable.
        let event_loop = glium::glutin::event_loop::EventLoopBuilder::new().build();

        // 2. Parameters for building the Window.
        let wb = glium::glutin::window::WindowBuilder::new()
            .with_inner_size(glium::glutin::dpi::LogicalSize::new(
                display_size.0 as f64,
                display_size.1 as f64,
            ))
            .with_title("Hello world");

        // 3. Parameters for building the OpenGL context.
        let cb = glium::glutin::ContextBuilder::new();

        // 4. Build the Display with the given window and OpenGL context parameters and register the
        //    window with the events_loop.
        let display = glium::Display::new(wb, cb, &event_loop).unwrap();

        let mut canvas = Self::new(&display);
        init(&mut canvas);

        let mut state = State {
            mouse_pos: (0.0, 0.0),
            display_size,
        };

        let mut prev_time = Instant::now();

        event_loop.run(move |ev, _, control_flow| {
            let time = Instant::now();
            let dt = (time - prev_time).as_secs_f64();
            prev_time = time;

            let mut stop = false;

            // println!("{:?}", ev);

            //events
            match &ev {
                glium::glutin::event::Event::WindowEvent { event, .. } => match event {
                    glium::glutin::event::WindowEvent::CloseRequested => {
                        stop = true;
                    }
                    glium::glutin::event::WindowEvent::CursorMoved { position, .. } => {
                        state.mouse_pos = (position.x, position.y);
                    }
                    glium::glutin::event::WindowEvent::Resized(size) => {
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

            //control flow
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
