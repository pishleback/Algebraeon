use std::{cell::OnceCell, sync::Arc};
use winit::{
    application::ApplicationHandler,
    error::EventLoopError,
    event::WindowEvent,
    event_loop::{ActiveEventLoop, ControlFlow, EventLoop},
    window::{Window, WindowId},
};

pub trait Canvas: Sized {
    type WindowState;
    fn new_state(&self, window: Arc<Window>) -> Self::WindowState;
    fn run(self) -> Result<(), EventLoopError> {
        let event_loop = EventLoop::new().unwrap();
        event_loop.set_control_flow(ControlFlow::Poll);
        let mut app_handler = CanvasApplicationHandler::new(self);
        event_loop.run_app(&mut app_handler)
    }
    fn window_event(
        &mut self,
        window_state: &mut Self::WindowState,
        event_loop: &ActiveEventLoop,
        _id: WindowId,
        event: &WindowEvent,
    );
}

struct CanvasApplicationHandler<C: Canvas> {
    canvas: C,
    window: OnceCell<Arc<Window>>,
    state: Option<C::WindowState>,
}

impl<C: Canvas> CanvasApplicationHandler<C> {
    pub fn new(canvas: C) -> Self {
        Self {
            canvas,
            window: OnceCell::default(),
            state: None,
        }
    }
}

impl<C: Canvas> ApplicationHandler for CanvasApplicationHandler<C> {
    fn resumed(&mut self, event_loop: &ActiveEventLoop) {
        let window = self
            .window
            .get_or_init(|| {
                Arc::new(
                    event_loop
                        .create_window(Window::default_attributes())
                        .unwrap(),
                )
            })
            .clone();
        self.state = Some(self.canvas.new_state(window));
    }

    fn window_event(&mut self, event_loop: &ActiveEventLoop, id: WindowId, event: WindowEvent) {
        match event {
            WindowEvent::CloseRequested => {
                event_loop.exit();
            }
            _ => (),
        }
        let window_state = self.state.as_mut().unwrap();
        self.canvas
            .window_event(window_state, event_loop, id, &event);
    }
}
