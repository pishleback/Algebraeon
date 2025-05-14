use std::sync::Arc;

use winit::event::WindowEvent;

use crate::{canvas::Canvas, wgpu_state::WgpuState};

pub struct ExampleCanvas {}

pub struct Example1WindowState {
    wgpu_state: WgpuState,
}

impl ExampleCanvas {
    pub fn new() -> Self {
        Self {}
    }

    fn render(&mut self, state: &mut WgpuState) -> Result<(), wgpu::SurfaceError> {
        let output = state.surface.get_current_texture()?;
        let view = output
            .texture
            .create_view(&wgpu::TextureViewDescriptor::default());
        let mut encoder = state
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("Render Encoder"),
            });

        state.queue.submit(std::iter::once(encoder.finish()));
        output.present();
        Ok(())
    }
}

impl Canvas for ExampleCanvas {
    type WindowState = Example1WindowState;

    fn new_state(&self, window: Arc<winit::window::Window>) -> Self::WindowState {
        Example1WindowState {
            wgpu_state: WgpuState::new(window),
        }
    }

    fn window_event(
        &mut self,
        window_state: &mut Self::WindowState,
        event_loop: &winit::event_loop::ActiveEventLoop,
        _id: winit::window::WindowId,
        event: &WindowEvent,
    ) {
        match event {
            WindowEvent::Resized(physical_size) => {
                window_state.wgpu_state.resize(*physical_size);
            }
            WindowEvent::RedrawRequested => {
                match self.render(&mut window_state.wgpu_state) {
                    Ok(()) => {}
                    // Reconfigure the surface if it's lost or outdated
                    Err(wgpu::SurfaceError::Lost | wgpu::SurfaceError::Outdated) => {
                        window_state.wgpu_state.reconfigure();
                    }
                    // The system is out of memory, we should probably quit
                    Err(wgpu::SurfaceError::OutOfMemory | wgpu::SurfaceError::Other) => {
                        log::error!("OutOfMemory");
                        event_loop.exit();
                    }

                    // This happens when the a frame takes too long to present
                    Err(wgpu::SurfaceError::Timeout) => {
                        log::warn!("Surface timeout")
                    }
                }

                window_state.wgpu_state.window.request_redraw();
            }
            _ => {}
        }
    }
}
