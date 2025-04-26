use cgmath::{Matrix2, Vector2};

use crate::canvas::{BasicCanvasState, Canvas, CanvasState};

pub trait Canvas2DState: CanvasState {
    fn camera_mat_and_slide() -> (Matrix2<f64>, Vector2<f64>);
}

pub trait Canvas2D: Canvas
where
    Self::State: Canvas2DState,
{
}

pub struct BasicCanvas2DState {
    state: BasicCanvasState,
}

impl CanvasState for BasicCanvas2DState {
    fn new(window: std::sync::Arc<winit::window::Window>) -> Self {
        Self {
            state: CanvasState::new(window),
        }
    }

    fn size(&self) -> winit::dpi::PhysicalSize<u32> {
        self.state.size()
    }

    fn device(&self) -> &wgpu::Device {
        self.state.device()
    }

    fn config(&self) -> &wgpu::SurfaceConfiguration {
        self.state.config()
    }

    fn resize(&mut self, new_size: winit::dpi::PhysicalSize<u32>) {
        self.state.resize(new_size)
    }

    fn render(
        &mut self,
        renderers: &Vec<Box<dyn crate::canvas::RendererInstance>>,
    ) -> Result<(), wgpu::SurfaceError> {
        self.state.render(renderers)
    }
}

impl Canvas2DState for BasicCanvas2DState {
    fn camera_mat_and_slide() -> (Matrix2<f64>, Vector2<f64>) {
        (Matrix2::new(1.0, 0.0, 0.0, 1.0), Vector2::new(0.0, 0.0))
    }
}
