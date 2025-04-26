use crate::canvas::{BasicCanvasState, Canvas, CanvasState};

pub trait Canvas2DState: CanvasState {
    fn camera_uniform(&self) -> CameraUniform;
}

pub trait Canvas2D: Canvas
where
    Self::State: Canvas2DState,
{
}

#[repr(C)]
#[derive(Debug, Copy, Clone, bytemuck::Pod, bytemuck::Zeroable)]
pub struct CameraUniform {
    mat: [[f32; 2]; 2],
    offset: [f32; 2],
}

struct BasicCamera {
    // the x coordinate at the centre of the screen
    mid_x: f32,
    // the y coordinate at the centre of the screen
    mid_y: f32,
    // the square root of the visible area
    sqrt_area: f32,
}

impl BasicCamera {
    fn camera_uniform(&self, display_size: winit::dpi::PhysicalSize<u32>) -> CameraUniform {
        let display_size = (display_size.width as f32, display_size.height as f32);
        let avg_side = (display_size.0 * display_size.1).sqrt();
        let x_mult = 2.0 * avg_side / (display_size.0 * self.sqrt_area);
        let y_mult = 2.0 * avg_side / (display_size.1 * self.sqrt_area);
        CameraUniform {
            mat: [[x_mult, 0.0], [0.0, y_mult]],
            offset: [-self.mid_x * x_mult, -self.mid_y * y_mult],
        }
    }
}

pub struct BasicCanvas2DState {
    state: BasicCanvasState,
    camera: BasicCamera,
}

impl CanvasState for BasicCanvas2DState {
    fn new(window: std::sync::Arc<winit::window::Window>) -> Self {
        let state: BasicCanvasState = CanvasState::new(window);

        let camera = BasicCamera {
            mid_x: 0.0,
            mid_y: 0.0,
            sqrt_area: 2.7,
        };

        Self { state, camera }
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
    fn camera_uniform(&self) -> CameraUniform {
        self.camera.camera_uniform(self.size())
    }
}
