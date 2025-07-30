use crate::canvas::*;
use std::sync::Arc;
use wgpu::{BindGroup, BindGroupLayout, CommandEncoder, TextureView, util::DeviceExt};
use winit::{
    dpi::{PhysicalPosition, PhysicalSize},
    event::WindowEvent,
    event_loop::ActiveEventLoop,
    window::{Window, WindowId},
};

pub mod complex_polynomial;
pub mod shapes;
pub mod test_pentagon;

#[repr(C)]
#[derive(Debug, Copy, Clone, bytemuck::Pod, bytemuck::Zeroable)]
pub struct CameraUniform {
    pub matrix: [[f32; 2]; 2],
    pub matrix_inv: [[f32; 2]; 2],
    pub shift: [f32; 2],
}

fn mat2x2inv([[a, b], [c, d]]: [[f32; 2]; 2]) -> [[f32; 2]; 2] {
    let det = a * d - b * c;
    [[d / det, -c / det], [-b / det, a / det]]
}

pub trait Camera {
    fn get_uniform(&self, display_size: PhysicalSize<u32>) -> CameraUniform;

    fn window_event(
        &mut self,
        display_size: PhysicalSize<u32>,
        mouse_pos: PhysicalPosition<f64>,
        event_loop: &ActiveEventLoop,
        id: WindowId,
        event: &WindowEvent,
    );

    fn pixel_to_wgpu(
        &self,
        display_size: PhysicalSize<u32>,
        pixels: PhysicalPosition<f64>,
    ) -> (f64, f64) {
        (
            2.0 * pixels.x / display_size.width as f64 - 1.0,
            -2.0 * pixels.y / display_size.height as f64 + 1.0,
        )
    }

    fn wgpu_to_pixel(
        &self,
        display_size: PhysicalSize<u32>,
        wgpu: (f64, f64),
    ) -> PhysicalPosition<f64> {
        PhysicalPosition::new(
            display_size.width as f64 * (wgpu.0 + 1.0) / 2.0,
            display_size.height as f64 * (-wgpu.1 + 1.0) / 2.0,
        )
    }

    fn wgpu_to_coord(&self, display_size: PhysicalSize<u32>, wgpu: (f64, f64)) -> (f64, f64) {
        let uniform = self.get_uniform(display_size);
        let [[a, b], [c, d]] = mat2x2inv(uniform.matrix);
        let v = (
            wgpu.0 - uniform.shift[0] as f64,
            wgpu.1 - uniform.shift[1] as f64,
        );
        (
            a as f64 * v.0 + b as f64 * v.1,
            c as f64 * v.0 + d as f64 * v.1,
        )
    }

    fn coord_to_wgpu(&self, display_size: PhysicalSize<u32>, coord: (f64, f64)) -> (f64, f64) {
        let uniform = self.get_uniform(display_size);
        let [[a, b], [c, d]] = uniform.matrix;
        let [x, y] = uniform.shift;
        let v = coord;
        (
            a as f64 * v.0 + b as f64 * v.1 + x as f64,
            c as f64 * v.0 + d as f64 * v.1 + y as f64,
        )
    }

    fn pixel_to_coord(
        &self,
        display_size: PhysicalSize<u32>,
        pixels: PhysicalPosition<f64>,
    ) -> (f64, f64) {
        self.wgpu_to_coord(display_size, self.pixel_to_wgpu(display_size, pixels))
    }

    fn coord_to_pixel(
        &self,
        display_size: PhysicalSize<u32>,
        coord: (f64, f64),
    ) -> PhysicalPosition<f64> {
        self.wgpu_to_pixel(display_size, self.coord_to_wgpu(display_size, coord))
    }
}

pub struct MouseWheelZoomCamera {
    // the x coordinate at the centre of the screen
    mid_x: f64,
    // the y coordinate at the centre of the screen
    mid_y: f64,
    // the square root of the visible area
    sqrt_area: f64,
}

#[allow(clippy::new_without_default)]
impl MouseWheelZoomCamera {
    pub fn new() -> Self {
        Self {
            mid_x: 0.0,
            mid_y: 0.0,
            sqrt_area: 6.0,
        }
    }
}

impl Camera for MouseWheelZoomCamera {
    fn get_uniform(&self, display_size: PhysicalSize<u32>) -> CameraUniform {
        let display_size = (display_size.width as f32, display_size.height as f32);
        let avg_side = (display_size.0 * display_size.1).sqrt();
        let x_mult = 2.0 * avg_side / (display_size.0 * self.sqrt_area as f32);
        let y_mult = 2.0 * avg_side / (display_size.1 * self.sqrt_area as f32);
        let matrix = [[x_mult, 0.0], [0.0, y_mult]];
        let matrix_inv = mat2x2inv(matrix);
        CameraUniform {
            matrix,
            matrix_inv,
            shift: [-self.mid_x as f32 * x_mult, -self.mid_y as f32 * y_mult],
        }
    }

    fn window_event(
        &mut self,
        display_size: PhysicalSize<u32>,
        mouse_pos: PhysicalPosition<f64>,
        _event_loop: &ActiveEventLoop,
        _id: WindowId,
        event: &WindowEvent,
    ) {
        #[allow(clippy::single_match)]
        match event {
            WindowEvent::MouseWheel { delta, .. } => {
                let dy = match delta {
                    winit::event::MouseScrollDelta::PixelDelta(pos) => pos.y,
                    winit::event::MouseScrollDelta::LineDelta(_x, y) => *y as f64,
                };
                self.zoom_event(display_size, mouse_pos, 0.8f64.powf(dy));
            }
            _ => {}
        }
    }
}

impl MouseWheelZoomCamera {
    fn zoom_event(
        &mut self,
        display_size: PhysicalSize<u32>,
        center: PhysicalPosition<f64>,
        mult: f64,
    ) {
        // println!("{:?} {:?}", center, mult);
        let center_before = self.pixel_to_coord(display_size, center);
        self.sqrt_area *= mult;
        let center_after = self.pixel_to_coord(display_size, center);
        self.mid_x += center_before.0 - center_after.0;
        self.mid_y += center_before.1 - center_after.1;
    }
}

pub struct Canvas2D {
    mouse_pos: PhysicalPosition<f64>,
    items: Vec<Box<dyn Canvas2DItem>>,
    camera: Box<dyn Camera>,
}

pub trait Canvas2DItemWgpu {
    fn render(
        &mut self,
        encoder: &mut CommandEncoder,
        view: &TextureView,
        camera_bind_group: &BindGroup,
    ) -> Result<(), wgpu::SurfaceError>;
}

#[allow(clippy::new_without_default)]
pub trait Canvas2DItem {
    fn new_wgpu(
        &self,
        wgpu_state: &WgpuState,
        camera_bind_group_layout: &BindGroupLayout,
    ) -> Box<dyn Canvas2DItemWgpu>;
}

pub struct Canvas2DWindowState {
    wgpu_state: WgpuState,

    camera_uniform: CameraUniform,
    camera_buffer: wgpu::Buffer,
    camera_bind_group_layout: wgpu::BindGroupLayout,
    camera_bind_group: wgpu::BindGroup,

    items: Vec<Box<dyn Canvas2DItemWgpu>>,
}

impl Canvas2DWindowState {
    fn new(window: Arc<Window>) -> Self {
        let wgpu_state = WgpuState::new(window);

        let camera_uniform = CameraUniform {
            matrix: [[1.0, 0.0], [0.0, 1.0]],
            matrix_inv: [[1.0, 0.0], [0.0, 1.0]],
            shift: [0.0, 0.0],
        };

        let camera_buffer =
            wgpu_state
                .device
                .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some("Camera Buffer"),
                    contents: bytemuck::cast_slice(&[camera_uniform]),
                    usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
                });

        let camera_bind_group_layout =
            wgpu_state
                .device
                .create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                    entries: &[wgpu::BindGroupLayoutEntry {
                        binding: 0,
                        visibility: wgpu::ShaderStages::VERTEX,
                        ty: wgpu::BindingType::Buffer {
                            ty: wgpu::BufferBindingType::Uniform,
                            has_dynamic_offset: false,
                            min_binding_size: None,
                        },
                        count: None,
                    }],
                    label: Some("camera_bind_group_layout"),
                });

        let camera_bind_group = wgpu_state
            .device
            .create_bind_group(&wgpu::BindGroupDescriptor {
                layout: &camera_bind_group_layout,
                entries: &[wgpu::BindGroupEntry {
                    binding: 0,
                    resource: camera_buffer.as_entire_binding(),
                }],
                label: Some("camera_bind_group"),
            });

        Self {
            items: vec![],
            wgpu_state,
            camera_uniform,
            camera_bind_group,
            camera_buffer,
            camera_bind_group_layout,
        }
    }

    fn render(&mut self) -> Result<(), wgpu::SurfaceError> {
        self.wgpu_state.queue.write_buffer(
            &self.camera_buffer,
            0,
            bytemuck::cast_slice(&[self.camera_uniform]),
        );

        let output = self.wgpu_state.surface.get_current_texture()?;
        let view = output
            .texture
            .create_view(&wgpu::TextureViewDescriptor::default());
        let mut encoder =
            self.wgpu_state
                .device
                .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                    label: Some("Render Encoder"),
                });

        for item in &mut self.items {
            item.render(&mut encoder, &view, &self.camera_bind_group)?;
        }

        // submit will accept anything that implements IntoIter
        self.wgpu_state
            .queue
            .submit(std::iter::once(encoder.finish()));
        output.present();

        Ok(())
    }
}

impl Canvas for Canvas2D {
    type WindowState = Canvas2DWindowState;

    fn new_window_state(&self, window: Arc<Window>) -> Self::WindowState {
        let mut state = Canvas2DWindowState::new(window);
        for item in &self.items {
            state
                .items
                .push(item.new_wgpu(&state.wgpu_state, &state.camera_bind_group_layout));
        }
        state
    }

    fn window_event(
        &mut self,
        window_state: &mut Self::WindowState,
        event_loop: &ActiveEventLoop,
        id: WindowId,
        event: WindowEvent,
    ) {
        self.camera.window_event(
            window_state.wgpu_state.size,
            self.mouse_pos,
            event_loop,
            id,
            &event,
        );

        #[allow(clippy::single_match)]
        match event {
            WindowEvent::CloseRequested => {
                event_loop.exit();
            }
            WindowEvent::Resized(physical_size) => {
                window_state.wgpu_state.resize(physical_size);
            }
            WindowEvent::RedrawRequested => {
                window_state.camera_uniform = self.camera.get_uniform(window_state.wgpu_state.size);
                match window_state.render() {
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
                        log::warn!("Surface timeout");
                    }
                }
                window_state.wgpu_state.window.request_redraw();
            }
            WindowEvent::CursorMoved { position, .. } => {
                self.mouse_pos = position;
            }
            WindowEvent::MouseInput { state, button, .. } => match (button, state) {
                (winit::event::MouseButton::Left, winit::event::ElementState::Pressed) => {
                    println!(
                        "{:?} -> {:?} -> {:?}",
                        self.mouse_pos,
                        self.camera
                            .pixel_to_coord(window_state.wgpu_state.size, self.mouse_pos),
                        self.camera.coord_to_pixel(
                            window_state.wgpu_state.size,
                            self.camera
                                .pixel_to_coord(window_state.wgpu_state.size, self.mouse_pos)
                        )
                    );
                }
                _ => {}
            },
            _ => (),
        }
    }
}

impl Canvas2D {
    pub fn new(camera: Box<dyn Camera>) -> Self {
        Self {
            mouse_pos: PhysicalPosition::new(0.0, 0.0),
            items: vec![],
            camera,
        }
    }

    pub fn add_item(&mut self, item: impl Canvas2DItem + 'static) {
        self.items.push(Box::new(item));
    }
}
