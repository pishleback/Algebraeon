use crate::canvas::*;
use std::sync::Arc;
use wgpu::util::DeviceExt;
use winit::{
    dpi::{PhysicalPosition, PhysicalSize},
    event::WindowEvent,
    event_loop::ActiveEventLoop,
    window::{Window, WindowId},
};

#[repr(C)]
#[derive(Debug, Copy, Clone, bytemuck::Pod, bytemuck::Zeroable)]
pub struct CameraUniform {
    pub matrix: [[f32; 2]; 2],
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
            2.0 * pixels.x as f64 / display_size.width as f64 - 1.0,
            -2.0 * pixels.y as f64 / display_size.height as f64 + 1.0,
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

impl MouseWheelZoomCamera {
    pub fn new() -> Self {
        Self {
            mid_x: -21.918549,
            mid_y: -44.939706,
            sqrt_area: 200.0,
        }
    }
}

impl Camera for MouseWheelZoomCamera {
    fn get_uniform(&self, display_size: PhysicalSize<u32>) -> CameraUniform {
        let display_size = (display_size.width as f32, display_size.height as f32);
        let avg_side = (display_size.0 * display_size.1).sqrt();
        let x_mult = 2.0 * avg_side / (display_size.0 * self.sqrt_area as f32);
        let y_mult = 2.0 * avg_side / (display_size.1 * self.sqrt_area as f32);
        CameraUniform {
            matrix: [[x_mult, 0.0], [0.0, y_mult]],
            shift: [-self.mid_x as f32 * x_mult, -self.mid_y as f32 * y_mult],
        }
    }

    fn window_event(
        &mut self,
        display_size: PhysicalSize<u32>,
        mouse_pos: PhysicalPosition<f64>,
        event_loop: &ActiveEventLoop,
        id: WindowId,
        event: &WindowEvent,
    ) {
        match event {
            WindowEvent::MouseWheel {
                device_id,
                delta,
                phase,
            } => {
                let dy = match delta {
                    winit::event::MouseScrollDelta::PixelDelta(pos) => pos.y,
                    winit::event::MouseScrollDelta::LineDelta(x, y) => *y as f64,
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
    camera: Box<dyn Camera>,
}

#[repr(C)]
#[derive(Copy, Clone, Debug, bytemuck::Pod, bytemuck::Zeroable)]
struct Vertex {
    position: [f32; 2],
    color: [f32; 3],
}

impl Vertex {
    const ATTRIBS: [wgpu::VertexAttribute; 2] =
        wgpu::vertex_attr_array![0 => Float32x2, 1 => Float32x3];

    fn desc() -> wgpu::VertexBufferLayout<'static> {
        wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<Vertex>() as wgpu::BufferAddress,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &Self::ATTRIBS,
        }
    }
}

const VERTICES: &[Vertex] = &[
    Vertex {
        position: [-8.68241, 49.240386],
        color: [0.5, 0.0, 0.5],
    }, // A
    Vertex {
        position: [-49.513406, 6.958647],
        color: [0.5, 0.5, 0.5],
    }, // B
    Vertex {
        position: [-21.918549, -44.939706],
        color: [0.5, 0.0, 0.5],
    }, // C
    Vertex {
        position: [35.966998, -34.73291],
        color: [0.5, 0.0, 0.5],
    }, // D
    Vertex {
        position: [44.147372, 23.47359],
        color: [0.5, 0.0, 0.5],
    }, // E
];

const INDICES: &[u16] = &[0, 1, 4, 1, 2, 4, 2, 3, 4];

pub struct Canvas2DWindowState {
    wgpu_state: WgpuState,

    vertex_buffer: wgpu::Buffer,
    num_vertices: u32,
    index_buffer: wgpu::Buffer,
    num_indices: u32,

    camera_uniform: CameraUniform,
    camera_buffer: wgpu::Buffer,
    camera_bind_group: wgpu::BindGroup,

    render_pipeline: wgpu::RenderPipeline,
}

impl Canvas2DWindowState {
    fn new(window: Arc<Window>) -> Self {
        let wgpu_state = WgpuState::new(window);

        let camera_uniform = CameraUniform {
            matrix: [[1.0, 0.0], [0.0, 1.0]],
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

        let vertex_buffer =
            wgpu_state
                .device
                .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some("Vertex Buffer"),
                    contents: bytemuck::cast_slice(VERTICES),
                    usage: wgpu::BufferUsages::VERTEX,
                });

        let num_vertices = VERTICES.len() as u32;

        let index_buffer =
            wgpu_state
                .device
                .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some("Index Buffer"),
                    contents: bytemuck::cast_slice(INDICES),
                    usage: wgpu::BufferUsages::INDEX,
                });

        let num_indices = INDICES.len() as u32;

        let shader = wgpu_state
            .device
            .create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some("Shader"),
                source: wgpu::ShaderSource::Wgsl(include_str!("shader.wgsl").into()),
            });

        let render_pipeline_layout =
            wgpu_state
                .device
                .create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                    label: Some("Render Pipeline Layout"),
                    bind_group_layouts: &[&camera_bind_group_layout],
                    push_constant_ranges: &[],
                });

        let render_pipeline =
            wgpu_state
                .device
                .create_render_pipeline(&wgpu::RenderPipelineDescriptor {
                    label: Some("Render Pipeline"),
                    layout: Some(&render_pipeline_layout),
                    vertex: wgpu::VertexState {
                        module: &shader,
                        entry_point: Some("vs_main"), // 1.
                        buffers: &[Vertex::desc()],   // 2.
                        compilation_options: wgpu::PipelineCompilationOptions::default(),
                    },
                    fragment: Some(wgpu::FragmentState {
                        // 3.
                        module: &shader,
                        entry_point: Some("fs_main"),
                        targets: &[Some(wgpu::ColorTargetState {
                            // 4.
                            format: wgpu_state.config.format,
                            blend: Some(wgpu::BlendState::REPLACE),
                            write_mask: wgpu::ColorWrites::ALL,
                        })],
                        compilation_options: wgpu::PipelineCompilationOptions::default(),
                    }),
                    primitive: wgpu::PrimitiveState {
                        topology: wgpu::PrimitiveTopology::TriangleList, // 1.
                        strip_index_format: None,
                        front_face: wgpu::FrontFace::Ccw, // 2.
                        cull_mode: Some(wgpu::Face::Back),
                        // Setting this to anything other than Fill requires Features::NON_FILL_POLYGON_MODE
                        polygon_mode: wgpu::PolygonMode::Fill,
                        // Requires Features::DEPTH_CLIP_CONTROL
                        unclipped_depth: false,
                        // Requires Features::CONSERVATIVE_RASTERIZATION
                        conservative: false,
                    },
                    depth_stencil: None, // 1.
                    multisample: wgpu::MultisampleState {
                        count: 1,                         // 2.
                        mask: !0,                         // 3.
                        alpha_to_coverage_enabled: false, // 4.
                    },
                    multiview: None, // 5.
                    cache: None,     // 6.
                });

        Self {
            wgpu_state,
            vertex_buffer,
            num_vertices,
            index_buffer,
            num_indices,
            camera_uniform,
            camera_bind_group,
            camera_buffer,
            render_pipeline,
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

        {
            let mut render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("Render Pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color {
                            r: 0.1,
                            g: 0.2,
                            b: 0.3,
                            a: 1.0,
                        }),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: None,
                occlusion_query_set: None,
                timestamp_writes: None,
            });

            render_pass.set_pipeline(&self.render_pipeline);
            render_pass.set_bind_group(0, &self.camera_bind_group, &[]);
            render_pass.set_vertex_buffer(0, self.vertex_buffer.slice(..));
            render_pass.set_index_buffer(self.index_buffer.slice(..), wgpu::IndexFormat::Uint16);

            render_pass.draw_indexed(0..self.num_indices, 0, 0..1);
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

    fn new_window_state(window: Arc<Window>) -> Self::WindowState {
        Canvas2DWindowState::new(window)
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
                        log::warn!("Surface timeout")
                    }
                }
                window_state.wgpu_state.window.request_redraw();
            }
            WindowEvent::CursorMoved {
                device_id,
                position,
            } => {
                self.mouse_pos = position;
            }
            WindowEvent::MouseInput {
                device_id,
                state,
                button,
            } => match (button, state) {
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
            camera,
        }
    }
}
