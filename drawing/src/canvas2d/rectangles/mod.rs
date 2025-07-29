use super::{Canvas2DItem, Canvas2DItemWgpu};
use crate::{canvas::WgpuState, canvas2d::Canvas2D};
use wgpu::{BindGroup, BindGroupLayout, CommandEncoder, TextureView, util::DeviceExt};

#[repr(C)]
#[derive(Copy, Clone, Debug, bytemuck::Pod, bytemuck::Zeroable)]
struct Vertex {
    position: [f32; 2],
}

impl Vertex {
    const ATTRIBS: [wgpu::VertexAttribute; 1] = wgpu::vertex_attr_array![0 => Float32x2];

    fn desc() -> wgpu::VertexBufferLayout<'static> {
        wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<Self>() as wgpu::BufferAddress,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &Self::ATTRIBS,
        }
    }
}

// [-1, 0] is a buffer zone
// [0, 1] is the rect
// [1, 2] is a buffer zone
const VERTICES: &[Vertex] = &[
    Vertex {
        position: [-1.0, -1.0],
    },
    Vertex {
        position: [2.0, -1.0],
    },
    Vertex {
        position: [-1.0, 2.0],
    },
    Vertex {
        position: [2.0, 2.0],
    },
];

const INDICES: &[u16] = &[0, 1, 2, 1, 3, 2];

#[repr(C)]
#[derive(Copy, Clone, bytemuck::Pod, bytemuck::Zeroable)]
struct Instance {
    a: f32,
    b: f32,
    c: f32,
    d: f32,
    thickness: f32,
}

impl Instance {
    const ATTRIBS: [wgpu::VertexAttribute; 5] = wgpu::vertex_attr_array![1 => Float32, 2 => Float32, 3 => Float32, 4 => Float32, 5 => Float32];

    fn desc() -> wgpu::VertexBufferLayout<'static> {
        wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<Self>() as wgpu::BufferAddress,
            step_mode: wgpu::VertexStepMode::Instance,
            attributes: &Self::ATTRIBS,
        }
    }
}

struct PointsCanvas2DWgpu {
    vertex_buffer: wgpu::Buffer,
    _num_vertices: u32,
    index_buffer: wgpu::Buffer,
    num_indices: u32,
    instance_buffer: wgpu::Buffer,
    num_instances: u32,
    render_pipeline: wgpu::RenderPipeline,
}

impl Canvas2DItemWgpu for PointsCanvas2DWgpu {
    fn render(
        &mut self,
        encoder: &mut CommandEncoder,
        view: &TextureView,
        camera_bind_group: &BindGroup,
    ) -> Result<(), wgpu::SurfaceError> {
        if self.num_instances == 0 {
            return Ok(());
        }

        let mut render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
            label: Some("Render Pass"),
            color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                view,
                resolve_target: None,
                ops: wgpu::Operations {
                    load: wgpu::LoadOp::Load,
                    store: wgpu::StoreOp::Store,
                },
            })],
            depth_stencil_attachment: None,
            occlusion_query_set: None,
            timestamp_writes: None,
        });

        render_pass.set_pipeline(&self.render_pipeline);
        render_pass.set_bind_group(0, camera_bind_group, &[]);
        render_pass.set_vertex_buffer(0, self.vertex_buffer.slice(..));
        render_pass.set_vertex_buffer(1, self.instance_buffer.slice(..));
        render_pass.set_index_buffer(self.index_buffer.slice(..), wgpu::IndexFormat::Uint16);
        render_pass.draw_indexed(0..self.num_indices, 0, 0..self.num_instances);

        Ok(())
    }
}

struct RectanglesCanvas2DItem {
    rectangles: Vec<Instance>,
}

impl Canvas2DItem for RectanglesCanvas2DItem {
    fn new_wgpu(
        &self,
        wgpu_state: &WgpuState,
        camera_bind_group_layout: &BindGroupLayout,
    ) -> Box<dyn Canvas2DItemWgpu> {
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
                    bind_group_layouts: &[camera_bind_group_layout],
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
                        entry_point: Some("vs_main"),
                        buffers: &[Vertex::desc(), Instance::desc()],
                        compilation_options: wgpu::PipelineCompilationOptions::default(),
                    },
                    fragment: Some(wgpu::FragmentState {
                        module: &shader,
                        entry_point: Some("fs_main"),
                        targets: &[Some(wgpu::ColorTargetState {
                            format: wgpu_state.config.format,
                            blend: Some(wgpu::BlendState::REPLACE),
                            write_mask: wgpu::ColorWrites::ALL,
                        })],
                        compilation_options: wgpu::PipelineCompilationOptions::default(),
                    }),
                    primitive: wgpu::PrimitiveState {
                        topology: wgpu::PrimitiveTopology::TriangleList,
                        strip_index_format: None,
                        front_face: wgpu::FrontFace::Ccw,
                        cull_mode: Some(wgpu::Face::Back),
                        // Setting this to anything other than Fill requires Features::NON_FILL_POLYGON_MODE
                        polygon_mode: wgpu::PolygonMode::Fill,
                        // Requires Features::DEPTH_CLIP_CONTROL
                        unclipped_depth: false,
                        // Requires Features::CONSERVATIVE_RASTERIZATION
                        conservative: false,
                    },
                    depth_stencil: None,
                    multisample: wgpu::MultisampleState {
                        count: 1,
                        mask: !0,
                        alpha_to_coverage_enabled: false,
                    },
                    multiview: None,
                    cache: None,
                });

        let instance_buffer =
            wgpu_state
                .device
                .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some("Instance Buffer"),
                    contents: bytemuck::cast_slice(&self.rectangles),
                    usage: wgpu::BufferUsages::VERTEX,
                });

        let num_instances = self.rectangles.len() as u32;

        Box::new(PointsCanvas2DWgpu {
            vertex_buffer,
            _num_vertices: num_vertices,
            index_buffer,
            num_indices,
            instance_buffer,
            num_instances,
            render_pipeline,
        })
    }
}

pub struct Rectangle {
    pub a: f32,
    pub b: f32,
    pub c: f32,
    pub d: f32,
}

impl Canvas2D {
    pub fn plot_rectangles(&mut self, rectangles: impl Iterator<Item = Rectangle>) {
        self.add_item(RectanglesCanvas2DItem {
            rectangles: rectangles
                .map(|r| Instance {
                    a: r.a,
                    b: r.b,
                    c: r.c,
                    d: r.d,
                    thickness: 0.01,
                })
                .collect(),
        });
    }
}
