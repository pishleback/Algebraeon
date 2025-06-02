use super::{Canvas2DItem, Canvas2DItemWgpu};
use crate::canvas::WgpuState;
use algebraeon_rings::{
    polynomial::Polynomial,
    structure::{ComplexSubsetSignature, MetaComplexSubset},
};
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
            array_stride: std::mem::size_of::<Vertex>() as wgpu::BufferAddress,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &Self::ATTRIBS,
        }
    }
}

const VERTICES: &[Vertex] = &[
    Vertex {
        position: [-1.0, -1.0],
    },
    Vertex {
        position: [1.0, -1.0],
    },
    Vertex {
        position: [1.0, 1.0],
    },
    Vertex {
        position: [-1.0, 1.0],
    },
];

const INDICES: &[u16] = &[0, 1, 2, 0, 2, 3];

struct PolynomialWgpu {
    vertex_buffer: wgpu::Buffer,
    index_buffer: wgpu::Buffer,
    num_indices: u32,
    render_pipeline: wgpu::RenderPipeline,
}

pub struct PolynomialPlot {
    coeffs: Vec<(f64, f64)>,
}

impl PolynomialPlot {
    pub fn new<MetaRing: MetaComplexSubset>(p: Polynomial<MetaRing>) -> Self
    where
        MetaRing::Signature: ComplexSubsetSignature,
    {
        let coeffs = p
            .into_coeffs()
            .into_iter()
            .map(|c| MetaRing::as_f64_real_and_imaginary_parts(&c))
            .collect();
        Self { coeffs }
    }

    pub fn make_shader(&self) -> String {
        let n = self.coeffs.len();
        String::from(include_str!("shader.wgsl")).replace(
            "// Generate eval_cfn HERE",
            format!(
                r#"fn eval_cfn(z: vec2<f32>) -> vec2<f32> {{{}}}"#,
                match n {
                    0 => {
                        format!("return vec2<f32>(0.0, 0.0);")
                    }
                    1 => {
                        format!(
                            "return vec2<f32>({}, {});",
                            self.coeffs[0].0, self.coeffs[0].1
                        )
                    }
                    n => {
                        let n2 = n - 2;
                        let n1 = n - 1;
                        let wgsl_coeffs = self
                            .coeffs
                            .iter()
                            .map(|(a, b)| format!("vec2<f32>({a}, {b}), "))
                            .collect::<Vec<_>>()
                            .join("");
                        println!("{:?}", wgsl_coeffs);
                        format!(
                            r#"
var coeffs = array<vec2<f32>, {n}>(
    {wgsl_coeffs}
);
var result = coeffs[{n1}];
for (var i = {n2}; i >= 0i; i = i - 1i) {{
    result = c_add(c_mul(result, z), coeffs[i]);
}}
return result;
                "#
                        )
                    }
                }
            )
            .as_str(),
        )
    }
}

impl Canvas2DItem for PolynomialPlot {
    fn new(
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
                source: wgpu::ShaderSource::Wgsl(self.make_shader().into()),
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
                        buffers: &[Vertex::desc()],
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
                        topology: wgpu::PrimitiveTopology::TriangleList,
                        strip_index_format: None,
                        front_face: wgpu::FrontFace::Ccw,
                        cull_mode: None,
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

        Box::new(PolynomialWgpu {
            vertex_buffer,
            index_buffer,
            num_indices,
            render_pipeline,
        })
    }
}

impl Canvas2DItemWgpu for PolynomialWgpu {
    fn render(
        &mut self,
        encoder: &mut CommandEncoder,
        view: &TextureView,
        camera_bind_group: &BindGroup,
    ) -> Result<(), wgpu::SurfaceError> {
        let mut render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
            label: Some("Render Pass"),
            color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                view: view,
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
        render_pass.set_bind_group(0, camera_bind_group, &[]);
        render_pass.set_vertex_buffer(0, self.vertex_buffer.slice(..));
        render_pass.set_index_buffer(self.index_buffer.slice(..), wgpu::IndexFormat::Uint16);

        render_pass.draw_indexed(0..self.num_indices, 0, 0..1);
        Ok(())
    }
}
