use std::time::Duration;
use algebraeon_drawing::canvas2d::*;
use glium::{
    DrawParameters, Program, Surface, glutin::surface::WindowSurface, implement_vertex, uniform,
    winit::event::Event,
};

struct TestElement {
    program: Option<Program>,
}

impl TestElement {
    pub fn new() -> Self {
        Self { program: None }
    }
}

impl EventHandler for TestElement {
    fn tick(&mut self, dt: &Duration) {}
    fn event(&mut self, ev: &Event<()>) {}
}

impl DrawElement for TestElement {
    fn draw(
        &mut self,
        display: &glium::Display<WindowSurface>,
        target: &mut glium::Frame,
        display_size: (u32, u32),
        camera: &Box<dyn Camera>,
    ) {
        if self.program.is_none() {
            let vertex_shader_src = r#"
                #version 330

                in vec2 position;
                in vec4 colour;
                in vec2 offset;

                out vec4 v_colour;

                uniform mat2 view_matrix;
                uniform vec2 view_offset;

                void main() {
                    gl_Position = vec4(
                        view_matrix * (offset + position) + view_offset,
                        0.0,
                        1.0
                    );
                    v_colour = colour;
                }
            "#;

            let fragment_shader_src = r#"
                #version 330

                in vec4 v_colour;

                out vec4 f_color;

                void main() {
                    f_color = v_colour;
                }
            "#;

            self.program = Some(
                glium::Program::from_source(display, vertex_shader_src, fragment_shader_src, None)
                    .unwrap(),
            );

        }
        let program = self.program.as_ref().unwrap();

        #[derive(Copy, Clone)]
        struct Vertex {
            position: [f32; 2],
            colour: [f32; 4],
        }
        implement_vertex!(Vertex, position, colour);

        let shape = vec![
            Vertex {
                position: [0.0, 0.0],
                colour: [1.0, 0.0, 0.0, 1.0],
            },
            Vertex {
                position: [0.0, 1.0],
                colour: [0.0, 1.0, 0.0, 1.0],
            },
            Vertex {
                position: [1.0, 0.0],
                colour: [0.0, 0.0, 1.0, 1.0],
            },
        ];

        #[derive(Copy, Clone)]
        struct Offset {
            offset: [f32; 2],
        }
        implement_vertex!(Offset, offset);
        let mut instances = vec![];
        for x in 0..10 {
            for y in 0..10 {
                instances.push(Offset {
                    offset: [x as f32, y as f32],
                });
            }
        }

        let vertex_buffer = glium::VertexBuffer::new(display, &shape).unwrap();
        let instance_buffer = glium::VertexBuffer::new(display, &instances).unwrap();
        let indices = glium::index::NoIndices(glium::index::PrimitiveType::TrianglesList);

        let (view_matrix, view_offset) = camera.view_matrix(display_size);

        target
            .draw(
                (&vertex_buffer, instance_buffer.per_instance().unwrap()),
                &indices,
                program,
                &uniform! {view_matrix : view_matrix,
                view_offset : view_offset },
                &DrawParameters::default(),
            )
            .unwrap();
    }
}

fn main() {
    let mut canvas = Canvas2D::new(Box::new(MouseWheelZoomCamera::new()));

    canvas.add_element(Box::new(TestElement::new()));

    canvas.run(1000, 600);
}
