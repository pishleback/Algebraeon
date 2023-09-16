use glium::{glutin::event::Event, Display, Program, Surface};

#[derive(Debug)]
struct Camera {
    center: (f64, f64),
    //big scale = more zoomed out
    //scale = 1 means area of screen = 1 coords squared
    scale: f64,
}

impl Camera {
    fn pos_to_pixel(&self, size_px: (u32, u32), pos: (f64, f64)) -> (f64, f64) {
        assert!(size_px.0 > 0);
        assert!(size_px.1 > 0);
        let size_px = (size_px.0 as f64, size_px.1 as f64);
        let side = (size_px.0 * size_px.1).sqrt();
        (
            size_px.0 * ((pos.0 - self.center.0) * self.scale * side / size_px.0 + 0.5),
            size_px.1 - size_px.1 * ((pos.1 - self.center.1) * self.scale * side / size_px.1 + 0.5),
        )
    }

    fn pixel_to_pos(&self, size_px: (u32, u32), pixels: (f64, f64)) -> (f64, f64) {
        assert!(size_px.0 > 0);
        assert!(size_px.1 > 0);
        let size_px = (size_px.0 as f64, size_px.1 as f64);
        let side = (size_px.0 * size_px.1).sqrt();
        (
            self.center.0 + size_px.0 * (pixels.0 / size_px.0 - 0.5) / (self.scale * side),
            self.center.1
                + size_px.1 * ((size_px.1 - pixels.1) / size_px.1 - 0.5) / (self.scale * side),
        )
    }

    fn change_zoom(&mut self, size_px: (u32, u32), center_px: (f64, f64), mult: f64) {
        assert!(mult > 0.0);
        let old_center = self.pixel_to_pos(size_px, center_px);
        self.scale *= mult;
        let new_center = self.pixel_to_pos(size_px, center_px);
        self.center.0 += old_center.0 - new_center.0;
        self.center.1 += old_center.1 - new_center.1;
    }
}

pub struct Canvas {
    camera: Camera,
    program: Option<Program>,
}

impl Canvas {
    pub fn new() -> Self {
        Self {
            camera: Camera {
                center: (0.0, 0.0),
                scale: 1.0,
            },
            program: None,
        }
    }

    fn make_program(&mut self, display: &Display) {
        if self.program.is_none() {
            let vertex_shader_src = r#"
                #version 460

                in vec2 position;
                in vec4 colour;
                in vec2 offset;

                out vec4 v_colour;

                uniform vec2 camera_center;
                uniform float camera_scale;
                uniform vec2 display_size;

                void main() {
                    float side = sqrt(display_size.x * display_size.y);

                    gl_Position = vec4(
                        2 * (position.x + offset.x - camera_center.x) * camera_scale * side / display_size.x,
                        2 * (position.y + offset.y - camera_center.y) * camera_scale * side / display_size.y,
                        0.0,
                        1.0);
                    v_colour = colour;
                }
            "#;

            let fragment_shader_src = r#"
                #version 460

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
    }
}

impl super::Canvas for Canvas {
    fn tick(&mut self, _state: &super::State, _dt: f64) {}

    fn draw(&mut self, state: &super::State, display: &Display) {
        self.make_program(display);

        let mut target = display.draw();

        target.clear_color(0.0, 0.0, 0.0, 1.0);

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
            // Vertex {
            //     position: [0.5, 0.0],
            //     colour: [1.0, 0.0, 0.0, 1.0],
            // },
            // Vertex {
            //     position: [0.5, 0.5],
            //     colour: [0.0, 1.0, 0.0, 1.0],
            // },
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

        target
            .draw(
                (&vertex_buffer, instance_buffer.per_instance().unwrap()),
                &indices,
                &self.program.as_ref().unwrap(),
                // &glium::uniforms::EmptyUniforms,
                &uniform! {
                    camera_scale : self.camera.scale as f32,
                    camera_center : (self.camera.center.0 as f32, self.camera.center.1 as f32),
                    display_size : (state.display_size.0 as f32, state.display_size.1 as f32)
                },
                &Default::default(),
            )
            .unwrap();

        target.finish().unwrap();
    }

    fn event(&mut self, state: &super::State, ev: &Event<'_, ()>) {
        // println!("{:?}", ev);
        match ev {
            Event::WindowEvent { event, .. } => match event {
                glium::glutin::event::WindowEvent::MouseWheel { delta, .. } => match delta {
                    glium::glutin::event::MouseScrollDelta::LineDelta(_x, y) => {
                        self.camera.change_zoom(
                            state.display_size,
                            state.mouse_pos,
                            (1.1 as f64).powf(*y as f64),
                        );
                        // println!("{:?}", self.camera);

                        let pos = self
                            .camera
                            .pixel_to_pos(state.display_size, state.mouse_pos);
                        let px_again = self.camera.pos_to_pixel(state.display_size, pos);
                        println!("{:?} {:?} {:?}", state.mouse_pos, pos, px_again);
                    }
                    _ => {}
                },
                _ => {}
            },
            _ => {}
        }
    }
}
