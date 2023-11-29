use glium::{glutin::event::Event, Display, Program, Surface, VertexBuffer};

use crate::{
    geometry::Shape,
    rings::{nzq::QQ, ring::Real},
};

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

#[derive(Debug, Copy, Clone)]
struct Vertex {
    position: [f32; 2],
    colour: [f32; 3],
}
implement_vertex!(Vertex, position, colour);

pub struct Shape2dCanvas {
    camera: Camera,
    point_program: Option<Program>,
    line_program: Option<Program>,
    triangle_program: Option<Program>,

    point_verts: Vec<Vertex>,
    line_verts: Vec<Vertex>,
    triangle_verts: Vec<Vertex>,

    points_vertex_buffer: Option<VertexBuffer<Vertex>>,
    lines_vertex_buffer: Option<VertexBuffer<Vertex>>,
    triangles_vertex_buffer: Option<VertexBuffer<Vertex>>,
}

impl Shape2dCanvas {
    pub fn new() -> Self {
        Self {
            camera: Camera {
                center: (0.0, 0.0),
                scale: 1.0,
            },
            point_program: None,
            line_program: None,
            triangle_program: None,
            point_verts: vec![],
            line_verts: vec![],
            triangle_verts: vec![],
            points_vertex_buffer: None,
            lines_vertex_buffer: None,
            triangles_vertex_buffer: None,
        }
    }

    fn make_program(&mut self, display: &Display) {
        if self.point_program.is_none() {
            self.point_program = Some(
                glium::Program::from_source(
                    display,
                    r#"
                        #version 460

                        in vec2 position;
                        in vec3 colour;

                        out vec4 g_colour;
                        out vec2 g_position;

                        void main() {
                            g_position = position;
                            g_colour = vec4(colour, 1.0);
                        }
                    "#,
                    r#"
                        #version 460

                        in vec4 v_colour;

                        out vec4 f_color;

                        void main() {
                            f_color = v_colour;
                        }
                    "#,
                    Some(r#"
                        #version 460

                        layout (points) in;
                        layout (triangle_strip, max_vertices = 36) out;

                        in vec4 g_colour[1];
                        in vec2 g_position[1];

                        out vec4 v_colour;

                        uniform vec2 camera_center;
                        uniform float camera_scale;
                        uniform vec2 display_size;

                        const float tau = 6.28318530718;

                        void main() {
                            float side = sqrt(display_size.x * display_size.y);

                            vec2 axis = 18.0 * vec2(1.0, 0.0) / (2.0 * camera_scale * side);

                            v_colour = g_colour[0];

                            for (int i = 0; i < 12; i ++) {
                                vec2 a = cos(tau * float(i) / 12.0) * axis + sin(tau * float(i) / 12.0) * vec2(axis.y, -axis.x);
                                vec2 b = cos(tau * float(i + 1) / 12.0) * axis + sin(tau * float(i + 1) / 12.0) * vec2(axis.y, -axis.x);

                                gl_Position = vec4(
                                    2 * (g_position[0] - camera_center) * camera_scale * side / display_size,
                                    0.0,
                                    1.0);
                                EmitVertex();
                                gl_Position = vec4(
                                    2 * (g_position[0] + a - camera_center) * camera_scale * side / display_size,
                                    0.0,
                                    1.0);
                                EmitVertex();
                                gl_Position = vec4(
                                    2 * (g_position[0] + b - camera_center) * camera_scale * side / display_size,
                                    0.0,
                                    1.0);
                                EmitVertex();
                                EndPrimitive();
                            }
                        }
                    "#),)
                    .unwrap(),
            );
        }

        if self.line_program.is_none() {
            self.line_program = Some(
                glium::Program::from_source(
                    display,
                    r#"
                        #version 460

                        in vec2 position;
                        in vec3 colour;

                        out vec4 g_colour;
                        out vec2 g_position;

                        void main() {
                            g_position = position;
                            g_colour = vec4(colour, 1.0);
                        }
                    "#,
                    r#"
                        #version 460

                        in vec4 v_colour;

                        out vec4 f_color;

                        void main() {
                            f_color = v_colour;
                        }
                    "#,
                    Some(r#"
                        #version 460

                        layout (lines) in;
                        layout (triangle_strip, max_vertices = 6) out;

                        in vec4 g_colour[2];
                        in vec2 g_position[2];

                        out vec4 v_colour;

                        uniform vec2 camera_center;
                        uniform float camera_scale;
                        uniform vec2 display_size;

                        void main() {
                            float side = sqrt(display_size.x * display_size.y);

                            vec2 perp = 5 * normalize(g_position[1] - g_position[0]) / (2 * camera_scale * side);
                            perp = vec2(perp.y, -perp.x);

                            v_colour = g_colour[0];
                            gl_Position = vec4(
                                2 * (g_position[0] + perp - camera_center) * camera_scale * side / display_size,
                                0.0,
                                 1.0);
                            EmitVertex();
                            gl_Position = vec4(
                                2 * (g_position[0] - perp - camera_center) * camera_scale * side / display_size,
                                0.0,
                                1.0);
                            EmitVertex();
                            v_colour = g_colour[1];
                            gl_Position = vec4(
                                2 * (g_position[1] + perp - camera_center) * camera_scale * side / display_size,
                                0.0,
                                1.0);
                            EmitVertex();
                            EndPrimitive();

                            v_colour = g_colour[0];
                            gl_Position = vec4(
                                2 * (g_position[0] - perp - camera_center) * camera_scale * side / display_size,
                                0.0,
                                1.0);
                            EmitVertex();
                            v_colour = g_colour[1];
                            gl_Position = vec4(
                                2 * (g_position[1] - perp - camera_center) * camera_scale * side / display_size,
                                0.0,
                                1.0);
                            EmitVertex();
                            gl_Position = vec4(
                                2 * (g_position[1] + perp - camera_center) * camera_scale * side / display_size,
                                0.0,
                                1.0);
                            EmitVertex();
                            EndPrimitive();
                        }
                    "#),
                )
                .unwrap(),
            );
        }

        if self.triangle_program.is_none() {
            self.triangle_program = Some(
                glium::Program::from_source(
                    display,
                    r#"
                        #version 460

                        in vec2 position;
                        in vec3 colour;

                        out vec4 v_colour;

                        uniform vec2 camera_center;
                        uniform float camera_scale;
                        uniform vec2 display_size;

                        void main() {
                            float side = sqrt(display_size.x * display_size.y);

                            gl_Position = vec4(
                                2 * (position - camera_center) * camera_scale * side / display_size,
                                0.0,
                                1.0);
                            v_colour = vec4(colour, 0.3);
                        }
                    "#,
                    r#"
                        #version 460

                        in vec4 v_colour;

                        out vec4 f_color;

                        void main() {
                            f_color = v_colour;
                        }
                    "#,
                    None,
                )
                .unwrap(),
            );
        }

        if self.points_vertex_buffer.is_none() {
            self.points_vertex_buffer =
                Some(glium::VertexBuffer::new(display, &self.point_verts).unwrap());
        }

        if self.lines_vertex_buffer.is_none() {
            self.lines_vertex_buffer =
                Some(glium::VertexBuffer::new(display, &self.line_verts).unwrap());
        }

        if self.triangles_vertex_buffer.is_none() {
            self.triangles_vertex_buffer =
                Some(glium::VertexBuffer::new(display, &self.triangle_verts).unwrap());
        }
    }

    pub fn draw_shape(&mut self, shape: &Shape, colour: (f32, f32, f32)) {
        assert_eq!(shape.dim(), 2);

        for simplex in shape.clone().simplices() {
            let points = simplex.points();

            match simplex.n() {
                1 => {
                    self.point_verts.push(Vertex {
                        position: [
                            QQ.as_f32(points[0].get_coord(0)),
                            QQ.as_f32(points[0].get_coord(1)),
                        ],
                        colour: [colour.0, colour.1, colour.2],
                    });
                    self.points_vertex_buffer = None;
                }
                2 => {
                    self.line_verts.push(Vertex {
                        position: [
                            QQ.as_f32(points[0].get_coord(0)),
                            QQ.as_f32(points[0].get_coord(1)),
                        ],
                        colour: [colour.0, colour.1, colour.2],
                    });
                    self.line_verts.push(Vertex {
                        position: [
                            QQ.as_f32(points[1].get_coord(0)),
                            QQ.as_f32(points[1].get_coord(1)),
                        ],
                        colour: [colour.0, colour.1, colour.2],
                    });
                    self.lines_vertex_buffer = None;
                }
                3 => {
                    self.triangle_verts.push(Vertex {
                        position: [
                            QQ.as_f32(points[0].get_coord(0)),
                            QQ.as_f32(points[0].get_coord(1)),
                        ],
                        colour: [colour.0, colour.1, colour.2],
                    });

                    self.triangle_verts.push(Vertex {
                        position: [
                            QQ.as_f32(points[1].get_coord(0)),
                            QQ.as_f32(points[1].get_coord(1)),
                        ],
                        colour: [colour.0, colour.1, colour.2],
                    });
                    self.triangle_verts.push(Vertex {
                        position: [
                            QQ.as_f32(points[2].get_coord(0)),
                            QQ.as_f32(points[2].get_coord(1)),
                        ],
                        colour: [colour.0, colour.1, colour.2],
                    });
                    self.triangles_vertex_buffer = None;
                }
                _ => panic!(),
            }
        }
    }

    // pub fn draw_partial_simplicial_complex(
    //     &mut self,
    //     psc: &PartialSimplicialComplex,
    //     colour: (f32, f32, f32),
    // ) {
    //     self.draw_shape(&psc.clone().as_shape(), colour)
    // }

    // pub fn draw_simplicial_complex(&mut self, sc: &SimplicialComplex, colour: (f32, f32, f32)) {
    //     self.draw_shape(&sc.clone().as_shape(), colour)
    // }
}

impl super::Canvas for Shape2dCanvas {
    fn tick(&mut self, _state: &super::State, _dt: f64) {}

    fn draw(&mut self, state: &super::State, display: &Display) {
        self.make_program(display);

        let mut target = display.draw();

        target.clear_color(0.0, 0.0, 0.0, 1.0);

        target
            .draw(
                self.triangles_vertex_buffer.as_ref().unwrap(),
                &glium::index::NoIndices(glium::index::PrimitiveType::TrianglesList),
                &self.triangle_program.as_ref().unwrap(),
                &uniform! {
                    camera_scale : self.camera.scale as f32,
                    camera_center : (self.camera.center.0 as f32, self.camera.center.1 as f32),
                    display_size : (state.display_size.0 as f32, state.display_size.1 as f32)
                },
                &glium::DrawParameters {
                    blend: glium::Blend::alpha_blending(),
                    ..Default::default()
                },
            )
            .unwrap();

        target
            .draw(
                self.lines_vertex_buffer.as_ref().unwrap(),
                &glium::index::NoIndices(glium::index::PrimitiveType::LinesList),
                &self.line_program.as_ref().unwrap(),
                &uniform! {
                    camera_scale : self.camera.scale as f32,
                    camera_center : (self.camera.center.0 as f32, self.camera.center.1 as f32),
                    display_size : (state.display_size.0 as f32, state.display_size.1 as f32)
                },
                &Default::default(),
            )
            .unwrap();

        target
            .draw(
                self.points_vertex_buffer.as_ref().unwrap(),
                &glium::index::NoIndices(glium::index::PrimitiveType::Points),
                &self.point_program.as_ref().unwrap(),
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
                    }
                    _ => {}
                },
                _ => {}
            },
            _ => {}
        }
    }
}

impl Shape {
    pub fn view2d(&self) -> ! {
        if self.dim() != 2 {
            panic!()
        }

        use super::Canvas;
        let mut canvas = Shape2dCanvas::new();
        canvas.draw_shape(self, (1.0, 1.0, 1.0));
        canvas.run()
    }
}
