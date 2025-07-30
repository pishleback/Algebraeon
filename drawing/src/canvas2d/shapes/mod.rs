use crate::canvas2d::Canvas2D;
use crate::colour::Colour;

mod lines;
mod points;
mod rectangles;
mod triangles;

#[derive(Debug, Clone)]
struct State {
    colour: Colour,
    thickness: f32,
}

pub enum Shape {
    /// Set the colour
    Colour(Colour),
    /// Set the thickness as a % of the screen
    Thickness(f32),
    /// Push the draw parameters
    Push,
    /// Pop the draw parameters
    Pop,
    /// A point at (x, y)
    Point { x: f32, y: f32 },
    /// A line connecting (x1, y1) to (x2, y2)
    Line { x1: f32, y1: f32, x2: f32, y2: f32 },
    /// A box with edges x=a, x=b, y=c, y=d
    Box { a: f32, b: f32, c: f32, d: f32 },
    /// A line connecting (x1 + x1s, y1 + y1s) to (x2 + x2s, y2 + y2s)
    /// where (x1, y1) and (x2, y2) are the endpoints and (x1s, y1s) and (x2s, y2s) are offsets for the endpoints in screen coordinates relative to the thickness
    LineRaw {
        x1: f32,
        x1s: f32,
        y1: f32,
        y1s: f32,
        x2: f32,
        x2s: f32,
        y2: f32,
        y2s: f32,
    },
    /// A solid triangle
    Triangle {
        x1: f32,
        y1: f32,
        x2: f32,
        y2: f32,
        x3: f32,
        y3: f32,
    },
}

impl Canvas2D {
    pub fn plot_shapes(&mut self, shapes: impl IntoIterator<Item = Shape>) {
        let mut points = vec![];
        let mut lines = vec![];
        let mut rectangles = vec![];
        let mut triangles = vec![];

        let mut state = State {
            colour: Colour::black(),
            thickness: 0.01,
        };
        let mut state_stack = vec![];

        for shape in shapes {
            match shape {
                Shape::Colour(colour) => {
                    state.colour = colour;
                }
                Shape::Thickness(thickness) => {
                    state.thickness = 2.0 * thickness / 100.0;
                }
                Shape::Push => state_stack.push(state.clone()),
                Shape::Pop => {
                    state = state_stack.pop().unwrap();
                }
                Shape::Point { x, y } => points.push(points::Instance {
                    pos: [x, y],
                    radius: state.thickness,
                    colour: state.colour.rgb,
                }),
                Shape::Line { x1, y1, x2, y2 } => lines.push(lines::Instance {
                    pos1: [x1, y1],
                    pos1_screen_offset: [0.0, 0.0],
                    pos2: [x2, y2],
                    pos2_screen_offset: [0.0, 0.0],
                    radius: state.thickness,
                    colour: state.colour.rgb,
                }),
                Shape::LineRaw {
                    x1,
                    x1s,
                    y1,
                    y1s,
                    x2,
                    x2s,
                    y2,
                    y2s,
                } => lines.push(lines::Instance {
                    pos1: [x1, y1],
                    pos1_screen_offset: [state.thickness * x1s, state.thickness * y1s],
                    pos2: [x2, y2],
                    pos2_screen_offset: [state.thickness * x2s, state.thickness * y2s],
                    radius: state.thickness,
                    colour: state.colour.rgb,
                }),
                Shape::Box {
                    mut a,
                    mut b,
                    mut c,
                    mut d,
                } => {
                    if b < a {
                        (a, b) = (b, a);
                    }
                    if d < c {
                        (c, d) = (d, c);
                    }
                    rectangles.push(rectangles::Instance {
                        a,
                        b,
                        c,
                        d,
                        thickness: state.thickness,
                        colour: state.colour.rgb,
                    });
                }
                Shape::Triangle {
                    x1,
                    y1,
                    x2,
                    y2,
                    x3,
                    y3,
                } => {
                    triangles.push(triangles::Instance {
                        pos1: [x1, y1],
                        pos2: [x2, y2],
                        pos3: [x3, y3],
                        colour: state.colour.rgb,
                    });
                }
            }
        }

        assert!(state_stack.is_empty());

        if !triangles.is_empty() {
            self.add_item(triangles::TrianglesCanvas2DItem { triangles });
        }
        if !rectangles.is_empty() {
            self.add_item(rectangles::RectanglesCanvas2DItem { rectangles });
        }
        if !lines.is_empty() {
            self.add_item(lines::LinesCanvas2DItem { lines });
        }
        if !points.is_empty() {
            self.add_item(points::PointsCanvas2DItem { points });
        }
    }
}
