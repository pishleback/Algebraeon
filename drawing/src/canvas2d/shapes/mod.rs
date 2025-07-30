use crate::canvas2d::Canvas2D;

mod lines;
mod points;
mod rectangles;

pub enum Shape {
    /// A point at (x, y)
    Point { x: f32, y: f32 },
    /// A line connecting (x1, y1) to (x2, y2)
    Line { x1: f32, y1: f32, x2: f32, y2: f32 },
    /// A box with edges x=a, x=b, y=c, y=d
    Rectangle { a: f32, b: f32, c: f32, d: f32 },
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
}

impl Canvas2D {
    pub fn plot_shapes(&mut self, shapes: impl IntoIterator<Item = Shape>) {
        let mut points = vec![];
        let mut lines = vec![];
        let mut rectangles = vec![];

        let thickness = 0.01;

        for shape in shapes {
            match shape {
                Shape::Point { x, y } => points.push(points::Instance {
                    pos: [x, y],
                    radius: thickness,
                }),
                Shape::Line { x1, y1, x2, y2 } => lines.push(lines::Instance {
                    pos1: [x1, y1],
                    pos1_screen_offset: [0.0, 0.0],
                    pos2: [x2, y2],
                    pos2_screen_offset: [0.0, 0.0],
                    radius: thickness,
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
                    pos1_screen_offset: [thickness * x1s, thickness * y1s],
                    pos2: [x2, y2],
                    pos2_screen_offset: [thickness * x2s, thickness * y2s],
                    radius: thickness,
                }),
                Shape::Rectangle {
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
                        thickness,
                    });
                }
            }
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
