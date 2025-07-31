use algebraeon_sets::structure::BorrowedStructure;

use crate::canvas2d::Canvas2D;
use crate::colour::Colour;

mod lines;
mod points;
mod rectangles;
mod triangles;

#[derive(Debug, Clone)]
struct State {
    colour: Colour,
    colour_alpha: f32,
    thickness: f32,
}

pub enum Shape {
    /// Set the colour
    SetColour(Colour),
    /// Set the colour alpha
    SetAlpha(f32),
    /// Set the thickness as a % of the screen
    SetThickness(f32),
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
            colour_alpha: 1.0,
            thickness: 0.01,
        };
        let mut state_stack = vec![];

        for shape in shapes {
            match shape {
                Shape::SetColour(colour) => {
                    state.colour = colour;
                }
                Shape::SetAlpha(alpha) => {
                    state.colour_alpha = alpha;
                }
                Shape::SetThickness(thickness) => {
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
                        colour: [
                            state.colour.rgb[0],
                            state.colour.rgb[1],
                            state.colour.rgb[2],
                            state.colour_alpha,
                        ],
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

pub fn simplicial_complex_shapes<
    FS: algebraeon_rings::structure::OrderedRingSignature
        + algebraeon_rings::structure::FieldSignature
        + algebraeon_rings::structure::RealToFloatSignature,
    FSB: BorrowedStructure<FS>,
    SP: std::borrow::Borrow<algebraeon_geometry::AffineSpace<FS, FSB>> + Clone,
    T: Eq + Clone,
>(
    line_colour: &Colour,
    fill_colour: &Colour,
    fill_alpha: f32,
    sc: &impl algebraeon_geometry::simplexes::LabelledSimplexCollection<FS, FSB, SP, T>,
) -> impl IntoIterator<Item = Shape>
where
    FS::Set: std::hash::Hash,
{
    let sp = sc.ambient_space();
    let field = sp.borrow().field();

    [Shape::Push, Shape::SetAlpha(1.0)]
        .into_iter()
        .chain(sc.simplexes().iter().flat_map(|s| {
            let pts = s.points();
            if pts.len() == 1 {
                let p = &pts[0];
                vec![
                    Shape::SetColour(line_colour.clone()),
                    Shape::Point {
                        x: field.as_f32(p.coordinate(0)),
                        y: field.as_f32(p.coordinate(1)),
                    },
                ]
            } else if pts.len() == 2 {
                let p1 = &pts[0];
                let p2 = &pts[1];
                vec![
                    Shape::SetColour(line_colour.clone()),
                    Shape::Line {
                        x1: field.as_f32(p1.coordinate(0)),
                        y1: field.as_f32(p1.coordinate(1)),
                        x2: field.as_f32(p2.coordinate(0)),
                        y2: field.as_f32(p2.coordinate(1)),
                    },
                ]
            } else if pts.len() == 3 {
                let p1 = &pts[0];
                let p2 = &pts[1];
                let p3 = &pts[2];
                vec![
                    Shape::SetColour(fill_colour.clone()),
                    Shape::Push,
                    Shape::SetAlpha(fill_alpha),
                    Shape::Triangle {
                        x1: field.as_f32(p1.coordinate(0)),
                        y1: field.as_f32(p1.coordinate(1)),
                        x2: field.as_f32(p2.coordinate(0)),
                        y2: field.as_f32(p2.coordinate(1)),
                        x3: field.as_f32(p3.coordinate(0)),
                        y3: field.as_f32(p3.coordinate(1)),
                    },
                    Shape::Pop,
                ]
            } else {
                vec![]
            }
        }))
        .chain([Shape::Pop])
        .collect::<Vec<_>>()
}
