use std::borrow::Borrow;

use crate::{
    drawing::canvas2d::Drawable,
    rings::ring_structure::structure::{
        FieldStructure, OrderedRingStructure, RealToFloatStructure,
    },
};

use super::{simplices::Simplex, AffineSpace};

impl<
        FS: OrderedRingStructure + FieldStructure + RealToFloatStructure,
        SP: Borrow<AffineSpace<FS>> + Clone,
    > Drawable for Simplex<FS, SP>
{
    fn draw(
        &self,
        canvas: &mut crate::drawing::canvas2d::Diagram2dCanvas,
        colour: (f32, f32, f32),
    ) {
        let space = self.ambient_space();
        let ordered_field = space.borrow().ordered_field();
        assert_eq!(space.borrow().linear_dimension(), Some(2));
        let points = self.points();
        match points.len() {
            0 => {}
            1 => {
                canvas.draw_point(
                    (
                        ordered_field.as_f32(points[0].coordinate(0)),
                        ordered_field.as_f32(points[0].coordinate(1)),
                    ),
                    colour,
                );
            }
            2 => {
                canvas.draw_line(
                    (
                        ordered_field.as_f32(points[0].coordinate(0)),
                        ordered_field.as_f32(points[0].coordinate(1)),
                    ),
                    (
                        ordered_field.as_f32(points[1].coordinate(0)),
                        ordered_field.as_f32(points[1].coordinate(1)),
                    ),
                    colour,
                );
            }
            3 => {
                canvas.draw_triangle(
                    (
                        ordered_field.as_f32(points[0].coordinate(0)),
                        ordered_field.as_f32(points[0].coordinate(1)),
                    ),
                    (
                        ordered_field.as_f32(points[1].coordinate(0)),
                        ordered_field.as_f32(points[1].coordinate(1)),
                    ),
                    (
                        ordered_field.as_f32(points[2].coordinate(0)),
                        ordered_field.as_f32(points[2].coordinate(1)),
                    ),
                    colour,
                );
            }
            _ => {
                unreachable!()
            }
        }
    }
}
