use std::borrow::Borrow;

use crate::canvas_old::canvas2d::Drawable;
use algebraeon_rings::structure::{FieldSignature, OrderedRingSignature, RealToFloatSignature};

use algebraeon_geometry::{
    AffineSpace,
    simplexes::{
        LabelledSimplexCollection, LabelledSimplicialComplex, LabelledSimplicialDisjointUnion,
        PartialSimplicialComplex, Simplex,
    },
};

impl<
    FS: OrderedRingSignature + FieldSignature + RealToFloatSignature,
    SP: Borrow<AffineSpace<FS>> + Clone,
> Drawable for Simplex<FS, SP>
{
    fn draw(&self, canvas: &mut crate::canvas_old::canvas2d::Diagram2dCanvas, colour: (f32, f32, f32)) {
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

impl<
    FS: OrderedRingSignature + FieldSignature + RealToFloatSignature,
    SP: Borrow<AffineSpace<FS>> + Clone,
    T: Eq + Clone,
> Drawable for LabelledSimplicialComplex<FS, SP, T>
where
    FS::Set: std::hash::Hash,
{
    fn draw(&self, canvas: &mut crate::canvas_old::canvas2d::Diagram2dCanvas, colour: (f32, f32, f32)) {
        for simplex in self.simplexes() {
            simplex.draw(canvas, colour);
        }
    }
}

impl<
    FS: OrderedRingSignature + FieldSignature + RealToFloatSignature,
    SP: Borrow<AffineSpace<FS>> + Clone,
> Drawable for PartialSimplicialComplex<FS, SP>
where
    FS::Set: std::hash::Hash,
{
    fn draw(&self, canvas: &mut crate::canvas_old::canvas2d::Diagram2dCanvas, colour: (f32, f32, f32)) {
        for simplex in self.simplexes() {
            simplex.draw(canvas, colour);
        }
    }
}

impl<
    FS: OrderedRingSignature + FieldSignature + RealToFloatSignature,
    SP: Borrow<AffineSpace<FS>> + Clone,
    T: Eq + Clone,
> Drawable for LabelledSimplicialDisjointUnion<FS, SP, T>
where
    FS::Set: std::hash::Hash,
{
    fn draw(&self, canvas: &mut crate::canvas_old::canvas2d::Diagram2dCanvas, colour: (f32, f32, f32)) {
        for simplex in self.simplexes() {
            simplex.draw(canvas, colour);
        }
    }
}

// impl<
//         FS: OrderedRingStructure + FieldStructure + RealToFloatStructure,
//         SP: Borrow<AffineSpace<FS>> + Clone,
//         SC: Borrow<SimplicialComplex<FS, SP>> + Clone,
//     > Drawable for FullSubSimplicialComplex<FS, SP, SC>
// where
//     FS::Set: std::hash::Hash,
// {
//     fn draw(
//         &self,
//         canvas: &mut crate::drawing::canvas2d::Diagram2dCanvas,
//         colour: (f32, f32, f32),
//     ) {
//         for simplex in self.simplexes() {
//             simplex.draw(canvas, colour);
//         }
//     }
// }
