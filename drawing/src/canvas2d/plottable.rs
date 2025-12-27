use crate::{
    canvas2d::{
        Canvas2D,
        shapes::{Shape, simplicial_complex_shapes},
    },
    colour::Colour,
};
use algebraeon_geometry::{
    boolean_operations::Difference, partial_simplicial_complex::PartialSimplicialComplex,
    simplex_collection::LabelledSimplexCollection, simplicial_complex::SimplicialComplex,
    simplicial_disjoint_union::SimplicialDisjointUnion,
};

pub trait Plottable {
    fn plot(&self, canvas: &mut Canvas2D);
}

impl Canvas2D {
    pub fn plot(&mut self, p: impl Plottable) {
        p.plot(self);
    }
}

fn plot_simplex_collection<
    'f,
    FS: algebraeon_rings::structure::OrderedRingSignature
        + algebraeon_rings::structure::FieldSignature
        + algebraeon_rings::structure::RealSubsetSignature
        + 'f,
>(
    canvas: &mut Canvas2D,
    shape: &impl LabelledSimplexCollection<'f, FS, ()>,
) where
    FS::Set: std::hash::Hash,
{
    let x = shape.to_partial_simplicial_complex();
    let y = x.closure().difference(&x);
    let z = y.closure().difference(&y);

    canvas.plot_shapes(
        [Shape::SetThickness(0.3)]
            .into_iter()
            .chain(simplicial_complex_shapes(
                &Colour::green(),
                &Colour::green(),
                0.5,
                &x,
            )),
    );
    canvas.plot_shapes(
        [Shape::SetThickness(0.3)]
            .into_iter()
            .chain(simplicial_complex_shapes(
                &Colour::red(),
                &Colour::red(),
                0.5,
                &y,
            )),
    );
    canvas.plot_shapes(
        [Shape::SetThickness(0.3)]
            .into_iter()
            .chain(simplicial_complex_shapes(
                &Colour::green(),
                &Colour::green(),
                0.5,
                &z,
            )),
    );
}

impl<
    'f,
    FS: algebraeon_rings::structure::OrderedRingSignature
        + algebraeon_rings::structure::FieldSignature
        + algebraeon_rings::structure::RealSubsetSignature
        + 'f,
> Plottable for SimplicialDisjointUnion<'f, FS>
where
    FS::Set: std::hash::Hash,
{
    fn plot(&self, canvas: &mut Canvas2D) {
        plot_simplex_collection(canvas, self);
    }
}

impl<
    'f,
    FS: algebraeon_rings::structure::OrderedRingSignature
        + algebraeon_rings::structure::FieldSignature
        + algebraeon_rings::structure::RealSubsetSignature
        + 'f,
> Plottable for PartialSimplicialComplex<'f, FS>
where
    FS::Set: std::hash::Hash,
{
    fn plot(&self, canvas: &mut Canvas2D) {
        plot_simplex_collection(canvas, self);
    }
}

impl<
    'f,
    FS: algebraeon_rings::structure::OrderedRingSignature
        + algebraeon_rings::structure::FieldSignature
        + algebraeon_rings::structure::RealSubsetSignature
        + 'f,
> Plottable for SimplicialComplex<'f, FS>
where
    FS::Set: std::hash::Hash,
{
    fn plot(&self, canvas: &mut Canvas2D) {
        plot_simplex_collection(canvas, self);
    }
}
