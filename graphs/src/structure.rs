use algebraeon_sets::structure::{SetSignature, UnorderedPair};

/// A directed graph permitting loops.
pub trait GraphSignature {
    type Vertices: SetSignature;

    fn has_directed_edge(
        &self,
        source: &<Self::Vertices as SetSignature>::Set,
        target: &<Self::Vertices as SetSignature>::Set,
    ) -> Result<(), String>;
}

/// A graph with no loops.
pub trait LooplessGraphSignature: GraphSignature {}

/// A graph such that `has_directed_edge(u, v) == has_directed_edge(v, u)`
pub trait UndirectedGraphSignature: GraphSignature {}

pub trait GraphWithEdgesSignature: GraphSignature {
    type Edges: SetSignature;

    /// Return the endpoints of an edge.
    fn endpoints(
        &self,
        edge: &<Self::Edges as SetSignature>::Set,
    ) -> UnorderedPair<<Self::Vertices as SetSignature>::Set>;
}
