use algebraeon_sets::structure::SetSignature;

/// A directed simple graph permitting loops
pub trait GraphSignature<Vertices: SetSignature> {
    fn vertices(&self) -> &Vertices;

    fn has_directed_edge(
        &self,
        source: &Vertices::Set,
        target: &Vertices::Set,
    ) -> Result<(), String>;
}

/// A graph without loops
pub trait LooplessGraphSignature<Vertices: SetSignature>: GraphSignature<Vertices> {}

/// A graph such that has_directed_edge(u, v) == has_directed_edge(v, u)
pub trait UndirectedGraphSignature<Vertices: SetSignature>: GraphSignature<Vertices> {}

pub trait GraphWithEdgesSignature<Vertices: SetSignature, Edges: SetSignature>:
    GraphSignature<Vertices>
{
    /// return an edge if has_directed_edge is true
    fn directed_edge(
        &self,
        source: &Vertices::Set,
        target: &Vertices::Set,
    ) -> Result<&Edges::Set, String>;

    fn source(&self, edge: &Edges::Set) -> &Vertices::Set;
    fn target(&self, edge: &Edges::Set) -> &Vertices::Set;
    fn opposite(&self, edge: &Edges::Set) -> &Edges::Set;
}
