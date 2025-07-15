use algebraeon_sets::structure::{EqSignature, SetSignature, UnorderedPair, UnorderedPairs};

use crate::structure::{
    GraphSignature, GraphWithEdgesSignature, LooplessGraphSignature, UndirectedGraphSignature,
};

#[allow(dead_code)]
pub struct CompleteUndirectedGraph<Vertices: SetSignature> {
    vertices: Vertices,
    pairs_of_vertices: UnorderedPairs<Vertices>,
}

impl<Vertices: SetSignature> CompleteUndirectedGraph<Vertices> {
    pub fn new(vertices: Vertices) -> Self {
        let pairs_of_vertices = UnorderedPairs::new(vertices.clone());
        Self {
            vertices,
            pairs_of_vertices,
        }
    }
}

impl<Vertices: SetSignature + EqSignature> GraphSignature for CompleteUndirectedGraph<Vertices> {
    type Vertices = Vertices;

    fn has_directed_edge(
        &self,
        source: &Vertices::Set,
        target: &Vertices::Set,
    ) -> Result<(), String> {
        if let Err(e) = self.vertices.is_element(source) {
            return Err(format!("Source is not an element of Vertices: {e}"));
        }
        if let Err(e) = self.vertices.is_element(target) {
            return Err(format!("Target is not an element of Vertices: {e}"));
        }
        if self.vertices.equal(source, target) {
            return Err("Complete graphs do not contain loops.".to_string());
        }
        // For undirected graphs, edge (u,v) exists iff edge (v,u) exists
        // Since this is a complete graph, any two distinct vertices are connected
        Ok(())
    }
}

impl<Vertices: SetSignature + EqSignature> LooplessGraphSignature
    for CompleteUndirectedGraph<Vertices>
{
}

impl<Vertices: SetSignature + EqSignature> UndirectedGraphSignature
    for CompleteUndirectedGraph<Vertices>
{
}

impl<Vertices: SetSignature + EqSignature> GraphWithEdgesSignature
    for CompleteUndirectedGraph<Vertices>
{
    type Edges = UnorderedPairs<Vertices>;

    fn endpoints(
        &self,
        edge: &<Self::Edges as SetSignature>::Set,
    ) -> UnorderedPair<<Self::Vertices as SetSignature>::Set> {
        edge.clone()
    }
}

#[cfg(test)]
mod tests {
    use algebraeon_sets::structure::{EnumeratedFiniteSetStructure, EqSignature};

    use crate::{
        examples::CompleteUndirectedGraph,
        structure::{GraphSignature, GraphWithEdgesSignature},
    };

    #[test]
    fn test_k5_undirected() {
        let fin5 = EnumeratedFiniteSetStructure::new(5);
        let k5 = CompleteUndirectedGraph::new(fin5.clone());

        // Test that edges exist in both directions
        assert!(k5.has_directed_edge(&1, &2).is_ok());
        assert!(k5.has_directed_edge(&2, &1).is_ok());

        // Test that loops are not allowed
        assert!(k5.has_directed_edge(&1, &1).is_err());

        // Test invalid vertices
        assert!(k5.has_directed_edge(&5, &1).is_err());
        assert!(k5.has_directed_edge(&1, &6).is_err());
    }

    #[test]
    fn test_k5_undirected_edges() {
        let fin5 = EnumeratedFiniteSetStructure::new(5);
        let k5 = CompleteUndirectedGraph::new(fin5.clone());

        // Create an edge between vertices 1 and 2
        let edge_12 = k5.pairs_of_vertices.new_pair(&1, &2).unwrap();
        let endpoints = k5.endpoints(&edge_12);

        // Verify the endpoints are the same as the original edge
        // (since for undirected graphs, the edge is just the unordered pair itself)
        assert!(k5.pairs_of_vertices.equal(&edge_12, &endpoints));
    }

    #[test]
    fn test_undirected_graph_symmetry() {
        let fin3 = EnumeratedFiniteSetStructure::new(3);
        let k3 = CompleteUndirectedGraph::new(fin3);

        // For an undirected graph, has_directed_edge should be symmetric
        for i in 0..3 {
            for j in 0..3 {
                if i != j {
                    assert_eq!(
                        k3.has_directed_edge(&i, &j).is_ok(),
                        k3.has_directed_edge(&j, &i).is_ok(),
                        "Edge ({i}, {j}) symmetry failed"
                    );
                }
            }
        }
    }
}
