use algebraeon_sets::structure::{SetSignature, UnorderedPair, UnorderedPairs};

use crate::structure::{
    GraphSignature, GraphWithEdgesSignature, LooplessGraphSignature, UndirectedGraphSignature,
};

/// An undirected cycle graph with n vertices arranged in a cycle.
/// Each vertex has exactly degree 2, forming a closed loop.
/// Also known as a circular graph or n-gon.
#[derive(Debug, Clone)]
pub struct UndirectedCycleGraph<Vertices: SetSignature> {
    vertices: Vertices,
    n: usize,
}

impl<Vertices: SetSignature> UndirectedCycleGraph<Vertices> {
    pub fn new(vertices: Vertices, n: usize) -> Result<Self, String> {
        if n < 3 {
            return Err("Cycle graphs must have at least 3 vertices for simple graphs".to_string());
        }

        Ok(Self { vertices, n })
    }

    pub fn is_bipartite(&self) -> bool {
        self.n.is_multiple_of(2)
    }

    pub fn girth(&self) -> usize {
        self.n
    }

    pub fn diameter(&self) -> usize {
        self.n / 2
    }

    pub fn chromatic_number(&self) -> usize {
        if self.n.is_multiple_of(2) { 2 } else { 3 }
    }
}

// Specialized implementation for EnumeratedFiniteSetStructure
use algebraeon_sets::structure::{EnumeratedFiniteSetStructure, EqSignature};

impl GraphSignature for UndirectedCycleGraph<EnumeratedFiniteSetStructure> {
    type Vertices = EnumeratedFiniteSetStructure;

    fn has_directed_edge(
        &self,
        source: &<Self::Vertices as SetSignature>::Set,
        target: &<Self::Vertices as SetSignature>::Set,
    ) -> Result<(), String> {
        if let Err(e) = self.vertices.validate_element(source) {
            return Err(format!("Source is not an element of Vertices: {e}"));
        }
        if let Err(e) = self.vertices.validate_element(target) {
            return Err(format!("Target is not an element of Vertices: {e}"));
        }

        if self.vertices.equal(source, target) {
            return Err("Cycle graphs do not contain self-loops".to_string());
        }

        let u = *source;
        let v = *target;

        if u >= self.n || v >= self.n {
            return Err("Vertex index out of range for cycle graph".to_string());
        }

        let adjacent = (u + 1) % self.n == v || (v + 1) % self.n == u;

        if adjacent {
            Ok(())
        } else {
            Err(format!(
                "No edge between vertices {u} and {v} in cycle graph"
            ))
        }
    }
}

impl LooplessGraphSignature for UndirectedCycleGraph<EnumeratedFiniteSetStructure> {}

impl UndirectedGraphSignature for UndirectedCycleGraph<EnumeratedFiniteSetStructure> {}

impl GraphWithEdgesSignature for UndirectedCycleGraph<EnumeratedFiniteSetStructure> {
    type Edges = UnorderedPairs<EnumeratedFiniteSetStructure>;

    fn endpoints(
        &self,
        edge: &<Self::Edges as SetSignature>::Set,
    ) -> UnorderedPair<<Self::Vertices as SetSignature>::Set> {
        edge.clone()
    }
}

/// A directed cycle graph with n vertices arranged in a cycle.
/// All edges are oriented in the same direction around the cycle.
/// Each vertex has in-degree 1 and out-degree 1.
#[derive(Debug, Clone)]
pub struct DirectedCycleGraph<Vertices: SetSignature> {
    vertices: Vertices,
    n: usize,
}

impl<Vertices: SetSignature> DirectedCycleGraph<Vertices> {
    pub fn new(vertices: Vertices, n: usize) -> Result<Self, String> {
        if n < 3 {
            return Err("Cycle graphs must have at least 3 vertices for simple graphs".to_string());
        }

        Ok(Self { vertices, n })
    }

    pub fn is_strongly_connected(&self) -> bool {
        true
    }

    pub fn girth(&self) -> usize {
        self.n
    }

    pub fn diameter(&self) -> usize {
        self.n - 1
    }
}

impl GraphSignature for DirectedCycleGraph<EnumeratedFiniteSetStructure> {
    type Vertices = EnumeratedFiniteSetStructure;

    fn has_directed_edge(
        &self,
        source: &<Self::Vertices as SetSignature>::Set,
        target: &<Self::Vertices as SetSignature>::Set,
    ) -> Result<(), String> {
        if let Err(e) = self.vertices.validate_element(source) {
            return Err(format!("Source is not an element of Vertices: {e}"));
        }
        if let Err(e) = self.vertices.validate_element(target) {
            return Err(format!("Target is not an element of Vertices: {e}"));
        }

        if self.vertices.equal(source, target) {
            return Err("Cycle graphs do not contain self-loops".to_string());
        }

        let u = *source;
        let v = *target;

        if u >= self.n || v >= self.n {
            return Err("Vertex index out of range for cycle graph".to_string());
        }

        let has_edge = (u + 1) % self.n == v;

        if has_edge {
            Ok(())
        } else {
            Err(format!(
                "No directed edge from vertex {u} to vertex {v} in directed cycle graph"
            ))
        }
    }
}

impl LooplessGraphSignature for DirectedCycleGraph<EnumeratedFiniteSetStructure> {}

#[cfg(test)]
mod tests {
    use algebraeon_sets::structure::EnumeratedFiniteSetStructure;

    use crate::{
        examples::cycle::{DirectedCycleGraph, UndirectedCycleGraph},
        structure::GraphSignature,
    };

    #[test]
    fn test_undirected_cycle_creation() {
        let fin5 = EnumeratedFiniteSetStructure::new(5);
        let c5 = UndirectedCycleGraph::new(fin5, 5);
        assert!(c5.is_ok());

        let c5 = c5.unwrap();
        assert_eq!(c5.n, 5);
        assert!(!c5.is_bipartite()); // C5 is odd, not bipartite
        assert_eq!(c5.girth(), 5);
        assert_eq!(c5.diameter(), 2);
        assert_eq!(c5.chromatic_number(), 3);
    }

    #[test]
    fn test_undirected_cycle_c4_properties() {
        let fin4 = EnumeratedFiniteSetStructure::new(4);
        let c4 = UndirectedCycleGraph::new(fin4, 4).unwrap();

        assert!(c4.is_bipartite()); // C4 is even, bipartite
        assert_eq!(c4.girth(), 4);
        assert_eq!(c4.diameter(), 2);
        assert_eq!(c4.chromatic_number(), 2);
    }

    #[test]
    fn test_undirected_cycle_c3_properties() {
        let fin3 = EnumeratedFiniteSetStructure::new(3);
        let c3 = UndirectedCycleGraph::new(fin3, 3).unwrap();

        assert!(!c3.is_bipartite()); // C3 is odd, not bipartite
        assert_eq!(c3.girth(), 3);
        assert_eq!(c3.diameter(), 1);
        assert_eq!(c3.chromatic_number(), 3);
    }

    #[test]
    fn test_undirected_cycle_invalid_size() {
        let fin2 = EnumeratedFiniteSetStructure::new(2);
        let c2 = UndirectedCycleGraph::new(fin2, 2);
        assert!(c2.is_err());

        let fin1 = EnumeratedFiniteSetStructure::new(1);
        let c1 = UndirectedCycleGraph::new(fin1, 1);
        assert!(c1.is_err());
    }

    #[test]
    fn test_undirected_cycle_adjacency() {
        let fin5 = EnumeratedFiniteSetStructure::new(5);
        let c5 = UndirectedCycleGraph::new(fin5.clone(), 5).unwrap();

        // In C5, vertices should be connected in cycle: 0-1-2-3-4-0
        // Test valid edges
        assert!(c5.has_directed_edge(&0, &1).is_ok()); // 0 -> 1
        assert!(c5.has_directed_edge(&1, &0).is_ok()); // 1 -> 0 (undirected)
        assert!(c5.has_directed_edge(&1, &2).is_ok()); // 1 -> 2
        assert!(c5.has_directed_edge(&2, &3).is_ok()); // 2 -> 3
        assert!(c5.has_directed_edge(&3, &4).is_ok()); // 3 -> 4
        assert!(c5.has_directed_edge(&4, &0).is_ok()); // 4 -> 0 (wraps around)
        assert!(c5.has_directed_edge(&0, &4).is_ok()); // 0 -> 4 (undirected)

        // Test invalid edges (non-adjacent vertices)
        assert!(c5.has_directed_edge(&0, &2).is_err()); // 0 and 2 not adjacent
        assert!(c5.has_directed_edge(&1, &3).is_err()); // 1 and 3 not adjacent
        assert!(c5.has_directed_edge(&2, &4).is_err()); // 2 and 4 not adjacent

        // Test self-loops (should be rejected)
        assert!(c5.has_directed_edge(&0, &0).is_err());
        assert!(c5.has_directed_edge(&1, &1).is_err());
    }

    #[test]
    fn test_undirected_cycle_c4_adjacency() {
        let fin4 = EnumeratedFiniteSetStructure::new(4);
        let c4 = UndirectedCycleGraph::new(fin4.clone(), 4).unwrap();

        // In C4, vertices should be connected: 0-1-2-3-0
        assert!(c4.has_directed_edge(&0, &1).is_ok());
        assert!(c4.has_directed_edge(&1, &2).is_ok());
        assert!(c4.has_directed_edge(&2, &3).is_ok());
        assert!(c4.has_directed_edge(&3, &0).is_ok());

        // Non-adjacent vertices should not have edges
        assert!(c4.has_directed_edge(&0, &2).is_err()); // Opposite vertices
        assert!(c4.has_directed_edge(&1, &3).is_err()); // Opposite vertices
    }

    #[test]
    fn test_directed_cycle_creation() {
        let fin5 = EnumeratedFiniteSetStructure::new(5);
        let dc5 = DirectedCycleGraph::new(fin5, 5);
        assert!(dc5.is_ok());

        let dc5 = dc5.unwrap();
        assert_eq!(dc5.n, 5);
        assert!(dc5.is_strongly_connected());
        assert_eq!(dc5.girth(), 5);
        assert_eq!(dc5.diameter(), 4);
    }

    #[test]
    fn test_directed_cycle_adjacency() {
        let fin5 = EnumeratedFiniteSetStructure::new(5);
        let dc5 = DirectedCycleGraph::new(fin5.clone(), 5).unwrap();

        // In directed C5, edges go: 0->1->2->3->4->0
        assert!(dc5.has_directed_edge(&0, &1).is_ok()); // 0 -> 1
        assert!(dc5.has_directed_edge(&1, &2).is_ok()); // 1 -> 2
        assert!(dc5.has_directed_edge(&2, &3).is_ok()); // 2 -> 3
        assert!(dc5.has_directed_edge(&3, &4).is_ok()); // 3 -> 4
        assert!(dc5.has_directed_edge(&4, &0).is_ok()); // 4 -> 0 (wraps around)

        // Test reverse directions (should fail in directed graph)
        assert!(dc5.has_directed_edge(&1, &0).is_err()); // 1 -> 0 (wrong direction)
        assert!(dc5.has_directed_edge(&2, &1).is_err()); // 2 -> 1 (wrong direction)
        assert!(dc5.has_directed_edge(&0, &4).is_err()); // 0 -> 4 (wrong direction)

        // Test non-adjacent vertices
        assert!(dc5.has_directed_edge(&0, &2).is_err()); // 0 -> 2 (skip vertex)
        assert!(dc5.has_directed_edge(&1, &3).is_err()); // 1 -> 3 (skip vertex)

        // Test self-loops
        assert!(dc5.has_directed_edge(&0, &0).is_err());
        assert!(dc5.has_directed_edge(&1, &1).is_err());
    }

    #[test]
    fn test_directed_cycle_properties() {
        let fin6 = EnumeratedFiniteSetStructure::new(6);
        let dc6 = DirectedCycleGraph::new(fin6, 6).unwrap();

        assert!(dc6.is_strongly_connected());
        assert_eq!(dc6.girth(), 6);
        assert_eq!(dc6.diameter(), 5);
    }

    #[test]
    fn test_directed_cycle_invalid_size() {
        let fin2 = EnumeratedFiniteSetStructure::new(2);
        let dc2 = DirectedCycleGraph::new(fin2, 2);
        assert!(dc2.is_err());

        let fin1 = EnumeratedFiniteSetStructure::new(1);
        let dc1 = DirectedCycleGraph::new(fin1, 1);
        assert!(dc1.is_err());
    }

    #[test]
    fn test_undirected_vs_directed_cycle_differences() {
        let fin4 = EnumeratedFiniteSetStructure::new(4);
        let uc4 = UndirectedCycleGraph::new(fin4.clone(), 4).unwrap();
        let dc4 = DirectedCycleGraph::new(fin4, 4).unwrap();

        // Both should have same girth
        assert_eq!(uc4.girth(), dc4.girth());

        // But different diameters
        assert_eq!(uc4.diameter(), 2); // Undirected: can go either way
        assert_eq!(dc4.diameter(), 3); // Directed: must follow direction

        // Test that undirected allows both directions
        assert!(uc4.has_directed_edge(&0, &1).is_ok());
        assert!(uc4.has_directed_edge(&1, &0).is_ok());

        // But directed only allows one direction
        assert!(dc4.has_directed_edge(&0, &1).is_ok());
        assert!(dc4.has_directed_edge(&1, &0).is_err());
    }

    #[test]
    fn test_cycle_graph_invalid_vertices() {
        let fin5 = EnumeratedFiniteSetStructure::new(5);
        let c5 = UndirectedCycleGraph::new(fin5.clone(), 5).unwrap();
        let dc5 = DirectedCycleGraph::new(fin5, 5).unwrap();

        // Test vertices out of range
        assert!(c5.has_directed_edge(&5, &0).is_err());
        assert!(c5.has_directed_edge(&0, &6).is_err());
        assert!(dc5.has_directed_edge(&5, &0).is_err());
        assert!(dc5.has_directed_edge(&0, &6).is_err());
    }

    #[test]
    fn test_cycle_properties_comprehensive() {
        // Test C3 (triangle)
        let fin3 = EnumeratedFiniteSetStructure::new(3);
        let c3 = UndirectedCycleGraph::new(fin3, 3).unwrap();
        assert_eq!(c3.diameter(), 1);
        assert_eq!(c3.chromatic_number(), 3);
        assert!(!c3.is_bipartite());

        // Test C6 (even cycle)
        let fin6 = EnumeratedFiniteSetStructure::new(6);
        let c6 = UndirectedCycleGraph::new(fin6.clone(), 6).unwrap();
        let dc6 = DirectedCycleGraph::new(fin6, 6).unwrap();

        assert_eq!(c6.diameter(), 3);
        assert_eq!(c6.chromatic_number(), 2);
        assert!(c6.is_bipartite());

        assert_eq!(dc6.diameter(), 5);
        assert!(dc6.is_strongly_connected());
    }

    #[test]
    fn test_even_vs_odd_cycle_properties() {
        // Even cycles are bipartite with chromatic number 2
        for n in [4, 6, 8, 10] {
            let fin_n = EnumeratedFiniteSetStructure::new(n);
            let c_n = UndirectedCycleGraph::new(fin_n, n).unwrap();
            assert!(c_n.is_bipartite(), "C{n} should be bipartite");
            assert_eq!(
                c_n.chromatic_number(),
                2,
                "C{n} should have chromatic number 2"
            );
        }

        // Odd cycles are not bipartite with chromatic number 3
        for n in [3, 5, 7, 9] {
            let fin_n = EnumeratedFiniteSetStructure::new(n);
            let c_n = UndirectedCycleGraph::new(fin_n, n).unwrap();
            assert!(!c_n.is_bipartite(), "C{n} should not be bipartite");
            assert_eq!(
                c_n.chromatic_number(),
                3,
                "C{n} should have chromatic number 3"
            );
        }
    }
}
