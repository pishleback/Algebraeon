use algebraeon_sets::structure::{
    EnumeratedFiniteSetStructure, EqSignature, SetSignature, UnorderedPair, UnorderedPairs,
};

use crate::structure::{
    GraphSignature, GraphWithEdgesSignature, LooplessGraphSignature, UndirectedGraphSignature,
};

/// A wheel graph `W_n` with a single universal center vertex
/// (index `0`) joined to all vertices of an outer cycle on `n - 1` vertices.
/// Requires `n >= 4`, where `W_4` is the complete graph `K_4`.
#[derive(Debug, Clone)]
pub struct WheelGraph<Vertices: SetSignature> {
    vertices: Vertices,
    n: usize,
}

impl<Vertices: SetSignature> WheelGraph<Vertices> {
    pub fn new(vertices: Vertices, n: usize) -> Result<Self, String> {
        if n < 4 {
            return Err("Wheel graphs require at least 4 vertices".to_string());
        }

        Ok(Self { vertices, n })
    }

    pub fn vertex_count(&self) -> usize {
        self.n
    }

    pub fn edge_count(&self) -> usize {
        2 * (self.n - 1)
    }

    pub fn diameter(&self) -> usize {
        if self.n == 4 { 1 } else { 2 }
    }

    pub fn girth(&self) -> usize {
        3
    }

    pub fn chromatic_number(&self) -> usize {
        if self.n.is_multiple_of(2) { 4 } else { 3 }
    }

    pub fn is_hamiltonian(&self) -> bool {
        true
    }
}

impl GraphSignature for WheelGraph<EnumeratedFiniteSetStructure> {
    type Vertices = EnumeratedFiniteSetStructure;

    fn has_directed_edge(
        &self,
        source: &<Self::Vertices as SetSignature>::Set,
        target: &<Self::Vertices as SetSignature>::Set,
    ) -> Result<(), String> {
        if let Err(e) = self.vertices.is_element(source) {
            return Err(format!("Source is not an element of Vertices: {e}"));
        }
        if let Err(e) = self.vertices.is_element(target) {
            return Err(format!("Target is not an element of Vertices: {e}"));
        }

        if self.vertices.equal(source, target) {
            return Err("Wheel graphs do not contain self-loops".to_string());
        }

        let u = *source;
        let v = *target;

        if u >= self.n || v >= self.n {
            return Err("Vertex index out of range for wheel graph".to_string());
        }

        let center = 0;
        if u == center || v == center {
            return Ok(());
        }

        // Both vertices lie on the outer cycle (indices 1 .. n-1)
        let cycle_len = self.n - 1;
        let u_cycle = u - 1;
        let v_cycle = v - 1;

        let adjacent = (u_cycle + 1) % cycle_len == v_cycle || (v_cycle + 1) % cycle_len == u_cycle;

        if adjacent {
            Ok(())
        } else {
            Err(format!(
                "No edge between vertices {u} and {v} in wheel graph"
            ))
        }
    }
}

impl LooplessGraphSignature for WheelGraph<EnumeratedFiniteSetStructure> {}

impl UndirectedGraphSignature for WheelGraph<EnumeratedFiniteSetStructure> {}

impl GraphWithEdgesSignature for WheelGraph<EnumeratedFiniteSetStructure> {
    type Edges = UnorderedPairs<EnumeratedFiniteSetStructure>;

    fn endpoints(
        &self,
        edge: &<Self::Edges as SetSignature>::Set,
    ) -> UnorderedPair<<Self::Vertices as SetSignature>::Set> {
        edge.clone()
    }
}

#[cfg(test)]
mod tests {
    use algebraeon_sets::structure::EnumeratedFiniteSetStructure;

    use crate::{examples::WheelGraph, structure::GraphSignature};

    #[test]
    fn test_wheel_creation() {
        let vertices = EnumeratedFiniteSetStructure::new(6);
        let wheel = WheelGraph::new(vertices, 6).unwrap();

        assert_eq!(wheel.vertex_count(), 6);
        assert_eq!(wheel.edge_count(), 10); // 2 * (6 - 1)
        assert_eq!(wheel.diameter(), 2);
        assert_eq!(wheel.girth(), 3);
        assert_eq!(wheel.chromatic_number(), 4); // n even
        assert!(wheel.is_hamiltonian());
    }

    #[test]
    fn test_wheel_w4_properties() {
        let vertices = EnumeratedFiniteSetStructure::new(4);
        let wheel = WheelGraph::new(vertices, 4).unwrap();

        assert_eq!(wheel.diameter(), 1); // W4 is K4
        assert_eq!(wheel.chromatic_number(), 4);
        assert_eq!(wheel.edge_count(), 6);
    }

    #[test]
    fn test_wheel_adjacency() {
        let vertices = EnumeratedFiniteSetStructure::new(7);
        let wheel = WheelGraph::new(vertices, 7).unwrap();

        for v in 1..7 {
            assert!(
                wheel.has_directed_edge(&0, &v).is_ok(),
                "Expected spoke 0-{v}"
            );
            assert!(
                wheel.has_directed_edge(&v, &0).is_ok(),
                "Expected spoke {v}-0"
            );
        }

        let rim = [1, 2, 3, 4, 5, 6, 1];
        for edge in rim.windows(2) {
            let u = edge[0];
            let v = edge[1];
            assert!(
                wheel.has_directed_edge(&u, &v).is_ok(),
                "Expected rim edge {u}-{v}"
            );
            assert!(
                wheel.has_directed_edge(&v, &u).is_ok(),
                "Expected rim edge {v}-{u}"
            );
        }

        assert!(wheel.has_directed_edge(&1, &3).is_err());
        assert!(wheel.has_directed_edge(&2, &5).is_err());
    }

    #[test]
    fn test_wheel_invalid_cases() {
        let vertices = EnumeratedFiniteSetStructure::new(5);
        let wheel = WheelGraph::new(vertices.clone(), 5).unwrap();

        assert!(wheel.has_directed_edge(&2, &2).is_err());
        assert!(wheel.has_directed_edge(&5, &1).is_err());
        assert!(wheel.has_directed_edge(&1, &7).is_err());

        assert!(WheelGraph::new(vertices, 3).is_err());
    }
}
