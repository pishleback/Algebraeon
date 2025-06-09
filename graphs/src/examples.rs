use algebraeon_sets::structure::{EqSignature, SetSignature};

use crate::structure::{GraphSignature, LooplessGraphSignature, UndirectedGraphSignature};

struct CompleteGraph<Vertices: SetSignature> {
    vertices: Vertices,
}

impl<Vertices: SetSignature + EqSignature> GraphSignature for CompleteGraph<Vertices> {
    type Vertices = Vertices;

    fn has_directed_edge(
        &self,
        source: &Vertices::Set,
        target: &Vertices::Set,
    ) -> Result<(), String> {
        if let Err(e) = self.vertices.is_element(source) {
            return Err(format!("Source is not an element of Vertices: {}", e));
        }
        if let Err(e) = self.vertices.is_element(target) {
            return Err(format!("Target is not an element of Vertices: {}", e));
        }
        if self.vertices.equal(source, target) {
            return Err("Complete graphs do not contain loops.".to_string());
        }
        Ok(())
    }
}

impl<Vertices: SetSignature + EqSignature> LooplessGraphSignature for CompleteGraph<Vertices> {}

impl<Vertices: SetSignature + EqSignature> UndirectedGraphSignature for CompleteGraph<Vertices> {}

#[cfg(test)]
mod tests {
    use algebraeon_sets::structure::EnumeratedFiniteSetStructure;

    use crate::{examples::CompleteGraph, structure::GraphSignature};

    #[test]
    fn test_k5() {
        let fin5 = EnumeratedFiniteSetStructure::new(5);
        let k5 = CompleteGraph { vertices: fin5 };

        assert!(k5.has_directed_edge(&1, &2).is_ok());
        assert!(k5.has_directed_edge(&1, &1).is_err());
        assert!(k5.has_directed_edge(&5, &1).is_err());
    }
}
