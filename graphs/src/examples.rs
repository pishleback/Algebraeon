use algebraeon_sets::structure::{EqSignature, Pairs, SetSignature, UnorderedPair, UnorderedPairs};

use crate::structure::{
    GraphSignature, GraphWithEdgesSignature, LooplessGraphSignature, UndirectedGraphSignature,
};

struct CompleteDirectedGraph<Vertices: SetSignature> {
    vertices: Vertices,

    pairs_of_vertices: UnorderedPairs<Vertices>,
}

impl<Vertices: SetSignature + EqSignature> GraphSignature for CompleteDirectedGraph<Vertices> {
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

impl<Vertices: SetSignature + EqSignature> LooplessGraphSignature
    for CompleteDirectedGraph<Vertices>
{
}

impl<Vertices: SetSignature + EqSignature> GraphWithEdgesSignature
    for CompleteDirectedGraph<Vertices>
{
    type Edges = Pairs<Vertices>;

    fn endpoints(
        &self,
        edge: &<Self::Edges as SetSignature>::Set,
    ) -> UnorderedPair<<Self::Vertices as SetSignature>::Set> {
        self.pairs_of_vertices.new_pair(&edge.0, &edge.1).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use algebraeon_sets::structure::{
        EnumeratedFiniteSetStructure, EqSignature, Pairs, UnorderedPairs,
    };

    use crate::{
        examples::CompleteDirectedGraph,
        structure::{GraphSignature, GraphWithEdgesSignature},
    };

    #[test]
    fn test_k5() {
        let fin5 = EnumeratedFiniteSetStructure::new(5);
        let fin5_unordered_pairs = UnorderedPairs::new(fin5.clone());
        let k5 = CompleteDirectedGraph {
            vertices: fin5.clone(),
            pairs_of_vertices: fin5_unordered_pairs.clone(),
        };

        assert!(k5.has_directed_edge(&1, &2).is_ok());
        assert!(k5.has_directed_edge(&1, &1).is_err());
        assert!(k5.has_directed_edge(&5, &1).is_err());

        let fin5_pairs = Pairs::new(fin5.clone());

        let e1 = fin5_pairs.clone().new_pair(1.clone(), 2.clone()).unwrap();
        let e2 = fin5_pairs.clone().new_pair(2.clone(), 1.clone()).unwrap();
        assert!(fin5_unordered_pairs.equal(&k5.endpoints(&e1), &k5.endpoints(&e2)));
    }
}
