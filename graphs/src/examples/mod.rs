pub mod complete_directed;
pub mod complete_undirected;
pub mod cycle;
pub mod wheel;

pub use complete_directed::CompleteDirectedGraph;
pub use complete_undirected::CompleteUndirectedGraph;
pub use cycle::{DirectedCycleGraph, UndirectedCycleGraph};
pub use wheel::WheelGraph;
