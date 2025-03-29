//! Contains combinatorial counting and enumeration algorithms.

mod combinations;
mod partitions;

pub use combinations::LexicographicCombinationsWithRemovals;
pub use combinations::combinations;
pub use combinations::subsets;
pub use partitions::partitions;
pub use partitions::partitions_sized;
pub use partitions::partitions_sized_zero;
pub use partitions::predicated_partitions;
pub use partitions::predicated_partitions_sized;
pub use partitions::predicated_partitions_sized_zero;
