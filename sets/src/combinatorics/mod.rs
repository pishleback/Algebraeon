//! Contains combinatorial counting and enumeration algorithms.

mod combinations;
mod number_partitions;
mod set_partitions;
// TODO: subsets

pub use combinations::LexicographicCombinationsWithRemovals;
pub use combinations::combinations;
pub use combinations::subsets;
pub use number_partitions::partitions;
pub use number_partitions::partitions_sized;
pub use number_partitions::partitions_sized_zero;
pub use number_partitions::predicated_partitions;
pub use number_partitions::predicated_partitions_sized;
pub use number_partitions::predicated_partitions_sized_zero;
pub use set_partitions::Partition;
