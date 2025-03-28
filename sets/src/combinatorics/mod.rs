//! Contains combinatorial counting and enumeration algorithms.

mod partition;
mod subset;

pub use partition::partitions;
pub use partition::partitions_sized;
pub use partition::partitions_sized_zero;
pub use partition::predicated_partitions;
pub use partition::predicated_partitions_sized;
pub use partition::predicated_partitions_sized_zero;
pub use subset::LexicographicSubsetsWithRemovals;
pub use subset::subsets;
pub use subset::subsets_vec;
