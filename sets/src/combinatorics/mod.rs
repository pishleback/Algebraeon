//! Contains combinatorial counting and enumeration algorithms.

mod number_compositions;
mod number_partitions;
mod set_partitions;
mod subsets;
mod twelvefold_way;

pub use number_compositions::compositions;
pub use number_compositions::compositions_sized;
pub use number_compositions::compositions_sized_zero;
pub use number_partitions::partitions;
pub use number_partitions::partitions_sized;
pub use number_partitions::partitions_sized_zero;
pub use number_partitions::predicated_partitions;
pub use number_partitions::predicated_partitions_sized;
pub use number_partitions::predicated_partitions_sized_zero;
pub use set_partitions::Partition;
pub use subsets::LexicographicSubsetsWithRemovals;
pub use subsets::subsets;
pub use subsets::subsets_of_vec;
pub use twelvefold_way::{FunctionType, TwelvefoldWay};
