//! Contains combinatorial counting and enumeration algorithms.

mod number_compositions;
mod number_partitions;
mod set_partitions;
mod subsets;
mod twelvefold_way;

pub use number_compositions::compositions;
pub use number_compositions::compositions_sized;
pub use number_compositions::compositions_sized_zero;
pub use number_partitions::num_partitions;
pub use number_partitions::num_partitions_part_pool;
pub use number_partitions::num_partitions_predicated;
pub use number_partitions::num_partitions_sized;
pub use number_partitions::num_partitions_sized_predicated;
pub use number_partitions::num_partitions_sized_zero;
pub use number_partitions::num_partitions_sized_zero_predicated;
pub use set_partitions::Partition;
pub use set_partitions::set_compositions_eq;
pub use set_partitions::set_partitions_eq;
pub use set_partitions::set_partitions_ge;
pub use set_partitions::set_partitions_le;
pub use set_partitions::set_partitions_range;
pub use subsets::LexicographicSubsetsWithRemovals;
pub use subsets::all_subsets;
pub use subsets::subsets;
pub use subsets::subsets_of_vec;
// pub use twelvefold_way::{FunctionType, TwelvefoldWay};
