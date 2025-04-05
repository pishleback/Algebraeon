// #[derive(Debug, Clone, Copy, PartialEq, Eq)]
// pub enum FunctionType {
//     Any,
//     Injective,
//     Surjective,
// }

// #[derive(Debug, Clone, Copy, PartialEq, Eq)]
// pub enum TwelvefoldWay {
//     // N-sequence in X
//     Sequence,
//     // N-permutation in X
//     Permutation,
//     // composition of N using X subsets
//     SetComposition,
//     // n-multisubset of X
//     Multisubset,
//     // n-subset of X
//     Subset,
//     // composition of n with X terms
//     NatComposition,
//     // partition of N into <= x subsets
//     SetPartitionLe,
//     // partition of N into <= x elements
//     SetLe,
//     // partition of N into x subsets
//     SetPartitionEq,
//     // partition of n into <= x parts
//     NatPartitionLe,
//     // partition of n into <= x parts of size 1
//     NatLe,
//     // partition of n into x parts
//     NatPartitionEq,
// }

// impl TwelvefoldWay {
//     pub fn new(function_type: FunctionType, permute_n: bool, permute_x: bool) -> Self {
//         match (function_type, permute_n, permute_x) {
//             (FunctionType::Any, false, false) => Self::Sequence,
//             (FunctionType::Injective, false, false) => Self::Permutation,
//             (FunctionType::Surjective, false, false) => Self::SetComposition,
//             (FunctionType::Any, true, false) => Self::Multisubset,
//             (FunctionType::Injective, true, false) => Self::Subset,
//             (FunctionType::Surjective, true, false) => Self::NatComposition,
//             (FunctionType::Any, false, true) => Self::SetPartitionLe,
//             (FunctionType::Injective, false, true) => Self::SetLe,
//             (FunctionType::Surjective, false, true) => Self::SetPartitionEq,
//             (FunctionType::Any, true, true) => Self::NatPartitionLe,
//             (FunctionType::Injective, true, true) => Self::NatLe,
//             (FunctionType::Surjective, true, true) => Self::NatPartitionEq,
//         }
//     }
// }
