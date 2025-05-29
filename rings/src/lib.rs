#![allow(dead_code)]
#![allow(
    clippy::uninlined_format_args,
    clippy::to_string_in_format_args,
    clippy::println_empty_string,
    clippy::manual_assert,
    clippy::manual_while_let_some,
    clippy::assertions_on_constants,
    clippy::must_use_candidate,
    clippy::return_self_not_must_use,
    clippy::missing_panics_doc,
    clippy::missing_errors_doc,
    clippy::needless_pass_by_value,
    clippy::ptr_arg,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::cast_possible_truncation,
    clippy::unnecessary_cast,
    clippy::needless_lifetimes,
    clippy::result_unit_err,
    clippy::needless_range_loop,
    clippy::too_many_lines,
    clippy::similar_names,
    clippy::many_single_char_names,
    clippy::wrong_self_convention,
    clippy::from_over_into,
    clippy::match_wild_err_arm,
    clippy::derived_hash_with_manual_eq,
    clippy::non_canonical_partial_ord_impl,
)]
#![allow(
    clippy::single_match,
    clippy::match_same_arms,
    clippy::wildcard_imports,
    clippy::needless_borrow,
)]

pub mod linear;
pub mod parsing;
pub mod polynomial;
pub mod rings;
pub mod structure;
