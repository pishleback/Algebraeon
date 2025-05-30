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
    clippy::doc_markdown,
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
    clippy::items_after_statements,
    clippy::type_complexity,
    clippy::too_many_arguments,
    clippy::module_inception,
    clippy::implicit_hasher,
    clippy::unnecessary_box_returns,
    clippy::no_effect_underscore_binding,
    clippy::multiple_bound_locations
)]
#![allow(
    clippy::wildcard_imports,
    clippy::redundant_closure,
    clippy::unnecessary_wraps,
    clippy::useless_conversion,
    clippy::nonminimal_bool,
    clippy::while_let_loop,
    clippy::let_and_return,
    clippy::op_ref,
    clippy::absurd_extreme_comparisons,
    clippy::assign_op_pattern,
    clippy::redundant_closure_for_method_calls,
    clippy::needless_return,
    clippy::needless_borrows_for_generic_args,
    clippy::filter_next,
    clippy::manual_find,
    clippy::if_not_else,
    clippy::default_constructed_unit_structs,
    clippy::into_iter_on_ref
)]

pub mod linear;
pub mod parsing;
pub mod polynomial;
pub mod rings;
pub mod structure;
