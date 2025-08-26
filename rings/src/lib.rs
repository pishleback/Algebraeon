#![allow(dead_code)]
#![allow(
    clippy::uninlined_format_args,
    clippy::to_string_in_format_args,
    clippy::println_empty_string,
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
    clippy::multiple_bound_locations,
    clippy::wildcard_imports,
    clippy::manual_assert
)]

pub mod algebraic_number_field;
pub mod finite_fields;
pub mod integer;
pub mod isolated_algebraic;
pub mod matrix;
pub mod module;
pub mod natural;
pub mod num_theory;
pub mod parsing;
pub mod polynomial;
pub mod quaternion_algebras;
pub mod rational;
pub mod structure;
pub mod valuation;
