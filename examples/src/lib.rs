#![allow(non_snake_case)]
#![allow(clippy::doc_lazy_continuation)] // because it parses the guide weirdly
#![doc = include_str!("../../README.md")]
include!(concat!(env!("OUT_DIR"), "/generated_docs.rs"));

#[cfg(test)]
mod tests {
    use algebraeon::{rings::finite_fields::modulo::Modulo, sets::structure::*};

    #[test]
    fn enumerate_functions_finite_sets() {
        let fns = Functions::new(Modulo::<3>::structure(), Modulo::<4>::structure());
        assert_eq!(fns.list_all_elements().len(), 64);
    }
}
