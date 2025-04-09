extern crate proc_macro;

use proc_macro::TokenStream;
use quote::quote;
use syn::{DeriveInput, Ident, parse_macro_input};

#[proc_macro_derive(CanonicalStructure)]
pub fn derive_newtype(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);

    let name = input.ident;
    let newtype_name = Ident::new(&format!("{}CanonicalStructure", name), name.span());

    let expanded = quote! {
        #[derive(Debug, Clone, PartialEq, Eq)]
        pub struct #newtype_name {}

        impl #newtype_name {
            fn new() -> Self {
                Self {}
            }
        }

        impl Structure for #newtype_name {}

        impl SetStructure for #newtype_name {
            type Set = #name;
        }

        impl PartialEqStructure for #newtype_name
            where #name: Eq
        {
            fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
                a == b
            }
        }

        impl EqStructure for #newtype_name {}

        impl MetaType for #name {
            type Structure = #newtype_name;

            fn structure() -> Self::Structure {
                #newtype_name::new()
            }
        }
    };

    TokenStream::from(expanded)
}
