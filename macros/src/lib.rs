extern crate proc_macro;

use proc_macro::TokenStream;
use quote::quote;
use syn::{DeriveInput, Ident, parse_macro_input};

#[proc_macro_derive(CanonicalStructure)]
pub fn derive_newtype(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);

    let name = input.ident;
    let vis = input.vis;
    let newtype_name = Ident::new(&format!("{name}CanonicalStructure"), name.span());

    let expanded = quote! {
        #[derive(Debug, Clone, PartialEq, Eq)]
        #vis struct #newtype_name {}

        impl #newtype_name {
            fn new() -> Self {
                Self {}
            }
        }

        impl Signature for #newtype_name {}

        impl SetSignature for #newtype_name {
            type Set = #name;

            fn is_element(&self, _x : &Self::Set) -> bool {
                true
            }
        }

        impl EqSignature for #newtype_name
            where #name: Eq
        {
            fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
                a == b
            }
        }

        impl MetaType for #name {
            type Signature = #newtype_name;

            fn structure() -> Self::Signature {
                #newtype_name::new()
            }
        }
    };

    TokenStream::from(expanded)
}
