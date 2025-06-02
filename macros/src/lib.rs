extern crate proc_macro;

use proc_macro::TokenStream;
use quote::quote;
use syn::{Attribute, DeriveInput, Ident, parse_macro_input};

fn has_option(attrs: &[Attribute], option_name: &str) -> bool {
    for attr in attrs
        .iter()
        .filter(|a| a.path().is_ident("canonical_structure"))
    {
        let mut found = false;

        // `parse_nested_meta` lets us walk through the arguments in #[canonical_structure(...)]
        let _ = attr.parse_nested_meta(|meta| {
            if meta.path.is_ident(option_name) {
                found = true;
            }
            // Continue parsing
            Ok(())
        });

        if found {
            return true;
        }
    }
    false
}

#[proc_macro_derive(CanonicalStructure, attributes(canonical_structure))]
pub fn derive_newtype(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);

    let name = input.ident;
    let vis = input.vis;
    let newtype_name = Ident::new(&format!("{name}CanonicalStructure"), name.span());

    let has_eq = has_option(&input.attrs, "eq");
    let has_ord = has_option(&input.attrs, "ord");

    let impl_eq_signature = if has_eq {
        quote! {
            impl EqSignature for #newtype_name
                where #name: Eq
            {
                fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
                    a == b
                }
            }
        }
    } else {
        quote! {}
    };

    let impl_ord_signature = if has_ord {
        quote! {
            impl OrdSignature for #newtype_name
                where #name: Ord
            {
                fn cmp(&self, a: &Self::Set, b: &Self::Set) -> std::cmp::Ordering {
                    a.cmp(b)
                }

                fn sort<S: std::borrow::Borrow<Self::Set>>(&self, mut a: Vec<S>) -> Vec<S> {
                    a.sort_unstable_by(|x, y| x.borrow().cmp(y.borrow()));
                    a
                }
            }
        }
    } else {
        quote! {}
    };

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

            fn is_element(&self, _x : &Self::Set) -> Result<(), String> {
                Ok(())
            }
        }

        #impl_eq_signature
        #impl_ord_signature

        impl MetaType for #name {
            type Signature = #newtype_name;

            fn structure() -> Self::Signature {
                #newtype_name::new()
            }
        }
    };

    TokenStream::from(expanded)
}
