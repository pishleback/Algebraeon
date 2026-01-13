extern crate proc_macro;

use proc_macro::TokenStream;
use proc_macro2::Span;
use quote::quote;
use syn::visit_mut::VisitMut;
use syn::{
    Attribute, DeriveInput, Error, FnArg, Ident, ItemTrait, TraitItem, TraitItemFn,
    parse_macro_input,
};

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
    let has_partial_ord = has_option(&input.attrs, "partial_ord");
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

    let impl_partial_ord_signature = if has_partial_ord {
        quote! {
            impl PartialOrdSignature for #newtype_name
                where #name: Ord
            {
                fn partial_cmp(&self, a: &Self::Set, b: &Self::Set) -> Option<std::cmp::Ordering> {
                    Some(a.cmp(b))
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
        #impl_partial_ord_signature
        #impl_ord_signature

        impl MetaType for #name {
            type Signature = #newtype_name;

            fn structure() -> Self::Signature {
                #newtype_name::new()
            }
        }

        impl #name {
            pub fn structure_ref() -> &'static #newtype_name{
                static CELL: std::sync::OnceLock<#newtype_name> = std::sync::OnceLock::new();
                CELL.get_or_init(|| #newtype_name::new())
            }
        }
    };

    TokenStream::from(expanded)
}

#[proc_macro_attribute]
pub fn skip_meta(_attr: TokenStream, item: TokenStream) -> TokenStream {
    item
}

#[proc_macro_attribute]
pub fn signature_meta_trait(_args: TokenStream, input: TokenStream) -> TokenStream {
    let trait_item = parse_macro_input!(input as ItemTrait);

    let expanded = expand_meta_trait(&trait_item);

    quote! {
        #trait_item
        #expanded
    }
    .into()
}

/// Expand MetaTrait + impl
fn expand_meta_trait(trait_item: &ItemTrait) -> proc_macro2::TokenStream {
    let sig_trait_ident = &trait_item.ident;
    let meta_trait_ident = Ident::new(&format!("Meta{}", sig_trait_ident), Span::call_site());

    let mut meta_methods = Vec::new();
    // let mut impl_methods = Vec::new();

    for item in &trait_item.items {
        if let TraitItem::Fn(TraitItemFn { attrs, sig, .. }) = item {
            if attrs.iter().any(|attr| attr.path().is_ident("skip_meta")) {
                continue;
            }

            let mut sig = sig.clone();
            // Check the first argument
            if let Some(first_arg) = sig.inputs.first() {
                match first_arg {
                    FnArg::Receiver(_) => {
                        sig.inputs = sig.inputs.into_iter().skip(1).collect();
                        ReplaceSelfSet.visit_signature_mut(&mut sig);

                        let ident = sig.ident.clone();

                        let mut arg_names = Vec::new();
                        #[allow(clippy::never_loop)]
                        for arg in &sig.inputs {
                            match arg {
                                FnArg::Typed(pat_type) => match pat_type.pat.as_ref() {
                                    syn::Pat::Ident(pat_ident) => {
                                        arg_names.push(pat_ident.clone());
                                    }
                                    _ => {
                                        return Error::new_spanned(
                                        trait_item,
                                        "Invalid pattern in argument list. Must be a plain Ident.",
                                    )
                                    .to_compile_error();
                                    }
                                },
                                FnArg::Receiver(_) => {
                                    panic!();
                                }
                            }
                        }

                        meta_methods.push(quote! {
                            #(#attrs)*
                            #sig {
                                Self::structure().#ident(#(#arg_names),*)
                            }
                        });
                    }
                    FnArg::Typed(_) => {
                        // Not a method receiver
                    }
                }
            }
        }
    }

    quote! {
        pub trait #meta_trait_ident: MetaType
        where
            Self::Signature: #sig_trait_ident {

            #(#meta_methods)*
        }

        impl<T> #meta_trait_ident for T
        where
            T: MetaType,
            T::Signature: #sig_trait_ident,
        {
        }
    }
}

struct ReplaceSelfSet;

impl VisitMut for ReplaceSelfSet {
    fn visit_type_path_mut(&mut self, ty: &mut syn::TypePath) {
        syn::visit_mut::visit_type_path_mut(self, ty);

        if ty.qself.is_none()
            && ty.path.segments.len() == 2
            && ty.path.segments[0].ident == "Self"
            && ty.path.segments[1].ident == "Set"
            && ty.path.segments[1].arguments.is_empty()
        {
            *ty = syn::parse_quote!(Self);
        }
    }
}
