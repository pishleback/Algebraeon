extern crate proc_macro;

use proc_macro::TokenStream;
use proc_macro2::Span;
use quote::quote;
use syn::visit_mut::VisitMut;
use syn::{
    Attribute, DeriveInput, Error, FnArg, Ident, ItemTrait, PatIdent, Receiver, TraitItem,
    TraitItemFn, parse_macro_input,
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

/// Decorate fns of a structure trait decorated with `signature_meta_trait` to omit the fn from the meta structure trait.
#[proc_macro_attribute]
pub fn skip_meta(_attr: TokenStream, item: TokenStream) -> TokenStream {
    item
}

/// Decorate a structure trait with this to auto-generate a meta structure trait.
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

    for item in &trait_item.items {
        if let TraitItem::Fn(TraitItemFn { attrs, sig, .. }) = item {
            if attrs.iter().any(|attr| attr.path().is_ident("skip_meta")) {
                continue;
            }

            let mut meta_sig = sig.clone();
            // Check the first argument is self, &self, or &mut self
            if let Some(first_arg) = meta_sig.inputs.first() {
                match first_arg {
                    FnArg::Receiver(_) => {
                        meta_sig.inputs = meta_sig.inputs.into_iter().skip(1).collect();
                        ReplaceSelfSetSignature {
                            sig_trait_ident: sig_trait_ident.clone(),
                        }
                        .visit_signature_mut(&mut meta_sig);

                        let ident = meta_sig.ident.clone();

                        let mut meta_args = Vec::new();
                        #[allow(clippy::never_loop)]
                        for arg in &mut meta_sig.inputs {
                            match arg {
                                FnArg::Typed(pat_type) => match pat_type.pat.as_mut() {
                                    syn::Pat::Ident(pat_ident) => {
                                        pat_ident.mutability = None;
                                        meta_args.push(pat_ident.clone());
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

                        if let Some(first) = sig.inputs.iter().nth(1) {
                            match first {
                                FnArg::Receiver(_) => {}
                                FnArg::Typed(pat_type) => match pat_type.ty.as_ref() {
                                    syn::Type::Reference(type_reference) => {
                                        if let syn::Type::Path(type_path) =
                                            type_reference.elem.as_ref()
                                            && is_self_set(type_path)
                                        {
                                            // if the first argument is `a: &Self::Set` then replace it with `&self` in the meta type
                                            // if the first argument is `a: &mut Self::Set` then replace it with `&mut self` in the meta type
                                            meta_args[0] = PatIdent {
                                                attrs: vec![],
                                                by_ref: None,
                                                mutability: None,
                                                ident: Ident::new("self", Span::call_site()),
                                                subpat: None,
                                            };
                                            meta_sig.inputs[0] = FnArg::Receiver(Receiver {
                                                attrs: vec![],
                                                reference: Some((
                                                    syn::token::And {
                                                        spans: [Span::call_site()],
                                                    },
                                                    None,
                                                )),
                                                mutability: type_reference.mutability,
                                                self_token: syn::token::SelfValue {
                                                    span: Span::call_site(),
                                                },
                                                colon_token: None,
                                                ty: Box::new(syn::Type::Reference(
                                                    syn::TypeReference {
                                                        and_token: syn::token::And {
                                                            spans: [Span::call_site()],
                                                        },
                                                        lifetime: None,
                                                        mutability: type_reference.mutability,
                                                        elem: Box::new(syn::Type::Path(
                                                            syn::TypePath {
                                                                qself: None,
                                                                path: syn::Path::from(Ident::new(
                                                                    "Self",
                                                                    Span::call_site(),
                                                                )),
                                                            },
                                                        )),
                                                    },
                                                )),
                                            });
                                        }
                                    }
                                    syn::Type::Path(type_path) => {
                                        // if the first argument is `a: Self::Set` then replace it with `self` in the meta type (TODO)
                                        if is_self_set(type_path) {
                                            meta_args[0] = PatIdent {
                                                attrs: vec![],
                                                by_ref: None,
                                                mutability: None,
                                                ident: Ident::new("self", Span::call_site()),
                                                subpat: None,
                                            };
                                            meta_sig.inputs[0] = FnArg::Receiver(Receiver {
                                                attrs: vec![],
                                                reference: None,
                                                mutability: None,
                                                self_token: syn::token::SelfValue {
                                                    span: Span::call_site(),
                                                },
                                                colon_token: None,
                                                ty: Box::new(syn::Type::Path(syn::TypePath {
                                                    qself: None,
                                                    path: syn::Path::from(Ident::new(
                                                        "Self",
                                                        Span::call_site(),
                                                    )),
                                                })),
                                            });
                                        }
                                    }
                                    _ => {}
                                },
                            }
                        }

                        meta_methods.push(quote! {
                            #(#attrs)*
                            #meta_sig {
                                Self::structure().#ident(#(#meta_args),*)
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

struct ReplaceSelfSetSignature {
    #[allow(unused)]
    sig_trait_ident: Ident,
}
impl VisitMut for ReplaceSelfSetSignature {
    fn visit_type_path_mut(&mut self, ty: &mut syn::TypePath) {
        syn::visit_mut::visit_type_path_mut(self, ty);
        if is_self_set(ty) {
            *ty = syn::parse_quote!(Self);
        } else if is_self(ty) {
            *ty = syn::parse_quote!(Self::Signature);
        }
    }
}

fn is_self_set(ty: &syn::TypePath) -> bool {
    ty.qself.is_none()
        && ty.path.segments.len() == 2
        && ty.path.segments[0].ident == "Self"
        && ty.path.segments[1].ident == "Set"
        && ty.path.segments[1].arguments.is_empty()
}

fn is_self(ty: &syn::TypePath) -> bool {
    ty.qself.is_none()
        && ty.path.segments.len() == 1
        && ty.path.segments[0].ident == "Self"
        && ty.path.segments[0].arguments.is_empty()
}
