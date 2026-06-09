mod sets;
mod signatures;

pub use sets::*;
pub use signatures::*;

use std::{borrow::Borrow, fmt::Debug};

#[macro_export]
macro_rules! make_maybe_trait {
    ($name:ident) => {
        paste! {
            pub trait [<Maybe $name Signature>]: $crate::SetSignature {
                type [<$name Structure>]: [<$name Signature>] <Elem = Self::Elem>;

                #[allow(clippy::result_unit_err)]
                fn [<$name:snake _structure>](
                    &self
                ) -> Result<Self::[<$name Structure>], ()>;
            }
        }
    };
}
