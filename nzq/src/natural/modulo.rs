// use super::*;

// #[derive(Debug)]
// pub enum InvertModuloError {
//     ZeroModulus,
//     NotInvertible,
// }

// pub trait InvertModulo<Modulus> {
//     type Output;

//     fn invert_modulo(self, modulus: Modulus) -> Result<Self::Output, InvertModuloError>;
// }

// impl InvertModulo<Natural> for Natural {
//     type Output = Natural;

//     fn invert_modulo(self, modulus: Natural) -> Result<Self::Output, InvertModuloError> {
//         if modulus == Natural::ZERO {
//             Err(InvertModuloError::ZeroModulus)
//         } else {
//             let
//         }
//     }
// }
