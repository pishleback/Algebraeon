use std::rc::Rc;

use super::{isolated_roots::complex::ComplexAlgebraic, number_field::*};

#[derive(Debug, Clone)]
pub struct EmbeddedAnf {
    //anf.modulus() == gen.min_poly()
    anf: Rc<ANFStructure>,
    gen: ComplexAlgebraic,
}

impl ANFStructure {
    pub fn all_complex_embeddings(&self) -> Vec<EmbeddedAnf> {
        self.modulus()
            .primitive_part_fof()
            .all_complex_roots()
            .into_iter()
            .map(|gen| EmbeddedAnf {
                anf: self.clone().into(),
                gen,
            })
            .collect()
    }
}

impl ComplexAlgebraic {
    pub fn abstract_generated_algebraic_number_field(&self) -> ANFStructure {
        new_anf(self.min_poly())
    }

    pub fn embedded_generated_algebraic_number_field(self) -> EmbeddedAnf {
        EmbeddedAnf {
            anf: self.abstract_generated_algebraic_number_field().into(),
            gen: self,
        }
    }
}

#[cfg(any())]
impl EmbeddedAnf {
    pub fn intersect_pair(field1: &Self, field2: &Self) -> Self {
        todo!()
    }

    pub fn intersect_list(fields: Vec<impl Borrow<Self>>) -> Self {
        todo!()
    }

    pub fn generated_pair(field1: &Self, field2: &Self) -> Self {
        todo!()
    }

    pub fn generated_list(fields: Vec<impl Borrow<Self>>) -> Self {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use malachite_nz::integer::Integer;

    use crate::{polynomial::polynomial::*, structure::*};

    #[test]
    fn test_embedded_anf_integral_basis() {
        let x = &Polynomial::<Integer>::var().into_ring();
        let f = (x.pow(2) + 7).into_set();
        for root in f.all_complex_roots() {
            println!(
                "{:?}",
                root.abstract_generated_algebraic_number_field()
                    .compute_integral_basis()
            );
        }
    }
}
