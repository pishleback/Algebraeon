use std::{borrow::Borrow, rc::Rc};

use super::{isolated_roots::*, number_field::*};

#[derive(Debug, Clone)]
struct EmbeddedAnf {
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
    use std::str::FromStr;

    use malachite_nz::integer::Integer;
    use malachite_q::Rational;

    use crate::rings::{polynomial::polynomial::*, structure::*};

    use super::*;

    #[test]
    fn test_embedded_anf_run() {
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
