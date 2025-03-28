use super::natural::*;

pub struct Rng {
    rng: malachite_base::num::random::RandomPrimitiveInts<u64>,
}

impl Rng {
    pub fn new() -> Self {
        Self {
            rng: malachite_base::num::random::random_primitive_ints(
                malachite_base::random::EXAMPLE_SEED,
            ),
        }
    }
}

impl Natural {
    pub fn random_below(&self, rng: &mut Rng) -> Natural {
        Natural::from_malachite(malachite_nz::natural::random::get_random_natural_less_than(
            &mut rng.rng,
            (self + Natural::ONE).to_malachite_ref(),
        ))
    }
}
