use orthoclase_groups::examples::symmetric::Permutation;
use orthoclase_rings::number::modulo::Modulo;

pub fn modular_permutation<const N: usize>(
    f: impl Fn(Modulo<N>) -> Modulo<N>,
) -> Result<Permutation<N>, &'static str> {
    let mut perm = [0; N];
    for i in 0..N {
        perm[i] = f(i.into()).into()
    }
    Permutation::new(perm)
}
