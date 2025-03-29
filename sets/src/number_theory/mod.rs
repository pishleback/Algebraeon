/// Return an iterator generating all prime numbers.
pub fn primes() -> impl Iterator<Item = usize> {
    use malachite_base::num::factorization::traits::Primes;
    usize::primes()
}
