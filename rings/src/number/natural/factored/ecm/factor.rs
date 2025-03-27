use std::collections::HashSet;

use algebraeon_nzq::{
    natural::{Natural, primes},
    random::Rng,
};

use crate::{
    number::natural::{functions::gcd, primes::is_prime},
    structure::structure::MetaSemiRing,
};

use super::point::Point;

/*
Returns one factor of n using
Lenstra's 2 Stage Elliptic curve Factorization
with Suyama's Parameterization. Here Montgomery
arithmetic is used for fast computation of addition
and doubling of points in elliptic curve.

Explanation
===========

This ECM method considers elliptic curves in Montgomery
form (E : b*y**2*z = x**3 + a*x**2*z + x*z**2) and involves
elliptic curve operations (mod N), where the elements in
Z are reduced (mod N). Since N is not a prime, E over FF(N)
is not really an elliptic curve but we can still do point additions
and doubling as if FF(N) was a field.

Stage 1 : The basic algorithm involves taking a random point (P) on an
elliptic curve in FF(N). The compute k*P using Montgomery ladder algorithm.
Let q be an unknown factor of N. Then the order of the curve E, |E(FF(q))|,
might be a smooth number that divides k. Then we have k = l * |E(FF(q))|
for some l. For any point belonging to the curve E, |E(FF(q))|*P = O,
hence k*P = l*|E(FF(q))|*P. Thus kP.z_cord = 0 (mod q), and the unknownn
factor of N (q) can be recovered by taking gcd(kP.z_cord, N).

Stage 2 : This is a continuation of Stage 1 if k*P != O. The idea utilize
the fact that even if kP != 0, the value of k might miss just one large
prime divisor of |E(FF(q))|. In this case we only need to compute the
scalar multiplication by p to get p*k*P = O. Here a second bound B2
restrict the size of possible values of p.

Parameters
==========

n : Number to be Factored
B1 : Stage 1 Bound. Must be an even number.
B2 : Stage 2 Bound. Must be an even number.
max_curve : Maximum number of curves generated

Returns
=======

integer | None : ``n`` (if it is prime) else a non-trivial divisor of ``n``. ``None`` if not found

References
==========

.. [1] Carl Pomerance, Richard Crandall, Prime Numbers: A Computational Perspective,
        2nd Edition (2005), page 344, ISBN:978-0387252827
*/
pub fn ecm_one_factor(
    n: &Natural,
    b1: usize,
    b2: usize,
    max_curve: usize,
    rng: &mut Rng,
) -> Result<Natural, ()> {
    debug_assert_eq!(b1 % 2, 0);
    debug_assert_eq!(b2 % 2, 0);

    debug_assert!(!is_prime(n));

    println!("n = {} b1 = {} b2 = {}", n, b1, b2);

    let b1: usize = 10000;
    let b2: usize = 100000;

    // When calculating T, if (B1 - 2*D) is negative, it cannot be calculated.
    let big_d = std::cmp::min(b2.isqrt(), b1 / 2 - 1);
    println!("D = {}", big_d);
    let mut beta = vec![Natural::default(); big_d];
    let mut s = vec![Point::default(); big_d];
    let mut k = Natural::ONE;
    for p in primes().take_while(|&p| p <= b1) {
        k *= Natural::from(p).nat_pow(&b1.ilog(p).into());
    }
    // Pre-calculate the prime numbers to be used in stage 2.
    // Using the fact that the x-coordinates of point P and its
    // inverse -P coincide, the number of primes to be checked
    // in stage 2 can be reduced.

    println!("got k");

    let mut deltas_list = vec![];
    for r in (b1 + 2 * big_d..b2 + 2 * big_d).step_by(4 * big_d) {
        let mut deltas = HashSet::new();
        for q in primes()
            .take_while(|&q| q < r + 2 * big_d)
            .filter(|&q| r - 2 * big_d < q)
        {
            deltas.insert((q.abs_diff(r) - 1) / 2);
        }
        // d in deltas iff r+(2d+1) and/or r-(2d+1) is prime
        deltas_list.push(deltas.into_iter().collect::<Vec<_>>())
    }

    println!("got deltas");

    debug_assert!(n >= &Natural::from(7u32));

    for c in 0..max_curve {
        println!("c = {}", c);

        // Suyama's Parametrization
        // random sigma in the range [6, n-1]
        let sigma = (n - Natural::from(7u32)).random_below(rng) + Natural::from(6u32);
        let u = (&sigma * &sigma - Natural::from(5u32)) % n;
        let v = (&sigma * Natural::from(4u32)) % n;
        let u3 = (&u * &u * &u) % n;
        let diff = &v + (-&u) % n;
        // We use the elliptic curve y**2 = x**3 + a*x**2 + x
        // where a = pow(v - u, 3, n)*(3*u + v)*invert(4*u_3*v, n) - 2
        // However, we do not declare a because it is more convenient
        // to use a24 = (a + 2)*invert(4, n) in the calculation.
        let u3_16_v = Natural::from(16u32) * &u3 * &v;
        let a24 = match u3_16_v.mod_inv(n) {
            Ok(u3_16_v_inv) => {
                (&diff * &diff * &diff * (Natural::from(3u32) * &u + &v) * u3_16_v_inv) % n
            }
            Err(()) => {
                let g = gcd(Natural::from(2u32) * u3 * v, n.clone());
                debug_assert_ne!(g, Natural::ONE);
                if &g == n {
                    continue;
                } else {
                    return Ok(g);
                }
            }
        };
        let v3 = (&v * &v * &v) % n;

        let q = Point::new(u3, v3, a24, n.clone());
        let q = q.mont_ladder(&k);
        let g = gcd(q.z_cord.clone(), n.clone());

        if g != Natural::ONE && g != *n {
            // Stage 1 factor found
            return Ok(g);
        } else if g == *n {
            // Stage 1 failure. Q.z = 0, Try another curve
            continue;
        }

        // Stage 2 - Improved Standard Continuation
        s[0] = q.clone();
        let q2 = q.double();
        s[1] = q2.add(&q, &q);
        beta[0] = (&s[0].x_cord * &s[0].z_cord) % n;
        beta[1] = (&s[1].x_cord * &s[1].z_cord) % n;
        for d in 2..big_d {
            s[d] = s[d - 1].add(&q2, &s[d - 2]);
            beta[d] = (&s[d].x_cord * &s[d].z_cord) % n;
        }
        // i.e., S[i] = Q.mont_ladder(2*i + 1)

        let mut g = Natural::ONE;
        let w = q.mont_ladder(&(4 * big_d).into());
        let mut t = q.mont_ladder(&(b1 - 2 * big_d).into());
        let mut r = q.mont_ladder(&(b1 + 2 * big_d).into());
        for deltas in &deltas_list {
            // R = Q.mont_ladder(r) where r in range(B1 + 2*D, B2 + 2*D, 4*D)
            let alpha = (&r.x_cord * &r.z_cord) % n;
            for delta in deltas {
                // We want to calculate
                // f = R.x_cord * S[delta].z_cord - S[delta].x_cord * R.z_cord
                let f = (&r.x_cord + (-&s[*delta].x_cord) % n) * (&r.z_cord + &s[*delta].z_cord)
                    + (-&alpha) % n
                    + &beta[*delta];
                g = (g * f) % n;
            }
            (t, r) = (r.clone(), r.add(&w, &t));
        }
        g = gcd(g, n.clone());
        if g != Natural::ONE && g != *n {
            // Stage 2 Factor found
            return Ok(g);
        }
    }

    Err(())
}

/// Optimal params retrieved from <https://gitlab.inria.fr/zimmerma/ecm>
fn optimal_params(digits: usize) -> (usize, usize, usize) {
    match digits {
        1..=10 => (2_000, 160_000, 35),
        11..=15 => (5_000, 500_000, 500),
        16..=20 => (11_000, 1_900_000, 74),
        21..=25 => (50_000, 13_000_000, 214),
        26..=30 => (250_000, 130_000_000, 430),
        31..=35 => (1_000_000, 1_000_000_000, 904),
        36..=40 => (3_000_000, 5_700_000_000, 2350),
        41..=45 => (11_000_000, 35_000_000_000, 4480),
        46..=50 => (44_000_000, 240_000_000_000, 7553),
        51..=55 => (110_000_000, 780_000_000_000, 17769),
        56..=60 => (260_000_000, 3_200_000_000_000, 42017),
        _ => (850_000_000, 16_000_000_000_000, 69408),
    }
}

/// Performs factorization using Lenstra's Elliptic curve method.
///
/// This function repeatedly calls `ecm_one_factor` to compute the factors
/// of n. First all the small factors are taken out using trial division.
/// Then `ecm_one_factor` is used to compute one factor at a time.
///
/// # Parameters
///
/// - `n`: Number to be factored.
pub fn ecm(n: &Natural) -> Natural {
    let mut optimal_params = optimal_params(n.to_string().len());
    let mut rgen = Rng::new();
    loop {
        match ecm_one_factor(
            n,
            optimal_params.0,
            optimal_params.1,
            optimal_params.2,
            &mut rgen,
        ) {
            Ok(d) => {
                return d;
            }
            Err(()) => {
                optimal_params.0 *= 2;
                optimal_params.1 *= 2;
                continue;
            }
        }
    }
}

// /// Performs factorization using Lenstra's Elliptic curve method.
// ///
// /// This function repeatedly calls `ecm_one_factor` to compute the factors
// /// of n. First all the small factors are taken out using trial division.
// /// Then `ecm_one_factor` is used to compute one factor at a time.
// ///
// /// # Parameters
// ///
// /// - `n`: Number to be factored.
// /// - `B1`: Stage 1 Bound.
// /// - `B2`: Stage 2 Bound.
// /// - `max_curve`: Maximum number of curves generated.
// /// - `seed`: Initialize pseudorandom generator.
// pub fn ecm_with_params(
//     n: &Integer,
//     b1: usize,
//     b2: usize,
//     max_curve: usize,
//     seed: usize,
// ) -> Result<HashMap<Integer, usize>, Error> {
//     let mut factors = HashMap::new();

//     let mut n: Integer = n.clone();
//     for prime in Primes::all().take(100_000) {
//         if n.is_divisible_u(prime as u32) {
//             let prime = Integer::from(prime);
//             while n.is_divisible(&prime) {
//                 n /= &prime;
//                 *factors.entry(prime.clone()).or_insert(0) += 1;
//             }
//         }
//     }

//     let mut rand_state = RandState::new();
//     rand_state.seed(&seed.into());

//     while n != 1 {
//         let factor = ecm_one_factor(&n, b1, b2, max_curve, &mut rand_state).unwrap_or(n.clone());

//         while n.is_divisible(&factor) {
//             n /= &factor;
//             *factors.entry(factor.clone()).or_insert(0) += 1;
//         }
//     }

//     Ok(factors)
// }

// #[cfg(test)]
// mod tests {
//     use std::str::FromStr;

//     use super::*;

//     fn ecm(n: &Integer) -> Result<HashMap<Integer, usize>, Error> {
//         super::ecm(n)
//     }

//     #[test]
//     fn sympy_1() {
//         assert_eq!(
//             ecm(&Integer::from_str("398883434337287").unwrap()).unwrap(),
//             HashMap::from([
//                 (Integer::from_str("99476569").unwrap(), 1),
//                 (Integer::from_str("4009823").unwrap(), 1),
//             ])
//         );
//     }

//     #[test]
//     fn sympy_2() {
//         assert_eq!(
//             ecm(&Integer::from_str("46167045131415113").unwrap()).unwrap(),
//             HashMap::from([
//                 (Integer::from_str("43").unwrap(), 1),
//                 (Integer::from_str("2634823").unwrap(), 1),
//                 (Integer::from_str("407485517").unwrap(), 1),
//             ])
//         );
//     }

//     #[test]
//     fn sympy_3() {
//         assert_eq!(
//             ecm(&Integer::from_str("64211816600515193").unwrap()).unwrap(),
//             HashMap::from([
//                 (Integer::from_str("281719").unwrap(), 1),
//                 (Integer::from_str("359641").unwrap(), 1),
//                 (Integer::from_str("633767").unwrap(), 1),
//             ])
//         );
//     }

//     #[test]
//     fn sympy_4() {
//         assert_eq!(
//             ecm(&Integer::from_str("168541512131094651323").unwrap()).unwrap(),
//             HashMap::from([
//                 (Integer::from_str("79").unwrap(), 1),
//                 (Integer::from_str("113").unwrap(), 1),
//                 (Integer::from_str("11011069").unwrap(), 1),
//                 (Integer::from_str("1714635721").unwrap(), 1),
//             ])
//         );
//     }

//     #[test]
//     fn sympy_5() {
//         assert_eq!(
//             ecm(&Integer::from_str("631211032315670776841").unwrap()).unwrap(),
//             HashMap::from([
//                 (Integer::from_str("9312934919").unwrap(), 1),
//                 (Integer::from_str("67777885039").unwrap(), 1),
//             ])
//         );
//     }

//     #[test]
//     fn sympy_6() {
//         assert_eq!(
//             ecm(&Integer::from_str("4132846513818654136451").unwrap()).unwrap(),
//             HashMap::from([
//                 (Integer::from_str("47").unwrap(), 1),
//                 (Integer::from_str("160343").unwrap(), 1),
//                 (Integer::from_str("2802377").unwrap(), 1),
//                 (Integer::from_str("195692803").unwrap(), 1),
//             ])
//         );
//     }

//     #[test]
//     fn sympy_7() {
//         assert_eq!(
//             ecm(&Integer::from_str("4516511326451341281684513").unwrap()).unwrap(),
//             HashMap::from([
//                 (Integer::from_str("3").unwrap(), 2),
//                 (Integer::from_str("39869").unwrap(), 1),
//                 (Integer::from_str("131743543").unwrap(), 1),
//                 (Integer::from_str("95542348571").unwrap(), 1),
//             ])
//         );
//     }

//     #[test]
//     fn sympy_8() {
//         assert_eq!(
//             ecm(&Integer::from_str("3146531246531241245132451321").unwrap(),).unwrap(),
//             HashMap::from([
//                 (Integer::from_str("3").unwrap(), 1),
//                 (Integer::from_str("100327907731").unwrap(), 1),
//                 (Integer::from_str("10454157497791297").unwrap(), 1),
//             ])
//         );
//     }

//     #[test]
//     fn sympy_9() {
//         assert_eq!(
//             ecm(&Integer::from_str("4269021180054189416198169786894227").unwrap()).unwrap(),
//             HashMap::from([
//                 (Integer::from_str("184039").unwrap(), 1),
//                 (Integer::from_str("241603").unwrap(), 1),
//                 (Integer::from_str("333331").unwrap(), 1),
//                 (Integer::from_str("477973").unwrap(), 1),
//                 (Integer::from_str("618619").unwrap(), 1),
//                 (Integer::from_str("974123").unwrap(), 1),
//             ])
//         );
//     }

//     #[test]
//     fn same_factors() {
//         assert_eq!(
//             ecm(&Integer::from_str("7853316850129").unwrap()).unwrap(),
//             HashMap::from([(Integer::from_str("2802377").unwrap(), 2)])
//         );
//     }

//     #[test]
//     fn small_prime() {
//         assert_eq!(
//             ecm(&Integer::from(17)).unwrap(),
//             HashMap::from([(Integer::from(17), 1)])
//         );
//     }

//     #[test]
//     fn big_prime() {
//         assert_eq!(
//             ecm(&Integer::from_str("21472883178031195225853317139").unwrap()).unwrap(),
//             HashMap::from([(
//                 Integer::from_str("21472883178031195225853317139").unwrap(),
//                 1
//             )])
//         );
//     }
// }
