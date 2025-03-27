use algebraeon_nzq::{integer::Integer, natural::Natural};

use crate::number::natural::primes::{PrimeGenerator, is_prime};

use super::point::Point;
use std::collections::HashMap;

/// Returns one factor of n using Lenstra's 2 Stage Elliptic curve Factorization
/// with Suyama's Parameterization. Here Montgomery arithmetic is used for fast
/// computation of addition and doubling of points in elliptic curve.
///
/// This ECM method considers elliptic curves in Montgomery form (E : b*y^2*z = x^3 + a*x^2*z + x*z^2)
/// and involves elliptic curve operations (mod N), where the elements in Z are reduced (mod N).
/// Since N is not a prime, E over FF(N) is not really an elliptic curve but we can still do point additions
/// and doubling as if FF(N) was a field.
///
/// Stage 1: The basic algorithm involves taking a random point (P) on an elliptic curve in FF(N).
/// The compute k*P using Montgomery ladder algorithm.
/// Let q be an unknown factor of N. Then the order of the curve E, |E(FF(q))|,
/// might be a smooth number that divides k. Then we have k = l * |E(FF(q))|
/// for some l. For any point belonging to the curve E, |E(FF(q))|*P = O,
/// hence k*P = l*|E(FF(q))|*P. Thus kP.z_cord = 0 (mod q), and the unknown factor of N (q)
/// can be recovered by taking gcd(kP.z_cord, N).
///
/// Stage 2: This is a continuation of Stage 1 if k*P != O. The idea is to utilize
/// the fact that even if kP != 0, the value of k might miss just one large prime divisor
/// of |E(FF(q))|. In this case, we only need to compute the scalar multiplication by p
/// to get p*k*P = O. Here a second bound B2 restricts the size of possible values of p.
///
/// Parameters:
///
/// - `n`: Number to be factored.
/// - `B1`: Stage 1 Bound.
/// - `B2`: Stage 2 Bound.
/// - `max_curve`: Maximum number of curves generated.
/// - `rgen`: Random number generator.
pub fn ecm_one_factor(n: &Natural, b1: usize, b2: usize, max_curve: usize) -> Result<Integer, ()> {
    debug_assert_ne!(b1 % 2, 0);
    debug_assert_ne!(b2 % 2, 0);

    debug_assert!(!is_prime(n));

    let mut curve = 0;
    let d = (b2 as f64).sqrt() as usize;
    let two_d = 2 * d;
    let mut beta: Vec<Integer> = vec![Integer::default(); d + 1];
    let mut s: Vec<Point> = vec![Point::default(); d + 1];
    let mut k = Integer::from(1);
    let three = Integer::from(3);

    for p in PrimeGenerator::new().take_while(|&p| p <= b1) {
        k *= p.pow(b1.ilog(p));
    }

    while curve <= max_curve {
        curve += 1;

        // Suyama's Parametrization
        let sigma = (n - Integer::from(1)).random_below(rgen);
        let u = (&sigma * &sigma - Integer::from(5)) % n;
        let v: Integer = (4 * sigma) % n;
        let diff = Integer::from(&v - &u);
        let u_3 = u.clone().pow_mod(&three, n).unwrap();
        let v_3 = v.clone().pow_mod(&three, n).unwrap();

        let c = match (Integer::from(4) * &u_3 * &v).invert(n) {
            Ok(c) => {
                (diff.pow_mod(&three, n).unwrap() * (Integer::from(4) * &u + &v) * c
                    - Integer::from(2))
                    % n
            }
            _ => return Ok((Integer::from(4) * u_3 * v).gcd(n)),
        };

        let a24 = (c + 2) * Integer::from(4).invert(n).unwrap() % n;
        let q = Point::new(u_3, v_3, a24, n.clone());
        let q = q.mont_ladder(&k);
        let g = q.z_cord.clone().gcd(n);

        // Stage 1 factor
        if &g != n && g != 1 {
            return Ok(g);
        }

        // Stage 1 failure. Q.z = 0, Try another curve
        if &g == n {
            continue;
        }

        // Stage 2 - Improved Standard Continuation
        s[1] = q.double();
        s[2] = s[1].double();
        beta[1] = Integer::from(&s[1].x_cord * &s[1].z_cord) % n;
        beta[2] = Integer::from(&s[2].x_cord * &s[2].z_cord) % n;

        for d in 3..=(d) {
            s[d] = s[d - 1].add(&s[1], &s[d - 2]);
            beta[d] = Integer::from(&s[d].x_cord * &s[d].z_cord) % n;
        }

        let mut g = Integer::from(1);
        let b = b1 - 1;
        let mut t = q.mont_ladder(&Integer::from(b - two_d));
        let mut r = q.mont_ladder(&Integer::from(b));

        let mut primes = Primes::all().skip_while(|&q| q < b);
        for rr in (b..b2).step_by(two_d) {
            let alpha = Integer::from(&r.x_cord * &r.z_cord) % n;
            for q in primes.by_ref().take_while(|&q| q <= rr + two_d) {
                let delta = (q - rr) / 2;
                let f = Integer::from(&r.x_cord - &s[d].x_cord)
                    * Integer::from(&r.z_cord + &s[d].z_cord)
                    - &alpha
                    + &beta[delta];
                g = (g * f) % n;
            }
            // Swap
            std::mem::swap(&mut t, &mut r);
            r = r.add(&s[d], &t);
        }
        g = g.gcd(n);

        // Stage 2 Factor found
        if &g != n && g != 1 {
            return Ok(g);
        }
    }

    // ECM failed, Increase the bounds
    Err(Error::ECMFailed)
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
pub fn ecm(n: &Natural) -> Result<HashMap<Integer, usize>, Error> {
    let optimal_params = optimal_params(n.to_string().len());

    ecm_with_params(
        n,
        optimal_params.0,
        optimal_params.1,
        optimal_params.2,
        1234,
    )
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
/// - `B1`: Stage 1 Bound.
/// - `B2`: Stage 2 Bound.
/// - `max_curve`: Maximum number of curves generated.
/// - `seed`: Initialize pseudorandom generator.
pub fn ecm_with_params(
    n: &Integer,
    b1: usize,
    b2: usize,
    max_curve: usize,
    seed: usize,
) -> Result<HashMap<Integer, usize>, Error> {
    let mut factors = HashMap::new();

    let mut n: Integer = n.clone();
    for prime in Primes::all().take(100_000) {
        if n.is_divisible_u(prime as u32) {
            let prime = Integer::from(prime);
            while n.is_divisible(&prime) {
                n /= &prime;
                *factors.entry(prime.clone()).or_insert(0) += 1;
            }
        }
    }

    let mut rand_state = RandState::new();
    rand_state.seed(&seed.into());

    while n != 1 {
        let factor = ecm_one_factor(&n, b1, b2, max_curve, &mut rand_state).unwrap_or(n.clone());

        while n.is_divisible(&factor) {
            n /= &factor;
            *factors.entry(factor.clone()).or_insert(0) += 1;
        }
    }

    Ok(factors)
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use super::*;

    fn ecm(n: &Integer) -> Result<HashMap<Integer, usize>, Error> {
        super::ecm(n)
    }

    #[test]
    fn sympy_1() {
        assert_eq!(
            ecm(&Integer::from_str("398883434337287").unwrap()).unwrap(),
            HashMap::from([
                (Integer::from_str("99476569").unwrap(), 1),
                (Integer::from_str("4009823").unwrap(), 1),
            ])
        );
    }

    #[test]
    fn sympy_2() {
        assert_eq!(
            ecm(&Integer::from_str("46167045131415113").unwrap()).unwrap(),
            HashMap::from([
                (Integer::from_str("43").unwrap(), 1),
                (Integer::from_str("2634823").unwrap(), 1),
                (Integer::from_str("407485517").unwrap(), 1),
            ])
        );
    }

    #[test]
    fn sympy_3() {
        assert_eq!(
            ecm(&Integer::from_str("64211816600515193").unwrap()).unwrap(),
            HashMap::from([
                (Integer::from_str("281719").unwrap(), 1),
                (Integer::from_str("359641").unwrap(), 1),
                (Integer::from_str("633767").unwrap(), 1),
            ])
        );
    }

    #[test]
    fn sympy_4() {
        assert_eq!(
            ecm(&Integer::from_str("168541512131094651323").unwrap()).unwrap(),
            HashMap::from([
                (Integer::from_str("79").unwrap(), 1),
                (Integer::from_str("113").unwrap(), 1),
                (Integer::from_str("11011069").unwrap(), 1),
                (Integer::from_str("1714635721").unwrap(), 1),
            ])
        );
    }

    #[test]
    fn sympy_5() {
        assert_eq!(
            ecm(&Integer::from_str("631211032315670776841").unwrap()).unwrap(),
            HashMap::from([
                (Integer::from_str("9312934919").unwrap(), 1),
                (Integer::from_str("67777885039").unwrap(), 1),
            ])
        );
    }

    #[test]
    fn sympy_6() {
        assert_eq!(
            ecm(&Integer::from_str("4132846513818654136451").unwrap()).unwrap(),
            HashMap::from([
                (Integer::from_str("47").unwrap(), 1),
                (Integer::from_str("160343").unwrap(), 1),
                (Integer::from_str("2802377").unwrap(), 1),
                (Integer::from_str("195692803").unwrap(), 1),
            ])
        );
    }

    #[test]
    fn sympy_7() {
        assert_eq!(
            ecm(&Integer::from_str("4516511326451341281684513").unwrap()).unwrap(),
            HashMap::from([
                (Integer::from_str("3").unwrap(), 2),
                (Integer::from_str("39869").unwrap(), 1),
                (Integer::from_str("131743543").unwrap(), 1),
                (Integer::from_str("95542348571").unwrap(), 1),
            ])
        );
    }

    #[test]
    fn sympy_8() {
        assert_eq!(
            ecm(&Integer::from_str("3146531246531241245132451321").unwrap(),).unwrap(),
            HashMap::from([
                (Integer::from_str("3").unwrap(), 1),
                (Integer::from_str("100327907731").unwrap(), 1),
                (Integer::from_str("10454157497791297").unwrap(), 1),
            ])
        );
    }

    #[test]
    fn sympy_9() {
        assert_eq!(
            ecm(&Integer::from_str("4269021180054189416198169786894227").unwrap()).unwrap(),
            HashMap::from([
                (Integer::from_str("184039").unwrap(), 1),
                (Integer::from_str("241603").unwrap(), 1),
                (Integer::from_str("333331").unwrap(), 1),
                (Integer::from_str("477973").unwrap(), 1),
                (Integer::from_str("618619").unwrap(), 1),
                (Integer::from_str("974123").unwrap(), 1),
            ])
        );
    }

    #[test]
    fn same_factors() {
        assert_eq!(
            ecm(&Integer::from_str("7853316850129").unwrap()).unwrap(),
            HashMap::from([(Integer::from_str("2802377").unwrap(), 2)])
        );
    }

    #[test]
    fn small_prime() {
        assert_eq!(
            ecm(&Integer::from(17)).unwrap(),
            HashMap::from([(Integer::from(17), 1)])
        );
    }

    #[test]
    fn big_prime() {
        assert_eq!(
            ecm(&Integer::from_str("21472883178031195225853317139").unwrap()).unwrap(),
            HashMap::from([(
                Integer::from_str("21472883178031195225853317139").unwrap(),
                1
            )])
        );
    }
}
