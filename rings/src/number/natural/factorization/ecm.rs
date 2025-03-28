/*
This module is derived from sympy https://www.sympy.org/en/index.html and is included under the following lisence

Copyright (c) 2006-2023 SymPy Development Team

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

  a. Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
  b. Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
  c. Neither the name of SymPy nor the names of its contributors
     may be used to endorse or promote products derived from this software
     without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

use crate::{
    number::natural::{factorization::primes::is_prime, functions::gcd},
    structure::MetaSemiRing,
};
use algebraeon_nzq::{
    natural::{Natural, primes},
    random::Rng,
};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::collections::HashSet;

/// Montgomery form of Points in an elliptic curve.
///
/// In this form, the addition and doubling of points
/// does not need any y-coordinate information thus
/// decreasing the number of operations.
/// Using Montgomery form we try to perform point addition
/// and doubling in least amount of multiplications.
///
/// The elliptic curve used here is of the form
/// `(E : b*y**2*z = x**3 + a*x**2*z + x*z**2)`.
/// The `a_24` parameter is equal to `(a + 2)/4`.
///
/// References
/// ----------
/// - http://www.hyperelliptic.org/tanja/SHARCS/talks06/Gaj.pdf
#[derive(Default, Debug, Clone)]
pub struct Point {
    /// X coordinate of the Point
    pub x_cord: Natural,
    /// Z coordinate of the Point
    pub z_cord: Natural,
    /// Parameter of the elliptic curve in Montgomery form
    pub a_24: Natural,
    /// modulus
    pub modulus: Natural,
}

impl Point {
    /// Initial parameters for the Point struct.
    ///
    /// # Parameters
    ///
    /// - `x_cord`: X coordinate of the Point
    /// - `z_cord`: Z coordinate of the Point
    /// - `a_24`: Parameter of the elliptic curve in Montgomery form
    /// - `mod`: modulus
    pub fn new(x_cord: Natural, z_cord: Natural, a_24: Natural, modulus: Natural) -> Point {
        Point {
            x_cord,
            z_cord,
            a_24,
            modulus,
        }
    }

    /// Adds two points `self` and `Q` where `diff = self - Q`.
    ///
    /// This algorithm requires 6 multiplications. The assumption is that `self.x_cord * Q.x_cord * (self.x_cord - Q.x_cord) != 0`.
    /// Using this algorithm speeds up the addition by reducing the number of multiplications required.
    ///
    /// The `mont_ladder` algorithm is constructed in a way that the difference between intermediate points is always equal to the initial point.
    /// So, we always know what the difference between the point is.
    ///
    /// # Parameters
    ///
    /// - `Q`: Point on the curve in Montgomery form.
    /// - `diff`: `self - Q`
    pub fn add(&self, q: &Point, diff: &Point) -> Point {
        let u = (&self.modulus + &self.x_cord - &self.z_cord) * (&q.x_cord + &q.z_cord);
        let v = (&self.x_cord + &self.z_cord) * (&self.modulus + &q.x_cord - &q.z_cord);
        let add = &u + &v;
        let subt = u + (-v) % &self.modulus;
        let x_cord = (&diff.z_cord * &add * &add) % &self.modulus;
        let z_cord = (&diff.x_cord * &subt * &subt) % &self.modulus;

        Point::new(x_cord, z_cord, self.a_24.clone(), self.modulus.clone())
    }

    /// Doubles a point in an elliptic curve in Montgomery form.
    pub fn double(&self) -> Point {
        let u = (&self.x_cord + &self.z_cord) * (&self.x_cord + &self.z_cord);
        let v = (&self.modulus + &self.x_cord - &self.z_cord)
            * (&self.modulus + &self.x_cord - &self.z_cord);
        let diff = &u + (-&v) % &self.modulus;
        let x_cord = (u * &v) % &self.modulus;
        let z_cord = ((v + &self.a_24 * &diff) * diff) % &self.modulus;

        Point::new(x_cord, z_cord, self.a_24.clone(), self.modulus.clone())
    }

    /// Scalar multiplication of a point in Montgomery form
    /// using Montgomery Ladder Algorithm.
    /// A total of 11 multiplications are required in each step of this
    /// algorithm.
    ///
    /// # Parameters
    ///
    /// - `k`: The positive integer multiplier
    pub fn mont_ladder(&self, k: &Natural) -> Point {
        let mut q = self.clone();
        let mut r = self.double();
        for i in k.bits().rev().skip(1) {
            if i {
                q = r.add(&q, self);
                r = r.double();
            } else {
                r = q.add(&r, self);
                q = q.double();
            }
        }
        q
    }
}

impl PartialEq for Point {
    /// Two points are equal if X/Z of both points are equal.
    fn eq(&self, other: &Self) -> bool {
        if self.a_24 != other.a_24 {
            return false;
        }
        let modulus = &self.modulus;
        if modulus != &other.modulus {
            return false;
        }

        (self.z_cord.mod_inv_ref(modulus).unwrap() * &self.x_cord) % modulus
            == (other.z_cord.mod_inv_ref(modulus).unwrap() * &other.x_cord) % modulus
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_point_add() {
        let p1 = Point::new(11u32.into(), 16u32.into(), 7u32.into(), 29u32.into());
        let p2 = Point::new(13u32.into(), 10u32.into(), 7u32.into(), 29u32.into());
        let p3 = p2.add(&p1, &p1);

        assert_eq!(p3.x_cord, Natural::from(23u32));
        assert_eq!(p3.z_cord, Natural::from(17u32));
    }

    #[test]
    fn test_point_double() {
        let p1 = Point::new(11u32.into(), 16u32.into(), 7u32.into(), 29u32.into());
        let p2 = p1.double();

        assert_eq!(p2.x_cord, Natural::from(13u32));
        assert_eq!(p2.z_cord, Natural::from(10u32));
    }

    #[test]
    fn test_point_mont_ladder() {
        let p1 = Point::new(11u32.into(), 16u32.into(), 7u32.into(), 29u32.into());
        let p3 = p1.mont_ladder(&3u32.into());

        assert_eq!(p3.x_cord, Natural::from(23u32));
        assert_eq!(p3.z_cord, Natural::from(17u32));
    }

    #[test]
    fn test_point() {
        let a: Natural = 10u32.into();
        let a_24: Natural =
            (a + Natural::from(2u32)) * Natural::from(4u32).mod_inv(&101u32.into()).unwrap();

        let p1 = Point::new(10u32.into(), 17u32.into(), a_24.clone(), 101u32.into());
        let p2 = p1.double();
        assert_eq!(
            p2,
            Point::new(68u32.into(), 56u32.into(), a_24.clone(), 101u32.into())
        );
        let p4 = p2.double();
        assert_eq!(
            p4,
            Point::new(22u32.into(), 64u32.into(), a_24.clone(), 101u32.into())
        );
        let p8 = p4.double();
        assert_eq!(
            p8,
            Point::new(71u32.into(), 95u32.into(), a_24.clone(), 101u32.into())
        );
        let p16 = p8.double();
        assert_eq!(
            p16,
            Point::new(5u32.into(), 16u32.into(), a_24.clone(), 101u32.into())
        );
        let p32 = p16.double();
        assert_eq!(
            p32,
            Point::new(33u32.into(), 96u32.into(), a_24.clone(), 101u32.into())
        );

        // p3 = p2 + p1
        let p3 = p2.add(&p1, &p1);
        assert_eq!(
            p3,
            Point::new(1u32.into(), 61u32.into(), a_24.clone(), 101u32.into())
        );
        // p5 = p3 + p2 or p4 + p1
        let p5 = p3.add(&p2, &p1);
        assert_eq!(
            p5,
            Point::new(49u32.into(), 90u32.into(), a_24.clone(), 101u32.into())
        );
        assert_eq!(p5, p4.add(&p1, &p3));
        // # p6 = 2*p3
        let p6 = p3.double();
        assert_eq!(
            p6,
            Point::new(87u32.into(), 43u32.into(), a_24.clone(), 101u32.into())
        );
        assert_eq!(p6, p4.add(&p2, &p2));
        // # p7 = p5 + p2
        let p7 = p5.add(&p2, &p3);
        assert_eq!(
            p7,
            Point::new(69u32.into(), 23u32.into(), a_24.clone(), 101u32.into())
        );
        assert_eq!(p7, p4.add(&p3, &p1));
        assert_eq!(p7, p6.add(&p1, &p5));
        // # p9 = p5 + p4
        let p9 = p5.add(&p4, &p1);
        assert_eq!(
            p9,
            Point::new(56u32.into(), 99u32.into(), a_24, 101u32.into())
        );
        assert_eq!(p9, p6.add(&p3, &p3));
        assert_eq!(p9, p7.add(&p2, &p5));
        assert_eq!(p9, p8.add(&p1, &p7));

        assert_eq!(p5, p1.mont_ladder(&5u32.into()));
        assert_eq!(p9, p1.mont_ladder(&9u32.into()));
        assert_eq!(p16, p1.mont_ladder(&16u32.into()));
        assert_eq!(p9, p3.mont_ladder(&3u32.into()));
    }
}

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
pub fn ecm_one_factor_raw(
    n: &Natural,
    b1: usize,
    b2: usize,
    max_curve: usize,
    rng: &mut Rng,
) -> Result<Natural, ()> {
    debug_assert_eq!(b1 % 2, 0);
    debug_assert_eq!(b2 % 2, 0);

    debug_assert!(!is_prime(n));

    // When calculating T, if (B1 - 2*D) is negative, it cannot be calculated.
    let big_d = std::cmp::min(b2.isqrt(), b1 / 2 - 1);
    let mut k = Natural::ONE;
    for p in primes().take_while(|&p| p <= b1) {
        k *= Natural::from(p).nat_pow(&b1.ilog(p).into());
    }
    // Pre-calculate the prime numbers to be used in stage 2.
    // Using the fact that the x-coordinates of point P and its
    // inverse -P coincide, the number of primes to be checked
    // in stage 2 can be reduced.

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

    debug_assert!(n >= &Natural::from(7u32));

    // Search over the curves in parallel
    let result = (0..max_curve)
        .map(|_| (n - Natural::from(7u32)).random_below(rng) + Natural::from(6u32))
        .collect::<Vec<_>>()
        .into_par_iter()
        .find_map_any(|sigma| {
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
                        return None;
                    } else {
                        return Some(g);
                    }
                }
            };
            let v3 = (&v * &v * &v) % n;

            let q = Point::new(u3, v3, a24, n.clone());
            let q = q.mont_ladder(&k);
            let g = gcd(q.z_cord.clone(), n.clone());

            if g != Natural::ONE && g != *n {
                // Stage 1 factor found
                return Some(g);
            } else if g == *n {
                // Stage 1 failure. Q.z = 0, Try another curve
                return None;
            }

            // Stage 2 - Improved Standard Continuation
            let mut beta = vec![Natural::default(); big_d];
            let mut s = vec![Point::default(); big_d];

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
                    let f = (&r.x_cord + (-&s[*delta].x_cord) % n)
                        * (&r.z_cord + &s[*delta].z_cord)
                        + (-&alpha) % n
                        + &beta[*delta];
                    g = (g * f) % n;
                }
                (t, r) = (r.clone(), r.add(&w, &t));
            }
            g = gcd(g, n.clone());
            if g != Natural::ONE && g != *n {
                // Stage 2 Factor found
                return Some(g);
            }

            None
        });

    match result {
        Some(d) => Ok(d),
        None => Err(()),
    }
}

// TODO: don't use this once I have qnfs
// Instead try ecm_one_factor_target_digits a on 20 digits to get small factors then qnfs
pub fn ecm_one_factor_target_digits(
    n: &Natural,
    fith_target_factor_digits: usize,
    rng: &mut Rng,
) -> Result<Natural, ()> {
    // The target factor has fith_target_factor_digits*5 many digits
    let (b1, b2, max_curve) = match fith_target_factor_digits {
        0 | 1 | 2 => (2_000, 160_000, 35),
        3 => (5_000, 500_000, 50),
        4 => (11_000, 1_900_000, 74),
        5 => (50_000, 13_000_000, 214),
        6 => (250_000, 130_000_000, 430),
        7 => (1_000_000, 1_000_000_000, 904),
        8 => (3_000_000, 5_700_000_000, 2350),
        9 => (11_000_000, 35_000_000_000, 4480),
        10 => (44_000_000, 240_000_000_000, 7553),
        11 => (110_000_000, 780_000_000_000, 17769),
        12 => (260_000_000, 3_200_000_000_000, 42017),
        _ => (850_000_000, 16_000_000_000_000, 69408),
    };
    ecm_one_factor_raw(n, b1, b2, max_curve, rng)
}
