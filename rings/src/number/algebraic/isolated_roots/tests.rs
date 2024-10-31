use std::str::FromStr;

use malachite_base::num::basic::traits::Two;
use malachite_nz::{integer::Integer, natural::Natural};
use malachite_q::Rational;

use crate::{number::algebraic::isolated_roots::{anf_multi_primitive_element_theorem, anf_pair_primitive_element_theorem, as_poly_expr, ComplexAlgebraic, RealAlgebraic}, polynomial::polynomial::Polynomial, ring_structure::cannonical::{PositiveRealNthRootOpps, Ring, UniqueFactorizationDomain}};

use super::{IntegralDomain, StructuredType};

#[test]
fn test_squarefree_polynomial_real_root_isolation() {
    let f = Polynomial::product(vec![
        &Polynomial::from_coeffs(vec![
            Integer::from(-2),
            Integer::from(-4),
            Integer::from(-2),
        ]),
        &Polynomial::from_coeffs(vec![Integer::from(6), Integer::from(0), Integer::from(-3)]),
        &Polynomial::from_coeffs(vec![Integer::from(1), Integer::from(-3), Integer::from(1)]),
        &Polynomial::from_coeffs(vec![
            Integer::from(2),
            Integer::from(-3),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(1),
        ]),
        &Polynomial::from_coeffs(vec![
            Integer::from(1),
            Integer::from(-3),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(1),
        ]),
        &Polynomial::from_coeffs(vec![
            Integer::from(-1),
            Integer::from(12),
            Integer::from(-4),
            Integer::from(-15),
            Integer::from(5),
            Integer::from(3),
            Integer::from(-1),
        ]),
    ]);
    let f = f.primitive_squarefree_part();
    //f is a squarefree polynomial with lots of roots
    println!("f = {:?}", f);
    let intervals = Polynomial::real_roots_squarefree(f, None, None, false, false);
    println!("intervals = {:?}", &intervals);
    intervals.check_invariants().unwrap();

    let f = Polynomial::from_coeffs(vec![Integer::from(1), Integer::from(-3), Integer::from(1)]);
    println!("f = {:?}", f);
    let mut intervals = Polynomial::real_roots_squarefree(f, None, None, false, false);
    intervals.check_invariants().unwrap();
    intervals.clone().to_real_roots();
    for root in intervals.clone().to_real_roots() {
        println!("root = {:?}", root);
        root.check_invariants().unwrap();
    }
    println!("refine");
    for _i in 0..10 {
        intervals.refine_all();
    }
    for root in intervals.clone().to_real_roots() {
        println!("root = {:?}", root);
        root.check_invariants().unwrap();
    }

    let f = Polynomial::product(vec![
        &Polynomial::from_coeffs(vec![Integer::from(-1), Integer::from(1)]),
        &Polynomial::from_coeffs(vec![Integer::from(-2), Integer::from(1)]),
        &Polynomial::from_coeffs(vec![Integer::from(-3), Integer::from(1)]),
        &Polynomial::from_coeffs(vec![Integer::from(-4), Integer::from(1)]),
    ]);
    assert_eq!(
        f.clone()
            .real_roots_squarefree(
                Some(&Rational::from(1)),
                Some(&Rational::from(4)),
                false,
                false
            )
            .intervals
            .len(),
        2
    );
    assert_eq!(
        f.clone()
            .real_roots_squarefree(
                Some(&Rational::from(1)),
                Some(&Rational::from(4)),
                true,
                false
            )
            .intervals
            .len(),
        3
    );
    assert_eq!(
        f.clone()
            .real_roots_squarefree(
                Some(&Rational::from(1)),
                Some(&Rational::from(4)),
                false,
                true
            )
            .intervals
            .len(),
        3
    );
    assert_eq!(
        f.clone()
            .real_roots_squarefree(
                Some(&Rational::from(1)),
                Some(&Rational::from(4)),
                true,
                true
            )
            .intervals
            .len(),
        4
    );
}

#[test]
fn test_real_roots_squarefree() {
    //a test where real_roots_squarefree find rational roots while refining the bounds on a real root such that no rational root lies on the end points
    let poly = Polynomial::from_coeffs(vec![
        Integer::from(0),
        Integer::from(9),
        Integer::from(0),
        Integer::from(-4),
    ]);
    let roots = poly.real_roots_squarefree(
        Some(&Rational::from(-2)),
        Some(&Rational::from(2)),
        false,
        true,
    );
    println!("{:?}", roots);
    roots.check_invariants().unwrap();
}

#[test]
fn test_real_root_irreducible_count() {
    assert_eq!(
        Polynomial::from_coeffs(vec![
            Integer::from(3),
            Integer::from(-3),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(1)
        ])
        .real_roots_irreducible(None, None, false, false)
        .len(),
        1
    );
    assert_eq!(
        Polynomial::from_coeffs(vec![
            Integer::from(1),
            Integer::from(-3),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(1)
        ])
        .real_roots_irreducible(None, None, false, false)
        .len(),
        3
    );
}

#[test]
fn test_real_algebraic_ordering() {
    let mut all_roots = vec![];
    for f in vec![
        Polynomial::from_coeffs(vec![
            Integer::from(-2),
            Integer::from(-4),
            Integer::from(-2),
        ]),
        Polynomial::from_coeffs(vec![Integer::from(6), Integer::from(0), Integer::from(-3)]),
        Polynomial::from_coeffs(vec![Integer::from(1), Integer::from(-3), Integer::from(1)]),
        Polynomial::from_coeffs(vec![
            Integer::from(2),
            Integer::from(-3),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(1),
        ]),
        Polynomial::from_coeffs(vec![
            Integer::from(1),
            Integer::from(-3),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(1),
        ]),
        Polynomial::from_coeffs(vec![
            Integer::from(-1),
            Integer::from(12),
            Integer::from(-4),
            Integer::from(-15),
            Integer::from(5),
            Integer::from(3),
            Integer::from(-1),
        ]),
    ] {
        for root in Polynomial::real_roots(&f, None, None, false, false) {
            all_roots.push(root.clone());
        }
    }

    all_roots.sort();

    for mut root in &mut all_roots {
        root.check_invariants().unwrap();
        match &mut root {
            RealAlgebraic::Rational(_a) => {}
            RealAlgebraic::Real(a) => a.refine_to_accuracy(&Rational::from_signeds(1, i64::MAX)),
        }
        println!("    {} {:?}", root.to_string(), root);
    }

    let mut all_roots_sorted_by_lower_tight_bound = all_roots.clone();
    all_roots_sorted_by_lower_tight_bound.sort_by_key(|root| match root {
        RealAlgebraic::Rational(a) => a.clone(),
        RealAlgebraic::Real(r) => r.tight_a.clone(),
    });
    assert_eq!(all_roots, all_roots_sorted_by_lower_tight_bound);
}

#[test]
fn test_at_fixed_re_and_im() {
    let f = Polynomial::from_coeffs(vec![
        Integer::from(-1),
        Integer::from(3),
        Integer::from(0),
        Integer::from(1),
    ]);

    println!("f = {}", Polynomial::to_string(&f));

    let (vert_re_f, vert_im_f) = Polynomial::at_fixed_re(&f, &Rational::TWO);
    println!("re = {}", Polynomial::to_string(&vert_re_f));
    println!("im = {}", Polynomial::to_string(&vert_im_f));
    // f(z) = z^3 + 3z - 1
    // f(2 + xi) = (2 + xi)^3 + 3(2 + xi) - 1
    //           = 8 + 12xi - 6x^2 - x^3i + 6 + 3xi - 1
    //           = 13 + 15ix - 6x^2 - ix^3
    debug_assert_eq!(
        vert_re_f,
        Polynomial::from_coeffs(vec![Integer::from(13), Integer::from(0), Integer::from(-6)])
    );
    debug_assert_eq!(
        vert_im_f,
        Polynomial::from_coeffs(vec![
            Integer::from(0),
            Integer::from(15),
            Integer::from(0),
            Integer::from(-1)
        ])
    );

    let (vert_re_f, vert_im_f) = Polynomial::at_fixed_re(&f, &Rational::from_signeds(1, 2));
    println!("re = {}", Polynomial::to_string(&vert_re_f));
    println!("im = {}", Polynomial::to_string(&vert_im_f));
    // f(z) = z^3 + 3z - 1
    // f(1/2 + xi) = 5 + 30ix - 12x^2 - 8ix^3
    debug_assert_eq!(
        vert_re_f,
        Polynomial::from_coeffs(vec![Integer::from(5), Integer::from(0), Integer::from(-12)])
    );
    debug_assert_eq!(
        vert_im_f,
        Polynomial::from_coeffs(vec![
            Integer::from(0),
            Integer::from(30),
            Integer::from(0),
            Integer::from(-8)
        ])
    );

    let (vert_re_f, vert_im_f) = Polynomial::at_fixed_im(&f, &Rational::TWO);
    println!("re = {}", Polynomial::to_string(&vert_re_f));
    println!("im = {}", Polynomial::to_string(&vert_im_f));
    // f(z) = z^3 + 3z - 1
    // f(x + 2i) = -1 -2i -9x + 6ix^2 + x^3
    debug_assert_eq!(
        &vert_re_f,
        &Polynomial::from_coeffs(vec![
            Integer::from(-1),
            Integer::from(-9),
            Integer::from(0),
            Integer::from(1)
        ])
    );
    debug_assert_eq!(
        &vert_im_f,
        &Polynomial::from_coeffs(vec![Integer::from(-2), Integer::from(0), Integer::from(6),])
    );

    let (vert_re_f, vert_im_f) = Polynomial::at_fixed_im(&f, &Rational::from_signeds(1, 2));
    println!("re = {}", Polynomial::to_string(&vert_re_f));
    println!("im = {}", Polynomial::to_string(&vert_im_f));
    // f(z) = z^3 + 3z - 1
    // f(x + 1/2i) = -8 +11i + 18x + 12ix^2 + 8x^3
    debug_assert_eq!(
        vert_re_f,
        Polynomial::from_coeffs(vec![
            Integer::from(-8),
            Integer::from(18),
            Integer::from(0),
            Integer::from(8)
        ])
    );
    debug_assert_eq!(
        vert_im_f,
        Polynomial::from_coeffs(vec![Integer::from(11), Integer::from(0), Integer::from(12)])
    );
}

#[test]
fn test_count_complex_roots() {
    //cyclotomic polynomials in a box of sidelength 4
    for k in 1..19 {
        let f = Polynomial::add(
            &Polynomial::var_pow(k),
            &Polynomial::neg(&Polynomial::one()),
        );
        let n = f
            .count_complex_roots(
                &Rational::from(-2),
                &Rational::from(2),
                &Rational::from(-2),
                &Rational::from(2),
            )
            .unwrap();
        assert_eq!(n, k);
    }

    //cyclotomic polynomials in a box with a boundary root iff k=0 mod 2
    for k in 1..19 {
        let f = Polynomial::add(
            &Polynomial::var_pow(k),
            &Polynomial::neg(&Polynomial::one()),
        );
        let n = Polynomial::count_complex_roots(
            &f,
            &Rational::from(-1),
            &Rational::from(3),
            &Rational::from(-3),
            &Rational::from(3),
        );
        if k % 2 == 0 {
            assert!(n.is_none());
        } else {
            assert_eq!(n.unwrap(), k);
        }
    }

    //cyclotomic polynomials in a box with a boundary root iff k=0 mod 4
    for k in 1..19 {
        let f = Polynomial::add(
            &Polynomial::var_pow(k),
            &Polynomial::neg(&Polynomial::one()),
        );
        let n = Polynomial::count_complex_roots(
            &f,
            &Rational::from(-2),
            &Rational::from(2),
            &Rational::from(-1),
            &Rational::from(1),
        );
        if k % 4 == 0 {
            assert!(n.is_none());
        } else {
            assert_eq!(n.unwrap(), k);
        }
    }

    //other test cases
    assert_eq!(
        Some(1),
        Polynomial::count_complex_roots(
            &Polynomial::from_coeffs(vec![
                Integer::from(2),
                Integer::from(-8),
                Integer::from(1),
                Integer::from(-4),
                Integer::from(0),
                Integer::from(1),
            ]),
            &Rational::from(-1),
            &Rational::from(1),
            &Rational::from(1),
            &Rational::from(2),
        )
    );

    assert_eq!(
        Some(3),
        Polynomial::count_complex_roots(
            &Polynomial::from_coeffs(vec![
                Integer::from(2),
                Integer::from(-8),
                Integer::from(1),
                Integer::from(-4),
                Integer::from(0),
                Integer::from(1),
            ]),
            &Rational::from(-3),
            &Rational::from(3),
            &Rational::from(-1),
            &Rational::from(1),
        )
    );

    //polynomial with roots 2+3i, 2-3i and counting box with 2+3i as a vertex
    assert_eq!(
        None,
        Polynomial::count_complex_roots(
            &Polynomial::from_coeffs(vec![Integer::from(13), Integer::from(-4), Integer::from(1),]),
            &Rational::from(2),
            &Rational::from(3),
            &Rational::from(3),
            &Rational::from(4),
        )
    );

    //x^2-x+1
    let f = Polynomial::from_coeffs(vec![Integer::from(1), Integer::from(-1), Integer::from(1)]);
    let n = f
        .count_complex_roots(
            &Rational::from(-1),
            &Rational::from(1),
            &Rational::from(-1),
            &Rational::from(1),
        )
        .unwrap();
    assert_eq!(n, 2);
}

#[test]
fn test_count_complex_roots_big_values() {
    {
        let a = Rational::from_str("667/19382").unwrap(); //0.034413373232896505
        let b = Rational::from_str("754/4405").unwrap(); //0.17116912599318956
        let c = Rational::from_str("899/9691").unwrap(); //0.09276648436693839
        let d = Rational::from_str("3161/19382").unwrap(); //0.16308946445155298

        //0.951343405 -6.707838852x + 27.141574009x^2 = 0
        // x = 0.123571 - 0.140646 i
        // x = 0.123571 + 0.140646 i

        let poly = Polynomial::from_coeffs(vec![
            Integer::from_str("951343405").unwrap(),
            Integer::from_str("-6707838852").unwrap(),
            Integer::from_str("27141574009").unwrap(),
        ]);
        assert_eq!(poly.count_complex_roots(&a, &b, &c, &d), Some(1));
    }

    // {
    //     let a = Rational::from_str("-163099/329494").unwrap(); //-0.49499839147298585
    //     let b = Rational::from_str("-26827/74885").unwrap(); //-0.3582426387126928
    //     let c = Rational::from_str("899/9691").unwrap(); //0.09276648436693839
    //     let d = Rational::from_str("3161/19382").unwrap(); //0.16308946445155298

    //     //2139048288466 + 8191288386270x + 7843914888601x^2 = 0
    //     //same as
    //     //2.139048288466 + 8.191288386270x + 7.843914888601x^2 = 0
    //     //roots are -0.52 +/- a tiny complex bit

    //     let poly = Polynomial::from_coeffs(vec![
    //         Integer::from_str("2139048288466").unwrap(),
    //         Integer::from_str("8191288386270").unwrap(),
    //         Integer::from_str("7843914888601").unwrap(),
    //     ]);
    //     assert_eq!(poly.count_complex_roots(&a, &b, &c, &d), Some(1));
    // }
}

#[test]
fn test_real_neg() {
    {
        let f =
            Polynomial::from_coeffs(vec![Integer::from(-2), Integer::from(0), Integer::from(1)]);
        let roots = Polynomial::all_real_roots(&f);

        assert_eq!(roots.len(), 2);
        let a = &roots[0];
        let b = &roots[1];

        let a_neg = RealAlgebraic::neg(a);
        let b_neg = RealAlgebraic::neg(b);

        a_neg.check_invariants().unwrap();
        b_neg.check_invariants().unwrap();

        println!("a = {}", a.to_string());
        println!("b = {}", b.to_string());
        println!("a_neg = {}", a_neg.to_string());
        println!("b_neg = {}", b_neg.to_string());

        assert_ne!(a, b);
        assert_eq!(a, &b_neg);
        assert_eq!(b, &a_neg);
    }
    {
        let f = Polynomial::from_coeffs(vec![
            Integer::from(-1),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(3),
            Integer::from(1),
        ]);
        let roots = Polynomial::all_real_roots(&f);

        assert_eq!(roots.len(), 3);
        for root in roots {
            RealAlgebraic::neg(&root).check_invariants().unwrap();
        }
    }
    {
        //example where f(g(x)) is not primitive even though f and g are
        let f =
            Polynomial::from_coeffs(vec![Integer::from(-4), Integer::from(-1), Integer::from(1)]);
        let roots = Polynomial::all_real_roots(&f);
        for root in roots {
            let root2 = RealAlgebraic::add(
                &root,
                &RealAlgebraic::from_rat(&Rational::from_signeds(1, 2)).unwrap(),
            );
            root2.check_invariants().unwrap();
        }
    }
}

#[test]
fn test_real_add() {
    let f = Polynomial::from_coeffs(vec![Integer::from(-2), Integer::from(0), Integer::from(3)]);
    let roots = Polynomial::all_real_roots(&f);
    let a = RealAlgebraic::sum(roots.iter().collect());
    assert_eq!(a, RealAlgebraic::zero());

    let f = Polynomial::from_coeffs(vec![
        Integer::from(-7),
        Integer::from(0),
        Integer::from(100),
    ]);
    let roots = Polynomial::all_real_roots(&f);
    let a = RealAlgebraic::sum(roots.iter().collect());
    assert_eq!(a, RealAlgebraic::zero());

    let f = Polynomial::from_coeffs(vec![
        Integer::from(-100),
        Integer::from(0),
        Integer::from(7),
    ]);
    let roots = Polynomial::all_real_roots(&f);
    let a = RealAlgebraic::sum(roots.iter().collect());
    assert_eq!(a, RealAlgebraic::zero());
}

#[test]
fn test_real_mul() {
    let f = Polynomial::from_coeffs(vec![
        Integer::from(-100),
        Integer::from(0),
        Integer::from(7),
    ]);
    // (x-a)(x-b) = x^2 - 100/7
    // so ab=-100/7
    let roots = Polynomial::all_real_roots(&f);
    let a = RealAlgebraic::product(roots.iter().collect());
    assert_eq!(
        a,
        RealAlgebraic::from_rat(&Rational::from_signeds(-100, 7)).unwrap()
    );
}

#[test]
fn test_real_nth_root() {
    let x = &Polynomial::<Integer>::var().into_ring();
    let f = ((4 * x.pow(5) - 12 * x.pow(3) + 8 * x + 1)
        * (x + 1)
        * (x)
        * (x - 1)
        * (x - 2)
        * (x - 3)
        * (x - 4)
        * (x - 5)
        * (x - 144)
        * (x.pow(2) - 3))
        .into_set();
    let n = 2;

    for root in f.all_real_roots() {
        println!();
        println!("root = {}", root);
        match root.nth_root(n) {
            Ok(nth_root) => {
                println!("YES {}-root = {}", n, nth_root);
                debug_assert!(RealAlgebraic::zero() <= root);
                debug_assert!(RealAlgebraic::zero() <= nth_root);
                debug_assert_eq!(nth_root.nat_pow(&Natural::from(n)), root);
            }
            Err(()) => {
                println!("NO {}-root", n);
                debug_assert!(RealAlgebraic::zero() > root);
            }
        }
    }
}

#[test]
fn test_all_complex_roots() {
    let f = Polynomial::from_coeffs(vec![
        Integer::from(-1),
        Integer::from(-1),
        Integer::from(0),
        Integer::from(0),
        Integer::from(0),
        Integer::from(1),
    ]);
    let roots = Polynomial::all_complex_roots(&f);
    assert_eq!(roots.len(), Polynomial::degree(&f).unwrap());
    for root in &roots {
        root.check_invariants().unwrap();
    }
}

#[test]
fn test_complex_equal() {
    let x = &Polynomial::<Integer>::var().into_ring();
    let f = (x.pow(5) - x + 1).into_set();
    assert!(f.is_irreducible());
    let mut count = 0;
    for a in f.all_complex_roots() {
        for b in f.all_complex_roots() {
            if a == b {
                count += 1;
            }
        }
    }
    assert_eq!(count, f.degree().unwrap());
}

#[test]
fn test_complex_root_sum() {
    let f = Polynomial::from_coeffs(vec![
        Integer::from(1),
        Integer::from(0),
        Integer::from(0),
        Integer::from(1),
    ]);
    let roots = f.all_complex_roots();
    let s = ComplexAlgebraic::sum(roots.iter().collect());
    println!("{:?}", s);
    assert_eq!(s, ComplexAlgebraic::zero());

    let f = Polynomial::from_coeffs(vec![
        Integer::from(7),
        Integer::from(-3),
        Integer::from(42),
        Integer::from(9),
    ]);
    println!("f = {}", f);
    let roots = f.all_complex_roots();
    for root in &roots {
        println!("root = {}", root);
    }
    let s = ComplexAlgebraic::sum(roots.iter().collect());
    println!("root sum = {:?}", s);
    assert_eq!(
        s,
        ComplexAlgebraic::Real(RealAlgebraic::Rational(Rational::from_integers(
            Integer::from(-14),
            Integer::from(3)
        )))
    );
}

#[test]
fn test_complex_mul() {
    let f = Polynomial::from_coeffs(vec![
        Integer::from(1),
        Integer::from(0),
        Integer::from(0),
        Integer::from(1),
    ]);
    let roots = f.all_complex_roots();
    let s = ComplexAlgebraic::product(roots.iter().collect());
    println!("{:?}", s);
    assert_eq!(s, ComplexAlgebraic::one().neg());

    let f = Polynomial::from_coeffs(vec![
        Integer::from(7),
        Integer::from(-3),
        Integer::from(42),
        Integer::from(9),
    ]);
    println!("f = {}", f);
    let roots = f.all_complex_roots();
    for root in &roots {
        println!("root = {}", root);
    }
    let s = ComplexAlgebraic::product(roots.iter().collect());
    println!("root prod = {:?}", s);
    assert_eq!(
        s,
        ComplexAlgebraic::Real(RealAlgebraic::Rational(Rational::from_integers(
            Integer::from(-7),
            Integer::from(9)
        )))
    );
}

#[test]
fn test_complex_inv() {
    let x = &Polynomial::<Integer>::var().into_ring();
    let f = (x.pow(4) - x + 1).into_set();

    for root in f.all_complex_roots() {
        assert_eq!(
            ComplexAlgebraic::mul(&root.inv().unwrap(), &root),
            ComplexAlgebraic::one()
        );
    }
}

#[test]
fn test_apply_poly() {
    let x = &Polynomial::<Rational>::var().into_ring();

    let f = (2 * x.pow(2) - 4 * x - 3).into_set();
    let g = (3 * x.pow(3) + 7 * x - 1).into_set();
    //h(x) = f(g(x))
    let h = Polynomial::compose(&f, &g);

    println!("f = {}", f);
    println!("g = {}", g);
    println!("h = {}", h);

    for x in h.primitive_part_fof().all_complex_roots() {
        println!("");
        println!("x = {} deg = {}", x, x.min_poly().degree().unwrap());
        let gx = x.clone().apply_poly(&g);
        println!("gx = {}", gx);
        let fgx = gx.clone().apply_poly(&f);
        println!("fgx = {}", fgx);
        debug_assert_eq!(fgx, ComplexAlgebraic::zero());
    }

    let i = &ComplexAlgebraic::i().into_ring();
    let a = (2 + 3 * i).into_set();
    let f = (x.pow(10)).into_set();
    assert_eq!(a.clone().apply_poly(&f), (-341525 - 145668 * i).into_set());
}

#[test]
fn test_as_poly_expr() {
    let x = &Polynomial::<Rational>::var().into_ring();

    let two = ComplexAlgebraic::Real(RealAlgebraic::Rational(Rational::from(2)));
    let sqrt_two = two.nth_root(2).unwrap();
    let fourrt_two = sqrt_two.nth_root(2).unwrap();
    assert_eq!(
        as_poly_expr(&sqrt_two, &fourrt_two),
        Some(x.pow(2).into_set())
    );
    assert_eq!(as_poly_expr(&fourrt_two, &sqrt_two), None);

    let f = (x.pow(3) - 2 * x.pow(2) - 2 * x - 1).into_set();
    for root in f.primitive_part_fof().all_complex_roots() {
        assert_eq!(as_poly_expr(&root, &root), Some(x.pow(1).into_set()));
    }
}

#[test]
fn test_pair_generated_anf() {
    // let x = &Polynomial::<Rational>::var().into_ring();

    let sqrt_two = ComplexAlgebraic::Real(RealAlgebraic::Rational(Rational::from(2)))
        .nth_root(2)
        .unwrap();
    let sqrt_three = ComplexAlgebraic::Real(RealAlgebraic::Rational(Rational::from(3)))
        .nth_root(2)
        .unwrap();
    let sqrt_six = ComplexAlgebraic::Real(RealAlgebraic::Rational(Rational::from(6)))
        .nth_root(2)
        .unwrap();

    println!("{}", sqrt_two);
    println!("{}", sqrt_three);
    println!("{}", sqrt_six);

    let (gen, _, _, x, y) = anf_pair_primitive_element_theorem(&sqrt_two, &sqrt_three);
    println!(
        "gen = {} min_poly = {} deg = {}",
        gen,
        gen.min_poly(),
        gen.min_poly().degree().unwrap()
    );
    assert_eq!(
        gen,
        (sqrt_two.clone().into_ring() + sqrt_three.clone().into_ring()).into_set()
    );
    let oof = (sqrt_two.clone().into_ring() + sqrt_three.clone().into_ring()).into_set();
    println!("{} {}", oof, oof.min_poly());
    println!("x = {}", x);
    println!("y = {}", y);
}

#[test]
fn test_multi_generated_anf() {
    // let x = &Polynomial::<Rational>::var().into_ring();

    let sqrt_two = ComplexAlgebraic::Real(RealAlgebraic::Rational(Rational::from(2)))
        .nth_root(2)
        .unwrap();
    let sqrt_three = ComplexAlgebraic::Real(RealAlgebraic::Rational(Rational::from(3)))
        .nth_root(2)
        .unwrap();
    let sqrt_six = ComplexAlgebraic::Real(RealAlgebraic::Rational(Rational::from(6)))
        .nth_root(2)
        .unwrap();

    println!("{}", sqrt_two);
    println!("{}", sqrt_three);
    println!("{}", sqrt_six);

    let (mut g, p) = anf_multi_primitive_element_theorem(vec![&sqrt_two, &sqrt_three, &sqrt_six]);

    println!("{} of degree {}", g, g.degree());
    println!("sqrt(2) = {} applied to {}", p[0], g);
    println!("sqrt(3) = {} applied to {}", p[1], g);
    println!("sqrt(6) = {} applied to {}", p[2], g);

    assert_eq!(sqrt_two, g.apply_poly(&p[0]));
    assert_eq!(sqrt_three, g.apply_poly(&p[1]));
    assert_eq!(sqrt_six, g.apply_poly(&p[2]));
}
