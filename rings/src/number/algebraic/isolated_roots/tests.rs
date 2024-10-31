use malachite_nz::{integer::Integer, natural::Natural};
use malachite_q::Rational;

use crate::{
    number::algebraic::isolated_roots::{
        anf_multi_primitive_element_theorem, anf_pair_primitive_element_theorem, as_poly_expr,
        ComplexAlgebraic, RealAlgebraic,
    },
    polynomial::polynomial::Polynomial,
    ring_structure::cannonical::{PositiveRealNthRootOpps, Ring, UniqueFactorizationDomain},
};

use super::{IntegralDomain, StructuredType};

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
