use algebraeon_nzq::traits::Fraction;

use super::*;

impl Polynomial<Integer> {
    fn at_fixed_re_or_im_impl<const RE_OR_IM: bool>(
        &self,
        a: &Rational,
    ) -> (Polynomial<Integer>, Polynomial<Integer>) {
        //find real and imag polys of
        //poly(a + xi) if RE_OR_IM = false
        //poly(x + ai) if RE_OR_IM = true
        //up to rational multiples (its the roots we care about)
        match self.degree() {
            Some(n) => {
                let (a_numer, a_denom) = a.numerator_and_denominator();
                let a_denom = Integer::from(a_denom);
                //multiply everything by a_d^n so that everything is integers

                //compute 1, a, a^2, a^3, ..., a^n (after multiplying everything by a_d)
                // a_d^n(a_n/a_d)^k = a_n^k a_d^{n-k}
                let mut a_numer_pow = vec![Integer::from(1)];
                let mut a_denom_pow = vec![Integer::from(1)];
                for k in 1..=n {
                    a_numer_pow.push(&a_numer * &a_numer_pow[k - 1]);
                    a_denom_pow.push(&a_denom * &a_denom_pow[k - 1]);
                }
                let mut a_pow = vec![];
                for k in 0..=n {
                    a_pow.push(&a_numer_pow[k] * &a_denom_pow[n - k]);
                }

                let mut re = Vec::with_capacity(n + 1);
                let mut im = Vec::with_capacity(n + 1);
                for _ in 0..=n {
                    re.push(Integer::from(0));
                    im.push(Integer::from(0));
                }
                let mut n_choose = vec![Integer::from(1)];
                for n in 0..=n {
                    if n == 0 {
                        debug_assert_eq!(n_choose, vec![Integer::from(1)]);
                    } else if n == 1 {
                        debug_assert_eq!(n_choose, vec![Integer::from(1), Integer::from(1)]);
                    } else if n == 2 {
                        debug_assert_eq!(
                            n_choose,
                            vec![Integer::from(1), Integer::from(2), Integer::from(1)]
                        );
                    } else if n == 3 {
                        debug_assert_eq!(
                            n_choose,
                            vec![
                                Integer::from(1),
                                Integer::from(3),
                                Integer::from(3),
                                Integer::from(1)
                            ]
                        );
                    }

                    //if fixed real add
                    //(a + xi)^n = \sum_{k=0,1,...,n} \binom{n}{k} a^{n-k} (xi)^k
                    //           = \sum_{k=0,1,...,n} \binom{n}{k} a^{n-k} x^k i^k
                    //           = \sum_{k=0,1,...,n} {
                    //               k = 0 mod 4        + \binom{n}{k} a^{n-k} x^k
                    //               k = 1 mod 4        + \binom{n}{k} a^{n-k} x^k i
                    //               k = 2 mod 4        - \binom{n}{k} a^{n-k} x^k
                    //               k = 3 mod 4        - \binom{n}{k} a^{n-k} x^k i
                    //                                }
                    //
                    //if fixed imag add
                    //(a + xi)^n = \sum_{k=0,1,...,n} \binom{n}{k} a^{n-k} (xi)^k
                    //           = \sum_{k=0,1,...,n} \binom{n}{k} a^{n-k} x^k i^k
                    //           = \sum_{k=0,1,...,n} {
                    //               k = 0 mod 4        + \binom{n}{k} a^{n-k} x^k
                    //               k = 1 mod 4        + \binom{n}{k} a^{n-k} x^k i
                    //               k = 2 mod 4        - \binom{n}{k} a^{n-k} x^k
                    //               k = 3 mod 4        - \binom{n}{k} a^{n-k} x^k i
                    //
                    if self.coeff(n).as_ref() != &Integer::ZERO {
                        let mut k = 0;
                        loop {
                            //k = 0 mod 4
                            re[if RE_OR_IM { n - k } else { k }] += self.coeff(n).as_ref()
                                * &n_choose[k]
                                * &a_pow[if RE_OR_IM { k } else { n - k }];
                            if k == n {
                                break;
                            }
                            k += 1;
                            //k = 1 mod 4
                            im[if RE_OR_IM { n - k } else { k }] += self.coeff(n).as_ref()
                                * &n_choose[k]
                                * &a_pow[if RE_OR_IM { k } else { n - k }];
                            if k == n {
                                break;
                            }
                            k += 1;
                            //k = 2 mod 4
                            re[if RE_OR_IM { n - k } else { k }] -= self.coeff(n).as_ref()
                                * &n_choose[k]
                                * &a_pow[if RE_OR_IM { k } else { n - k }];
                            if k == n {
                                break;
                            }
                            k += 1;
                            //k = 3 mod 4
                            im[if RE_OR_IM { n - k } else { k }] -= self.coeff(n).as_ref()
                                * &n_choose[k]
                                * &a_pow[if RE_OR_IM { k } else { n - k }];
                            if k == n {
                                break;
                            }
                            k += 1;
                        }
                    }
                    //update n choose k
                    //e.g. for n=3 do
                    //[1, 3, 3, 1]
                    //[1, 3, 3, 1, 1]
                    //[1, 3, 3, 4, 1]
                    //[1, 3, 6, 4, 1]
                    //[1, 4, 6, 4, 1]
                    n_choose.push(Integer::from(1));
                    for i in (1..=n).rev() {
                        n_choose[i] = &n_choose[i] + &n_choose[i - 1];
                    }
                }
                (Polynomial::from_coeffs(re), Polynomial::from_coeffs(im))
            }
            None => (Self::zero(), Self::zero()),
        }
    }

    fn at_fixed_re(&self, a: &Rational) -> (Polynomial<Integer>, Polynomial<Integer>) {
        self.at_fixed_re_or_im_impl::<false>(a)
    }

    fn at_fixed_im(&self, a: &Rational) -> (Polynomial<Integer>, Polynomial<Integer>) {
        self.at_fixed_re_or_im_impl::<true>(a)
    }

    //count how many complex roots are in the box a < re < b, c < im < d
    //or return None if there is a root on the boundary
    pub fn count_complex_roots(
        &self,
        a: &Rational,
        b: &Rational,
        c: &Rational,
        d: &Rational,
    ) -> Option<usize> {
        assert!(a < b);
        assert!(c < d);

        //the idea is to compute the winding number of the path around the boundary of the box
        //this is done by computing where the value of the polynomial crosses the real and imaginary axes as the input traces the path
        //the crossing points and their order is done using the exact total ordering of real polynomial roots
        let (a_vert_re, a_vert_im) = self.at_fixed_re(a);
        let (b_vert_re, b_vert_im) = self.at_fixed_re(b);
        let (c_horz_re, c_horz_im) = self.at_fixed_im(c);
        let (d_horz_re, d_horz_im) = self.at_fixed_im(d);

        // println!("poly = {} abcd = {} {} {} {}", &self, a, b, c, d);
        // println!(
        //     "a_vert_re = {}, a_vert_im = {}",
        //     Polynomial::to_string(&a_vert_re),
        //     Polynomial::to_string(&a_vert_im)
        // );
        // println!(
        //     "b_vert_re = {}, b_vert_im = {}",
        //     Polynomial::to_string(&b_vert_re),
        //     Polynomial::to_string(&b_vert_im)
        // );
        // println!(
        //     "c_horz_re = {}, c_horz_im = {}",
        //     Polynomial::to_string(&c_horz_re),
        //     Polynomial::to_string(&c_horz_im)
        // );
        // println!(
        //     "d_horz_re = {}, d_horz_im = {}",
        //     Polynomial::to_string(&d_horz_re),
        //     Polynomial::to_string(&d_horz_im)
        // );

        // //checks will fail - the real and imaginary parts are only up to scalar multiples
        // debug_assert_eq!(
        //     evaluate_at_rational(&a_vert_re, c),
        //     evaluate_at_rational(&c_horz_re, a)
        // );
        // debug_assert_eq!(
        //     evaluate_at_rational(&a_vert_re, d),
        //     evaluate_at_rational(&d_horz_re, a)
        // );
        // debug_assert_eq!(
        //     evaluate_at_rational(&b_vert_re, c),
        //     evaluate_at_rational(&c_horz_re, b)
        // );
        // debug_assert_eq!(
        //     evaluate_at_rational(&b_vert_re, d),
        //     evaluate_at_rational(&d_horz_re, b)
        // );
        // debug_assert_eq!(
        //     evaluate_at_rational(&a_vert_im, c),
        //     evaluate_at_rational(&c_horz_im, a)
        // );
        // debug_assert_eq!(
        //     evaluate_at_rational(&a_vert_im, d),
        //     evaluate_at_rational(&d_horz_im, a)
        // );
        // debug_assert_eq!(
        //     evaluate_at_rational(&b_vert_im, c),
        //     evaluate_at_rational(&c_horz_im, b)
        // );
        // debug_assert_eq!(
        //     evaluate_at_rational(&b_vert_im, d),
        //     evaluate_at_rational(&d_horz_im, b)
        // );

        //compute squarefree versions because when only care about the roots without multiplicity
        let a_vert_re_sqfr = a_vert_re.clone().primitive_squarefree_part();
        let a_vert_im_sqfr = a_vert_im.clone().primitive_squarefree_part();
        let b_vert_re_sqfr = b_vert_re.clone().primitive_squarefree_part();
        let b_vert_im_sqfr = b_vert_im.clone().primitive_squarefree_part();
        let c_horz_re_sqfr = c_horz_re.clone().primitive_squarefree_part();
        let c_horz_im_sqfr = c_horz_im.clone().primitive_squarefree_part();
        let d_horz_re_sqfr = d_horz_re.clone().primitive_squarefree_part();
        let d_horz_im_sqfr = d_horz_im.clone().primitive_squarefree_part();

        //trace an anticlockwise path around the box and create a list of crossings which encode what happens to the value of the polynomial
        #[derive(Debug)]
        enum Crossing {
            PosRe,
            PosIm,
            NegRe,
            NegIm,
        }

        fn crossings<const REVERSE: bool>(
            re: &Polynomial<Integer>,
            mut re_sqfr: Polynomial<Integer>,
            im: &Polynomial<Integer>,
            mut im_sqfr: Polynomial<Integer>,
            s: &Rational,
            t: &Rational,
        ) -> Option<Vec<Crossing>> {
            use super::super::real::polynomial::*;
            // println!(
            //     "REVERSE={} re={}, re_sqfr={}, im={}, im_sqfr={}",
            //     REVERSE,
            //     Polynomial::to_string(&re),
            //     Polynomial::to_string(&re_sqfr),
            //     Polynomial::to_string(&im),
            //     Polynomial::to_string(&im_sqfr)
            // );
            debug_assert_eq!(re == &Polynomial::zero(), re_sqfr == Polynomial::zero());
            debug_assert_eq!(im == &Polynomial::zero(), im_sqfr == Polynomial::zero());
            //because if the real and imaginary part are both constant at 0 then poly has infinitely many complex zeros which is not possible
            debug_assert!(re_sqfr != Polynomial::zero() || im_sqfr != Polynomial::zero());
            if re_sqfr == Polynomial::zero() {
                //the image is doing a path confied to the imaginary axis
                let roots_im = im.real_roots(Some(s), Some(t), true, true);
                if roots_im.is_empty() {
                    //the image stays once side of the real axis
                    let val = evaluate_at_rational(im, s);
                    debug_assert_eq!(
                        val > Rational::ZERO,
                        evaluate_at_rational(im, t) > Rational::ZERO
                    );
                    if val > Rational::ZERO {
                        Some(vec![Crossing::PosIm]) //this whole line segment is a positive imaginary crossing
                    } else {
                        Some(vec![Crossing::NegIm]) //this whole line segment is a negative imaginary crossing
                    }
                } else {
                    //the image crosses the real axis and hence passes through 0
                    None
                }
            } else if im_sqfr == Polynomial::zero() {
                //the image is doing a path confied to the real axis
                let roots_re = re.real_roots(Some(s), Some(t), true, true);
                if roots_re.is_empty() {
                    //the image stays one side of the imaginary axis
                    let val = evaluate_at_rational(re, s);
                    debug_assert_eq!(
                        val > Rational::ZERO,
                        evaluate_at_rational(re, t) > Rational::ZERO
                    );
                    if val > Rational::ZERO {
                        Some(vec![Crossing::PosRe]) //this whole line segment is a positive real crossing
                    } else {
                        Some(vec![Crossing::NegRe]) //this whole line segment is a negative real crossing
                    }
                } else {
                    //the image crosses the imaginary axis and hence passes through 0
                    None
                }
            } else {
                //want to isolate roots of squarefree polynomials without factoring
                //get ordered real roots in some structure
                //get ordered imag roots in some structure
                //do a merge sort pass to interleave the real and imaginary roots in the correct order
                //    if a real root equals an imaginary root then there is a root on the boundary
                //for each root of one type, compute the sign of the other part when evaluated at the root

                let mut crossings = vec![];

                //check the value of the real and imaginary part at the vertex at the start of this path
                let v = { if REVERSE { t } else { s } };
                match (
                    evaluate_at_rational(re, v).cmp(&Rational::from(0)),
                    evaluate_at_rational(im, v).cmp(&Rational::from(0)),
                ) {
                    (std::cmp::Ordering::Equal, std::cmp::Ordering::Equal) => {
                        //the polynomial is zero at vertex v
                        return None;
                    }
                    (std::cmp::Ordering::Equal, std::cmp::Ordering::Less) => {
                        crossings.push(Crossing::NegIm);
                    }
                    (std::cmp::Ordering::Equal, std::cmp::Ordering::Greater) => {
                        crossings.push(Crossing::PosIm);
                    }
                    (std::cmp::Ordering::Less, std::cmp::Ordering::Equal) => {
                        crossings.push(Crossing::NegRe);
                    }
                    (std::cmp::Ordering::Greater, std::cmp::Ordering::Equal) => {
                        crossings.push(Crossing::PosRe);
                    }
                    (_, _) => {}
                }

                if evaluate_at_rational(&re_sqfr, s) == Rational::from(0) {
                    re_sqfr = Polynomial::try_div(
                        &re_sqfr,
                        &Polynomial::from_coeffs(vec![-s.numerator(), s.denominator().into()]),
                    )
                    .unwrap();
                }

                if evaluate_at_rational(&re_sqfr, t) == Rational::from(0) {
                    re_sqfr = Polynomial::try_div(
                        &re_sqfr,
                        &Polynomial::from_coeffs(vec![-t.numerator(), t.denominator().into()]),
                    )
                    .unwrap();
                }

                if evaluate_at_rational(&im_sqfr, s) == Rational::from(0) {
                    im_sqfr = Polynomial::try_div(
                        &im_sqfr,
                        &Polynomial::from_coeffs(vec![-s.numerator(), s.denominator().into()]),
                    )
                    .unwrap();
                }
                if evaluate_at_rational(&im_sqfr, t) == Rational::from(0) {
                    im_sqfr = Polynomial::try_div(
                        &im_sqfr,
                        &Polynomial::from_coeffs(vec![-t.numerator(), t.denominator().into()]),
                    )
                    .unwrap();
                }
                debug_assert_ne!(evaluate_at_rational(&re_sqfr, s), Rational::from(0));
                debug_assert_ne!(evaluate_at_rational(&re_sqfr, t), Rational::from(0));
                debug_assert_ne!(evaluate_at_rational(&im_sqfr, s), Rational::from(0));
                debug_assert_ne!(evaluate_at_rational(&im_sqfr, t), Rational::from(0));

                let mut re_roots =
                    re_sqfr
                        .clone()
                        .real_roots_squarefree(Some(s), Some(t), REVERSE, !REVERSE);
                let mut im_roots =
                    im_sqfr
                        .clone()
                        .real_roots_squarefree(Some(s), Some(t), REVERSE, !REVERSE);

                // println!("re_roots = {:?}", re_roots);
                // println!("im_roots = {:?}", im_roots);

                debug_assert!(re_roots.check_invariants().is_ok());
                debug_assert!(im_roots.check_invariants().is_ok());

                match SquarefreePolyRealRoots::separate(&mut re_roots, &mut im_roots) {
                    Ok(mut all_roots) => {
                        //the isolating intervals for re_roots and im_roots no longer overlap
                        //we can use this to our advantage...

                        for (interleave, root_idx) in {
                            if REVERSE {
                                all_roots.reverse();
                                all_roots
                            } else {
                                all_roots
                            }
                        } {
                            // println!("interleave = {:?} root_idx = {:?}", interleave, root_idx);
                            match interleave {
                                Interleave::First => {
                                    //a real root
                                    loop {
                                        let re_root = re_roots.interval(root_idx);
                                        match re_root {
                                            SquarefreePolyRealRootInterval::Rational(x) => {
                                                match evaluate_at_rational(im, x)
                                                    .cmp(&Rational::from(0))
                                                {
                                                    std::cmp::Ordering::Less => {
                                                        crossings.push(Crossing::NegIm);
                                                        break;
                                                    }
                                                    std::cmp::Ordering::Equal => panic!(),
                                                    std::cmp::Ordering::Greater => {
                                                        crossings.push(Crossing::PosIm);
                                                        break;
                                                    }
                                                }
                                            }
                                            SquarefreePolyRealRootInterval::Real(a, b, _) => {
                                                match evaluate_at_rational(im, a)
                                                    .cmp(&Rational::from(0))
                                                {
                                                    std::cmp::Ordering::Less => {
                                                        crossings.push(Crossing::NegIm);
                                                        break;
                                                    }
                                                    std::cmp::Ordering::Equal => {
                                                        //need to refine
                                                    }
                                                    std::cmp::Ordering::Greater => {
                                                        crossings.push(Crossing::PosIm);
                                                        break;
                                                    }
                                                }
                                                match evaluate_at_rational(im, b)
                                                    .cmp(&Rational::from(0))
                                                {
                                                    std::cmp::Ordering::Less => {
                                                        crossings.push(Crossing::NegIm);
                                                        break;
                                                    }
                                                    std::cmp::Ordering::Equal => {
                                                        //need to refine
                                                    }
                                                    std::cmp::Ordering::Greater => {
                                                        crossings.push(Crossing::PosIm);
                                                        break;
                                                    }
                                                }
                                            }
                                        }
                                        re_roots.refine(root_idx);
                                    }
                                }
                                Interleave::Second => {
                                    //an imaginary root
                                    loop {
                                        let im_root = im_roots.interval(root_idx);
                                        match im_root {
                                            SquarefreePolyRealRootInterval::Rational(x) => {
                                                match evaluate_at_rational(re, x)
                                                    .cmp(&Rational::from(0))
                                                {
                                                    std::cmp::Ordering::Less => {
                                                        crossings.push(Crossing::NegRe);
                                                        break;
                                                    }
                                                    std::cmp::Ordering::Equal => panic!(),
                                                    std::cmp::Ordering::Greater => {
                                                        crossings.push(Crossing::PosRe);
                                                        break;
                                                    }
                                                }
                                            }
                                            SquarefreePolyRealRootInterval::Real(a, b, _) => {
                                                match evaluate_at_rational(re, a)
                                                    .cmp(&Rational::from(0))
                                                {
                                                    std::cmp::Ordering::Less => {
                                                        crossings.push(Crossing::NegRe);
                                                        break;
                                                    }
                                                    std::cmp::Ordering::Equal => {
                                                        //need to refine
                                                    }
                                                    std::cmp::Ordering::Greater => {
                                                        crossings.push(Crossing::PosRe);
                                                        break;
                                                    }
                                                }
                                                match evaluate_at_rational(re, b)
                                                    .cmp(&Rational::from(0))
                                                {
                                                    std::cmp::Ordering::Less => {
                                                        crossings.push(Crossing::NegRe);
                                                        break;
                                                    }
                                                    std::cmp::Ordering::Equal => {
                                                        //need to refine
                                                    }
                                                    std::cmp::Ordering::Greater => {
                                                        crossings.push(Crossing::PosRe);
                                                        break;
                                                    }
                                                }
                                            }
                                        }
                                        im_roots.refine(root_idx);
                                    }
                                }
                            }
                        }
                        Some(crossings)
                    }
                    Err(()) => None,
                }
            }
        }

        /*
             a           b

         d   +-----------+
             |           |
             |           |
             |           |
             |           |
             |           |
         c   +-----------+

        */

        // println!("c = {:?}", crossings::<false>(&c_horz_re, c_horz_re_sqfr.clone(), &c_horz_im, c_horz_im_sqfr.clone(), a, b));
        // println!("b = {:?}", crossings::<false>(&b_vert_re, b_vert_re_sqfr.clone(), &b_vert_im, b_vert_im_sqfr.clone(), c, d));
        // println!("d = {:?}", crossings::<true>(&d_horz_re, d_horz_re_sqfr.clone(), &d_horz_im, d_horz_im_sqfr.clone(), a, b));
        // println!("a = {:?}", crossings::<true>(&a_vert_re, a_vert_re_sqfr.clone(), &a_vert_im, a_vert_im_sqfr.clone(), c, d));

        let mut winding = vec![];
        for cr in [
            crossings::<false>(&c_horz_re, c_horz_re_sqfr, &c_horz_im, c_horz_im_sqfr, a, b),
            crossings::<false>(&b_vert_re, b_vert_re_sqfr, &b_vert_im, b_vert_im_sqfr, c, d),
            crossings::<true>(&d_horz_re, d_horz_re_sqfr, &d_horz_im, d_horz_im_sqfr, a, b),
            crossings::<true>(&a_vert_re, a_vert_re_sqfr, &a_vert_im, a_vert_im_sqfr, c, d),
        ] {
            match cr {
                Some(mut w) => winding.append(&mut w),
                None => {
                    return None;
                }
            }
        }

        // println!("winding = {:?}", winding);

        //compute the winding number = number of roots
        if winding.is_empty() {
            Some(0)
        } else {
            fn axis_pair_to_num_offset(ax1: &Crossing, ax2: &Crossing) -> isize {
                match (ax1, ax2) {
                    (Crossing::PosRe, Crossing::PosRe)
                    | (Crossing::PosIm, Crossing::PosIm)
                    | (Crossing::NegRe, Crossing::NegRe)
                    | (Crossing::NegIm, Crossing::NegIm) => 0,
                    (Crossing::PosRe, Crossing::NegRe)
                    | (Crossing::NegRe, Crossing::PosRe)
                    | (Crossing::PosIm, Crossing::NegIm)
                    | (Crossing::NegIm, Crossing::PosIm) => panic!(),
                    (Crossing::PosRe, Crossing::PosIm)
                    | (Crossing::PosIm, Crossing::NegRe)
                    | (Crossing::NegRe, Crossing::NegIm)
                    | (Crossing::NegIm, Crossing::PosRe) => 1,
                    (Crossing::PosRe, Crossing::NegIm)
                    | (Crossing::NegIm, Crossing::NegRe)
                    | (Crossing::NegRe, Crossing::PosIm)
                    | (Crossing::PosIm, Crossing::PosRe) => -1,
                }
            }

            let mut num: isize = 0;
            num += axis_pair_to_num_offset(&winding[winding.len() - 1], &winding[0]);
            for i in 0..winding.len() - 1 {
                num += axis_pair_to_num_offset(&winding[i], &winding[i + 1]);
            }

            assert!(num >= 0, "winding should always be overall anti-clockwise");
            let num = num as usize;
            match num % 4 {
                0 => Some(num / 4),
                _ => panic!("invalid remainder modulo four"),
            }
        }
    }

    fn uhp_complex_roots_irreducible_impl(
        &self,
        num_real_roots: usize,
    ) -> Vec<ComplexAlgebraicRoot> {
        debug_assert!(self.is_irreducible());
        debug_assert_eq!(self.all_real_roots().len(), num_real_roots);
        let deg = self.degree().unwrap();
        debug_assert!(num_real_roots <= deg);
        if num_real_roots == deg {
            vec![]
        } else {
            //search the upper half plane for the complete roots with positive imaginary part
            debug_assert_eq!((deg - num_real_roots) % 2, 0);
            let target_uhp_num = (deg - num_real_roots) / 2;

            let mut a = Rational::from(-1);
            let mut b = Rational::from(1);
            let mut c = Rational::from_integers(1, 2);
            let mut d = Rational::from(2);

            loop {
                if let Some(n) = self.count_complex_roots(&a, &b, &c, &d) {
                    debug_assert!(n <= target_uhp_num);
                    if n == target_uhp_num {
                        break;
                    }
                } else {
                    // boundary root
                }
                a *= Rational::from(2);
                b *= Rational::from(2);
                c *= Rational::from_integers(1, 2);
                d *= Rational::from(2);
            }

            fn bisect(
                poly: &Polynomial<Integer>,
                n: usize,
                a: &Rational,
                b: &Rational,
                c: &Rational,
                d: &Rational,
            ) -> Vec<ComplexAlgebraicRoot> {
                debug_assert!(a < b);
                debug_assert!(c < d);
                debug_assert_eq!(poly.count_complex_roots(a, b, c, d).unwrap(), n);
                if n == 0 {
                    vec![]
                } else if n == 1 {
                    vec![ComplexAlgebraicRoot {
                        poly: poly.clone(),
                        tight_a: a.clone(),
                        tight_b: b.clone(),
                        tight_c: c.clone(),
                        tight_d: d.clone(),
                    }]
                } else {
                    let ((n1, a1, b1, c1, d1), (n2, a2, b2, c2, d2)) =
                        bisect_box(poly, n, a, b, c, d);

                    let mut roots = bisect(poly, n1, &a1, &b1, &c1, &d1);
                    roots.append(&mut bisect(poly, n2, &a2, &b2, &c2, &d2));
                    roots
                }
            }

            bisect(self, target_uhp_num, &a, &b, &c, &d)
        }
    }

    fn lhp_complex_roots_irreducible_impl(
        &self,
        num_real_roots: usize,
    ) -> Vec<ComplexAlgebraicRoot> {
        #[allow(clippy::redundant_closure_for_method_calls)]
        self.uhp_complex_roots_irreducible_impl(num_real_roots)
            .into_iter()
            .map(|root| root.conj())
            .collect()
    }

    pub fn all_complex_roots_irreducible(&self) -> Vec<ComplexAlgebraic> {
        debug_assert!(self.is_irreducible());
        let deg = self.degree().unwrap();

        let mut all_roots = vec![];
        for real_root in self.all_real_roots() {
            all_roots.push(ComplexAlgebraic::Real(real_root));
        }
        let num_real_roots = all_roots.len();

        debug_assert!(num_real_roots <= deg);
        if num_real_roots == deg {
            return all_roots;
        }

        for complex_root in self.uhp_complex_roots_irreducible_impl(num_real_roots) {
            all_roots.push(ComplexAlgebraic::Complex(complex_root.clone().conj()));
            all_roots.push(ComplexAlgebraic::Complex(complex_root));
        }

        debug_assert_eq!(all_roots.len(), deg);
        all_roots
    }

    pub fn all_complex_roots(&self) -> Vec<ComplexAlgebraic> {
        assert_ne!(self, &Self::zero());
        let factors = self.factor();
        let mut roots = vec![];
        for (factor, k) in factors.into_powers().unwrap() {
            for root in factor.all_complex_roots_irreducible() {
                let mut i = Natural::from(0u8);
                while i < k {
                    roots.push(root.clone());
                    i += Natural::from(1u8);
                }
            }
        }
        roots
    }
}

impl Polynomial<Rational> {
    pub fn all_complex_roots(&self) -> Vec<ComplexAlgebraic> {
        assert_ne!(self, &Self::zero());
        self.primitive_part_fof().all_complex_roots()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

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

        let (vert_re_f, vert_im_f) = Polynomial::at_fixed_re(&f, &Rational::from_integers(1, 2));
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

        let (vert_re_f, vert_im_f) = Polynomial::at_fixed_im(&f, &Rational::from_integers(1, 2));
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
                &Polynomial::from_coeffs(vec![
                    Integer::from(13),
                    Integer::from(-4),
                    Integer::from(1),
                ]),
                &Rational::from(2),
                &Rational::from(3),
                &Rational::from(3),
                &Rational::from(4),
            )
        );

        //x^2-x+1
        let f =
            Polynomial::from_coeffs(vec![Integer::from(1), Integer::from(-1), Integer::from(1)]);
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
    }
}
