use malachite_q::Rational;

pub mod bisection_gen;
pub mod complex;
pub mod padic;
pub mod poly_tools;
pub mod real;

fn rat_to_string(a: Rational) -> String {
    if a == 0 {
        return "0".into();
    }
    let neg = a < Rational::from(0);
    let (mant, exp, _): (f64, _, _) = a
        .sci_mantissa_and_exponent_round(malachite_base::rounding_modes::RoundingMode::Nearest)
        .unwrap();
    let mut b = (2.0 as f64).powf(exp as f64) * mant;
    if neg {
        b = -b;
    }
    b = (1000.0 * b).round() / 1000.0;

    b.to_string()
}
