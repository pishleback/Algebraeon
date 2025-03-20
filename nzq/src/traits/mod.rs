pub trait DivMod<T> {
    type DivOutput;
    type ModOutput;
    fn div_mod(self, other: T) -> (Self::DivOutput, Self::ModOutput);
}

pub trait Abs {
    type Output;
    fn abs(self) -> Self::Output;
}

pub trait AbsDiff<T> {
    type Output;
    fn abs_diff(self, rhs: T) -> Self::Output;
}

pub trait Fraction {
    type NumeratorOutput;
    type DenominatorOutput;
    fn numerator(self) -> Self::NumeratorOutput;
    fn denominator(self) -> Self::DenominatorOutput;
}

pub trait Floor {
    type Output;
    fn floor(self) -> Self::Output;
}

pub trait Ceil {
    type Output;
    fn ceil(self) -> Self::Output;
}
