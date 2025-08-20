use crate::partial_simplicial_complex::PartialSimplicialComplex;
use algebraeon_nzq::RationalCanonicalStructure;
use lalrpop_util::lalrpop_mod;

lalrpop_mod!(shape_parser, "/parse/shape_grammar.rs");
mod ast;

pub fn parse_shape(string: &str) -> PartialSimplicialComplex<'static, RationalCanonicalStructure> {
    shape_parser::ShapeParser::new()
        .parse(string)
        .unwrap()
        .to_partial_simplicial_complex_root()
}

#[cfg(test)]
mod tests {
    use crate::parse::ast::{Point, SignedValue};

    use super::*;

    #[test]
    fn test() {
        assert_eq!(
            shape_parser::PointParser::new()
                .parse("(6, +7, -8, + 9, - 10)")
                .unwrap(),
            Point {
                coordinates: vec![
                    SignedValue {
                        sign: ast::Sign::Positive,
                        value: String::from("6")
                    },
                    SignedValue {
                        sign: ast::Sign::Positive,
                        value: String::from("7")
                    },
                    SignedValue {
                        sign: ast::Sign::Negative,
                        value: String::from("8")
                    },
                    SignedValue {
                        sign: ast::Sign::Positive,
                        value: String::from("9")
                    },
                    SignedValue {
                        sign: ast::Sign::Negative,
                        value: String::from("10")
                    }
                ]
            }
        );

        assert!(
            shape_parser::ShapeParser::new()
                .parse("ConvexHull((),(),)")
                .is_ok()
        );
        assert!(
            shape_parser::ShapeParser::new()
                .parse("ConvexHull((),())")
                .is_ok()
        );
        assert!(
            shape_parser::ShapeParser::new()
                .parse("ConvexHull((),)")
                .is_ok()
        );
        assert!(
            shape_parser::ShapeParser::new()
                .parse("ConvexHull(())")
                .is_ok()
        );
        assert!(
            shape_parser::ShapeParser::new()
                .parse("ConvexHull()")
                .is_ok()
        );
        assert!(
            shape_parser::ShapeParser::new()
                .parse("ConvexHull() | ConvexHull()")
                .is_ok()
        );
        assert!(
            shape_parser::ShapeParser::new()
                .parse("ConvexHull() | ConvexHull() |)")
                .is_err()
        );
        assert!(
            shape_parser::ShapeParser::new()
                .parse("ConvexHull() | ConvexHull() | ConvexHull()")
                .is_ok()
        );
        assert!(
            shape_parser::ShapeParser::new()
                .parse("ConvexHull() & ConvexHull()")
                .is_ok()
        );
        assert!(
            shape_parser::ShapeParser::new()
                .parse("ConvexHull() & ConvexHull() &")
                .is_err()
        );
        assert!(
            shape_parser::ShapeParser::new()
                .parse("ConvexHull() & ConvexHull() & ConvexHull()")
                .is_ok()
        );
        assert!(
            shape_parser::ShapeParser::new()
                .parse("ConvexHull() | ConvexHull() & ConvexHull()")
                .is_err()
        );
        assert!(
            shape_parser::ShapeParser::new()
                .parse("ConvexHull() & ConvexHull() | ConvexHull()")
                .is_err()
        );
        assert!(
            shape_parser::ShapeParser::new()
                .parse("ConvexHull() \\ ConvexHull()")
                .is_ok()
        );
        assert!(
            shape_parser::ShapeParser::new()
                .parse("ConvexHull() \\ ConvexHull() \\ ConvexHull()")
                .is_err()
        );
        assert!(
            shape_parser::ShapeParser::new()
                .parse("ConvexHull() \\ ConvexHull() | ConvexHull()")
                .is_err()
        );
        assert!(
            shape_parser::ShapeParser::new()
                .parse("(ConvexHull() \\ ConvexHull()) | ConvexHull()")
                .is_ok()
        );
        assert!(
            shape_parser::ShapeParser::new()
                .parse("ConvexHull() \\ ConvexHull() & ConvexHull()")
                .is_err()
        );
        assert!(
            shape_parser::ShapeParser::new()
                .parse("(ConvexHull() \\ ConvexHull()) & ConvexHull()")
                .is_ok()
        );
    }
}
