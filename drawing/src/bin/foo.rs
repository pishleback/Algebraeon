use algebraeon_drawing::{
    canvas::Canvas,
    canvas2d::{Canvas2D, MouseWheelZoomCamera, shapes::*},
};
use algebraeon_nzq::Integer;
use algebraeon_rings::{
    polynomial::{Polynomial, PolynomialFromStr},
    structure::MetaRealToFloat,
};

fn main() {
    simplelog::CombinedLogger::init(vec![simplelog::TermLogger::new(
        simplelog::LevelFilter::Trace,
        simplelog::Config::default(),
        simplelog::TerminalMode::Mixed,
        simplelog::ColorChoice::Auto,
    )])
    .unwrap();
    log::set_max_level(log::LevelFilter::Warn);

    let p = Polynomial::<Integer>::from_str(
        "(x^7 - 3 * x^3 - 5 + x) * (x + 1) * (x^3 - 3 * x + 1)",
        "x",
    )
    .unwrap();
    // let p = Polynomial::<Integer>::from_str("x^2 - 2", "x").unwrap();

    let mut canvas = Canvas2D::new(Box::new(MouseWheelZoomCamera::new()));

    canvas.plot_complex_polynomial(p.clone());

    let mut roots = p.all_complex_roots();
    for _ in 0..6 {
        for root in &mut roots {
            root.refine();
        }
    }

    canvas.plot_shapes(roots.iter().flat_map(|root| match root.isolate() {
        algebraeon_rings::isolated_algebraic::ComplexIsolatingRegion::Rational(r) => {
            vec![Shape::Point {
                x: r.as_f32(),
                y: 0.0,
            }]
        }
        algebraeon_rings::isolated_algebraic::ComplexIsolatingRegion::RealInterval(a, b) => vec![
            Shape::Line {
                x1: a.as_f32(),
                y1: 0.0,
                x2: b.as_f32(),
                y2: 0.0,
            },
            Shape::LineRaw {
                x1: a.as_f32(),
                x1s: 0.0,
                y1: 0.0,
                y1s: -1.0,
                x2: a.as_f32(),
                x2s: 0.0,
                y2: 0.0,
                y2s: 1.0,
            },
            Shape::LineRaw {
                x1: b.as_f32(),
                x1s: 0.0,
                y1: 0.0,
                y1s: -1.0,
                x2: b.as_f32(),
                x2s: 0.0,
                y2: 0.0,
                y2s: 1.0,
            },
        ],
        algebraeon_rings::isolated_algebraic::ComplexIsolatingRegion::Box(a, b, c, d) => {
            vec![Shape::Rectangle {
                a: a.as_f32(),
                b: b.as_f32(),
                c: c.as_f32(),
                d: d.as_f32(),
            }]
        }
    }));

    canvas.run();
}
