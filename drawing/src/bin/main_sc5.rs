use algebraeon_drawing::canvas::Canvas;
use algebraeon_drawing::canvas2d::Canvas2D;
use algebraeon_drawing::canvas2d::MouseWheelZoomCamera;
use algebraeon_geometry::parse::parse_shape;

fn main() {
    // let shape = parse_shape(
    //     "ConvexHull((0, 0), (1, 0), (0, 1), (1, 1)) \\ ((1/2, 1/2) | (1/3, 1/3) | (1/4, 1/4) | (1/5, 1/5) | (1/6, 1/6) | (1/3, 4/5))",
    // );

    // let shape = parse_shape(
    //     "ConvexHull((0, 0), (3, 0), (0, 3), (3, 3)) \\ ConvexHull((1, 1), (2, 1), (1, 2), (2, 2))",
    // );

    // let shape = parse_shape("Polygon((0, 0), (0, 3), (1/2, -1), (3, 3), (3, 0), (5/2, 4))");

    // let shape = parse_shape("ConvexHull((0, 0), (1, 1/10)) + Loop((0, 0), (10, 0), (10, 10), (0, 10))");

    // let shape = parse_shape(
    //     "((Polygon((0, 0), (4, 0), (4, 1), (3, 1), (3, -1), (2, -1), (2, 2), (1, 2)) | Polygon((5/2, 0/2), (4/2, 1/2), (3/2, 0/2), (4/2, -1/2))) \\ (Loop((0, 0), (1, 1), (11/4, -1/4)) \\ (1/3, 1/3))) + Lines((0, 0), (1/10, 1/10)) + Lines((0, 0), (1/10, -1/10))",
    // );

    let shape = parse_shape(
        "Polygon((0, 0), (6, 0), (6, 6), (2, 3), (0, 4)) \\ Polygon((1, 1), (2, 1), (2, 2), (1, 3))",
    );

    let mut canvas = Canvas2D::new(Box::new(MouseWheelZoomCamera::new()));
    canvas.plot(shape);
    canvas.run();
}
