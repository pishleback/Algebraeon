// Vertex shader

struct VertexInput {
    @location(0) position: vec2<f32>,
};

struct VertexOutput {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) coord: vec2<f32>,
};

struct CameraUniform {
    matrix: mat2x2<f32>,
    matrix_inv: mat2x2<f32>,
    shift: vec2<f32>,
};
@group(0) @binding(0)
var<uniform> camera: CameraUniform;

@vertex
fn vs_main(
    model: VertexInput,
) -> VertexOutput {
    var out: VertexOutput;
    out.clip_position = vec4<f32>(model.position, 0.0, 1.0);
    out.coord = camera.matrix_inv * (model.position - camera.shift);
    return out;
}

// Complex addition
fn c_add(a: vec2<f32>, b: vec2<f32>) -> vec2<f32> {
    return vec2<f32>(a.x + b.x, a.y + b.y);
}

// Complex multiplication
fn c_mul(a: vec2<f32>, b: vec2<f32>) -> vec2<f32> {
    return vec2<f32>(
        a.x * b.x - a.y * b.y,
        a.x * b.y + a.y * b.x
    );
}

// Complex magnitude
fn c_abs(z: vec2<f32>) -> f32 {
    return length(z);
}

// Complex argument (angle)
fn c_arg(z: vec2<f32>) -> f32 {
    return atan2(z.y, z.x); // in radians
}

// Convert HSL to RGB (H in radians, L in [0,1], S fixed at 1.0)
fn hsl_to_rgb(h: f32, s: f32, l: f32) -> vec3<f32> {
    let c = (1.0 - abs(2.0 * l - 1.0)) * s;
    let h_ = h / (3.1415926 * 2.0); // Normalize radians to [0, 1]
    let x = c * (1.0 - abs((h_ * 6.0) % 2.0 - 1.0));
    var r = 0.0;
    var g = 0.0;
    var b = 0.0;

    if (h_ < 1.0/6.0) {
        r = c; g = x; b = 0.0;
    } else if (h_ < 2.0/6.0) {
        r = x; g = c; b = 0.0;
    } else if (h_ < 3.0/6.0) {
        r = 0.0; g = c; b = x;
    } else if (h_ < 4.0/6.0) {
        r = 0.0; g = x; b = c;
    } else if (h_ < 5.0/6.0) {
        r = x; g = 0.0; b = c;
    } else {
        r = c; g = 0.0; b = x;
    }

    let m = l - c * 0.5;
    return vec3<f32>(r + m, g + m, b + m);
}

@fragment
fn fs_main(in: VertexOutput) -> @location(0) vec4<f32> {
    let z = in.coord;
    let w = eval_cfn(z);

    let hue = c_arg(w) + 3.1415926;
    let mag = c_abs(w);
    let lightness = (0.5 * sin(1.44269504089 * log2(mag)) + 1.0) * 0.5;
    let saturation = 1.0;

    let rgb = hsl_to_rgb(hue, saturation, lightness);
    return vec4<f32>(rgb * rgb, 1.0);
}

// Evaluate a complex polynomial using Horner's method
fn eval_cfn(z: vec2<f32>) -> vec2<f32> {
    var coeffs = array<vec2<f32>, 4>(
        vec2<f32>(-1.0, 0.0),
        vec2<f32>( 0.0, 0.0),
        vec2<f32>( 0.0, 0.0),
        vec2<f32>( 1.0, 0.0)
    );

    var result = coeffs[3];
    for (var i = 2i; i >= 0i; i = i - 1i) {
        result = c_add(c_mul(result, z), coeffs[i]);
    }
    return result;
}

 

 