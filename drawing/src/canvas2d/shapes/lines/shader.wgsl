// Vertex shader

struct VertexInput {
    @location(0) offset: vec2<f32>,
    @location(1) which_end: f32
};

struct InstanceInput {
    @location(2) pos1: vec2<f32>,
    @location(3) pos1_screen_offset: vec2<f32>,
    @location(4) pos2: vec2<f32>,
    @location(5) pos2_screen_offset: vec2<f32>,
    @location(6) radius: f32,
    @location(7) colour: vec3<f32>,
};

struct VertexOutput {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) colour: vec3<f32>,
};

struct CameraUniform {
    matrix: mat2x2<f32>,
    matrix_inv: mat2x2<f32>,
    shift: vec2<f32>,
};
@group(0) @binding(0)
var<uniform> camera: CameraUniform;

fn rotate2D(angle: f32, v: vec2<f32>) -> vec2<f32> {
    let cosA = cos(angle);
    let sinA = sin(angle);
    return vec2<f32>(
        v.x * cosA - v.y * sinA,
        v.x * sinA + v.y * cosA
    );
}

@vertex
fn vs_main(
    vertex: VertexInput,
    instance: InstanceInput,
) -> VertexOutput {
    var out: VertexOutput;
    out.colour = instance.colour;
    var scale = sqrt(abs(camera.matrix[0][0] * camera.matrix[1][1] - camera.matrix[0][1] * camera.matrix[1][0]));

    var total_pos1 = (instance.pos1_screen_offset) / scale + instance.pos1;
    var total_pos2 = (instance.pos2_screen_offset) / scale + instance.pos2;

    var vec = total_pos2 - total_pos1;
    var angle = atan2(vec.y, vec.x);

    if vertex.which_end < 0.5 {
        // == 0.0
        out.clip_position = vec4<f32>(camera.matrix * (instance.radius * rotate2D(angle, 0.5 * vertex.offset) / scale + total_pos1) + camera.shift, 0.0, 1.0);
    } else {
        // == 1.0
        out.clip_position = vec4<f32>(camera.matrix * (instance.radius * rotate2D(angle, 0.5 * vertex.offset) / scale + total_pos2) + camera.shift, 0.0, 1.0);
    }

    return out;
}

// Fragment shader
@fragment
fn fs_main(in: VertexOutput) -> @location(0) vec4<f32> {
    return vec4<f32>(in.colour, 1.0);
}

 

 