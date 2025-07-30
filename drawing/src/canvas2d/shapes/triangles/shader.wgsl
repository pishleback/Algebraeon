// Vertex shader

struct VertexInput {
    @location(0) n: u32,
};

struct InstanceInput {
    @location(1) pos1: vec2<f32>,
    @location(2) pos2: vec2<f32>,
    @location(3) pos3: vec2<f32>,
    @location(4) colour : vec3<f32>,
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

@vertex
fn vs_main(
    vertex: VertexInput,
    instance: InstanceInput,
) -> VertexOutput {
    var out: VertexOutput;
    out.colour = instance.colour;

    if vertex.n == 0 {
        out.clip_position = vec4<f32>(camera.matrix * instance.pos1 + camera.shift, 0.0, 1.0);
    } else if vertex.n == 1 {
        out.clip_position = vec4<f32>(camera.matrix * instance.pos2 + camera.shift, 0.0, 1.0);
    } else {
        // n == 2
        out.clip_position = vec4<f32>(camera.matrix * instance.pos3 + camera.shift, 0.0, 1.0);
    }
    return out;
}

// Fragment shader
@fragment
fn fs_main(in: VertexOutput) -> @location(0) vec4<f32> {
    return vec4<f32>(in.colour, 1.0);
}

 

 