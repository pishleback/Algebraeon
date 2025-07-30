// Vertex shader

struct VertexInput {
    @location(0) offset: vec2<f32>,
};

struct InstanceInput {
    @location(1) pos: vec2<f32>,
    @location(2) radius: f32,
};

struct VertexOutput {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) color: vec3<f32>,
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
    out.color = vec3<f32>(0.0, 0.0, 0.0);
    var scale = sqrt(abs(camera.matrix[0][0] * camera.matrix[1][1] - camera.matrix[0][1] * camera.matrix[1][0]));
    out.clip_position = vec4<f32>(camera.matrix * (instance.radius * vertex.offset / scale + instance.pos) + camera.shift, 0.0, 1.0);
    return out;
}

// Fragment shader
@fragment
fn fs_main(in: VertexOutput) -> @location(0) vec4<f32> {
    return vec4<f32>(in.color, 1.0);
}

 

 