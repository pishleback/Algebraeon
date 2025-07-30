// Vertex shader

struct VertexInput {
    @location(0) offset: vec2<f32>,
};

struct InstanceInput {
    @location(1) a: f32,
    @location(2) b: f32,
    @location(3) c: f32,
    @location(4) d: f32,
    @location(5) thickness: f32,
    @location(6) colour: vec3<f32>,
};

struct VertexOutput {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) colour: vec3<f32>,
    @location(1) pos: vec2<f32>,
    @location(2) a: f32,
    @location(3) b: f32,
    @location(4) c: f32,
    @location(5) d: f32,
    @location(6) thickness: f32,
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
    var scale = sqrt(abs(camera.matrix[0][0] * camera.matrix[1][1] - camera.matrix[0][1] * camera.matrix[1][0]));

    var model_pos = vec2<f32>(
        instance.a + (instance.b - instance.a) * vertex.offset.x,
        instance.c + (instance.d - instance.c) * vertex.offset.y
    );

    out.clip_position = vec4<f32>(camera.matrix * model_pos + camera.shift, 0.0, 1.0);
    out.pos = model_pos;
    out.a = instance.a;
    out.b = instance.b;
    out.c = instance.c;
    out.d = instance.d;
    out.thickness = 0.5 * instance.thickness / scale;
    return out;
}

// Fragment shader
@fragment
fn fs_main(in: VertexOutput) -> @location(0) vec4<f32> {
    if in.pos.x < in.a - in.thickness {
        discard;
    }
    if in.pos.x > in.b + in.thickness {
        discard;
    }
    if in.pos.y < in.c - in.thickness {
        discard;
    }
    if in.pos.y > in.d + in.thickness {
        discard;
    }
    if abs(in.pos.x - in.a) > in.thickness && abs(in.pos.x - in.b) > in.thickness && abs(in.pos.y - in.c) > in.thickness && abs(in.pos.y - in.d) > in.thickness {
        discard;
    }
    if length(in.pos - vec2<f32>(in.a + in.thickness, in.c + in.thickness)) > 2 * in.thickness && in.pos.x < in.a + in.thickness && in.pos.y < in.c + in.thickness {
        discard;
    }
    if length(in.pos - vec2<f32>(in.b - in.thickness, in.c + in.thickness)) > 2 * in.thickness && in.pos.x > in.b - in.thickness && in.pos.y < in.c + in.thickness {
        discard;
    }
    if length(in.pos - vec2<f32>(in.a + in.thickness, in.d - in.thickness)) > 2 * in.thickness && in.pos.x < in.a + in.thickness && in.pos.y > in.d - in.thickness {
        discard;
    }
    if length(in.pos - vec2<f32>(in.b - in.thickness, in.d - in.thickness)) > 2 * in.thickness && in.pos.x > in.b - in.thickness && in.pos.y > in.d - in.thickness {
        discard;
    }
    return vec4<f32>(in.colour, 1.0);
}

 

 