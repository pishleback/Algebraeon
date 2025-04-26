use std::{
    cell::OnceCell,
    sync::Arc,
    time::{Duration, Instant},
};
use winit::{
    application::ApplicationHandler,
    event::{Event, WindowEvent},
    event_loop::{ActiveEventLoop, ControlFlow, EventLoop},
    window::{Window, WindowId},
};

pub trait EventHandler {
    fn tick(&mut self, dt: &Duration);
    fn event(&mut self, ev: &Event<()>);
}

pub trait Camera: EventHandler {
    fn view_matrix_and_shift(&self, display_size: (u32, u32)) -> ([[f64; 2]; 2], [f64; 2]);
}

// pub trait DrawElement: EventHandler {
//     fn draw(
//         &mut self,
//         display: &Display<WindowSurface>,
//         target: &mut Frame,
//         display_size: (u32, u32),
//         camera: &Box<dyn Camera>,
//     );
// }

pub struct MouseWheelZoomCamera {
    // the x coordinate at the centre of the screen
    mid_x: f64,
    // the y coordinate at the centre of the screen
    mid_y: f64,
    // the square root of the visible area
    sqrt_area: f64,
}

impl MouseWheelZoomCamera {
    pub fn new() -> Self {
        Self {
            mid_x: 1.0,
            mid_y: 1.0,
            sqrt_area: 6.0,
        }
    }
}

impl EventHandler for MouseWheelZoomCamera {
    fn tick(&mut self, dt: &Duration) {}
    fn event(&mut self, ev: &Event<()>) {}
}

impl Camera for MouseWheelZoomCamera {
    fn view_matrix_and_shift(&self, display_size: (u32, u32)) -> ([[f64; 2]; 2], [f64; 2]) {
        let display_size = (display_size.0 as f64, display_size.1 as f64);
        let avg_side = (display_size.0 * display_size.1).sqrt();
        let x_mult = 2.0 * avg_side / (display_size.0 * self.sqrt_area);
        let y_mult = 2.0 * avg_side / (display_size.1 * self.sqrt_area);
        (
            [[x_mult, 0.0], [0.0, y_mult]],
            [-self.mid_x * x_mult, -self.mid_y * y_mult],
        )
    }
}

pub struct Canvas2D {
    camera: Box<dyn Camera>,
    // elements: Vec<Box<dyn DrawElement>>,
}

struct WindowState {
    window: Arc<Window>,
    surface: wgpu::Surface<'static>,
    device: wgpu::Device,
    queue: wgpu::Queue,
    config: wgpu::SurfaceConfiguration,
    size: winit::dpi::PhysicalSize<u32>,
    // rendering
    render_pipeline: wgpu::RenderPipeline,
}

impl WindowState {
    fn new(window: Arc<Window>) -> Self {
        let size = window.inner_size();

        // The instance is a handle to our GPU
        // Backends::all => Vulkan + Metal + DX12 + Browser WebGPU
        let instance = wgpu::Instance::new(&wgpu::InstanceDescriptor {
            backends: wgpu::Backends::PRIMARY,
            ..Default::default()
        });

        let surface = instance.create_surface(window.clone()).unwrap();

        // Use instance.enumerate_adapters(wgpu::Backends::all()) for a full list of adapters
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::default(),
            compatible_surface: Some(&surface),
            force_fallback_adapter: false,
        }))
        .unwrap();

        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
            required_features: wgpu::Features::empty(),
            required_limits: wgpu::Limits::default(),
            label: None,
            memory_hints: Default::default(),
            trace: wgpu::Trace::Off,
        }))
        .unwrap();

        let surface_caps = surface.get_capabilities(&adapter);
        let surface_format = surface_caps
            .formats
            .iter()
            .find(|f| f.is_srgb())
            .copied()
            .unwrap_or(surface_caps.formats[0]);
        let config = wgpu::SurfaceConfiguration {
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            format: surface_format,
            width: size.width,
            height: size.height,
            present_mode: wgpu::PresentMode::Fifo,
            // present_mode: surface_caps.present_modes[0],
            alpha_mode: surface_caps.alpha_modes[0],
            view_formats: vec![],
            desired_maximum_frame_latency: 2,
        };

        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("Shader"),
            source: wgpu::ShaderSource::Wgsl(include_str!("shader.wgsl").into()),
        });

        let render_pipeline_layout =
            device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                label: Some("Render Pipeline Layout"),
                bind_group_layouts: &[],
                push_constant_ranges: &[],
            });

        let render_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("Render Pipeline"),
            layout: Some(&render_pipeline_layout),
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: Some("vs_main"), // 1.
                buffers: &[],                 // 2.
                compilation_options: wgpu::PipelineCompilationOptions::default(),
            },
            fragment: Some(wgpu::FragmentState {
                // 3.
                module: &shader,
                entry_point: Some("fs_main"),
                targets: &[Some(wgpu::ColorTargetState {
                    // 4.
                    format: config.format,
                    blend: Some(wgpu::BlendState::REPLACE),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
                compilation_options: wgpu::PipelineCompilationOptions::default(),
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList, // 1.
                strip_index_format: None,
                front_face: wgpu::FrontFace::Ccw, // 2.
                cull_mode: Some(wgpu::Face::Back),
                // Setting this to anything other than Fill requires Features::NON_FILL_POLYGON_MODE
                polygon_mode: wgpu::PolygonMode::Fill,
                // Requires Features::DEPTH_CLIP_CONTROL
                unclipped_depth: false,
                // Requires Features::CONSERVATIVE_RASTERIZATION
                conservative: false,
            },
            depth_stencil: None, // 1.
            multisample: wgpu::MultisampleState {
                count: 1,                         // 2.
                mask: !0,                         // 3.
                alpha_to_coverage_enabled: false, // 4.
            },
            multiview: None, // 5.
            cache: None,     // 6.
        });

        Self {
            window,
            surface,
            device,
            queue,
            config,
            size,
            render_pipeline,
        }
    }

    pub fn resize(&mut self, new_size: winit::dpi::PhysicalSize<u32>) {
        if new_size.width > 0 && new_size.height > 0 {
            self.size = new_size;
            self.config.width = new_size.width;
            self.config.height = new_size.height;
            self.surface.configure(&self.device, &self.config);
        }
    }

    pub fn reconfigure(&mut self) {
        self.resize(self.size);
    }

    fn render(&mut self) -> Result<(), wgpu::SurfaceError> {
        let output = self.surface.get_current_texture()?;
        let view = output
            .texture
            .create_view(&wgpu::TextureViewDescriptor::default());
        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("Render Encoder"),
            });

        {
            let mut render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("Render Pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color {
                            r: 0.1,
                            g: 0.2,
                            b: 0.3,
                            a: 1.0,
                        }),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: None,
                occlusion_query_set: None,
                timestamp_writes: None,
            });

            render_pass.set_pipeline(&self.render_pipeline);
            render_pass.draw(0..3, 0..1);
        }

        // submit will accept anything that implements IntoIter
        self.queue.submit(std::iter::once(encoder.finish()));
        output.present();

        Ok(())
    }
}

struct Canvas2DApp {
    canvas: Canvas2D,
    window: OnceCell<Arc<Window>>,
    window_state: Option<WindowState>,
}

impl ApplicationHandler for Canvas2DApp {
    fn resumed(&mut self, event_loop: &ActiveEventLoop) {
        let window = self
            .window
            .get_or_init(|| {
                Arc::new(
                    event_loop
                        .create_window(Window::default_attributes())
                        .unwrap(),
                )
            })
            .clone();
        self.window_state = Some(WindowState::new(window));
    }

    fn window_event(&mut self, event_loop: &ActiveEventLoop, id: WindowId, event: WindowEvent) {
        match event {
            WindowEvent::CloseRequested => {
                event_loop.exit();
            }
            WindowEvent::Resized(physical_size) => {
                self.window_state.as_mut().unwrap().resize(physical_size);
            }
            WindowEvent::RedrawRequested => {
                match self.window_state.as_mut().unwrap().render() {
                    Ok(()) => {}
                    // Reconfigure the surface if it's lost or outdated
                    Err(wgpu::SurfaceError::Lost | wgpu::SurfaceError::Outdated) => {
                        self.window_state.as_mut().unwrap().reconfigure();
                    }
                    // The system is out of memory, we should probably quit
                    Err(wgpu::SurfaceError::OutOfMemory | wgpu::SurfaceError::Other) => {
                        log::error!("OutOfMemory");
                        event_loop.exit();
                    }

                    // This happens when the a frame takes too long to present
                    Err(wgpu::SurfaceError::Timeout) => {
                        log::warn!("Surface timeout")
                    }
                }

                self.window.get().unwrap().request_redraw();
            }
            _ => (),
        }
    }
}

impl Canvas2D {
    pub fn new(camera: Box<dyn Camera>) -> Self {
        Self {
            camera,
            // elements: vec![],
        }
    }

    // pub fn add_element(&mut self, element: Box<dyn DrawElement>) {
    //     self.elements.push(element);
    // }

    fn tick(&mut self, dt: &Duration) {
        // for element in &mut self.elements {
        //     element.tick(dt);
        // }
    }

    // fn draw(&mut self, display: &Display<WindowSurface>, display_size: (u32, u32)) {
    //     let mut target = display.draw();
    //     target.clear_color(0.0, 0.0, 0.0, 1.0);
    //     for element in &mut self.elements {
    //         element.draw(display, &mut target, display_size, &self.camera);
    //     }
    //     target.finish().unwrap();
    // }

    fn event(&mut self, ev: &Event<()>) {
        // for element in &mut self.elements {
        //     element.event(ev);
        // }
    }

    pub fn run(mut self, display_width: u32, display_height: u32) {
        let event_loop = EventLoop::new().unwrap();
        event_loop.set_control_flow(ControlFlow::Poll);
        let mut app = Canvas2DApp {
            canvas: self,
            window: OnceCell::default(),
            window_state: None,
        };
        event_loop.run_app(&mut app).unwrap();
    }
}
