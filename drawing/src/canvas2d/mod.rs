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

        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::default(),
            compatible_surface: Some(&surface),
            force_fallback_adapter: false,
        }))
        .unwrap();

        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
            required_features: wgpu::Features::empty(),
            // WebGL doesn't support all of wgpu's features, so if
            // we're building for the web, we'll have to disable some.
            required_limits: if cfg!(target_arch = "wasm32") {
                wgpu::Limits::downlevel_webgl2_defaults()
            } else {
                wgpu::Limits::default()
            },
            label: None,
            memory_hints: Default::default(),
            trace: wgpu::Trace::Off,
        }))
        .unwrap();

        let surface_caps = surface.get_capabilities(&adapter);
        // Shader code in this tutorial assumes an sRGB surface texture. Using a different
        // one will result in all the colors coming out darker. If you want to support non
        // sRGB surfaces, you'll need to account for that when drawing to the frame.
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
            present_mode: surface_caps.present_modes[0],
            alpha_mode: surface_caps.alpha_modes[0],
            view_formats: vec![],
            desired_maximum_frame_latency: 2,
        };

        Self {
            window,
            surface,
            device,
            queue,
            config,
            size,
        }
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
        static mut count: usize = 0;
        match event {
            WindowEvent::CloseRequested => {
                event_loop.exit();
            }
            WindowEvent::RedrawRequested => {
                #[allow(static_mut_refs)]
                unsafe {
                    count += 1;

                    println!("frame# = {:?}", count);
                }
                // Redraw the application.
                //
                // It's preferable for applications that do not render continuously to render in
                // this event rather than in AboutToWait, since rendering in here allows
                // the program to gracefully handle redraws requested by the OS.

                // Draw.

                // Queue a RedrawRequested event.
                //
                // You only need to call this if you've determined that you need to redraw in
                // applications which do not always need to. Applications that redraw continuously
                // can render here instead.
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
