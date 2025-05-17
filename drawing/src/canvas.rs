use std::{cell::OnceCell, sync::Arc};
use winit::{
    application::ApplicationHandler,
    dpi::PhysicalSize,
    event::WindowEvent,
    event_loop::{ActiveEventLoop, ControlFlow, EventLoop},
    window::{Window, WindowId},
};

pub struct WgpuState {
    pub window: Arc<Window>,
    pub surface: wgpu::Surface<'static>,
    pub device: wgpu::Device,
    pub queue: wgpu::Queue,
    pub config: wgpu::SurfaceConfiguration,
    pub size: winit::dpi::PhysicalSize<u32>,
}

impl WgpuState {
    pub fn new(window: Arc<Window>) -> Self {
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
            present_mode: wgpu::PresentMode::Immediate,
            // present_mode: surface_caps.present_modes[0],
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
}

pub trait Canvas: Sized {
    type WindowState;

    fn new_window_state(&self, window: Arc<Window>) -> Self::WindowState;

    fn run(self) {
        let event_loop = EventLoop::new().unwrap();
        event_loop.set_control_flow(ControlFlow::Poll);
        let mut app = CanvasApplicationHandler {
            canvas: self,
            window: OnceCell::default(),
            window_state: None,
        };
        event_loop.run_app(&mut app).unwrap();
    }

    fn window_event(
        &mut self,
        window_state: &mut Self::WindowState,
        event_loop: &ActiveEventLoop,
        id: WindowId,
        event: WindowEvent,
    );
}

struct CanvasApplicationHandler<C: Canvas> {
    canvas: C,
    window: OnceCell<Arc<Window>>,
    window_state: Option<C::WindowState>,
}

impl<C: Canvas> ApplicationHandler for CanvasApplicationHandler<C> {
    fn resumed(&mut self, event_loop: &ActiveEventLoop) {
        let window = self
            .window
            .get_or_init(|| {
                Arc::new(
                    event_loop
                        .create_window(
                            Window::default_attributes()
                                .with_inner_size(PhysicalSize::new(2000, 1600)),
                        )
                        .unwrap(),
                )
            })
            .clone();
        self.window_state = Some(self.canvas.new_window_state(window));
    }

    fn window_event(&mut self, event_loop: &ActiveEventLoop, id: WindowId, event: WindowEvent) {
        self.canvas
            .window_event(self.window_state.as_mut().unwrap(), event_loop, id, event);
    }
}
