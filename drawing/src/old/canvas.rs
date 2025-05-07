use std::{cell::OnceCell, sync::Arc};

use wgpu::{CommandEncoder, SurfaceTexture, TextureView};
use winit::{
    application::ApplicationHandler,
    error::EventLoopError,
    event::WindowEvent,
    event_loop::{ActiveEventLoop, ControlFlow, EventLoop},
    window::{Window, WindowId},
};

pub struct RenderTarget<'r> {
    pub output: &'r SurfaceTexture,
    pub view: &'r TextureView,
    pub encoder: &'r mut CommandEncoder,
}

pub trait Renderer<State: CanvasState> {
    fn init(&self, window_state: &State) -> Box<dyn RendererInstance>;
}

pub trait RendererInstance {
    fn render(&self, state: &WgpuState, target: RenderTarget);
}

pub struct WgpuState {
    surface: wgpu::Surface<'static>,
    pub device: wgpu::Device,
    pub queue: wgpu::Queue,
    pub config: wgpu::SurfaceConfiguration,
    size: winit::dpi::PhysicalSize<u32>,
}

pub struct BasicCanvasState {
    wgpu: WgpuState,
}

pub trait CanvasState {
    fn new(window: Arc<Window>) -> Self;

    fn wgpu_state(&self) -> &WgpuState;

    fn size(&self) -> winit::dpi::PhysicalSize<u32> {
        self.wgpu_state().size
    }

    fn resize(&mut self, new_size: winit::dpi::PhysicalSize<u32>);

    fn reconfigure(&mut self) {
        self.resize(self.size());
    }

    fn render(
        &mut self,
        renderers: &Vec<Box<dyn RendererInstance>>,
    ) -> Result<(), wgpu::SurfaceError>;
}

impl CanvasState for BasicCanvasState {
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

        Self {
            wgpu: WgpuState {
                surface,
                device,
                queue,
                config,
                size,
            },
        }
    }

    fn resize(&mut self, new_size: winit::dpi::PhysicalSize<u32>) {
        if new_size.width > 0 && new_size.height > 0 {
            self.wgpu.size = new_size;
            self.wgpu.config.width = new_size.width;
            self.wgpu.config.height = new_size.height;
            self.wgpu.surface.configure(&self.wgpu.device, &self.wgpu.config);
        }
    }

    fn render(
        &mut self,
        renderers: &Vec<Box<dyn RendererInstance>>,
    ) -> Result<(), wgpu::SurfaceError> {
        let output = self.wgpu.surface.get_current_texture()?;
        let view = output
            .texture
            .create_view(&wgpu::TextureViewDescriptor::default());
        let mut encoder = self.wgpu
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("Render Encoder"),
            });
        for renderer in renderers {
            renderer.render(
                &self.wgpu,
                RenderTarget {
                    output: &output,
                    view: &view,
                    encoder: &mut encoder,
                },
            );
        }
        self.wgpu.queue.submit(std::iter::once(encoder.finish()));
        output.present();
        Ok(())
    }
    
    fn wgpu_state(&self) -> &WgpuState {
        &self.wgpu
    }
}

pub trait Canvas: Sized {
    type State: CanvasState;

    fn run(&self) -> Result<(), EventLoopError> {
        let event_loop = EventLoop::new().unwrap();
        event_loop.set_control_flow(ControlFlow::Poll);
        let mut app_handler = CanvasApplicationHandler::new(self);
        event_loop.run_app(&mut app_handler)
    }

    fn renderers(&self) -> impl Iterator<Item = &dyn Renderer<Self::State>>;
}

struct CanvasApplicationHandler<'c, C: Canvas> {
    app: &'c C,
    window: OnceCell<Arc<Window>>,
    window_state: Option<C::State>,
    renderers: Vec<Box<dyn RendererInstance>>,
}

impl<'c, C: Canvas> CanvasApplicationHandler<'c, C> {
    pub fn new(app: &'c C) -> Self {
        Self {
            app,
            window: OnceCell::default(),
            window_state: None,
            renderers: vec![],
        }
    }
}

impl<'c, C: Canvas> ApplicationHandler for CanvasApplicationHandler<'c, C> {
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
        self.window_state = Some(C::State::new(window));
        self.renderers = self
            .app
            .renderers()
            .map(|renderer| renderer.init(self.window_state.as_ref().unwrap()))
            .collect();
    }

    fn window_event(&mut self, event_loop: &ActiveEventLoop, _id: WindowId, event: WindowEvent) {
        match event {
            WindowEvent::CloseRequested => {
                event_loop.exit();
            }
            WindowEvent::Resized(physical_size) => {
                self.window_state.as_mut().unwrap().resize(physical_size);
            }
            WindowEvent::RedrawRequested => {
                match self.window_state.as_mut().unwrap().render(&self.renderers) {
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
