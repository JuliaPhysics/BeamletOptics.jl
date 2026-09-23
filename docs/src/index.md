```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: A digital optics laboratory
  tagline: 3D ray tracing and Gaussian beamlet propagation for optical setups in Julia
  image:
    light: /logo.svg
    dark: /logo-dark.svg
    alt: BeamletOptics
  actions:
    - theme: brand
      text: Get started
      link: /tutorials/
    - theme: alt
      text: Basics
      link: /basics/intro
    - theme: alt
      text: View on Github
      link: https://github.com/JuliaPhysics/BeamletOptics.jl

features:
  - icon: "🔦"
    title: 3D ray tracing
    details: Hybrid sequential and non-sequential ray tracing without paraxial approximation.
    link: /basics/rays
    linkText: Learn more
  - icon: "🌊"
    title: Gaussian beamlets
    details: TEM₀₀ and astigmatic Gaussian beamlets, including coherent diffraction.
    link: /basics/beams/stigmatic_beam
    linkText: Learn more
  - icon: "🔭"
    title: Optical components
    details: Mirrors, lenses, beamsplitters, detectors and polarizing optics.
    link: /basics/components/components
    linkText: Learn more
  - icon: "🎛️"
    title: Kinematic API
    details: Translate, rotate and group elements to model moving or vibrating setups.
    link: /basics/components/components#Moving-optical-elements
    linkText: Learn more
  - icon: "📊"
    title: Makie visualization
    details: Render setups and beams in 2D and 3D with CairoMakie or GLMakie.
    link: /basics/render
    linkText: Learn more
  - icon: "🧩"
    title: Extensible
    details: Implement your own optical interactions via the API.
    link: /api/api
    linkText: Learn more
---
```

```@raw html
<p style="margin-bottom:2cm"></p>

<div class="vp-doc" style="width:80%; margin:auto">
```

Building optical setups in a laboratory environment -- for instance a laser interferometer -- is a common task for optical engineers and physicists. This package is intended to provide a simulation environment in which the user can quickly analyze and layout simple optical components like lenses or beamsplitters before committing to a breadboard setup.  

## What is the purpose of this package

This package mainly tries to provide a simple Gaussian beamlet propagation tool for coherent, monochromatic and directed light sources. It also offers a convenient kinematic API that allows for the easy placing of optical elements and straight-forward simulation of moving or vibrating components. 

For this purpose, the package implements a traditional ray tracing solver. This forms the backbone of the Gaussian beamlet tracing scheme that has been implemented to model the propagation of laser beams.

!!! info "What this package is not"
    This package does not include tools for optimizing optical systems, such as fine-tuning lens surfaces to minimize specific aberrations in multi-lens setups. Instead, the package is designed as a digital laboratory where you can play around with stuff before buying it.

## BMO in 30 seconds

White light through a dense flint prism -- dispersion comes for free with a Sellmeier glass model.

```@example quickstart
using GLMakie, BeamletOptics

# dense flint glass N-SF11, Sellmeier coefficients in µm²
SF11 = SellmeierEquation(1.73759695, 0.313747346, 1.89878101, 0.013188707, 0.0623068142, 155.23629)
prism = RightAnglePrism(40e-3, 20e-3, SF11)
system = System([prism])

fig = Figure(size=(800, 400))
ax = Axis3(fig[1,1], aspect=:data, azimuth=-π/2, elevation=π/2, limits=(-0.1, 0.2, -0.06, 0.06, -0.02, 0.02))
hidedecorations!(ax); hidespines!(ax)
render!(ax, system)

θ = deg2rad(45)                                     # angle of incidence
for (λ, c) in zip(420e-9:40e-9:660e-9, cgrad(:rainbow, 7, categorical=true))
    beam = Beam([-0.1, -0.088, 0], [cos(θ), sin(θ), 0], λ)
    solve_system!(system, beam)
    render!(ax, beam, color=c, flen=0.25)
    render!(ax, first(rays(beam)), color=:black)    # incoming white light
end
save("quickstart.png", fig, px_per_unit=3); nothing # hide
```

![Prism dispersion](quickstart.png)

## Showcase

```@raw html
<div class="bmo-gallery">
<div class="bmo-tile">
```

![Michelson interferometer](tutorials/mi_intro_fig.png)

```@raw html
<p><Badge type="warning" text="Intermediate" /></p>
```

[Michelson interferometer](@ref)

```@raw html
<p class="bmo-teaser">Simulate fringes of a Michelson interferometer with moving mirrors.</p>
</div>
<div class="bmo-tile">
```

![Raman spectroscopy](tutorials/or_intro_fig.png)

```@raw html
<p><Badge type="danger" text="Advanced" /></p>
```

[Raman spectroscopy](@ref)

```@raw html
<p class="bmo-teaser">Model the OpenRAMAN spectrometer with custom components.</p>
</div>
<div class="bmo-tile">
```

![Miniature microscope](tutorials/ucla_intro_fig.png)

```@raw html
<p><Badge type="warning" text="Intermediate" /></p>
```

[Miniature microscope](@ref)

```@raw html
<p class="bmo-teaser">Rebuild the imaging path of the UCLA 2P miniscope.</p>
</div>
</div>
```

[Browse all tutorials and examples](@ref "Tutorials and examples")

## Installation

!!! warning
    This package requires Julia ≥ 1.12

You can add this package to your project by entering the package manager (press `]` in the REPL) and typing `add BeamletOptics`. It is also recommended that you `add GLMakie`. You can include this package into your current scope via `using BeamletOptics`. If a Makie version is loaded before or after the inclusion of this package, the extension provided as part of this package will enable additional visualization functions. 

## Citation and license

The BeamletOptics package is made available under the MIT license. If you use this package for your research, we encourage you to cite it. For your convenience, a BibTeX entry is provided as part of the package (CITATION.bib) or on [Zenodo](https://zenodo.org/records/15090784).

## Similar packages

A variety of packages and tools exist that implement similar approaches or offer optics modeling capabilities.

```@raw html
<details class="details custom-block">
<summary>Show similar packages</summary>
```

Within the Julia ecosystem, the following packages need to be mentioned:

- [OpticSim.jl](https://github.com/brianguenter/OpticSim.jl)
- [FluxOptics.jl](https://github.com/anscoil/FluxOptics.jl)
- [ABCDMatrixOptics.jl](https://github.com/JuliaPhysics/ABCDMatrixOptics.jl)
- [WaveOpticsPropagation.jl](https://github.com/JuliaPhysics/WaveOpticsPropagation.jl)

More broadly speaking, have a look at these packages as well:

- [DynamicalBilliards.jl](https://github.com/JuliaDynamics/DynamicalBilliards.jl)
- [RayTracer.jl](https://github.com/avik-pal/RayTracer.jl)

```@raw html
</details>
```

There also exists a plethora of commercial and non-commercial simulation frameworks outside of the Julia ecosystem. For specific examples regarding the beamlet method used in this package, refer to the [Complex ray tracing](@ref) section. 

```@raw html
</div>
```
