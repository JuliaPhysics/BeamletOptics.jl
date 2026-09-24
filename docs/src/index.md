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
    link: /basics/visualization/overview
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

A circularly polarized Gaussian laser beam through a Keplerian beam expander -- beam radius, focus and polarization are traced along with the rays.

```@example quickstart
using GLMakie, BeamletOptics

const mm = 1e-3

# Keplerian beam expander: f₁ = 15 mm, f₂ = 45 mm → 3× magnification
lens1 = ThinLens(15mm, 15mm, 12mm, 1.5)      # biconvex, f = R for n = 1.5
lens2 = ThinLens(45mm, 45mm, 25mm, 1.5)
translate3d!(lens1, [0, 20mm, 0])
translate3d!(lens2, [0, 80mm, 0])              # spacing f₁ + f₂
system = System([lens1, lens2])

# 532 nm laser, 1.5 mm waist radius, circularly polarized
laser = AstigmaticGaussianBeamlet([0, 0, 0], [0, 1, 0], 532e-9, 1.5mm; E0=[1, 0, im]/√2)
solve_system!(system, laser)

# plot lenses and beam
fig = Figure(size=(800, 300))
ax = LScene(fig[1, 1], show_axis=false)
render!(ax, system)
render!(ax, laser; color=RGBAf(0.1, 0.8, 0.1, 0.25), flen=0.04,
        show_polarization=true, pol_λ=4mm, pol_gain_max=4, pol_scale=1.5)
render_lcs!(ax, [-5mm, 50mm, -25mm]; scale=3, show_labels=true)
set_orthographic(ax)

cview = [                                               # hide
 -0.558947    0.829203   -5.72459e-17  -0.0528982       # hide
 -0.0248834  -0.0167734   0.99955       0.00802568      # hide
  0.82883     0.558696    0.0300089    -2.26736         # hide
  0.0         0.0         0.0           1.0             # hide     
]                                                       # hide
set_view(ax, cview) # hide
save("quickstart.png", fig, px_per_unit=3, update=false) # hide

# beam radius behind the expander relative to the input beam
gauss_parameters(laser, 0.12)[1] / gauss_parameters(laser, 0.0)[1]
```

![Circularly polarized beam in a Keplerian beam expander](quickstart.png)

## Resources to get you started

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
