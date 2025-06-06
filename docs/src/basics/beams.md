```@setup beams
beam_showcase_dir = joinpath(@__DIR__, "..", "assets", "beam_renders")

include(joinpath(beam_showcase_dir, "gb_showcase.jl"))
```

# Beams

As mentioned in the [Rays](@ref) section, a beam within the context of this package serves as a data structure for storing collections of rays, forming the backbone of the simulation framework. Beams are intended to be designed as [AbstractTrees](https://github.com/JuliaCollections/AbstractTrees.jl) to allow for ray bifurcations, e.g. in the case of optical elements such as beamsplitters. The [`solve_system!`](@ref) function relies on this data structure to perform ray tracing computations within optical systems. 

To ensure compatibility and extensibility, beam types must adhere to the [`BeamletOptics.AbstractBeam`](@ref) interface. Refer to its documentation for more information.

## Basic beam

A minimal implementation of the [`BeamletOptics.AbstractBeam`](@ref) type is provided by the [`Beam`](@ref). It can be used to store a light path through an optical system. If the beam is split, its children will be recursively traced until all paths are solved.

```@docs; canonical=false
Beam
```

The propagation of multiple parallel beams through an imaging system is illustrated below. The [Beam expander](@ref) tutorial covers the use of the [`Beam`](@ref) in more detail. A ray tracing example using [`Beam`](@ref)s is shown below.

```@example telescope_with_beams
using CairoMakie, BeamletOptics

# setup system
LA1353 = SphericalLens(103e-3, Inf, 10.1e-3, 75e-3, λ->1.5007)
LA1131 = SphericalLens(25.8e-3, Inf, 5.3e-3, BeamletOptics.inch, λ->1.5007)
LA1805 = SphericalLens(15.5e-3, Inf, 8.6e-3, BeamletOptics.inch, λ->1.5007)

zrotate3d!(LA1353, deg2rad(180))
zrotate3d!(LA1805, deg2rad(180))

translate3d!(LA1131, [0, 0.25, 0])
translate3d!(LA1805, [0, 0.29, 0])

system = System([LA1353, LA1131, LA1805])

# create render
fig = Figure(size=(600,220))
rend = Axis3(fig[1,1], aspect=(1,5,1), limits=(-0.05, 0.05, -0.1, 0.4, -0.05, 0.05), azimuth=0, elevation=1e-3)
hidexdecorations!(rend)
hidezdecorations!(rend)

render!(rend, system)

zs = LinRange(-0.025, 0.025, 7)

for z in zs
    beam = Beam(Ray([0, -0.1, z], [0, 1, 0.0], 1550e-9))
    solve_system!(system, beam)
    render!(rend, beam, flen=0.1)
end

save("telescope_with_beams.png", fig, px_per_unit=4); nothing # hide
```

![Telescope with beams](telescope_with_beams.png)

## Gaussian beamlet

Lasers are common devices in modern optical laboratories. Modeling their propagation through an optical setup can be of interest when planning new experiments. Geometrical ray tracing struggles to capture the propagation of a laser beam correctly, since it can not inherently capture the wave nature of e.g. the [Gaussian beam](https://www.rp-photonics.com/gaussian_beams.html).

The electric field of the ``\text{TEM}_{00}`` spatial Gaussian mode can be calculated analytically using the [`BeamletOptics.electric_field`](@ref) function: 

```@docs; canonical=false
BeamletOptics.electric_field(r::Real, z::Real, E0, w0, w, k, ψ, R)
```

The evolution of this field through an optical system can be modeled e.g. by the ray transfer matrix formalism using the complex ``q``-factor [Saleh2019; pp. 27](@cite). A Julia-based implementation of this approach can be found in [ABCDMatrixOptics.jl](https://github.com/JuliaPhysics/ABCDMatrixOptics.jl). However, in the case of this package another approach will be used.

### Complex ray tracing

In 1968 an internal publication at Bell Labs by J. Arnaud introduced the concept of complex rays wherein three geometrical beams can be used to model the propagation of a Gaussian in fundamental mode through a symmetric optical system, i.e. without the Gaussian obtaining astigmatism and/or higher-order abberations. This method is analoguos to the ray transfer matrix based ``q``-method [Arnaud1968](@cite).

Without extensions of the original method, the following key assumptions must be met such that this method can be applied

- all (complex) beams of the Gaussian in question must intersect the same optical elements
- the optical elements are large compared to the beam (waist)
- the paraxial approximation must hold for each beam
- the Gaussian may not be clipped by hard apertures
- Lagrange invariant must be fulfilled

Various versions of this approach have been implemented under different names in commercial software, most notably [FRED](https://photonengr.com/fred-software/) and [Code V](https://www.synopsys.com/optical-solutions/codev.html), as well as in open source software, e.g. 

- [Raypier](https://github.com/bryancole/raypier_optics) - based on Cython, maintenance status not known
- [Poke](https://github.com/Jashcraf/poke) - based on Zemax API and Python, maintained by J. Ashcraft et al. [Ashcraft:2022](@cite)
- [IfoCAD](https://www.aei.mpg.de/ifocad-de) - maintenance status not known, refer to Wanner et al. [Wanner:2017](@cite)

This package implements the above method via the [`GaussianBeamlet`](@ref) and the `BeamletOptics.AstigmaticGaussianBeamlet` (**Work in progress**).

### Stigmatic Beamlets

The [`GaussianBeamlet`](@ref) implements the [`BeamletOptics.AbstractBeam`](@ref) interface and can be used to model the propagation of a monochromatic Gaussian (``\text{TEM}_{00}``-mode) through optical system where all optics lie on the optical axis, e.g. no tip and/or tilt dealignment, and abberations can be neglected. It is represented by a `chief` (red), `waist` (blue) and `divergence` (green) beam. See below how these beans are placed in relation to the envelope of the Gaussian beam.

![Complex ray tracing I](gbtest1.png)

```@docs; canonical=false
GaussianBeamlet
```

A [`GaussianBeamlet`](@ref) can be constructed via:

```@docs; canonical=false
GaussianBeamlet(::AbstractArray{<:Real}, ::AbstractArray{<:Real}, ::Real, ::Real)
```

#### Obtaining the beam parameters 

Once a `GaussianBeamlet` has been traced through an optical system, several parameters might be of interest for further analysis. In order to relate the traced geometrical beams/rays to the Gaussian parameters, the publications of Arnaud, Herloski et al. and DeJager et al. are used [Arnaud1968, Herloski1983, DeJager1992](@cite). Consider the following system where a Gaussian beam with arbitrary parameters has been traced through a lens using the approach outlined in the [Complex ray tracing](@ref) section.

![Complex ray tracing II](gbtest2.png)

The user can obtain parameters such as the beam waist radius, the radius of curvature and more using the [`BeamletOptics.gauss_parameters`](@ref) function. Below the local waist radius and curvature ``R = r^{-1}`` have been calculated for the example above.

![Gauss beam parameters](gauss_parameters.png)

### Astigmatic Polarized Beamlets

- Assumptions
    - homogeneous polarization distribution across waist
    - Lagrange invariant

!!! info
    WORK IN PROGRESS