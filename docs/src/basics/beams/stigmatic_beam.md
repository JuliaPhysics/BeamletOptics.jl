```@setup stigmatic_beamlet
beam_showcase_dir = joinpath(@__DIR__, "..", "..", "assets", "beam_renders")

Main.DocUtils.conditional_include(joinpath(beam_showcase_dir, "gb_showcase.jl"))
```

# Gaussian beamlet

Lasers are common devices in modern optical laboratories. Modeling their propagation through an optical setup can be of interest when planning new experiments. Geometrical ray tracing struggles to capture the propagation of a laser beam correctly, since it can not inherently capture the wave nature of e.g. the [Gaussian beam](https://www.rp-photonics.com/gaussian_beams.html).

The electric field of the ``\text{TEM}_{00}`` spatial Gaussian mode can be calculated analytically using the [`BeamletOptics.electric_field`](@ref) function: 

```@docs; canonical=false
BeamletOptics.electric_field(r::Real, z::Real, E0, w0, w, k, ψ, R)
```

The evolution of this field through an optical system can be modeled e.g. by the ray transfer matrix formalism using the complex ``q``-factor [Saleh2019; pp. 27](@cite). A Julia-based implementation of this approach can be found in [ABCDMatrixOptics.jl](https://github.com/JuliaPhysics/ABCDMatrixOptics.jl). However, in the case of this package another approach will be used.

## Complex ray tracing

In 1968 an internal publication at Bell Labs by J. Arnaud introduced the concept of complex rays wherein three geometrical beams can be used to model the propagation of a Gaussian in fundamental mode through a symmetric optical system, i.e. without the Gaussian obtaining astigmatism and/or higher-order abberations. This method is analoguos to the ray transfer matrix based ``q``-method [Arnaud1968](@cite).

Without extensions of the original method, the following key assumptions must be met such that this method can be applied

- all (complex) beams of the Gaussian in question must intersect the same optical elements
- the optical elements are large compared to the beam (waist)
- the paraxial approximation must hold for each beam
- the Gaussian may not be clipped by hard apertures
- Lagrange invariant must be fulfilled

Various versions of this approach have been implemented under different names in commercial software, most notably [FRED](https://photonengr.com/fred-software/), [Code V](https://www.synopsys.com/optical-solutions/codev.html) and [QUADOA](https://www.quadoa.com), as well as in open source software, e.g. 

- [Raypier](https://github.com/bryancole/raypier_optics) - based on Cython, maintenance status not known
- [Poke](https://github.com/Jashcraf/poke) - based on Zemax API and Python, maintained by J. Ashcraft et al. [Ashcraft:2022](@cite)
- [IfoCAD](https://www.aei.mpg.de/ifocad-de) - maintenance status not known, refer to Wanner et al. [Wanner:2017](@cite)

This package implements the above method via the stigmatic [`GaussianBeamlet`](@ref) and the [`AstigmaticGaussianBeamlet`](@ref), the latter of which is presented in more detail in the [Astigmatic polarized beamlets](@ref) chapter.

## Stigmatic beamlets

The [`GaussianBeamlet`](@ref) implements the [`BeamletOptics.AbstractBeam`](@ref) interface and can be used to model the propagation of a monochromatic Gaussian (``\text{TEM}_{00}``-mode) through optical system where all optics lie on the optical axis, e.g. no tip and/or tilt dealignment, and abberations can be neglected. It is represented by a `chief` (red), `waist` (blue) and `divergence` (green) beam. See below how these beams are placed in relation to the envelope of the Gaussian beam.

![Complex ray tracing I](gbtest1.png)

```@docs; canonical=false
GaussianBeamlet
```

A [`GaussianBeamlet`](@ref) can be constructed via:

```@docs; canonical=false
GaussianBeamlet(::AbstractArray{<:Real}, ::AbstractArray{<:Real}, ::Real, ::Real)
```

### Obtaining the beam parameters 

Once a `GaussianBeamlet` has been traced through an optical system, several parameters might be of interest for further analysis. In order to relate the traced geometrical beams/rays to the Gaussian parameters, the publications of Arnaud, Herloski et al. and DeJager et al. are used [Arnaud1968, Herloski1983, DeJager1992](@cite). Consider the following system where a Gaussian beam with arbitrary parameters has been traced through a lens using the approach outlined in the [Complex ray tracing](@ref) section.

![Complex ray tracing II](gbtest2.png)

The user can obtain parameters such as the beam waist radius, the radius of curvature and more using the [`gauss_parameters`](@ref) and/or [`waist_parameters`](@ref) functions. Below the local waist radius and curvature ``R = r^{-1}`` have been calculated for the example above.

![Gauss beam parameters](gauss_parameters.png)