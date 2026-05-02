```@setup beams
include(joinpath(@__DIR__, "..", "assets", "cond_save.jl"))

beam_showcase_dir = joinpath(@__DIR__, "..", "assets", "beam_renders")

conditional_include(joinpath(beam_showcase_dir, "beam_showcase.jl"))
conditional_include(joinpath(beam_showcase_dir, "gb_showcase.jl"))
conditional_include(joinpath(beam_showcase_dir, "collimated_sc.jl"))
conditional_include(joinpath(beam_showcase_dir, "pointsource_sc.jl"))
```

# Beams

As mentioned in the [Rays](@ref) section, a beam within the context of this package serves as a data structure for storing collections of rays, forming the backbone of the simulation framework. Beams are intended to be designed as [AbstractTrees](https://github.com/JuliaCollections/AbstractTrees.jl) to allow for ray bifurcations, e.g. in the case of optical elements such as beamsplitters. The [`solve_system!`](@ref) function relies on this data structure to perform ray tracing computations within optical systems. 

To ensure compatibility and extensibility, beam types must adhere to the [`BeamletOptics.AbstractBeam`](@ref) interface. Refer to its documentation for more information.

## Basic beam

A minimal implementation of the [`BeamletOptics.AbstractBeam`](@ref) type is provided by the [`Beam`](@ref). It can be used to store a light path through an optical system. If the beam is split, its children will be recursively traced until all paths are solved.

```@docs; canonical=false
Beam
```

A ray tracing example through an arbitrary system using a [`Beam`](@ref) is shown below. Individual [`Ray`](@ref) segments are marked by their starting position and direction. The [Beam expander](@ref) and [Miniature microscope](@ref) tutorial covers the use of the [`Beam`](@ref) in more detail. 

![Beam structure](beam_showcase.png)

## Beam groups

For convenience, the [`BeamletOptics.AbstractBeamGroup`](@ref) offers a container-like interface for groups of [`Beam`](@ref)s as commonly used in other software packages. The following concrete implementations are currently provided:

```@repl
using BeamletOptics # hide
BeamletOptics.list_subtypes(BeamletOptics.AbstractBeamGroup);
```

Refer to the following sections for convenience constructors to generate the sources listed above.

### Collimated beam source

The collimated beam source is ideal to model light coming from a focal plane at infinity. This is useful for simulating plane wavefronts. You can define a collimated monochromatic [`Beam`](@ref) source as follows:

```@docs; canonical=false
CollimatedSource(::AbstractArray{<:Real}, ::AbstractArray{<:Real}, ::Real, ::Real)
```

![Collimated group of beams](collimated_beam_source.png)

A special constructor called [`UniformDiscSource`](@ref) is available, which offers an equal-area
sampling (Fibonnaci-pattern) sampling and is thus favorable in situations where the weighting of the
individual beams becomes important, e.g. for calculating a point spread function using the [`intensity`](@ref) function.

```@docs; canonical=false
UniformDiscSource
```

![Collimated uniform group of beams](collimated_uniform_beam_source.png)

### Point beam source

The `PointSource` type is used to model emission from a spatially localized source that radiates [`Beam`](@ref)s in a range of directions. This is commonly used to simulate conical emission patterns, such as light emerging from a fiber tip or a light source for a lens objective with a known focal distance. You can specify the origin and a propagation direction, which are then used to construct the monochromatic `PointSource`.

```@docs; canonical=false
PointSource(::AbstractArray{<:Real}, ::AbstractArray{<:Real}, ::Real, ::Real)
```

Below you can find an exemplary illustration of a `PointSource`.

![Point source of beams](point_beam_source.png)

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

## Astigmatic Polarized Beamlets

Standard stigmatic beamlets are limited to symmetric systems. For modeling arbitrary optical systems—including off-axis components, tilted lenses, or complex media like atmospheric turbulence—the package provides the [`AstigmaticGaussianBeamlet`](@ref). This type implements the **Parabasal Gaussian Beamlet** formalism TODO: Add citations

### Formalism

An astigmatic beamlet is represented by a cluster of **9 rays**:
1.  **Chief Ray**: A [`PolarizedRay`](@ref) defining the central path and polarization state.
2.  **Eight Auxiliary Rays**: Four "positional" rays (waist $x+, x-, y+, y-$) and four "directional" rays (divergence $x+, x-, y+, y-$) that track the complex curvature matrix $\mathbf{Q}(z)$ through the system.

The complex scalar field at a local transverse position $\mathbf{r}$ is given by:
```math
\psi(\mathbf{r}, z) = \sqrt{\frac{A_{ref}}{A(z)}} \exp\left[ i k \left( z + \Delta l + \frac{1}{2} \mathbf{r}^T \mathbf{Q}(z) \mathbf{r} \right) \right]
```
where $A(z)$ is the complex beam area and $\Delta l$ is the accumulated optical path length difference. Note that the quadratic phase term uses the **analytic bilinear form** $\mathbf{r}^T \mathbf{u}$ rather than the conjugated dot product to ensure consistency in the complex domain.

### API and Usage

#### Construction

You can construct a single astigmatic beamlet with a specified waist and polarization:
```julia
using BeamletOptics
# A 10mm waist beamlet propagating along Y, polarized along X
agb = AstigmaticGaussianBeamlet([0,0,0], [0,1,0], 1064e-9, 10e-3; E0 = [1,0,0])
```

#### Astigmatic Beam Groups

For complex sources, the package provides the [`AstigmaticBeamGroup`](@ref) container. Several constructors are available for different scenarios:

*   [`GaussianBeamletDecomposition`](@ref): Tiling a large Gaussian beam into many small stable beamlets.
*   [`WavefrontBeamletDecomposition`](@ref): Importing an arbitrary complex field (e.g. from a phase screen or camera data).
*   [`CollimatedGaussianBeamletSource`](@ref): A square grid of parallel beamlets (ideal for aperture diffraction).
*   [`SphericalGaussianBeamletSource`](@ref): A point-like source emitting a cone of beamlets (ideal for focused/divergent beams).

#### Field Calculation

To obtain the physical vector electric field [V/m] at any point:
```julia
# Field at world position (x, y, z)
E_vec = polarized_field(agb, [0.01, 0, 0], 10.0)
```

```julia
# Decompose a field from a detector into a new AstigmaticBeamGroup
beams = WavefrontBeamletDecomposition(x, z, amplitude, phase, dir, λ)

# Solve system (propagate all beamlets to next segment)
solve_system!(system, beams)
```

### Key Assumptions
- **Paraxiality**: The auxiliary rays must remain within the paraxial regime relative to the chief ray.
- **Analytic Continuity**: Quadratic phase calculations use the bilinear dot product to avoid unphysical conjugation of complex ray parameters.
- **Power Conservation**: The decomposition methods automatically scale amplitudes to preserve total field energy regardless of grid resolution.