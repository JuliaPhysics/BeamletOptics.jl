# Rays

Individual monochromatic rays form the basic building blocks to describe the propagation of light through an optical system using geometrical optics. In general, the ray path in the context of this package is described by 

```math
\vec{x}(t) = \vec{p} + t \cdot \vec{d}
```

where ``\vec{p}`` and ``\vec{d}`` are the position and direction ``\mathbb{R}^3``-vectors, respectively. The ray length ``t`` is used to describe the geometrical length of the ray. This assumes that the [`SCDI.RefractiveIndex](@ref) along the ray path is constant. If after solving an optical system a ray intersection is determined, a new ray must be spawned to model an arbitrary light path. This data is stored, e.g., in a [`SCDI.Beam`](@ref). More on this can be found in the [Beams](@ref) chapter. 

## Basic `Ray`

The generic type that describes geometrical rays is [`SCDI.AbstractRay`](@ref). Refer to its documentation for more information about what data is used to model light propagation. A minimal implementation of this API (i.e. subtype) is provided by [`SCDI.Ray`](@ref):

```@docs; canonical=false
SCDI.Ray
```

This ray type is able to model reflection and [refraction](https://www.rp-photonics.com/refraction.html), using Snell's law, within the paraxial approximation. A single ray only ever describes the optical path between its starting point and the closest intersection, e.g. the surface of a [`SCDI.Lens`](@ref).

## Polarized Rays

In order to model the effect of polarizing elements, e.g. a ``\lambda/4``-plate, the polarization ray tracing calculus of Yun et. al is used [Yun2011_1, Yun2011_2](@cite). This formalism allows to model the effects of said elements on the electric field vector ``E_0`` using the [Jones formalism](https://www.rp-photonics.com/polarization_of_light.html) in global coordinates:

```@docs; canonical=false
SCDI.PolarizedRay
```

### Fresnel coefficients

This package uses the equations of Fowles [Fowles1989; p. 44](@cite) and Peatross [Peatross2015; p. 78](@cite) to determine the [Fresnel coefficents](https://www.rp-photonics.com/fresnel_equations.html) at a given surface where the [`SCDI.PolarizedRay`](@ref) enters from a medium with (complex) refractive index ``n_1`` into a medium with ``n_2``. The ability to define coatings is currently not included. 

!!! warning
    When a [`SCDI.PolarizedRay`](@ref) interacts with a refractive medium, e.g. an [`SCDI.AbstractRefractiveOptic`](@ref), the default tracing behaviour is to only trace the refracted and ignore the reflected ray, unless [Total Internal Reflection (TIR)](https://www.rp-photonics.com/total_internal_reflection.html) occurs. 

Below the Fresnel coefficients for different ``n_1 \rightarrow n_2`` interfaces are listed to highlight the definition of signs used for this package.

#### Vacuum to glass

First, the Fresnel coefficients for ``n_1 = 1.0`` to ``n_2 = 1.5`` will be calculated. The angle of incidence ``\theta`` refers to the plane of incidence in the `s`enkrecht and `p`arallel coordinate system. Note that the imaginary part of the coefficents is shown by the dash-dotted lines.

```@example fresnel_vacuum_glass
using CairoMakie #hide
using SCDI
include("fresnel.jl") # hide

# Angle of incidence
θ = deg2rad.(0:.01:90)

# Define refractive indices - vacuum to glass
n1 = 1.0
n2 = 1.5

# Calculate complex Fresnel coefficients
rs, rp, ts, tp = SCDI.fresnel_coefficients(θ, n2/n1)

plot_and_save_fresnel_coeffs(n1, n2, save_fig=false) # hide
```

Note that for this example, the imaginary part of the coefficients is zero for all considered `θ`s.

#### Glass to vacuum

For a glass-vacuum interface with ``n_1 = 1.5`` to ``n_2 = 1.0`` the coefficients are calculated likewise. Note the unsteadiness of the coefficients at around 40°. This is the critical angle where [TIR](https://www.rp-photonics.com/total_internal_reflection.html) occurs.

```@example fresnel_vacuum_glass
# Define refractive indices - glass to vacuum
n1 = 1.5
n2 = 1.0

# Calculate complex Fresnel coefficients
rs, rp, ts, tp = SCDI.fresnel_coefficients(θ, n2/n1)

plot_and_save_fresnel_coeffs(n1, n2, save_fig=false) # hide
```