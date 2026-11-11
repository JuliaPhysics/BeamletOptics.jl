# Polarizers

Polarizers in the context of this package are optical elements that select or modify the R³ polarization vector of a [`PolarizedRay`](@ref), e.g. filters. This is mainly done by two approaches:

1. 3D polarization ray-tracing calculus
2. 3D modified Jones matrix calculus

For more information on the first method refer to the section: [Polarized rays](@ref). For the second approach, elements fall under the category of the [`BeamletOptics.AbstractJonesPolarizer`](@ref).

```@docs; canonical=false
BeamletOptics.AbstractJonesPolarizer
```

## Jones matrix element representation

Fundamentally, the approach used here to simulate the effect of polarizers is referred to as [Jones calculus](https://en.wikipedia.org/wiki/Jones_calculus) and gives a "0th order" approximation of the physical effect. Out of plane tilts with respect to the optical axis of an incoming ray are currently only considered via a projection into the transverse plane of the incoming ray [Korger2013](@cite).

In a nutshell, elements are characterized by a 2x2 matrix that determines how the E-field components in the transverse plane to the optical axis are passed through in a global coordinate system where a ray of polarized light propagates along the z-axis. For instance, the entries for a linear filter that blocks in the y-direction are

```math
J = 
\begin{pmatrix}
1 & 0\\
0 & 0
\end{pmatrix} \,.
```

While this allows to simply model a specific set of polarizing elements, its important to note that more complex phenomena need more extensive implementations. For 3D-calculations, the ``J``-matrix representation is embedded into 

```math
P =
\begin{bmatrix}
  J & \begin{matrix}0 \\ 0\end{matrix} \\
  \begin{matrix}0 & 0\end{matrix} & 1
\end{bmatrix}
```

in order to calculate ``\vec{E}_1 = P \cdot \vec{E}_0`` for normal incidence. Additional details are provided in the docs above.

## Polarisation filter

A polarisation filter or linear polarizer is the simplest practical polarizer and is commonly used to select a desired polarization state. This package provides the [`PolarizationFilter`](@ref) as an idealized implementation for a zero-thickness filter.

The keyword `cutoff_strength` (default `eps()`) is a threshold on the norm of the polarization vector behind the filter: if the norm is at most `cutoff_strength`, the ray is blocked and its beam ends at the filter. A ray that is not a [`PolarizedRay`](@ref), e.g. an unpolarized [`Ray`](@ref), passes the filter unchanged.

```@docs; canonical=false
PolarizationFilter(::Real)
RoundPolarizationFilter
transmission_axis(::PolarizationFilter)
```

The transmission axis of a filter can be queried via [`transmission_axis`](@ref). When rendered, both the [`PolarizationFilter`](@ref) and the [`LinearPolarizer`](@ref) show this axis as a line across the film; the [`LinearPolarizer`](@ref) additionally shows rim marks akin to the engraved line of real components. This can be disabled via `render!(ax, filter; show_transmission_axis=false)`, and styled with `axis_color` and `axis_linewidth`.

## Linear polarizer

The [`LinearPolarizer`](@ref) models a real film polarizer (e.g. Thorlabs [LPNIRE100-B](https://www.thorlabs.com/item/LPNIRE100-B?aID=af1e1f34f538197c483fc7db11fa2b22&aC=1)): a 2D [`PolarizationFilter`](@ref) cemented between two flush glass plates. The film acts at the inner glass interface, the uncoated outer surfaces cause Fresnel losses, and the film thickness itself is not modeled. An example is shown in the [Polarized rays](@ref) section.

```@docs; canonical=false
LinearPolarizer
RoundLinearPolarizer
transmission_axis(::LinearPolarizer)
```

## Wave plates

A [`QuarterWavePlate`](@ref) or [`HalfWavePlate`](@ref) is an infinitesimally thin retarder film, a
[`CoatingLayer`](@ref) with a [`Retarder`](@ref). Its fast axis lies along x; as in the lab, the plate is
turned in its mount about its normal to set the angle of the fast axis.

```@example waveplates
using BeamletOptics
const mm = 1e-3

hwp = HalfWavePlate(25.4mm)
yrotate3d!(hwp, deg2rad(22.5))   # turns linear x by 45°
translate3d!(hwp, [0, 20mm, 0])

beam = Beam(PolarizedRay([0, 0, 0], [0, 1, 0], 633e-9, [1, 0, 0]))
solve_system!(System([hwp]), beam)
BeamletOptics.polarization(last(rays(beam)))
```

A zero-order plate is a quarter or half wave at its design wavelength `λ0` only. Films that act on a
face of a component, e.g. a retarder under a reflective polarizer, are described on the
[Coatings](@ref) page.

```@docs; canonical=false
QuarterWavePlate
HalfWavePlate
fast_axis
```
