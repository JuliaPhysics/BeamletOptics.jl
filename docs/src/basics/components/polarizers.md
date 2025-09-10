# Polarizers

Polarizers in the context of this package are optical elements that select or modify the R³ polarization vector of a [`PolarizedRay`](@ref), e.g. filters or λ/2 waveplates. This is mainly done by two approaches:

1. 3D polarization ray-tracing calculus
2. 3D modified Jones matrix calculus

For more information on the first method refer to the section: [Polarized Rays](@ref). For the second approach, elements fall under the category of the [`BeamletOptics.AbstractJonesPolarizer`](@ref).

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

```@docs; canonical=false
PolarizationFilter(::Real)
```