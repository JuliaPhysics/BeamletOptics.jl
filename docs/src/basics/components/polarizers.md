# Polarizers

Polarizers in the context of this package are optical elements that select or modify the R³ polarization vector of a [`PolarizedRay`](@ref), e.g. filters or λ/2 waveplates. This is mainly done by two approaches:

1. 3D polarization ray-tracing calculus
2. 3D modified Jones matrix calculus

For more information on the first method refer to the section: [Polarized rays](@ref). For the second approach, elements use a [`Coating`](@ref) wrapping a [`JonesCoating`](@ref) model.

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

---

## Waveplates

Waveplates (or retarders) introduce a phase shift $\Gamma$ between two orthogonal polarization components. The fast axis is aligned along the local $x$-axis and the slow axis along the local $z$-axis (with propagation along the local $y$-axis), resulting in the local retardance matrix:

```math
J_{\text{retarder}} =
\begin{pmatrix}
e^{-i\Gamma/2} & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & e^{i\Gamma/2}
\end{pmatrix} \,.
```

This package supports both **flat** (zero-thickness) and **thick/plate** waveplates (which include a bulk glass substrate of finite thickness $d$ and index $n$). The constructors automatically distinguish between **round** shapes (by passing a single `diameter` value) and **rectangular** shapes (by passing separate `width` and `height` values).

```@docs; canonical=false
Waveplate
HalfWavePlate
QuarterWavePlate
RectangularPlateWaveplate
RoundPlateWaveplate
```

---

## Polarizing Cube Beamsplitter

A Polarizing Cube Beamsplitter (PBS) consists of two joined right-angle prisms, separating the s- and p-polarization components of light at the hypotenuse interface: s-polarized light is completely reflected at 90°, while p-polarized light is completely transmitted.

```@docs; canonical=false
PolarizingCubeBeamsplitter
```