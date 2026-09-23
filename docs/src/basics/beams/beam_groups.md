```@setup beam_groups
beam_showcase_dir = joinpath(@__DIR__, "..", "..", "assets", "beam_renders")

Main.DocUtils.conditional_include(joinpath(beam_showcase_dir, "collimated_sc.jl"))
Main.DocUtils.conditional_include(joinpath(beam_showcase_dir, "pointsource_sc.jl"))
```

# Beam groups

For convenience, the [`BeamletOptics.AbstractBeamGroup`](@ref) offers a container-like interface for groups of [`Beam`](@ref)s as commonly used in other software packages. The following concrete implementations are currently provided:

```@repl
using BeamletOptics # hide
BeamletOptics.list_subtypes(BeamletOptics.AbstractBeamGroup);
```

Refer to the following sections for convenience constructors to generate the sources listed above.

!!! tip "Moving a beam group"
    Every beam group has a source position (`center`) and an [`orientation`](@ref) matrix, whose second column is the central source [`direction`](@ref) and whose first column is the sampling reference vector (`basis`). Both are updated when the group is moved. Refer to [`BeamletOptics.AbstractBeamGroup`](@ref) for the convention and to [Moving sources](@ref) for the kinematic API.

## Collimated beam source

The collimated beam source is ideal to model light coming from a focal plane at infinity. This is useful for simulating plane wavefronts. You can define a collimated monochromatic [`Beam`](@ref) source as follows:

```@docs; canonical=false
CollimatedSource(::AbstractArray{<:Real}, ::AbstractArray{<:Real}, ::Real, ::Real)
```

![Collimated group of beams](collimated_beam_source.png)

Already existing beams, e.g. [`Beam`](@ref)s of [`PolarizedRay`](@ref)s as in the [Vectorial focusing at high NA](@ref) example, can be wrapped into a group with `CollimatedSource(beams, diameter, pos, dir)`. `PointSource(beams, NA, pos, dir)` and `AstigmaticBeamGroup(beams, pos, dir)` work analogously; instead of `dir`, a full `orientation` matrix can be passed.

A special constructor called [`UniformDiscSource`](@ref) is available, which offers an equal-area
sampling (Fibonacci pattern) and is thus favorable in situations where the weighting of the
individual beams becomes important, e.g. for calculating a point spread function using the [`intensity`](@ref) function.

```@docs; canonical=false
UniformDiscSource
```

![Collimated uniform group of beams](collimated_uniform_beam_source.png)

## Point beam source

The `PointSource` type is used to model emission from a spatially localized source that radiates [`Beam`](@ref)s in a range of directions. This is commonly used to simulate conical emission patterns, such as light emerging from a fiber tip or a light source for a lens objective with a known focal distance. You can specify the origin and a propagation direction, which are then used to construct the monochromatic `PointSource`.

```@docs; canonical=false
PointSource(::AbstractArray{<:Real}, ::AbstractArray{<:Real}, ::Real, ::Real)
```

Below you can find an exemplary illustration of a `PointSource`.

![Point source of beams](point_beam_source.png)

Analogous to [`UniformDiscSource`](@ref), a special constructor called [`UniformPointSource`](@ref)
is available, which samples the spherical cap around `dir` with an equal solid angle per ray
(Fibonacci/sunflower pattern), rather than in concentric rings. It returns a plain [`PointSource`](@ref)
and, unlike the ring-based constructor, has no dedicated center beam along `dir`.

```@docs; canonical=false
UniformPointSource
```

## Astigmatic Beam Groups

For complex sources, the package provides the [`AstigmaticBeamGroup`](@ref) container. Several constructors are available for different scenarios:

- [`GaussianBeamletDecomposition`](@ref): Tiling a large Gaussian beam into many small stable beamlets.
- [`WavefrontBeamletDecomposition`](@ref): Importing an arbitrary complex field (e.g. from a phase screen or camera data).
- [`CollimatedGaussianBeamletSource`](@ref): A square grid of parallel beamlets (ideal for aperture diffraction).
- [`SphericalGaussianBeamletSource`](@ref): A point-like source emitting a cone of beamlets (ideal for focused/divergent beams).
- [`EllipticalGaussianBeamletSource`](@ref): A point-like source emitting an elliptical cone of beamlets (ideal for sources with different fast/slow axis divergence).