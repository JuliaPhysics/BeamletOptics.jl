```@setup beam_groups
beam_showcase_dir = joinpath(@__DIR__, "..", "..", "assets", "beam_renders")

Main.DocUtils.conditional_include(joinpath(beam_showcase_dir, "collimated_sc.jl"), use_placeholder=true)
Main.DocUtils.conditional_include(joinpath(beam_showcase_dir, "pointsource_sc.jl"), use_placeholder=true)
```

# Beam groups

For convenience, the [`BeamletOptics.AbstractBeamGroup`](@ref) offers a container-like interface for groups of [`Beam`](@ref)s as commonly used in other software packages. The following concrete implementations are currently provided:

```@repl
using BeamletOptics # hide
BeamletOptics.list_subtypes(BeamletOptics.AbstractBeamGroup);
```

Refer to the following sections for convenience constructors to generate the sources listed above.

## Collimated beam source

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

## Point beam source

The `PointSource` type is used to model emission from a spatially localized source that radiates [`Beam`](@ref)s in a range of directions. This is commonly used to simulate conical emission patterns, such as light emerging from a fiber tip or a light source for a lens objective with a known focal distance. You can specify the origin and a propagation direction, which are then used to construct the monochromatic `PointSource`.

```@docs; canonical=false
PointSource(::AbstractArray{<:Real}, ::AbstractArray{<:Real}, ::Real, ::Real)
```

Below you can find an exemplary illustration of a `PointSource`.

![Point source of beams](point_beam_source.png)

## Astigmatic Beam Groups

For complex sources, the package provides the [`AstigmaticBeamGroup`](@ref) container. Several constructors are available for different scenarios:

- [`GaussianBeamletDecomposition`](@ref): Tiling a large Gaussian beam into many small stable beamlets.
- [`WavefrontBeamletDecomposition`](@ref): Importing an arbitrary complex field (e.g. from a phase screen or camera data).
- [`CollimatedGaussianBeamletSource`](@ref): A square grid of parallel beamlets (ideal for aperture diffraction).
- [`SphericalGaussianBeamletSource`](@ref): A point-like source emitting a cone of beamlets (ideal for focused/divergent beams).
- [`EllipticalGaussianBeamletSource`](@ref): A point-like source emitting an elliptical cone of beamlets (ideal for sources with different fast/slow axis divergence).