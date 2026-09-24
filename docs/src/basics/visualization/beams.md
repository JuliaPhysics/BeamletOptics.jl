# Rendering rays and beams

Rays and beams are drawn as 3D lines. A ray without an intersection is drawn with the finite length `flen`, since its actual length is infinite. Example renderings can be found in the [Basic rays](@ref), [Basic beam](@ref) and [Beam groups](../beams/beam_groups.md) sections.

## Single rays

```@docs
render!(::Union{GLMakie.Axis3, GLMakie.LScene}, ::BeamletOptics.AbstractRay)
```

## Beams of rays

A [`Beam`](@ref) is rendered by drawing every ray of the beam tree, including all child beams created at beamsplitters. Keyword arguments such as `color` or `linewidth` are passed on to each ray.

```@docs
render!(::Union{GLMakie.Axis3, GLMakie.LScene}, ::Beam)
```

## Groups of beams

```@docs
render!(::Union{GLMakie.Axis3, GLMakie.LScene}, ::BeamletOptics.AbstractBeamGroup)
```

## Polarization overlay

For a [`PolarizedRay`](@ref) or a `Beam` of polarized rays, `render!(ax, beam; show_polarization=true)` overlays a curve that shows the electric field component perpendicular to the ray direction along the optical path. The appearance of the curve is controlled via the `pol_*` keyword arguments listed above. Passing `show_polarization=true` for non-polarized rays throws an `ArgumentError`. An example is shown in the [Polarized rays](@ref) section.

```julia
render!(ax, beam; show_polarization=true, pol_color=:orange)
```
