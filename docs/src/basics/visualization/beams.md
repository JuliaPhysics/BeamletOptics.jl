# Rendering rays and beams

Rays and beams are drawn as 3D lines: all ray segments of a ray, a beam (including its child beams) or a beam group are drawn by a single `linesegments` plot. A ray without an intersection is drawn with the finite length `flen`, since its actual length is infinite. [`live_render!`](@ref) draws the same plot and updates it after a solve, see [Live rendering](@ref). Example renderings can be found in the [Basic rays](@ref), [Basic beam](@ref) and [Beam groups](../beams/beam_groups.md) sections.

## Single rays

```@docs
render!(::Union{GLMakie.Axis3, GLMakie.LScene}, ::BeamletOptics.AbstractRay)
```

## Beams of rays

A [`Beam`](@ref) is rendered by drawing every ray of the beam tree, including all child beams created at beamsplitters. Keyword arguments such as `color` or `linewidth` are passed on to the plot of the segments.

```@docs
render!(::Union{GLMakie.Axis3, GLMakie.LScene}, ::Beam)
```

## Groups of beams

```@docs
render!(::Union{GLMakie.Axis3, GLMakie.LScene}, ::BeamletOptics.AbstractBeamGroup)
```

## Coloring by wavelength

`color = :wavelength` draws each ray in the display color of its own wavelength, e.g. to show white light that is dispersed by a prism: child beams that leave a dispersive element keep the color of their wavelength, and every beam of a beam group gets its own color. The opacity is set by `color = (:wavelength, 0.3)`. The same option colors the envelope of Gaussian beamlets, segment by segment. [`live_render!`](@ref) takes it as well.

```julia
render!(ax, beam_group; color = :wavelength)
render!(ax, gauss; color = (:wavelength, 0.5))
```

The color of a wavelength is available without a plotting backend via [`wavelength_color`](@ref), e.g. to color other plots consistently:

```@docs
wavelength_color
```

## Polarization overlay

For a [`PolarizedRay`](@ref) or a `Beam` of polarized rays, `render!(ax, beam; show_polarization=true)` overlays a curve that shows the electric field component perpendicular to the ray direction along the optical path. The appearance of the curve is controlled via the `pol_*` keyword arguments listed above. Passing `show_polarization=true` for non-polarized rays throws an `ArgumentError`. An example is shown in the [Polarized rays](@ref) section.

```julia
render!(ax, beam; show_polarization=true, pol_color=:orange)
```
