# Rendering Gaussian beamlets

Gaussian beamlets are drawn as a mesh of their 1/e² intensity envelope, one mesh for all segments of a beamlet (or of a group of astigmatic beamlets). [`live_render!`](@ref) draws the same mesh and updates it after a solve, see [Live rendering](@ref). With `show_beams = true`, the chief, waist and divergence rays that generate the beamlet are drawn as well. Example renderings can be found in the [Stigmatic beamlets](@ref) and [Astigmatic beamlets](@ref) sections.

## Stigmatic Gaussian beamlets

```@docs
render!(::Union{GLMakie.Axis3, GLMakie.LScene}, ::GaussianBeamlet)
```

## Astigmatic Gaussian beamlets

```@docs
render!(::Union{GLMakie.Axis3, GLMakie.LScene}, ::AstigmaticGaussianBeamlet)
```

## Groups of astigmatic beamlets

```@docs
render!(::Union{GLMakie.Axis3, GLMakie.LScene}, ::BeamletOptics.AstigmaticBeamGroup)
```
