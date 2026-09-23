# Rendering Gaussian beamlets

Gaussian beamlets are drawn as the surface of their 1/e² intensity envelope. With `show_beams = true`, the chief, waist and divergence rays that generate the beamlet are drawn as well. Example renderings can be found in the [Stigmatic beamlets](@ref) and [Astigmatic beamlets](@ref) sections.

## Stigmatic Gaussian beamlets

```@docs
render!(::Union{GLMakie.Axis3, GLMakie.LScene}, ::GaussianBeamlet{T}) where T
```

## Astigmatic Gaussian beamlets

```@docs
render!(::Union{GLMakie.Axis3, GLMakie.LScene}, ::AstigmaticGaussianBeamlet{T}) where T
```
