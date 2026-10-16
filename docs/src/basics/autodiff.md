# Automatic differentiation

Results of [`solve_system!`](@ref) can be differentiated with respect to the parameters of a setup
via forward-mode automatic differentiation with
[ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl), e.g. for gradient-based optimization
or sensitivity analysis. The derivative is exact up to the accuracy of the ray tracing, typically
better than ``10^{-5}`` relative, and it is free of the step size problem of finite differences.

## Usage

BMO is generic in its number type. Write the setup as a function of the parameter and build every
component, beam and detector in the number type `T` of the parameter, which is a
`ForwardDiff.Dual` while the derivative is computed:

```julia
using BeamletOptics, ForwardDiff

# x coordinate of the spot behind a lens that is shifted laterally by `dx`
function spot_x(dx)
    T = typeof(dx)
    lens = SphericalLens(T(50e-3), T(-50e-3), T(5e-3), T(25e-3), _ -> T(1.5))
    translate3d!(lens, [dx, zero(T), zero(T)])
    pd = Detector(T(40e-3))
    translate3d!(pd, T[0, 0.1, 0])
    beam = Beam(T[0, -0.05, 0], T[0, 1, 0], T(1e-6))
    solve_system!(System([lens, pd]), beam)
    return BeamletOptics.hit_point(only(BeamletOptics.hits(pd)))[1]
end

ForwardDiff.derivative(spot_x, 0.0)   # ≈ 1.90
```

The parameter can enter anywhere: the pose of a component (`translate3d!`, `rotate3d!`), its
geometry (radii, thicknesses), a refractive index (e.g. `_ -> n0` for the parameter `n0`) or a
source. Many parameters at once are differentiated with `ForwardDiff.gradient` or
`ForwardDiff.jacobian`.

The test suite (`test/TestAutoDiff.jl`) compares the derivatives with finite differences for rays
and beamlets through spherical, aspheric, meniscus and cylindrical lenses, spherical, parabolic
and conic mirrors, prisms (including total internal reflection), thin, plate and cube
beamsplitters, polarizing filters, for polarized rays and for the optical power on a detector,
e.g. of a Michelson interferometer vs. the position of a mirror.

## Limitations

- **Discrete events are not differentiable.** Which surface a ray hits, the onset of total internal
  reflection, a ray that misses an aperture or a detector, and the cutoff by `r_max` change the
  result in steps. The derivative is that of the current path; at a step it is not defined, also
  not for finite differences.
- **Edges and corners.** At the edge of a prism or the rim of a lens the normal changes abruptly, a
  ray that hits such an edge has no meaningful derivative.
- **Number type of object groups.** An [`ObjectGroup`](@ref) is `Float64` unless its number type is
  passed: `ObjectGroup(objects, T)`.
- **Forward mode only.** The cost grows with the number of parameters. Reverse-mode AD (Zygote,
  Enzyme) is not supported, since the tracing mutates its beams and objects.
- **Compile time.** Each new combination of dual number and setup is compiled on its first call,
  which takes several seconds.
