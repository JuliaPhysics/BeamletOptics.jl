# Optical systems

**FIXME** Rewrite with focus on `System` and `StaticSystem`

A collection of optical elements is referred to as a system in the context of this package. Optical systems are used together with beams for the [`solve_system!`](@ref) function. In general, optical systems must fulfill the [`BeamletOptics.AbstractSystem`](@ref) interface in order to be compatible with the standard solvers in this package:

```@docs; canonical=false
BeamletOptics.AbstractSystem
```

As is, the package provides two basic system types: [`System`](@ref) and [`StaticSystem`](@ref). Refer to the [Beam expander](@ref) tutorial for an example on how to define a simple optical system.

## Solving systems

In order to solve optical systems, this package uses a hybrid sequential and non-sequential mode. Which mode is being used is determined automatically by the [`solve_system!`](@ref) function. This is explained in more detail in the section: [Tracing logic](@ref).

```@docs; canonical=false
solve_system!
```