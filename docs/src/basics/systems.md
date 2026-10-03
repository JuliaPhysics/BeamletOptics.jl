# Optical systems

A collection of optical elements forms an optical system. Optical systems are used together with beams for the [`solve_system!`](@ref) function. Depending on the system type, different solver implementations can be activated. Currently the following types are available:

```@repl
using BeamletOptics # hide
BeamletOptics.list_subtypes(BeamletOptics.AbstractSystem); # hide
```

Refer to the **Tutorials** section for examples on how to define optical systems.

```@docs; canonical=false
System
StaticSystem
```

## Changing a system

The objects of a [`System`](@ref) can be changed after its construction, e.g. to build up a setup
step by step, starting from an empty system:

```julia
system = System()
push!(system, lens, mirror)   # objects or object groups
delete!(system, mirror)       # removes this object, compared by identity
pop!(system)                  # removes the last top-level object and returns it
popat!(system, 1)             # removes the first top-level object and returns it
```

A [`StaticSystem`](@ref) can not be changed.

```@docs; canonical=false
Base.push!(::System, ::Vararg{BeamletOptics.AbstractObject})
Base.pop!(::System)
Base.popat!(::System, ::Integer)
Base.delete!(::System, ::BeamletOptics.AbstractObject)
```

## Solving systems

In order to solve optical systems, this package uses a hybrid sequential and non-sequential mode. Which mode is being used is determined automatically by the [`solve_system!`](@ref) function. This is explained in more detail in the section: [Tracing logic](@ref).

```@docs; canonical=false
solve_system!
```

## Clear aperture and vignetting

To trace e.g. a zoom lens without vignetting, the largest bundle diameter that passes the system completely is needed. [`clear_aperture`](@ref) finds it by bisection. A ray counts as vignetted if the sequence of objects it hits differs from that of the axial ray or if it leaves the system early; this is derived from the traced geometry, so no order of the objects in the [`System`](@ref) is assumed. For an already traced source, [`vignetted`](@ref) returns the indices of the vignetted beams.

```julia
D = clear_aperture(system, [0, -0.1, 0], [0, 1, 0]; λ = 633e-9)
src = UniformDiscSource([0, -0.1, 0], [0, 1, 0], 0.9D)
```

The result is sampling-based and for a collimated bundle, see the limitations in the docstring.

```@docs; canonical=false
clear_aperture
vignetted
```
