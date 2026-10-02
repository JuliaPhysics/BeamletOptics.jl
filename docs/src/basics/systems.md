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

`push!` rejects an object that is already part of the system, `delete!` an object that is not, or
that is part of an object group (remove the group instead). The index of `popat!` counts the objects
and groups as they were added, a group counts as one. Beams that were solved before the change
must be solved again from their start: call `empty!(beam)` and then
`solve_system!(system, beam; retrace = false)`. A [`StaticSystem`](@ref) can not be changed.

## Solving systems

In order to solve optical systems, this package uses a hybrid sequential and non-sequential mode. Which mode is being used is determined automatically by the [`solve_system!`](@ref) function. This is explained in more detail in the section: [Tracing logic](@ref).

```@docs; canonical=false
solve_system!
```