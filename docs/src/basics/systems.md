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
Base.push!(::System, ::Vararg{BeamletOptics.ObjectOrGroup})
Base.pop!(::System)
Base.popat!(::System, ::Integer)
Base.delete!(::System, ::BeamletOptics.ObjectOrGroup)
```

## Solving systems

Optical systems are solved non-sequentially by the [`solve_system!`](@ref) function: for every ray the solver searches the system for the object that is hit next. A component can shortcut this search for the following ray by returning a [`BeamletOptics.Hint`](@ref). This is explained in more detail in the section: [Tracing logic](@ref).

```@docs; canonical=false
solve_system!
```

Some objects store data during a solve, e.g. a [`Detector`](@ref) its hits. This data is kept between solves unless the system is initialized:

```@docs; canonical=false
initialize!
```