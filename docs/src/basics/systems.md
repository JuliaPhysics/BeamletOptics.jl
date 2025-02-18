# Optical systems

A collection of optical elements is referred to as a system in the context of this package. Optical systems are used together with beams for the [`SCDI.solve_system!`](@ref) function. In general, optical systems must fulfill the [`SCDI.AbstractSystem`](@ref) interface in order to be compatible with the standard solvers in this package:

```@docs; canonical=false
SCDI.AbstractSystem
```

As is, the package provides two basic system types: [`SCDI.System`](@ref) and [`SCDI.StaticSystem`](@ref). Refer to the [Beam expander](@ref) tutorial for an example on how to define a simple optical system.

## Tracing logic

In order to solve optical systems, this package uses a hybrid sequential and non-sequential mode. Which mode is being used is determined automatically by the [`SCDI.solve_system!`](@ref) function. This will be explained in more detail below.

```@docs; canonical=false
SCDI.solve_system!
```

### Tracing systems

In the initial state, is is assumed that the problem consists of `objects` <: [`SCDI.AbstractObject`](@ref)s (in a system) and a `beam` <: [`SCDI.AbstractBeam`](@ref) with a defined starting position and direction. No additional information is provided, and the specific path of the beam is not known beforehand. Consequently, brute force tracing of the optical system is required, involving testing against each individual element to determine the trajectory of the beam.

This non-sequential mode is comparatively safe in determining the "true" beam path, but will scale suboptimally in time-complexity with the amount of optical elements. After solving the system, the beam path is known and can be potentially reused in the future.

!!! info "Object order"
    Unlike with classic, surface-based ray tracers, the order in which objects are listed in the [`SCDI.System`](@ref) object vector/tuple is not considered for the purpose of tracing or retracing.

### Retracing systems

Once a system has been traced for the first time, the system and beam can be solved again. However, this time the solver will try to reuse as much information from the previous run as possible by testing if the previous beam trajectory is still valid in a sequential tracing mode. Retracing systems assumes that the kinematic changes (e.g. optomechanical aligment) between the current tracing procedure and the previous one are small. If an intersection along the beam trajectory becomes invalid, the solver will perform a non-sequential trace for all invalidated parts of the beam.

!!! warning "Retracing blocked beam paths"
    The  implemented standard retracing procedure can handle beam path invalidations under certain conditions. However, one case that will lead to a **silent error** is if an element in the system is moved such that it **blocks the beam path between two other elements**. The retracer will not be able to detect this, since the testing of the previous intersection will return a valid intersection.

    If this kind of situation must be modeled, e.g. in the case of an optical chopper wheel, retracing should be disabled.