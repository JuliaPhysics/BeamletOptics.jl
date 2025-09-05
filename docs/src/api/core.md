# Core design

The BMO package is intended to provide optical simulation capabilites with as much "out-of-the-box" comfort as possible, meaning that users should not have to worry about e.g. providing the correct optical sequence of objects. The following four design principles are core assumptions of the underlying API:

1. Optical interactions are decoupled from the underlying geometry representation
2. Optical elements are closed volumes or must mimic as such
3. Without additional knowledge, tracing is performed non-sequentially
4. With additional knowledge, tracing is performed sequentially

## Intersect-Interact-Repeat-Loop

The first two principles will be elaborated upon in more detail in the [Geometry representation](@ref) section. For the latter two design decisions, the following high-level solver schematic is used to explain the steps that are performed when calling [`solve_system!`](@ref) with an input system and beam:

```@raw html
<img src="iir_loop.svg" alt="my figure" style="width:100%; height:auto;"/>
```

This scheme is loosely referred to as the **Intersect-Interact-Repeat-Loop** and consists of the following steps:

1. Calculate the closest **Intersection** between a ray/beam and the objects within the system
2. Calculate the optical **Interaction** that occurs at the surface or within the volume of the element
3. Attach or overwrite the next part of the ray chain
4. Use the new information to repeat 1.

Once this procedure has been completed, the alignment of the system or other time-dependent optical properties (e.g. the phase of a Gaussian beamlet) can be updated. When rerunning the solver, the algorithm will try to reuse information about the previously intersected objects to speed up the calculation of the next simulation step. This is described in more detail in the sections: [Tracing systems](@ref) and [Retracing systems](@ref).

The following sections will focus on the **Intersection** and **Interaction** steps.

## Intersections

## Interactions

## Hints



## CPU and GPU support

!!! info
    GPU processing (tracing) of optical systems is  not supported at the moment.



