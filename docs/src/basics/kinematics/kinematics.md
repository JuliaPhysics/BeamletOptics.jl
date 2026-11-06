# Kinematics

Optical elements, rays, beams and beam groups can move around freely in three-dimensional space, which enables the modeling of kinematics within optical setups, e.g. the alignment of a mirror or the displacement of a light source. All of them share one common kinematic API, which is explained on this page. Particularities of the individual entities are covered in the following sections.

```@contents
Pages = ["objects.md", "sources.md"]
Depth = 2
```

!!! important "Optical system kinematics"
    Elements and sources can be moved freely between each call of [`solve_system!`](@ref). However, during tracing it is assumed that all elements remain static.

## Position, orientation and direction

Every movable entity has a [`position`](@ref), represented as a ``\mathbb{R}^3``-vector, which serves as its reference point for all movements. For optical elements, this is the origin of their local coordinate system, which the constructor of the element defines (e.g. the front vertex of a lens); for rays and beams, it is the start of the light path; for beam groups, it is the pivot `center` of the group, not the start of a light path. Entities with a full local coordinate system additionally have an [`orientation`](@ref), represented by an orthonormal matrix in ``\mathbb{R}^3``. Its columns are the local axes in global coordinates, so the matrix transforms vectors from local to global coordinates and its transpose transforms from global to local coordinates. [`direction`](@ref) returns the local y-axis of such an entity, i.e. its optical axis, as `orientation(x)[:, 2]`. Rays and beams only have a [`direction`](@ref), but no orientation.

## Movement commands

The following movement commands are provided:

- Translation
    - [`translate3d!`](@ref)
    - [`translate_to3d!`](@ref)
- Rotation
    - [`rotate3d!`](@ref)
    - [`xrotate3d!`](@ref)
    - [`yrotate3d!`](@ref)
    - [`zrotate3d!`](@ref)
    - [`align3d!`](@ref)
- Reset commands
    - [`reset_translation3d!`](@ref)
    - [`reset_rotation3d!`](@ref) (not available for rays and beams)
- Pivot commands (groups only)
    - [`set_pivot3d!`](@ref)

!!! important "Relative motion"
    Unless specified otherwise, the translation and rotation commands result in relative motions to the current position and orientation. This must be taken into account when trying to model a specific set of movements.

Rotations are performed about the [`position`](@ref) of the moved entity, unless a `pivot` is passed to [`rotate3d!`](@ref) explicitly. Angles are given in radians, and the axes of [`xrotate3d!`](@ref), [`yrotate3d!`](@ref) and [`zrotate3d!`](@ref) are the global coordinate axes (see [Conventions](@ref)).

!!! note "Alignment of oriented entities"
    [`align3d!`](@ref) rotates an entity with an [`orientation`](@ref) by the shortest rotation that brings its [`direction`](@ref) onto the target. How the entity is twisted about the new direction therefore depends on its previous orientation.

```@docs; canonical=false
translate3d!(::Any, ::Any)
translate_to3d!(::Any, ::Any)
rotate3d!(::Any, ::AbstractMatrix)
xrotate3d!(::Any, ::Any)
yrotate3d!(::Any, ::Any)
zrotate3d!(::Any, ::Any)
align3d!(::Any, ::Any)
reset_translation3d!(::Any)
reset_rotation3d!(::Any)
```

## Static and movable entities

All shapes, optical elements, object groups, rays, beams and beam groups provided by BMO can be moved, whereas a [`System`](@ref) is static. A custom type can also be declared static, in which case every movement command throws an `ArgumentError`. [`set_pivot3d!`](@ref) is the exception, it does not check the trait. How the kinematic API is implemented, and how a new element or beam type opts in or out of it, is explained in the [Kinematic system](@ref) section of the developer documentation.
