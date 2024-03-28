#=
Types:
    AbstractEntity
        AbstractSystem
            System
        AbstractBeam
            Beam
            GaussianBeamlet
        AbstractRay
            Ray
            PolarizedRay
        AbstractInteraction
            BeamInteraction
            GaussianBeamletInteraction
        AbstractShape
            AbstractSDF
                SphereSDF
                CylinderSDF
                CutSphereSDF
                AbstractSphericalLensSDF
                    ThinLensSDF
                    BiConvexLensSDF
                    BiConcaveLensSDF
                    PlanoConvexLensSDF
                    PlanoConcaveLensSDF
                    ConvexConcaveLensSDF
            AbstractMesh
                Mesh
        AbstractObject
            AbstractReflectiveOptic
                Mirror
            AbstractRefractiveOptic
                Lens
                Prism
            AbstractDetector
                Photodetector
            AbstractBeamSplitter
                ThinBeamSplitter
            AbstractObjectGroup
                ObjectGroup
        Intersection

Core Functions:
    intersect3d(AbstractObject, AbstractRay)
    interact3d(AbstractSystem, AbstractObject, AbstractBeam, AbstractRay)
    rotate3d!(AbstractObject, axis, angle)
    translate3d!(AbstractObject, offset)
    trace_system!(system, Beam)
    trace_system!(system, GaussianBeamlet)
    retrace_system!(system, Beam)
    retrace_system!(system, GaussianBeamlet)
=#

"""
    AbstractEntity

A generic type for something that exists independently.

# Implementation reqs.
Subtypes of `AbstractEntity` should implement the following:

# Fields
- `id`: a unique identifier (UUID 4)
"""
abstract type AbstractEntity end

id(entity::AbstractEntity)::UUID = entity.id
id(::Any) = nothing

"""
    AbstractSystem <: AbstractEntity

A generic type for a container type which holds objects, beams, etc.
"""
abstract type AbstractSystem <: AbstractEntity end

refractive_index(::AbstractSystem) = 1.0

"""
    Intersection{T}

Stores some data resulting from ray tracing a `System`. This information can be used, i.e. for retracing.

# Fields:
- `t`: length of the ray parametrization in [m]
- `n`: normal vector at the point of intersection
- `id`: index of the intersected object in the `System` object-vector
"""
mutable struct Intersection{T} <: AbstractEntity
    t::T
    n::Point3{T}
    id::Nullable{UUID}
end

function Intersection(t::T, n::AbstractVector{T}, id::Nullable{UUID}) where {T}
    length(n) == 3 || throw(ArgumentError("`n` has to be of length 3"))
    Intersection(t, Point3{T}(n), id)
end

Base.length(intersection::Intersection) = intersection.t

normal3d(intersection::Intersection) = intersection.n

"""
    AbstractRay{T<:Real} <: AbstractEntity

An implementation for a geometrical optics ray in R³. In general, a `AbstractRay` is described by ``\\vec{p} + t\\cdot\\vec{d}`` with ``t\\in(0,\\infty)``. 
`AbstractRay`s are intended to model the propagation of light between optical interactions according to the laws of geometrical optics.
To store the result of a ray tracing solution, refer to [`AbstractBeam`](@ref).

# Intersections:
Since the length of a ray can not be known before solving an optical system, the [`Intersection`](@ref)-type is used.
This [`Nullable`](@ref) type can represent the intersection with an optical element, or lack thereof.

# Implementation reqs.
Subtypes of `AbstractBeam` must implement the following:

## Fields:
- `pos`: a R³-vector that stores the current position ``\\vec{p}``
- `dir`: a R³-vector that stores the current direction ``\\vec{d}``
- `intersection`: a `Nullable` field that stores the current intersection or `nothing`
- `λ`: wavelength in [m]
- `n`: refractive index along the ray length

# Additional information

!!! info "Ray length"
    Base.`length`: this function is used to return the length of the `AbstractRay` (if no intersection exists, the ray length is `Inf`)

!!! warning "Ray direction"
    Many functions assume that the `dir`ection vector has unit length (i.e. ``|\\vec{p}| = 1``). 
    Violating this assumption might lead to spurious results.
"""
abstract type AbstractRay{T <: Real} <: AbstractEntity end

position(ray::AbstractRay) = ray.pos
position!(ray::AbstractRay, pos) = (ray.pos = pos)

direction(ray::AbstractRay) = ray.dir
function direction!(ray::AbstractRay, dir)
    ray.dir = normalize(dir)
    return nothing
end

wavelength(ray::AbstractRay) = ray.λ
wavelength!(ray::AbstractRay, λ) = (ray.λ = λ)

refractive_index(ray::AbstractRay) = ray.n
refractive_index!(ray::AbstractRay, n) = (ray.n = n)

intersection(ray::AbstractRay) = ray.intersection
function intersection!(ray::AbstractRay{T}, _intersection::Nullable{Intersection{T}}) where {T}
     ray.intersection = _intersection
     return nothing
end

"""
    length(ray::AbstractRay; opl::Bool=false)

Returns the geometric length of a `ray` between its start and intersection point. If no intersection exists, `Inf` is returned.
The `opl` keyword can be used to calculate the `o`ptical `p`ath `l`ength instead, i.e. ``OPL = n \\cdot l``.
Default is the geometrical length.
"""
function Base.length(ray::AbstractRay; opl::Bool=false)
    if isnothing(intersection(ray))
        return Inf
    end
    if opl        
        return length(intersection(ray)) * refractive_index(ray)
    else
        return length(intersection(ray))
    end
end

"""
    intersect3d(plane_position, plane_normal, ray)

Returns the intersection between a `ray` and an infinitely large plane which is characterized by its `position` and `normal`.
"""
function intersect3d(plane_position::AbstractArray,
        plane_normal::AbstractArray,
        ray::AbstractRay{T}) where {T}
    t = line_plane_distance3d(plane_position, plane_normal, position(ray), direction(ray))
    isnothing(t) && return nothing
    return Intersection{T}(t, plane_normal, nothing)
end

"""
    line_point_distance3d(ray, point)

Returns value for the shortest distance between the `ray` (extended to ∞) and `point`.
"""
line_point_distance3d(ray::AbstractRay, point) = line_point_distance3d(position(ray),
    direction(ray),
    point)

"""
    angle3d(ray::AbstractRay, intersect::Intersection=intersection(ray))

Calculates the angle between a `ray` and its or some other `intersection`.
"""
function angle3d(ray::AbstractRay, intersect::Intersection = intersection(ray))
    return angle3d(direction(ray), normal3d(intersect))
end

"""
    AbstractBeam{T <: Real, R <: AbstractRay{T}} <: AbstractEntity

A generic type for a container type which holds rays, beams etc.

# Implementation reqs.
Subtypes of `AbstractBeam` must implement the following:

## Parametrization:
A subtype of `AbstractBeam` is parameterized by its main data type `T <: Real`, as well as the underlying ray representation `R <: AbstractRay{T}`.
If a beam is to be compatible with different [`AbstractRay`](@ref) implementations, it must be parameterized by `T` and `R`.
However, it can also be set to a fixed type for `T` and `R`, i.e. `MyBeam <: AbstractBeam{Float32, MyRay}`.

## Fields:
- `parent`: a [`Nullable`](@ref) field that holds the same type as the subtype, used for tree navigation
- `children`: a vector that holds the same type as the subtype, used for sub-beam tracking, i.e. beam splitting

## Functions:
- `_modify_beam_head!`: modifies the beam path for retracing purposes
- `_last_beam_intersection`: returns the last `Beam` intersection
"""
abstract type AbstractBeam{T <: Real, R <: AbstractRay{T}} <: AbstractEntity end

AbstractTrees.NodeType(::Type{<:AbstractBeam{T, R}}) where {T, R} = HasNodeType()
AbstractTrees.nodetype(beamtype::Type{<:AbstractBeam{T, R}}) where {T, R} = beamtype

AbstractTrees.parent(beam::AbstractBeam) = beam.parent
parent!(beam::B, parent::B) where {B <: AbstractBeam} = (beam.parent = parent)

AbstractTrees.children(b::AbstractBeam) = b.children

AbstractTrees.printnode(io::IO, node::B; kw...) where {B <: AbstractBeam} = show(io, B)

"""
    children!(beam::B, child::B) where {B<:AbstractBeam}

Handles the inclusion of adding a single `child` to an existing `beam`. The function behaves as follows:

1. If no previous children exist, add child
2. If `beam` already has a single child, modify child beam starting ray (retracing)
3. Else throw error
"""
function children!(beam::B, child::B) where {B <: AbstractBeam}
    if isempty(children(beam))
        # Link parent and add child to tree
        parent!(child, beam)
        push!(children(beam), child)
        return nothing
    end
    if length(children(beam)) == 1
        _modify_beam_head!(first(children(beam)), child)
        return nothing
    end
    return error("Adding child to beam failed")
end

function children!(beam::B, _children::AbstractVector{B}) where {B <: AbstractBeam}
    if isempty(children(beam))
        # Link parent and add children to tree
        parent!.(_children, Ref(beam))
        append!(children(beam), _children)
        return nothing
    end
    if length(children(beam)) == length(_children)
        for (i, child) in enumerate(children(beam))
            _modify_beam_head!(child, _children[i])
        end
        return nothing
    end
    return error("Adding children to beam failed")
end

_drop_beams!(b::B) where {B <: AbstractBeam} = (b.children = Vector{B}())

function _modify_beam_head!(::B, ::B) where {B <: AbstractBeam}
    error("_modify_beam_head not implemented for $B")
end

function _last_beam_intersection(::B) where {B <: AbstractBeam}
    error("_last_beam_intersection not implemented for $B")
end

"""
    AbstractShape{T<:Real} <: AbstractEntity

A generic type for a shape that exists in 3D-space. Must have a `pos`ition and `dir`ection.\\
Types used to describe an `object` should be subtypes of `Real`.\\

# Implementation reqs.
Subtypes of `AbstractShape` should implement the following:

## Fields:
- `pos`: a 3D-vector that stores the current `position` of the object-specific coordinate system
- `dir`: a 3x3-matrix that represents the orthonormal basis of the object and therefore, the `orientation`

## Kinematic:
- `translate3d!`: the object is moved by a translation vector relative to its current position
- `translate_to3d!`: the object is moved towards the target position
- `rotate3d!`: the object is rotated by an angle around a reference vector
- `xrotate3d!`: rotation around the x-axis
- `yrotate3d!`: rotation around the y-axis
- `zrotate3d!`: rotation around the z-axis
- `reset_translation3d!`: return the `object` to the global origin
- `reset_rotation3d!`: rotate the `object` back into its original state

## Ray Tracing:
- `intersect3d`: returns the intersection between an `AbstractShape` and `AbstractRay`, or lack thereof. See also [`Intersection`](@ref)

## Rendering (with GLMakie):
- `render_shape!`: plot the `shape` into an `Axis3` environment
- `render_shape_normals!`: plot the `shape` surface normals into an `Axis3` environment (optional)
"""
abstract type AbstractShape{T <: Real} <: AbstractEntity end

"Defines the intersection between a `shape` and a ray, defaults to no intersection."
function intersect3d(shape::AbstractShape, ::AbstractRay)
    @warn lazy"No intersect3d method defined for:" typeof(shape)
    return nothing
end

"Enforces that `object` has to have the field `pos` or implement `position()`."
position(shape::AbstractShape) = shape.pos
position!(shape::AbstractShape, pos) = (shape.pos = pos)

"Enforces that `object` has to have the field `dir` or implement `orientation()`."
orientation(shape::AbstractShape) = shape.dir
orientation!(shape::AbstractShape, dir) = (shape.dir = dir)

"""
    translate3d!(shape::AbstractShape, offset)

Translates the `pos`ition of `shape` by the `offset`-vector.
"""
function translate3d!(shape::AbstractShape, offset)
    position!(shape, position(shape) + offset)
    return nothing
end

"""
    translate_to3d!(shape::AbstractShape, target)

Translates the `shape` to the `target` position. 
"""
function translate_to3d!(shape::AbstractShape, target)
    current = position(shape)
    translate3d!(shape, target - current)
    return nothing
end

"""
    rotate3d!(shape::AbstractShape, axis, θ)

Rotates the `dir`-matrix of `shape` around the reference-`axis` by an angle of `θ`.
"""
function rotate3d!(shape::AbstractShape, axis, θ)
    R = rotate3d(axis, θ)
    orientation!(shape, orientation(shape) * R)
    return nothing
end

function xrotate3d!(shape::AbstractShape{T}, θ) where {T}
    rotate3d!(shape, Point3(one(T), zero(T), zero(T)), θ)
end
function yrotate3d!(shape::AbstractShape{T}, θ) where {T}
    rotate3d!(shape, Point3(zero(T), one(T), zero(T)), θ)
end
function zrotate3d!(shape::AbstractShape{T}, θ) where {T}
    rotate3d!(shape, Point3(zero(T), zero(T), one(T)), θ)
end

function reset_translation3d!(shape::AbstractShape{T}) where {T}
    shape.pos = Point3(zero(T))
    return nothing
end

reset_rotation3d!(shape::AbstractShape{T}) where {T} = (shape.dir = Matrix{T}(I, 3, 3))

render_shape!(::Any, ::AbstractShape) = nothing
render_shape_normals!(::Any, ::AbstractShape) = nothing

"""
    isinfrontof(shape::AbstractShape, ray::AbstractRay)

A simple test to check if a `shape` lies "in front of" a `ray`.
The forward direction is here defined as the ray `orientation`.
Only works well if `ray` is **outside** of the volume of `shape`.
Can be dispatched to return more accurate results for subtypes of `AbstractShape`.
"""
isinfrontof(shape::AbstractShape, ray::AbstractRay) = isinfrontof(position(shape),
    position(ray),
    direction(ray))

"""
    AbstractObject <: AbstractEntity

A generic type for 2D/3D objects. Optical elements are supposed to fall under this type.

# Implementation reqs.
Subtypes of `AbstractObject` must implement the following:

## Fields:
- `shape`: stores an [`AbstractShape`](@ref) that represents the object geometry, called via `shape(object)`

## Functions:
- [`intersect3d`](@ref): returns the intersection between the object and an incoming ray, defaults to the intersection with `shape(object)`
- [`interact3d`](@ref): defines the optical interaction, should return `nothing` or an [`AbstractInteraction`](@ref)
- for the kinematic API, all corresponding functions should be forwarded to the underlying [`AbstractShape`](@ref) i.e. `rotate3d!(shape(object))`
"""
abstract type AbstractObject <: AbstractEntity end

shape(object::AbstractObject) = object.shape

position(object::AbstractObject) = position(shape(object))
orientation(object::AbstractObject) = orientation(shape(object))

"""
    intersect3d(object::AbstractObject, ray::AbstractRay)

In general, the intersection between an `object` and a `ray` is defined as the intersection with `shape(object)`.
"""
function intersect3d(object::AbstractObject, ray::AbstractRay)
    intersection = intersect3d(shape(object), ray)
    # Ensure that the intersection knows about the object id
    if !isnothing(intersection)
        intersection.id = id(object)
    end
    return intersection
end

"""
    interact3d(::AbstractSystem, object::AbstractObject, ::AbstractBeam)

Defines optical interactions between beams/rays and objects, defaults to `nothing` which stops the ray tracer.
"""
function interact3d(::AbstractSystem, ::ObjectType, ::AbstractBeam, ::RayType) where {ObjectType<:AbstractObject, RayType<:AbstractRay}
    @warn lazy"No interact3d method defined for:" ObjectType RayType
    return nothing
end

translate3d!(object::AbstractObject, offset) = translate3d!(shape(object), offset)

translate_to3d!(object::AbstractObject, target) = translate_to3d!(shape(object), target)

rotate3d!(object::AbstractObject, axis, θ) = rotate3d!(shape(object), axis, θ)

xrotate3d!(object::AbstractObject, θ) = xrotate3d!(shape(object), θ)
yrotate3d!(object::AbstractObject, θ) = yrotate3d!(shape(object), θ)
zrotate3d!(object::AbstractObject, θ) = zrotate3d!(shape(object), θ)

reset_translation3d!(object::AbstractObject) = reset_translation3d!(shape(object))

reset_rotation3d!(object::AbstractObject) = reset_rotation3d!(shape(object))

"""
    AbstractObjectGroup <: AbstractObject

Container type for groups of optical elements, based on a tree-like data structure. Intended for easier kinematic handling of connected elements.

# Implementation reqs.
Subtypes of `AbstractObjectGroup` must implement the following:

## Fields:
- `objects`: stores objects or additional subgroups of objects, allows for hierarchical structures

## Functions:
- for the kinematic API, all corresponding functions of [`AbstractObject`](@ref) must be implemented
"""
abstract type AbstractObjectGroup <: AbstractObject end

AbstractTrees.children(group::AbstractObjectGroup) = group.objects

"""
    AbstractInteraction <: AbstractEntity

Describes how an `AbstractBeam` and an `AbstractSystem` interact with each other.
Contains the parameters for the continuation of the `AbstractBeam`.
Can store an optional `id` which hints at the next object intersection.
"""
abstract type AbstractInteraction <: AbstractEntity end

hint(interaction::AbstractInteraction) = interaction.id
hint(::Nothing) = nothing
