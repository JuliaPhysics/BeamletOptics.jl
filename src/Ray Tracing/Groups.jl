"""
    ObjectGroup{T} <: AbstractObjectGroup

A tree-like storage container for groups of objects. Can store individual objects and subgroups.
Main purpose is handling of, i.e., groups of lenses.

## Fields
- `center`: a point in 3D space which is regarded as the reference origin of the group
- `dir`: a 3x3 matrix that describes the common `orientation` of the group
- `objects`: stores [`AbstractObject`](@ref), can also store subgroups of type [`AbstractObjectGroup`](@ref) 

## Kinematic
A `ObjectGroup` implements the kinematic functions of [`AbstractObject`](@ref). The following logic is applied to

- [`translate3d!`](@ref): all objects in the group are translated by the offset vector
- [`translate_to3d!`](@ref): all objects are moved in parallel such that the group `center` is equal to the target position
- [`rotate3d!`](@ref): all objects are rotated around the `center` point with respect to their relative position
"""
mutable struct ObjectGroup{T} <: AbstractObjectGroup
    const dir::Matrix{T}
    center::Point3{T}
    const objects::Vector{AbstractObject}
end

position(group::ObjectGroup) = group.center
position!(group::ObjectGroup, pos) = (group.center = pos)

orientation(group::ObjectGroup) = group.dir
orientation!(group::ObjectGroup, dir) = (group.dir .= dir)

function ObjectGroup(v::AbstractArray{<:AbstractObject}, T = Float64)
    ObjectGroup{T}(Matrix{T}(I, 3, 3), Point3{T}(0), v)
end

"""
    objects(group::ObjectGroup)

Exposes all objects/subgroups stored within the group.
"""
objects(group::ObjectGroup) = group.objects

"""
    translate3d!(group::ObjectGroup, offset)

Moves all objects in the group by the specified `offset` vector.
"""
function translate3d!(group::ObjectGroup, offset)
    # Translate tracking vector
    position!(group, position(group) .+ offset)
    # Recursively translate all subgroups
    for object in objects(group)
        translate3d!(object, offset)
    end
    return nothing
end

"""
    translate_to3d!(group::ObjectGroup, target)

Translates all objects in parallel to the specified `target` position.
The `group` center point will be equal to the `target`.
"""
function translate_to3d!(group::ObjectGroup, target)
    current = position(group)
    translate3d!(group, target - current)
    return nothing
end

"""
    rotate3d!(group::ObjectGroup, axis, θ)

All objects in the group are rotated around the group center via the specified angle `θ` and `axis`. 
"""
function rotate3d!(group::ObjectGroup, axis, θ)
    R = rotate3d(axis, θ)
    # Update group orientation
    orientation!(group, orientation(group) * R)
    # Recursively rotate all subgroups and objects
    for object in objects(group)
        rotate3d!(object, axis, θ)
        v = position(object) - position(group)
        # Translate group around pivot point
        v = (R * v) - v
        translate3d!(object, v)
    end
    return nothing
end

function xrotate3d!(group::ObjectGroup{T}, θ) where {T}
    rotate3d!(group, @SArray(T[one(T), zero(T), zero(T)]), θ)
end
function yrotate3d!(group::ObjectGroup{T}, θ) where {T}
    rotate3d!(group, @SArray(T[zero(T), one(T), zero(T)]), θ)
end
function zrotate3d!(group::ObjectGroup{T}, θ) where {T}
    rotate3d!(group, @SArray(T[zero(T), zero(T), one(T)]), θ)
end

function reset_translation3d!(group::ObjectGroup{T}) where T
    position!(group, Point3{T}(0))
    # Recursively reset all subgroups and objects
    for object in objects(group)
        reset_translation3d!(object)
    end
    return nothing
end

function reset_rotation3d!(group::ObjectGroup{T}) where T
    orientation!(group, Matrix{T}(I, 3, 3))
    # Recursively reset all subgroups and objects
    for object in objects(group)
        reset_rotation3d!(object)
    end
    return nothing
end

function render_object!(axis::Any, group::AbstractObjectGroup)
    # Recursively render objects
    for object in objects(group)
        render_object!(axis, object)
    end
    return nothing
end

Base.show(::IO, ::MIME"text/plain", group::ObjectGroup) = print_tree(group)
