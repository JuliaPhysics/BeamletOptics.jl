"""
    AbstractShapeTrait

The shape trait defines how many shapes an [`AbstractObject`](@ref) consists of. 
Two different traits are defined:

1. [`SingleShape`](@ref): the `AbstractObject` consists of a single [`AbstractShape`](@ref)
2. [`MultiShape`](@ref): the `AbstractObject` consists of two or more [`AbstractShape`](@ref)s

Refer to the respective documentation for more information
"""
abstract type AbstractShapeTrait end

function shape_getter_not_implemented_error(T::Type, O::Type)
    throw(ErrorException("Shape getter of $O not implemented! Refer to $T documentation"))
end

"""
    SingleShape <: AbstractShapeTrait

Represents that the [`AbstractObject`](@ref) consists of a single underlying shape.

# AbstractObject implementation reqs.

If `shape_trait_of(::Foo) = SingleShape()` is defined, `Foo` must implement the following: 

## Fields

- `shape`: a single concrete [`AbstractShape`](@ref), e.g. a [`CylinderSDF`](@ref)
"""
struct SingleShape <: AbstractShapeTrait end

shape(::SingleShape, object::AbstractObject) = object.shape

Base.position(::SingleShape, object::AbstractObject) = position(shape(object))
position!(::SingleShape, object::AbstractObject, pos) = position!(shape(object), pos)

orientation(::SingleShape, object::AbstractObject) = orientation(shape(object))
orientation!(::SingleShape, object::AbstractObject, dir) = orientation!(shape(object), dir)

translate3d!(::SingleShape, object::AbstractObject, offset) = translate3d!(shape(object), offset)

translate_to3d!(::SingleShape, object::AbstractObject, target) = translate_to3d!(shape(object), target)

rotate3d!(::SingleShape, object::AbstractObject, axis, θ) = rotate3d!(shape(object), axis, θ)

align3d!(::SingleShape, object::AbstractObject, axis) = align3d!(shape(object), axis)

reset_translation3d!(::SingleShape, object::AbstractObject) = reset_translation3d!(shape(object))

reset_rotation3d!(::SingleShape, object::AbstractObject) = reset_rotation3d!(shape(object))

"""
    MultiShape <: AbstractShapeTrait

Represents that the [`AbstractObject`](@ref) consists of a two or more [`AbstractShape`](@ref)s.

# AbstractObject implementation reqs.

If `shape_trait_of(::Foo) = MultiShape()` is defined, `Foo` must implement the following: 

## Functions

- `shape(::Foo)`: a getter function that returns a `Tuple` of all relevant shapes, e.g. `(foo.front, foo.middle, foo.back)`

# Additional information

!!! warn "Kinematic center"
    Unless specified otherwise by dispatching `position` / `position!` and `orientation` / `orientation!`
    onto custom `pos` and `dir` data fields, the **position and orientation of the first element** returned
    by `shape(object)` will be used as the **kinematic center** for e.g. `translate3d!`.
""" 
struct MultiShape <: AbstractShapeTrait end

shape(::MultiShape, ::O) where O <: AbstractObject = shape_getter_not_implemented_error(MultiShape, O)

Base.position(::MultiShape, object::AbstractObject) = position(first(shape(object)))
position!(::MultiShape, object::AbstractObject, ::Any) = nothing

orientation(::MultiShape, object::AbstractObject) = orientation(first(shape(object)))
orientation!(::MultiShape, object::AbstractObject, ::Any) = nothing

"""
    translate3d!(::MultiShape, object, offset)

Moves all parts of the [`MultiShape`](@ref) `object` along the specified `offset` vector.
"""
function translate3d!(::MultiShape, object::AbstractObject, offset)
    # Translate tracking vector
    position!(object, position(object) .+ offset)
    # Recursively translate all subparts
    for subpart in shape(object)
        translate3d!(subpart, offset)
    end
    return nothing
end

"""
    translate_to3d!(::MultiShape, object, target)

Translates all parts of the [`MultiShape`](@ref) `object` in parallel to the specified `target` position.
The `object` center point will be equal to the `target`.
"""
function translate_to3d!(::MultiShape, object::AbstractObject, target)
    current = position(object)
    translate3d!(object, target .- current)
    return nothing
end

"""
    rotate3d!(::MultiShape, object, axis, θ)

All parts of the [`MultiShape`](@ref) `object` are rotated around the pivot center via the specified angle `θ` and `axis`.
"""
function rotate3d!(::MultiShape, object::AbstractObject, axis, θ)
    R = rotate3d(axis, θ)
    # Update group orientation
    orientation!(object, R * orientation(object))
    # Recursively rotate all subgroups and objects
    for subpart in shape(object)
        rotate3d!(subpart, axis, θ)
        v = position(subpart) .- position(object)
        # Translate group around pivot point
        v = (R * v) - v
        translate3d!(subpart, v)
    end
    return nothing
end

function align3d!(::MultiShape, object::AbstractObject, target_vec)
    # FIXME
    @warn "align3d! not yet implemented for MultiShape"
    return nothing
end

"""
    reset_translation3d!(::MultiShape, object)

Resets all applied translations of the `object`, i.e. moves the center back to the origin.

!!! info "Parts within parts"
    Sub-part relative translations are not reset!
"""
function reset_translation3d!(::MultiShape, object::AbstractObject{T}) where T
    # Reset object back to origin
    translate3d!(object, -position(object))
    # Reset center of kinematics (removes precision artifacts)
    position!(object, Point3{T}(0))
    return nothing
end

"""
    reset_rotation3d!(::MultiShape, object)

Reset all applied rotations of the `object`, i.e. resets the local coordinate system to the standard base.

!!! info "Parts within parts"
    Sub-part relative rotations are not reset!
"""
function reset_rotation3d!(::MultiShape, object::AbstractObject{T}) where T
    # Calculate rotation reset angle (thx LLMs)
    R = orientation(object)
    θ = acos(clamp((tr(R)-1)/2, -1, 1))
    if iszero(θ)
        return nothing
    end
    # Calculate rotation reset axis
    axis = 1/(2*sin(θ)) * [R[3,2]-R[2,3], R[1,3]-R[3,1], R[2,1]-R[1,2]]
    # Reset object rotation
    rotate3d!(object, axis, -θ)
    # Reset center of kinematics (removes precision artifacts)
    orientation!(object, Matrix{T}(I, 3, 3))
    return nothing
end