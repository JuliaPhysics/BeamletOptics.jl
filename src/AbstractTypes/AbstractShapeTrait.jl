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

- `shape(::Foo)`: a getter function that returns a `Tuple` of all relevant parts, e.g. `(foo.front, foo.middle, foo.back)`.
  Parts may be [`AbstractShape`](@ref)s or nested [`AbstractObject`](@ref)s (e.g. a coating modeled as its own object) —
  `intersect3d`/kinematic functions dispatch on whichever is returned.

# Additional information

!!! warning "Kinematic center"
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

@inline _in_tuple(x, ::Tuple{}) = false
@inline _in_tuple(x, t::Tuple) = (x === first(t)) || _in_tuple(x, Base.tail(t))

"""
    unique_shapes(object::Union{AbstractObject, AbstractShape}) -> Tuple

Collects all unique underlying leaf shapes within `object` without duplicate references (stack-allocated, zero-allocation Tuple).
"""
@inline unique_shapes(object::Union{AbstractObject, AbstractShape}) = _collect_unique_shapes(object, ())

@inline function _collect_unique_shapes(item::AbstractShape, visited::Tuple)
    if _in_tuple(item, visited)
        return visited
    else
        return (visited..., item)
    end
end

@inline function _collect_unique_shapes(item::AbstractObject, visited::Tuple)
    if shape_trait_of(item) isa SingleShape
        return _collect_unique_shapes(shape(item), visited)
    elseif shape_trait_of(item) isa MultiShape
        return _collect_unique_shapes_tuple(shape(item), visited)
    else
        return visited
    end
end

@inline _collect_unique_shapes_tuple(::Tuple{}, visited::Tuple) = visited
@inline function _collect_unique_shapes_tuple(parts::Tuple, visited::Tuple)
    head = first(parts)
    tail = Base.tail(parts)
    v_next = _collect_unique_shapes(head, visited)
    return _collect_unique_shapes_tuple(tail, v_next)
end

@inline function _collect_unique_shapes_tuple(parts::AbstractArray, visited::Tuple)
    v = visited
    for part in parts
        v = _collect_unique_shapes(part, v)
    end
    return v
end

"""
    translate3d!(::MultiShape, object, offset)

Moves all parts of the [`MultiShape`](@ref) `object` along the specified `offset` vector (zero allocations).
Also backs [`AbstractObjectGroup`](@ref)'s `translate3d!`, whose `shape(group) = objects(group)`
lets it reuse this same implementation.
"""
@inline function translate3d!(::MultiShape, object::Union{AbstractObject, AbstractObjectGroup}, offset)
    position!(object, position(object) .+ offset)
    _translate3d_multishape_tuple(shape(object), offset, (object,))
    return nothing
end

@inline _translate3d_multishape_tuple(::Tuple{}, offset, visited::Tuple) = visited
@inline function _translate3d_multishape_tuple(parts::Tuple, offset, visited::Tuple)
    head = first(parts)
    tail = Base.tail(parts)
    v_next = _translate3d_part!(head, offset, visited)
    return _translate3d_multishape_tuple(tail, offset, v_next)
end

@inline function _translate3d_multishape_tuple(parts::AbstractArray, offset, visited::Tuple)
    v = visited
    for part in parts
        v = _translate3d_part!(part, offset, v)
    end
    return v
end

@inline function _translate3d_part!(shape_part::AbstractShape, offset, visited::Tuple)
    if !_in_tuple(shape_part, visited)
        translate3d!(shape_part, offset)
        return (visited..., shape_part)
    end
    return visited
end

@inline function _translate3d_part!(obj_part::AbstractObject, offset, visited::Tuple)
    if shape_trait_of(obj_part) isa SingleShape
        s = shape(obj_part)
        if !_in_tuple(s, visited)
            translate3d!(s, offset)
            return (visited..., obj_part, s)
        else
            return (visited..., obj_part)
        end
    elseif shape_trait_of(obj_part) isa MultiShape
        if !_in_tuple(obj_part, visited)
            position!(obj_part, position(obj_part) .+ offset)
            v1 = (visited..., obj_part)
            return _translate3d_multishape_tuple(shape(obj_part), offset, v1)
        end
    else
        translate3d!(obj_part, offset)
    end
    return visited
end

@inline function _translate3d_part!(part::Any, offset, visited::Tuple)
    if !_in_tuple(part, visited)
        translate3d!(part, offset)
        return (visited..., part)
    end
    return visited
end

"""
    translate_to3d!(::MultiShape, object, target)

Translates all parts of the [`MultiShape`](@ref) `object` in parallel to the specified `target` position.
The `object` center point will be equal to the `target`.
"""
@inline function translate_to3d!(::MultiShape, object::Union{AbstractObject, AbstractObjectGroup}, target)
    current = position(object)
    translate3d!(object, target .- current)
    return nothing
end

"""
    rotate3d!(::MultiShape, object, axis, θ)

All parts of the [`MultiShape`](@ref) `object` are rotated around the pivot center via the specified angle `θ` and `axis` (zero allocations).
"""
@inline function rotate3d!(::MultiShape, object::Union{AbstractObject, AbstractObjectGroup}, axis, θ)
    R = rotate3d(axis, θ)
    orig_pos = position(object)
    orientation!(object, R * orientation(object))
    _rotate3d_multishape_tuple(shape(object), axis, θ, R, orig_pos, (object,))
    return nothing
end

@inline _rotate3d_multishape_tuple(::Tuple{}, axis, θ, R, orig_pos, visited::Tuple) = visited
@inline function _rotate3d_multishape_tuple(parts::Tuple, axis, θ, R, orig_pos, visited::Tuple)
    head = first(parts)
    tail = Base.tail(parts)
    v_next = _rotate3d_part!(head, axis, θ, R, orig_pos, visited)
    return _rotate3d_multishape_tuple(tail, axis, θ, R, orig_pos, v_next)
end

@inline function _rotate3d_multishape_tuple(parts::AbstractArray, axis, θ, R, orig_pos, visited::Tuple)
    v = visited
    for part in parts
        v = _rotate3d_part!(part, axis, θ, R, orig_pos, v)
    end
    return v
end

@inline function _rotate3d_part!(shape_part::AbstractShape, axis, θ, R, orig_pos, visited::Tuple)
    if !_in_tuple(shape_part, visited)
        rotate3d!(shape_part, axis, θ)
        v = position(shape_part) .- orig_pos
        v_shifted = (R * v) - v
        translate3d!(shape_part, v_shifted)
        return (visited..., shape_part)
    end
    return visited
end

@inline function _rotate3d_part!(obj_part::AbstractObject, axis, θ, R, orig_pos, visited::Tuple)
    if shape_trait_of(obj_part) isa SingleShape
        s = shape(obj_part)
        if !_in_tuple(s, visited)
            rotate3d!(obj_part, axis, θ)
            v = position(obj_part) .- orig_pos
            v_shifted = (R * v) - v
            translate3d!(obj_part, v_shifted)
            return (visited..., obj_part, s)
        else
            return (visited..., obj_part)
        end
    elseif shape_trait_of(obj_part) isa MultiShape
        if !_in_tuple(obj_part, visited)
            orientation!(obj_part, R * orientation(obj_part))
            v1 = (visited..., obj_part)
            return _rotate3d_multishape_tuple(shape(obj_part), axis, θ, R, orig_pos, v1)
        end
    else
        rotate3d!(obj_part, axis, θ)
    end
    return visited
end

@inline function _rotate3d_part!(part::Any, axis, θ, R, orig_pos, visited::Tuple)
    if !_in_tuple(part, visited)
        rotate3d!(part, axis, θ)
        return (visited..., part)
    end
    return visited
end

"""
    align3d!(::MultiShape, object, target_vec)

Rotates all constituent parts of the [`MultiShape`](@ref) `object` around its kinematic origin
such that the primary local y-axis aligns with `target_vec`.
"""
@inline function align3d!(::MultiShape, object::Union{AbstractObject, AbstractObjectGroup}, target_vec)
    start_vec = normalize(orientation(object)[:, 2])
    tgt = normalize(target_vec)
    c = dot(start_vec, tgt)
    if c ≈ 1
        return nothing
    elseif c ≈ -1
        ax = normal3d(start_vec)
        rotate3d!(MultiShape(), object, ax, π)
    else
        ax = normalize(cross(start_vec, tgt))
        θ = acos(clamp(c, -1.0, 1.0))
        rotate3d!(MultiShape(), object, ax, θ)
    end
    return nothing
end

"""
    reset_translation3d!(::MultiShape, object)

Resets all applied translations of the `object`, i.e. moves the center back to the origin.

!!! info "Parts within parts"
    Sub-part relative translations are not reset!
"""
function reset_translation3d!(::MultiShape, object::Union{AbstractObject{T}, AbstractObjectGroup{T}}) where T
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
function reset_rotation3d!(::MultiShape, object::Union{AbstractObject{T}, AbstractObjectGroup{T}}) where T
    # Calculate rotation reset angle (thx LLMs)
    R = orientation(object)
    θ = acos(clamp((tr(R)-1)/2, -1, 1))
    if iszero(θ)
        return nothing
    end
    # Calculate rotation reset axis
    if θ > π - 1e-5
        # Since θ ≈ π, sin(θ) ≈ 0, R is symmetric: R + I = 2 * axis * axisᵀ
        # Find the column of R + I with the largest norm to extract the axis robustly
        M = R + I
        col_idx = argmax(vec(sum(abs2, M, dims=1)))
        axis = normalize(M[:, col_idx])
    else
        axis = 1/(2*sin(θ)) * [R[3,2]-R[2,3], R[1,3]-R[3,1], R[2,1]-R[1,2]]
    end
    # Reset object rotation
    rotate3d!(object, axis, -θ)
    # Reset center of kinematics (removes precision artifacts)
    orientation!(object, Matrix{T}(I, 3, 3))
    return nothing
end