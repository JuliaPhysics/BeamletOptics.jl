"""
    AbstractCompositeSDF{T} <: AbstractSDF{T}

Generic supertype for boolean composite SDFs, i.e. SDFs that combine two or more
[`AbstractSDF`](@ref) operands into a single shape. Owns the field layout and
pivot-aware kinematics shared by all boolean composites; concrete subtypes only
need to define the boolean combination semantics (`sdf`, `normal3d`, `thickness`).

# Implementation reqs.

Subtypes of `AbstractCompositeSDF` should implement all reqs. of `AbstractSDF` as well as
the following:

# Fields

- `dir::SMatrix{3, 3, T, 9}`: the composite's own orientation matrix
- `transposed_dir::SMatrix{3, 3, T, 9}`: transpose of `dir`
- `pos::Point3{T}`: the composite's own position, in that field order, and the struct
  must be `mutable` since every setter in `AbstractShape.jl`/`AbstractSDF.jl` assigns
  fields directly

# Functions

- [`operands`](@ref): returns a `Tuple` of every child `AbstractSDF`, in any order
- `sdf`, `normal3d` and `thickness` specific to the boolean combination
"""
abstract type AbstractCompositeSDF{T} <: AbstractSDF{T} end

"""
    operands(c::AbstractCompositeSDF)

Returns a `Tuple` of every child [`AbstractSDF`](@ref) that the composite `c` combines, in
any order. This is the accessor the shared kinematics of [`AbstractCompositeSDF`](@ref)
iterate over, so every concrete composite must implement it; the order carries no meaning
and must not be relied upon to identify an operand's role in the boolean expression.
"""
function operands end

"""
    translate3d!(c::AbstractCompositeSDF, offset)

Translates `c` and all of its [`operands`](@ref) by `offset`.
"""
function translate3d!(c::AbstractCompositeSDF, offset)
    position!(c, position(c) .+ offset)
    for s in operands(c)
        translate3d!(s, offset)
    end
    return nothing
end

"""
    rotate3d!(c::AbstractCompositeSDF, R::AbstractMatrix)

Rotates `c` and all of its [`operands`](@ref) around `c`'s own origin (pivot), by the
rotation matrix `R`.
"""
function rotate3d!(c::AbstractCompositeSDF, R::AbstractMatrix)
    # Update group orientation
    orientation!(c, R * orientation(c))
    # Rotate all operands around the composite center
    for s in operands(c)
        rotate3d!(s, R)
        v = position(s) - position(c)
        # Translate group around pivot point
        v = (R * v) - v
        translate3d!(s, v)
    end
    return nothing
end

"""
    rotate3d!(c::AbstractCompositeSDF, axis, θ)

Rotates `c` and all of its [`operands`](@ref) around `c`'s own origin (pivot), by an
angle `θ` around `axis`.
"""
function rotate3d!(c::AbstractCompositeSDF, axis, θ)
    R = rotate3d(axis, θ)
    return rotate3d!(c, R)
end

"""
    align3d!(c::AbstractCompositeSDF, target_axis)

Rotates `c` and all of its [`operands`](@ref) around `c`'s own origin (pivot) such that
its local y-axis aligns with `target_axis`.
"""
function align3d!(c::AbstractCompositeSDF, target_axis)
    R = align3d(orientation(c)[:, 2], target_axis)
    rotate3d!(c, R)
    return nothing
end

"""
    reset_translation3d!(c::AbstractCompositeSDF)

Resets the translation of `c` and all of its [`operands`](@ref), returning the composite
origin to `(0, 0, 0)` while preserving relative operand positions.
"""
function reset_translation3d!(c::AbstractCompositeSDF{T}) where {T}
    translate3d!(c, -position(c))
    position!(c, Point3{T}(0))
    return nothing
end

"""
    reset_rotation3d!(c::AbstractCompositeSDF)

Resets the orientation of `c` and all of its [`operands`](@ref) back to the standard basis,
preserving relative operand positions and orientations.
"""
function reset_rotation3d!(c::AbstractCompositeSDF{T}) where {T}
    R = orientation(c)'
    rotate3d!(c, R)
    orientation!(c, Matrix{T}(I, 3, 3))
    return nothing
end
