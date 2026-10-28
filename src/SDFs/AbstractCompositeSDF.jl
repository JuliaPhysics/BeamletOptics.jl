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
    kinematic_trait_of(c::AbstractCompositeSDF)

A composite takes the kinematic class of its [`operands`](@ref), which the constructor ensures
to be either all static or all movable, see [`BeamletOptics.AbstractKinematicTrait`](@ref).
"""
kinematic_trait_of(c::AbstractCompositeSDF) = _container_trait(operands(c))

"""
    translate3d!(::Movable, c::AbstractCompositeSDF, offset)

Translates `c` and all of its [`operands`](@ref) by `offset`, a shared operand once.
"""
function translate3d!(::Movable, c::AbstractCompositeSDF, offset)
    translate3d!(c, offset, IdSet{Any}())
    return nothing
end

"""
    rotate3d!(::Movable, c::AbstractCompositeSDF, R::AbstractMatrix)

Rotates `c` and all of its [`operands`](@ref) around `c`'s own origin (pivot), by the
rotation matrix `R`, a shared operand once.
"""
function rotate3d!(::Movable, c::AbstractCompositeSDF, R::AbstractMatrix)
    rotate3d!(c, R, position(c), IdSet{Any}())
    return nothing
end

function translate3d!(c::AbstractCompositeSDF, offset, carried::IdSet)
    _mark_carried!(carried, c) || return nothing
    position!(c, position(c) .+ offset)
    for s in operands(c)
        translate3d!(s, offset, carried)
    end
    return nothing
end

function rotate3d!(c::AbstractCompositeSDF, R::AbstractMatrix, pivot::AbstractVector, carried::IdSet)
    _mark_carried!(carried, c) || return nothing
    p = position(c)
    # Update group orientation and move the composite center about the pivot
    orientation!(c, R * orientation(c))
    position!(c, pivot + R * (p - pivot))
    # Rotate all operands around the pivot
    for s in operands(c)
        rotate3d!(s, R, pivot, carried)
    end
    return nothing
end
