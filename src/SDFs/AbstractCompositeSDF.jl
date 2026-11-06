"""
    AbstractCompositeSDF{T} <: AbstractSDF{T}

Generic supertype for boolean composite SDFs, i.e. SDFs that combine two or more
[`AbstractSDF`](@ref) operands into a single shape. Owns the field layout and
pivot-aware kinematics shared by all boolean composites; concrete subtypes only
need to define the boolean combination semantics (`sdf`, `normal3d`, `thickness`).

# Implementation reqs.

Subtypes of `AbstractCompositeSDF` should implement all reqs. of `AbstractSDF` as well as
the following:

## Fields

- `dir::SMatrix{3, 3, T, 9}`: the composite's own orientation matrix
- `transposed_dir::SMatrix{3, 3, T, 9}`: transpose of `dir`
- `pos::Point3{T}`: the composite's own position, in that field order, and the struct
  must be `mutable` since every setter in `AbstractShape.jl`/`AbstractSDF.jl` assigns
  fields directly

## Getters/setters

The [`AbstractShape`](@ref) getters and setters read these fields, so the composite needs no getters of its own.

## Functions

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
operands(c::AbstractCompositeSDF) = error(lazy"operands not defined for $(typeof(c))")

"""
    kinematic_trait_of(c::AbstractCompositeSDF)

A composite takes the kinematic class of its [`operands`](@ref), which the constructor ensures
to be either all static or all movable, see [`BeamletOptics.AbstractKinematicTrait`](@ref).
"""
kinematic_trait_of(c::AbstractCompositeSDF) = _container_trait(operands(c))

# Claim composite and its operands
function _claim_shapes!(owner::IdDict, c::AbstractCompositeSDF, part)
    invoke(_claim_shapes!, Tuple{IdDict, AbstractShape, Any}, owner, c, part)
    foreach(s -> _claim_shapes!(owner, s, part), operands(c))
    return nothing
end

"""
    translate3d!(::Movable, c::AbstractCompositeSDF, offset)

Translates `c` and all of its [`operands`](@ref) by `offset`, an operand listed several times once.
"""
function translate3d!(::Movable, c::AbstractCompositeSDF, offset)
    _move_rigidly!((c,), _movable_parts(operands(c)), offset)
    return nothing
end

"""
    rotate3d!(::Movable, c::AbstractCompositeSDF, R::AbstractMatrix)

Rotates `c` and all of its [`operands`](@ref) around `c`'s own origin (pivot), by the
rotation matrix `R`, an operand listed several times once.
"""
function rotate3d!(::Movable, c::AbstractCompositeSDF, R::AbstractMatrix)
    _move_rigidly!((c,), _movable_parts(operands(c)), R, position(c))
    return nothing
end
