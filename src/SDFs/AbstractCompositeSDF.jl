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
    _enclosing_sphere(a, b)

Returns the smallest sphere that encloses the spheres `a` and `b`, each a tuple `(center, radius)`
in the same frame, or `nothing` if one of them is `nothing`.
"""
function _enclosing_sphere(a::Tuple, b::Tuple)
    (c1, r1), (c2, r2) = a, b
    Δ = c2 - c1
    d = norm(Δ)
    # one sphere contains the other
    d + r2 ≤ r1 && return (c1, r1)
    d + r1 ≤ r2 && return (c2, r2)
    r = (d + r1 + r2) / 2
    return (c1 + (r - r1) / d * Δ, r)
end
_enclosing_sphere(::Nothing, ::Tuple) = nothing
_enclosing_sphere(::Tuple, ::Nothing) = nothing
_enclosing_sphere(::Nothing, ::Nothing) = nothing

"""
    _enclosing_sphere(shapes::Tuple)

Returns the sphere that encloses the [`world_bounding_sphere`](@ref) of every entry of `shapes`, in
world coordinates, or `nothing` if one of the `shapes` has none.
"""
_enclosing_sphere(shapes::Tuple{Any}) = world_bounding_sphere(shapes[1])
function _enclosing_sphere(shapes::Tuple{Any, Any, Vararg})
    return _enclosing_sphere(world_bounding_sphere(shapes[1]), _enclosing_sphere(Base.tail(shapes)))
end

"""
    _local_sphere(c::AbstractSDF, sphere)

Converts a `sphere` given in world coordinates into the local frame of `c`, as returned by
[`bounding_sphere`](@ref). Passes on `nothing`.
"""
function _local_sphere(c::AbstractSDF, sphere::Tuple)
    center, r = sphere
    return (transposed_orientation(c) * (center - position(c)), r)
end
_local_sphere(::AbstractSDF, ::Nothing) = nothing

"""
    kinematic_trait_of(c::AbstractCompositeSDF)

A composite takes the kinematic class of its [`operands`](@ref), which the constructor ensures
to be either all static or all movable, see [`BeamletOptics.AbstractKinematicTrait`](@ref).
"""
kinematic_trait_of(c::AbstractCompositeSDF) = _container_trait(operands(c))

"""
    translate3d!(::Movable, c::AbstractCompositeSDF, offset)

Translates `c` and all of its [`operands`](@ref) by `offset`.
"""
function translate3d!(::Movable, c::AbstractCompositeSDF, offset)
    position!(c, position(c) .+ offset)
    for s in operands(c)
        translate3d!(s, offset)
    end
    return nothing
end

"""
    rotate3d!(::Movable, c::AbstractCompositeSDF, R::AbstractMatrix)

Rotates `c` and all of its [`operands`](@ref) around `c`'s own origin (pivot), by the
rotation matrix `R`.
"""
function rotate3d!(::Movable, c::AbstractCompositeSDF, R::AbstractMatrix)
    # Update group orientation
    orientation!(c, R * orientation(c))
    # Rotate all operands around the composite center
    pivot = position(c)
    for s in operands(c)
        rotate3d!(s, R, pivot)
    end
    return nothing
end
