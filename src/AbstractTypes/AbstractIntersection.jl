"""
    AbstractIntersection{T}

Stores data calculated by the [`intersect3d`](@ref) method.

# Required fields

- `t`: length of the ray parametrization in [m]
- `n`: normal vector at the point of intersection

# Interface

Subtypes must have the `t` and `n` fields above, which back the generic [`Base.length`](@ref)/[`normal3d`](@ref)
accessors defined for `AbstractIntersection`. Where applicable, subtypes should also implement `shape`
and/or `object`, see [`ShapeIntersection`](@ref) and [`ObjectIntersection`](@ref).
"""
abstract type AbstractIntersection{T} end

Base.length(i::AbstractIntersection) = i.t

normal3d(i::AbstractIntersection) = i.n

"""
    ShapeIntersection{T} <: AbstractIntersection

Stores the result of a shape-level [`intersect3d`](@ref) call, i.e. before the enclosing [`AbstractObject`](@ref)
is known.

# Fields

- `shape`: a reference to the [`AbstractShape`](@ref) of the `shape` that has been hit
"""
struct ShapeIntersection{T <: Real, S <: AbstractShape{T}} <: AbstractIntersection{T}
    shape::S
    t::T
    n::Point3{T}
end

function ShapeIntersection(t::Real, n::AbstractArray{<:Real}, shape::AbstractShape{S}) where {S <: Real}
    return ShapeIntersection(shape, S(t), Point3{S}(n))
end

shape(i::ShapeIntersection) = i.shape

function Base.show(io::IO, ::MIME"text/plain", i::ShapeIntersection)
    println(io, "Intersected shape: $(typeof(shape(i)))")
    println(io, "Normal vector at intersection: $(normal3d(i))")
    println(io, "Length to intersection: $(length(i))")
end

"""
    ObjectIntersection{T} <: AbstractIntersection

Stores the result of an object-level [`intersect3d`](@ref) call, i.e. once the enclosing [`AbstractObject`](@ref)
has been attached to a [`ShapeIntersection`](@ref).

# Fields

- `object`: a reference to the [`AbstractObject`](@ref) of the `object` that has been hit
- `hit`: a reference to the underlying [`ShapeIntersection`](@ref)
"""
mutable struct ObjectIntersection{T <: Real, S <: AbstractShape{T}, O <: AbstractObject{T}} <: AbstractIntersection{T}
    object::O
    hit::ShapeIntersection{T, S}
end

shape(i::ObjectIntersection) = shape(i.hit)
object(i::ObjectIntersection) = i.object

Base.length(i::ObjectIntersection) = length(i.hit)

normal3d(i::ObjectIntersection) = normal3d(i.hit)

function Base.show(io::IO, ::MIME"text/plain", i::ObjectIntersection)
    println(io, "Intersected object: $(typeof(object(i)))")
    println(io, "Intersected shape: $(typeof(shape(i)))")
    println(io, "Normal vector at intersection: $(normal3d(i))")
    println(io, "Length to intersection: $(length(i))")
end