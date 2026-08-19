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
    PlaneIntersection{T} <: AbstractIntersection{T}

Stores the result of a ray/plane intersection test performed by the free-standing
`intersect3d(plane_position, plane_normal, ray)` helper. Unlike [`ShapeIntersection`](@ref),
no concrete [`AbstractShape`](@ref) can be attached, since the plane is defined only by a
`position`/`normal` pair rather than backed by a real shape instance.

# Fields

- `t`: length of the ray parametrization in [m]
- `n`: normal vector at the point of intersection
"""
struct PlaneIntersection{T <: Real} <: AbstractIntersection{T}
    t::T
    n::Point3{T}
end

function PlaneIntersection(t::Real, n::AbstractArray{<:Real})
    T = promote_type(typeof(t), eltype(n))
    return PlaneIntersection{T}(T(t), Point3{T}(n))
end

function Base.show(io::IO, ::MIME"text/plain", i::PlaneIntersection)
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
struct ObjectIntersection{T <: Real, S <: AbstractShape{T}, O <: AbstractObject{T}} <: AbstractIntersection{T}
    object::O
    hit::ShapeIntersection{T, S}
end

"""
    ObjectIntersection(object, i::ObjectIntersection)

Re-tags an existing [`ObjectIntersection`](@ref) with a different `object`, keeping the
underlying [`ShapeIntersection`](@ref). Used e.g. by composite [`MultiShape`](@ref) objects
to attach themselves as the hit `object` once the closest sub-hit has been determined,
reproducing what the old `object!` mutation used to do.
"""
function ObjectIntersection(object::O, i::ObjectIntersection{T, S}) where {T, S, O <: AbstractObject{T}}
    return ObjectIntersection{T, S, O}(object, i.hit)
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

"""
    MultiIntersection{T} <: AbstractIntersection{T}

Stores a primary [`ObjectIntersection`](@ref) together with up to two coincident
`ObjectIntersection`s sharing (numerically, within [`Config.get_coincident_boundary_tolerance`](@ref))
the same boundary — representing a ray landing on a boundary shared by multiple objects.

# Fields

- `hit`: the primary [`ObjectIntersection`](@ref), i.e. the surface actually reported by `intersect3d` (e.g. a coating)
- `exiting`: [`Nullable`](@ref) reference to the coincident [`ObjectIntersection`](@ref) of the object being exited
- `entering`: [`Nullable`](@ref) reference to the coincident [`ObjectIntersection`](@ref) of the object being entered

# Additional information

Which fields are set depends on the physical case:

- Cemented doublet interface (no coating): `hit` is whichever lens surface `intersect3d`
  found; the other lens is recorded as `exiting` **or** `entering` (only one of the two).
- Coating sandwiched between two substrates: `hit` is the coating; `exiting` is the
  substrate being left, `entering` is the substrate being entered.

!!! info "Not yet consumed by the generic solver"
    This type is scaffolding for future work on coincident multi-surface boundaries.
    `intersect3d` and component `interact3d` implementations currently resolve such
    ambiguity via [`Hint`](@ref)-directed re-tracing and shape-identity comparison instead.
"""
struct MultiIntersection{T <: Real} <: AbstractIntersection{T}
    hit::ObjectIntersection{T}
    exiting::Nullable{ObjectIntersection{T}}
    entering::Nullable{ObjectIntersection{T}}
end

function MultiIntersection(hit::ObjectIntersection{T};
        exiting::Nullable{ObjectIntersection{T}} = nothing,
        entering::Nullable{ObjectIntersection{T}} = nothing) where {T}
    return MultiIntersection{T}(hit, exiting, entering)
end

exiting(mi::MultiIntersection) = mi.exiting
entering(mi::MultiIntersection) = mi.entering

Base.length(mi::MultiIntersection) = length(mi.hit)

normal3d(mi::MultiIntersection) = normal3d(mi.hit)

object(mi::MultiIntersection) = object(mi.hit)

shape(mi::MultiIntersection) = shape(mi.hit)

function Base.show(io::IO, ::MIME"text/plain", mi::MultiIntersection)
    println(io, "Intersected object: $(typeof(object(mi)))")
    println(io, "Intersected shape: $(typeof(shape(mi)))")
    println(io, "Exiting object: $(isnothing(exiting(mi)) ? nothing : typeof(object(exiting(mi))))")
    println(io, "Entering object: $(isnothing(entering(mi)) ? nothing : typeof(object(entering(mi))))")
    println(io, "Normal vector at intersection: $(normal3d(mi))")
    println(io, "Length to intersection: $(length(mi))")
end