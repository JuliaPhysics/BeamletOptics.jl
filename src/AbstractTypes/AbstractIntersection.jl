"""
    AbstractIntersection{T}

Stores data calculated by the [`intersect3d`](@ref) method.
"""
abstract type AbstractIntersection{T} end

Base.length(i::AbstractIntersection) = i.t
normal3d(i::AbstractIntersection) = i.n
Base.position(i::AbstractIntersection) = i.p

"""
    Intersection{T<:Real} <: AbstractIntersection{T}

Stores the geometric details of a hit interface:
- `t`: distance parameter along ray parametrization in [m]
- `p`: 3D position of the hit point
- `n`: outward unit normal vector at the point of intersection
"""
struct Intersection{T <: Real} <: AbstractIntersection{T}
    t::T
    p::Point3{T}
    n::Point3{T}
end

function Intersection(t::Real, p::AbstractArray{<:Real}, n::AbstractArray{<:Real})
    T = promote_type(typeof(t), eltype(p), eltype(n))
    return Intersection{T}(T(t), Point3{T}(p), Point3{T}(n))
end

function Base.show(io::IO, ::MIME"text/plain", i::Intersection)
    println(io, "Intersection distance: $(length(i))")
    println(io, "Intersection position: $(position(i))")
    println(io, "Normal vector: $(normal3d(i))")
end

"""
    MultiIntersection{T<:Real, O1<:Union{Nothing, AbstractObject}, O2<:Union{Nothing, AbstractObject}} <: AbstractIntersection{T}

Stores coincident intersections sharing the same boundary within tolerance.
"""
struct MultiIntersection{T <: Real, O1 <: Union{Nothing, AbstractObject}, O2 <: Union{Nothing, AbstractObject}} <: AbstractIntersection{T}
    hit::Intersection{T}
    exiting::O1
    entering::O2
end

function MultiIntersection(hit::Intersection{T};
        exiting::O1 = nothing,
        entering::O2 = nothing) where {T, O1 <: Union{Nothing, AbstractObject}, O2 <: Union{Nothing, AbstractObject}}
    return MultiIntersection{T, O1, O2}(hit, exiting, entering)
end

exiting(mi::MultiIntersection) = mi.exiting
entering(mi::MultiIntersection) = mi.entering
Base.length(mi::MultiIntersection) = length(mi.hit)
normal3d(mi::MultiIntersection) = normal3d(mi.hit)
Base.position(mi::MultiIntersection) = position(mi.hit)

function Base.show(io::IO, ::MIME"text/plain", mi::MultiIntersection)
    println(io, "Intersection distance: $(length(mi))")
    println(io, "Intersection position: $(position(mi))")
    println(io, "Exiting object: $(isnothing(exiting(mi)) ? nothing : typeof(exiting(mi)))")
    println(io, "Entering object: $(isnothing(entering(mi)) ? nothing : typeof(entering(mi)))")
    println(io, "Normal vector: $(normal3d(mi))")
end