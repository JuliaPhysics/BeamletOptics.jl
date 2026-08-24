"""
    AbstractSegment{T<:Real}

Abstract supertype representing the spatial trajectory of a ray between interactions.

# Implementation requirements
Subtypes `S <: AbstractSegment{T}` must implement:
- `position(seg::S)`: start position of the segment as `Point3{T}`
- `direction(seg::S)`: normalized propagation direction as `Point3{T}`
- `intersection(seg::S)`: the boundary intersection or `nothing`
- `normal3d(seg::S)`: surface normal at the end of the segment or `nothing`
- `accumulated_opl(seg::S)`: cumulative optical path length up to this segment
- `Base.length(seg::S)`: geometrical length of the segment (or `Inf` if unhit)
- `propagate(seg::S, t)`: evaluate position along the segment at distance `t`
"""
abstract type AbstractSegment{T <: Real} end

position(s::AbstractSegment) = s.pos
direction(s::AbstractSegment) = s.dir
accumulated_opl(s::AbstractSegment) = s.accumulated_opl

"""
    OpenSegment{T<:Real} <: AbstractSegment{T}

Stores an unbounded / semi-infinite ray trajectory starting at `pos` in direction `dir`.
Used for initial rays before tracing and rays escaping to infinity without hitting an optical surface.

# Fields
- `pos::Point3{T}`: ray origin
- `dir::Point3{T}`: normalized propagation direction
- `accumulated_opl::T`: cumulative optical path length up to this segment
"""
struct OpenSegment{T <: Real} <: AbstractSegment{T}
    pos::Point3{T}
    dir::Point3{T}
    accumulated_opl::T
end

function OpenSegment(pos::AbstractArray{P}, dir::AbstractArray{D}, opl::Real = 0.0) where {P <: Real, D <: Real}
    T = promote_type(P, D, typeof(opl))
    if isapprox(norm(dir), 0, atol = 1e-14)
        throw(ErrorException("Direction vector too short for normalization."))
    end
    return OpenSegment{T}(Point3{T}(pos), normalize(Point3{T}(dir)), T(opl))
end

function OpenSegment(pos::Point3{T}, dir::Point3{T}, opl::Real = 0.0) where {T <: Real}
    if isapprox(norm(dir), 0, atol = 1e-14)
        throw(ErrorException("Direction vector too short for normalization."))
    end
    return OpenSegment{T}(pos, normalize(dir), T(opl))
end

position(s::OpenSegment) = s.pos
direction(s::OpenSegment) = s.dir
intersection(::OpenSegment) = nothing
normal3d(::OpenSegment) = nothing
accumulated_opl(s::OpenSegment) = s.accumulated_opl
Base.length(::OpenSegment{T}) where {T} = T(Inf)
@inline propagate(s::OpenSegment{T}, t::Real) where {T} = position(s) + T(t) * direction(s)


"""
    LineSegment{T<:Real} <: AbstractSegment{T}

Stores a bounded straight line segment starting at `p0` in direction `dir` and ending at `intersection`.

# Fields
- `p0::Point3{T}`: starting position of the segment
- `dir::Point3{T}`: normalized propagation direction
- `intersection::Intersection{T}`: surface intersection at the end of the segment
- `accumulated_opl::T`: cumulative optical path length up to the end of this segment
"""
struct LineSegment{T <: Real} <: AbstractSegment{T}
    p0::Point3{T}
    dir::Point3{T}
    intersection::Intersection{T}
    accumulated_opl::T
end

function LineSegment(p0::AbstractArray{P}, dir::AbstractArray{D}, int::Intersection{I}, opl::Real = 0.0) where {P <: Real, D <: Real, I <: Real}
    T = promote_type(P, D, I, typeof(opl))
    if isapprox(norm(dir), 0, atol = 1e-14)
        throw(ErrorException("Direction vector too short for normalization."))
    end
    return LineSegment{T}(Point3{T}(p0), normalize(Point3{T}(dir)), Intersection{T}(int.t, Point3{T}(int.p), Point3{T}(int.n)), T(opl))
end

function LineSegment(p0::Point3{T}, dir::Point3{T}, int::Intersection{T}, opl::Real = 0.0) where {T <: Real}
    if isapprox(norm(dir), 0, atol = 1e-14)
        throw(ErrorException("Direction vector too short for normalization."))
    end
    return LineSegment{T}(p0, normalize(dir), int, T(opl))
end

position(s::LineSegment) = s.p0
direction(s::LineSegment) = s.dir
intersection(s::LineSegment) = s.intersection
normal3d(s::LineSegment) = normal3d(s.intersection)
accumulated_opl(s::LineSegment) = s.accumulated_opl
Base.length(s::LineSegment) = length(s.intersection)
@inline propagate(s::LineSegment{T}, t::Real) where {T} = position(s) + T(t) * direction(s)
