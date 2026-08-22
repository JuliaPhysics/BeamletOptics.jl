"""
    AbstractRay{T<:Real}

An implementation for a geometrical optics ray in R³. In general, an `AbstractRay` is described by
an immutable position ``\\vec{p}``, normalized direction ``\\vec{d}``, wavelength ``\\lambda``, and refractive index ``n``.
`AbstractRay`s model the propagation of light between optical interactions according to the laws of geometrical optics.
To store the solved sequence of ray propagation segments through an optical system, refer to [`AbstractBeam`](@ref) and [`RaySegment`](@ref).

# Fields
Subtypes of `AbstractRay` must implement:
- `pos`: a point in R³ that stores the ray origin ``\\vec{p}``
- `dir`: a unit vector in R³ that stores the propagation direction ``\\vec{d}``
- `λ`: wavelength in [m]
- `n`: refractive index along the ray path
"""
abstract type AbstractRay{T <: Real} end

Base.position(ray::AbstractRay) = ray.pos
position!(ray::AbstractRay, pos) = (ray.pos = pos)

direction(ray::AbstractRay) = ray.dir
function direction!(ray::AbstractRay, dir)
    ray.dir = normalize(dir)
    return nothing
end

wavelength(ray::AbstractRay) = ray.λ
wavelength!(ray::AbstractRay, λ) = (ray.λ = λ)

wavenumber(ray::AbstractRay) = 2π / wavelength(ray)

refractive_index(ray::AbstractRay) = ray.n
refractive_index!(ray::AbstractRay, n) = (ray.n = n)



"""
    propagate(ray::AbstractRay, t)

Evaluates the 3D point along the ray at distance parameter `t`: ``\\vec{p}(t) = \\vec{p} + t \\cdot \\vec{d}``.
"""
@inline propagate(ray::AbstractRay{T}, t::Real) where {T} = position(ray) + T(t) * direction(ray)

"""
    RaySegment{T<:Real, R<:AbstractRay{T}}

Stores the trace state of a ray propagating through a medium until its intersection
with an optical boundary (or escape to infinity).

# Fields
- `ray::R`: the ray propagating along this segment
- `intersection::Nullable{AbstractIntersection{T}}`: the boundary intersection point, or `nothing` if unhit
- `hit_object::Nullable{AbstractObject}`: the object hit at the end of this segment, or `nothing`
- `accumulated_opl::T`: cumulative optical path length up to the end of this segment
"""
struct RaySegment{T <: Real, R <: AbstractRay{T}}
    ray::R
    intersection::Nullable{AbstractIntersection{T}}
    hit_object::Nullable{AbstractObject}
    accumulated_opl::T
end

function RaySegment(ray::R, int::Nullable{AbstractIntersection{T}} = nothing, obj::Nullable{AbstractObject} = nothing, opl::Real = 0.0) where {T <: Real, R <: AbstractRay{T}}
    return RaySegment{T, R}(ray, int, obj, T(opl))
end

ray(seg::RaySegment) = seg.ray
intersection(seg::RaySegment) = seg.intersection
hit_object(seg::RaySegment) = seg.hit_object
accumulated_opl(seg::RaySegment) = seg.accumulated_opl

Base.position(seg::RaySegment) = position(seg.ray)
direction(seg::RaySegment) = direction(seg.ray)
wavelength(seg::RaySegment) = wavelength(seg.ray)
refractive_index(seg::RaySegment) = refractive_index(seg.ray)
normal3d(seg::RaySegment) = isnothing(seg.intersection) ? nothing : normal3d(seg.intersection)

geom_length(seg::RaySegment{T}) where {T} = isnothing(seg.intersection) ? T(Inf) : length(seg.intersection)
optical_path_length(seg::RaySegment{T}) where {T} = isnothing(seg.intersection) ? T(Inf) : length(seg.intersection) * real(refractive_index(seg.ray))
Base.length(seg::RaySegment) = geom_length(seg)

isentering(seg::RaySegment) = isnothing(seg.intersection) ? false : isentering(direction(seg.ray), normal3d(seg.intersection))


"""
    bounding_sphere(obj)

Returns `(c_local, r)` where `c_local` is the center of the bounding sphere in local coordinates and `r` is its radius, or `nothing` if not implemented.
"""
bounding_sphere(::Any) = nothing

"""
    intersect3d(shape::AbstractShape, ::AbstractRay)

Defines the intersection between an [`AbstractShape`](@ref) and an [`AbstractRay`](@ref), must return an [`Intersection`](@ref) or `nothing`.
The default behavior for concrete `shape`s and rays is to indicate no intersection, that is `nothing`, which will inform the tracing algorithm to stop.
"""
function intersect3d(shape::AbstractShape, ::AbstractRay)
    @warn lazy"No intersect3d method defined for:" typeof(shape)
    return nothing
end

"""
    intersect3d(object::AbstractObject, ray::AbstractRay)

In general, the intersection logic between an [`AbstractObject`](@ref) and an [`AbstractRay`](@ref) depends on the [`AbstractShapeTrait`](@ref).
"""
intersect3d(object::AbstractObject, ray::AbstractRay) = intersect3d(
    shape_trait_of(object), object, ray)

function intersect3d(::SingleShape, object::AbstractObject, ray::AbstractRay)
    return intersect3d(shape(object), ray)
end

function intersect3d(::MultiShape, object::AbstractObject, ray::AbstractRay{R}) where {R}
    closest::Nullable{Intersection{R}} = nothing
    for part in shape(object)
        temp = intersect3d(part, ray)
        isnothing(temp) && continue
        if isnothing(closest) || length(temp) < length(closest)
            closest = temp
        end
    end
    return closest
end

"""
    intersect3d(plane_position, plane_normal, ray)

Returns the intersection between a `ray` and an infinitely large plane which is characterized by its `position` and `normal`.
"""
function intersect3d(plane_position::AbstractArray,
        plane_normal::AbstractArray,
        ray::AbstractRay{T}) where {T}
    t = line_plane_distance3d(plane_position, plane_normal, position(ray), direction(ray))
    if isnothing(t)
        return nothing
    else
        hit_pos = position(ray) + T(t) * direction(ray)
        return Intersection(T(t), hit_pos, Point3{T}(plane_normal))
    end
end

"""
    line_point_distance3d(ray, point)

Returns value for the shortest distance between the `ray` (extended to ∞) and `point`.
"""
line_point_distance3d(ray::AbstractRay, point) = line_point_distance3d(position(ray),
    direction(ray),
    point)

"""
    angle3d(ray::AbstractRay, intersect::AbstractIntersection)
    angle3d(seg::RaySegment)

Calculates the angle between a `ray` and its surface `intersection`.
"""
angle3d(ray::AbstractRay, intersect::AbstractIntersection) = angle3d(direction(ray), normal3d(intersect))
angle3d(seg::RaySegment) = isnothing(seg.intersection) ? nothing : angle3d(direction(seg.ray), normal3d(seg.intersection))

"""
    isinfrontof(shape::AbstractShape, ray::AbstractRay)

A simple test to check if a `shape` lies "in front of" a `ray`.
The forward direction is here defined as the ray `orientation`.
Only works well if `ray` is **outside** of the volume of `shape`.
Can be dispatched to return more accurate results for subtypes of `AbstractShape`.
"""
isinfrontof(shape::AbstractShape, ray::AbstractRay) = isinfrontof(position(shape),
    position(ray),
    direction(ray))

"""
    isentering(ray, int)
    isentering(seg::RaySegment)

Tests whether the ray is entering a shape based on the orientation of the `ray` direction and surface normal.
"""
isentering(r::AbstractRay, i::AbstractIntersection) = isentering(direction(r), normal3d(i))
isentering(d::AbstractArray, n::AbstractArray) = dot(d, n) < 0
isentering(::AbstractRay, ::Nothing) = false


"""
    refraction3d(ray, int, n2)
    refraction3d(seg::RaySegment, n2)

Calculates the new direction of a `ray` entering into a new medium with ref. index `n2`.
"""
function refraction3d(seg::RaySegment, n2)
    isnothing(seg.intersection) && error("Cannot calculate refraction for segment without intersection")
    dir = direction(seg.ray)
    nml = normal3d(seg.intersection)
    if !isentering(seg)
        nml *= -1
    end
    n1 = refractive_index(seg.ray)
    return refraction3d(dir, nml, n1, n2)
end

function refraction3d(ray::AbstractRay, int::AbstractIntersection, n2)
    dir = direction(ray)
    nml = normal3d(int)
    if !isentering(dir, nml)
        nml *= -1
    end
    n1 = refractive_index(ray)
    return refraction3d(dir, nml, n1, n2)
end

function refraction3d(ray::AbstractRay, n2)
    int = intersection(ray)
    isnothing(int) && error("Cannot calculate refraction for ray without intersection")
    return refraction3d(ray, int, n2)
end

