"""
    AbstractRay{T<:Real, S<:AbstractSegment{T}}

An implementation for a geometrical optics ray in R³. An `AbstractRay` is described by
its optical payload (wavelength `λ`, refractive index `n`, polarization state, etc.) and
wraps a spatial trajectory `segment::S` where `S <: AbstractSegment{T}`.

`AbstractRay`s model the propagation of light between optical interactions according to the laws of geometrical optics.
To store the solved sequence of ray propagation segments through an optical system, refer to [`AbstractBeam`](@ref).

# Implementation requirements
Subtypes of `AbstractRay` must implement:
- `segment`: an [`AbstractSegment`](@ref) describing the ray's spatial trajectory
- `λ`: wavelength in \\[m\\]
- `n`: refractive index along the ray path
"""
abstract type AbstractRay{T <: Real, S <: AbstractSegment{T}} end

segment(ray::AbstractRay) = ray.segment

Base.position(ray::AbstractRay) = position(segment(ray))
direction(ray::AbstractRay) = direction(segment(ray))
intersection(ray::AbstractRay) = intersection(segment(ray))
normal3d(ray::AbstractRay) = normal3d(segment(ray))
accumulated_opl(ray::AbstractRay) = accumulated_opl(segment(ray))

Base.length(ray::AbstractRay) = length(segment(ray))
geom_length(ray::AbstractRay) = length(segment(ray))
optical_path_length(ray::AbstractRay{T}) where {T} = isinf(length(ray)) ? T(Inf) : length(ray) * real(refractive_index(ray))

wavelength(ray::AbstractRay) = ray.λ
wavenumber(ray::AbstractRay) = 2π / wavelength(ray)
refractive_index(ray::AbstractRay) = ray.n

"""
    propagate(ray::AbstractRay, t)

Evaluates the 3D point along the ray at distance parameter `t`: ``\\vec{p}(t) = \\vec{p} + t \\cdot \\vec{d}``.
"""
@inline propagate(ray::AbstractRay, t::Real) = propagate(segment(ray), t)

"""
    with_segment(ray::AbstractRay, seg::AbstractSegment) -> AbstractRay

Replaces the trajectory segment of `ray` with `seg` while preserving optical properties.
"""
function with_segment end

"""
    with_accumulated_opl(ray::AbstractRay, opl::Real) -> AbstractRay

Returns a copy of `ray` with updated `accumulated_opl` on its segment.
"""
@inline with_accumulated_opl(ray::AbstractRay, opl::Real) = with_segment(ray, with_accumulated_opl(segment(ray), opl))


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
    angle3d(ray::AbstractRay)

Calculates the angle between a `ray` and its surface `intersection`.
"""
angle3d(ray::AbstractRay, intersect::AbstractIntersection) = angle3d(direction(ray), normal3d(intersect))
angle3d(ray::AbstractRay) = isnothing(intersection(ray)) ? nothing : angle3d(direction(ray), normal3d(intersection(ray)))

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
    isentering(ray)

Tests whether the ray is entering a shape based on the orientation of the `ray` direction and surface normal.
"""
isentering(r::AbstractRay, i::AbstractIntersection) = isentering(direction(r), normal3d(i))
isentering(d::AbstractArray, n::AbstractArray) = dot(d, n) < 0
isentering(::AbstractRay, ::Nothing) = false
isentering(r::AbstractRay) = isentering(r, intersection(r))


"""
    refraction3d(ray, int, n2)
    refraction3d(ray, n2)

Calculates the new direction of a `ray` entering into a new medium with ref. index `n2`.
"""
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
