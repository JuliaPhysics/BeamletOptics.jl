"""
    AbstractRay{T<:Real}

An implementation for a geometrical optics ray in R³. In general, a `AbstractRay` is described by ``\\vec{p} + t\\cdot\\vec{d}`` with ``t\\in(0,\\infty)``.
`AbstractRay`s are intended to model the propagation of light between optical interactions according to the laws of geometrical optics.
To store the result of a ray tracing solution, refer to [`AbstractBeam`](@ref).

# Intersections:

Since the length of a ray can not be known before solving an optical system, the [`AbstractIntersection`](@ref)-type is used.
This [`Nullable`](@ref) type can represent the intersection with an optical element, or lack thereof.

# Implementation reqs.

Subtypes of `AbstractBeam` must implement the following:

## Fields:

- `pos`: a R³-vector that stores the current position ``\\vec{p}``
- `dir`: a R³-vector that stores the current direction ``\\vec{d}``
- `intersection`: a `Nullable` field that stores the current [`AbstractIntersection`](@ref) or `nothing`
- `λ`: wavelength in [m]
- `n`: refractive index along the ray path

# Additional information

!!! info "Ray length"
    Base.`length`: this function is used to return the length of the `AbstractRay` (if no intersection exists, the ray length is `Inf`).
    The `opl` keyword can be used to obtain the optical path length instead.

!!! warning "Ray direction"
    Many functions assume that the `dir`ection vector has unit length (i.e. ``|\\vec{p}| = 1``).
    Violating this assumption might lead to spurious results.
"""
abstract type AbstractRay{T <: Real} end

Base.position(ray::AbstractRay) = ray.pos
position!(ray::AbstractRay, pos) = (ray.pos = pos)

"""
    direction(ray::AbstractRay)

Returns the direction vector of the `ray`.
"""
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

intersection(ray::AbstractRay) = ray.intersection
function intersection!(ray::AbstractRay, _intersection::Nullable{<:AbstractIntersection})
    ray.intersection = _intersection
    return nothing
end

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
    # QUESTION returns only shape information upstream, Coating handling when?
    return intersect3d(shape(object), ray)
end

function intersect3d(::MultiShape, object::AbstractObject, ray::AbstractRay{R}) where {R}
    closest::Nullable{Intersection{R}} = nothing
    # QUESTION proper recursion handling (unwrapping shapes of sub-objects) tested?
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

function geom_length(ray::AbstractRay{T}) where {T}
    isnothing(intersection(ray)) && return T(Inf)
    return length(intersection(ray))
end

"""
    optical_path_length(ray::AbstractRay{T}) where {T}

Calculate the optical path length of the `ray`, i.e. ``\\mathrm{OPL} = n \\cdot l``.
"""
function optical_path_length(ray::AbstractRay{T}) where {T}
    isnothing(intersection(ray)) && return T(Inf)
    return length(intersection(ray)) * real(refractive_index(ray))
end

"""
    Base.length(ray::AbstractRay)

Returns the geometric length of a `ray` between its start and intersection point. If no intersection exists, `Inf` is returned.

!!! tip
    Use [`optical_path_length`](@ref) to get the optical path length instead.
"""
Base.length(ray::AbstractRay{T}) where {T} = geom_length(ray)

"""
    line_point_distance3d(ray, point)

Returns value for the shortest distance between the `ray` (extended to ∞) and `point`.
"""
line_point_distance3d(ray::AbstractRay, point) = line_point_distance3d(position(ray),
    direction(ray),
    point)

"""
    angle3d(ray::AbstractRay, intersect::AbstractIntersection=intersection(ray))

Calculates the angle between a `ray` and its or some other `intersection`.
"""
function angle3d(ray::AbstractRay, intersect::AbstractIntersection = intersection(ray))
    return angle3d(direction(ray), normal3d(intersect))
end

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
    isentering(ray)

Tests whether the ray is entering a shape based on the orientation of the `ray` direction and surface normal.
If no intersection is present, default behavior is to return `false`.
"""
isentering(r::AbstractRay) = isentering(r, intersection(r))
isentering(r::AbstractRay, i::AbstractIntersection) = isentering(direction(r), normal3d(i))
isentering(d::AbstractArray, n::AbstractArray) = dot(d, n) < 0
isentering(::AbstractRay, ::Nothing) = false

"""
    refraction3d(ray, n2)

Calculates the new direction of a `ray` entering into a new medium with ref. index `n2`.
"""
function refraction3d(ray::AbstractRay, n2)
    dir = direction(ray)
    nml = normal3d(intersection(ray))
    # if beam is leaving substrate, flip normal
    if !isentering(ray)
        nml *= -1
    end
    n1 = refractive_index(ray)
    return refraction3d(dir, nml, n1, n2)
end
