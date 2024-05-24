"""
    Intersection{T}

Stores data calculated by the [`intersect3d`](@ref) method. This information can be reused, i.e. for retracing.

# Fields:

- `shape`: a reference to the [`AbstractShape`](@ref) that has been hit
- `object`: a nullable reference to the [`AbstractObject`](@ref) that has been hit (optional but recommended)
- `t`: length of the ray parametrization in [m]
- `n`: normal vector at the point of intersection
"""
mutable struct Intersection
    object::Nullable{AbstractObject}
    shape::Nullable{AbstractShape}
    t::Real
    n::Point3{Real}
end

function Intersection(t::T, n::AbstractArray{T}) where {T}
    return Intersection(nothing, nothing, t, Point3{T}(n))
end

function Intersection(t::T, n::AbstractArray{T}, shape::Nullable{AbstractShape{T}}) where {T}
    return Intersection(nothing, shape, t, Point3{T}(n))
end

shape(i::Intersection) = i.shape
object(i::Intersection) = i.object
object!(i::Intersection, new::AbstractObject) = (i.object = new) 

Base.length(i::Intersection) = i.t

normal3d(i::Intersection) = i.n

function Base.show(io::IO, ::MIME"text/plain", _intersection::Intersection)
    println(io, "Intersected object: $(typeof(object(_intersection)))")
    println(io, "Intersected shape: $(typeof(shape(_intersection)))")
    println(io, "Normal vector at intersection: $(normal3d(_intersection))")
    println(io, "Length to intersection: $(length(_intersection))")
end

"""
    AbstractRay{T<:Real}

An implementation for a geometrical optics ray in R³. In general, a `AbstractRay` is described by ``\\vec{p} + t\\cdot\\vec{d}`` with ``t\\in(0,\\infty)``. 
`AbstractRay`s are intended to model the propagation of light between optical interactions according to the laws of geometrical optics.
To store the result of a ray tracing solution, refer to [`AbstractBeam`](@ref).

# Intersections:

Since the length of a ray can not be known before solving an optical system, the [`Intersection`](@ref)-type is used.
This [`Nullable`](@ref) type can represent the intersection with an optical element, or lack thereof.

# Implementation reqs.

Subtypes of `AbstractBeam` must implement the following:

## Fields:

- `pos`: a R³-vector that stores the current position ``\\vec{p}``
- `dir`: a R³-vector that stores the current direction ``\\vec{d}``
- `intersection`: a `Nullable` field that stores the current [`Intersection`] or `nothing`
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

position(ray::AbstractRay) = ray.pos
position!(ray::AbstractRay, pos) = (ray.pos = pos)

direction(ray::AbstractRay) = ray.dir
function direction!(ray::AbstractRay, dir)
    ray.dir = normalize(dir)
    return nothing
end

wavelength(ray::AbstractRay) = ray.λ
wavelength!(ray::AbstractRay, λ) = (ray.λ = λ)

refractive_index(ray::AbstractRay) = ray.n
refractive_index!(ray::AbstractRay, n) = (ray.n = n)

intersection(ray::AbstractRay) = ray.intersection
function intersection!(ray::AbstractRay, _intersection::Nullable{Intersection})
     ray.intersection = _intersection
     return nothing
end

"""
    intersect3d(shape::AbstractShape, ::AbstractRay)

Defines the intersection between an [`AbstractShape`](@ref) and an [`AbstractRay`], must return an [`Intersection`](@ref) or `nothing`.
The default behavior for concrete shapes and rays is to indicate no intersection, that is `nothing`, which will inform the tracing algorithm to stop.
Refer to the [`Intersection`](@ref) typedocs for more information on the return type value.
"""
function intersect3d(shape::AbstractShape, ::AbstractRay)
    @warn lazy"No intersect3d method defined for:" typeof(shape)
    return nothing
end

"""
    intersect3d(object::AbstractObject, ray::AbstractRay)

In general, the intersection between an `object` and a `ray` is defined as the intersection with `shape(object)`.
"""
function intersect3d(object::AbstractObject, ray::AbstractRay)
    # FIXME: insinfrontof check?
    intersection = intersect3d(shape(object), ray)
    # Ensure that the intersection knows about the object if intersected
    if !isnothing(intersection)
        object!(intersection, object)
    end
    return intersection
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
        return Intersection(T(t), T.(plane_normal))
    end
end

"""
    length(ray::AbstractRay; opl::Bool=false)

Returns the geometric length of a `ray` between its start and intersection point. If no intersection exists, `Inf` is returned.
The `opl` keyword can be used to calculate the `o`ptical `p`ath `l`ength instead, i.e. ``OPL = n \\cdot l``.
Default is the geometrical length.
"""
function Base.length(ray::AbstractRay{T}; opl::Bool=false)::T where {T}
    if isnothing(intersection(ray))
        return T(Inf)
    end
    if opl        
        return length(intersection(ray)) * refractive_index(ray)
    else
        return length(intersection(ray))
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
    angle3d(ray::AbstractRay, intersect::Intersection=intersection(ray))

Calculates the angle between a `ray` and its or some other `intersection`.
"""
function angle3d(ray::AbstractRay, intersect::Intersection = intersection(ray))
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