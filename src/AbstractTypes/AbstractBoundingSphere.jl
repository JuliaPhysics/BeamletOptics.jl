using Base.ScopedValues: ScopedValue, with

"""
    AbstractBoundingSphere

Supertype of the results of [`bounding_sphere_of`](@ref): a [`SingleBoundingSphere`](@ref), or
[`NoBoundingSphere`](@ref) for anything without one. Like a trait, `bounding_sphere_of(x)` selects
by the type of `x` how `x` is tested against a ray, see `intersect3d(sphere, shape, ray)`; in
contrast to a trait, a `SingleBoundingSphere` also holds data of `x`, its center and radius.
"""
abstract type AbstractBoundingSphere end

"""
    NoBoundingSphere <: AbstractBoundingSphere

States that a shape or object has no bounding sphere. It is then tested against every ray via
[`intersect3d`](@ref), see [`bounding_sphere_of`](@ref).
"""
struct NoBoundingSphere <: AbstractBoundingSphere end

"""
    SingleBoundingSphere{T} <: AbstractBoundingSphere

A sphere that encloses a shape, an object or a group, see [`bounding_sphere_of`](@ref).

# Fields

- `pos`: the center in world coordinates [m]
- `radius`: the radius [m]
"""
struct SingleBoundingSphere{T <: Real} <: AbstractBoundingSphere
    pos::Point3{T}
    radius::T
end

function SingleBoundingSphere(pos::AbstractVector{P}, radius::R) where {P <: Real, R <: Real}
    T = promote_type(P, R)
    return SingleBoundingSphere{T}(Point3{T}(pos), T(radius))
end

"""
    SingleBoundingSphere(shape::AbstractShape, center, radius)

The bounding sphere of the `shape` from a `center` given in the local frame of the `shape`, i.e.
relative to its [`position`](@ref) and along the axes of its [`orientation`](@ref), and the
`radius` in [m]. The center is converted to world coordinates for the current pose of the `shape`.
This is the constructor for a method of [`bounding_sphere_of`](@ref).
"""
function SingleBoundingSphere(shape::AbstractShape, center::AbstractVector, radius::Real)
    return SingleBoundingSphere(position(shape) + orientation(shape) * center, radius)
end

"""
    SingleBoundingSphere(shape::AbstractShape, sphere::AbstractBoundingSphere)

The `sphere`, whose center is given in the local frame of the `shape`, in world coordinates, e.g.
the sphere around the parts of a shape that are positioned in its local frame. `NoBoundingSphere()`
is passed on.
"""
SingleBoundingSphere(shape::AbstractShape, sphere::SingleBoundingSphere) = SingleBoundingSphere(shape, sphere.pos, sphere.radius)
SingleBoundingSphere(::AbstractShape, sphere::NoBoundingSphere) = sphere

Base.position(sphere::SingleBoundingSphere) = sphere.pos

"""
    bounding_sphere_of(shape::AbstractShape)
    bounding_sphere_of(object::AbstractObject)

Computes a sphere that encloses the `shape` or `object` in its current pose and returns it as a
[`SingleBoundingSphere`](@ref), or returns [`NoBoundingSphere`](@ref)`()` if there is none.

For a shape, the method is an optional part of the [`AbstractShape`](@ref) interface. With the
default `NoBoundingSphere()` every ray is tested via [`intersect3d`](@ref). With a sphere, a
solve tests a ray against the sphere first and skips the shape if it misses, which pays off for
shapes with a costly `intersect3d`. The sphere is computed once per [`solve_system!`](@ref) and
not stored in the shape, hence the method may take some time, e.g. look at all vertices of a mesh.

The sphere must enclose every point of the shape, for all its parameters: a hit outside of the
sphere is lost without a warning. It should also be tight, since a ray that hits the sphere is
tested as before. A safety margin is not needed, the solver adds its own tolerance. Use
[`render_bounding_sphere!`](@ref) to look at the result.

```julia
# a cylinder of the radius `r` and the height `h`, with the origin at the center of its base
function BeamletOptics.bounding_sphere_of(c::MyCylinder)
    return BeamletOptics.SingleBoundingSphere(c, Point3(0, c.h / 2, 0), sqrt(c.r^2 + (c.h / 2)^2))
end
```

For an object, the result follows from its shapes, see [`AbstractShapeTrait`](@ref): the sphere of
the shape of a [`SingleShape`](@ref) object, and for a [`MultiShape`](@ref) object or an object
group the sphere that encloses the spheres of all its parts, or `NoBoundingSphere()` if one part has
none.

A solve does not call this method per ray. It computes the spheres once into a table and looks
them up via `bounding_sphere_of(table, x)`.
"""
bounding_sphere_of(::AbstractShape) = NoBoundingSphere()

bounding_sphere_of(object::AbstractObject) = bounding_sphere_of(shape_trait_of(object), object)

bounding_sphere_of(::SingleShape, object::AbstractObject) = bounding_sphere_of(shape(object))

function bounding_sphere_of(::MultiShape, object::AbstractObject)
    return _enclosing_sphere(map(bounding_sphere_of, shape(object)))
end

"""
    _enclosing_sphere(a, b)

Returns the smallest sphere that encloses the bounding spheres `a` and `b`, or
`NoBoundingSphere()` if one of them is none.
"""
function _enclosing_sphere(a::SingleBoundingSphere, b::SingleBoundingSphere)
    Δ = b.pos - a.pos
    d = norm(Δ)
    # one sphere contains the other
    d + b.radius ≤ a.radius && return SingleBoundingSphere(a.pos + zero(Δ), a.radius + zero(d))
    d + a.radius ≤ b.radius && return SingleBoundingSphere(b.pos + zero(Δ), b.radius + zero(d))
    r = (d + a.radius + b.radius) / 2
    return SingleBoundingSphere(a.pos + (r - a.radius) / d * Δ, r)
end
_enclosing_sphere(::AbstractBoundingSphere, ::AbstractBoundingSphere) = NoBoundingSphere()

# The sphere that encloses all `spheres`, none for an object without parts
_enclosing_sphere(spheres::Tuple{Vararg{AbstractBoundingSphere}}) = foldl(_enclosing_sphere, spheres)
_enclosing_sphere(::Tuple{}) = NoBoundingSphere()

"""
    bounding_box(sphere::SingleBoundingSphere)

Returns the limits `(xmin, xmax, ymin, ymax, zmin, zmax)` of the axis-aligned box around the `sphere`.
"""
function bounding_box(sphere::SingleBoundingSphere)
    c, r = sphere.pos, sphere.radius
    return c[1] - r, c[1] + r, c[2] - r, c[2] + r, c[3] - r, c[3] + r
end

#=
Table of the bounding spheres of a solve
=#

"""
The table of the bounding spheres of a solve, see [`bounding_spheres`](@ref): the sphere of a shape
under the shape, the sphere of a [`MultiShape`](@ref) object or group under the object. The keys are
compared by identity.
"""
const BoundingSphereTable = IdDict{Any, SingleBoundingSphere{Float64}}

"""
The [`BoundingSphereTable`](@ref) of the solve that is running in the current task, or `nothing`
outside of a solve. It is set by [`with_bounding_spheres`](@ref), and the intersection code passes
it to `bounding_sphere_of(table, x)`.
"""
const BOUNDING_SPHERES = ScopedValue{Union{Nothing, BoundingSphereTable}}(nothing)

"""
    bounding_sphere_of(table::BoundingSphereTable, x)
    bounding_sphere_of(::Nothing, x)

Returns the bounding sphere that the `table` holds for the shape or object `x`, or
[`NoBoundingSphere`](@ref)`()` if it holds none. Nothing is computed, in contrast to
`bounding_sphere_of(x)`, which fills the table, see [`bounding_spheres`](@ref).

The intersection code calls this method per ray with the table of the solve that is running in the
current task, `BOUNDING_SPHERES[]`. Outside of a solve this table is `nothing` and the result is
always `NoBoundingSphere()`, i.e. `intersect3d(object, ray)` then tests the object without a sphere.
"""
bounding_sphere_of(table::BoundingSphereTable, x) = get(table, x, NoBoundingSphere())
bounding_sphere_of(::Nothing, x) = NoBoundingSphere()

"""
    bounding_spheres(system::AbstractSystem) -> BoundingSphereTable

Computes the bounding spheres of all objects of the `system` in their current poses, via
[`bounding_sphere_of`](@ref): of the shape of every [`SingleShape`](@ref) object, and of every
[`MultiShape`](@ref) object and object group, which encloses the spheres of its parts. Shapes and
objects without a sphere have no entry.

The spheres are stored as `Float64`, with the radius enlarged by the rounding error of the number
type of the shape. A sphere of a number type that is no `AbstractFloat` is not stored, i.e. this
shape is tested without a sphere.
"""
function bounding_spheres(system::AbstractSystem)
    table = BoundingSphereTable()
    for object in _top_level(system)
        _register!(table, object)
    end
    return table
end

# Stores the sphere of `x` and of its parts in the `table`, returns the stored sphere of `x`
_register!(table::BoundingSphereTable, object::AbstractObject) = _register!(table, shape_trait_of(object), object)
_register!(table::BoundingSphereTable, shape::AbstractShape) = _store!(table, shape, bounding_sphere_of(shape))
# anything else that a custom system type holds has no sphere
_register!(::BoundingSphereTable, ::Any) = NoBoundingSphere()

function _register!(table::BoundingSphereTable, ::SingleShape, object::AbstractObject)
    # the object decides, e.g. a `NonInteractableObject` has no sphere whatever its shape
    return _store!(table, shape(object), bounding_sphere_of(object))
end

function _register!(table::BoundingSphereTable, ::MultiShape, object::AbstractObject)
    # the spheres of the parts are computed once and reused for the sphere around them
    parts = map(part -> _register!(table, part), shape(object))
    return _store!(table, object, _enclosing_sphere(parts))
end

_store!(::BoundingSphereTable, key, ::AbstractBoundingSphere) = NoBoundingSphere()
function _store!(table::BoundingSphereTable, key, sphere::SingleBoundingSphere{T}) where {T <: AbstractFloat}
    pos = Point3{Float64}(sphere.pos)
    r = Float64(sphere.radius)
    # rounding of the pose and the parameters of the shape in its own number type
    stored = SingleBoundingSphere(pos, r + sqrt(eps(T)) * (r + norm(pos)))
    table[key] = stored
    return stored
end

"""
    with_bounding_spheres(f, system::AbstractSystem)

Calls `f()` with the [`bounding_spheres`](@ref) of the `system` as the table of the current task
and of the tasks it spawns, such that the intersection code finds them. If a table is
already set, e.g. by the solve of the beam group that a beam belongs to, it is kept and nothing is
computed.

The `system` must not be changed while `f` runs: the spheres are those of the poses at the call.
"""
function with_bounding_spheres(f, system::AbstractSystem)
    isnothing(BOUNDING_SPHERES[]) || return f()
    return with(f, BOUNDING_SPHERES => bounding_spheres(system))
end

#=
Test of a ray against a bounding sphere
=#

"""
    _sphere_exit(sphere::SingleBoundingSphere, ray, margin = 0)

Returns the path length in [m] at which the `ray` leaves the `sphere`, or `nothing` if the `ray`
does not pass through it, i.e. if its line misses the sphere or the sphere lies behind its start.
The radius is enlarged by the `margin` in [m] and by the rounding error of the test, such that a ray
towards a point on the sphere still passes through it.
"""
function _sphere_exit(sphere::SingleBoundingSphere, ray::AbstractRay, margin = 0)
    center, radius = sphere.pos, sphere.radius
    dir = direction(ray)
    oc = center - position(ray)
    b = dot(oc, dir)
    # Distance of the center from the line of the ray. This form does not cancel for a far start,
    # in contrast to |oc|² - b².
    perp = oc - b * dir
    q = dot(perp, perp)
    tol = sqrt(eps(float(typeof(radius))))
    ρ = radius + margin + tol * (radius + norm(center) + norm(oc))
    q > ρ^2 && return nothing
    t_out = b + sqrt(ρ^2 - q)
    t_out < 0 && return nothing
    return t_out
end

# Whether the `ray` can not hit anything within the bounding sphere
_misses(::NoBoundingSphere, ::AbstractRay) = false
_misses(sphere::SingleBoundingSphere, ray::AbstractRay) = isnothing(_sphere_exit(sphere, ray))

"""
    intersect3d(sphere::AbstractBoundingSphere, shape::AbstractShape, ray::AbstractRay)

Returns the intersection between the `shape` and the `ray` like `intersect3d(shape, ray)`, but tests
the `ray` against the bounding `sphere` of the `shape` first:

- [`NoBoundingSphere`](@ref): the shape is tested via `intersect3d(shape, ray)`
- the `ray` misses the [`SingleBoundingSphere`](@ref), or the sphere lies behind it: returns
  `nothing` without a test of the `shape`
- otherwise the shape is tested via `intersect3d(shape, ray)`

The intersection code of the objects calls this method with the sphere that the table of the
running solve holds for the `shape`, `bounding_sphere_of(BOUNDING_SPHERES[], shape)`. A shape type
implements `intersect3d(shape, ray)` and, optionally, `bounding_sphere_of(shape)`, not this method.
"""
intersect3d(::NoBoundingSphere, shape::AbstractShape, ray::AbstractRay) = intersect3d(shape, ray)

function intersect3d(sphere::SingleBoundingSphere, shape::AbstractShape, ray::AbstractRay)
    _misses(sphere, ray) && return nothing
    return intersect3d(shape, ray)
end
