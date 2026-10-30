"""
    AbstractBoundingSphere

Supertype of the results of [`bounding_sphere_of`](@ref): a [`SingleBoundingSphere`](@ref) for a
shape, a [`MultiBoundingSphere`](@ref) for several shapes or objects, or [`NoBoundingSphere`](@ref)
for anything without one. Like a trait, `bounding_sphere_of(x)` selects by the type of `x` how `x`
is tested against a ray, see `intersect3d(sphere, x, ray)`; in contrast to a trait, a sphere also
holds data of `x`, its center and radius.

# Implementation reqs.

## Fields

- `pos`: the center in world coordinates [m]
- `radius`: the radius [m]

## Getters/setters

- `position(sphere)`: reads `pos`
- `radius(sphere)`: reads `radius`

## Functions

None required. The generic `intersect3d(sphere, x, ray)` and `bounding_box` use the getters.
[`NoBoundingSphere`](@ref) has no fields. It has methods of its own wherever no sphere is a valid
input, and its `position` and `radius` throw an `ArgumentError`.
"""
abstract type AbstractBoundingSphere end

Base.position(sphere::AbstractBoundingSphere) = sphere.pos
radius(sphere::AbstractBoundingSphere) = sphere.radius

"""
    NoBoundingSphere <: AbstractBoundingSphere

States that a shape or object has no bounding sphere. It is then tested against every ray via
[`intersect3d`](@ref), see [`bounding_sphere_of`](@ref).

It has no center and no radius: `position` and `radius` throw an `ArgumentError`. Code that needs
them handles `NoBoundingSphere` by a method of its own, or tests for it first.
"""
struct NoBoundingSphere <: AbstractBoundingSphere end

function no_sphere_error(what::AbstractString)
    throw(ArgumentError(lazy"a NoBoundingSphere has no $what, handle the case of no bounding sphere before the $what is used"))
end

Base.position(::NoBoundingSphere) = no_sphere_error("center")
radius(::NoBoundingSphere) = no_sphere_error("radius")

"""
    SingleBoundingSphere{T} <: AbstractBoundingSphere

A sphere that encloses one shape, see [`bounding_sphere_of`](@ref).

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
    MultiBoundingSphere{T} <: AbstractBoundingSphere

A sphere that encloses several bounding spheres, i.e. all parts of a [`MultiShape`](@ref) object or
all objects of an object group, see [`bounding_sphere_of`](@ref). A ray that misses it hits none of
the parts, hence a solve tests it before the parts. It is tested against a ray like a
[`SingleBoundingSphere`](@ref) and differs only in what it states: the sphere of several things.

# Fields

- `pos`: the center in world coordinates [m]
- `radius`: the radius [m]
"""
struct MultiBoundingSphere{T <: Real} <: AbstractBoundingSphere
    pos::Point3{T}
    radius::T
end

function MultiBoundingSphere(pos::AbstractVector{P}, radius::R) where {P <: Real, R <: Real}
    T = promote_type(P, R)
    return MultiBoundingSphere{T}(Point3{T}(pos), T(radius))
end

"""
    MultiBoundingSphere(a::AbstractBoundingSphere, b::AbstractBoundingSphere)
    MultiBoundingSphere(spheres::Tuple)
    MultiBoundingSphere(spheres::AbstractVector)

Returns the smallest sphere that encloses the bounding spheres `a` and `b`, or all `spheres`, e.g.
those of the parts of an object or of the objects of a group. The result is `NoBoundingSphere()` if
one of them is none, since then no sphere encloses all parts, and for no spheres at all.
"""
function MultiBoundingSphere(a::AbstractBoundingSphere, b::AbstractBoundingSphere)
    # Assert inferred types
    pa, pb = position(a)::Point3, position(b)::Point3
    ra, rb = radius(a)::Real, radius(b)::Real
    Δ = pb - pa
    d = norm(Δ)
    # One sphere contains the other
    d + rb ≤ ra && return MultiBoundingSphere(pa + zero(Δ), ra + zero(d))
    d + ra ≤ rb && return MultiBoundingSphere(pb + zero(Δ), rb + zero(d))
    r = (d + ra + rb) / 2
    return MultiBoundingSphere(pa + (r - ra) / d * Δ, r)
end
MultiBoundingSphere(::NoBoundingSphere, ::AbstractBoundingSphere) = NoBoundingSphere()
MultiBoundingSphere(::AbstractBoundingSphere, ::NoBoundingSphere) = NoBoundingSphere()
MultiBoundingSphere(::NoBoundingSphere, ::NoBoundingSphere) = NoBoundingSphere()

function MultiBoundingSphere(spheres::Tuple{AbstractBoundingSphere, Vararg{AbstractBoundingSphere}})
    # Closure, a type as `foldl` operator is not inferred
    return MultiBoundingSphere(foldl((a, b) -> MultiBoundingSphere(a, b), spheres))
end
MultiBoundingSphere(::Tuple{}) = NoBoundingSphere()
# Vector of spheres, element type not inferred
function MultiBoundingSphere(spheres::AbstractVector)
    isempty(spheres) && return NoBoundingSphere()
    return MultiBoundingSphere(foldl((a, b) -> MultiBoundingSphere(a, b), spheres))
end

"""
    SingleBoundingSphere(sphere::AbstractBoundingSphere)
    MultiBoundingSphere(sphere::AbstractBoundingSphere)

The `sphere` with the same center and radius as a sphere of the other kind, e.g. the sphere around
the parts of a composite shape as the sphere of this one shape. `NoBoundingSphere()` is passed on.
"""
SingleBoundingSphere(sphere::AbstractBoundingSphere) = SingleBoundingSphere(position(sphere), radius(sphere))
SingleBoundingSphere(sphere::NoBoundingSphere) = sphere
MultiBoundingSphere(sphere::AbstractBoundingSphere) = MultiBoundingSphere(position(sphere), radius(sphere))
MultiBoundingSphere(sphere::NoBoundingSphere) = sphere

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
function SingleBoundingSphere(shape::AbstractShape, sphere::AbstractBoundingSphere)
    return SingleBoundingSphere(shape, position(sphere), radius(sphere))
end
SingleBoundingSphere(::AbstractShape, sphere::NoBoundingSphere) = sphere

"""
    bounding_sphere_of(shape::AbstractShape)
    bounding_sphere_of(object::AbstractObject)
    bounding_sphere_of(group::AbstractObjectGroup)

Computes a sphere that encloses the `shape`, `object` or `group` in its current pose and returns it as a
[`SingleBoundingSphere`](@ref) or [`MultiBoundingSphere`](@ref), or returns
[`NoBoundingSphere`](@ref)`()` if there is none.

For a shape, the method is an optional part of the [`AbstractShape`](@ref) interface. With the
default `NoBoundingSphere()` every ray is tested via [`intersect3d`](@ref). With a sphere, a
solve tests a ray against the sphere first and skips the shape if it misses, which pays off for
shapes with a costly `intersect3d`. The sphere is computed once per [`solve_system!`](@ref) and
not stored in the shape, hence the method may take some time, e.g. look at all vertices of a mesh.

The sphere must enclose every point of the shape, for all its parameters: a hit outside of the
sphere is lost without a warning. It should also be tight, since a ray that hits the sphere is
tested as before. A safety margin is not needed, the solver adds its own tolerance.
Draw the result with [`render!`](@ref)`(ax, sphere)` to look at it.

```julia
# a cylinder of the radius `r` and the height `h`, with the origin at the center of its base
function BeamletOptics.bounding_sphere_of(c::MyCylinder)
    return BeamletOptics.SingleBoundingSphere(c, Point3(0, c.h / 2, 0), sqrt(c.r^2 + (c.h / 2)^2))
end
```

For an object, the result follows from its shapes, see [`AbstractShapeTrait`](@ref): the sphere of
the shape of a [`SingleShape`](@ref) object, and for a [`MultiShape`](@ref) object or an object
group the [`MultiBoundingSphere`](@ref) that encloses the spheres of all its parts, or
`NoBoundingSphere()` if one part has none.
"""
bounding_sphere_of(::AbstractShape) = NoBoundingSphere()

bounding_sphere_of(object::AbstractObject) = bounding_sphere_of(shape_trait_of(object), object)

bounding_sphere_of(::SingleShape, object::AbstractObject) = bounding_sphere_of(shape(object))

function bounding_sphere_of(::MultiShape, object::AbstractObject)
    return MultiBoundingSphere(map(bounding_sphere_of, shape(object)))
end

bounding_sphere_of(group::AbstractObjectGroup) = MultiBoundingSphere(map(bounding_sphere_of, objects(group)))

"""
    bounding_box(sphere::AbstractBoundingSphere)

Returns the limits `(xmin, xmax, ymin, ymax, zmin, zmax)` of the axis-aligned box around the `sphere`.
"""
function bounding_box(sphere::AbstractBoundingSphere)
    c, r = position(sphere), radius(sphere)
    return c[1] - r, c[1] + r, c[2] - r, c[2] + r, c[3] - r, c[3] + r
end

#=
Table of the bounding spheres of a solve
=#

"""
The table of the bounding spheres of a solve, keyed by shape, object or group (by identity).
The values have concrete types, such that the lookup per ray needs no dynamic dispatch.
"""
const BoundingSphereTable = IdDict{Any, Union{SingleBoundingSphere{Float64}, MultiBoundingSphere{Float64}}}

# Key of the table in the task-local storage
const bounding_spheres_key = :BeamletOptics_bounding_spheres

"""
    current_bounding_spheres() -> Union{Nothing, BoundingSphereTable}

Returns the table of the solve that is running in the current task, or `nothing` outside of a solve.
"""
function current_bounding_spheres()
    return get(task_local_storage(), bounding_spheres_key, nothing)::Union{Nothing, BoundingSphereTable}
end

"""
    bounding_sphere_of(table::BoundingSphereTable, x)
    bounding_sphere_of(::Nothing, x)

Returns the sphere that the `table` holds for `x`, or `NoBoundingSphere()` if there is none. Nothing
is computed.
"""
bounding_sphere_of(table::BoundingSphereTable, x) = get(table, x, NoBoundingSphere())
bounding_sphere_of(::Nothing, x) = NoBoundingSphere()

# The table holds the sphere of a single-shape object under its shape, see `_register!`
bounding_sphere_of(table::BoundingSphereTable, object::AbstractObject) =
    bounding_sphere_of(table, shape_trait_of(object), object)
bounding_sphere_of(table::BoundingSphereTable, ::SingleShape, object::AbstractObject) =
    get(table, shape(object), NoBoundingSphere())
bounding_sphere_of(table::BoundingSphereTable, ::MultiShape, object::AbstractObject) =
    get(table, object, NoBoundingSphere())

"""
    bounding_spheres(system::AbstractSystem) -> BoundingSphereTable

Computes the bounding spheres of all objects of the `system` in their current poses. The spheres are
stored as `Float64`, with the radius enlarged by the rounding error of the shape's number type.
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
    # The object decides, not its shape
    return _store!(table, shape(object), bounding_sphere_of(object))
end

function _register!(table::BoundingSphereTable, ::MultiShape, object::AbstractObject)
    # Reuse the spheres of the parts
    parts = map(part -> _register!(table, part), shape(object))
    return _store!(table, object, MultiBoundingSphere(parts))
end

function _register!(table::BoundingSphereTable, group::AbstractObjectGroup)
    parts = map(part -> _register!(table, part), objects(group))
    return _store!(table, group, MultiBoundingSphere(parts))
end

# Center and radius as the table holds them
function table_values(sphere::AbstractBoundingSphere)
    T = typeof(radius(sphere))
    pos = Point3{Float64}(position(sphere))
    r = Float64(radius(sphere))
    return pos, r + sqrt(eps(T)) * (r + norm(pos))
end

# Store the `sphere` under the `key`, non-float spheres are skipped
_store!(::BoundingSphereTable, key, ::AbstractBoundingSphere) = NoBoundingSphere()
function _store!(table::BoundingSphereTable, key, sphere::SingleBoundingSphere{<:AbstractFloat})
    return table[key] = SingleBoundingSphere(table_values(sphere)...)
end
function _store!(table::BoundingSphereTable, key, sphere::MultiBoundingSphere{<:AbstractFloat})
    return table[key] = MultiBoundingSphere(table_values(sphere)...)
end

"""
    with_bounding_spheres(f, system::AbstractSystem)
    with_bounding_spheres(f, table::BoundingSphereTable)

Calls `f()` with the bounding spheres of the `system`, or the given `table`, as the table of the
current task. A table that is already set is kept when a `system` is given.
"""
function with_bounding_spheres(f, system::AbstractSystem)
    isnothing(current_bounding_spheres()) || return f()
    return with_bounding_spheres(f, bounding_spheres(system))
end

with_bounding_spheres(f, table::BoundingSphereTable) = task_local_storage(f, bounding_spheres_key, table)

#=
Test of a ray against a bounding sphere
=#

"""
    sphere_chord(sphere::AbstractBoundingSphere, ray, margin = 0)

Returns `(b, s)` such that the line of the `ray` passes through the `sphere` between the path lengths
`b - s` and `b + s` in [m], or `nothing` if the line misses the sphere. The radius is enlarged by the
`margin` in [m] and by the rounding error of the test, such that a ray towards a point on the sphere
still passes through it.
"""
function sphere_chord(sphere::AbstractBoundingSphere, ray::AbstractRay, margin = 0)
    # Assert inferred types
    center, r = position(sphere)::Point3, radius(sphere)::Real
    dir = direction(ray)::Point3
    oc = center - position(ray)::Point3
    b = dot(oc, dir)
    # Distance of the center from the line, no cancellation
    perp = oc - b * dir
    q = dot(perp, perp)
    tol = sqrt(eps(float(typeof(r))))
    ρ = r + margin + tol * (r + norm(center) + norm(oc))
    q > ρ^2 && return nothing
    return b, sqrt(ρ^2 - q)
end

"""
    _sphere_exit(sphere::AbstractBoundingSphere, ray, margin = 0)

Returns the path length in [m] at which the `ray` leaves the `sphere`, or `nothing` if the `ray`
does not pass through it, i.e. if its line misses the sphere or the sphere lies behind its start.
The radius is enlarged by the `margin` in [m] and by the rounding error of the test, such that a ray
towards a point on the sphere still passes through it.
"""
function _sphere_exit(sphere::AbstractBoundingSphere, ray::AbstractRay, margin = 0)
    chord = sphere_chord(sphere, ray, margin)
    isnothing(chord) && return nothing
    b, s = chord
    t_out = b + s
    t_out < 0 && return nothing
    return t_out
end

"""
    beyond(sphere::AbstractBoundingSphere, ray, t)

Returns `true` if nothing within the `sphere` can be hit by the `ray` at a path length up to `t` in
[m]: the `ray` misses the `sphere`, the sphere lies behind its start, or the `ray` enters it only
beyond `t`. Without a sphere ([`NoBoundingSphere`](@ref)) the result is `false`. The solver skips
such parts once it has found a hit at `t`.
"""
beyond(::NoBoundingSphere, ::AbstractRay, t) = false

function beyond(sphere::AbstractBoundingSphere, ray::AbstractRay, t)
    chord = sphere_chord(sphere, ray)
    isnothing(chord) && return true
    b, s = chord
    return b + s < 0 || b - s > t
end

"""
    intersect3d(sphere::AbstractBoundingSphere, x::Union{AbstractShape, Tuple, AbstractVector, Leaves, LeafObjects}, ray::AbstractRay)

Returns the intersection between `x` and the `ray` like `intersect3d(x, ray)`, but tests the `ray`
against the bounding `sphere` of `x` first. `x` is a shape, or the parts of a [`MultiShape`](@ref)
object or of an object group, see `intersect3d(parts, ray)`:

- [`NoBoundingSphere`](@ref): `x` is tested via `intersect3d(x, ray)`
- the `ray` misses the [`SingleBoundingSphere`](@ref) or [`MultiBoundingSphere`](@ref), or the
  sphere lies behind it: returns
  `nothing` without a test of `x`, i.e. of any of the parts
- otherwise `x` is tested via `intersect3d(x, ray)`
"""
intersect3d(::NoBoundingSphere, shape::AbstractShape, ray::AbstractRay) = intersect3d(shape, ray)

function intersect3d(sphere::AbstractBoundingSphere, shape::AbstractShape, ray::AbstractRay)
    isnothing(_sphere_exit(sphere, ray)) && return nothing
    return intersect3d(shape, ray)
end

# Parts have methods of their own, separate from the shapes
intersect3d(::NoBoundingSphere, parts::Union{Tuple, AbstractVector, Leaves, LeafObjects}, ray::AbstractRay) = intersect3d(parts, ray)

function intersect3d(sphere::AbstractBoundingSphere, parts::Union{Tuple, AbstractVector, Leaves, LeafObjects}, ray::AbstractRay)
    isnothing(_sphere_exit(sphere, ray)) && return nothing
    return intersect3d(parts, ray)
end
