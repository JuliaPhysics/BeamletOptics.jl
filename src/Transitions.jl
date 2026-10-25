#=
Transitions: what a ray crosses at a boundary

The solver collects all surfaces that a ray hits at one point, see `coincident_hits`. The bodies
that the ray leaves and enters there define the media on both sides of the boundary, and the thin
layers at it (coatings, films, beamsplitters) act between them, see `transition`. Hence touching
elements, e.g. cemented lenses or a coating on its substrate, need no knowledge of each other.
=#

"""
    find_part(x, s)

Returns the object within `x` whose shape is `s`, or `nothing`. `x` is an object, which may consist
of objects (a [`MultiShape`](@ref) object), or a shape, which holds no object.
"""
find_part(::AbstractShape, s) = nothing
find_part(x::AbstractObject, s) = find_part(shape_trait_of(x), x, s)
find_part(::SingleShape, x::AbstractObject, s) = shape(x) === s ? x : nothing

function find_part(::MultiShape, x::AbstractObject, s)
    for part in shape(x)
        found = find_part(part, s)
        isnothing(found) || return found
    end
    return nothing
end

"""
    part_of(object::AbstractObject, s)

Returns the part of the `object` whose shape is `s`: the object itself for a single shape or a
[`MultiShape`](@ref) object of shapes, else the object among its parts, e.g. the back lens of a
[`DoubletLens`](@ref) for the shape of the back lens. This part interacts with a ray that hits `s`.
"""
part_of(object::AbstractObject, s) = something(find_part(object, s), object)

"""Returns the part of the object that the hit `hit` belongs to, see [`part_of`](@ref)."""
part_of(hit::Intersection) = isnothing(hit.part) ? part_of(object(hit), shape(hit)) : hit.part

"""
    interaction_rank(hit, ray, medium)

Returns the rank of the `hit` of the `ray` among the [`coincident_hits`](@ref): the hit of the lowest
rank interacts with the `ray`. `medium` is the [`medium`](@ref) of the part that was hit. A surface
that is no body (`medium` is `nothing`) ranks first, e.g. a thin beamsplitter, a mirror or a detector
on a lens: it acts at the boundary, see [`Transition`](@ref). The surface of a body that the `ray`
enters ranks before the surface of a body that it leaves.
"""
interaction_rank(hit, ray, ::Nothing) = 0
interaction_rank(hit, ray, ::AbstractMedium) = isentering(direction(ray), normal3d(hit)) ? 1 : 2

"""
    first_to_interact(hit, ray)

Returns the chain of [`coincident_hits`](@ref) `hit` of the `ray` with the hit that interacts first,
see [`interaction_rank`](@ref). Hits of the same rank keep their order. No hit (`nothing`) stays
`nothing`.
"""
first_to_interact(::Nothing, ray) = nothing

function first_to_interact(hit::Intersection, ray::AbstractRay)
    isnothing(hit.next) && return hit
    best, best_rank = hit, _rank(hit, ray)
    for h in coincident_hits(hit)
        r = _rank(h, ray)
        if r < best_rank
            best, best_rank = h, r
        end
    end
    return first_of(hit, best)
end

# the part is not inferred: the rank is computed behind a function barrier
_rank(hit::Intersection, ray::AbstractRay) = _rank(part_of(hit), hit, ray)
_rank(part, hit, ray) = interaction_rank(hit, ray, medium(part))

"""
    first_of(hit, chosen)

Returns the chain of [`coincident_hits`](@ref) `hit` with the hit `chosen` of the chain first.
"""
function first_of(hit::Intersection, chosen::Intersection)
    chosen === hit && return hit
    prev = hit
    while prev.next !== chosen
        prev = prev.next
    end
    prev.next = chosen.next
    chosen.next = hit
    return chosen
end

"""
    Transition{T}

The boundary that a ray crosses at its [`Intersection`](@ref): the medium on each side and the thin
layer between them, see [`transition`](@ref).

# Fields

- `system`: the system of the boundary
- `λ`: the wavelength of the ray
- `n_in`, `n_out`: the refractive indices before and behind the boundary at `λ`
- `normal`: the unit normal of the boundary at the hit, against the direction of the ray
- `exiting`: the hit of the body that the ray leaves, or `nothing`
- `entering`: the hit of the body that the ray enters, or `nothing`
- `layer`: the hit of the thin layer at the boundary (a surface that is no body, e.g. a thin
  beamsplitter or a polarizing film), or `nothing`
- `enclosure_in`, `enclosure_out`: the [`Enclosure`](@ref) of a ray before and behind the boundary

# Media

The media on both sides are [`medium_in`](@ref) and [`medium_out`](@ref), each an
[`AbstractMedium`](@ref). They follow from the [`Enclosure`](@ref) of the ray, the bodies that
enclose it: before the boundary, the innermost body, or the body that the ray leaves there. Behind
the boundary, the body that the ray enters, else the body around the one that it leaves (e.g. the
water of a cuvette around a lens), else the medium of the system. A ray that a component created
without an enclosure is taken to be in the medium of its refractive index.

# Layers

At most one thin layer acts at a boundary. If several coincide, e.g. two films on top of each other,
the closest one acts and a warning is issued once.
"""
struct Transition{T}
    system::AbstractSystem
    λ::T
    n_in::T
    n_out::T
    normal::Point3{T}
    exiting::Nullable{Intersection{T}}
    entering::Nullable{Intersection{T}}
    layer::Nullable{Intersection{T}}
    enclosure_in::Union{Nothing, Enclosure}
    enclosure_out::Union{Nothing, Enclosure}
end

"""
    medium_in(t::Transition)

The [`AbstractMedium`](@ref) before the boundary, in which the ray propagates, see [`Transition`](@ref).
"""
medium_in(t::Transition) = _medium(t, innermost(t.enclosure_in), t.n_in)

"""
    medium_out(t::Transition)

The [`AbstractMedium`](@ref) behind the boundary, into which the ray is transmitted, see
[`Transition`](@ref).
"""
medium_out(t::Transition) = _medium(t, innermost(t.enclosure_out), t.n_out)

# outside of all bodies: the medium of the system, or of the index of a ray without an enclosure
_medium(::Transition, body::AbstractObject, n) = medium(body)
_medium(t::Transition, ::Nothing, n) =
    n == refractive_index(t.system, t.λ) ? medium(t.system) : IsotropicMedium(n)

"""The [`Enclosure`](@ref) of a ray that stays before the boundary, e.g. after a reflection."""
enclosure_in(t::Transition) = t.enclosure_in
"""The [`Enclosure`](@ref) of a ray that is transmitted through the boundary."""
enclosure_out(t::Transition) = t.enclosure_out

"""The refractive index before the boundary at the wavelength of the ray, see [`Transition`](@ref)."""
refractive_index_in(t::Transition) = t.n_in
"""The refractive index behind the boundary at the wavelength of the ray, see [`Transition`](@ref)."""
refractive_index_out(t::Transition) = t.n_out
normal3d(t::Transition) = t.normal

"""
    hint_in(t::Transition)

The [`Hint`](@ref) for a ray that stays in the medium before the boundary, e.g. after a reflection:
the body that the ray would have left, or `nothing`.
"""
hint_in(t::Transition) = Hint(t.exiting)

"""
    hint_out(t::Transition)

The [`Hint`](@ref) for a ray that is transmitted through the boundary: the body that it enters, or
`nothing`.
"""
hint_out(t::Transition) = Hint(t.entering)

Hint(hit::Intersection) = Hint(object(hit), shape(hit))
Hint(::Nothing) = nothing

"""
    transition(system::AbstractSystem, ray::AbstractRay)

Returns the [`Transition`](@ref) at the [`Intersection`](@ref) of the `ray`. The parts whose surfaces
the `ray` hits at this point (see [`coincident_hits`](@ref)) determine it: a part with a
[`medium`](@ref) is a body, which the `ray` leaves or enters depending on the normal of its hit, a
part without one is a thin layer. The [`Enclosure`](@ref) of the `ray` gives the bodies around it. The
normal of the transition is that of the hit that interacts with the `ray` (see
[`first_to_interact`](@ref)), turned against the `ray`.
"""
function transition(system::AbstractSystem, ray::AbstractRay{T}) where {T}
    hit = intersection(ray)::Intersection{T}
    λ = wavelength(ray)
    dir = direction(ray)
    exiting::Nullable{Intersection{T}} = nothing
    entering::Nullable{Intersection{T}} = nothing
    layer::Nullable{Intersection{T}} = nothing
    layers = 0
    for h in coincident_hits(hit)
        # the part is not inferred: whether it is a body is read behind a function barrier
        if !_is_body(part_of(h))
            layers += 1
            isnothing(layer) && (layer = h)
        elseif isentering(dir, normal3d(h))
            isnothing(entering) && (entering = h)
        else
            isnothing(exiting) && (exiting = h)
        end
    end
    layers > 1 && @warn "Several thin layers coincide at a boundary, only the closest one acts" maxlog = 1
    # the bodies around the ray: the body that it leaves encloses it, also for a ray that a component
    # created without an enclosure
    left = isnothing(exiting) ? nothing : part_of(exiting)
    e_in = _enclosing(enclosure(ray), left)
    e_out = leave(e_in, left)
    isnothing(entering) || (e_out = enter(e_out, part_of(entering)))
    n_in = _index(innermost(e_in), λ, T(refractive_index(ray)))
    n_out = isnothing(entering) && isnothing(exiting) ? n_in :
            _index(innermost(e_out), λ, T(refractive_index(system, λ)))
    normal = normal3d(hit)
    normal = isentering(dir, normal) ? normal : -normal
    return Transition{T}(system, T(λ), n_in, n_out, normal, exiting, entering, layer, e_in, e_out)
end

_is_body(part) = !isnothing(medium(part))

_enclosing(e, ::Nothing) = e
_enclosing(e, body::AbstractObject) = encloses(e, body) ? e : enter(e, body)

# the index of the innermost body (read behind a function barrier), else the index `n` outside of all
# bodies
_index(::Nothing, λ, n::T) where {T} = n
_index(body, λ, n::T) where {T} = T(_medium_index(body, λ)::Float64)

_medium_index(part, λ) = _medium_index(medium(part), λ)
_medium_index(::Nothing, λ) = NaN
_medium_index(m::AbstractMedium, λ) = Float64(refractive_index(m, λ))

"""
    AbstractSurfaceModel

The optical model of a boundary between two media, e.g. an uncoated surface or a coating. A model
`M` implements `interact3d(system, model::M, t::Transition, beam, ray)`, which returns the
interaction of the `ray` at the boundary of the transition `t`, see [`surface_model`](@ref).
"""
abstract type AbstractSurfaceModel end

"""
    Uncoated <: AbstractSurfaceModel

A bare boundary between two media: refraction or total internal reflection, see
[`cross_boundary`](@ref).
"""
struct Uncoated <: AbstractSurfaceModel end

"""
    surface_model(t::Transition)

Returns the [`AbstractSurfaceModel`](@ref) of the boundary of the transition `t`. All surfaces are
[`Uncoated`](@ref) so far; coatings on the faces of the bodies on both sides will enter here.
"""
surface_model(::Transition) = Uncoated()

"""
    interact3d(system, model::AbstractSurfaceModel, t::Transition, beam, ray)
    interact3d(system, layer::AbstractObject, t::Transition, beam, ray)

The interaction of the `ray` at the boundary of the transition `t`. A body, e.g. a lens, calls the
first method with the [`surface_model`](@ref) of the boundary. A thin layer, e.g. a
[`ThinBeamsplitter`](@ref), implements the second method for itself: it acts between the media of the
transition, e.g. it refracts the transmitted beam from one into the other.
"""
interact3d(::AbstractSystem, ::Uncoated, t::Transition, beam::Beam, ray::AbstractRay) =
    cross_boundary(t, beam, ray)

"""
    cross_boundary(t::Transition, beam::Beam, ray)

Returns the [`BeamInteraction`](@ref) of the `ray` with the bare boundary of the transition `t`:
refraction from the medium before into the medium behind the boundary via [`refraction3d`](@ref), or
total internal reflection beyond the critical angle. A [`PolarizedRay`](@ref) is transmitted or
reflected with the [`fresnel_coefficients`](@ref), only the stronger part is traced. Between equal
indices the ray passes unchanged. The hint points into the medium in which the new ray propagates.
"""
function cross_boundary(t::Transition, ::Beam{T, R}, ray::R) where {T <: Real, R <: Ray{T}}
    n1, n2 = refractive_index_in(t), refractive_index_out(t)
    ndir, TIR = n1 == n2 ? (direction(ray), false) : refraction3d(direction(ray), normal3d(t), n1, n2)
    n, hint, e = TIR ? (n1, hint_in(t), enclosure_in(t)) : (n2, hint_out(t), enclosure_out(t))
    return BeamInteraction{T, R}(hint, Ray{T}(hit_point(ray), ndir, nothing, wavelength(ray), n, e))
end

function cross_boundary(t::Transition, ::Beam{T, R}, ray::R) where {T <: Real, R <: PolarizedRay{T}}
    normal = normal3d(t)
    n1, n2 = refractive_index_in(t), refractive_index_out(t)
    if n1 == n2
        return BeamInteraction{T, R}(hint_out(t), PolarizedRay{T}(hit_point(ray), direction(ray),
            nothing, wavelength(ray), n2, polarization(ray), enclosure_out(t)))
    end
    # angle of incidence in the first quadrant
    θi = angle3d(direction(ray), -normal)
    rs, rp, ts, tp = fresnel_coefficients(θi, n2 / n1)
    if is_internally_reflected(rp, rs)
        new_dir = reflection3d(direction(ray), normal)
        J = SPBasis(-rs, 0, 0, rp)
        n, hint, e = n1, hint_in(t), enclosure_in(t)
    else
        new_dir, _ = refraction3d(direction(ray), normal, n1, n2)
        J = SPBasis(ts, 0, 0, tp)
        n, hint, e = n2, hint_out(t), enclosure_out(t)
    end
    P = _calculate_global_E0(direction(ray), new_dir, normal3d(intersection(ray)), J)
    E0 = _transverse(P * polarization(ray), new_dir)
    return BeamInteraction{T, R}(
        hint, PolarizedRay{T}(hit_point(ray), new_dir, nothing, wavelength(ray), n, E0, e))
end

"""Returns `true` if the `ray` is totally reflected at the bare boundary of the transition `t`."""
is_totally_reflected(ray::AbstractRay, t::Transition) =
    last(refraction3d(direction(ray), normal3d(t), refractive_index_in(t), refractive_index_out(t)))

"""
    cross_boundary(system::AbstractSystem, beam::Beam, ray)

The interaction of the `ray` with the bare boundary at its hit, see
`cross_boundary(transition(system, ray), beam, ray)`.
"""
cross_boundary(system::AbstractSystem, beam::Beam, ray::AbstractRay) =
    cross_boundary(transition(system, ray), beam, ray)
