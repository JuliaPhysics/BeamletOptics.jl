#=
Transitions: what a ray crosses at a boundary

The solver collects all surfaces that a ray hits at one point, see `coincident_hits`. Of these, one
hit interacts with the ray (`first_to_interact`), and the bodies that the ray leaves and enters
there define the media on both sides of the boundary (`transition`). Hence touching elements, e.g.
cemented lenses or a coating on its substrate, need no knowledge of each other.
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
    body_refractive_index(part, λ)

Returns the refractive index of the body of the `part` at the wavelength `λ`, or `nothing` if the
`part` is no refracting body. A ray that enters such a body propagates within it, see
[`Transition`](@ref).
"""
body_refractive_index(::AbstractObject, λ) = nothing

"""
    interaction_rank(hit, ray, n)

Returns the rank of the `hit` of the `ray` among the [`coincident_hits`](@ref): the hit of the lowest
rank interacts with the `ray`. `n` is the [`body_refractive_index`](@ref) of the part that was hit. A
surface that is no refracting body (`n` is `nothing`) ranks first, e.g. a coating, a mirror or a
detector on a lens. The surface of a body that the `ray` enters ranks before the surface of a body
that it leaves.
"""
interaction_rank(hit, ray, ::Nothing) = 0
interaction_rank(hit, ray, ::Real) = isentering(direction(ray), normal3d(hit)) ? 1 : 2

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
_rank(part, hit, ray) = interaction_rank(hit, ray, body_refractive_index(part, wavelength(ray)))

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

The boundary that a ray crosses at its [`Intersection`](@ref), see [`transition`](@ref).

# Fields

- `medium_in`: the medium before the boundary, in which the ray propagates
- `medium_out`: the medium behind the boundary, into which the ray is transmitted
- `normal`: the unit normal of the boundary at the hit, against the direction of the ray
- `exiting`: the hit of the body that the ray leaves, or `nothing`
- `entering`: the hit of the body that the ray enters, or `nothing`

# Media

A medium is represented by its refractive index at the wavelength of the ray. Behind a boundary at
which the ray leaves a body and enters none, the ray propagates in the medium of the system, see
[`refractive_index(::AbstractSystem, λ)`](@ref). Bodies within bodies (e.g. a lens immersed in a
cuvette) are not resolved: the medium behind the boundary of the inner body is the medium of the
system.
"""
struct Transition{T}
    medium_in::T
    medium_out::T
    normal::Point3{T}
    exiting::Nullable{Intersection{T}}
    entering::Nullable{Intersection{T}}
end

"""Refractive index of the medium before the boundary, see [`Transition`](@ref)."""
medium_in(t::Transition) = t.medium_in
"""Refractive index of the medium behind the boundary, see [`Transition`](@ref)."""
medium_out(t::Transition) = t.medium_out
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

Returns the [`Transition`](@ref) at the [`Intersection`](@ref) of the `ray`. The bodies whose surfaces
the `ray` hits at this point (see [`coincident_hits`](@ref)) determine the media:

- the `ray` leaves a body: the medium before the boundary is the body, else the refractive index of
  the `ray`
- the `ray` enters a body: the medium behind the boundary is the body; it leaves a body and enters
  none: the medium of the `system`; else the medium before the boundary, e.g. at a coating within a
  body or in air

A body is a part for which [`body_refractive_index`](@ref) returns a value, see
[`AbstractRefractiveOptic`](@ref). The normal is the normal of the hit that interacts with the `ray`,
see [`first_to_interact`](@ref), turned against the `ray`.
"""
_body_index(part, λ) = _as_index(body_refractive_index(part, λ))
_as_index(::Nothing) = NaN
_as_index(n::Real) = Float64(n)

function transition(system::AbstractSystem, ray::AbstractRay{T}) where {T}
    hit = intersection(ray)::Intersection{T}
    λ = wavelength(ray)
    dir = direction(ray)
    exiting::Nullable{Intersection{T}} = nothing
    entering::Nullable{Intersection{T}} = nothing
    n_exiting = n_entering = T(refractive_index(ray))
    for h in coincident_hits(hit)
        # the part is not inferred: its index is read behind a function barrier, NaN for no body
        n = _body_index(part_of(h), λ)::Float64
        isnan(n) && continue
        if isentering(dir, normal3d(h))
            isnothing(entering) && ((entering, n_entering) = (h, T(n)))
        else
            isnothing(exiting) && ((exiting, n_exiting) = (h, T(n)))
        end
    end
    n_in = n_exiting
    n_out = if !isnothing(entering)
        n_entering
    elseif !isnothing(exiting)
        T(refractive_index(system, λ))
    else
        n_in
    end
    normal = normal3d(hit)
    normal = isentering(dir, normal) ? normal : -normal
    return Transition{T}(n_in, n_out, normal, exiting, entering)
end

"""
    cross_boundary(t::Transition, beam::Beam, ray)

Returns the [`BeamInteraction`](@ref) of the `ray` with the bare boundary of the transition `t`:
refraction from `medium_in` into `medium_out` via [`refraction3d`](@ref), or total internal reflection
beyond the critical angle. A [`PolarizedRay`](@ref) is transmitted or reflected with the
[`fresnel_coefficients`](@ref), only the stronger part is traced. The hint points into the medium in
which the new ray propagates.
"""
function cross_boundary(t::Transition, ::Beam{T, R}, ray::R) where {T <: Real, R <: Ray{T}}
    ndir, TIR = refraction3d(direction(ray), normal3d(t), medium_in(t), medium_out(t))
    n, hint = TIR ? (medium_in(t), hint_in(t)) : (medium_out(t), hint_out(t))
    return BeamInteraction{T, R}(hint, Ray{T}(hit_point(ray), ndir, nothing, wavelength(ray), n))
end

function cross_boundary(t::Transition, ::Beam{T, R}, ray::R) where {T <: Real, R <: PolarizedRay{T}}
    normal = normal3d(t)
    n1, n2 = medium_in(t), medium_out(t)
    # angle of incidence in the first quadrant
    θi = angle3d(direction(ray), -normal)
    rs, rp, ts, tp = fresnel_coefficients(θi, n2 / n1)
    if is_internally_reflected(rp, rs)
        new_dir = reflection3d(direction(ray), normal)
        J = SPBasis(-rs, 0, 0, rp)
        n, hint = n1, hint_in(t)
    else
        new_dir, _ = refraction3d(direction(ray), normal, n1, n2)
        J = SPBasis(ts, 0, 0, tp)
        n, hint = n2, hint_out(t)
    end
    P = _calculate_global_E0(direction(ray), new_dir, normal3d(intersection(ray)), J)
    E0 = _transverse(P * polarization(ray), new_dir)
    return BeamInteraction{T, R}(
        hint, PolarizedRay{T}(hit_point(ray), new_dir, nothing, wavelength(ray), n, E0))
end

"""Returns `true` if the `ray` is totally reflected at the bare boundary of the transition `t`."""
is_totally_reflected(ray::AbstractRay, t::Transition) =
    last(refraction3d(direction(ray), normal3d(t), medium_in(t), medium_out(t)))

"""
    cross_boundary(system::AbstractSystem, beam::Beam, ray)

The interaction of the `ray` with the bare boundary at its hit, see
`cross_boundary(transition(system, ray), beam, ray)`.
"""
cross_boundary(system::AbstractSystem, beam::Beam, ray::AbstractRay) =
    cross_boundary(transition(system, ray), beam, ray)
