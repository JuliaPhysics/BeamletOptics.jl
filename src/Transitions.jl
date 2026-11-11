#=
Transitions: what a ray crosses at a boundary

The solver collects all surfaces that a ray hits at one point, see `coincident_hits`. The bodies
that the ray leaves and enters there define the media on both sides of the boundary, and the thin
layers at it (coatings, films, beamsplitters) act between them, see `transition`. Hence touching
elements, e.g. cemented lenses or a coating on its substrate, need no knowledge of each other.
=#

"""
    find_part(object, target)

Returns the object within `object` whose shape is `target`, or `nothing`. A part that is a plain
shape holds no object.
"""
find_part(::AbstractShape, target) = nothing
find_part(object::AbstractObject, target) = find_part(shape_trait_of(object), object, target)
find_part(::SingleShape, object::AbstractObject, target) = shape(object) === target ? object : nothing

function find_part(::MultiShape, object::AbstractObject, target)
    for part in shape(object)
        found = find_part(part, target)
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

"""Returns the part of the object that the hit `hit` belongs to, see [`part_of`](@ref). No hit
(`nothing`) belongs to no part."""
part_of(hit::Intersection) = isnothing(hit.part) ? part_of(object(hit), shape(hit)) : hit.part
part_of(::Nothing) = nothing

"""
    isbody(hit::Intersection)

Returns `true` if the part that the `hit` belongs to is a body, i.e. has a [`medium`](@ref). The solver
stores this with the part of the hit; a hit of an object with its own `intersect3d` is looked up.
"""
isbody(hit::Intersection) = isnothing(hit.part) ? !isnothing(medium(part_of(hit))) : hit.isbody

"""
    interaction_rank(hit, ray)

Returns the rank of the `hit` of the `ray` among the [`coincident_hits`](@ref): the hit of the lowest
rank interacts with the `ray`. A surface that is no body (see [`isbody`](@ref)) ranks first, e.g. a
thin beamsplitter, a mirror or a detector on a lens: it acts at the boundary, see
[`Transition`](@ref). The surface of a body that the `ray` enters ranks before the surface of a body
that it leaves.
"""
interaction_rank(hit::Intersection, ray::AbstractRay) =
    !isbody(hit) ? 0 : isentering(direction(ray), normal3d(hit)) ? 1 : 2

"""
    first_to_interact(hit, ray)

Returns the chain of [`coincident_hits`](@ref) `hit` of the `ray` with the hit that interacts first,
see [`interaction_rank`](@ref). Hits of the same rank keep their order. No hit (`nothing`) stays
`nothing`.
"""
first_to_interact(::Nothing, ray) = nothing

function first_to_interact(hit::Intersection, ray::AbstractRay)
    isnothing(hit.next) && return hit
    best, best_rank = hit, interaction_rank(hit, ray)
    for h in coincident_hits(hit)
        r = interaction_rank(h, ray)
        if r < best_rank
            best, best_rank = h, r
        end
    end
    return first_of(hit, best)
end

"""
    Transition{T, S}

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
- `enclosure`: the [`Enclosure`](@ref) of the ray before the boundary, see [`enclosure_in`](@ref) and
  [`enclosure_out`](@ref)

# Media

The media on both sides are [`medium_in`](@ref) and [`medium_out`](@ref), each an
[`AbstractMedium`](@ref). They follow from the [`Enclosure`](@ref) of the ray, the bodies that
enclose it: before the boundary, the innermost body, or the body that the ray leaves there. Behind
the boundary, the body that the ray enters, else the body around the one that it leaves (e.g. the
water of a cuvette around a lens), else the medium of the system. The solver assigns the enclosure
of every new ray that starts at the boundary, see [`enclose!`](@ref).

# Layers

At most one thin layer acts at a boundary. If several coincide, e.g. two films on top of each other,
the first one in the chain of [`coincident_hits`](@ref) acts, which the geometry does not define, and a
warning is issued once. Films on top of each other belong into the model of one layer instead.
"""
struct Transition{T, S <: AbstractSystem}
    system::S
    λ::T
    n_in::T
    n_out::T
    normal::Point3{T}
    exiting::Nullable{Intersection{T}}
    entering::Nullable{Intersection{T}}
    layer::Nullable{Intersection{T}}
    enclosure::Union{Nothing, Enclosure}
end

"""
    medium_in(t::Transition)

The [`AbstractMedium`](@ref) before the boundary, in which the ray propagates, see [`Transition`](@ref).
"""
medium_in(t::Transition) = side_medium(t, innermost(enclosure_in(t)), t.n_in)

"""
    medium_out(t::Transition)

The [`AbstractMedium`](@ref) behind the boundary, into which the ray is transmitted, see
[`Transition`](@ref).
"""
medium_out(t::Transition) = side_medium(t, innermost(enclosure_out(t)), t.n_out)

# Outside of all bodies: the medium of the system, or of the index of a ray without an enclosure
side_medium(::Transition, body::AbstractObject, n) = medium(body)
side_medium(t::Transition, ::Nothing, n) =
    n == refractive_index(t.system, t.λ) ? medium(t.system) : IsotropicMedium(ConstantRefractiveIndex(n))

"""The refractive index before the boundary at the wavelength of the ray, see [`Transition`](@ref)."""
refractive_index_in(t::Transition) = t.n_in
"""The refractive index behind the boundary at the wavelength of the ray, see [`Transition`](@ref)."""
refractive_index_out(t::Transition) = t.n_out
normal3d(t::Transition) = t.normal

"""The [`Enclosure`](@ref) of a ray that stays before the boundary, e.g. after a reflection."""
enclosure_in(t::Transition) = t.enclosure

"""The [`Enclosure`](@ref) of a ray that is transmitted through the boundary."""
enclosure_out(t::Transition) = enclosure_beyond(t.enclosure, t.exiting, t.entering)

# The enclosure `e` without the body of the hit `exiting` and within the body of the hit `entering`
enclosure_beyond(e::Nullable{Enclosure}, exiting, entering) =
    enter(leave(e, part_of(exiting)), part_of(entering))

"""
    enclose!(system::AbstractSystem, ray::AbstractRay, new_rays)

Sets the [`Enclosure`](@ref) of each of the `new_rays` that starts at the hit of the `ray`: the
continued ray and the first rays of child beams. The side of the boundary into which a new ray
propagates decides: back into the medium before the boundary ([`enclosure_in`](@ref)), or into
the medium behind it ([`enclosure_out`](@ref)). The solver calls it after each interaction, hence a
component needs no knowledge of the bodies around the ray. A new ray that starts elsewhere, e.g. one
that a component re-emits at another place, keeps its enclosure.
"""
function enclose!(::AbstractSystem, ray::AbstractRay{T}, new_rays) where {T}
    hit = intersection(ray)::Intersection{T}
    p = hit_point(ray)
    exiting, entering, _ = boundary_hits(hit, ray)
    e_in = enter(enclosure(ray), part_of(exiting))
    normal = normal_against(normal3d(hit), ray)
    e_out = nothing
    for r in new_rays
        norm(position(r) - p) ≤ get_coincidence_tolerance() || continue
        if isentering(direction(r), normal)
            isnothing(e_out) && (e_out = enclosure_beyond(e_in, exiting, entering))
            enclosure!(r, e_out)
        else
            enclosure!(r, e_in)
        end
    end
    return nothing
end

# The `normal` turned against the direction of the `ray`
normal_against(normal::AbstractVector, ray::AbstractRay) =
    isentering(direction(ray), normal) ? normal : -normal

"""
    boundary_hits(hit::Intersection, ray) -> (exiting, entering, layer)

Sorts the [`coincident_hits`](@ref) `hit` of the `ray` by their role at the boundary, see
[`Transition`](@ref): the first hit of a body that the `ray` leaves, of a body that it enters and of a
thin layer, each `nothing` if there is none. Several thin layers warn once, only the first one acts.
"""
function boundary_hits(hit::Intersection{T}, ray::AbstractRay) where {T}
    dir = direction(ray)
    exiting::Nullable{Intersection{T}} = nothing
    entering::Nullable{Intersection{T}} = nothing
    layer::Nullable{Intersection{T}} = nothing
    for h in coincident_hits(hit)
        if !isbody(h)
            if isnothing(layer)
                layer = h
            else
                @warn "Several thin layers coincide at a boundary, only one of them acts" maxlog = 1
            end
        elseif isentering(dir, normal3d(h))
            isnothing(entering) && (entering = h)
        else
            isnothing(exiting) && (exiting = h)
        end
    end
    return exiting, entering, layer
end

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

Hint(hit::Intersection) = Hint(object(hit)::AbstractObject, shape(hit)::AbstractShape)
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
    exiting, entering, layer = boundary_hits(hit, ray)
    # Bodies around the ray, including the body it leaves
    e_in = enter(enclosure(ray), part_of(exiting))
    # Index of the ray, unless it leaves a body outside of its enclosure
    n_ray = T(refractive_index(ray))
    n_in = e_in === enclosure(ray) ? n_ray : _index(innermost(e_in), λ, n_ray)
    n_out = if !isnothing(entering)
        _index(part_of(entering), λ, n_in)
    elseif !isnothing(exiting)
        _index(innermost(leave(e_in, part_of(exiting))), λ, T(refractive_index(system, λ)))
    else
        n_in
    end
    normal = normal_against(normal3d(hit), ray)
    return Transition{T, typeof(system)}(system, T(λ), n_in, n_out, normal, exiting, entering, layer,
        e_in)
end

"""
    transition(system::AbstractSystem, beam::AbstractBeam, at)

The [`Transition`](@ref) at the hit of the `beam` at its ray `at`: of the ray itself for a
[`Beam`](@ref), of the chief ray of the segment `at` for a beamlet, see `intersection(beam, at)`.
"""
transition(system::AbstractSystem, ::Beam, ray::AbstractRay) = transition(system, ray)
transition(system::AbstractSystem, beamlet::GaussianBeamlets, ray_id::Int) =
    transition(system, rays(beamlet)[ray_id])

# The index of the innermost body in the number type `T` of the ray, e.g. a `ForwardDiff.Dual` within
# automatic differentiation (read behind a function barrier), else the index `n` outside of all bodies
_index(::Nothing, λ, n::T) where {T} = n
_index(body::AbstractObject, λ, n::T) where {T} = convert(T, medium_index(body, λ))::T

medium_index(part::AbstractObject, λ) = refractive_index(medium(part), λ)

"""
    enclose_new!(system, beam, ray, interaction)
    enclose_new!(system, beamlet, ray_id, interaction)

Sets the [`Enclosure`](@ref) of the rays that an interaction at the hit of the `ray` created: the
continued ray of the `interaction`, or, without one, the first rays of the children of the `beam`.
For a beamlet, each component ray at the segment `ray_id` is handled on its own. See
[`enclose!`](@ref).
"""
enclose_new!(system::AbstractSystem, ::Beam, ray::AbstractRay, interaction::BeamInteraction) =
    enclose!(system, ray, (interaction.ray,))

enclose_new!(system::AbstractSystem, beam::Beam, ray::AbstractRay, ::Nothing) =
    enclose!(system, ray, (first(rays(c)) for c in children(beam)))

enclose_new!(system::AbstractSystem, beam::GaussianBeamlets, ray_id::Int,
    interaction::AbstractInteraction) =
    enclose_components!(system, ray_id, _component_beams(beam),
        map(bi -> bi.ray, component_interactions(interaction)))

function enclose_new!(system::AbstractSystem, beam::GaussianBeamlets, ray_id::Int, ::Nothing)
    for c in children(beam)
        enclose_components!(system, ray_id, _component_beams(beam),
            map(cb -> first(rays(cb)), _component_beams(c)))
    end
    return nothing
end

# The component beams of a beamlet differ in their ray type, hence they are handled by recursion over
# the tuples, which is unrolled, not by a loop; a closure over the system would be allocated
enclose_components!(::AbstractSystem, ::Int, ::Tuple{}, ::Tuple{}) = nothing

@inline function enclose_components!(system::AbstractSystem, ray_id::Int, beams::Tuple, new_rays::Tuple)
    enclose!(system, rays(first(beams))[ray_id], (first(new_rays),))
    return enclose_components!(system, ray_id, Base.tail(beams), Base.tail(new_rays))
end

component_interactions(interaction::GaussianBeamletInteraction) =
    (interaction.chief, interaction.waist, interaction.divergence)
component_interactions(interaction::AstigmaticGaussianBeamletInteraction) =
    (interaction.chief, interaction.wxp, interaction.wxm, interaction.wyp, interaction.wym,
        interaction.dxp, interaction.dxm, interaction.dyp, interaction.dym)

"""
    AbstractSurfaceModel

The optical model of a boundary between two media, e.g. an uncoated surface or a coating, see
[`surface_model`](@ref).

# Implementation reqs.

## Fields

None required.

## Getters/setters

None required.

## Functions

- `interact3d(system, model, t::Transition, beam::Beam, ray)`: the interaction of the `ray` at the
  boundary of the transition `t`; without a method, the boundary warns and blocks the beam
- `interact3d(system, model, t::Transition, beamlet, ray_id)`: optional, the interaction of a
  [`GaussianBeamlet`](@ref) or an [`AstigmaticGaussianBeamlet`](@ref) at the segment `ray_id`. By
  default, each component ray interacts with the model at its own hit; a model that splits the beamlet,
  e.g. a [`BeamsplitterCoating`](@ref), implements it.

A thin layer object carries its model via [`surface_model`](@ref)`(object)`.
"""
abstract type AbstractSurfaceModel end

function interact3d(::AbstractSystem, model::AbstractSurfaceModel, ::Transition, ::Beam,
        ray::AbstractRay)
    @warn lazy"No interact3d method defined for:" ModelType = typeof(model) RayType = typeof(ray)
    return nothing
end

# Each component ray of the beamlet interacts with the model at its own hit
function interact3d(system::AbstractSystem, model::AbstractSurfaceModel, ::Transition,
        beamlet::GaussianBeamlets, ray_id::Int)
    return beamlet_interaction(beamlet) do b
        ray = rays(b)[ray_id]
        interact3d(system, model, transition(system, ray), b, ray)
    end
end

"""
    interact3d(system, ::SingleShape, object, beamlet, ray_id)

The interaction of a [`GaussianBeamlet`](@ref) or an [`AstigmaticGaussianBeamlet`](@ref) with an
`object` without a beamlet method of its own. A thin layer with a model (see [`surface_model`](@ref))
hands the beamlet to its model, which may split it, e.g. a [`BeamsplitterCoating`](@ref). Any other
object interacts with each component ray, e.g. a body refracts each at its own hit. A
[`MultiShape`](@ref) object hands the beamlet to the part that the chief ray hits.
"""
function interact3d(system::AbstractSystem, ::SingleShape, object::AbstractObject,
        beamlet::GaussianBeamlets, ray_id::Int)
    isnothing(surface_model(object)) ||
        return interact3d(system, transition(system, beamlet, ray_id), beamlet, ray_id)
    return beamlet_interaction(b -> interact3d(system, object, b, rays(b)[ray_id]), beamlet)
end

# The interaction of a beamlet from the interaction `f(b)` of each component beam `b`, or `nothing`
# if one of them ends; map unrolls the heterogeneous tuple of component beams. `F` forces
# specialization on `f`, which is only passed on to `map`
function beamlet_interaction(f::F, beamlet::GaussianBeamlets) where {F}
    interactions = map(f, _component_beams(beamlet))
    any(isnothing, interactions) && return nothing
    return beamlet_interaction(beamlet, interactions)
end

beamlet_interaction(::GaussianBeamlet{R}, interactions::Tuple) where {R} =
    GaussianBeamletInteraction{R}(interactions...)
beamlet_interaction(::AstigmaticGaussianBeamlet{R}, interactions::Tuple) where {R} =
    AstigmaticGaussianBeamletInteraction{R}(interactions...)

"""
    Uncoated <: AbstractSurfaceModel

A bare boundary between two media: refraction or total internal reflection, see
[`cross_boundary`](@ref).
"""
struct Uncoated <: AbstractSurfaceModel end

"""
    surface_model(t::Transition)

Returns the [`AbstractSurfaceModel`](@ref) that acts at the boundary of the transition `t`: the model
of its thin layer, e.g. the [`BeamsplitterCoating`](@ref) of a [`ThinBeamsplitter`](@ref), else
[`Uncoated`](@ref). A thin layer without a model, e.g. a [`Detector`](@ref) that does not stop the
ray, leaves the boundary as it is. Coatings on the faces of bodies will enter here.
"""
surface_model(t::Transition) = something(surface_model(t.layer), Uncoated())

"""
    surface_model(object::AbstractObject)

Returns the [`AbstractSurfaceModel`](@ref) of the thin layer `object`, or `nothing` (default). It acts
at every boundary at which the `object` lies, between the media on both sides, see
[`surface_model(::Transition)`](@ref).
"""
surface_model(::AbstractObject) = nothing
surface_model(hit::Intersection) = surface_model(part_of(hit))
surface_model(::Nothing) = nothing

"""
    interact3d(system, t::Transition, beam, at)

The interaction of the `beam` at the boundary of the transition `t`, at its ray `at` (see
[`transition`](@ref)): the [`surface_model`](@ref) of the boundary acts. Bodies and thin layers without
an `interact3d` method of their own interact here, an own method may call it, e.g. a [`Detector`](@ref)
that does not stop the ray.
"""
interact3d(system::AbstractSystem, t::Transition, beam::AbstractBeam, at::Union{AbstractRay, Int}) =
    interact3d(system, surface_model(t), t, beam, at)

"""
    interact3d(system, model::AbstractSurfaceModel, t::Transition, beam, ray)

The interaction of the `ray` with the `model` at the boundary of the transition `t`. The model acts
between the media of the transition, e.g. an [`Uncoated`](@ref) boundary refracts from one into the
other, see [`cross_boundary`](@ref).
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

For a [`PolarizedRay`](@ref), `cross_boundary(t, beam, ray, E)` crosses the boundary with the field `E`
instead of the field of the `ray`, e.g. the field behind a polarizing film.
"""
function cross_boundary(t::Transition, ::Beam{T, R}, ray::R) where {T <: Real, R <: Ray{T}}
    n1, n2 = refractive_index_in(t), refractive_index_out(t)
    ndir, TIR = n1 == n2 ? (direction(ray), false) : refraction3d(direction(ray), normal3d(t), n1, n2)
    n, hint = TIR ? (n1, hint_in(t)) : (n2, hint_out(t))
    return BeamInteraction{T, R}(hint, Ray{T}(hit_point(ray), ndir, nothing, wavelength(ray), n))
end

cross_boundary(t::Transition, beam::Beam{T, R}, ray::R) where {T <: Real, R <: PolarizedRay{T}} =
    cross_boundary(t, beam, ray, polarization(ray))

function cross_boundary(t::Transition, ::Beam{T, R}, ray::R, E) where {T <: Real, R <: PolarizedRay{T}}
    normal = normal3d(t)
    n1, n2 = refractive_index_in(t), refractive_index_out(t)
    if n1 == n2
        return BeamInteraction{T, R}(hint_out(t), PolarizedRay{T}(hit_point(ray), direction(ray),
            nothing, wavelength(ray), n2, E))
    end
    # Fresnel coefficients from the cosine, smooth for AD at normal incidence
    rs, rp, ts, tp = _fresnel_coefficients(cos_angle3d(direction(ray), -normal), n2 / n1)
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
    E0 = _transverse(P * E, new_dir)
    return BeamInteraction{T, R}(
        hint, PolarizedRay{T}(hit_point(ray), new_dir, nothing, wavelength(ray), n, E0))
end

"""Returns `true` if the `ray` is totally reflected at the bare boundary of the transition `t`."""
is_totally_reflected(ray::AbstractRay, t::Transition) =
    last(refraction3d(direction(ray), normal3d(t), refractive_index_in(t), refractive_index_out(t)))
