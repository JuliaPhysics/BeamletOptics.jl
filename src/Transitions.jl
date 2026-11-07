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
    iscoated(hit::Intersection)

Returns `true` if faces of the part that was hit carry surface models, see [`surface_model!`](@ref).
"""
iscoated(hit::Intersection) = isnothing(hit.part) ? !isempty(face_models(part_of(hit))) : hit.iscoated

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
the closest one acts and a warning is issued once.
"""
struct Transition{T, S <: AbstractSystem}
    system::S
    λ::T
    n_in::T
    n_out::T
    normal::Point3{T}
    point::Point3{T}
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
                @warn "Several thin layers coincide at a boundary, only the closest one acts" maxlog = 1
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
    exiting, entering, layer = boundary_hits(hit, ray)
    # Bodies around the ray, including the body it leaves
    e_in = enter(enclosure(ray), part_of(exiting))
    n_in = _index(innermost(e_in), λ, T(refractive_index(ray)))
    n_out = if !isnothing(entering)
        _index(part_of(entering), λ, n_in)
    elseif !isnothing(exiting)
        _index(innermost(leave(e_in, part_of(exiting))), λ, T(refractive_index(system, λ)))
    else
        n_in
    end
    normal = normal_against(normal3d(hit), ray)
    return Transition{T, typeof(system)}(system, T(λ), n_in, n_out, normal, Point3{T}(hit_point(ray)),
        exiting, entering, layer, e_in)
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

# The component beams of a beamlet differ in their ray type, hence they are handled by `map` over the
# tuples, which is unrolled, not by a loop
function enclose_new!(system::AbstractSystem, beam::GaussianBeamlets, ray_id::Int,
        interaction::AbstractInteraction)
    map(_component_beams(beam), component_interactions(interaction)) do b, bi
        enclose!(system, rays(b)[ray_id], (bi.ray,))
    end
    return nothing
end

function enclose_new!(system::AbstractSystem, beam::GaussianBeamlets, ray_id::Int, ::Nothing)
    for c in children(beam)
        map(_component_beams(beam), _component_beams(c)) do b, cb
            enclose!(system, rays(b)[ray_id], (first(rays(cb)),))
        end
    end
    return nothing
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

- `channels(model, t::Transition, ray, cs)`: the [`BoundaryChannel`](@ref)s of the boundary of the
  transition `t` with the model, from the channels `cs` of the bare boundary, see [`channels`](@ref)
- `splitting_of(model)`: optional, which channels the solver traces, [`Dominant`](@ref) (default) or
  [`AllChannels`](@ref)
- `interact3d(system, model, t::Transition, beam::Beam, ray)`: optional, the interaction of the `ray`
  at the boundary of the transition `t`; by default the solver builds it from the channels, see
  [`realize`](@ref)

The solver builds the new rays, child beams and child beamlets of every beam type from the channels,
see [`beamlet_interaction`](@ref); a model needs no code per beam type.

A thin layer object carries its model via [`surface_model`](@ref)`(object)`.
"""
abstract type AbstractSurfaceModel end

"""
    interact3d(system, model::AbstractSurfaceModel, t::Transition, beam::Beam, ray)

The interaction of the `ray` with the `model` at the boundary of the transition `t`: the
[`channels`](@ref) of the bare boundary, transformed by the `model`, of which the solver traces those
that [`splitting_of`](@ref)`(model)` selects.
"""
function interact3d(::AbstractSystem, model::AbstractSurfaceModel, t::Transition, beam::Beam,
        ray::AbstractRay)
    cs = channels(model, t, ray, channels(t, ray))
    return realize(t, beam, ray, select_channels(splitting_of(model), cs))
end

# The splitting of the model decides how the beamlet interacts, see beamlet_interaction
interact3d(system::AbstractSystem, model::AbstractSurfaceModel, ::Transition,
    beamlet::GaussianBeamlets, ray_id::Int) =
    beamlet_interaction(system, splitting_of(model), model, beamlet, ray_id)

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

Returns the [`AbstractSurfaceModel`](@ref) that acts at the boundary of the transition `t`, the first
of
1. the model of its thin layer, e.g. the [`BeamsplitterCoating`](@ref) of a [`ThinBeamsplitter`](@ref);
   a thin layer without a model, e.g. a [`Detector`](@ref) that does not stop the ray, has none,
2. the model of the coated face of the body that the ray leaves, then of the body that it enters, see
   [`surface_model!`](@ref),
3. the model of the system, see [`surface_model(::AbstractSystem)`](@ref), [`Uncoated`](@ref) by
   default.
"""
surface_model(t::Transition) = something(boundary_model(t), surface_model(t.system))

# The model of the thin layer or of a coated face at the boundary of `t`, or `nothing`
function boundary_model(t::Transition)
    model = surface_model(t.layer)
    isnothing(model) || return model
    model = face_model(t, t.exiting)
    isnothing(model) || return model
    return face_model(t, t.entering)
end

# The model of the coated face of a body at the boundary of `t`
face_model(::Transition, ::Nothing) = nothing
face_model(t::Transition, hit::Intersection) =
    iscoated(hit) ? face_model(part_of(hit), hit, t.point) : nothing

"""
    surface_model(system::AbstractSystem)

Returns the surface model of every boundary of the `system` without a thin layer or a coated face:
the keyword `surface_model` of [`System`](@ref) and [`StaticSystem`](@ref), e.g. [`FresnelGhosts`](@ref)
to trace the ghost reflections of all uncoated surfaces, [`Uncoated`](@ref) by default.
"""
surface_model(::AbstractSystem) = Uncoated()
surface_model(system::Union{System, StaticSystem}) = system.surface_model

# The hit of the surface that carries the model of `t`, in the order of `surface_model(t)`
function carrier(t::Transition)
    isnothing(surface_model(t.layer)) || return t.layer
    isnothing(face_model(t, t.exiting)) || return t.exiting
    isnothing(face_model(t, t.entering)) || return t.entering
    return something(t.layer, t.entering, t.exiting)
end

"""
    isfacing(t::Transition, ray)

Returns `true` if the `ray` comes from the side into which the surface that carries the model of the
transition `t` faces: the normal of a thin layer, or the outside of a coated body. A coating lies on
this side of its substrate.
"""
isfacing(t::Transition, ray::AbstractRay) = isentering(direction(ray), normal3d(carrier(t)))

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
function interact3d(system::AbstractSystem, t::Transition, beam::AbstractBeam,
        at::Union{AbstractRay, Int})
    model = boundary_model(t)
    # The model of the system has the type that the system knows
    isnothing(model) && return interact3d(system, surface_model(t.system), t, beam, at)
    return interact3d(system, model, t, beam, at)
end

"""
    cross_boundary(t::Transition, beam::Beam, ray)

Returns the [`BeamInteraction`](@ref) of the `ray` with the bare boundary of the transition `t`, the
[`Uncoated`](@ref) boundary: refraction from the medium before into the medium behind the boundary, or
total internal reflection, see [`channels`](@ref). The hint points into the medium in which the new ray
propagates.
"""
cross_boundary(t::Transition, beam::Beam, ray::AbstractRay) =
    realize(t, beam, ray, select_channels(Dominant(), channels(t, ray)))

"""
    start_enclosure(system, point, dir)

Returns the [`Enclosure`](@ref) of a ray that starts at the `point` along `dir`: the bodies of the
`system` that contain it (see [`isinside`](@ref)), innermost first, or `nothing` outside of all bodies.
A ray that starts on the surface of a body is within it if it runs into it: the bodies are tested at
the distance `get_coincidence_tolerance()` along `dir`. Nested bodies
are ordered by the size of their bounding spheres, the larger one outside. The bounding spheres of the
solve skip the bodies far from the `point`.
"""
function start_enclosure(system::AbstractSystem, point::AbstractVector, dir::AbstractVector)
    table = current_bounding_spheres()
    # Just ahead of the start, which decides for a start on a surface
    point = point + get_coincidence_tolerance() * dir
    bodies = AbstractObject[]
    foreach(x -> collect_bodies!(bodies, table, x, point), _top_level(system))
    isempty(bodies) && return nothing
    sort!(bodies; by = body -> enclosing_radius(bounding_sphere_of(table, body)), rev = true)
    return foldl(enter, bodies; init = nothing)
end

# The bodies within `x` that contain the `point`, pushed onto `bodies`
collect_bodies!(bodies, table, object::AbstractObject, point) =
    collect_bodies!(bodies, table, shape_trait_of(object), object, point)
collect_bodies!(bodies, table, group::AbstractObjectGroup, point) =
    foreach(o -> collect_bodies!(bodies, table, o, point), objects(group))
# A plain shape of a composite, or anything else a custom system holds, is no body
collect_bodies!(bodies, table, ::Any, point) = bodies

function collect_bodies!(bodies, table, ::SingleShape, object::AbstractObject, point)
    isnothing(medium(object)) && return bodies
    isinside(bounding_sphere_of(table, object), point) || return bodies
    isinside(shape(object), point) && push!(bodies, object)
    return bodies
end

function collect_bodies!(bodies, table, ::MultiShape, object::AbstractObject, point)
    isinside(bounding_sphere_of(table, object), point) || return bodies
    foreach(part -> collect_bodies!(bodies, table, part, point), shape(object))
    return bodies
end

enclosing_radius(sphere::AbstractBoundingSphere) = radius(sphere)
enclosing_radius(::NoBoundingSphere) = Inf

"""
    enclose_start!(system, beam)

Sets the [`Enclosure`](@ref) of the first ray of each component of the root `beam`, see
[`start_enclosure`](@ref), and its refractive index: that of the innermost body that contains its
start, else that of the system. A source that starts within a lens or a cuvette starts in its medium.
"""
function enclose_start!(system::AbstractSystem, beam::AbstractBeam)
    foreach(_component_beams(beam)) do b
        ray = first(rays(b))
        λ = wavelength(ray)
        enclosure = start_enclosure(system, position(ray), direction(ray))
        enclosure!(ray, enclosure)
        body = innermost(enclosure)
        refractive_index!(ray, isnothing(body) ? refractive_index(system, λ) : medium_index(body, λ))
    end
    return nothing
end
