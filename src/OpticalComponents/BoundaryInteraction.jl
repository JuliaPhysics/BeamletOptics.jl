"""
    current_intersection(beam::AbstractBeam, obj::AbstractObject, ray::AbstractRay) -> Nullable{AbstractIntersection}

Returns the surface intersection for the current interaction step: the intersection
already carried by `ray`, or by the most recently traced ray in `beam`, or a freshly computed one via
[`intersect3d`](@ref) against `obj` if `ray` has not been resolved yet (e.g. when [`interact3d`](@ref) is called
directly with an unresolved ray instead of via [`trace_system!`](@ref)).

# Arguments

- `beam`: The beam containing previously traced ray segments.
- `obj`: The optical component being hit.
- `ray`: The incoming ray segment.
"""
function current_intersection(beam::AbstractBeam, obj::AbstractObject, ray::AbstractRay)
    ray_int = intersection(ray)
    !isnothing(ray_int) && return ray_int
    r_vec = rays(beam)
    last_int = isempty(r_vec) ? nothing : intersection(last(r_vec))
    return isnothing(last_int) ? intersect3d(obj, ray) : last_int
end

"""
    transition_for(optic::AbstractObject, system::AbstractSystem, ray::AbstractRay, int::AbstractIntersection) -> Transition

Resolves the media [`Transition`](@ref) across `optic`'s boundary for the universal
boundary interaction pipeline (see [`interact3d`](@ref)).
Subtypes of [`AbstractRefractiveOptic`](@ref), [`AbstractReflectiveOptic`](@ref) and
[`AbstractCoating`](@ref) may specialize this to override the default per-family behavior.

# Arguments

- `optic`: The optical component being traversed.
- `system`: The optical system providing the ambient medium.
- `ray`: The incident ray.
- `int`: The intersection on the optic's boundary.
"""
transition_for(optic::AbstractRefractiveOptic, system::AbstractSystem, ray::AbstractRay, int::AbstractIntersection) = resolve_transition(
    medium_from(optic), ambient_medium(system), ray, normal3d(int))

function transition_for(::AbstractReflectiveOptic, system::AbstractSystem,
        ::AbstractRay, ::AbstractIntersection)
    amb = ambient_medium(system)
    return Transition(amb, amb, false)
end

function transition_for(::AbstractCoating, system::AbstractSystem,
        ray::AbstractRay, int::AbstractIntersection)
    amb = ambient_medium(system)
    return resolve_transition(amb, amb, ray, normal3d(int))
end

function transition_for(::AbstractObject, system::AbstractSystem,
        ::AbstractRay, ::AbstractIntersection)
    amb = ambient_medium(system)
    return Transition(amb, amb, false)
end

"""
    interaction_hint(optic::AbstractObject, trans::Transition, out_ray::AbstractRay, ray::AbstractRay) -> Nullable{Hint}

Returns an optional [`Hint`](@ref) to attach to the [`BeamInteraction`](@ref) produced by the
universal boundary interaction pipeline. Defaults to `nothing`; [`AbstractRefractiveOptic`](@ref)
overrides this to flag total internal reflection / re-entry ambiguities for the solver.

# Arguments

- `optic`: The optical component involved in the interaction.
- `trans`: The media transition across the boundary.
- `out_ray`: The outgoing ray after surface interaction.
- `ray`: The incoming ray before surface interaction.
"""
interaction_hint(::AbstractObject, ::Transition, ::AbstractRay, ::AbstractRay) = nothing

function interaction_hint(optic::AbstractRefractiveOptic, trans::Transition,
        out_ray::AbstractRay{T}, ray::AbstractRay{T}) where {T <: Real}
    entering = trans.is_entering
    optic_n = refractive_index(optic, wavelength(ray))
    tol = sqrt(eps(T))
    tir = (!entering && isapprox(refractive_index(out_ray), optic_n; atol = tol))
    return (entering || tir) ? Hint(optic) : nothing
end

"""
    interact3d(system::AbstractSystem, optic::Union{AbstractRefractiveOptic, AbstractReflectiveOptic, AbstractCoating}, beam::AbstractBeam{T}, ray::AbstractRay{T}) where {T <: Real}

Universal boundary interaction pipeline shared by refractive optics, reflective optics and
coatings: resolves the current intersection and media [`Transition`](@ref) (see
[`transition_for`](@ref)), evaluates the optic's `surface_model`, and packages the result as a
[`BeamInteraction`](@ref) (with an optional [`Hint`](@ref), see [`interaction_hint`](@ref)).

# Arguments

- `system`: The optical system enclosing the scene.
- `optic`: The target optical element.
- `beam`: The beam being traced.
- `ray`: The active ray segment.
"""
function interact3d(
        system::AbstractSystem,
        optic::Union{AbstractRefractiveOptic, AbstractReflectiveOptic, AbstractCoating},
        beam::AbstractBeam{T},
        ray::AbstractRay{T}
) where {T <: Real}
    int = current_intersection(beam, optic, ray)
    isnothing(int) && return nothing
    trans = transition_for(optic, system, ray, int)
    out_ray = interact3d(surface_model(optic), trans, int, ray)
    isnothing(out_ray) && return nothing
    return BeamInteraction(interaction_hint(optic, trans, out_ray, ray), out_ray)
end
