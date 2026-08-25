"""
    current_intersection(beam, obj, ray) -> Nullable{AbstractIntersection}

Returns the surface intersection for the current interaction step: the intersection
already carried by `ray`, or by the most recently traced ray in `beam`, or a freshly computed one via
`intersect3d` against `obj`'s shape if `ray` has not been resolved yet (e.g. when `interact3d` is called
directly with an unresolved ray instead of via [`trace_system!`](@ref)).
"""
function current_intersection(beam::AbstractBeam, obj, ray::AbstractRay)
    ray_int = intersection(ray)
    !isnothing(ray_int) && return ray_int
    r_vec = rays(beam)
    last_int = isempty(r_vec) ? nothing : intersection(last(r_vec))
    return isnothing(last_int) ? intersect3d(shape(obj), ray) : last_int
end

"""
    transition_for(optic, system, ray, int) -> Transition

Resolves the media [`Transition`](@ref) across `optic`'s boundary for the universal
boundary interaction pipeline (see [`interact3d(::AbstractSystem, ::Union, ::AbstractBeam, ::AbstractRay)`](@ref)).
Subtypes of [`AbstractRefractiveOptic`](@ref), [`AbstractReflectiveOptic`](@ref) and
[`AbstractCoating`](@ref) may specialize this to override the default per-family behavior.
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

function transition_for(optic::AbstractObject, system::AbstractSystem,
        ray::AbstractRay, int::AbstractIntersection)
    is_refractive(optic) ?
    resolve_transition(medium_from(optic), ambient_medium(system), ray, normal3d(int)) :
    Transition(ambient_medium(system), ambient_medium(system), false)
end

"""
    interaction_hint(optic, trans, out_ray, ray) -> Nullable{Hint}

Returns an optional [`Hint`](@ref) to attach to the [`BeamInteraction`](@ref) produced by the
universal boundary interaction pipeline. Defaults to `nothing`; [`AbstractRefractiveOptic`](@ref)
overrides this to flag total internal reflection / re-entry ambiguities for the solver.
"""
interaction_hint(::AbstractObject, ::Transition, ::AbstractRay, ::AbstractRay) = nothing

function interaction_hint(optic::AbstractRefractiveOptic, trans::Transition,
        out_ray::AbstractRay, ray::AbstractRay)
    entering = trans.is_entering
    optic_n = refractive_index(optic, wavelength(ray))
    tir = (!entering && isapprox(refractive_index(out_ray), optic_n; atol = 1e-10))
    return (entering || tir) ? Hint(optic) : nothing
end

"""
    interact3d(system, optic::Union{AbstractRefractiveOptic,AbstractReflectiveOptic,AbstractCoating}, beam, ray)

Universal boundary interaction pipeline shared by refractive optics, reflective optics and
coatings: resolves the current intersection and media [`Transition`](@ref) (see
[`transition_for`](@ref)), evaluates the optic's `surface_model`, and packages the result as a
[`BeamInteraction`](@ref) (with an optional [`Hint`](@ref), see [`interaction_hint`](@ref)).
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
