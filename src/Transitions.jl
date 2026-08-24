"""
    Transition{Min<:AbstractMedium, Mout<:AbstractMedium}

Carries media transition state across an optical boundary.

# Fields
- `medium_in::Min`: Incident medium
- `medium_out::Mout`: Transmitted medium
- `is_entering::Bool`: `true` if entering the component, `false` if exiting
"""
struct Transition{Min <: AbstractMedium, Mout <: AbstractMedium}
    medium_in::Min
    medium_out::Mout
    is_entering::Bool
end

"""
    resolve_transition(component_medium, ambient, ray, normal) -> Transition

Resolves the incident and transmitted media across a boundary given the component medium,
system ambient medium, ray propagation direction, and boundary outward surface normal.
"""
@inline function resolve_transition(
        component_medium::AbstractMedium,
        ambient::AbstractMedium,
        ray::AbstractRay,
        normal::Point3
)
    # Media transition for an isolated component immersed in ambient medium
    is_entering = dot(direction(ray), normal) < 0
    if is_entering
        return Transition(ambient, component_medium, true)
    else
        return Transition(component_medium, ambient, false)
    end
end

"""
    resolve_transition(mi::MultiIntersection, ambient, ray) -> Transition

Resolves the media transition across a coincident boundary between touching optics.
"""
@inline function resolve_transition(
        mi::MultiIntersection,
        ambient::AbstractMedium,
        ray::AbstractRay
)
    ex_obj = exiting(mi)
    en_obj = entering(mi)
    if ex_obj !== nothing && en_obj !== nothing
        # Direct interface between two touching optical components (e.g. cemented doublet)
        med_in = medium_from(ex_obj)
        med_out = medium_from(en_obj)
        return Transition(med_in, med_out, true)
    elseif ex_obj !== nothing
        # Ray emerging from an optical component into ambient medium
        return Transition(medium_from(ex_obj), ambient, false)
    elseif en_obj !== nothing
        # Ray entering an optical component from ambient medium
        return Transition(ambient, medium_from(en_obj), true)
    else
        # Standalone virtual interface or detector plane in ambient medium
        return resolve_transition(Ambient(), ambient, ray, normal3d(mi))
    end
end

# Fallback for rays storing scalar refractive index instead of AbstractMedium
current_medium(ray::AbstractRay) = IsotropicMedium(refractive_index(ray))
