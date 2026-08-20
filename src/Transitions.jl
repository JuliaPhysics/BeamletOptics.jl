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
    # QUESTION is this method intended primarily for medium -> air and vice versa?
    is_entering = dot(direction(ray), normal) < 0
    med_in  = is_entering ? current_medium(ray) : component_medium
    med_out = is_entering ? component_medium : ambient
    return Transition(med_in, med_out, is_entering)
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
    # QUESTION two touching objects, entry + exit + possibly "coating" 
    if ex_obj !== nothing && en_obj !== nothing
        # QUESTION `medium_from` called for obj with no defined medium, e.g. Mirror?
        med_in = medium_from(ex_obj)
        med_out = medium_from(en_obj)
        return Transition(med_in, med_out, true)
    # QUESTION only exiting object + possibly "coating" 
    elseif ex_obj !== nothing
        return Transition(medium_from(ex_obj), ambient, false)
    # QUESTIOn only entering object + possibly "coating" 
    elseif en_obj !== nothing
        return Transition(ambient, medium_from(en_obj), true)
    # QUESTION which case is represented here?
    else
        return resolve_transition(Ambient(), ambient, ray, normal3d(mi))
    end
end

# QUESTION this is the default behavior, correct? extensions considered in the future?
# Fallback for rays storing refractive index instead of AbstractMedium
current_medium(ray::AbstractRay) = IsotropicMedium(refractive_index(ray))
