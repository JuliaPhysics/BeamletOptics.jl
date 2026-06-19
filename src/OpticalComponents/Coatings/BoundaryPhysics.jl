# Shared Refractive and Reflective Boundary Physics Helpers

# Shared Refractive Boundary Interaction Helper
# Dispatched helpers for interact_refractive_boundary by CoatingBehavior
# For Ray{T}:
function interact_refractive_boundary(
        ::Transmissive,
        system::AbstractSystem,
        substrate_obj::AbstractObject{T},
        coating_model,
        beam::AbstractBeam{T, R},
        ray::R,
        n_incident,
        n_transmitted,
        hint,
        normal,
        λ,
        from_front
) where {T, R <: Ray{T}}
    ndir, TIR = refraction3d(direction(ray), normal, n_incident, n_transmitted)
    npos = position(ray) + length(ray) * direction(ray)
    if TIR
        ndir = reflection3d(direction(ray), normal)
        n_out = n_incident
        hint = Hint(substrate_obj)
        T_coeff = 1.0
    else
        n_out = n_transmitted
        J_r = get_jones_matrix(coating_model, angle3d(direction(ray), -normal), λ, n_incident, n_transmitted, true; from_front=from_front)
        R_coeff = 0.5 * (abs2(J_r[1,1]) + abs2(J_r[1,2]) + abs2(J_r[2,1]) + abs2(J_r[2,2]))
        T_coeff = 1.0 - R_coeff
    end
    return BeamInteraction{T, R}(hint, Ray{T}(npos, ndir, nothing, λ, n_out, weight(ray) * T_coeff))
end

function interact_refractive_boundary(
        ::Reflective,
        system::AbstractSystem,
        substrate_obj::AbstractObject{T},
        coating_model,
        beam::AbstractBeam{T, R},
        ray::R,
        n_incident,
        n_transmitted,
        hint,
        normal,
        λ,
        from_front
) where {T, R <: Ray{T}}
    ndir = reflection3d(direction(ray), normal)
    npos = position(ray) + length(ray) * direction(ray)
    n_out = n_incident
    hint_out = isentering(ray) ? nothing : Hint(substrate_obj)
    J_r = get_jones_matrix(coating_model, angle3d(direction(ray), -normal), λ, n_incident, n_transmitted, true; from_front=from_front)
    R_coeff = 0.5 * (abs2(J_r[1,1]) + abs2(J_r[1,2]) + abs2(J_r[2,1]) + abs2(J_r[2,2]))
    return BeamInteraction{T, R}(hint_out, Ray{T}(npos, ndir, nothing, λ, n_out, weight(ray) * R_coeff))
end

function interact_refractive_boundary(
        ::Splitting,
        system::AbstractSystem,
        substrate_obj::AbstractObject{T},
        coating_model,
        beam::AbstractBeam{T, R},
        ray::R,
        n_incident,
        n_transmitted,
        hint,
        normal,
        λ,
        from_front
) where {T, R <: Ray{T}}
    ndir_t, TIR = refraction3d(direction(ray), normal, n_incident, n_transmitted)
    if TIR
        return interact_refractive_boundary(
            Reflective(), system, substrate_obj, coating_model, beam,
            ray, n_incident, n_transmitted, hint, normal, λ, from_front)
    end
    ndir_r = reflection3d(direction(ray), normal)
    npos = position(ray) + length(ray) * direction(ray)

    J_r = get_jones_matrix(coating_model, angle3d(direction(ray), -normal), λ, n_incident, n_transmitted, true; from_front=from_front)
    R_coeff = 0.5 * (abs2(J_r[1,1]) + abs2(J_r[1,2]) + abs2(J_r[2,1]) + abs2(J_r[2,2]))
    T_coeff = 1.0 - R_coeff

    beam_t = Beam(Ray{T}(npos, ndir_t, nothing, λ, n_transmitted, weight(ray) * T_coeff))
    beam_r = Beam(Ray{T}(npos, ndir_r, nothing, λ, n_incident, weight(ray) * R_coeff))
    children!(beam, (beam_t, beam_r))
    return nothing
end

# For PolarizedRay{T}:
function interact_refractive_boundary(
        ::Transmissive,
        system::AbstractSystem,
        substrate_obj::AbstractObject{T},
        coating_model,
        beam::AbstractBeam{T, R},
        ray::R,
        n_incident,
        n_transmitted,
        hint,
        normal,
        λ,
        from_front
) where {T, R <: PolarizedRay{T}}
    ndir, TIR = refraction3d(direction(ray), normal, n_incident, n_transmitted)
    npos = position(ray) + length(ray) * direction(ray)
    if TIR
        ndir = reflection3d(direction(ray), normal)
        n_out = n_incident
        hint = Hint(substrate_obj)
        J = get_jones_matrix(coating_model, angle3d(direction(ray), -normal),
            λ, n_incident, n_transmitted, true; from_front=from_front)
    else
        n_out = n_transmitted
        J = get_jones_matrix(coating_model, angle3d(direction(ray), -normal),
            λ, n_incident, n_transmitted, false; from_front=from_front)
    end
    E0 = _calculate_global_E0(substrate_obj, ray, ndir, J)
    return BeamInteraction{T, R}(hint, PolarizedRay{T}(npos, ndir, nothing, λ, n_out, E0))
end

function interact_refractive_boundary(
        ::Reflective,
        system::AbstractSystem,
        substrate_obj::AbstractObject{T},
        coating_model,
        beam::AbstractBeam{T, R},
        ray::R,
        n_incident,
        n_transmitted,
        hint,
        normal,
        λ,
        from_front
) where {T, R <: PolarizedRay{T}}
    ndir = reflection3d(direction(ray), normal)
    npos = position(ray) + length(ray) * direction(ray)
    n_out = n_incident
    hint_out = isentering(ray) ? nothing : Hint(substrate_obj)
    J = get_jones_matrix(
        coating_model, angle3d(direction(ray), -normal), λ, n_incident, n_transmitted, true; from_front=from_front)
    E0 = _calculate_global_E0(substrate_obj, ray, ndir, J)
    return BeamInteraction{T, R}(
        hint_out, PolarizedRay{T}(npos, ndir, nothing, λ, n_out, E0))
end

function interact_refractive_boundary(
        ::Splitting,
        system::AbstractSystem,
        substrate_obj::AbstractObject{T},
        coating_model,
        beam::AbstractBeam{T, R},
        ray::R,
        n_incident,
        n_transmitted,
        hint,
        normal,
        λ,
        from_front
) where {T, R <: PolarizedRay{T}}
    ndir_t, TIR = refraction3d(direction(ray), normal, n_incident, n_transmitted)
    if TIR
        return interact_refractive_boundary(
            Reflective(), system, substrate_obj, coating_model, beam,
            ray, n_incident, n_transmitted, hint, normal, λ, from_front)
    end
    ndir_r = reflection3d(direction(ray), normal)
    npos = position(ray) + length(ray) * direction(ray)

    J_t = get_jones_matrix(coating_model, angle3d(direction(ray), -normal),
        λ, n_incident, n_transmitted, false; from_front=from_front)
    J_r = get_jones_matrix(
        coating_model, angle3d(direction(ray), -normal), λ, n_incident, n_transmitted, true; from_front=from_front)

    E0_t = _calculate_global_E0(substrate_obj, ray, ndir_t, J_t)
    E0_r = _calculate_global_E0(substrate_obj, ray, ndir_r, J_r)

    beam_t = Beam(PolarizedRay{T}(npos, ndir_t, nothing, λ, n_transmitted, E0_t))
    beam_r = Beam(PolarizedRay{T}(npos, ndir_r, nothing, λ, n_incident, E0_r))
    children!(beam, (beam_t, beam_r))
    return nothing
end

"""
    _base_optic(obj)

Helper to extract the underlying optic from an object wrapper for coincident boundary evaluation.
Overridden by `CoatedComponent` wrappers.
"""
_base_optic(obj) = obj

# Main entry points for interact_refractive_boundary
function interact_refractive_boundary(
        system::AbstractSystem,
        substrate_obj::AbstractObject{T},
        coating_model,
        beam::AbstractBeam{T, R},
        ray::R
) where {T, R <: AbstractRay{T}}
    int = intersection(ray)
    λ = wavelength(ray)
    n_substrate = refractive_index(substrate_obj, λ)
    
    obj = object(int)
    is_substrate = (_base_optic(obj) === substrate_obj)
    
    normal = normal3d(int)
    if !is_substrate
        normal = -normal
    end
    
    entering_substrate = dot(direction(ray), normal) < 0
    from_front = entering_substrate

    if entering_substrate
        n_incident = refractive_index(ray)
        n_transmitted = n_substrate
        hint = Hint(substrate_obj)
    else
        n_incident = n_substrate
        coin_obj = is_substrate ? int.coincident_object : obj
        if !isnothing(coin_obj) && is_refractive(coin_obj)
            n_transmitted = refractive_index(coin_obj, λ)
            hint = Hint(coin_obj)
        else
            n_transmitted = refractive_index(system, λ)
            hint = nothing
        end
        normal = -normal
    end

    behavior = coating_behavior(coating_model, ray)
    return interact_refractive_boundary(
        behavior, system, substrate_obj, coating_model, beam,
        ray, n_incident, n_transmitted, hint, normal, λ, from_front)
end

# Shared Reflective Boundary Interaction Helper
function interact_reflective_boundary(
        system::AbstractSystem,
        substrate_obj::AbstractObject{T},
        coating_model,
        beam::AbstractBeam{T, R},
        ray::R
) where {T, R <: Ray{T}}
    int = intersection(ray)
    obj = object(int)
    is_substrate = (_base_optic(obj) === substrate_obj)
    
    normal = normal3d(int)
    if !is_substrate
        normal = -normal
    end
    
    entering_substrate = dot(direction(ray), normal) < 0
    from_front = entering_substrate
    
    npos = position(ray) + length(ray) * direction(ray)
    ndir = reflection3d(direction(ray), normal)

    λ = wavelength(ray)
    n_incident = refractive_index(ray)
    n_transmitted = n_incident
    
    coin_obj = is_substrate ? int.coincident_object : obj
    if !isnothing(coin_obj) && is_refractive(coin_obj)
        n_transmitted = refractive_index(coin_obj, λ)
    else
        n_transmitted = refractive_index(system, λ)
    end

    J_r = if coating_model isa Uncoated
        SPBasis(-1, 0, 0, 1)
    else
        get_jones_matrix(coating_model, angle3d(direction(ray), -normal), λ, n_incident, n_transmitted, true; from_front=from_front)
    end
    R_coeff = 0.5 * (abs2(J_r[1,1]) + abs2(J_r[1,2]) + abs2(J_r[2,1]) + abs2(J_r[2,2]))

    return BeamInteraction{T, R}(
        nothing, Ray{T}(npos, ndir, nothing, λ, n_incident, weight(ray) * R_coeff))
end

function interact_reflective_boundary(
        system::AbstractSystem,
        substrate_obj::AbstractObject{T},
        coating_model,
        beam::AbstractBeam{T, R},
        ray::R
) where {T, R <: PolarizedRay{T}}
    int = intersection(ray)
    obj = object(int)
    is_substrate = (_base_optic(obj) === substrate_obj)
    
    normal = normal3d(int)
    if !is_substrate
        normal = -normal
    end
    
    entering_substrate = dot(direction(ray), normal) < 0
    from_front = entering_substrate
    
    npos = position(ray) + length(ray) * direction(ray)
    ndir = reflection3d(direction(ray), normal)

    λ = wavelength(ray)
    n_incident = refractive_index(ray)
    n_transmitted = n_incident
    
    coin_obj = is_substrate ? int.coincident_object : obj
    if !isnothing(coin_obj) && is_refractive(coin_obj)
        n_transmitted = refractive_index(coin_obj, λ)
    else
        n_transmitted = refractive_index(system, λ)
    end

    J = if coating_model isa Uncoated
        SPBasis(-1, 0, 0, 1)
    else
        get_jones_matrix(coating_model, angle3d(direction(ray), -normal),
            λ, n_incident, n_transmitted, true; from_front=from_front)
    end
    E0 = _calculate_global_E0(substrate_obj, ray, ndir, J)
    return BeamInteraction{T, R}(
        nothing, PolarizedRay{T}(npos, ndir, nothing, λ, refractive_index(ray), E0))
end
