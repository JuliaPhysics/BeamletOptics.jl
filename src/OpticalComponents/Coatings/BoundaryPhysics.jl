# Shared Refractive and Reflective Boundary Interaction Helpers
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
    att_power, _ = bulk_attenuation_factor(n_incident, λ, length(ray))
    if TIR
        ndir = reflection3d(direction(ray), normal)
        n_out = n_incident
        hint = Hint(substrate_obj)
        T_coeff = 1.0
    else
        n_out = n_transmitted
        θi = angle3d(direction(ray), -normal)
        T_coeff = coating_transmittance(coating_model, θi, λ, n_incident, n_transmitted; from_front=from_front)
    end
    return BeamInteraction{T, R}(hint, Ray{T}(npos, ndir, nothing, λ, n_out, weight(ray) * att_power * T_coeff))
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
    θi = angle3d(direction(ray), -normal)
    R_coeff = coating_reflectance(coating_model, θi, λ, n_incident, n_transmitted; from_front=from_front)
    att_power, _ = bulk_attenuation_factor(n_incident, λ, length(ray))
    return BeamInteraction{T, R}(hint_out, Ray{T}(npos, ndir, nothing, λ, n_out, weight(ray) * att_power * R_coeff))
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

    θi = angle3d(direction(ray), -normal)
    R_coeff = coating_reflectance(coating_model, θi, λ, n_incident, n_transmitted; from_front=from_front)
    T_coeff = coating_transmittance(coating_model, θi, λ, n_incident, n_transmitted; from_front=from_front)
    att_power, _ = bulk_attenuation_factor(n_incident, λ, length(ray))

    beam_t = Beam(Ray{T}(npos, ndir_t, nothing, λ, n_transmitted, weight(ray) * att_power * T_coeff))
    beam_r = Beam(Ray{T}(npos, ndir_r, nothing, λ, n_incident, weight(ray) * att_power * R_coeff))
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
    _, att_field = bulk_attenuation_factor(n_incident, λ, length(ray))
    E0 = _calculate_global_E0(substrate_obj, ray, ndir, J) * att_field
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
    _, att_field = bulk_attenuation_factor(n_incident, λ, length(ray))
    E0 = _calculate_global_E0(substrate_obj, ray, ndir, J) * att_field
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

    _, att_field = bulk_attenuation_factor(n_incident, λ, length(ray))
    E0_t = _calculate_global_E0(substrate_obj, ray, ndir_t, J_t) * att_field
    E0_r = _calculate_global_E0(substrate_obj, ray, ndir_r, J_r) * att_field

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

@inline function _resolve_interface_context(
        system::AbstractSystem, substrate_obj::AbstractObject, ray::AbstractRay)
    int = intersection(ray)
    λ = wavelength(ray)
    n_substrate = is_refractive(substrate_obj) ? refractive_index(substrate_obj, λ) : refractive_index(system, λ)
    
    obj = object(int)
    is_substrate = (shape(obj) === shape(substrate_obj))
    
    normal = normal3d(int)
    if !is_substrate
        normal = -normal
    end
    
    entering_substrate = dot(direction(ray), normal) < 0
    from_front = entering_substrate

    if entering_substrate
        n_incident = refractive_index(ray)
        n_transmitted = n_substrate
        hint = is_refractive(substrate_obj) ? Hint(substrate_obj) : nothing
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
    return (normal, n_incident, n_transmitted, hint, λ, from_front)
end

# Main entry points for interact_refractive_boundary
function interact_refractive_boundary(
        system::AbstractSystem,
        substrate_obj::AbstractObject{T},
        coating_model,
        beam::AbstractBeam{T, R},
        ray::R
) where {T, R <: AbstractRay{T}}
    normal, n_incident, n_transmitted, hint, λ, from_front =
        _resolve_interface_context(system, substrate_obj, ray)
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
    normal, n_incident, n_transmitted, _, λ, from_front =
        _resolve_interface_context(system, substrate_obj, ray)
    
    npos = position(ray) + length(ray) * direction(ray)
    ndir = reflection3d(direction(ray), normal)

    J_r = if coating_model isa Uncoated
        SPBasis(-1, 0, 0, 1)
    else
        get_jones_matrix(coating_model, angle3d(direction(ray), -normal), λ, n_incident, n_transmitted, true; from_front=from_front)
    end
    R_coeff = clamp(unpolarized_reflectance(J_r), 0.0, 1.0)
    att_power, _ = bulk_attenuation_factor(n_incident, λ, length(ray))

    return BeamInteraction{T, R}(
        nothing, Ray{T}(npos, ndir, nothing, λ, n_incident, weight(ray) * att_power * R_coeff))
end

function interact_reflective_boundary(
        system::AbstractSystem,
        substrate_obj::AbstractObject{T},
        coating_model,
        beam::AbstractBeam{T, R},
        ray::R
) where {T, R <: PolarizedRay{T}}
    normal, n_incident, n_transmitted, _, λ, from_front =
        _resolve_interface_context(system, substrate_obj, ray)
    
    npos = position(ray) + length(ray) * direction(ray)
    ndir = reflection3d(direction(ray), normal)

    J = if coating_model isa Uncoated
        SPBasis(-1, 0, 0, 1)
    else
        get_jones_matrix(coating_model, angle3d(direction(ray), -normal),
            λ, n_incident, n_transmitted, true; from_front=from_front)
    end
    _, att_field = bulk_attenuation_factor(n_incident, λ, length(ray))
    E0 = _calculate_global_E0(substrate_obj, ray, ndir, J) * att_field
    return BeamInteraction{T, R}(
        nothing, PolarizedRay{T}(npos, ndir, nothing, λ, refractive_index(ray), E0))
end

# Absorptive behavior
function interact_refractive_boundary(
        ::Absorptive,
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
) where {T, R <: AbstractRay{T}}
    return nothing
end

# Shared splitting propagation helpers to eliminate code duplication
function _propagate_splitting_gaussian_beamlet(
    system::AbstractSystem,
    coated_or_coating::AbstractObject,
    coating_model,
    gauss::GaussianBeamlet{T},
    ray_id::Int,
    n_incident,
    n_transmitted,
    normal,
    from_front,
    tir_callback
) where {T}
    c_ray = gauss.chief.rays[ray_id]
    w_ray = gauss.waist.rays[ray_id]
    d_ray = gauss.divergence.rays[ray_id]

    c_dir_t, TIR_c = refraction3d(direction(c_ray), normal, n_incident, n_transmitted)
    w_dir_t, TIR_w = refraction3d(direction(w_ray), normal, n_incident, n_transmitted)
    d_dir_t, TIR_d = refraction3d(direction(d_ray), normal, n_incident, n_transmitted)
    if TIR_c || TIR_w || TIR_d
        return tir_callback()
    end

    c_dir_r = reflection3d(direction(c_ray), normal)
    w_dir_r = reflection3d(direction(w_ray), normal)
    d_dir_r = reflection3d(direction(d_ray), normal)

    pos = position(c_ray) + length(c_ray) * direction(c_ray)
    λ = wavelength(c_ray)

    J_t = get_jones_matrix(coating_model, angle3d(direction(c_ray), -normal),
        λ, n_incident, n_transmitted, false; from_front=from_front)
    J_r = get_jones_matrix(coating_model, angle3d(direction(c_ray), -normal),
        λ, n_incident, n_transmitted, true; from_front=from_front)

    if !iszero(J_t[1, 2]) || !iszero(J_t[2, 1]) || !iszero(J_r[1, 2]) || !iszero(J_r[2, 1])
        @warn "Jones matrix contains off-diagonal cross-polarization terms. Scalar GaussianBeamlet only tracks diagonal transmission J[1,1]. Use AstigmaticGaussianBeamlet for full vector polarization." maxlog = 1
    end

    _, att_field = bulk_attenuation_factor(n_incident, λ, length(c_ray))

    # Spawn transmitted
    chief_t = Beam(Ray{T}(pos, c_dir_t, nothing, λ, n_transmitted))
    waist_t = Beam(Ray{T}(
        position(w_ray) + length(w_ray) * direction(w_ray), w_dir_t, nothing, λ, n_transmitted))
    divergent_t = Beam(Ray{T}(
        position(d_ray) + length(d_ray) * direction(d_ray), d_dir_t, nothing, λ, n_transmitted))

    w0_t = gauss_parameters(gauss, length(gauss))[4]
    # NOTE: Using J_t[1,1] / J_r[1,1] is a polarization-independent approximation for plain GaussianBeamlets.
    # It ignores the full Jones matrix (e.g. cross-polarization and off-diagonal elements).
    # For full polarization support, AstigmaticGaussianBeamlet should be used instead.
    t_val = J_t[1, 1]
    E0_t = t_val * electric_field(gauss) * att_field * (beam_waist(gauss) / w0_t)
    t = GaussianBeamlet(chief_t, waist_t, divergent_t, wavelength(gauss), w0_t, E0_t)

    # Spawn reflected
    chief_r = Beam(Ray{T}(pos, c_dir_r, nothing, λ, n_incident))
    waist_r = Beam(Ray{T}(
        position(w_ray) + length(w_ray) * direction(w_ray), w_dir_r, nothing, λ, n_incident))
    divergent_r = Beam(Ray{T}(
        position(d_ray) + length(d_ray) * direction(d_ray), d_dir_r, nothing, λ, n_incident))

    w0_r = w0_t
    r_val = -J_r[1, 1]
    E0_r = r_val * electric_field(gauss) * att_field * (beam_waist(gauss) / w0_r)
    r = GaussianBeamlet(chief_r, waist_r, divergent_r, wavelength(gauss), w0_r, E0_r)

    children!(gauss, (t, r))
    return nothing
end

function _propagate_splitting_astigmatic_beamlet(
    system::AbstractSystem,
    coated_or_coating::AbstractObject,
    coating_model,
    agb::AstigmaticGaussianBeamlet{T},
    ray_id::Int,
    n_incident,
    n_transmitted,
    normal,
    from_front,
    tir_callback
) where {T}
    c_ray = rays(agb.c)[ray_id]
    λ = wavelength(c_ray)

    res_c   = refraction3d(direction(c_ray), normal, n_incident, n_transmitted)
    res_wxp = refraction3d(direction(rays(agb.wxp)[ray_id]), normal, n_incident, n_transmitted)
    res_wxm = refraction3d(direction(rays(agb.wxm)[ray_id]), normal, n_incident, n_transmitted)
    res_wyp = refraction3d(direction(rays(agb.wyp)[ray_id]), normal, n_incident, n_transmitted)
    res_wym = refraction3d(direction(rays(agb.wym)[ray_id]), normal, n_incident, n_transmitted)
    res_dxp = refraction3d(direction(rays(agb.dxp)[ray_id]), normal, n_incident, n_transmitted)
    res_dxm = refraction3d(direction(rays(agb.dxm)[ray_id]), normal, n_incident, n_transmitted)
    res_dyp = refraction3d(direction(rays(agb.dyp)[ray_id]), normal, n_incident, n_transmitted)
    res_dym = refraction3d(direction(rays(agb.dym)[ray_id]), normal, n_incident, n_transmitted)

    if res_c[2] || res_wxp[2] || res_wxm[2] || res_wyp[2] || res_wym[2] ||
       res_dxp[2] || res_dxm[2] || res_dyp[2] || res_dym[2]
        return tir_callback()
    end

    dirs_t = (
        res_c[1],
        res_wxp[1],
        res_wxm[1],
        res_wyp[1],
        res_wym[1],
        res_dxp[1],
        res_dxm[1],
        res_dyp[1],
        res_dym[1]
    )
    dirs_r = (
        reflection3d(
            direction(rays(agb.c)[ray_id]), normal3d(intersection(rays(agb.c)[ray_id]))),
        reflection3d(direction(rays(agb.wxp)[ray_id]),
            normal3d(intersection(rays(agb.wxp)[ray_id]))),
        reflection3d(direction(rays(agb.wxm)[ray_id]),
            normal3d(intersection(rays(agb.wxm)[ray_id]))),
        reflection3d(direction(rays(agb.wyp)[ray_id]),
            normal3d(intersection(rays(agb.wyp)[ray_id]))),
        reflection3d(direction(rays(agb.wym)[ray_id]),
            normal3d(intersection(rays(agb.wym)[ray_id]))),
        reflection3d(direction(rays(agb.dxp)[ray_id]),
            normal3d(intersection(rays(agb.dxp)[ray_id]))),
        reflection3d(direction(rays(agb.dxm)[ray_id]),
            normal3d(intersection(rays(agb.dxm)[ray_id]))),
        reflection3d(direction(rays(agb.dyp)[ray_id]),
            normal3d(intersection(rays(agb.dyp)[ray_id]))),
        reflection3d(direction(rays(agb.dym)[ray_id]),
            normal3d(intersection(rays(agb.dym)[ray_id])))
    )

    J_t = get_jones_matrix(coating_model, angle3d(direction(c_ray), -normal),
        λ, n_incident, n_transmitted, false; from_front=from_front)
    J_r = get_jones_matrix(coating_model, angle3d(direction(c_ray), -normal),
        λ, n_incident, n_transmitted, true; from_front=from_front)

    _, att_field = bulk_attenuation_factor(n_incident, λ, length(c_ray))
    E0_t = _calculate_global_E0(coated_or_coating, c_ray, dirs_t[1], J_t) * att_field
    E0_r = _calculate_global_E0(coated_or_coating, c_ray, dirs_r[1], J_r) * att_field

    # Spawn transmitted
    chief_t = Beam(PolarizedRay{T}(
        position(c_ray) + length(c_ray) * direction(c_ray), dirs_t[1], nothing, λ, n_transmitted, E0_t))
    wxp_t = Beam(Ray{T}(
        position(rays(agb.wxp)[ray_id]) +
        length(rays(agb.wxp)[ray_id]) * direction(rays(agb.wxp)[ray_id]),
        dirs_t[2], nothing, λ, n_transmitted))
    wxm_t = Beam(Ray{T}(
        position(rays(agb.wxm)[ray_id]) +
        length(rays(agb.wxm)[ray_id]) * direction(rays(agb.wxm)[ray_id]),
        dirs_t[3], nothing, λ, n_transmitted))
    wyp_t = Beam(Ray{T}(
        position(rays(agb.wyp)[ray_id]) +
        length(rays(agb.wyp)[ray_id]) * direction(rays(agb.wyp)[ray_id]),
        dirs_t[4], nothing, λ, n_transmitted))
    wym_t = Beam(Ray{T}(
        position(rays(agb.wym)[ray_id]) +
        length(rays(agb.wym)[ray_id]) * direction(rays(agb.wym)[ray_id]),
        dirs_t[5], nothing, λ, n_transmitted))
    dxp_t = Beam(Ray{T}(
        position(rays(agb.dxp)[ray_id]) +
        length(rays(agb.dxp)[ray_id]) * direction(rays(agb.dxp)[ray_id]),
        dirs_t[6], nothing, λ, n_transmitted))
    dxm_t = Beam(Ray{T}(
        position(rays(agb.dxm)[ray_id]) +
        length(rays(agb.dxm)[ray_id]) * direction(rays(agb.dxm)[ray_id]),
        dirs_t[7], nothing, λ, n_transmitted))
    dyp_t = Beam(Ray{T}(
        position(rays(agb.dyp)[ray_id]) +
        length(rays(agb.dyp)[ray_id]) * direction(rays(agb.dyp)[ray_id]),
        dirs_t[8], nothing, λ, n_transmitted))
    dym_t = Beam(Ray{T}(
        position(rays(agb.dym)[ray_id]) +
        length(rays(agb.dym)[ray_id]) * direction(rays(agb.dym)[ray_id]),
        dirs_t[9], nothing, λ, n_transmitted))

    t = AstigmaticGaussianBeamlet(
        chief_t, wxp_t, wxm_t, wyp_t, wym_t, dxp_t, dxm_t, dyp_t, dym_t)

    # Spawn reflected
    chief_r = Beam(PolarizedRay{T}(
        position(c_ray) + length(c_ray) * direction(c_ray), dirs_r[1], nothing, λ, n_incident, E0_r))
    wxp_r = Beam(Ray{T}(
        position(rays(agb.wxp)[ray_id]) +
        length(rays(agb.wxp)[ray_id]) * direction(rays(agb.wxp)[ray_id]),
        dirs_r[2], nothing, λ, n_incident))
    wxm_r = Beam(Ray{T}(
        position(rays(agb.wxm)[ray_id]) +
        length(rays(agb.wxm)[ray_id]) * direction(rays(agb.wxm)[ray_id]),
        dirs_r[3], nothing, λ, n_incident))
    wyp_r = Beam(Ray{T}(
        position(rays(agb.wyp)[ray_id]) +
        length(rays(agb.wyp)[ray_id]) * direction(rays(agb.wyp)[ray_id]),
        dirs_r[4], nothing, λ, n_incident))
    wym_r = Beam(Ray{T}(
        position(rays(agb.wym)[ray_id]) +
        length(rays(agb.wym)[ray_id]) * direction(rays(agb.wym)[ray_id]),
        dirs_r[5], nothing, λ, n_incident))
    dxp_r = Beam(Ray{T}(
        position(rays(agb.dxp)[ray_id]) +
        length(rays(agb.dxp)[ray_id]) * direction(rays(agb.dxp)[ray_id]),
        dirs_r[6], nothing, λ, n_incident))
    dxm_r = Beam(Ray{T}(
        position(rays(agb.dxm)[ray_id]) +
        length(rays(agb.dxm)[ray_id]) * direction(rays(agb.dxm)[ray_id]),
        dirs_r[7], nothing, λ, n_incident))
    dyp_r = Beam(Ray{T}(
        position(rays(agb.dyp)[ray_id]) +
        length(rays(agb.dyp)[ray_id]) * direction(rays(agb.dyp)[ray_id]),
        dirs_r[8], nothing, λ, n_incident))
    dym_r = Beam(Ray{T}(
        position(rays(agb.dym)[ray_id]) +
        length(rays(agb.dym)[ray_id]) * direction(rays(agb.dym)[ray_id]),
        dirs_r[9], nothing, λ, n_incident))

    r = AstigmaticGaussianBeamlet(
        chief_r, wxp_r, wxm_r, wyp_r, wym_r, dxp_r, dxm_r, dyp_r, dym_r)

    children!(agb, (t, r))
    return nothing
end
