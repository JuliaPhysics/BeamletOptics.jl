# Standalone Coating struct and ray/beamlet interaction methods

# Coating object wrapping a shape and a model
"""
    Coating(shape, model; side = :either)
    Coating{T, S, M} <: AbstractCoating{T}

A standalone optical component representing an infinitesimally thin coated boundary floating in space.
Wraps an `AbstractShape` and a coating model `M` (e.g. `ThinFilmCoating`, `SimpleBeamsplitterCoating`).
Can be configured to only interact with light hitting a specific `side` (e.g. `:front` or `:back`), passing straight through from the other direction without effect.
"""
struct Coating{T <: Real, S <: AbstractShape{T}, M} <: AbstractCoating{T}
    shape::S
    model::M
    normal_filter::Union{Nothing, SVector{3, Float64}}
end

function Coating(shape::AbstractShape{T}, model, side::Symbol) where {T}
    normal_filter = if side == :front
        SVector{3, Float64}(0.0, -1.0, 0.0)
    elseif side == :back
        SVector{3, Float64}(0.0, 1.0, 0.0)
    else
        nothing
    end
    return Coating{T, typeof(shape), typeof(model)}(shape, model, normal_filter)
end

function Coating(shape::AbstractShape{T}, model,
        normal_filter::Union{Nothing, AbstractVector}) where {T}
    nf = isnothing(normal_filter) ? nothing : SVector{3, Float64}(normal_filter)
    return Coating{T, typeof(shape), typeof(model)}(shape, model, nf)
end

function Coating(shape::AbstractShape{T}, model; side::Symbol = :either,
        normal_filter = nothing) where {T}
    nf = if !isnothing(normal_filter)
        SVector{3, Float64}(normal_filter)
    elseif side == :front
        SVector{3, Float64}(0.0, -1.0, 0.0)
    elseif side == :back
        SVector{3, Float64}(0.0, 1.0, 0.0)
    else
        nothing
    end
    return Coating{T, typeof(shape), typeof(model)}(shape, model, nf)
end

shape_trait_of(::Coating) = SingleShape()
shape(c::Coating) = c.shape

function intersect3d(::SingleShape, coating::Coating, ray::AbstractRay)
    int = intersect3d(shape(coating), ray)
    isnothing(int) && return nothing

    if !isnothing(coating.normal_filter)
        wn = normal3d(int)
        parent_shape = shape(coating)
        local_n = transposed_orientation(parent_shape) * wn
        if dot(local_n, coating.normal_filter) <= 0.5
            return nothing
        end
    end

    object!(int, coating)
    return int
end

@inline function _coating_media(system::AbstractSystem, ray::AbstractRay, int::Intersection)
    λ = wavelength(ray)
    n_incident = refractive_index(ray)
    n_sys = refractive_index(system, λ)

    n_exit = (!isnothing(int.coincident_object) && is_refractive(int.coincident_object)) ?
             refractive_index(int.coincident_object, λ) : n_sys
    n_enter = (!isnothing(int.coincident_object_2) &&
               is_refractive(int.coincident_object_2)) ?
              refractive_index(int.coincident_object_2, λ) : n_sys

    tol = get_index_matching_tolerance()
    if abs(n_incident - n_exit) < tol
        n_transmitted = n_enter
        n_reflected = n_exit
    else
        n_transmitted = n_exit
        n_reflected = n_enter
    end
    return n_transmitted, n_reflected
end

# Route interact3d through coating_behavior trait
function interact3d(system::AbstractSystem, coating::Coating{T},
        beam::AbstractBeam, ray::AbstractRay) where {T}
    return interact3d(coating_behavior(coating.model, ray), system, coating, beam, ray)
end

# Transmissive behavior
function interact3d(::Transmissive, system::AbstractSystem, coating::Coating{T},
        ::Beam{T, R}, ray::R) where {T, R <: Ray{T}}
    int = intersection(ray)
    n_transmitted, _ = _coating_media(system, ray, int)

    normal = normal3d(int)
    if dot(direction(ray), normal) > 0
        normal = -normal
    end

    ndir, TIR = refraction3d(direction(ray), normal, refractive_index(ray), n_transmitted)
    npos = position(ray) + length(ray) * direction(ray)

    if TIR
        ndir = reflection3d(direction(ray), normal)
        n_out = refractive_index(ray)

        tol = get_index_matching_tolerance()
        n_exit = isnothing(int.coincident_object) ? nothing : int.coincident_object
        n_enter = isnothing(int.coincident_object_2) ? nothing : int.coincident_object_2

        λ = wavelength(ray)
        n_sys = refractive_index(system, λ)
        n_exit_val = isnothing(n_exit) ? n_sys :
                     (is_refractive(n_exit) ? refractive_index(n_exit, λ) : n_sys)
        incident_obj = abs(refractive_index(ray) - n_exit_val) < tol ? n_exit : n_enter
        hint = isnothing(incident_obj) ? nothing : Hint(incident_obj)
        T_coeff = 1.0
    else
        n_out = n_transmitted

        tol = get_index_matching_tolerance()
        n_exit = isnothing(int.coincident_object) ? nothing : int.coincident_object
        n_enter = isnothing(int.coincident_object_2) ? nothing : int.coincident_object_2

        λ = wavelength(ray)
        n_sys = refractive_index(system, λ)
        n_exit_val = isnothing(n_exit) ? n_sys :
                     (is_refractive(n_exit) ? refractive_index(n_exit, λ) : n_sys)
        entered_obj = abs(refractive_index(ray) - n_exit_val) < tol ? n_enter : n_exit
        hint = isnothing(entered_obj) ? nothing : Hint(entered_obj)

        J_r = get_jones_matrix(coating.model, angle3d(direction(ray), -normal),
            wavelength(ray), refractive_index(ray), n_transmitted, true)
        R_coeff = 0.5 *
                  (abs2(J_r[1, 1]) + abs2(J_r[1, 2]) + abs2(J_r[2, 1]) + abs2(J_r[2, 2]))
        T_coeff = 1.0 - R_coeff
    end

    return BeamInteraction{T, R}(
        hint, Ray{T}(npos, ndir, nothing, wavelength(ray), n_out, weight(ray) * T_coeff))
end

function interact3d(::Transmissive, system::AbstractSystem, coating::Coating{T},
        ::Beam{T, R}, ray::R) where {T, R <: PolarizedRay{T}}
    int = intersection(ray)
    n_transmitted, _ = _coating_media(system, ray, int)

    normal = normal3d(int)
    if dot(direction(ray), normal) > 0
        normal = -normal
    end

    ndir, TIR = refraction3d(direction(ray), normal, refractive_index(ray), n_transmitted)
    npos = position(ray) + length(ray) * direction(ray)

    if TIR
        ndir = reflection3d(direction(ray), normal)
        n_out = refractive_index(ray)

        tol = get_index_matching_tolerance()
        n_exit = isnothing(int.coincident_object) ? nothing : int.coincident_object
        n_enter = isnothing(int.coincident_object_2) ? nothing : int.coincident_object_2

        λ = wavelength(ray)
        n_sys = refractive_index(system, λ)
        n_exit_val = isnothing(n_exit) ? n_sys :
                     (is_refractive(n_exit) ? refractive_index(n_exit, λ) : n_sys)
        incident_obj = abs(refractive_index(ray) - n_exit_val) < tol ? n_exit : n_enter
        hint = isnothing(incident_obj) ? nothing : Hint(incident_obj)

        J = get_jones_matrix(coating.model, angle3d(direction(ray), -normal),
            wavelength(ray), refractive_index(ray), n_transmitted, true)
    else
        n_out = n_transmitted

        tol = get_index_matching_tolerance()
        n_exit = isnothing(int.coincident_object) ? nothing : int.coincident_object
        n_enter = isnothing(int.coincident_object_2) ? nothing : int.coincident_object_2

        λ = wavelength(ray)
        n_sys = refractive_index(system, λ)
        n_exit_val = isnothing(n_exit) ? n_sys :
                     (is_refractive(n_exit) ? refractive_index(n_exit, λ) : n_sys)
        entered_obj = abs(refractive_index(ray) - n_exit_val) < tol ? n_enter : n_exit
        hint = isnothing(entered_obj) ? nothing : Hint(entered_obj)

        J = get_jones_matrix(coating.model, angle3d(direction(ray), -normal),
            wavelength(ray), refractive_index(ray), n_transmitted, false)
    end

    E0 = _calculate_global_E0(coating, ray, ndir, J)
    return BeamInteraction{T, R}(
        hint, PolarizedRay{T}(npos, ndir, nothing, wavelength(ray), n_out, E0))
end

# Reflective behavior
function interact3d(::Reflective, system::AbstractSystem, coating::Coating{T},
        ::Beam{T, R}, ray::R) where {T, R <: Ray{T}}
    int = intersection(ray)
    normal = normal3d(int)
    if dot(direction(ray), normal) > 0
        normal = -normal
    end

    ndir = reflection3d(direction(ray), normal)
    npos = position(ray) + length(ray) * direction(ray)

    tol = get_index_matching_tolerance()
    n_exit = isnothing(int.coincident_object) ? nothing : int.coincident_object
    n_enter = isnothing(int.coincident_object_2) ? nothing : int.coincident_object_2

    λ = wavelength(ray)
    n_sys = refractive_index(system, λ)
    n_exit_val = isnothing(n_exit) ? n_sys :
                 (is_refractive(n_exit) ? refractive_index(n_exit, λ) : n_sys)
    incident_obj = abs(refractive_index(ray) - n_exit_val) < tol ? n_exit : n_enter
    hint = isnothing(incident_obj) ? nothing : Hint(incident_obj)

    n_transmitted, _ = _coating_media(system, ray, int)
    J_r = get_jones_matrix(coating.model, angle3d(direction(ray), -normal),
        wavelength(ray), refractive_index(ray), n_transmitted, true)
    R_coeff = 0.5 * (abs2(J_r[1, 1]) + abs2(J_r[1, 2]) + abs2(J_r[2, 1]) + abs2(J_r[2, 2]))

    return BeamInteraction{T, R}(
        hint, Ray{T}(npos, ndir, nothing, wavelength(ray),
            refractive_index(ray), weight(ray) * R_coeff))
end

function interact3d(::Reflective, system::AbstractSystem, coating::Coating{T},
        ::Beam{T, R}, ray::R) where {T, R <: PolarizedRay{T}}
    int = intersection(ray)
    normal = normal3d(int)
    if dot(direction(ray), normal) > 0
        normal = -normal
    end

    ndir = reflection3d(direction(ray), normal)
    npos = position(ray) + length(ray) * direction(ray)

    tol = get_index_matching_tolerance()
    n_exit = isnothing(int.coincident_object) ? nothing : int.coincident_object
    n_enter = isnothing(int.coincident_object_2) ? nothing : int.coincident_object_2

    λ = wavelength(ray)
    n_sys = refractive_index(system, λ)
    n_exit_val = isnothing(n_exit) ? n_sys :
                 (is_refractive(n_exit) ? refractive_index(n_exit, λ) : n_sys)
    incident_obj = abs(refractive_index(ray) - n_exit_val) < tol ? n_exit : n_enter
    hint = isnothing(incident_obj) ? nothing : Hint(incident_obj)

    n_transmitted, _ = _coating_media(system, ray, int)

    J = get_jones_matrix(coating.model, angle3d(direction(ray), -normal),
        wavelength(ray), refractive_index(ray), n_transmitted, true)
    E0 = _calculate_global_E0(coating, ray, ndir, J)
    return BeamInteraction{T, R}(hint,
        PolarizedRay{T}(npos, ndir, nothing, wavelength(ray), refractive_index(ray), E0))
end

# Splitting behavior
function interact3d(::Splitting, system::AbstractSystem, coating::Coating{T},
        beam::Beam{T, R}, ray::R) where {T, R <: Ray{T}}
    int = intersection(ray)
    n_transmitted, n_reflected = _coating_media(system, ray, int)

    normal = normal3d(int)
    if dot(direction(ray), normal) > 0
        normal = -normal
    end
    dir_t, TIR = refraction3d(direction(ray), normal, refractive_index(ray), n_transmitted)
    if TIR
        return interact3d(Reflective(), system, coating, beam, ray)
    end
    dir_r = reflection3d(direction(ray), normal)

    pos = position(ray) + length(ray) * direction(ray)

    J_r = get_jones_matrix(coating.model, angle3d(direction(ray), -normal),
        wavelength(ray), refractive_index(ray), n_transmitted, true)
    R_coeff = 0.5 * (abs2(J_r[1, 1]) + abs2(J_r[1, 2]) + abs2(J_r[2, 1]) + abs2(J_r[2, 2]))
    T_coeff = 1.0 - R_coeff

    beam_t = Beam(Ray{T}(
        Point3{T}(pos), Point3{T}(dir_t), nothing, wavelength(ray), n_transmitted, weight(ray) *
                                                                                   T_coeff))
    beam_r = Beam(Ray{T}(
        Point3{T}(pos), Point3{T}(dir_r), nothing, wavelength(ray), n_reflected, weight(ray) *
                                                                                 R_coeff))

    children!(beam, (beam_t, beam_r))
    return nothing
end

function interact3d(::Splitting, system::AbstractSystem, coating::Coating{T},
        beam::Beam{T, R}, ray::R) where {T, R <: PolarizedRay{T}}
    int = intersection(ray)
    n_transmitted, n_reflected = _coating_media(system, ray, int)

    normal = normal3d(int)
    if dot(direction(ray), normal) > 0
        normal = -normal
    end
    dir_t, TIR = refraction3d(direction(ray), normal, refractive_index(ray), n_transmitted)
    if TIR
        return interact3d(Reflective(), system, coating, beam, ray)
    end
    dir_r = reflection3d(direction(ray), normal)

    pos = position(ray) + length(ray) * direction(ray)

    J_t = get_jones_matrix(coating.model, angle3d(direction(ray), -normal),
        wavelength(ray), refractive_index(ray), n_transmitted, false)
    J_r = get_jones_matrix(coating.model, angle3d(direction(ray), -normal),
        wavelength(ray), refractive_index(ray), n_transmitted, true)

    E0_t = _calculate_global_E0(coating, ray, dir_t, J_t)
    E0_r = _calculate_global_E0(coating, ray, dir_r, J_r)

    beam_t = Beam(PolarizedRay{T}(
        pos, dir_t, nothing, wavelength(ray), n_transmitted, E0_t))
    beam_r = Beam(PolarizedRay{T}(pos, dir_r, nothing, wavelength(ray), n_reflected, E0_r))

    children!(beam, (beam_t, beam_r))
    return nothing
end

# GaussianBeamlet and AstigmaticGaussianBeamlet splitting methods
function interact3d(system::AbstractSystem, coating::Coating{T},
        gauss::GaussianBeamlet, ray_id::Int) where {T}
    chief_ray = rays(gauss.chief)[ray_id]
    return interact3d(
        coating_behavior(coating.model, chief_ray), system, coating, gauss, ray_id)
end

function interact3d(::CoatingBehavior, system::AbstractSystem, coating::Coating{T},
        gauss::GaussianBeamlet, ray_id::Int) where {T}
    # For transmissive or reflective coatings, they act as single-beam interactions (no splitting),
    # so we fall back to the generic GaussianBeamlet interaction.
    i_c = interact3d(system, coating, gauss.chief, rays(gauss.chief)[ray_id])
    i_w = interact3d(system, coating, gauss.waist, rays(gauss.waist)[ray_id])
    i_d = interact3d(system, coating, gauss.divergence, rays(gauss.divergence)[ray_id])
    if any(isnothing, (i_c, i_w, i_d))
        return nothing
    end
    return GaussianBeamletInteraction{T}(i_c, i_w, i_d)
end

function interact3d(::Splitting, system::AbstractSystem, coating::Coating{T},
        gauss::GaussianBeamlet, ray_id::Int) where {T}
    # Ensure all component rays have valid intersections to prevent clipping crashes
    if isnothing(intersection(gauss.chief.rays[ray_id])) ||
       isnothing(intersection(gauss.waist.rays[ray_id])) ||
       isnothing(intersection(gauss.divergence.rays[ray_id]))
        return nothing
    end

    c_ray = gauss.chief.rays[ray_id]
    w_ray = gauss.waist.rays[ray_id]
    d_ray = gauss.divergence.rays[ray_id]

    int = intersection(c_ray)
    n_transmitted, n_reflected = _coating_media(system, c_ray, int)

    normal = normal3d(int)
    from_front = dot(direction(c_ray), normal) < 0
    if !from_front
        normal = -normal
    end

    c_dir_t, TIR = refraction3d(
        direction(c_ray), normal, refractive_index(c_ray), n_transmitted)
    if TIR
        i_c = interact3d(Reflective(), system, coating, gauss.chief, rays(gauss.chief)[ray_id])
        i_w = interact3d(Reflective(), system, coating, gauss.waist, rays(gauss.waist)[ray_id])
        i_d = interact3d(Reflective(), system, coating, gauss.divergence, rays(gauss.divergence)[ray_id])
        if any(isnothing, (i_c, i_w, i_d))
            return nothing
        end
        return GaussianBeamletInteraction{T}(i_c, i_w, i_d)
    end
    w_dir_t, _ = refraction3d(
        direction(w_ray), normal, refractive_index(w_ray), n_transmitted)
    d_dir_t, _ = refraction3d(
        direction(d_ray), normal, refractive_index(d_ray), n_transmitted)

    c_dir_r = reflection3d(direction(c_ray), normal)
    w_dir_r = reflection3d(direction(w_ray), normal)
    d_dir_r = reflection3d(direction(d_ray), normal)

    pos = position(c_ray) + length(c_ray) * direction(c_ray)

    J_t = get_jones_matrix(coating.model, angle3d(direction(c_ray), -normal),
        wavelength(c_ray), refractive_index(c_ray), n_transmitted, false; from_front = from_front)
    J_r = get_jones_matrix(coating.model, angle3d(direction(c_ray), -normal),
        wavelength(c_ray), refractive_index(c_ray), n_transmitted, true; from_front = from_front)

    # Spawn transmitted
    λ = wavelength(c_ray)
    chief_t = Beam(Ray{T}(pos, c_dir_t, nothing, λ, n_transmitted))
    waist_t = Beam(Ray{T}(
        position(w_ray) + length(w_ray) * direction(w_ray), w_dir_t, nothing, λ, n_transmitted))
    divergent_t = Beam(Ray{T}(
        position(d_ray) + length(d_ray) * direction(d_ray), d_dir_t, nothing, λ, n_transmitted))

    w0_t = gauss_parameters(gauss, length(gauss))[4]
    t_val = J_t[1, 1]
    E0_t = t_val * electric_field(gauss) * (beam_waist(gauss) / w0_t)
    t = GaussianBeamlet(chief_t, waist_t, divergent_t, wavelength(gauss), w0_t, E0_t)

    # Spawn reflected
    chief_r = Beam(Ray{T}(pos, c_dir_r, nothing, λ, n_reflected))
    waist_r = Beam(Ray{T}(
        position(w_ray) + length(w_ray) * direction(w_ray), w_dir_r, nothing, λ, n_reflected))
    divergent_r = Beam(Ray{T}(
        position(d_ray) + length(d_ray) * direction(d_ray), d_dir_r, nothing, λ, n_reflected))
    w0_r = w0_t
    r_val = -J_r[1, 1]
    E0_r = r_val * electric_field(gauss) * (beam_waist(gauss) / w0_r)
    r = GaussianBeamlet(chief_r, waist_r, divergent_r, wavelength(gauss), w0_r, E0_r)

    children!(gauss, (t, r))
    return nothing
end

function interact3d(system::AbstractSystem, coating::Coating{T},
        agb::AstigmaticGaussianBeamlet, ray_id::Int) where {T}
    chief_ray = rays(agb.c)[ray_id]
    return interact3d(
        coating_behavior(coating.model, chief_ray), system, coating, agb, ray_id)
end

function interact3d(::CoatingBehavior, system::AbstractSystem, coating::Coating{T},
        agb::AstigmaticGaussianBeamlet, ray_id::Int) where {T}
    # Non-splitting behaviors fall back to individual ray tracing
    i_c = interact3d(system, coating, agb.c, rays(agb.c)[ray_id])
    i_wxp = interact3d(system, coating, agb.wxp, rays(agb.wxp)[ray_id])
    i_wxm = interact3d(system, coating, agb.wxm, rays(agb.wxm)[ray_id])
    i_wyp = interact3d(system, coating, agb.wyp, rays(agb.wyp)[ray_id])
    i_wym = interact3d(system, coating, agb.wym, rays(agb.wym)[ray_id])
    i_dxp = interact3d(system, coating, agb.dxp, rays(agb.dxp)[ray_id])
    i_dxm = interact3d(system, coating, agb.dxm, rays(agb.dxm)[ray_id])
    i_dyp = interact3d(system, coating, agb.dyp, rays(agb.dyp)[ray_id])
    i_dym = interact3d(system, coating, agb.dym, rays(agb.dym)[ray_id])
    if any(isnothing, (i_c, i_wxp, i_wxm, i_wyp, i_wym, i_dxp, i_dxm, i_dyp, i_dym))
        return nothing
    end
    return AstigmaticGaussianBeamletInteraction{T}(
        i_c, i_wxp, i_wxm, i_wyp, i_wym, i_dxp, i_dxm, i_dyp, i_dym)
end

function interact3d(::Splitting, system::AbstractSystem, coating::Coating{T},
        agb::AstigmaticGaussianBeamlet, ray_id::Int) where {T}
    # Ensure all component rays have valid intersections to prevent clipping crashes
    if isnothing(intersection(rays(agb.c)[ray_id])) ||
       isnothing(intersection(rays(agb.wxp)[ray_id])) ||
       isnothing(intersection(rays(agb.wxm)[ray_id])) ||
       isnothing(intersection(rays(agb.wyp)[ray_id])) ||
       isnothing(intersection(rays(agb.wym)[ray_id])) ||
       isnothing(intersection(rays(agb.dxp)[ray_id])) ||
       isnothing(intersection(rays(agb.dxm)[ray_id])) ||
       isnothing(intersection(rays(agb.dyp)[ray_id])) ||
       isnothing(intersection(rays(agb.dym)[ray_id]))
        return nothing
    end

    c_ray = rays(agb.c)[ray_id]
    int = intersection(c_ray)
    n_transmitted, n_reflected = _coating_media(system, c_ray, int)

    normal = normal3d(int)
    from_front = dot(direction(c_ray), normal) < 0
    if !from_front
        normal = -normal
    end

    _, TIR = refraction3d(direction(c_ray), normal, refractive_index(c_ray), n_transmitted)
    if TIR
        i_c = interact3d(Reflective(), system, coating, agb.c, rays(agb.c)[ray_id])
        i_wxp = interact3d(Reflective(), system, coating, agb.wxp, rays(agb.wxp)[ray_id])
        i_wxm = interact3d(Reflective(), system, coating, agb.wxm, rays(agb.wxm)[ray_id])
        i_wyp = interact3d(Reflective(), system, coating, agb.wyp, rays(agb.wyp)[ray_id])
        i_wym = interact3d(Reflective(), system, coating, agb.wym, rays(agb.wym)[ray_id])
        i_dxp = interact3d(Reflective(), system, coating, agb.dxp, rays(agb.dxp)[ray_id])
        i_dxm = interact3d(Reflective(), system, coating, agb.dxm, rays(agb.dxm)[ray_id])
        i_dyp = interact3d(Reflective(), system, coating, agb.dyp, rays(agb.dyp)[ray_id])
        i_dym = interact3d(Reflective(), system, coating, agb.dym, rays(agb.dym)[ray_id])
        if any(isnothing, (i_c, i_wxp, i_wxm, i_wyp, i_wym, i_dxp, i_dxm, i_dyp, i_dym))
            return nothing
        end
        return AstigmaticGaussianBeamletInteraction{T}(
            i_c, i_wxp, i_wxm, i_wyp, i_wym, i_dxp, i_dxm, i_dyp, i_dym)
    end

    dirs_t = (
        refraction3d(rays(agb.c)[ray_id], n_transmitted)[1],
        refraction3d(rays(agb.wxp)[ray_id], n_transmitted)[1],
        refraction3d(rays(agb.wxm)[ray_id], n_transmitted)[1],
        refraction3d(rays(agb.wyp)[ray_id], n_transmitted)[1],
        refraction3d(rays(agb.wym)[ray_id], n_transmitted)[1],
        refraction3d(rays(agb.dxp)[ray_id], n_transmitted)[1],
        refraction3d(rays(agb.dxm)[ray_id], n_transmitted)[1],
        refraction3d(rays(agb.dyp)[ray_id], n_transmitted)[1],
        refraction3d(rays(agb.dym)[ray_id], n_transmitted)[1]
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
    # Transmitted
    λ = wavelength(c_ray)
    J_t = get_jones_matrix(coating.model, angle3d(direction(c_ray), -normal),
        λ, refractive_index(c_ray), n_transmitted, false; from_front = from_front)
    J_r = get_jones_matrix(coating.model, angle3d(direction(c_ray), -normal),
        λ, refractive_index(c_ray), n_transmitted, true; from_front = from_front)

    E0_t = _calculate_global_E0(coating, c_ray, dirs_t[1], J_t)
    E0_r = _calculate_global_E0(coating, c_ray, dirs_r[1], J_r)

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

    # Reflected
    chief_r = Beam(PolarizedRay{T}(
        position(c_ray) + length(c_ray) * direction(c_ray), dirs_r[1], nothing, λ, n_reflected, E0_r))
    wxp_r = Beam(Ray{T}(
        position(rays(agb.wxp)[ray_id]) +
        length(rays(agb.wxp)[ray_id]) * direction(rays(agb.wxp)[ray_id]),
        dirs_r[2], nothing, λ, n_reflected))
    wxm_r = Beam(Ray{T}(
        position(rays(agb.wxm)[ray_id]) +
        length(rays(agb.wxm)[ray_id]) * direction(rays(agb.wxm)[ray_id]),
        dirs_r[3], nothing, λ, n_reflected))
    wyp_r = Beam(Ray{T}(
        position(rays(agb.wyp)[ray_id]) +
        length(rays(agb.wyp)[ray_id]) * direction(rays(agb.wyp)[ray_id]),
        dirs_r[4], nothing, λ, n_reflected))
    wym_r = Beam(Ray{T}(
        position(rays(agb.wym)[ray_id]) +
        length(rays(agb.wym)[ray_id]) * direction(rays(agb.wym)[ray_id]),
        dirs_r[5], nothing, λ, n_reflected))
    dxp_r = Beam(Ray{T}(
        position(rays(agb.dxp)[ray_id]) +
        length(rays(agb.dxp)[ray_id]) * direction(rays(agb.dxp)[ray_id]),
        dirs_r[6], nothing, λ, n_reflected))
    dxm_r = Beam(Ray{T}(
        position(rays(agb.dxm)[ray_id]) +
        length(rays(agb.dxm)[ray_id]) * direction(rays(agb.dxm)[ray_id]),
        dirs_r[7], nothing, λ, n_reflected))
    dyp_r = Beam(Ray{T}(
        position(rays(agb.dyp)[ray_id]) +
        length(rays(agb.dyp)[ray_id]) * direction(rays(agb.dyp)[ray_id]),
        dirs_r[8], nothing, λ, n_reflected))
    dym_r = Beam(Ray{T}(
        position(rays(agb.dym)[ray_id]) +
        length(rays(agb.dym)[ray_id]) * direction(rays(agb.dym)[ray_id]),
        dirs_r[9], nothing, λ, n_reflected))

    r = AstigmaticGaussianBeamlet(
        chief_r, wxp_r, wxm_r, wyp_r, wym_r, dxp_r, dxm_r, dyp_r, dym_r)

    children!(agb, (t, r))
    return nothing
end
