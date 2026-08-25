"""
    interact3d(model::AbstractSurfaceModel, trans::Transition, int::Intersection, ray)

Universal boundary interaction function evaluating boundary physics on the given surface model.
"""
function interact3d(
    ::FresnelSurface,
    trans::Transition,
    int::Intersection{T},
    ray::Ray{T}
) where {T <: Real}
    λ = wavelength(ray)
    n1 = refractive_index(trans.medium_in, λ)
    n2 = refractive_index(trans.medium_out, λ)
    
    normal = dot(normal3d(int), direction(ray)) < 0 ? normal3d(int) : -normal3d(int)
    new_dir, tir = refraction3d(direction(ray), normal, n1, n2)
    
    new_pos = position(int)
    new_n = tir ? n1 : n2
    
    return Ray{T}(new_pos, new_dir, λ, new_n)
end

function interact3d(
    ::FresnelSurface,
    trans::Transition,
    int::Intersection{T},
    ray::PolarizedRay{T}
) where {T <: Real}
    λ = wavelength(ray)
    n1 = refractive_index(trans.medium_in, λ)
    n2 = refractive_index(trans.medium_out, λ)
    
    normal = dot(normal3d(int), direction(ray)) < 0 ? normal3d(int) : -normal3d(int)
    θi = angle3d(direction(ray), -normal)
    rs, rp, ts, tp = fresnel_coefficients(θi, n2 / n1)
    
    if is_internally_reflected(rp, rs)
        new_dir = reflection3d(direction(ray), normal)
        J = SPBasis(-rs, 0, 0, rp)
        new_n = n1
    else
        new_dir, _ = refraction3d(direction(ray), normal, n1, n2)
        J = SPBasis(ts, 0, 0, tp)
        new_n = n2
    end
    
    P = _calculate_global_E0(direction(ray), new_dir, normal, J)
    new_E0 = P * polarization(ray)
    new_pos = position(int)
    
    return PolarizedRay{T}(new_pos, new_dir, λ, new_n, new_E0)
end

function interact3d(
    ::IdealMirrorSurface,
    trans::Transition,
    int::Intersection{T},
    ray::Ray{T}
) where {T <: Real}
    normal = normal3d(int)
    new_dir = reflection3d(direction(ray), normal)
    new_pos = position(int)
    return Ray{T}(new_pos, new_dir, wavelength(ray), refractive_index(ray))
end

function interact3d(
    ::IdealMirrorSurface,
    trans::Transition,
    int::Intersection{T},
    ray::PolarizedRay{T}
) where {T <: Real}
    normal = normal3d(int)
    new_dir = reflection3d(direction(ray), normal)
    new_pos = position(int)
    J = SPBasis(-1, 0, 0, 1)
    P = _calculate_global_E0(direction(ray), new_dir, normal, J)
    new_E0 = P * polarization(ray)
    return PolarizedRay{T}(new_pos, new_dir, wavelength(ray), refractive_index(ray), new_E0)
end

function interact3d(
    ::AbsorbingSurface,
    trans::Transition,
    int::Intersection{T},
    ray::AbstractRay{T}
) where {T <: Real}
    return nothing
end

"""
    record!(data, d::DetectorSurface, int::Intersection, ray)

Extensible data recording interface for detector surfaces.
"""
function record!(data::Any, ::DetectorSurface, int::Intersection{T}, ray::AbstractRay{T}) where {T <: Real}
    push!(data, (position(int), ray))
    return nothing
end

function interact3d(
    d::DetectorSurface,
    trans::Transition,
    int::Intersection{T},
    ray::AbstractRay{T}
) where {T <: Real}
    if d.detector_data !== nothing
        record!(d.detector_data, d, int, ray)
    end
    return nothing
end

@inline _eval_coating_model(m::AbstractSurfaceModel, trans, int, ray) = interact3d(m, trans, int, ray)
@inline _eval_coating_model(f::Function, trans, int, ray) = f(trans, int, ray)
@inline _eval_coating_model(::Any, trans, int, ray) = interact3d(FresnelSurface(), trans, int, ray)

function interact3d(
    c::CoatedSurface,
    trans::Transition,
    int::Intersection{T},
    ray::AbstractRay{T}
) where {T <: Real}
    return _eval_coating_model(c.coating_model, trans, int, ray)
end

function interact3d(
    g::GratingSurface,
    trans::Transition,
    int::Intersection{T},
    ray::AbstractRay{T}
) where {T <: Real}
    normal = normal3d(int)
    dir = direction(ray)
    λ = wavelength(ray)
    n_med = refractive_index(trans.medium_in, λ)

    # Determine groove vector in the tangent plane
    g_vec = if g.groove_vector !== nothing
        normalize(Point3{T}(g.groove_vector))
    else
        # Default groove axis perpendicular to incidence plane or along local tangent
        t_ref = isparallel3d(dir, normal) ? normal3d(dir) : cross(dir, normal)
        normalize(Point3{T}(t_ref))
    end
    p_vec = normalize(cross(g_vec, normal))

    # Grating equation projection components
    d_g = dot(dir, g_vec)
    d_p_in = dot(dir, p_vec)
    m = g.order
    d_p_out = d_p_in + (m * λ * g.groove_density) / n_med

    # Unit normalization constraint
    d_n_sq = 1 - d_g^2 - d_p_out^2
    if d_n_sq < 0
        return nothing # Order is evanescent / non-propagating
    end

    # Reflective diffracted order directed back into incident half-space
    d_n_out = (dot(dir, normal) < 0) ? sqrt(d_n_sq) : -sqrt(d_n_sq)
    new_dir = normalize(d_g * g_vec + d_p_out * p_vec + d_n_out * normal)
    new_pos = position(int)

    return Ray{T}(new_pos, new_dir, λ, refractive_index(ray))
end

"""
    interact3d(system, mi::MultiIntersection, beam, ray)

Handles coincident interface interaction across touching optical components.
"""
function interact3d(
    system::AbstractSystem,
    mi::MultiIntersection,
    beam::AbstractBeam{T},
    ray::R
) where {T <: Real, R <: AbstractRay{T}}
    normal = normal3d(mi)
    int = mi.hit
    
    ambient = ambient_medium(system)
    trans = resolve_transition(mi, ambient, ray)
    
    bs_obj = _first_beamsplitter(coating(mi), entering(mi), exiting(mi))
    if bs_obj !== nothing
        return interact3d(system, bs_obj, trans, int, beam, ray)
    end
    
    target_obj = something(entering(mi), exiting(mi))
    model = surface_model(target_obj)
    
    out_ray = interact3d(model, trans, int, ray)
    if out_ray === nothing
        return nothing
    end
    
    hint = entering(mi) !== nothing ? Hint(entering(mi)) : nothing
    return BeamInteraction(hint, out_ray)
end

"""
    interact3d(system, c::AbstractCoating, beam, ray)

Handles standalone optical boundary coating interactions.
"""
function interact3d(
    system::AbstractSystem,
    c::AbstractCoating,
    beam::AbstractBeam{T},
    ray::AbstractRay{T}
) where {T <: Real}
    int = isempty(beam.segments) ? intersect3d(shape(c), ray) : intersection(last(beam.segments))
    isnothing(int) && return nothing
    ambient = ambient_medium(system)
    trans = resolve_transition(ambient, ambient, ray, normal3d(int))
    out_ray = interact3d(surface_model(c), trans, int, ray)
    if out_ray === nothing
        return nothing
    end
    return BeamInteraction(nothing, out_ray)
end



