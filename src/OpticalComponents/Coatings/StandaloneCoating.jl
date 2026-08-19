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
        # A threshold of 0.5 corresponds to cos(60°), restricting the interaction
        # to rays hitting the surface within a ±60° angle of the target normal direction.
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

"""
    _resolve_coincident_hint(system, ray, int, is_entering)

Resolves the target object hint (`Hint`) at a standalone coating boundary by comparing the ray's
current refractive index against the coincident adjacent objects (`coincident_object` and `coincident_object_2`).
If `is_entering` is `true`, returns the object entered upon transmission; otherwise returns the object remaining in upon reflection.
"""
@inline function _resolve_coincident_hint(
        system::AbstractSystem, ray::AbstractRay, int::Intersection, is_entering::Bool)
    n_exit = int.coincident_object
    n_enter = int.coincident_object_2
    (isnothing(n_exit) && isnothing(n_enter)) && return nothing

    λ = wavelength(ray)
    n_sys = refractive_index(system, λ)
    n_exit_val = isnothing(n_exit) ? n_sys :
                 (is_refractive(n_exit) ? refractive_index(n_exit, λ) : n_sys)

    tol = get_index_matching_tolerance()
    is_matching = abs(refractive_index(ray) - n_exit_val) < tol
    target_obj = is_entering ? (is_matching ? n_enter : n_exit) : (is_matching ? n_exit : n_enter)
    return isnothing(target_obj) ? nothing : Hint(target_obj)
end

"""
    _resolve_coating_normal(ray, int)

Computes the interface normal vector aligned against the incoming ray direction (`-normal`), along with
the boolean `from_front` flag indicating whether the ray hit the interface from the front side.
"""
@inline function _resolve_coating_normal(ray::AbstractRay, int::Intersection)
    normal = normal3d(int)
    from_front = dot(direction(ray), normal) < 0
    return from_front ? normal : -normal, from_front
end

@inline function _resolve_standalone_coating_context(
        system::AbstractSystem, coating::Coating, ray::AbstractRay)
    int = intersection(ray)
    λ = wavelength(ray)
    n_transmitted, _ = _coating_media(system, ray, int)
    normal, from_front = _resolve_coating_normal(ray, int)
    hint_trans = _resolve_coincident_hint(system, ray, int, true)
    n_incident = refractive_index(ray)
    return (normal, n_incident, n_transmitted, hint_trans, λ, from_front)
end

# Route interact3d through coating_behavior trait into unified boundary physics engine
function interact3d(system::AbstractSystem, coating::Coating{T},
        beam::AbstractBeam, ray::AbstractRay) where {T}
    normal, n_incident, n_transmitted, hint, λ, from_front =
        _resolve_standalone_coating_context(system, coating, ray)
    behavior = coating_behavior(coating.model, ray)
    return interact_refractive_boundary(
        behavior, system, coating, coating.model, beam,
        ray, n_incident, n_transmitted, hint, normal, λ, from_front)
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
    int = intersection(c_ray)
    n_transmitted, _ = _coating_media(system, c_ray, int)
    normal, from_front = _resolve_coating_normal(c_ray, int)

    return _propagate_splitting_gaussian_beamlet(
        system,
        coating,
        coating.model,
        gauss,
        ray_id,
        refractive_index(c_ray),
        n_transmitted,
        normal,
        from_front,
        () -> begin
            ints = _interact3d_reflective_component_beams(system, coating, coating.model, gauss, ray_id)
            isnothing(ints) && return nothing
            return GaussianBeamletInteraction{T}(ints...)
        end
    )
end

function interact3d(system::AbstractSystem, coating::Coating{T},
        agb::AstigmaticGaussianBeamlet, ray_id::Int) where {T}
    chief_ray = rays(agb.c)[ray_id]
    return interact3d(
        coating_behavior(coating.model, chief_ray), system, coating, agb, ray_id)
end

function interact3d(::CoatingBehavior, system::AbstractSystem, coating::Coating{T},
        agb::AstigmaticGaussianBeamlet, ray_id::Int) where {T}
    ints = _interact3d_component_beams(system, coating, agb, ray_id)
    isnothing(ints) && return nothing
    return AstigmaticGaussianBeamletInteraction{T}(ints...)
end

function interact3d(::Splitting, system::AbstractSystem, coating::Coating{T},
        agb::AstigmaticGaussianBeamlet, ray_id::Int) where {T}
    # Ensure all component rays have valid intersections to prevent clipping crashes
    beams = _component_beams(agb)
    if any(b -> isnothing(intersection(rays(b)[ray_id])), beams)
        return nothing
    end

    c_ray = rays(agb.c)[ray_id]
    int = intersection(c_ray)
    n_transmitted, _ = _coating_media(system, c_ray, int)
    normal, from_front = _resolve_coating_normal(c_ray, int)

    return _propagate_splitting_astigmatic_beamlet(
        system,
        coating,
        coating.model,
        agb,
        ray_id,
        refractive_index(c_ray),
        n_transmitted,
        normal,
        from_front,
        () -> begin
            ints = _interact3d_reflective_component_beams(system, coating, coating.model, agb, ray_id)
            isnothing(ints) && return nothing
            return AstigmaticGaussianBeamletInteraction{T}(ints...)
        end
    )
end

# Absorptive behavior
function interact3d(::Absorptive, ::AbstractSystem, ::Coating{T},
        ::AbstractBeam, ::AbstractRay) where {T}
    nothing
end
function interact3d(
        ::Absorptive, ::AbstractSystem, ::Coating{T}, ::GaussianBeamlet, ::Int) where {T}
    nothing
end
function interact3d(::Absorptive, ::AbstractSystem, ::Coating{T},
        ::AstigmaticGaussianBeamlet, ::Int) where {T}
    nothing
end
