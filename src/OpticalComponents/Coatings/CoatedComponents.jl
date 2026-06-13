# Coated unified components (CoatedLens, CoatedMirror) and fluent API pipeline helpers

# Matching Helpers for Coated Components
function eval_filter(face::Symbol, shape::AbstractShape, local_p, local_n)
    if face === :either || face === :all
        return true
    end
    return face_id(shape, local_n) === face
end

function eval_filter(filter_vec::AbstractVector, shape::AbstractShape, local_p, local_n)
    return dot(local_n, filter_vec) > 0.0
end

function eval_filter(pred::Function, shape::AbstractShape, local_p, local_n)
    return pred(local_p, local_n)
end

function eval_filter(::Nothing, shape::AbstractShape, local_p, local_n)
    return true
end

function get_matching_coating(coatings, shape::AbstractShape, local_p, local_n)
    for (filt, model) in coatings
        if eval_filter(filt, shape, local_p, local_n)
            return model
        end
    end
    return Uncoated()
end

# Helper to determine if a shape is a 2D flat shape (thickness = 0)
is_flat_shape(::AbstractShape) = false
function is_flat_shape(mesh::Mesh)
    n = size(mesh.vertices, 1)
    n < 3 && return true
    v1 = @view mesh.vertices[1, :]
    v2 = @view mesh.vertices[2, :]
    v3 = @view mesh.vertices[3, :]
    n_vec = cross(v2 - v1, v3 - v1)
    n_norm = norm(n_vec)
    if n_norm < 1e-12
        return true
    end
    n_vec = n_vec / n_norm
    for i in 4:n
        vi = @view mesh.vertices[i, :]
        if abs(dot(vi - v1, n_vec)) > 1e-6
            return false
        end
    end
    return true
end

# Composite Groups / Unified wrappers
"""
    CoatedRefractive{T, L, C} <: AbstractRefractiveOptic{T, RefractiveIndex}

Wrapper type that attaches a tuple of coating models to an `AbstractRefractiveOptic`.
"""
struct CoatedRefractive{T, L <: AbstractRefractiveOptic{T}, C} <:
       AbstractRefractiveOptic{T, RefractiveIndex}
    optic::L
    coatings::C
end

"""
    CoatedMirror{T, M, C} <: AbstractReflectiveOptic{T}

Wrapper type that attaches a tuple of coating models to an `AbstractReflectiveOptic`.
"""
struct CoatedMirror{T, M <: AbstractReflectiveOptic{T}, C} <: AbstractReflectiveOptic{T}
    optic::M
    coatings::C
end

const CoatedComponent = Union{CoatedRefractive, CoatedMirror}
const CoatedLens = CoatedRefractive

# Dispatch helpers for extracting coating models
_coating_model(c::AbstractCoating) = c.model
_coating_model(c) = c

# Shared coatings tuple constructor helper
function _build_coatings_tuple(parent_shape::AbstractShape, front, back)
    f_side = is_flat_shape(parent_shape) ? :either : :front
    b_side = is_flat_shape(parent_shape) ? :either : :back

    p1 = !isnothing(front) ? (f_side => _coating_model(front)) : nothing
    p2 = !isnothing(back) ? (b_side => _coating_model(back)) : nothing

    if isnothing(p1) && isnothing(p2)
        return ()
    elseif isnothing(p1)
        return (p2,)
    elseif isnothing(p2)
        return (p1,)
    else
        return (p1, p2)
    end
end



function CoatedRefractive(
        optic::AbstractRefractiveOptic; front = nothing, back = nothing)
    coatings = _build_coatings_tuple(shape(optic), front, back)
    return CoatedRefractive(optic, coatings)
end

function CoatedMirror(
        mirror::AbstractReflectiveOptic; front = nothing, back = nothing)
    coatings = _build_coatings_tuple(shape(mirror), front, back)
    return CoatedMirror(mirror, coatings)
end

function Base.getproperty(cr::CoatedRefractive, sym::Symbol)
    if sym === :lens
        return getfield(cr, :optic)
    else
        return getfield(cr, sym)
    end
end

function Base.getproperty(cm::CoatedMirror, sym::Symbol)
    if sym === :mirror
        return getfield(cm, :optic)
    else
        return getfield(cm, sym)
    end
end

# Unified forwarding methods
shape_trait_of(::CoatedComponent) = SingleShape()
shape(cl::CoatedComponent) = shape(cl.optic)

Base.position(cl::CoatedComponent) = position(cl.optic)
position!(cl::CoatedComponent, pos) = position!(cl.optic, pos)
orientation(cl::CoatedComponent) = orientation(cl.optic)
orientation!(cl::CoatedComponent, dir) = orientation!(cl.optic, dir)
translate3d!(cl::CoatedComponent, offset) = translate3d!(cl.optic, offset)
translate_to3d!(cl::CoatedComponent, target) = translate_to3d!(cl.optic, target)
rotate3d!(cl::CoatedComponent, axis, θ) = rotate3d!(cl.optic, axis, θ)

# CoatedRefractive specific forwards
refractive_index(cl::CoatedRefractive) = refractive_index(cl.optic)
refractive_index(cl::CoatedRefractive, λ::Real) = refractive_index(cl.optic, λ)
is_refractive(::CoatedRefractive) = true
thickness(cl::CoatedRefractive) = thickness(cl.optic)

# Get matching coating for a hit
function get_coating_model_at_hit(coated::CoatedComponent, ray::AbstractRay)
    int = intersection(ray)
    isnothing(int) && return Uncoated()
    p_hit = position(ray) + length(ray) * direction(ray)
    local_p = world_to_local(shape(coated), p_hit)
    n_world = normal3d(int)
    if object(int) !== coated && object(int) !== coated.optic
        n_world = -n_world
    end
    local_n = transposed_orientation(shape(coated)) * n_world
    return get_matching_coating(coated.coatings, shape(coated), local_p, local_n)
end

function get_coating_model_at_hit(coated::CoatedComponent, gauss::GaussianBeamlet, ray_id::Int)
    return get_coating_model_at_hit(coated, gauss.chief.rays[ray_id])
end

function get_coating_model_at_hit(coated::CoatedComponent, agb::AstigmaticGaussianBeamlet, ray_id::Int)
    return get_coating_model_at_hit(coated, agb.c.rays[ray_id])
end

# Coincident boundary resolution dispatch helpers
resolve_coated_boundary(system::AbstractSystem, obj::AbstractObject, ray::AbstractRay) = resolve_coated_boundary_dispatch(obj, system, ray)

function resolve_coated_boundary_dispatch(obj::CoatedComponent, system::AbstractSystem, ray::AbstractRay)
    coating = get_coating_model_at_hit(obj, ray)
    if !(coating isa Uncoated)
        return obj, coating
    end
    return resolve_coincident_coatings(intersection(ray), system, ray)
end

resolve_coated_boundary_dispatch(::AbstractObject, system::AbstractSystem, ray::AbstractRay) = resolve_coincident_coatings(intersection(ray), system, ray)

resolve_coincident_coatings(::Nothing, system::AbstractSystem, ray::AbstractRay) = (nothing, Uncoated())

function resolve_coincident_coatings(int::Intersection, system::AbstractSystem, ray::AbstractRay)
    res1 = check_coincident_coating(int.coincident_object, ray)
    res1 !== nothing && return res1
    res2 = check_coincident_coating(int.coincident_object_2, ray)
    res2 !== nothing && return res2
    return nothing, Uncoated()
end

check_coincident_coating(::Nothing, ray::AbstractRay) = nothing
check_coincident_coating(::AbstractObject, ray::AbstractRay) = nothing

function check_coincident_coating(obj::CoatedComponent, ray::AbstractRay)
    coating = get_coating_model_at_hit(obj, ray)
    if !(coating isa Uncoated)
        return (obj, coating)
    end
    return nothing
end

# interact3d definitions for CoatedComponent single-ray propagation
function interact3d(system::AbstractSystem,
        cl::CoatedRefractive,
        beam::Beam{T, R},
        ray::R) where {T <: Real, R <: AbstractRay{T}}
    coating_model = get_coating_model_at_hit(cl, ray)
    return interact_refractive_boundary(system, cl.optic, coating_model, beam, ray)
end

function interact3d(system::AbstractSystem,
        cm::CoatedMirror,
        beam::Beam{T, R},
        ray::R) where {T <: Real, R <: AbstractRay{T}}
    coating_model = get_coating_model_at_hit(cm, ray)
    return interact_reflective_boundary(system, cm.optic, coating_model, beam, ray)
end

# interact3d definitions for CoatedComponent beamlet propagation using trait dispatch
function interact3d(
        system::AbstractSystem, cl::CoatedComponent, gauss::GaussianBeamlet, ray_id::Int)
    coating_model = get_coating_model_at_hit(cl, gauss, ray_id)
    chief_ray = rays(gauss.chief)[ray_id]
    return interact3d_behavior(coating_behavior(coating_model, chief_ray), system, cl, coating_model, gauss, ray_id)
end

function interact3d_behavior(::Splitting, system::AbstractSystem, cl::CoatedComponent, coating_model, gauss::GaussianBeamlet, ray_id::Int)
    return interact_splitting_boundary(system, cl, coating_model, gauss, ray_id)
end

function interact3d_behavior(::CoatingBehavior, system::AbstractSystem, cl::CoatedComponent, coating_model, gauss::GaussianBeamlet, ray_id::Int)
    i_c = interact3d(system, cl, gauss.chief, rays(gauss.chief)[ray_id])
    i_w = interact3d(system, cl, gauss.waist, rays(gauss.waist)[ray_id])
    i_d = interact3d(system, cl, gauss.divergence, rays(gauss.divergence)[ray_id])
    if any(isnothing, (i_c, i_w, i_d))
        return nothing
    end
    T = eltype(position(gauss.chief.rays[ray_id]))
    return GaussianBeamletInteraction{T}(i_c, i_w, i_d)
end

function interact3d(
        system::AbstractSystem, cl::CoatedComponent, agb::AstigmaticGaussianBeamlet, ray_id::Int)
    coating_model = get_coating_model_at_hit(cl, agb, ray_id)
    chief_ray = rays(agb.c)[ray_id]
    return interact3d_behavior(coating_behavior(coating_model, chief_ray), system, cl, coating_model, agb, ray_id)
end

function interact3d_behavior(::Splitting, system::AbstractSystem, cl::CoatedComponent, coating_model, agb::AstigmaticGaussianBeamlet, ray_id::Int)
    return interact_splitting_boundary(system, cl, coating_model, agb, ray_id)
end

function interact3d_behavior(::CoatingBehavior, system::AbstractSystem, cl::CoatedComponent, coating_model, agb::AstigmaticGaussianBeamlet, ray_id::Int)
    i_c = interact3d(system, cl, agb.c, rays(agb.c)[ray_id])
    i_wxp = interact3d(system, cl, agb.wxp, rays(agb.wxp)[ray_id])
    i_wxm = interact3d(system, cl, agb.wxm, rays(agb.wxm)[ray_id])
    i_wyp = interact3d(system, cl, agb.wyp, rays(agb.wyp)[ray_id])
    i_wym = interact3d(system, cl, agb.wym, rays(agb.wym)[ray_id])
    i_dxp = interact3d(system, cl, agb.dxp, rays(agb.dxp)[ray_id])
    i_dxm = interact3d(system, cl, agb.dxm, rays(agb.dxm)[ray_id])
    i_dyp = interact3d(system, cl, agb.dyp, rays(agb.dyp)[ray_id])
    i_dym = interact3d(system, cl, agb.dym, rays(agb.dym)[ray_id])
    if any(isnothing, (i_c, i_wxp, i_wxm, i_wyp, i_wym, i_dxp, i_dxm, i_dyp, i_dym))
        return nothing
    end
    T = eltype(position(agb.c.rays[ray_id]))
    return AstigmaticGaussianBeamletInteraction{T}(
        i_c, i_wxp, i_wxm, i_wyp, i_wym, i_dxp, i_dxm, i_dyp, i_dym)
end

function interact_splitting_boundary(
        system::AbstractSystem,
        coated::CoatedComponent,
        coating_model,
        gauss::GaussianBeamlet,
        ray_id::Int
)
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
    λ = wavelength(c_ray)
    substrate_obj = coated.optic
    n_substrate = refractive_index(substrate_obj, λ)

    is_coated_obj = (object(int) === coated || object(int) === coated.optic)
    normal_coated = normal3d(int)
    if !is_coated_obj
        normal_coated = -normal_coated
    end

    entering_coated = dot(direction(c_ray), normal_coated) < 0

    if entering_coated
        n_incident = refractive_index(c_ray)
        n_transmitted = n_substrate
        normal = normal_coated
    else
        n_incident = n_substrate
        other_obj = is_coated_obj ? int.coincident_object : object(int)
        if !isnothing(other_obj) && is_refractive(other_obj)
            n_transmitted = refractive_index(other_obj, λ)
        else
            n_transmitted = refractive_index(system, λ)
        end
        normal = -normal_coated
    end

    from_front = entering_coated

    c_dir_t, _ = refraction3d(direction(c_ray), normal, n_incident, n_transmitted)
    w_dir_t, _ = refraction3d(direction(w_ray), normal, n_incident, n_transmitted)
    d_dir_t, _ = refraction3d(direction(d_ray), normal, n_incident, n_transmitted)

    c_dir_r = reflection3d(direction(c_ray), normal)
    w_dir_r = reflection3d(direction(w_ray), normal)
    d_dir_r = reflection3d(direction(d_ray), normal)

    pos = position(c_ray) + length(c_ray) * direction(c_ray)

    J_t = get_jones_matrix(coating_model, angle3d(direction(c_ray), -normal),
        λ, n_incident, n_transmitted, false; from_front=from_front)
    J_r = get_jones_matrix(coating_model, angle3d(direction(c_ray), -normal),
        λ, n_incident, n_transmitted, true; from_front=from_front)

    T = eltype(pos)

    # Spawn transmitted
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
    chief_r = Beam(Ray{T}(pos, c_dir_r, nothing, λ, n_incident))
    waist_r = Beam(Ray{T}(
        position(w_ray) + length(w_ray) * direction(w_ray), w_dir_r, nothing, λ, n_incident))
    divergent_r = Beam(Ray{T}(
        position(d_ray) + length(d_ray) * direction(d_ray), d_dir_r, nothing, λ, n_incident))

    w0_r = w0_t
    r_val = -J_r[1, 1]
    E0_r = r_val * electric_field(gauss) * (beam_waist(gauss) / w0_r)
    r = GaussianBeamlet(chief_r, waist_r, divergent_r, wavelength(gauss), w0_r, E0_r)

    children!(gauss, (t, r))
    return nothing
end

function interact_splitting_boundary(
        system::AbstractSystem,
        coated::CoatedComponent,
        coating_model,
        agb::AstigmaticGaussianBeamlet,
        ray_id::Int
)
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
    λ = wavelength(c_ray)
    substrate_obj = coated.optic
    n_substrate = refractive_index(substrate_obj, λ)

    is_coated_obj = (object(int) === coated || object(int) === coated.optic)
    normal_coated = normal3d(int)
    if !is_coated_obj
        normal_coated = -normal_coated
    end

    entering_coated = dot(direction(c_ray), normal_coated) < 0

    if entering_coated
        n_incident = refractive_index(c_ray)
        n_transmitted = n_substrate
        normal = normal_coated
    else
        n_incident = n_substrate
        other_obj = is_coated_obj ? int.coincident_object : object(int)
        if !isnothing(other_obj) && is_refractive(other_obj)
            n_transmitted = refractive_index(other_obj, λ)
        else
            n_transmitted = refractive_index(system, λ)
        end
        normal = -normal_coated
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

    from_front = entering_coated

    J_t = get_jones_matrix(coating_model, angle3d(direction(c_ray), -normal),
        λ, n_incident, n_transmitted, false; from_front=from_front)
    J_r = get_jones_matrix(coating_model, angle3d(direction(c_ray), -normal),
        λ, n_incident, n_transmitted, true; from_front=from_front)

    E0_t = _calculate_global_E0(coated, c_ray, dirs_t[1], J_t)
    E0_r = _calculate_global_E0(coated, c_ray, dirs_r[1], J_r)

    T = eltype(position(c_ray))

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

    # Reflected
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

# Fluent API / Pipeline helper
"""
    with_coatings(coatings::Pair...)
    with_coatings(; front=nothing, back=nothing)
    with_coatings(optic, coatings::Pair...)
    with_coatings(optic; front=nothing, back=nothing)

Fluent API helper to attach coatings to an optical component. 
Returns a `CoatedRefractive` or `CoatedMirror` depending on the type of `optic`.

# Examples
```julia
with_coatings(lens, :front => AR, :back => HR)
with_coatings(lens, :side => SimpleARCoating())
with_coatings(mirror, :either => HR)
```
"""
function with_coatings(coatings::Pair...)
    obj -> with_coatings(obj, coatings...)
end

function with_coatings(lens::AbstractRefractiveOptic, coatings::Pair...)
    mapped = map(c -> (c.first => _coating_model(c.second)), coatings)
    return CoatedRefractive(lens, mapped)
end

function with_coatings(mirror::AbstractReflectiveOptic, coatings::Pair...)
    mapped = map(c -> (c.first => _coating_model(c.second)), coatings)
    return CoatedMirror(mirror, mapped)
end

# Keyword variants
function with_coatings(; front = nothing, back = nothing)
    obj -> with_coatings(obj; front = front, back = back)
end
function with_coatings(lens::AbstractRefractiveOptic; front = nothing, back = nothing)
    CoatedRefractive(lens; front = front, back = back)
end
function with_coatings(mirror::AbstractReflectiveOptic; front = nothing, back = nothing)
    CoatedMirror(mirror; front = front, back = back)
end
