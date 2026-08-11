# Coated unified components and fluent API pipeline helpers using pure composition

# Generic getter for coatings property on any AbstractObject
coatings(obj::AbstractObject) = hasproperty(obj, :coatings) ? getproperty(obj, :coatings) : ()

# # Matching Helpers for Coated Components
function eval_filter(face::Symbol, shape::AbstractShape, local_p, local_n)
    if face === :either || face === :all
        return true
    end
    tag = surface_tag(shape, local_p, local_n)
    if tag !== :unknown
        return tag === face
    end
    return face_id(shape, local_n) === face
end

function eval_filter(faces::Union{Tuple{Vararg{Symbol}}, AbstractSet{Symbol}}, shape::AbstractShape, local_p, local_n)
    tag = surface_tag(shape, local_p, local_n)
    if tag !== :unknown
        return tag in faces
    end
    return face_id(shape, local_n) in faces
end

function eval_filter(faces::AbstractVector{Symbol}, shape::AbstractShape, local_p, local_n)
    tag = surface_tag(shape, local_p, local_n)
    if tag !== :unknown
        return tag in faces
    end
    return face_id(shape, local_n) in faces
end

function eval_filter(filter_vec::AbstractVector{<:Real}, shape::AbstractShape, local_p, local_n)
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

    # Find a non-collinear triplet of vertices to define the plane
    v1 = @view mesh.vertices[1, :]

    # Find the first vertex vi that is not coincident with v1
    i = 2
    while i <= n && norm(mesh.vertices[i, :] - v1) < 1e-12
        i += 1
    end
    i > n && return true  # All vertices are coincident (degenerate, flat)
    vi = @view mesh.vertices[i, :]

    # Find the first vertex vj such that v1, vi, vj are non-collinear
    j = i + 1
    n_vec = cross(vi - v1, vi - v1)  # Initial dummy value
    n_norm = 0.0
    found = false
    while j <= n
        vj = @view mesh.vertices[j, :]
        n_vec = cross(vi - v1, vj - v1)
        n_norm = norm(n_vec)
        if n_norm > 1e-12
            found = true
            break
        end
        j += 1
    end
    !found && return true

    # Normal vector of plane
    n_flat = n_vec / n_norm

    # Check if all remaining vertices lie on the plane
    for k in 1:n
        vk = @view mesh.vertices[k, :]
        if abs(dot(vk - v1, n_flat)) > 1e-6
            return false
        end
    end
    return true
end

# Interaction Dispatch Strategy for Coated Components
function _resolve_coated_splitting_context(system, obj, ray)
    shape_obj = shape(obj)
    p_hit = world_to_local(shape_obj, position(ray))
    n_hit = world_to_local(shape_obj, normal3d(intersection(ray)))
    c_model = get_matching_coating(coatings(obj), shape_obj, p_hit, n_hit)

    # Detect if ray is entering or exiting the shape
    from_front = dot(direction(ray), normal3d(intersection(ray))) < 0.0
    n_incident = from_front ? refractive_index(system, ray) : refractive_index(obj, wavelength(ray))
    n_transmitted = from_front ? refractive_index(obj, wavelength(ray)) : refractive_index(system, ray)

    return c_model, n_incident, n_transmitted, normal3d(intersection(ray)), from_front
end

function resolve_coated_boundary(system, obj, ray)
    shape_obj = shape(obj)
    p_hit = world_to_local(shape_obj, position(ray))
    n_hit = world_to_local(shape_obj, normal3d(intersection(ray)))
    coating_model = get_matching_coating(coatings(obj), shape_obj, p_hit, n_hit)
    behavior = coating_behavior(coating_model, ray)

    return coating_model, behavior
end

function _interact3d_reflective_component_beams(system, obj, coating_model, agb, ray_id)
    c_ray = agb.c.rays[ray_id]
    shape_obj = shape(obj)
    p_hit = world_to_local(shape_obj, position(c_ray))
    n_hit = world_to_local(shape_obj, normal3d(intersection(c_ray)))
    from_front = dot(direction(c_ray), normal3d(intersection(c_ray))) < 0.0
    n1 = from_front ? refractive_index(system, c_ray) : refractive_index(obj, wavelength(c_ray))
    n2 = from_front ? refractive_index(obj, wavelength(c_ray)) : refractive_index(system, c_ray)
    normal = normal3d(intersection(c_ray))

    refl_agb = copy(agb)
    refl_agb.c.rays[ray_id] = reflect(c_ray, normal)
    return (agb, refl_agb)
end

function interact3d_splitting_polarized_ray(system, obj, ray)
    c_model, n_incident, n_transmitted, normal, from_front =
        _resolve_coated_splitting_context(system, obj, ray)

    return _propagate_splitting_polarized_ray(
        system, obj, c_model, ray, n_incident, n_transmitted, normal, from_front
    )
end

function interact3d_splitting_astigmatic_beamlet(system, obj, agb, ray_id)
    c_ray = agb.c.rays[ray_id]
    coating_model, n_incident, n_transmitted, normal, from_front =
        _resolve_coated_splitting_context(system, obj, c_ray)

    return _propagate_splitting_astigmatic_beamlet(
        system,
        obj,
        coating_model,
        agb,
        ray_id,
        n_incident,
        n_transmitted,
        normal,
        from_front,
        () -> begin
            ints = _interact3d_reflective_component_beams(system, obj, coating_model, agb, ray_id)
            isnothing(ints) && return nothing
            T_type = eltype(position(c_ray))
            return AstigmaticGaussianBeamletInteraction{T_type}(ints...)
        end
    )
end

# Fluent API / Pipeline helper
"""
    with_coatings(coatings::Pair...; deepcopy_shape=false)
    with_coatings(; front=nothing, back=nothing, deepcopy_shape=false)
    with_coatings(optic, coatings::Pair...; deepcopy_shape=false)
    with_coatings(optic; front=nothing, back=nothing, deepcopy_shape=false)

Fluent API helper to attach coatings to an optical component natively.
Returns a new instance of the same optic type with coatings attached.
When `deepcopy_shape=true`, deepcopies the inner shape before attaching coatings to avoid shared mutation.

# Examples
```julia
with_coatings(lens, :front => AR, :back => HR)
with_coatings(lens, (:front, :back) => SimpleARCoating())
with_coatings(mirror, :either => HR; deepcopy_shape=true)
```
"""
function with_coatings(coatings::Pair...; deepcopy_shape::Bool = false)
    obj -> with_coatings(obj, coatings...; deepcopy_shape = deepcopy_shape)
end

function with_coatings(obj::AbstractObject, coatings::Pair...; deepcopy_shape::Bool = false)
    mapped = map(c -> (c.first => _coating_model(c.second)), coatings)
    return _attach_coatings(obj, mapped; deepcopy_shape = deepcopy_shape)
end

function with_coatings(; front = nothing, back = nothing, deepcopy_shape::Bool = false)
    obj -> with_coatings(obj; front = front, back = back, deepcopy_shape = deepcopy_shape)
end

function with_coatings(obj::AbstractObject; front = nothing, back = nothing, deepcopy_shape::Bool = false)
    mapped = _build_coatings_tuple(shape(obj), front, back)
    return _attach_coatings(obj, mapped; deepcopy_shape = deepcopy_shape)
end

_attach_coatings(lens::Lens, c_tuple; deepcopy_shape::Bool = false) =
    Lens(deepcopy_shape ? deepcopy(lens.shape) : lens.shape, lens.n, c_tuple)
_attach_coatings(prism::Prism, c_tuple; deepcopy_shape::Bool = false) =
    Prism(deepcopy_shape ? deepcopy(prism.shape) : prism.shape, prism.n, c_tuple)
_attach_coatings(mirror::Mirror, c_tuple; deepcopy_shape::Bool = false) =
    Mirror(deepcopy_shape ? deepcopy(mirror.shape) : mirror.shape, c_tuple)

# Dispatch helpers for extracting coating models
_coating_model(c::AbstractCoating) = c.model
_coating_model(c) = c

# Shared coatings tuple constructor helper
function _build_coatings_tuple(parent_shape::AbstractShape, front, back)
    if is_flat_shape(parent_shape)
        if !isnothing(front) && !isnothing(back)
            throw(ArgumentError("Cannot specify both front and back coatings for a flat (2D) shape."))
        end
    end

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

# Active constituent SDF helper for CSG / composite shapes (allocation-free)
active_constituent_sdf(s::AbstractShape, p_hit) = s

function active_constituent_sdf(s::UnionSDF, p_hit)
    min_d = typemax(Float64)
    best_idx = 1
    for i in 1:length(s.sdfs)
        d = abs(sdf(s.sdfs[i], p_hit))
        if d < min_d
            min_d = d
            best_idx = i
        end
    end
    return active_constituent_sdf(s.sdfs[best_idx], p_hit)
end

# Get matching coating for a hit
function get_coating_model_at_hit(obj::AbstractObject, ray::AbstractRay)
    obj_coatings = coatings(obj)
    int = intersection(ray)
    isnothing(int) && return Uncoated()
    p_hit = position(ray) + length(ray) * direction(ray)

    parent_shape = shape(obj)
    local_p = world_to_local(parent_shape, p_hit)
    n_world = normal3d(int)
    if object(int) !== obj
        n_world = -n_world
    end
    local_n = transposed_orientation(parent_shape) * n_world

    hit_sub_sdf = isnothing(int.shape) ? active_constituent_sdf(parent_shape, p_hit) : int.shape
    coat_list = hasproperty(hit_sub_sdf, :coatings) ? getproperty(hit_sub_sdf, :coatings) : obj_coatings
    return get_matching_coating(coat_list, parent_shape, local_p, local_n)
end

function get_coating_model_at_hit(
        obj::AbstractObject, gauss::GaussianBeamlet, ray_id::Int)
    return get_coating_model_at_hit(obj, gauss.chief.rays[ray_id])
end

function get_coating_model_at_hit(
        obj::AbstractObject, agb::AstigmaticGaussianBeamlet, ray_id::Int)
    return get_coating_model_at_hit(obj, agb.c.rays[ray_id])
end

# Coincident boundary resolution dispatch helpers
function resolve_coated_boundary(
        system::AbstractSystem, obj::AbstractObject, ray::AbstractRay)
    coating = get_coating_model_at_hit(obj, ray)
    if !(coating isa Uncoated)
        return obj, coating
    end
    return resolve_coincident_coatings(intersection(ray), system, ray)
end

function resolve_coincident_coatings(::Nothing, system::AbstractSystem, ray::AbstractRay)
    (nothing, Uncoated())
end

function resolve_coincident_coatings(
        int::Intersection, system::AbstractSystem, ray::AbstractRay)
    res1 = check_coincident_coating(int.coincident_object, ray)
    res2 = check_coincident_coating(int.coincident_object_2, ray)
    if res1 !== nothing && res2 !== nothing
        obj1, coat1 = res1
        obj2, coat2 = res2
        if typeof(coat1) != typeof(coat2) || coat1 != coat2
            @warn "Conflicting coatings detected at coincident boundary between $(typeof(obj1)) and $(typeof(obj2)): $(typeof(coat1)) vs $(typeof(coat2)). Using $(typeof(coat1))."
        end
    end
    res1 !== nothing && return res1
    res2 !== nothing && return res2
    return nothing, Uncoated()
end

check_coincident_coating(::Nothing, ray::AbstractRay) = nothing

function check_coincident_coating(obj::AbstractObject, ray::AbstractRay)
    coating = get_coating_model_at_hit(obj, ray)
    if !(coating isa Uncoated)
        return (obj, coating)
    end
    return nothing
end

@inline function _interact3d_component_beams(system::AbstractSystem, target::AbstractObject, beamlet, ray_id::Int)
    beams = _component_beams(beamlet)
    ints = map(b -> interact3d(system, target, b, rays(b)[ray_id]), beams)
    return any(isnothing, ints) ? nothing : ints
end

@inline function _interact3d_reflective_component_beams(system::AbstractSystem, target::AbstractObject, coating_model, beamlet, ray_id::Int)
    beams = _component_beams(beamlet)
    ints = map(b -> interact3d_reflective(system, target, coating_model, b, rays(b)[ray_id]), beams)
    return any(isnothing, ints) ? nothing : ints
end

# interact3d definitions for beamlet propagation using trait dispatch
function interact3d_behavior(::Splitting, system::AbstractSystem, obj::AbstractObject,
        coating_model, gauss::GaussianBeamlet, ray_id::Int)
    return interact_splitting_boundary(system, obj, coating_model, gauss, ray_id)
end

function interact3d_behavior(
        ::CoatingBehavior, system::AbstractSystem, obj::AbstractObject,
        coating_model, gauss::GaussianBeamlet, ray_id::Int)
    ints = _interact3d_component_beams(system, obj, gauss, ray_id)
    isnothing(ints) && return nothing
    T = eltype(position(gauss.chief.rays[ray_id]))
    return GaussianBeamletInteraction{T}(ints...)
end

function interact3d_behavior(::Splitting, system::AbstractSystem, obj::AbstractObject,
        coating_model, agb::AstigmaticGaussianBeamlet, ray_id::Int)
    return interact_splitting_boundary(system, obj, coating_model, agb, ray_id)
end

function interact3d_behavior(
        ::CoatingBehavior, system::AbstractSystem, obj::AbstractObject,
        coating_model, agb::AstigmaticGaussianBeamlet, ray_id::Int)
    ints = _interact3d_component_beams(system, obj, agb, ray_id)
    isnothing(ints) && return nothing
    T = eltype(position(agb.c.rays[ray_id]))
    return AstigmaticGaussianBeamletInteraction{T}(ints...)
end

function interact3d_reflective(system::AbstractSystem, obj::AbstractRefractiveOptic,
        coating_model, beam::AbstractBeam, ray::AbstractRay)
    int = intersection(ray)
    λ = wavelength(ray)
    n_substrate = refractive_index(obj, λ)

    is_substrate = (_base_optic(object(int)) === obj)

    normal = normal3d(int)
    if !is_substrate
        normal = -normal
    end

    entering_substrate = dot(direction(ray), normal) < 0
    from_front = entering_substrate

    if entering_substrate
        n_incident = refractive_index(ray)
        n_transmitted = n_substrate
        hint = Hint(obj)
    else
        n_incident = n_substrate
        coin_obj = is_substrate ? int.coincident_object : object(int)
        if !isnothing(coin_obj) && is_refractive(coin_obj)
            n_transmitted = refractive_index(coin_obj, λ)
            hint = Hint(coin_obj)
        else
            n_transmitted = refractive_index(system, λ)
            hint = nothing
        end
        normal = -normal
    end
    return interact_refractive_boundary(
        Reflective(), system, obj, coating_model, beam, ray,
        n_incident, n_transmitted, hint, normal, λ, from_front)
end

function interact3d_reflective(system::AbstractSystem, obj::AbstractReflectiveOptic,
        coating_model, beam::AbstractBeam, ray::AbstractRay)
    return interact_reflective_boundary(system, obj, coating_model, beam, ray)
end

@inline function _resolve_coated_splitting_context(
        system::AbstractSystem, obj::AbstractObject, c_ray::AbstractRay)
    int = intersection(c_ray)
    λ = wavelength(c_ray)
    n_substrate = is_refractive(obj) ? refractive_index(obj, λ) : refractive_index(system, λ)

    is_target_obj = (_base_optic(object(int)) === obj)
    normal_coated = normal3d(int)
    if !is_target_obj
        normal_coated = -normal_coated
    end

    entering_coated = dot(direction(c_ray), normal_coated) < 0

    if entering_coated
        n_incident = refractive_index(c_ray)
        n_transmitted = n_substrate
        normal = normal_coated
    else
        n_incident = n_substrate
        other_obj = is_target_obj ? int.coincident_object : object(int)
        if !isnothing(other_obj) && is_refractive(other_obj)
            n_transmitted = refractive_index(other_obj, λ)
        else
            n_transmitted = refractive_index(system, λ)
        end
        normal = -normal_coated
    end

    from_front = entering_coated
    return (n_incident, n_transmitted, normal, from_front)
end

function interact_splitting_boundary(
        system::AbstractSystem,
        obj::AbstractObject,
        coating_model,
        gauss::GaussianBeamlet,
        ray_id::Int
)
    c_ray = rays(gauss.chief)[ray_id]
    n_incident, n_transmitted, normal, from_front =
        _resolve_coated_splitting_context(system, obj, c_ray)

    return _propagate_splitting_gaussian_beamlet(
        system,
        obj,
        coating_model,
        gauss,
        ray_id,
        n_incident,
        n_transmitted,
        normal,
        from_front,
        () -> begin
            ints = _interact3d_reflective_component_beams(system, obj, coating_model, gauss, ray_id)
            isnothing(ints) && return nothing
            T_type = eltype(position(c_ray))
            return GaussianBeamletInteraction{T_type}(ints...)
        end
    )
end

function interact_splitting_boundary(
        system::AbstractSystem,
        obj::AbstractObject,
        coating_model,
        agb::AstigmaticGaussianBeamlet,
        ray_id::Int
)
    # Ensure all component rays have valid intersections to prevent clipping crashes
    beams = _component_beams(agb)
    if any(b -> isnothing(intersection(rays(b)[ray_id])), beams)
        return nothing
    end

    c_ray = rays(agb.c)[ray_id]
    n_incident, n_transmitted, normal, from_front =
        _resolve_coated_splitting_context(system, obj, c_ray)

    return _propagate_splitting_astigmatic_beamlet(
        system,
        obj,
        coating_model,
        agb,
        ray_id,
        n_incident,
        n_transmitted,
        normal,
        from_front,
        () -> begin
            ints = _interact3d_reflective_component_beams(system, obj, coating_model, agb, ray_id)
            isnothing(ints) && return nothing
            T_type = eltype(position(c_ray))
            return AstigmaticGaussianBeamletInteraction{T_type}(ints...)
        end
    )
end
