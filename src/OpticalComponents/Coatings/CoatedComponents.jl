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
        if n_norm >= 1e-12
            found = true
            break
        end
        j += 1
    end
    !found && return true  # All vertices are collinear (line of points, flat)

    # Check if all other vertices lie in the defined plane
    n_vec = n_vec / n_norm
    for k in 1:n
        vk = @view mesh.vertices[k, :]
        if abs(dot(vk - v1, n_vec)) > 1e-6
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

_base_optic(cl::CoatedComponent) = cl.optic

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
is_refractive(::CoatedMirror) = false
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

function get_coating_model_at_hit(
        coated::CoatedComponent, gauss::GaussianBeamlet, ray_id::Int)
    return get_coating_model_at_hit(coated, gauss.chief.rays[ray_id])
end

function get_coating_model_at_hit(
        coated::CoatedComponent, agb::AstigmaticGaussianBeamlet, ray_id::Int)
    return get_coating_model_at_hit(coated, agb.c.rays[ray_id])
end

# Coincident boundary resolution dispatch helpers
function resolve_coated_boundary(
        system::AbstractSystem, obj::AbstractObject, ray::AbstractRay)
    resolve_coated_boundary_dispatch(obj, system, ray)
end

function resolve_coated_boundary_dispatch(
        obj::CoatedComponent, system::AbstractSystem, ray::AbstractRay)
    coating = get_coating_model_at_hit(obj, ray)
    if !(coating isa Uncoated)
        return obj, coating
    end
    return resolve_coincident_coatings(intersection(ray), system, ray)
end

function resolve_coated_boundary_dispatch(
        ::AbstractObject, system::AbstractSystem, ray::AbstractRay)
    resolve_coincident_coatings(intersection(ray), system, ray)
end

function resolve_coincident_coatings(::Nothing, system::AbstractSystem, ray::AbstractRay)
    (nothing, Uncoated())
end

function resolve_coincident_coatings(
        int::Intersection, system::AbstractSystem, ray::AbstractRay)
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
    if coating_behavior(coating_model, ray) isa Absorptive
        return nothing
    end
    return interact_refractive_boundary(system, cl.optic, coating_model, beam, ray)
end

function interact3d(system::AbstractSystem,
        cm::CoatedMirror,
        beam::Beam{T, R},
        ray::R) where {T <: Real, R <: AbstractRay{T}}
    coating_model = get_coating_model_at_hit(cm, ray)
    if coating_behavior(coating_model, ray) isa Absorptive
        return nothing
    end
    return interact_reflective_boundary(system, cm.optic, coating_model, beam, ray)
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

# interact3d definitions for CoatedComponent beamlet propagation using trait dispatch
function interact3d(
        system::AbstractSystem, cl::CoatedComponent, gauss::GaussianBeamlet, ray_id::Int)
    coating_model = get_coating_model_at_hit(cl, gauss, ray_id)
    chief_ray = rays(gauss.chief)[ray_id]
    return interact3d_behavior(coating_behavior(coating_model, chief_ray),
        system, cl, coating_model, gauss, ray_id)
end

function interact3d_behavior(::Splitting, system::AbstractSystem, cl::CoatedComponent,
        coating_model, gauss::GaussianBeamlet, ray_id::Int)
    return interact_splitting_boundary(system, cl, coating_model, gauss, ray_id)
end

function interact3d_behavior(
        ::CoatingBehavior, system::AbstractSystem, cl::CoatedComponent,
        coating_model, gauss::GaussianBeamlet, ray_id::Int)
    ints = _interact3d_component_beams(system, cl, gauss, ray_id)
    isnothing(ints) && return nothing
    T = eltype(position(gauss.chief.rays[ray_id]))
    return GaussianBeamletInteraction{T}(ints...)
end

function interact3d(
        system::AbstractSystem, cl::CoatedComponent, agb::AstigmaticGaussianBeamlet, ray_id::Int)
    coating_model = get_coating_model_at_hit(cl, agb, ray_id)
    chief_ray = rays(agb.c)[ray_id]
    return interact3d_behavior(
        coating_behavior(coating_model, chief_ray), system, cl, coating_model, agb, ray_id)
end

function interact3d_behavior(::Splitting, system::AbstractSystem, cl::CoatedComponent,
        coating_model, agb::AstigmaticGaussianBeamlet, ray_id::Int)
    return interact_splitting_boundary(system, cl, coating_model, agb, ray_id)
end

function interact3d_behavior(
        ::CoatingBehavior, system::AbstractSystem, cl::CoatedComponent,
        coating_model, agb::AstigmaticGaussianBeamlet, ray_id::Int)
    ints = _interact3d_component_beams(system, cl, agb, ray_id)
    isnothing(ints) && return nothing
    T = eltype(position(agb.c.rays[ray_id]))
    return AstigmaticGaussianBeamletInteraction{T}(ints...)
end

function interact3d_reflective(system::AbstractSystem, cl::CoatedRefractive,
        coating_model, beam::AbstractBeam, ray::AbstractRay)
    int = intersection(ray)
    λ = wavelength(ray)
    n_substrate = refractive_index(cl.optic, λ)

    obj = object(int)
    is_substrate = (_base_optic(obj) === cl.optic)

    normal = normal3d(int)
    if !is_substrate
        normal = -normal
    end

    entering_substrate = dot(direction(ray), normal) < 0
    from_front = entering_substrate

    if entering_substrate
        n_incident = refractive_index(ray)
        n_transmitted = n_substrate
        hint = Hint(cl.optic)
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
    return interact_refractive_boundary(
        Reflective(), system, cl.optic, coating_model, beam, ray,
        n_incident, n_transmitted, hint, normal, λ, from_front)
end

function interact3d_reflective(system::AbstractSystem, cl::CoatedMirror,
        coating_model, beam::AbstractBeam, ray::AbstractRay)
    return interact_reflective_boundary(system, cl.optic, coating_model, beam, ray)
end

function interact_splitting_boundary(
        system::AbstractSystem,
        coated::CoatedComponent,
        coating_model,
        gauss::GaussianBeamlet,
        ray_id::Int
)
    c_ray = rays(gauss.chief)[ray_id]
    int = intersection(c_ray)
    λ = wavelength(c_ray)
    substrate_obj = coated.optic
    n_substrate = refractive_index(substrate_obj, λ)

    is_coated_obj = (_base_optic(object(int)) === coated.optic)
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

    return _propagate_splitting_gaussian_beamlet(
        system,
        coated,
        coating_model,
        gauss,
        ray_id,
        n_incident,
        n_transmitted,
        normal,
        from_front,
        () -> begin
            ints = _interact3d_reflective_component_beams(system, coated, coating_model, gauss, ray_id)
            isnothing(ints) && return nothing
            T_type = eltype(position(c_ray))
            return GaussianBeamletInteraction{T_type}(ints...)
        end
    )
end

function interact_splitting_boundary(
        system::AbstractSystem,
        coated::CoatedComponent,
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
    int = intersection(c_ray)
    λ = wavelength(c_ray)
    substrate_obj = coated.optic
    n_substrate = refractive_index(substrate_obj, λ)

    is_coated_obj = (_base_optic(object(int)) === coated.optic)
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

    return _propagate_splitting_astigmatic_beamlet(
        system,
        coated,
        coating_model,
        agb,
        ray_id,
        n_incident,
        n_transmitted,
        normal,
        from_front,
        () -> begin
            ints = _interact3d_reflective_component_beams(system, coated, coating_model, agb, ray_id)
            isnothing(ints) && return nothing
            T_type = eltype(position(c_ray))
            return AstigmaticGaussianBeamletInteraction{T_type}(ints...)
        end
    )
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

# Base.show methods for coated components
function Base.show(io::IO, ::MIME"text/plain", cr::CoatedRefractive)
    print(io, "CoatedRefractive(", cr.optic, ", coatings = ", cr.coatings, ")")
end
function Base.show(io::IO, ::MIME"text/plain", cm::CoatedMirror)
    print(io, "CoatedMirror(", cm.optic, ", coatings = ", cm.coatings, ")")
end
