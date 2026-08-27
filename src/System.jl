"""
    _ObjectOrGroup

Internal union of everything that can be passed into a [`System`](@ref)/[`StaticSystem`](@ref):
a regular [`AbstractObject`](@ref), or an [`AbstractObjectGroup`](@ref) bundle thereof.
"""
const _ObjectOrGroup = Union{AbstractObject, AbstractObjectGroup}

"""
    _flatten_system_objects(objs) -> Vector{AbstractObject}

Recursively expands `objs` (a single [`_ObjectOrGroup`](@ref) or a collection thereof) into
a flat `Vector{AbstractObject}` suitable for tracing: [`AbstractObjectGroup`](@ref)s are
expanded via `objects(group)`, and any [`MultiShape`](@ref) `AbstractObject` whose parts
(`shape(obj)`) are themselves all `AbstractObject`/`AbstractObjectGroup` is expanded the
same way (e.g. a [`CubeBeamsplitter`](@ref) becomes its `front`/`back`/`coating`). Any
other object is kept as-is.
"""
function _flatten_system_objects(objs)
    acc = AbstractObject[]
    for o in objs
        _flatten_system_objects!(acc, o)
    end
    return acc
end
_flatten_system_objects(obj::_ObjectOrGroup) = _flatten_system_objects!(AbstractObject[], obj)

function _flatten_system_objects!(acc::Vector{AbstractObject}, obj::AbstractObject)
    if shape_trait_of(obj) isa MultiShape
        parts = shape(obj)
        if all(p -> p isa _ObjectOrGroup, parts)
            for p in parts
                _flatten_system_objects!(acc, p)
            end
            return acc
        end
    end
    push!(acc, obj)
    return acc
end
function _flatten_system_objects!(acc::Vector{AbstractObject}, group::AbstractObjectGroup)
    for o in objects(group)
        _flatten_system_objects!(acc, o)
    end
    return acc
end

"""
    System{M<:AbstractMedium} <: AbstractSystem

A container storing the optical elements of, i.e. a camera lens or lab setup, and the surrounding ambient medium.

# Fields

- `objects`: vector containing the objects that are part of the system — [`AbstractObject`](@ref)s
  and/or [`AbstractObjectGroup`](@ref) bundles thereof, exactly as passed to the constructor
- `ambient_medium`: surrounding medium of the optical system (default: [`Ambient`](@ref))
"""
struct System{M <: AbstractMedium} <: AbstractSystem
    objects::Vector{_ObjectOrGroup}
    ambient_medium::M
    flat_objects::Vector{AbstractObject}
end

function System(objects::Vector{_ObjectOrGroup}, ambient::M) where {M <: AbstractMedium}
    flat = _flatten_system_objects(objects)
    return System{M}(objects, ambient, flat)
end
System(objects::Vector{_ObjectOrGroup}) = System(objects, Ambient())
System(objects::AbstractVector) = System(Vector{_ObjectOrGroup}(objects), Ambient())
System(object::_ObjectOrGroup) = System([object], Ambient())
System(objects::AbstractVector, ambient::AbstractMedium) = System(Vector{_ObjectOrGroup}(objects), ambient)
System(objects::AbstractVector, n::RefractiveIndex) = System(objects, medium_from(n))
System(object::_ObjectOrGroup, ambient::AbstractMedium) = System([object], ambient)
System(object::_ObjectOrGroup, n::RefractiveIndex) = System([object], medium_from(n))

ambient_medium(system::System) = system.ambient_medium

"""
    objects(system::System)

Exposes all objects stored within the system as a flat, trace-ready list: any
[`AbstractObjectGroup`](@ref) bundle or [`MultiShape`](@ref) composite (e.g. a
[`CubeBeamsplitter`](@ref) or [`DoubletLens`](@ref)) is recursively expanded into its
constituent [`AbstractObject`](@ref)s, see [`_flatten_system_objects`](@ref).
"""
objects(system::System) = system.flat_objects

function Base.push!(system::System, obj::_ObjectOrGroup)
    push!(system.objects, obj)
    _flatten_system_objects!(system.flat_objects, obj)
    return system
end

function Base.empty!(system::System)
    empty!(system.objects)
    empty!(system.flat_objects)
    return system
end

"""
    StaticSystem{T<:Tuple, M<:AbstractMedium} <: AbstractSystem

A static container storing the optical elements of, i.e. a camera lens or lab setup.
Compared to `System` this way defining the system is less flexible, i.e. no elements
can be added or removed after construction but it allows for more performant ray-tracing.

!!! warning
    This type uses long tuples for storing the elements. This container should not be used
    for very large optical systems as it puts a lot of stress onto the compiler.

# Fields

- `objects`: tuple containing the different objects that are part of the system (subtypes of [`AbstractObject`](@ref))
- `ambient_medium`: surrounding medium of the optical system (default: [`Ambient`](@ref))
"""
struct StaticSystem{T <: Tuple, M <: AbstractMedium} <: AbstractSystem
    objects::T
    ambient_medium::M
end

StaticSystem(objects::T) where {T <: Tuple} = StaticSystem(objects, Ambient())
StaticSystem(object::AbstractObject) = StaticSystem((object,), Ambient())
StaticSystem(object::AbstractObjectGroup) = StaticSystem(tuple(_flatten_system_objects(object)...), Ambient())
function StaticSystem(objects::AbstractArray)
    StaticSystem(tuple(_flatten_system_objects(objects)...), Ambient())
end
function StaticSystem(objects::AbstractArray, ambient::AbstractMedium)
    StaticSystem(tuple(_flatten_system_objects(objects)...), ambient)
end

ambient_medium(system::StaticSystem) = system.ambient_medium

objects(system::StaticSystem) = system.objects

function trace_system!(::AbstractSystem, beam::B; r_max = 0) where {B <: AbstractBeam}
    @warn "Tracing for $B not implemented"
    return nothing
end

@noinline function _throw_too_many_coincident_hits()
    throw(ErrorException("Too many coincident hits detected at interface (> 3 objects)"))
end

@noinline function _throw_multiple_exiting_objects(o1, o2)
    throw(ErrorException("Multiple exiting objects detected at coincident interface: $(typeof(o1)) and $(typeof(o2))"))
end

@noinline function _throw_multiple_entering_objects(o1, o2)
    throw(ErrorException("Multiple entering objects detected at coincident interface: $(typeof(o1)) and $(typeof(o2))"))
end

"""
    trace_all(system::AbstractSystem, ray::AbstractRay{R}) where {R<:Real}

Performs non-sequential ray tracing across all objects in `system` to find the closest
geometrical intersection. Detects and groups coincident boundaries within tolerance.

# Returns
- `(intersection, object)`: A tuple of the resolved [`Intersection`](@ref) or
  [`MultiIntersection`](@ref) and the associated [`AbstractObject`](@ref).
- `nothing`: If no object in the system was hit.
"""
@inline function trace_all(system::AbstractSystem, ray::AbstractRay{R}) where {R <: Real}
    tol = get_coincident_boundary_tolerance()
    best_t = R(Inf)
    hits = ()
    
    # Iterate over all optical objects/interfaces in system
    for obj in objects(system)
        temp = intersect3d(obj, ray)
        temp === nothing && continue
        t_hit = length(temp)

        if isempty(hits)
            best_t = t_hit
            hits = ((temp, obj),)
        elseif abs(t_hit - best_t) <= tol
            length(hits) >= 3 && _throw_too_many_coincident_hits()
            hits = (hits..., (temp, obj))
            if t_hit < best_t
                best_t = t_hit
            end
        elseif t_hit < best_t - tol
            best_t = t_hit
            hits = ((temp, obj),)
        end
    end
    
    return _resolve_coincident_hits(hits, ray)
end

"""
    _resolve_coincident_hits(hits::Union{Tuple, AbstractVector}, ray::AbstractRay{R}) where {R<:Real}

Resolves multiple coincident geometric intersections into a unified [`MultiIntersection`](@ref),
classifying entering, exiting, and coating objects via ray propagation direction and outward surface normals.
"""
@inline function _resolve_coincident_hits(hits::Union{Tuple, AbstractVector}, ray::AbstractRay{R}) where {R <: Real}
    # No intersection detected
    if isempty(hits)
        return nothing
    end
    # Single intersection detected (fast path)
    if length(hits) == 1
        return hits[1]
    end
    # Only up to three objects allowed at coincident interface (e.g. 2x substrate, 1x coating)
    if length(hits) > 3
        _throw_too_many_coincident_hits()
    end

    best_hit = hits[1][1]
    best_obj = hits[1][2]

    exiting_obj = nothing
    entering_obj = nothing
    coating_obj = nothing
    dir = direction(ray)

    # Classify coating, entering object, and exiting object via ray direction and outward normal
    for (hit, obj) in hits
        if is_thin_interface(obj)
            # Zero-thickness interface, not a solid volume: cannot be entered/exited
            coating_obj = obj
        else
            # Ray exits a volume through a surface whose outward normal points in the same
            # direction as the ray (dot > 0); otherwise it enters
            is_exit = dot(dir, normal3d(hit)) > 0
            if is_exit
                if exiting_obj !== nothing
                    _throw_multiple_exiting_objects(exiting_obj, obj)
                end
                exiting_obj = obj
            else
                if entering_obj !== nothing
                    _throw_multiple_entering_objects(entering_obj, obj)
                end
                entering_obj = obj
            end
        end
    end

    # Fallback: All coincident objects classified as surface interfaces (or no paired exit/entry)
    if exiting_obj === nothing && entering_obj === nothing
        if length(hits) >= 2
            h1, o1 = hits[1]
            h2, o2 = hits[2]
            h1_is_exit = dot(dir, normal3d(h1)) > 0
            exiting_obj = h1_is_exit ? o1 : o2
            entering_obj = h1_is_exit ? o2 : o1
        end
    end

    # Construct MultiIntersection bundle for downstream interact3d dispatch
    mi = MultiIntersection(best_hit; exiting = exiting_obj, entering = entering_obj, coating = coating_obj)
    return (mi, best_obj)
end

"""
    trace_one(system::AbstractSystem, ray::AbstractRay{R}, hint::Hint) where {R<:Real}

Traces `ray` using an optional [`Hint`](@ref) for fast shape-intersection. If the hinted
object is hit, checks for coincident interfaces; otherwise falls back to [`trace_all`](@ref).
"""
@inline function trace_one(
        system::AbstractSystem, ray::AbstractRay{R}, hint::Hint) where {R <: Real}
    # Fast path: Trace against hinted shape of object
    shape_intersection = intersect3d(shape(hint)::AbstractShape{R}, ray)
    if isnothing(shape_intersection)
        return trace_all(system, ray)
    end

    # Check if this hit is coincident with other objects in the system
    tol = get_coincident_boundary_tolerance()
    best_t = length(shape_intersection)
    hits = ((shape_intersection, object(hint)),)

    for obj in objects(system)
        obj === object(hint) && continue
        temp = intersect3d(obj, ray)
        if temp !== nothing && abs(length(temp) - best_t) <= tol
            length(hits) >= 3 && _throw_too_many_coincident_hits()
            hits = (hits..., (temp, obj))
        end
    end

    return _resolve_coincident_hits(hits, ray)
end

@inline trace_one(system::AbstractSystem, ray::AbstractRay, ::Nothing) = trace_all(system, ray)


"""
    tracing_step(system::AbstractSystem, ray::AbstractRay{R}, hint::Nullable{Hint} = nothing, opl_accum::R = zero(R)) where {R <: Real}

Tests if the `ray` intersects an `object` in the optical `system`. Returns `(resolved_ray, hit_intersection, hit_object)` or `(resolved_ray, nothing, nothing)` if unhit.
Delegates volumetric propagation and segment construction to `propagate_volume`.
"""
@inline function tracing_step(
        system::AbstractSystem, ray::AbstractRay{T}, hint::Nullable{Hint} = nothing, opl_accum::T = zero(T)) where {T <: Real}
    med = current_medium(system, hint)
    return propagate_volume(system, med, ray, hint, opl_accum)
end

@inline function tracing_step!(
        system::AbstractSystem, ray::AbstractRay{R}, hint::Nullable{Hint} = nothing) where {R <: Real}
    res = trace_one(system, ray, hint)
    if res === nothing
        return nothing
    else
        hit, obj = res
        return obj
    end
end

# Dispatched interaction helper for hit vs multi-intersection
@inline interact_hit(system::AbstractSystem, obj::AbstractObject, hit::Intersection, beam, ray) = interact3d(system, obj, beam, ray)
@inline interact_hit(system::AbstractSystem, ::AbstractObject, mi::MultiIntersection, beam, ray) = interact3d(system, mi, beam, ray)

@inline interact_hit(system::AbstractSystem, obj::AbstractObject, hit::Intersection, beam, seg_counter::Int) = interact3d(system, obj, beam, seg_counter)
@inline interact_hit(system::AbstractSystem, ::AbstractObject, mi::MultiIntersection, beam, seg_counter::Int) = interact3d(system, mi, beam, seg_counter)

"""
    trace_system!(system::AbstractSystem, beam::Beam{T}; r_max = get_default_r_max()) where {T <: Real}

Trace a [`Beam`](@ref) through an optical `system`. Maximum number of tracing steps can be capped by `r_max`.
"""
function trace_system!(
        system::AbstractSystem,
        beam::Beam{T};
        # kwargs
        r_max::Int = get_default_r_max(),
        kwargs...
) where {T <: Real}
    _reset_beam!(beam)
    interaction::Nullable{BeamInteraction} = nothing
    current_ray::AbstractRay{T} = first(beam.rays)
    current_opl::T = zero(T)
    empty!(beam.rays)

    while length(beam.rays) < r_max
        resolved_ray, hit, obj = tracing_step(system, current_ray, hint(interaction), current_opl)
        push!(beam, resolved_ray)

        if isnothing(hit) || isnothing(obj)
            break
        end

        current_opl = accumulated_opl(resolved_ray)
        # Dispatch to interface interaction depending on intersection type
        interaction = interact_hit(system, obj, hit, beam, resolved_ray)::Union{Nothing, BeamInteraction}
        if isnothing(interaction)
            break
        end
        current_ray = interaction.ray
    end
    return nothing
end

"""
    trace_system!(system::System, gauss::GaussianBeamlet{T}; r_max = get_default_r_max()) where {T <: Real}

Trace a [`GaussianBeamlet`](@ref) through an optical `system`. Maximum number of tracing steps can be capped by `r_max`.
"""
function trace_system!(
        system::AbstractSystem,
        gauss::GaussianBeamlet{T};
        # kwargs...
        r_max::Int = get_default_r_max(),
        kwargs...
) where {T <: Real}
    _reset_beam!(gauss)
    interaction::Nullable{GaussianBeamletInteraction{T}} = nothing
    seg_counter::Int = 1

    current_c::Ray{T} = first(gauss.chief.rays)
    current_w::Ray{T} = first(gauss.waist.rays)
    current_d::Ray{T} = first(gauss.divergence.rays)

    opl_c::T = zero(T)
    opl_w::T = zero(T)
    opl_d::T = zero(T)

    empty!(gauss.chief.rays)
    empty!(gauss.waist.rays)
    empty!(gauss.divergence.rays)

    while length(gauss.chief.rays) < r_max
        resolved_c, hit_c, obj_c = tracing_step(system, current_c, hint(interaction), opl_c)
        resolved_w, hit_w, obj_w = tracing_step(system, current_w, hint(interaction), opl_w)
        resolved_d, hit_d, obj_d = tracing_step(system, current_d, hint(interaction), opl_d)

        push!(gauss.chief, resolved_c)
        push!(gauss.waist, resolved_w)
        push!(gauss.divergence, resolved_d)

        int_c, int_w, int_d = hit_c, hit_w, hit_d

        if any(isnothing, (int_c, int_w, int_d)) || (obj_c !== obj_w || obj_c !== obj_d) || !_beams_hits_same_shape(gauss, seg_counter)
            # Finish rays as OpenSegments
            gauss.chief.rays[end] = with_segment(current_c, OpenSegment(position(current_c), direction(current_c), opl_c))
            gauss.waist.rays[end] = with_segment(current_w, OpenSegment(position(current_w), direction(current_w), opl_w))
            gauss.divergence.rays[end] = with_segment(current_d, OpenSegment(position(current_d), direction(current_d), opl_d))
            break
        end

        opl_c = accumulated_opl(resolved_c)
        opl_w = accumulated_opl(resolved_w)
        opl_d = accumulated_opl(resolved_d)

        # Calculate interaction
        interaction = interact_hit(system, obj_c, int_c, gauss, seg_counter)
        if isnothing(interaction)
            break
        end

        current_c = interaction.chief.ray
        current_w = interaction.waist.ray
        current_d = interaction.divergence.ray
        seg_counter += 1
    end
    return nothing
end

"""
    trace_system!(system, agb::AstigmaticGaussianBeamlet; r_max = get_default_r_max(), check_invariant = true, threshold = get_invariant_threshold())

Trace an [`AstigmaticGaussianBeamlet`](@ref) through an optical `system`.
"""
function trace_system!(
        system::AbstractSystem,
        agb::AstigmaticGaussianBeamlet{T};
        # kwargs...
        r_max::Int = get_default_r_max(),
        check_invariant::Bool = true,
        threshold::Real = get_invariant_threshold(),
        kwargs...
) where {T <: Real}
    _reset_beam!(agb)
    interaction::Nullable{AstigmaticGaussianBeamletInteraction{T}} = nothing
    seg_counter::Int = 1

    current_c::PolarizedRay{T} = first(agb.c.rays)
    current_aux = Vector{Ray{T}}(undef, 8)
    aux = _aux_beams(agb)
    for idx in 1:length(aux)
        current_aux[idx] = first(aux[idx].rays)
    end

    opl_c::T = zero(T)
    opl_aux = zeros(T, 8)

    empty!(agb.c.rays)
    for beam in aux
        empty!(beam.rays)
    end

    while length(agb.c.rays) < r_max
        resolved_c, hit_c, obj_c = tracing_step(system, current_c, hint(interaction), opl_c)
        push!(agb.c, resolved_c)

        missed = isnothing(hit_c)
        for idx in 1:length(aux)
            resolved_aux, hit_a, obj_a = tracing_step(system, current_aux[idx], hint(interaction), opl_aux[idx])
            push!(aux[idx], resolved_aux)
            if isnothing(hit_a) || obj_a !== obj_c
                missed = true
            end
        end

        if missed || !_beams_hits_same_shape(agb, seg_counter)
            agb.c.rays[end] = with_segment(current_c, OpenSegment(position(current_c), direction(current_c), opl_c))
            for idx in 1:length(aux)
                aux[idx].rays[end] = with_segment(current_aux[idx], OpenSegment(position(current_aux[idx]), direction(current_aux[idx]), opl_aux[idx]))
            end
            break
        end

        opl_c = accumulated_opl(resolved_c)
        for idx in 1:length(aux)
            opl_aux[idx] = accumulated_opl(aux[idx].rays[seg_counter])
        end

        int_c = hit_c
        interaction = interact_hit(system, obj_c, int_c, agb, seg_counter)
        if isnothing(interaction)
            break
        end

        current_c = interaction.chief.ray
        current_aux[1] = interaction.wxp.ray
        current_aux[2] = interaction.wxm.ray
        current_aux[3] = interaction.wyp.ray
        current_aux[4] = interaction.wym.ray
        current_aux[5] = interaction.dxp.ray
        current_aux[6] = interaction.dxm.ray
        current_aux[7] = interaction.dyp.ray
        current_aux[8] = interaction.dym.ray

        # Verify that the paraxial assumptions still hold for the current segment
        if check_invariant &&
           !check_optical_invariant(agb, seg_counter; threshold = threshold)
            break
        end

        seg_counter += 1
    end
    return nothing
end



"""
    solve_system!(system::System, beam::AbstractBeam; r_max=100, depth_max=typemax(Int))

Manage the tracing of an `AbstractBeam` through an optical `system`. Each beam in the tree
is first reset back to just its head ray(s) (see [`BeamletOptics._reset_beam!`](@ref)) and
then traced fresh, so calling this again on an already-solved `beam` (e.g. after moving an
object, or mutating a beam/ray property) performs a full non-sequential trace.
The condition to stop ray tracing is that the last `beam` intersection is `nothing` or the beam interaction is `nothing`. Then, the system is considered to be solved.
A maximum number of rays per `beam` (`r_max`) can be specified in order to avoid infinite calculations under resonant conditions, i.e. two facing mirrors. Likewise, `depth_max` limits how many branching levels are explored when new sub-beams are generated (for example, by beamsplitters) so that the tree cannot grow without bound.

# Arguments

- `system::System`: The optical system in which the beam will be traced.
- `beam::AbstractBeam`: The beam object to be traced through the system.

## Keyword Arguments

- `r_max = get_default_r_max()`: Maximum number of tracing iterations for each leaf.
- `depth_max = get_default_depth_max()`: Maximum number of branching levels explored from the root beam.
- `check_invariant = true`: enables or disables optical invariant checks where applicable
- `threshold = get_invariant_threshold()`: threshold for paraxial invariant checks
"""
function solve_system!(
        system::AbstractSystem,
        beam::B;
        # kwargs
        r_max::Int = get_default_r_max(),
        depth_max::Int = get_default_depth_max(),
        check_invariant::Bool = true,
        threshold::Real = get_invariant_threshold()
) where {B <: AbstractBeam}
    queue = Tuple{B, Int}[(beam, 1)]
    while !isempty(queue)
        # Process beams in FIFO order.
        current, depth = popfirst!(queue)
        # Process the current leaf beam. `trace_system!` resets the beam back to its
        # head ray(s) before (re-)tracing (see `BeamletOptics._reset_beam!`).
        solve_leaf!(system, current; r_max, check_invariant, threshold)
        # Check if the maximum branching depth has been reached.
        if depth <= depth_max
            # Enqueue all child beams for subsequent processing.
            for child in children(current) # 'children' returns an iterable of sub-beams.
                push!(queue, (child, depth + 1))
            end
        else
            # Maximum braching depth is reached. Remove childrens of the current beam because they will not be solved.
            _drop_beams!(current)
            @debug lazy"Maximum branching depth of $depth_max levels reached."
        end
    end
    return nothing
end

function solve_system!(system::AbstractSystem, bg::AbstractBeamGroup; kwargs...)
    Threads.@threads for _beam in beams(bg)
        solve_system!(system, _beam; kwargs...)
    end
    return nothing
end

@inline function solve_leaf!(
        system::AbstractSystem, beam::AbstractBeam; r_max = get_default_r_max(), kwargs...)
    trace_system!(system, beam; r_max, kwargs...)
    return nothing
end

function AbstractTrees.printnode(io::IO, node::B; kw...) where {B <: AbstractObject}
    show(io, node)
end
function AbstractTrees.printnode(io::IO,
        node::B;
        kw...) where {B <: AbstractObjectGroup}
    show(io, node)
end

Base.show(::IO, ::MIME"text/plain", system::System) =
    for obj in system.objects
        print_tree(obj)
    end
