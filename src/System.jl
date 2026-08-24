"""
    System{M<:AbstractMedium} <: AbstractSystem

A container storing the optical elements of, i.e. a camera lens or lab setup, and the surrounding ambient medium.

# Fields

- `objects`: vector containing the different objects that are part of the system (subtypes of [`AbstractObject`](@ref))
- `ambient_medium`: surrounding medium of the optical system (default: [`Ambient`](@ref))
"""
struct System{M <: AbstractMedium} <: AbstractSystem
    objects::Vector{AbstractObject}
    ambient_medium::M
end

System(objects::Vector{AbstractObject}) = System(objects, Ambient())
System(objects::Vector{<:AbstractObject}) = System(Vector{AbstractObject}(objects), Ambient())
System(object::AbstractObject) = System([object], Ambient())
System(objects::Vector{<:AbstractObject}, ambient::AbstractMedium) = System(Vector{AbstractObject}(objects), ambient)
System(objects::Vector{<:AbstractObject}, n::RefractiveIndex) = System(objects, medium_from(n))
System(object::AbstractObject, ambient::AbstractMedium) = System([object], ambient)
System(object::AbstractObject, n::RefractiveIndex) = System([object], medium_from(n))

ambient_medium(system::System) = system.ambient_medium

"""
    objects(system::System)

Exposes all objects stored within the system. By exposing the `Leaves` of the tree only, it is ensured that `AbstractObjectGroup`s are flattened into a regular vector.
"""
objects(system::System) = Leaves(system.objects)

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
StaticSystem(object::AbstractObjectGroup) = StaticSystem(tuple(collect(Leaves([object]))...), Ambient())
function StaticSystem(objects::AbstractArray{<:AbstractObject})
    StaticSystem(tuple(collect(Leaves(objects))...), Ambient())
end
function StaticSystem(objects::AbstractArray{<:AbstractObject}, ambient::AbstractMedium)
    StaticSystem(tuple(collect(Leaves(objects))...), ambient)
end

ambient_medium(system::StaticSystem) = system.ambient_medium

objects(system::StaticSystem) = system.objects

function trace_system!(::AbstractSystem, beam::B; r_max = 0) where {B <: AbstractBeam}
    @warn "Tracing for $B not implemented"
    return nothing
end

@inline function trace_all(system::AbstractSystem, ray::AbstractRay{R}) where {R}
    tol = get_coincident_boundary_tolerance()
    best_t = R(Inf)
    coincident_hits = Tuple{Intersection{R}, AbstractObject}[]
    
    # Iterate over all optical objects/interfaces in system
    for obj in objects(system)
        temp = intersect3d(obj, ray)
        if temp === nothing
            continue
        end
        t_hit = length(temp)

        if isempty(coincident_hits)
            best_t = t_hit
            push!(coincident_hits, (temp, obj))
        elseif abs(t_hit - best_t) <= tol
            push!(coincident_hits, (temp, obj))
            if t_hit < best_t
                best_t = t_hit
            end
        elseif t_hit < best_t - tol
            best_t = t_hit
            empty!(coincident_hits)
            push!(coincident_hits, (temp, obj))
        end
    end
    
    return _resolve_coincident_hits(coincident_hits, ray)
end

@inline function _resolve_coincident_hits(hits::Vector{Tuple{Intersection{R}, AbstractObject}}, ray::AbstractRay{R}) where {R}
    if isempty(hits)
        return nothing
    elseif length(hits) == 1
        return hits[1]
    end

    best_hit = hits[1][1]
    best_obj = hits[1][2]

    exiting_obj = nothing
    entering_obj = nothing
    coating_obj = nothing

    for (hit, obj) in hits
        if obj isa AbstractCoating || obj isa AbstractBeamsplitter
            coating_obj = obj
        else
            is_exit = dot(direction(ray), normal3d(hit)) > 0
            if is_exit
                exiting_obj = obj
            else
                entering_obj = obj
            end
        end
    end

    # Fallback if both/all are classified or ambiguous
    if exiting_obj === nothing && entering_obj === nothing
        if length(hits) >= 2
            h1, o1 = hits[1]
            h2, o2 = hits[2]
            h1_is_exit = dot(direction(ray), normal3d(h1)) > 0
            exiting_obj = h1_is_exit ? o1 : o2
            entering_obj = h1_is_exit ? o2 : o1
        end
    end

    mi = MultiIntersection(best_hit; exiting = exiting_obj, entering = entering_obj, coating = coating_obj)
    return (mi, best_obj)
end

@inline function trace_one(
        system::AbstractSystem, ray::AbstractRay{R}, hint::Hint) where {R}
    # Trace against hinted shape of object
    shape_intersection = intersect3d(shape(hint)::AbstractShape{R}, ray)
    if isnothing(shape_intersection)
        # If hinted object is not intersected, trace the entire system
        return trace_all(system, ray)
    end

    # Check if this hit is coincident with another object in system
    tol = get_coincident_boundary_tolerance()
    best_t = length(shape_intersection)
    coincident_hits = Tuple{Intersection{R}, AbstractObject}[(shape_intersection, object(hint))]

    for obj in objects(system)
        obj === object(hint) && continue
        temp = intersect3d(obj, ray)
        if temp !== nothing && abs(length(temp) - best_t) <= tol
            push!(coincident_hits, (temp, obj))
        end
    end

    return _resolve_coincident_hits(coincident_hits, ray)
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
    current_ray::AbstractRay{T} = first(beam.segments)
    current_opl::T = zero(T)
    empty!(beam.segments)

    while length(beam.segments) < r_max
        resolved_ray, hit, obj = tracing_step(system, current_ray, hint(interaction), current_opl)
        push!(beam, resolved_ray)

        if isnothing(hit) || isnothing(obj)
            break
        end

        current_opl = accumulated_opl(resolved_ray)
        # Dispatch to interface interaction depending on intersection type
        interaction = (hit isa MultiIntersection ?
            interact3d(system, hit, beam, resolved_ray) :
            interact3d(system, obj, beam, resolved_ray))::Union{Nothing, BeamInteraction}
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

    current_c::Ray{T} = first(gauss.chief.segments)
    current_w::Ray{T} = first(gauss.waist.segments)
    current_d::Ray{T} = first(gauss.divergence.segments)

    empty!(gauss.chief.segments)
    empty!(gauss.waist.segments)
    empty!(gauss.divergence.segments)

    opl_c::T = zero(T)
    opl_w::T = zero(T)
    opl_d::T = zero(T)

    while seg_counter <= r_max
        resolved_c, hit_c, obj_c = tracing_step(system, current_c, hint(interaction), opl_c)
        resolved_w, hit_w, obj_w = tracing_step(system, current_w, hint(interaction), opl_w)
        resolved_d, hit_d, obj_d = tracing_step(system, current_d, hint(interaction), opl_d)

        push!(gauss.chief, resolved_c)
        push!(gauss.waist, resolved_w)
        push!(gauss.divergence, resolved_d)

        int_c = hit_c
        int_w = hit_w
        int_d = hit_d

        # If any beam misses or beams do not hit same target stop tracing
        if any(isnothing, (int_c, int_w, int_d)) || (obj_c !== obj_w || obj_c !== obj_d) || !_beams_hits_same_shape(gauss, seg_counter)
            gauss.chief.segments[end] = Ray(current_c, OpenSegment(position(current_c), direction(current_c), opl_c))
            gauss.waist.segments[end] = Ray(current_w, OpenSegment(position(current_w), direction(current_w), opl_w))
            gauss.divergence.segments[end] = Ray(current_d, OpenSegment(position(current_d), direction(current_d), opl_d))
            break
        end

        opl_c = accumulated_opl(resolved_c)
        opl_w = accumulated_opl(resolved_w)
        opl_d = accumulated_opl(resolved_d)

        # Calculate interaction
        interaction = if int_c isa MultiIntersection
            interact3d(system, int_c, gauss, seg_counter)
        else
            interact3d(system, obj_c, gauss, seg_counter)
        end
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
        # kwargs
        r_max::Int = get_default_r_max(),
        check_invariant::Bool = true,
        threshold::Real = get_invariant_threshold()
) where {T <: Real}
    _reset_beam!(agb)
    interaction::Nullable{AstigmaticGaussianBeamletInteraction{T}} = nothing
    seg_counter::Int = 1
    aux = _aux_beams(agb)

    current_c::PolarizedRay{T} = first(agb.c.segments)
    current_aux = Ray{T}[first(b.segments) for b in aux]

    empty!(agb.c.segments)
    for b in aux
        empty!(b.segments)
    end

    opl_c::T = zero(T)
    opl_aux = zeros(T, length(aux))

    while seg_counter <= r_max
        resolved_c, hit_c, obj_c = tracing_step(system, current_c, hint(interaction), opl_c)
        push!(agb.c, resolved_c)

        missed = isnothing(hit_c)
        for (idx, b) in enumerate(aux)
            resolved_a, hit_a, obj_a = tracing_step(system, current_aux[idx], hint(interaction), opl_aux[idx])
            push!(b, resolved_a)
            if isnothing(hit_a) || obj_a !== obj_c
                missed = true
            end
        end

        if missed || !_beams_hits_same_shape(agb, seg_counter)
            agb.c.segments[end] = PolarizedRay(current_c, OpenSegment(position(current_c), direction(current_c), opl_c))
            for (idx, b) in enumerate(aux)
                b.segments[end] = Ray(current_aux[idx], OpenSegment(position(current_aux[idx]), direction(current_aux[idx]), opl_aux[idx]))
            end
            break
        end

        opl_c = accumulated_opl(resolved_c)
        for idx in 1:length(aux)
            opl_aux[idx] = accumulated_opl(aux[idx].segments[seg_counter])
        end

        int_c = hit_c
        interaction = if int_c isa MultiIntersection
            interact3d(system, int_c, agb, seg_counter)
        else
            interact3d(system, obj_c, agb, seg_counter)
        end
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
        # Reset the current beam back to its head ray(s) before (re-)tracing.
        _reset_beam!(current)
        # Process the current leaf beam.
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
