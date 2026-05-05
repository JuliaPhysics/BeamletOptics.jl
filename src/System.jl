"""
    System <: AbstractSystem

A container storing the optical elements of, i.e. a camera lens or lab setup.

# Fields

- `objects`: vector containing the different objects that are part of the system (subtypes of [`AbstractObject`](@ref))
"""
struct System <: AbstractSystem
    objects::Vector{AbstractObject}
end

System(object::AbstractObject) = System([object])

"""
    objects(system::System)

Exposes all objects stored within the system. By exposing the `Leaves` of the tree only, it is ensured that `AbstractObjectGroup`s are flattened into a regular vector.
"""
objects(system::System) = Leaves(system.objects)

"""
    StaticSystem <: AbstractSystem

A static container storing the optical elements of, i.e. a camera lens or lab setup.
Compared to `System` this way defining the system is less flexible, i.e. no elements
can be added or removed after construction but it allows for more performant ray-tracing.

!!! warning
    This type uses long tuples for storing the elements. This container should not be used
    for very large optical systems as it puts a lot of stress onto the compiler.

# Fields

- `objects`: vector containing the different objects that are part of the system (subtypes of [`AbstractObject`](@ref))
"""
struct StaticSystem{T <: Tuple} <: AbstractSystem
    objects::T
end
StaticSystem(object::AbstractObject) = StaticSystem((object))
StaticSystem(object::AbstractObjectGroup) = StaticSystem([object])
function StaticSystem(objects::AbstractArray{<:AbstractObject})
    StaticSystem(tuple(collect(Leaves(objects))...))
end

objects(system::StaticSystem) = system.objects

function trace_system!(::AbstractSystem, beam::B; r_max = 0) where {B <: AbstractBeam}
    @warn "Tracing for $B not implemented"
    return nothing
end

function retrace_system!(::AbstractSystem, beam::B) where {B <: AbstractBeam}
    @warn "Retracing for $B not implemented"
    return nothing
end

@inline function trace_all(system::AbstractSystem, ray::AbstractRay{R}) where {R}
    result::Union{Nothing, Intersection{R}} = nothing
    for obj in objects(system)
        # Find shortest intersection
        temp::Union{Nothing, Intersection{R}} = intersect3d(obj, ray)
        if temp === nothing
            continue
        end

        # Catch first valid intersection and replace current with closer intersection
        if result === nothing || length(temp) < length(result)
            result = temp
        end
    end
    return result
end

@inline function trace_one(
        system::AbstractSystem, ray::AbstractRay{R}, hint::Hint) where {R}
    # Trace against hinted shape of object
    intersection::Nullable{Intersection{R}} = intersect3d(
        shape(hint)::AbstractShape{R}, ray)
    if isnothing(intersection)
        # If hinted object is not intersected, trace the entire system
        intersection = trace_all(system, ray)
    else
        # If hinted object is intersected, update intersection
        object!(intersection, object(hint))
    end
    return intersection
end

"""
    tracing_step!(system::AbstractSystem, ray::AbstractRay{R}, hint::Hint)

Tests if the `ray` intersects an `object` in the optical `system`. Returns the closest intersection.

# Hint

An optional [`Hint`](@ref) can be provided to test against a specific object (and shape) in the `system` first.

!!! warning
    If a hint is provided and the object intersection is valid, the intersection will be returned immediately.
    However, it is not guaranteed that this is the true closest intersection.
"""
@inline function tracing_step!(
        system::AbstractSystem, ray::AbstractRay{R}, hint::Hint) where {R <: Real}
    # Test against hinted object
    intersection!(ray, trace_one(system, ray, hint))
    return nothing
end

@inline function tracing_step!(
        system::AbstractSystem, ray::AbstractRay{R}, ::Nothing) where {R <: Real}
    # Test against all objects in system
    intersection!(ray, trace_all(system, ray))
    return nothing
end

"""
    trace_system!(system::AbstractSystem, beam::Beam{T}; r_max::Int = 20) where {T <: Real}

Trace a [`Beam`](@ref) through an optical `system`. Maximum number of tracing steps can be capped by `r_max`.

# Tracing logic

The intersection of the last ray of the `beam` with any objects contained within the `system` is tested.
If an object is hit, the optical interaction is calculated. If no interaction occurs or no
further objects are hit, the tracing procedure is stopped.

# Arguments

- `system:`: The optical system through which the [`Beam`](@ref) is traced.
- `beam`: The [`Beam`](@ref) object to be traced.
- `r_max`: Maximum number of tracing iterations. Default is 20.
"""
function trace_system!(system::AbstractSystem, beam::Beam{T, R};
        r_max::Int = 20) where {T <: Real, R <: AbstractRay{T}}
    # Test until max. number of rays in beam reached
    interaction::Nullable{BeamInteraction{T, R}} = nothing
    while length(rays(beam)) < r_max
        ray = last(rays(beam))
        if interaction === nothing
            tracing_step!(system, ray, nothing)
        else
            tracing_step!(system, ray, hint(interaction))
        end
        # Test if intersection is valid
        ray_intersection = intersection(ray)
        if isnothing(ray_intersection)
            break
        end
        obj = object(ray_intersection)
        interaction = interact3d(
            system, obj, beam, ray)::Union{Nothing, BeamInteraction{T, R}}
        if isnothing(interaction)
            break
        end
        # Append ray to beam tail
        push!(beam, interaction)
    end
    return nothing
end

"""
    retrace_system!(system, beam)

This function tries to reuse data from a previously solved `beam` in order to solve the `system` againg using a sequential approach.

# Retracing

The retracing logic for an already solved `beam` loops over the rays and children and is as follows:

`Begin`
  1. Test if current `ray` has a valid `intersection`
      - If not, mark beam tail for cleanup and go to `End`
  2. Recalculate the `intersection`
      - If a hint was provided by a previous interaction, use hinted object
      - Else, test against previous `intersection`
  3. Test if the `ray` still has a valid `intersection` after recalculation
      - If no object is hit, mark beam tail for cleanup and go to `End`
`Interact`
  1. Recalculate the optical `interaction`
      - Catch hints provided for next `ray`
      - If no `interaction` occurs, mark beam tail for conditional cleanup and go to `End`
  2. Add the interaction to the current `beam`
      - If another `ray` follows, modify the next starting position
            - Go to `Begin`
      - Else mark children for cleanup, push new ray to `beam` tail
            - Go to `End`
`End`
  1. If cleanup is required, do conditionally
      - remove all beam tail rays after current `ray`
      - remove all beam children
      - reset beam tail ray intersection to nothing

!!! warning "Retracing blocked beam paths"
    The  implemented standard retracing procedure can handle beam path invalidations under certain conditions. However, one case that will lead to a **silent error** is if an element in the system is moved such that it **blocks the beam path between two other elements**. The retracer will not be able to detect this, since the testing of the previous intersection will return a valid intersection.

    If this kind of situation must be modeled, e.g. in the case of an optical chopper wheel, retracing should be disabled.
"""
function retrace_system!(
        system::AbstractSystem, beam::Beam{T, R}) where {T <: Real, R <: AbstractRay{T}}
    # Cleanup flags
    cleanup_children = false
    cleanup_tail = false
    reset_intersection = false
    cutoff = 0
    # Buffer variables
    _interaction::Nullable{BeamInteraction{T, R}} = nothing
    _hint::Nullable{Hint} = nothing
    for (i, ray) in enumerate(rays(beam))
        # Test if intersection is valid
        _intersection = intersection(ray)
        if isnothing(_intersection)
            cleanup_children = true
            cleanup_tail = true
            reset_intersection = true
            cutoff = i
            break
        end
        # Recalculate current intersection
        if isnothing(_hint)
            # Retrace intersection
            intersection!(ray, intersect3d(object(_intersection), ray))
        else
            # Retrace hint
            intersection!(ray, intersect3d(shape(_hint), ray))
            if !isnothing(intersection(ray))
                object!(intersection(ray), object(_hint))
            end
        end
        # Test if intersection is valid
        if isnothing(intersection(ray))
            cleanup_children = true
            cleanup_tail = true
            reset_intersection = true
            cutoff = i
            break
        end
        # Test if interaction is still valid
        _interaction = interact3d(system, object(intersection(ray)), beam, ray)
        # Catch hint
        _hint = hint(_interaction)
        if isnothing(_interaction)
            # Do not set cleanup since nothing is a valid interaction
            if length(beam.rays) > i
                cleanup_tail = true
                cutoff = i
            end
            break
        end
        if i < length(rays(beam))
            # Modify following ray (NOT THREAD-SAFE)
            replace!(beam, _interaction, i + 1)
        else
            # Valid new interaction, drop children and add new ray
            cleanup_children = true
            push!(beam, _interaction)
            break
        end
    end
    # Drop current branch since path has been altered
    cleanup_children && _drop_beams!(beam)
    # Drop all disconnected rays after last valid intersection, reset tail intersection to nothing
    cleanup_tail && deleteat!(rays(beam), (cutoff + 1):length(rays(beam)))
    reset_intersection && intersection!(last(rays(beam)), nothing)

    return nothing
end

"""
    trace_system!(system::System, gauss::GaussianBeamlet{T}; r_max::Int = 20) where {T <: Real}

Trace a [`GaussianBeamlet`](@ref) through an optical `system`. Maximum number of tracing steps can be capped by `r_max`.

# Tracing logic

The chief, waist and divergence beams are traced step-by-step through the `system`.
For each intersection after a [`tracing_step!`](@ref), the intersections are compared.
If all rays hit the same target, the optical interaction is analyzed, else the tracing stops.

# Arguments

- `system`: The optical system through which the [`GaussianBeamlet`](@ref) is traced.
- `gauss`: The [`GaussianBeamlet`](@ref) object to be traced.
- `r_max`: Maximum number of tracing iterations. Default is 20.
"""
function trace_system!(system::AbstractSystem,
        gauss::GaussianBeamlet{T};
        r_max::Int = 20) where {T <: Real}
    # Test until bundle is stopped
    interaction::Nullable{GaussianBeamletInteraction{T}} = nothing
    # Buffer variable
    seg_counter::Int = length(rays(gauss.chief))
    while seg_counter < r_max
        # Trace chief ray first
        ray = last(rays(gauss.chief))
        tracing_step!(system, ray, hint(interaction))
        isnothing(intersection(ray)) && break
        _object = object(intersection(ray))
        # Follow up with waist ray
        ray = last(rays(gauss.waist))
        tracing_step!(system, ray, hint(interaction))
        # if the waist ray is no longer hitting the same object as the chief ray stop here
        isnothing(intersection(ray)) && break
        # Follow up with divergence ray
        ray = last(rays(gauss.divergence))
        tracing_step!(system, ray, hint(interaction))
        # if the divergence ray is no longer hitting the same object as the chief ray stop here
        isnothing(intersection(ray)) && break
        # If beams do not hit same target stop tracing
        if !_beams_hits_same_shape(gauss, seg_counter)
            # Ensure that no intersection artifacts remain
            intersection!(last(rays(gauss.chief)), nothing)
            intersection!(last(rays(gauss.waist)), nothing)
            intersection!(last(rays(gauss.divergence)), nothing)
            break
        end
        # Calculate interaction
        interaction = interact3d(system,
            _object,
            gauss,
            seg_counter)
        if isnothing(interaction)
            break
        end
        # Add rays to gauss beam
        push!(gauss, interaction)
        seg_counter += 1
    end
    return nothing
end

"""
    retrace_system!(system::System, gauss::GaussianBeamlet{T}) where {T <: Real}

Retrace the beam stored in `GaussianBeamlet` through the optical `system`. Chief, waist and divergence ray intersections and interactions are recalculated.
All rays must hit the same object, or the retracing step is aborted. If retracing is stopped before the end of the beam is reached, further rays are dropped.
"""
function retrace_system!(
        system::AbstractSystem, gauss::GaussianBeamlet{T}) where {T <: Real}
    # Cleanup flags
    cleanup_children = false
    cleanup_tail = false
    reset_intersection = false
    cutoff = 0
    # Buffer variables
    _interaction::Nullable{GaussianBeamletInteraction{T}} = nothing
    _hint::Nullable{Hint} = nothing
    # Test if gauss beam is healthy
    n_c = length(rays(gauss.chief))
    n_w = length(rays(gauss.waist))
    n_d = length(rays(gauss.divergence))
    if !(n_c == n_w == n_d)
        error("Gaussian beamlet is broken")
    end
    # Iterate over chief rays
    for (i, c_ray) in enumerate(rays(gauss.chief))
        # Test if intersection is valid
        _intersection = intersection(c_ray)
        if isnothing(_intersection)
            cleanup_children = true
            cleanup_tail = true
            reset_intersection = true
            cutoff = i
            break
        end
        # Recalculate current intersection (use shape of hint if available)
        w_ray = rays(gauss.waist)[i]
        d_ray = rays(gauss.divergence)[i]
        if isnothing(_hint)
            _object = object(_intersection)
            _shape = shape(_intersection)
            intersection!(c_ray, intersect3d(_object, c_ray))
            intersection!(w_ray, intersect3d(_object, w_ray))
            intersection!(d_ray, intersect3d(_object, d_ray))
        else
            _object = object(_hint)
            _shape = BeamletOptics.shape(_hint)
            intersection!(c_ray, intersect3d(_shape, c_ray))
            intersection!(w_ray, intersect3d(_shape, w_ray))
            intersection!(d_ray, intersect3d(_shape, d_ray))
        end
        # Test if all beams still hit the same target
        if !_beams_hits_same_shape(gauss, i)
            cleanup_children = true
            cleanup_tail = true
            reset_intersection = true
            cutoff = i
            break
        end
        # Test if valid intersection
        if isnothing(intersection(c_ray))
            cleanup_children = true
            cleanup_tail = true
            reset_intersection = true
            cutoff = i
            break
        end
        # Update object field
        object!(intersection(c_ray), _object)
        object!(intersection(w_ray), _object)
        object!(intersection(d_ray), _object)
        # Test if interaction is still valid
        _interaction = interact3d(system, _object, gauss, i)
        # Catch hint
        _hint = hint(_interaction)
        if isnothing(_interaction)
            # Do not set cleanup since nothing is a valid interaction
            if n_c > i
                cleanup_tail = true
                cutoff = i
            end
            break
        end
        # Update next beamlet section
        if i < n_c
            # NOT THREAD-SAFE
            replace!(gauss, _interaction, i + 1)
        else
            # Valid new interaction, drop children and add new ray
            cleanup_children = true
            push!(gauss, _interaction)
            break
        end
    end
    # Drop current branch since path has been altered
    if cleanup_children
        _drop_beams!(gauss)
    end
    # Drop all disconnected rays after last valid intersection, reset tail intersection to nothing
    if cleanup_tail
        deleteat!(rays(gauss.chief), (cutoff + 1):length(rays(gauss.chief)))
        deleteat!(rays(gauss.waist), (cutoff + 1):length(rays(gauss.waist)))
        deleteat!(rays(gauss.divergence), (cutoff + 1):length(rays(gauss.divergence)))
    end
    if reset_intersection
        intersection!(last(rays(gauss.chief)), nothing)
        intersection!(last(rays(gauss.waist)), nothing)
        intersection!(last(rays(gauss.divergence)), nothing)
    end
    return nothing
end

"""
    trace_system!(system, agb::AstigmaticGaussianBeamlet; r_max=20)

Trace an [`AstigmaticGaussianBeamlet`](@ref) through an optical `system`.
All 9 component beams (chief + 8 parabasal) are traced in lockstep:
the chief ray is traced first for each intersection, followed by the auxiliary rays.
All rays must hit the same shape; otherwise tracing stops.
"""
function trace_system!(system::AbstractSystem,
        agb::AstigmaticGaussianBeamlet{T};
        r_max::Int = 20) where {T <: Real}
    interaction::Nullable{AstigmaticGaussianBeamletInteraction{T}} = nothing
    seg_counter::Int = length(rays(agb.c))
    aux = _aux_beams(agb)
    while seg_counter < r_max
        # Trace chief ray first
        ray = last(rays(agb.c))
        tracing_step!(system, ray, hint(interaction))
        isnothing(intersection(ray)) && break
        _object = object(intersection(ray))
        # Follow up with all auxiliary rays
        stopped = false
        for beam in aux
            ray = last(rays(beam))
            tracing_step!(system, ray, hint(interaction))
            if isnothing(intersection(ray))
                stopped = true
                break
            end
        end
        stopped && break
        # If beams do not hit same target stop tracing
        if !_beams_hits_same_shape(agb, seg_counter)
            # Ensure that no intersection artifacts remain
            intersection!(last(rays(agb.c)), nothing)
            for beam in aux
                intersection!(last(rays(beam)), nothing)
            end
            break
        end
        # Calculate interaction
        interaction = interact3d(system, _object, agb, seg_counter)
        if isnothing(interaction)
            break
        end
        # Add rays to beamlet
        push!(agb, interaction)
        seg_counter += 1
        
        # Verify that the paraxial assumptions still hold for the new segment
        if !check_optical_invariant(agb, seg_counter)
            break
        end
    end
    return nothing
end

"""
    retrace_system!(system, agb::AstigmaticGaussianBeamlet)

Retrace the beam stored in [`AstigmaticGaussianBeamlet`](@ref) through the optical
`system`. All 9 component ray intersections and interactions are recalculated.
All rays must hit the same object, or the retracing step is aborted.
"""
function retrace_system!(system::AbstractSystem,
        agb::AstigmaticGaussianBeamlet{T}) where {T <: Real}
    # Cleanup flags
    cleanup_children = false
    cleanup_tail = false
    reset_intersection = false
    cutoff = 0
    # Buffer variables
    _interaction::Nullable{AstigmaticGaussianBeamletInteraction{T}} = nothing
    _hint::Nullable{Hint} = nothing

    all_beams = _component_beams(agb)
    aux = _aux_beams(agb)

    # Test if beam is healthy
    n_c = length(rays(agb.c))
    n_lens = map(b -> length(rays(b)), all_beams)
    if !all(==(n_c), n_lens)
        error("Astigmatic Gaussian beamlet is broken")
    end

    # Iterate over chief rays
    for (i, c_ray) in enumerate(rays(agb.c))
        # Test if intersection is valid
        _intersection = intersection(c_ray)
        if isnothing(_intersection)
            cleanup_children = true
            cleanup_tail = true
            reset_intersection = true
            cutoff = i
            break
        end
        # Recalculate current intersection (use shape of hint if available)
        if isnothing(_hint)
            _object = object(_intersection)
            intersection!(c_ray, intersect3d(_object, c_ray))
            for beam in aux
                intersection!(rays(beam)[i], intersect3d(_object, rays(beam)[i]))
            end
        else
            _object = object(_hint)
            _shape = BeamletOptics.shape(_hint)
            intersection!(c_ray, intersect3d(_shape, c_ray))
            for beam in aux
                intersection!(rays(beam)[i], intersect3d(_shape, rays(beam)[i]))
            end
        end
        # Test if all beams still hit the same target
        if !_beams_hits_same_shape(agb, i)
            cleanup_children = true
            cleanup_tail = true
            reset_intersection = true
            cutoff = i
            break
        end
        # Test if valid intersection
        if isnothing(intersection(c_ray))
            cleanup_children = true
            cleanup_tail = true
            reset_intersection = true
            cutoff = i
            break
        end
        # Update object field
        object!(intersection(c_ray), _object)
        for beam in aux
            r_int = intersection(rays(beam)[i])
            !isnothing(r_int) && object!(r_int, _object)
        end
        # Test if interaction is still valid
        _interaction = interact3d(system, _object, agb, i)
        # Catch hint
        _hint = hint(_interaction)
        if isnothing(_interaction)
            # Do not set cleanup since nothing is a valid interaction
            if n_c > i
                cleanup_tail = true
                cutoff = i
            end
            break
        end
        # Update next beamlet section
        if i < n_c
            # NOT THREAD-SAFE
            replace!(agb, _interaction, i + 1)
        else
            # Valid new interaction, drop children and add new ray
            cleanup_children = true
            push!(agb, _interaction)
        end

        # Verify that the paraxial assumptions still hold for the new segment
        if !check_optical_invariant(agb, i + 1)
            cleanup_tail = true
            cutoff = i
            break
        end
        
        if i >= n_c
            break
        end
    end
    # Drop current branch since path has been altered
    if cleanup_children
        _drop_beams!(agb)
    end
    # Drop all disconnected rays after last valid intersection, reset tail intersection
    if cleanup_tail
        deleteat!(rays(agb.c), (cutoff + 1):length(rays(agb.c)))
        for beam in aux
            deleteat!(rays(beam), (cutoff + 1):length(rays(beam)))
        end
    end
    if reset_intersection
        intersection!(last(rays(agb.c)), nothing)
        for beam in aux
            intersection!(last(rays(beam)), nothing)
        end
    end
    return nothing
end

"""
    solve_system!(system::System, beam::AbstractBeam; r_max=100, retrace=true, depth_max=typemax(Int))

Manage the tracing of an `AbstractBeam` through an optical `system`. The function retraces the `beam` if possible and then proceeds to trace each leaf of the beam tree through the system.
The condition to stop ray tracing is that the last `beam` intersection is `nothing` or the beam interaction is `nothing`. Then, the system is considered to be solved.
A maximum number of rays per `beam` (`r_max`) can be specified in order to avoid infinite calculations under resonant conditions, i.e. two facing mirrors. Likewise, `depth_max` limits how many branching levels are explored when new sub-beams are generated (for example, by beamsplitters) so that the tree cannot grow without bound.

# Arguments

- `system::System`: The optical system in which the beam will be traced.
- `beam::AbstractBeam`: The beam object to be traced through the system.
- `r_max::Int=100` (optional): Maximum number of tracing iterations for each leaf. Default is 100.
- `retrace::Bool=true` (optional): Flag to indicate if the system should be retraced. Default is true.
- `depth_max::Int=typemax(Int)` (optional): Maximum number of branching levels explored from the root beam. Default allows unlimited depth.
"""
function solve_system!(
        system::AbstractSystem,
        beam::B;
        r_max::Int = 100,
        retrace::Bool = true,
        depth_max::Int = typemax(Int)
) where {B <: AbstractBeam}
    queue = Tuple{B, Int}[(beam, 1)]
    while !isempty(queue)
        # Process beams in FIFO order.
        current, depth = popfirst!(queue)
        # Optionally retrace the current beam.
        if retrace
            retrace_system!(system, current)
        end
        # Process the current leaf beam.
        solve_leaf!(system, current; r_max = r_max)
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
    for _beam in beams(bg)
        solve_system!(system, _beam; kwargs...)
    end
    return nothing
end

@inline function solve_leaf!(system::AbstractSystem, beam::AbstractBeam; r_max = 100)
    if isnothing(_last_beam_intersection(beam))
        trace_system!(system, beam, r_max = r_max)
    end
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
