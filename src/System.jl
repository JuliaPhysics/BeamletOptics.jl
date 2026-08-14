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

@inline resolve_coincident_boundary(exiting_int, entering_int) = resolve_coincident_boundary(
    exiting_int, entering_int, object(exiting_int), object(entering_int))

@inline function resolve_coincident_boundary(
        exiting_int, entering_int, ::AbstractObject, ::AbstractObject)
    exiting_int.coincident_object = object(entering_int)
    return exiting_int
end

"""
    _resolve_coincident_group(thin_primary, exiting_int, entering_int, fallback_int)

Resolves the primary intersection and assigns adjacent coincident objects for thin-interface coatings
and flush bulk optical boundaries (e.g., cemented doublet lenses).
"""
@inline function _resolve_coincident_group(
        thin_primary, exiting_int, entering_int, fallback_int)
    if thin_primary !== nothing
        thin_primary.coincident_object = exiting_int !== nothing ? object(exiting_int) :
                                         nothing
        thin_primary.coincident_object_2 = entering_int !== nothing ? object(entering_int) :
                                           nothing
        return thin_primary
    elseif exiting_int !== nothing && entering_int !== nothing
        return resolve_coincident_boundary(exiting_int, entering_int)
    elseif exiting_int !== nothing
        return exiting_int
    elseif entering_int !== nothing
        return entering_int
    else
        return fallback_int
    end
end

@inline function trace_all(system::AbstractSystem, ray::AbstractRay{R}) where {R}
    tol = get_coincident_boundary_tolerance()
    best_t = R(Inf)
    best_int = nothing
    thin_primary = nothing
    exiting_int = nothing
    entering_int = nothing

    # Find the closest intersection and collect any coincident ones at the same distance
    for obj in objects(system)
        temp = intersect3d(obj, ray)
        (temp === nothing || length(temp) <= tol) && continue
        t = length(temp)

        if t < best_t - tol
            best_t = t
            best_int = temp
            thin_primary = exiting_int = entering_int = nothing
        elseif abs(t - best_t) > tol
            continue
        end

        if is_thin_interface(object(temp))
            thin_primary = temp
        elseif dot(direction(ray), normal3d(temp)) > 0
            exiting_int = temp
        else
            entering_int = temp
        end
    end

    best_int === nothing && return nothing
    return _resolve_coincident_group(thin_primary, exiting_int, entering_int, best_int)
end

@inline function trace_one(
        system::AbstractSystem, ray::AbstractRay{R}, hint::Hint) where {R}
    tol = get_coincident_boundary_tolerance()
    # Trace against hinted shape of object first
    primary = intersect3d(
        shape(hint)::AbstractShape{R}, ray)
    if primary === nothing
        # If hinted object is not intersected, trace the entire system
        return trace_all(system, ray)
    end

    object!(primary, object(hint))
    t_best = length(primary)
    if t_best <= tol
        # Ignore self-intersection/coincident boundary re-entry
        return trace_all(system, ray)
    end

    thin_primary = nothing
    exiting_int = nothing
    entering_int = nothing

    if is_thin_interface(object(primary))
        thin_primary = primary
    elseif dot(direction(ray), normal3d(primary)) > 0
        exiting_int = primary
    else
        entering_int = primary
    end

    # Gather other intersections at the same distance to resolve cement/contact transitions
    for obj in objects(system)
        if obj === object(hint)
            continue
        end
        temp = intersect3d(obj, ray)
        if temp !== nothing && abs(length(temp) - t_best) <= tol
            if is_thin_interface(object(temp))
                thin_primary = temp
            elseif dot(direction(ray), normal3d(temp)) > 0
                exiting_int = temp
            else
                entering_int = temp
            end
        end
    end

    return _resolve_coincident_group(thin_primary, exiting_int, entering_int, primary)
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
    trace_system!(system::AbstractSystem, beam::Beam{T}; r_max = get_default_r_max()) where {T <: Real}

Trace a [`Beam`](@ref) through an optical `system`. Maximum number of tracing steps can be capped by `r_max`.

# Tracing logic

The intersection of the last ray of the `beam` with any objects contained within the `system` is tested.
If an object is hit, the optical interaction is calculated. If no interaction occurs or no
further objects are hit, the tracing procedure is stopped.

# Arguments

- `system:`: The optical system through which the [`Beam`](@ref) is traced.
- `beam`: The [`Beam`](@ref) object to be traced.
- `r_max`: Maximum number of tracing iterations.
"""
function trace_system!(
        system::AbstractSystem,
        beam::Beam{T, R};
        # kwargs
        r_max::Int = get_default_r_max(),
        kwargs...
) where {T <: Real, R <: AbstractRay{T}}
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
        system::AbstractSystem,
        beam::Beam{T, R};
        # kwargs
        kwargs...
) where {T <: Real, R <: AbstractRay{T}}
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
            intersection!(ray, trace_all(system, ray))
        else
            intersection!(ray, trace_one(system, ray, _hint))
        end
        # Test if intersection and object are valid
        int = intersection(ray)
        _obj = isnothing(int) ? nothing : object(int)
        if isnothing(_obj)
            cleanup_children = true
            cleanup_tail = true
            reset_intersection = true
            cutoff = i
            break
        end
        _interaction = interact3d(system, _obj, beam, ray)
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
    trace_system!(system::System, gauss::GaussianBeamlet{T}; r_max = get_default_r_max()) where {T <: Real}

Trace a [`GaussianBeamlet`](@ref) through an optical `system`. Maximum number of tracing steps can be capped by `r_max`.

# Tracing logic

The chief, waist and divergence beams are traced step-by-step through the `system`.
For each intersection after a [`tracing_step!`](@ref), the intersections are compared.
If all rays hit the same target, the optical interaction is analyzed, else the tracing stops.

# Arguments

- `system`: The optical system through which the [`GaussianBeamlet`](@ref) is traced.
- `gauss`: The [`GaussianBeamlet`](@ref) object to be traced.
- `r_max`: Maximum number of tracing iterations.
"""
function trace_system!(
        system::AbstractSystem,
        gauss::GaussianBeamlet{T};
        # kwargs...
        r_max::Int = get_default_r_max(),
        kwargs...
) where {T <: Real}
    # Test until bundle is stopped
    interaction::Nullable{GaussianBeamletInteraction{T}} = nothing
    # Buffer variable
    seg_counter::Int = length(rays(gauss.chief))
    beams = _component_beams(gauss)
    while seg_counter < r_max
        h = hint(interaction)
        stopped = false
        for b in beams
            r = last(rays(b))
            tracing_step!(system, r, h)
            if isnothing(intersection(r))
                stopped = true
                break
            end
        end
        if stopped || !_beams_hits_same_shape(gauss, seg_counter)
            for b in beams
                intersection!(last(rays(b)), nothing)
            end
            break
        end
        _object = object(intersection(last(rays(gauss.chief))))
        if isnothing(_object)
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
        system::AbstractSystem,
        gauss::GaussianBeamlet{T};
        # kwargs
        kwargs...
) where {T <: Real}
    # Cleanup flags
    cleanup_children = false
    cleanup_tail = false
    reset_intersection = false
    cutoff = 0
    # Buffer variables
    _interaction::Nullable{GaussianBeamletInteraction{T}} = nothing
    _hint::Nullable{Hint} = nothing
    beams = _component_beams(gauss)
    # Test if gauss beam is healthy
    n_c = length(rays(gauss.chief))
    if !all(b -> length(rays(b)) == n_c, beams)
        error("Gaussian beamlet is broken")
    end
    # Iterate over chief rays
    for (i, c_ray) in enumerate(rays(gauss.chief))
        # Recalculate current intersection
        if isnothing(_hint)
            for b in beams
                r = rays(b)[i]
                intersection!(r, trace_all(system, r))
            end
        else
            for b in beams
                r = rays(b)[i]
                intersection!(r, trace_one(system, r, _hint))
            end
        end
        # Test if all beams still hit the same target
        if !_beams_hits_same_shape(gauss, i)
            cleanup_children = true
            cleanup_tail = true
            reset_intersection = true
            cutoff = i
            break
        end
        # Test if valid intersection and object
        c_int = intersection(c_ray)
        _object = isnothing(c_int) ? nothing : object(c_int)
        if isnothing(_object)
            cleanup_children = true
            cleanup_tail = true
            reset_intersection = true
            cutoff = i
            break
        end
        for b in beams
            r_int = intersection(rays(b)[i])
            !isnothing(r_int) && object!(r_int, _object)
        end
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
    trace_system!(system, agb::AstigmaticGaussianBeamlet; r_max = get_default_r_max(), check_invariant = true, threshold = get_invariant_threshold())

Trace an [`AstigmaticGaussianBeamlet`](@ref) through an optical `system`.
All 9 component beams (chief + 8 parabasal) are traced in lockstep:
the chief ray is traced first for each intersection, followed by the auxiliary rays.
All rays must hit the same shape; otherwise tracing stops.
"""
function trace_system!(
        system::AbstractSystem,
        agb::AstigmaticGaussianBeamlet{T};
        # kwargs
        r_max::Int = get_default_r_max(),
        check_invariant::Bool = true,
        threshold::Real = get_invariant_threshold()
) where {T <: Real}
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
        if isnothing(_object)
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
        if check_invariant &&
           !check_optical_invariant(agb, seg_counter; threshold = threshold)
            break
        end
    end
    return nothing
end

"""
    retrace_system!(system, agb::AstigmaticGaussianBeamlet; check_invariant=true, threshold = get_invariant_threshold())

Retrace the beam stored in [`AstigmaticGaussianBeamlet`](@ref) through the optical
`system`. All 9 component ray intersections and interactions are recalculated.
All rays must hit the same object, or the retracing step is aborted.
"""
function retrace_system!(
        system::AbstractSystem,
        agb::AstigmaticGaussianBeamlet{T};
        # kwargs
        check_invariant::Bool = true,
        threshold::Real = get_invariant_threshold()
) where {T <: Real}
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
        # Recalculate current intersection
        if isnothing(_hint)
            intersection!(c_ray, trace_all(system, c_ray))
            for beam in aux
                intersection!(rays(beam)[i], trace_all(system, rays(beam)[i]))
            end
        else
            intersection!(c_ray, trace_one(system, c_ray, _hint))
            for beam in aux
                intersection!(rays(beam)[i], trace_one(system, rays(beam)[i], _hint))
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
        # Test if valid intersection and object
        c_int = intersection(c_ray)
        _object = isnothing(c_int) ? nothing : object(c_int)
        if isnothing(_object)
            cleanup_children = true
            cleanup_tail = true
            reset_intersection = true
            cutoff = i
            break
        end
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
        if check_invariant && !check_optical_invariant(agb, i + 1; threshold = threshold)
            cleanup_children = true
            cleanup_tail = true
            reset_intersection = true
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

## Keyword Arguments

- `r_max = get_default_r_max()`: Maximum number of tracing iterations for each leaf.
- `retrace = true`: Flag to indicate if the system should be retraced. Default is true.
- `depth_max = get_default_depth_max()`: Maximum number of branching levels explored from the root beam.
- `check_invariant = true`: enables or disables optical invariant checks where applicable
- `threshold = get_invariant_threshold()`: threshold for paraxial invariant checks
- `power_cutoff = get_default_power_cutoff()`: power threshold below which sub-beams/split paths are dropped to prevent infinite branching (default: 0.0)
"""
function solve_system!(
        system::AbstractSystem,
        beam::B;
        # kwargs
        r_max::Int = get_default_r_max(),
        retrace::Bool = true,
        depth_max::Int = get_default_depth_max(),
        check_invariant::Bool = true,
        threshold::Real = get_invariant_threshold(),
        power_cutoff::Real = get_default_power_cutoff()
) where {B <: AbstractBeam}
    queue = Tuple{B, Int}[(beam, 1)]
    head = 1
    while head <= length(queue)
        # Process beams in FIFO order without O(K) array shifts.
        current, depth = queue[head]
        head += 1
        # Check power cutoff before tracing the current beam.
        if optical_power(current) < power_cutoff
            _drop_beams!(current)
            continue
        end
        # Optionally retrace the current beam.
        if retrace
            retrace_system!(system, current; check_invariant, threshold)
        end
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
    @sync for _beam in beams(bg)
        Threads.@spawn solve_system!(system, _beam; kwargs...)
    end
    return nothing
end

@inline function solve_leaf!(
        system::AbstractSystem, beam::AbstractBeam; r_max = get_default_r_max(), kwargs...)
    if isnothing(_last_beam_intersection(beam))
        trace_system!(system, beam; r_max, kwargs...)
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
