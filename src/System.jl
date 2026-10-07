"""
    System <: AbstractSystem

A container storing the optical elements of, i.e. a camera lens or lab setup.

# Fields

- `objects`: vector containing the different objects that are part of the system (subtypes of [`AbstractObject`](@ref)).
  The constructor copies the vector it is given, such that `push!` and `delete!` of the system do not change it.
"""
struct System <: AbstractSystem
    objects::Vector{AbstractObject}
    System(objects::AbstractVector) = new(collect(AbstractObject, objects))
end

System(object::AbstractObject) = System([object])

"""
    System()

Creates an empty system, to which objects are added via `push!`.
"""
System() = System(AbstractObject[])

# `Leaves` of an empty vector is the vector itself, which is no object
"""
    objects(system::System)

Exposes all objects stored within the system. By exposing the `Leaves` of the tree only, it is ensured that `AbstractObjectGroup`s are flattened into a regular vector.
An empty system exposes no objects.
"""
objects(system::System) = isempty(system.objects) ? () : Leaves(system.objects)

"""Returns `true` if `obj` is one of the `leaves`, compared by identity."""
_is_leaf_of(obj, leaves) = any(leaf -> leaf === obj, leaves)

"""
    push!(system::System, objects::AbstractObject...) -> system

Adds the `objects` (or object groups) at the top level of the `system`, such that they are traced
by the following calls of [`solve_system!`](@ref). An object that is already part of the `system`,
directly or within a group, throws an `ArgumentError`; in this case nothing is added.

Beams and beam groups that were solved before do not know the new object: solve them again with
`solve_system!(system, beam)`. The `system` must not be changed while it is being solved.
"""
function Base.push!(system::System, objs::AbstractObject...)
    leaves = AbstractObject[o for o in objects(system)]
    for obj in objs, leaf in Leaves(obj)
        _is_leaf_of(leaf, leaves) && throw(ArgumentError(
            "the $(nameof(typeof(leaf))) is already an object of the system"))
        push!(leaves, leaf)
    end
    append!(system.objects, objs)
    return system
end

"""
    pop!(system::System) -> AbstractObject

Removes the last top-level object (or object group) of the `system` and returns it. Beams that
were solved before must be solved again, see `delete!(system, object)`.
"""
Base.pop!(system::System) = pop!(system.objects)

"""
    popat!(system::System, i::Integer) -> AbstractObject

Removes the `i`-th top-level object (or object group) of the `system` and returns it. The index
counts the objects and groups as they were added, a group counts as one. An index out of bounds
throws a `BoundsError`. Beams that were solved before must be solved again, see
`delete!(system, object)`.
"""
Base.popat!(system::System, i::Integer) = popat!(system.objects, i)

"""
    delete!(system::System, object::AbstractObject) -> system

Removes the top-level `object` (or object group) from the `system`, compared by identity. Nothing
happens for an `object` that is not part of the `system`. An object within a group can not be
removed on its own and throws an `ArgumentError`: remove the group instead.

Beams and beam groups that were solved before still end on the removed object: solve them again
with `solve_system!(system, beam)`. The `system` must not be changed while it is being solved.
"""
function Base.delete!(system::System, obj::AbstractObject)
    i = findfirst(o -> o === obj, system.objects)
    if isnothing(i)
        for top in system.objects
            any(o -> o === obj, PreOrderDFS(top)) && throw(ArgumentError(
                "the $(nameof(typeof(obj))) is part of a $(nameof(typeof(top))) of the system, remove the group instead"))
        end
        return system
    end
    deleteat!(system.objects, i)
    return system
end

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

function trace_system!(::AbstractSystem, beam::B; r_max = 0, kwargs...) where {B <: AbstractBeam}
    @warn "Tracing for $B not implemented"
    return nothing
end

# The objects and object groups at the top level of the system. A group tests its objects itself, see
# `intersect3d(::MultiShape, group, ray)`, hence the tracing does not use the flattened `objects(system)`.
_top_level(system::AbstractSystem) = objects(system)
_top_level(system::System) = system.objects

# Find the shortest intersection among all objects of the system
@inline trace_all(system::AbstractSystem, ray::AbstractRay) = intersect3d(_top_level(system), ray)

@inline function trace_one(
        system::AbstractSystem, ray::AbstractRay{R}, hint::Hint) where {R}
    # Trace against hinted shape of object
    _shape = shape(hint)::AbstractShape{R}
    intersection::Nullable{Intersection{R}} = intersect3d(
        bounding_sphere_of(current_bounding_spheres(), _shape), _shape, ray)
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
    solve_system!(system::AbstractSystem, beam::AbstractBeam; r_max=get_default_r_max(), depth_max=get_default_depth_max(), check_invariant=true, threshold=get_invariant_threshold())

Trace an `AbstractBeam` and all sub-beams it generates through an optical `system`. Every call solves the `beam` from its start: the `beam` is first reset to its untraced state with `empty!`, so nothing of an earlier solve is reused. A call after objects were moved, added or removed therefore gives the same result as solving a newly constructed beam.
The condition to stop ray tracing is that the last `beam` intersection is `nothing` or the beam interaction is `nothing`. Then, the system is considered to be solved.
A maximum number of rays per `beam` (`r_max`) can be specified in order to avoid infinite calculations under resonant conditions, i.e. two facing mirrors. Likewise, `depth_max` limits how many branching levels are explored when new sub-beams are generated (for example, by beamsplitters) so that the tree cannot grow without bound. Sub-beams beyond the depth limit are dropped from the tree.

Between two calls only the `beam` itself and its first ray keep their identity. The child beams and all later rays are new objects after every call: read `children(beam)` again instead of keeping a child.
A [`Detector`](@ref) is not reset by this function. Call `empty!(detector)` before solving again, otherwise the hits of both solves add up.

# Arguments

- `system::AbstractSystem`: The optical system in which the beam will be traced.
- `beam::AbstractBeam`: The beam object to be traced through the system.

## Keyword Arguments

- `r_max = get_default_r_max()`: Maximum number of tracing iterations for each leaf.
- `depth_max = get_default_depth_max()`: Maximum number of branching levels explored from the root beam
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
    # The bounding spheres of the objects are computed once, unless the solve of a beam group did so
    with_bounding_spheres(system) do
        _solve_beam!(system, beam; r_max, depth_max, check_invariant, threshold)
    end
    return nothing
end

function _solve_beam!(system::AbstractSystem, beam::B; r_max, depth_max, check_invariant,
        threshold) where {B <: AbstractBeam}
    queue = Tuple{B, Int}[(beam, 1)]
    while !isempty(queue)
        # Process beams in FIFO order.
        current, depth = popfirst!(queue)
        # Nothing of an earlier solve survives, also for a child that a component attached again.
        empty!(current)
        trace_system!(system, current; r_max, check_invariant, threshold)
        # Check if the maximum branching depth has been reached.
        if depth <= depth_max
            # Enqueue all child beams for subsequent processing.
            for child in children(current) # 'children' returns an iterable of sub-beams.
                push!(queue, (child, depth + 1))
            end
        else
            # Maximum braching depth is reached. Remove childrens of the current beam because they will not be solved.
            _drop_beams!(current)
        end
    end
    return nothing
end

"""
    solve_system!(system::AbstractSystem, bg::AbstractBeamGroup; progress=true, kwargs...)

Trace every beam of the beam group `bg` through the `system`, multithreaded over the member
beams. All other `kwargs` are passed on to [`solve_system!`](@ref) for each beam. Every beam
of the group is reset with `empty!` and solved from its start. The reset happens before the
first beam is traced, such that a cancelled or failed solve leaves no beam with an old path.

## Keyword Arguments

- `progress = true`: show a progress bar once tracing has run for `get_progress_threshold()`
  seconds (default 5 s). It is only drawn if `stderr` is a terminal, so documentation builds,
  CI logs and piped output stay clean.
"""
function solve_system!(
        system::AbstractSystem, bg::AbstractBeamGroup; progress::Bool = true, kwargs...)
    empty!(bg)
    # One table of bounding spheres for all beams. It belongs to a task, hence each task of the
    # threads sets it for its beams.
    with_bounding_spheres(system) do
        table = current_bounding_spheres()::BoundingSphereTable
        _with_progress(progress, length(bg), "Tracing beams: ") do prog
            Threads.@threads for _beam in beams(bg)
                with_bounding_spheres(table) do
                    solve_system!(system, _beam; kwargs...)
                end
                _tick!(prog)
            end
        end
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
