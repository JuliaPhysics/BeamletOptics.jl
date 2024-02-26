"""
    System <: AbstractSystem

A container storing the optical elements of, i.e. a camera lens or lab setup.

# Fields

- `id::UUID`: system ID (uuid4)
- `objects`: vector containing the different objects that are part of the system (subtypes of [`AbstractObject`](@ref))
"""
struct System <: AbstractSystem
    id::UUID
    objects::Vector{AbstractObject}
end

System(object::AbstractObject) = System(uuid4(), [object])
System(objects::AbstractArray{<:AbstractObject}) = System(uuid4(), objects)

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

- `id::UUID`: system ID (uuid4)
- `objects`: vector containing the different objects that are part of the system (subtypes of [`AbstractObject`](@ref))
"""
struct StaticSystem{T <: Tuple} <: AbstractSystem
    id::UUID
    objects::T
end
StaticSystem(object::AbstractObject) = StaticSystem(uuid4(), (object))
StaticSystem(object::AbstractObjectGroup) = StaticSystem([object])
StaticSystem(objects::AbstractArray{<:AbstractObject}) = StaticSystem(uuid4(), tuple(collect(Leaves(objects))...))

objects(system::StaticSystem) = system.objects

"""
    object(system::System, obj_id::UUID)

Find a specific object in the `system` based on its unique `obj_id`.
"""
function object(system::System, obj_id::UUID)::Union{AbstractObject, Nothing}
    for object in objects(system)
        if id(object) === obj_id
            return object
        end
    end
    # If no match
    return nothing
end

"""
    object(system::StaticSystem, obj_id::UUID)

Find a specific object in the `system` based on its unique `obj_id`.
"""
function object(system::StaticSystem, obj_id::UUID)
    index = object_index(system, obj_id)
    isnothing(index) && return nothing

    return system.objects[index]
end

function object_index(system::StaticSystem, obj_id::UUID)
    return findfirst(object -> id(object) === obj_id, system.objects)
end

object(::AbstractSystem, ::Nothing) = nothing

function trace_system!(::AbstractSystem, beam::B; r_max = 0) where {B <: AbstractBeam}
    @warn "Tracing for $B not implemented"
    return nothing
end

function retrace_system!(::AbstractSystem, beam::B) where {B <: AbstractBeam}
    @warn "Retracing for $B not implemented"
    return nothing
end

@inline function trace_all(system::AbstractSystem, ray::AbstractRay{R}) where {R}
    intersection::Nullable{Intersection{R}} = nothing
    for object in objects(system)
        # Find shortest intersection
        temp::Nullable{Intersection{R}} = intersect3d(object, ray)
        # Ignore miss
        if isnothing(temp)
            continue
        end
        # Catch first valid intersection
        if isnothing(intersection)
            intersection = temp
            continue
        end
        # Replace current with closer intersection
        if length(temp) < length(intersection)
            intersection = temp
        end
    end
    return intersection
end

@inline function trace_one(system::AbstractSystem, ray::AbstractRay{R}, hint::UUID) where {R}
    # Trace against hinted object
    intersection = intersect3d(object(system, hint), ray)
    # If hinted object is not intersected, trace the entire system
    if isnothing(intersection)
        intersection = trace_all(system, ray)
    end
    return intersection
end

"""
    tracing_step!(system::AbstractSystem, ray::AbstractRay{R}, hint::Nullable{UUID})

Tests if the `ray` intersects an `object` in the optical `system`. Returns the closest intersection.

# Hint UUID

An optional `hint` in the form of an UUID can be provided to test against a specific object in the `system` first.

!!! warning
    If a hint is provided and the object intersection is valid, the intersection will be returned immediately.
    However, it is not guaranteed that this is the true closest intersection. 
"""
@inline function tracing_step!(system::AbstractSystem,
        ray::AbstractRay{R},
        hint::Nullable{UUID}) where {R <: Real}
    if isnothing(hint)
        # Test against all objects in system
        intersection!(ray, trace_all(system, ray))
    else
        # Test against hinted object
        intersection!(ray, trace_one(system, ray, hint))
    end

    return nothing
end

"""
    trace_system!(system::AbstractSystem, beam::Beam{T}; r_max::Int = 20) where {T <: Real}

Trace a [`Beam`](@ref) through an optical `system`. Maximum number of tracing steps can be capped by `r_max`.

# Tracing logic

The intersection of the last ray of the `beam` with any objects in the `system` is tested.
If an object is hit, the optical interaction is analyzed and tracing continues. 
Else the tracing procedure is stopped.

# Arguments

- `system:`: The optical system through which the GaussianBeamlet is traced.
- `beam`: The GaussianBeamlet object to be traced.
- `r_max`: Maximum number of tracing iterations. Default is 20.
"""
function trace_system!(system::AbstractSystem, beam::Beam{T}; r_max::Int = 20) where {T <: Real}
    # Test until max. number of rays in beam reached
    interaction::Nullable{BeamInteraction{T}} = nothing
    while length(rays(beam)) < r_max
        ray = last(rays(beam))
        tracing_step!(system, ray, hint(interaction))
        # Test if intersection is valid
        if isnothing(intersection(ray))
            break
        end
        interaction = interact3d(system, object(system, id(intersection(ray))), beam, ray)
        if isnothing(interaction)
            break
        end
        # Append ray to beam tail
        push!(beam, interaction)
    end
    return nothing
end

function retrace_system!(system::AbstractSystem, beam::Beam{T}) where {T <: Real}
    # Retrace existing beams (NOT THREAD-SAFE)
    cutoff::Nullable{Int} = nothing
    interaction::Nullable{BeamInteraction{T}} = nothing
    for (i, ray) in enumerate(rays(beam))
        # Test if intersection is valid
        isect = intersection(ray)
        if isnothing(isect)
            cutoff = i
            break
        end

        # Recalculate current intersection
        intersection!(ray, intersect3d(object(system, id(isect)), ray))
        # Test if intersection is valid
        if isnothing(intersection(ray))
            cutoff = i
            break
        end
        # Test if interaction is still valid
        interaction = interact3d(system, object(system, id(intersection(ray))), beam, ray)
        if isnothing(interaction)
            # Do not set cutoff since nothing is a valid interaction
            break
        end
        # Modify following ray (NOT THREAD-SAFE)
        if i < length(rays(beam))
            replace!(beam, interaction, i + 1)
        end
    end
    # Drop no beams / branches
    if isnothing(cutoff)
        return nothing
    end
    # Drop current branch since path has been altered
    _drop_beams!(beam)
    # Drop all disconnected rays after last valid intersection, reset tail intersection to nothing
    if cutoff < length(rays(beam))
        deleteat!(rays(beam), (cutoff + 1):length(rays(beam)))
        intersection!(last(rays(beam)), nothing)
        return nothing
    end
end

"""
    trace_system!(system::System, gauss::GaussianBeamlet{T}; r_max::Int = 20) where {T <: Real}

Trace a [`GaussianBeamlet`](@ref) through an optical `system`. Maximum number of tracing steps can be capped by `r_max`.

# Tracing logic 

The chief, waist and divergence beams are traced step-by-step through the `system`. 
For each intersection after a [`tracing_step!`](@ref), the intersections are compared.
If all rays hit the same target, the optical interaction is analyzed, else the tracing stops.

# Arguments

- `system`: The optical system through which the GaussianBeamlet is traced.
- `gauss`: The GaussianBeamlet object to be traced.
- `r_max`: Maximum number of tracing iterations. Default is 20.
"""
function trace_system!(system::AbstractSystem,
        gauss::GaussianBeamlet{T};
        r_max::Int = 20) where {T <: Real}
    # Test until bundle is stopped
    interaction::Nullable{GaussianBeamletInteraction{T}} = nothing
    while length(rays(gauss.chief)) < r_max
        # Trace chief ray first
        ray = last(rays(gauss.chief))
        tracing_step!(system, ray, hint(interaction))
        if isnothing(intersection(ray))
            break
        end
        # Follow up with waist and divergence ray
        id_c = id(intersection(ray))
        ray = last(rays(gauss.waist))
        tracing_step!(system, ray, hint(interaction))
        id_w = id(intersection(ray))
        ray = last(rays(gauss.divergence))
        tracing_step!(system, ray, hint(interaction))
        id_d = id(intersection(ray))
        # If beams do not hit same target stop tracing
        if !(id_c == id_w == id_d)
            break
        end
        interaction = interact3d(system,
            object(system, id_c),
            gauss,
            length(rays(gauss.chief)))
        if isnothing(interaction)
            break
        end
        # Add rays to gauss beam
        push!(gauss, interaction)
    end
    return nothing
end

"""
    retrace_system!(system::System, gauss::GaussianBeamlet{T}) where {T <: Real}

Retrace the beam stored in `GaussianBeamlet` through the optical `system`. Chief, waist and divergence ray intersections and interactions are recalculated.
All rays must hit the same object, or the retracing step is aborted. If retracing is stopped before the end of the beam is reached, further rays are dropped.
"""
function retrace_system!(system::AbstractSystem, gauss::GaussianBeamlet{T}) where {T <: Real}
    cutoff::Nullable{Int} = nothing
    interaction::Nullable{GaussianBeamletInteraction{T}} = nothing
    # Test if gauss beam is healthy
    n_c = length(rays(gauss.chief))
    n_w = length(rays(gauss.waist))
    n_d = length(rays(gauss.divergence))
    if !(n_c == n_w == n_d)
        error("Gaussian beamlet is broken")
    end
    # Iterate over chief rays
    for (i, c_ray) in enumerate(rays(gauss.chief))
        if isnothing(intersection(c_ray))
            cutoff = i
            break
        end
        w_ray = rays(gauss.waist)[i]
        d_ray = rays(gauss.divergence)[i]
        obj = object(system, id(intersection(c_ray)))
        intersect_c = intersect3d(obj, c_ray)
        intersect_w = intersect3d(obj, w_ray)
        intersect_d = intersect3d(obj, d_ray)
        # Test if intersections are still valid
        if any(isnothing, (intersect_c, intersect_w, intersect_d))
            cutoff = i
            break
        end
        # Test if all beams still hit the same target
        if !(id(intersect_c) == id(intersect_w) == id(intersect_d))
            cutoff = i
            break
        end
        # Set intersection
        intersection!(c_ray, intersect_c)
        intersection!(w_ray, intersect_w)
        intersection!(d_ray, intersect_d)
        # Test if interaction is still valid
        interaction = interact3d(system, object(system, id(intersect_c)), gauss, i)
        if !isnothing(interaction) && i < n_c
            replace!(gauss, interaction, i + 1) # NOT THREAD-SAFE
        end
    end
    # Drop no beams / branches
    if isnothing(cutoff)
        return nothing
    end
    # Drop current branch since path has been altered
    _drop_beams!(gauss)
    # Drop all disconnected rays after last valid intersection, reset tail intersection to nothing
    if cutoff < n_c
        deleteat!(rays(gauss.chief), (cutoff + 1):n_c)
        deleteat!(rays(gauss.waist), (cutoff + 1):n_w)
        deleteat!(rays(gauss.divergence), (cutoff + 1):n_d)
        intersection!(last(rays(gauss.chief)), nothing)
        intersection!(last(rays(gauss.waist)), nothing)
        intersection!(last(rays(gauss.divergence)), nothing)
        return nothing
    end
    return nothing
end

function trace_system!(system::AbstractSystem, agb::AstigmaticGaussianBeamlet; r_max::Int = 20)
    # placeholder solver
    trace_system!(system, agb.c, r_max=r_max)
    trace_system!(system, agb.wxp, r_max=r_max)
    trace_system!(system, agb.wxm, r_max=r_max)
    trace_system!(system, agb.wyp, r_max=r_max)
    trace_system!(system, agb.wym, r_max=r_max)
    trace_system!(system, agb.dxp, r_max=r_max)
    trace_system!(system, agb.dxm, r_max=r_max)
    trace_system!(system, agb.dyp, r_max=r_max)
    trace_system!(system, agb.dym, r_max=r_max)
end

"""
    solve_system!(system::System, beam::AbstractBeam; r_max=100, retrace=true)

Manage the tracing of an `AbstractBeam` through an optical `system`. The function retraces the `beam` if possible and then proceeds to trace each leaf of the beam tree through the system.
The condition to stop ray tracing is that the last `beam` intersection is `nothing` or the beam interaction is `nothing`. Then, the system is considered to be solved.
A maximum number of rays per `beam` (`r_max`) can be specified in order to avoid infinite calculations under resonant conditions, i.e. two facing mirrors.

# Arguments

- `system::System`: The optical system in which the beam will be traced.
- `beam::AbstractBeam`: The beam object to be traced through the system.
- `r_max::Int=20` (optional): Maximum number of tracing iterations for each leaf. Default is 100.
- `retrace::Bool=true` (optional): Flag to indicate if the system should be retraced. Default is true.
"""
function solve_system!(system::AbstractSystem, beam::AbstractBeam; r_max = 100, retrace = true)
    B = nodetype(beam)
    # Retrace system, use stateless iterator for appendability
    if retrace
        for node::B in StatelessBFS(beam)
            retrace_system!(system, node::B)
        end
    end
    # Trace starting of each leaf in beam tree
    for leaf::B in Leaves(beam)
        for beam::B in StatelessBFS(leaf)
            if isnothing(_last_beam_intersection(beam))
                trace_system!(system, beam, r_max = r_max)
            end
        end
    end
    return nothing
end

function AbstractTrees.printnode(io::IO, node::B; kw...) where {B <: AbstractObject}
    show(io, B)
end
function AbstractTrees.printnode(io::IO,
        node::B;
        kw...) where {B <: AbstractObjectGroup}
    show(io, B)
end

Base.show(::IO, ::MIME"text/plain", system::System) =
    for obj in system.objects
        print_tree(obj)
    end
