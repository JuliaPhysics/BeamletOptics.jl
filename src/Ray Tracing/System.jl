struct System <: AbstractSystem
    id::UUID
    objects::Vector{AbstractObject}
end

System(object::AbstractObject) = System(uuid4(), [object])
System(objects::AbstractArray{<:AbstractObject}) = System(uuid4(), objects)

"""
    objects(system::System)

Exposes all objects stored within the system. By exposing the [`AbstractTrees.Leaves`](@ref) only, it is ensured that `AbstractObjectGroup`s are flattened into a regular vector.
"""
objects(system::System) = Leaves(system.objects)

"""
    object(system::System, obj_id::UUID)

Find a specific object in the `system` based on its unique `obj_id`.
"""
function object(system::System, obj_id::UUID)
    for object in objects(system)
        if id(object) == obj_id
            return object
        end
    end
    # If no match
    return error("Object ID not in system")
end

function trace_system!(::System, beam::B; r_max=0) where B<:AbstractBeam
    @warn "Tracing for $B not implemented"
    return nothing
end

function retrace_system!(::System, beam::B) where B<:AbstractBeam
    @warn "Retracing for $B not implemented"
    return nothing
end

@inline function trace_all(system::System, ray::Ray{R}) where {R}
    intersection::Nullable{Intersection{R}} = nothing
    for object in objects(system)
        # Find shortest intersection
        temp = intersect3d(object, ray)
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

@inline function trace_one(system::System, ray::Ray{R}, hint::UUID) where {R}
    # Trace against hinted object
    intersection = intersect3d(object(system, hint), ray)
    # If hinted object is not intersected, trace the entire system
    if isnothing(intersection)
        intersection = trace_all(system, ray)
    end
    return intersection
end

@inline function tracing_step!(system::System, ray::Ray{R}, hint::Nullable{UUID}) where R<:Real
    if isnothing(hint)
        # Test against all objects in system
        intersection!(ray, trace_all(system, ray))
    else
        # Test against hinted object
        intersection!(ray, trace_one(system, ray, hint))
    end
end

function trace_system!(system::System, beam::Beam{T}; r_max::Int=20) where T<:Real
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

function retrace_system!(system::System, beam::Beam{T}) where T<:Real
    # Retrace existing beams (NOT THREAD-SAFE)
    cutoff::Nullable{Int} = nothing
    interaction::Nullable{BeamInteraction{T}} = nothing
    for (i, ray) in enumerate(rays(beam))
        # Test if intersection is valid
        if isnothing(intersection(ray))
            cutoff = i
            break
        end
        # Recalculate current intersection
        intersection!(ray, intersect3d(object(system, intersection(ray).id), ray))
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
            replace!(beam, interaction, i+1)
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
        deleteat!(rays(beam), cutoff+1:length(rays(beam)))
        intersection!(last(rays(beam)), nothing)
        return nothing
    end
end

function trace_system!(system::System, gauss::GaussianBeamlet{T}; r_max::Int=20) where T<:Real
    # Test until bundle is stopped
    interaction::Nullable{GaussianBeamletInteraction{T}} = nothing
    while length(rays(gauss.chief)) < r_max
        # Trace chief ray first
        ray = last(rays(gauss.chief))
        tracing_step!(system, ray, hint(interaction))
        if isnothing(intersection(ray))
            break
        end
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
        interaction = interact3d(system, object(system, id_c), gauss, length(rays(gauss.chief)))
        if isnothing(interaction)
            break
        end
        # Add rays to gauss beam
        push!(gauss, interaction)
    end
    return nothing
end

function retrace_system!(system::System, gauss::GaussianBeamlet{T}) where T<:Real
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
        if any(isnothing.((intersect_c, intersect_w, intersect_d)))
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
        if i < n_c
            replace!(gauss, interaction, i+1) # NOT THREAD-SAFE
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
        deleteat!(rays(gauss.chief), cutoff+1:n_c)
        deleteat!(rays(gauss.waist), cutoff+1:n_w)
        deleteat!(rays(gauss.divergence), cutoff+1:n_d)
        intersection!(last(rays(gauss.chief)), nothing)
        intersection!(last(rays(gauss.waist)), nothing)
        intersection!(last(rays(gauss.divergence)), nothing)
        return nothing
    end
    return nothing
end

function solve_system!(system::System, beam::AbstractBeam; r_max=20, retrace=true)
    # Retrace system, use stateless iterator for appendability
    if retrace
        for node in StatelessBFS(beam)
            retrace_system!(system, node)
        end
    end
    # Trace starting of each leaf in beam tree
    for leaf in Leaves(beam)
        for beam in StatelessBFS(leaf)
            if isnothing(_last_beam_intersection(beam))
                trace_system!(system, beam, r_max=r_max)
            end
        end
    end
    return nothing
end

AbstractTrees.printnode(io::IO, node::B; kw...) where B<:SCDI.AbstractObject = show(io, B)
AbstractTrees.printnode(io::IO, node::B; kw...) where B<:SCDI.AbstractObjectGroup = show(io, B)

Base.show(::IO, ::MIME"text/plain", system::System) = for obj in system.objects; print_tree(obj); end