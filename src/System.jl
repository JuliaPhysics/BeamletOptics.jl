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
    best_hit::Union{Nothing, Intersection{R}} = nothing
    best_obj = nothing
    second_hit::Union{Nothing, Intersection{R}} = nothing
    second_obj = nothing
    
    # Iterate over all optical objects/interfaces in system
    for obj in objects(system)
        temp = intersect3d(obj, ray)
        if temp === nothing
            continue
        end

        # Check for coincident hits within tolerance
        if best_hit === nothing
            best_hit = temp
            best_obj = obj
        elseif abs(length(temp) - length(best_hit)) <= tol
            second_hit = temp
            second_obj = obj
            if length(temp) < length(best_hit)
                second_hit = best_hit
                second_obj = best_obj
                best_hit = temp
                best_obj = obj
            end
        elseif length(temp) < length(best_hit) - tol
            best_hit = temp
            best_obj = obj
            second_hit = nothing
            second_obj = nothing
        end
    end
    
    if best_hit === nothing
        return nothing
    elseif second_hit === nothing
        return (best_hit, best_obj)
    else
        h1_is_exit = dot(direction(ray), normal3d(best_hit)) > 0
        exiting_obj = h1_is_exit ? best_obj : second_obj
        entering_obj = h1_is_exit ? second_obj : best_obj
        mi = MultiIntersection(best_hit; exiting = exiting_obj, entering = entering_obj)
        return (mi, best_obj)
    end
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
    best_hit = shape_intersection
    best_obj = object(hint)
    second_hit = nothing
    second_obj = nothing

    for obj in objects(system)
        obj === best_obj && continue
        temp = intersect3d(obj, ray)
        if temp !== nothing && abs(length(temp) - length(best_hit)) <= tol
            second_hit = temp
            second_obj = obj
            break
        end
    end

    if second_hit === nothing
        return (best_hit, best_obj)
    else
        h1_is_exit = dot(direction(ray), normal3d(best_hit)) > 0
        exiting_obj = h1_is_exit ? best_obj : second_obj
        entering_obj = h1_is_exit ? second_obj : best_obj
        mi = MultiIntersection(best_hit; exiting = exiting_obj, entering = entering_obj)
        return (mi, best_obj)
    end
end

"""
    tracing_step!(system::AbstractSystem, ray::AbstractRay{R}, hint::Hint)

Tests if the `ray` intersects an `object` in the optical `system`. Returns the hit object, or `nothing`.
"""
@inline function tracing_step!(
        system::AbstractSystem, ray::AbstractRay{R}, hint::Hint) where {R <: Real}
    res = trace_one(system, ray, hint)
    if res === nothing
        intersection!(ray, nothing)
        return nothing
    else
        hit, obj = res
        intersection!(ray, hit)
        return obj
    end
end

@inline function tracing_step!(
        system::AbstractSystem, ray::AbstractRay{R}, ::Nothing) where {R <: Real}
    res = trace_all(system, ray)
    if res === nothing
        intersection!(ray, nothing)
        return nothing
    else
        hit, obj = res
        intersection!(ray, hit)
        return obj
    end
end

"""
    trace_system!(system::AbstractSystem, beam::Beam{T}; r_max = get_default_r_max()) where {T <: Real}

Trace a [`Beam`](@ref) through an optical `system`. Maximum number of tracing steps can be capped by `r_max`.
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
        obj = if interaction === nothing
            tracing_step!(system, ray, nothing)
        else
            tracing_step!(system, ray, hint(interaction))
        end
        # Test if intersection is valid
        ray_intersection = intersection(ray)
        if isnothing(ray_intersection) || isnothing(obj)
            break
        end
        # Dispatch to interface interaction depending on intersection type
        interaction = (ray_intersection isa MultiIntersection ?
            interact3d(system, ray_intersection, beam, ray) :
            interact3d(system, obj, beam, ray))::Union{Nothing, BeamInteraction{T, R}}
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
        ray_c = last(rays(gauss.chief))
        obj_c = tracing_step!(system, ray_c, hint(interaction))
        # Follow up with waist ray
        ray_w = last(rays(gauss.waist))
        obj_w = tracing_step!(system, ray_w, hint(interaction))
        # Follow up with divergence ray
        ray_d = last(rays(gauss.divergence))
        obj_d = tracing_step!(system, ray_d, hint(interaction))

        int_c = intersection(ray_c)
        int_w = intersection(ray_w)
        int_d = intersection(ray_d)

        # If any beam misses or beams do not hit same target stop tracing
        if any(isnothing, (int_c, int_w, int_d)) || (obj_c !== obj_w || obj_c !== obj_d) || !_beams_hits_same_shape(gauss, seg_counter)
            # Ensure that no intersection artifacts remain
            intersection!(ray_c, nothing)
            intersection!(ray_w, nothing)
            intersection!(ray_d, nothing)
            break
        end
        # Calculate interaction
        ray_intersection = intersection(ray_c)
        interaction = if ray_intersection isa MultiIntersection
            interact3d(system, ray_intersection, gauss, seg_counter)
        else
            interact3d(system, obj_c, gauss, seg_counter)
        end
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
        ray_c = last(rays(agb.c))
        obj_c = tracing_step!(system, ray_c, hint(interaction))
        
        # Follow up with all auxiliary rays
        missed = isnothing(intersection(ray_c))
        for beam in aux
            ray_aux = last(rays(beam))
            obj_aux = tracing_step!(system, ray_aux, hint(interaction))
            if isnothing(intersection(ray_aux)) || obj_aux !== obj_c
                missed = true
            end
        end

        if missed || !_beams_hits_same_shape(agb, seg_counter)
            # Ensure that no intersection artifacts remain
            intersection!(last(rays(agb.c)), nothing)
            for beam in aux
                intersection!(last(rays(beam)), nothing)
            end
            break
        end
        # Calculate interaction
        ray_intersection = intersection(ray_c)
        interaction = if ray_intersection isa MultiIntersection
            interact3d(system, ray_intersection, agb, seg_counter)
        else
            interact3d(system, obj_c, agb, seg_counter)
        end
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
