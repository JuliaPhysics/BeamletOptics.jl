struct System
    object::Vector
end

function object(system::System, obj_id::UUID)
    for object in system.object
        if id(object) == obj_id
            return object
        end
    end
    # If no match
    return error("Object ID not in system")
end

function trace_all(system::System, ray::Ray{R}) where {R}
    intersection::Nullable{Intersection{R}} = nothing
    for object in system.object
        temp = intersect3d(object, ray)
        # Find shortest intersection
        if isnothing(temp)
            continue
        elseif isnothing(intersection) || length(temp) < length(intersection)
            intersection = temp
        end
    end
    return intersection
end

function trace_one(system::System, ray::Ray{R}, hint::UUID) where {R}
    # Trace against hinted object
    intersection = intersect3d(object(system, hint), ray)
    # If hinted object is not intersected, trace the entire system
    if isnothing(intersection)
        intersection = trace_all(system, ray)
    end
    return intersection
end

function trace_system!(system::System, beam::Beam{B}; r_max::Int=20) where {B}
    interaction::Nullable{Interaction{B}} = nothing
    # Test until max. number of rays in beam reached
    while length(beam.rays) < r_max
        ray = beam.rays[end]
        if isnothing(interaction) || isnothing(hint(interaction))
            # Test against all objects in system
            intersection!(ray, trace_all(system, ray))
        else
            # Test against hinted object
            intersection!(ray, trace_one(system, ray, hint(interaction)))
        end
        # Test if intersection is valid
        if isnothing(intersection(ray))
            break
        end
        interaction = interact3d(object(system, intersection(ray).id), ray)
        if isnothing(interaction)
            break
        end
        # Append ray to beam tail
        append!(beam.rays, [Ray(uuid4(), interaction.pos, interaction.dir, nothing, interaction.parameters)])
    end
    return nothing
end

function retrace_system!(system::System, beam::Beam{B}) where {B}
    # Retrace existing beams (NOT THREAD-SAFE)
    ctr = 1
    numEl = length(beam.rays)
    for ray in beam.rays
        # Test if intersection is valid
        if isnothing(intersection(ray))
            break
        end
        # Recalculate current intersection)
        intersection!(ray, intersect3d(object(system, intersection(ray).id), ray))
        # Test if intersection is valid
        if isnothing(intersection(ray))
            break
        end
        # Test if interaction is still valid
        interaction = interact3d(object(system, intersection(ray).id), ray)
        if isnothing(interaction)
            break
        end
        # Modify following beam (NOT THREAD-SAFE)
        if ctr < numEl
            next_ray = beam.rays[ctr+1]
            position!(next_ray, position(interaction))
            direction!(next_ray, direction(interaction))
            parameters!(next_ray, parameters(interaction))
        end
        ctr += 1
    end
    # Drop all disconnected rays after last valid intersection
    if ctr < numEl+1
        deleteat!(beam.rays, ctr+1:numEl)
    end
    return nothing
end

function solve_system!(system::System, beam::Beam; r_max=20)
    # Test if anything can be retraced
    retrace_system!(system, beam)
    # Trace freely
    trace_system!(system, beam, r_max=r_max)
    return nothing
end

# there is a smarter way to do this, i.e. trace the waist and div. ray based on chief ray path
function solve_system!(system::System, gauss::GaussianBeamlet; r_max=20)
    # Test if anything can be retraced
    retrace_system!(system, gauss.chief)
    retrace_system!(system, gauss.waist)
    retrace_system!(system, gauss.divergence)
    # Trace freely
    trace_system!(system, gauss.chief, r_max=r_max)
    trace_system!(system, gauss.waist, r_max=r_max)
    trace_system!(system, gauss.divergence, r_max=r_max)
    return nothing
end
