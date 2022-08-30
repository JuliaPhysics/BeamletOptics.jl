struct System{T}
    objects::Vector{T}
    function System(objects)
        # Iteratively create union of all object types
        types = :(Union{})
        for obj in objects
            append!(types.args, [typeof(obj)])
        end
        new{eval(types)}(objects)
    end
end

function trace_system!(system::System, beam::Beam{T}; r_max::Int=20) where T
    # Test until max. number of rays in beam reached
    while length(beam.rays) < r_max
        intersection = NoIntersection(T)
        ray = beam.rays[end]
        # Test against all objects in system
        for (i, object) in enumerate(system.objects)
            temp = intersect3d(object, ray)
            if temp.t < intersection.t
                intersection = temp
                intersection.oID = i
            end
        end
        # Test if intersection is valid
        if intersection === NoIntersection(T)
            break
        end
        ray.intersection = intersection
        # Test if interaction is valid
        interaction = interact3d(system.objects[intersection.oID], ray)
        if interaction === NoInteraction(T)
            break
        end
        # Append ray to beam tail
        append!(beam.rays, [Ray(interaction.pos, interaction.dir, NoIntersection(T), interaction.information)])
    end
    return nothing
end

function retrace_system!(system::SCDI.System, beam::SCDI.Beam{T}) where T
    # Retrace existing beams (NOT THREAD-SAFE)
    ctr = 1
    numEl = length(beam.rays)
    for ray in beam.rays
        # Test if object ID is valid
        oID = ray.intersection.oID
        if ismissing(oID)
            break
        end
        # Test if intersection is still valid (modifys current intersection)
        ray.intersection = intersect3d(system.objects[oID], ray)
        if ray.intersection === NoIntersection(T)
            break
        end
        ray.intersection.oID = oID
        # Test if interaction is still valid
        interaction = interact3d(system.objects[oID], ray)        
        if interaction === NoInteraction(T)
            break
        end
        # Modify following beam (NOT THREAD-SAFE)
        if ctr <= numEl
            next_ray = beam.rays[ctr+1]
            position!(next_ray, interaction.pos)
            direction!(next_ray, interaction.dir)
            information!(next_ray, interaction.information)
        end
        ctr += 1
    end
    # Drop all disconnected rays after last valid intersection
    if ctr < numEl
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