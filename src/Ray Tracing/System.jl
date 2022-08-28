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

function trace_system(system::System, ray::Ray)
    intersection::Intersection = NoIntersection(Float64)
    for (i, object) in enumerate(system.objects)
        temp = intersect3d(object, ray)
        if temp.t < intersection.t
            intersection = temp
            intersection.oID = i
        end
    end
    return intersection
end

function solve_system!(system::System, beam::Beam; i_max=20)
    flag = true
    iter = 1
    @debug "Ray tracing routine started."
    while flag && iter <= i_max
        intersection = trace_system(system, beam.rays[end])
        if isinf(intersection.t)
            @debug "No further intersection found!"
            flag = false
        else
            @debug "Intersection found with:" intersection.oID
            beam.rays[end].intersection = intersection
            flag = interact(system.objects[intersection.oID], beam)
        end
        iter += 1
    end
    @debug "Ray tracing routine ended."
    return nothing
end
