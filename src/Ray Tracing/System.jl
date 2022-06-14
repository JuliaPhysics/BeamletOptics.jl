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
    t0 = Inf
    oID::Int = 0 # object index
    fID::Int = 0 # face index of object
    temp::Int = 0 # buffer for face index
    for (i, object) in enumerate(system.objects)
        t, temp = intersect3d(object, ray)
        if t < t0
            t0 = t
            oID = i
            fID = temp
        end
    end
    return t0, oID, fID
end

function solve_system!(system::System, beam::Beam; i_max=20)
    flag = true
    iter = 0
    # find intersect
    @debug "Ray tracing routine started."
    while flag && iter <= i_max
        t0, oID, fID = trace_system(system, beam.rays[end])
        if (t0 == Inf) || (oID == 0)
            @debug "No further intersection found!"
            flag = false
        else
            @debug "Intersection found with:" oID fID
            beam.rays[end].len = t0
            flag = interact(system.objects[oID], beam, fID)
        end
        iter += 1
    end
    @debug "Ray tracing routine ended."
end