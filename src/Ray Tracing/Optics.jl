@Geometry struct Mirror end

function interact(element::Mirror, beam::Beam, fID)
    normal = orthogonal3d(element, fID)
    oray = beam.rays[end]
    npos = oray.pos + oray.len*oray.dir
    append!(beam.rays, [Ray(npos, Float64.(reflection(oray.dir, normal)))])
    return true
end

function reflection(dir, normal)
    return dir - 2*dot(dir,normal) * normal
end

@Geometry struct Lens
    ref_index::Float64
    function Lens(ref_index)
        @assert ref_index >= 1 "It is assumed that n>=1!"
        new(Float64(ref_index))
    end
end 

function interact(element::Lens, beam::Beam, fID)
    @debug "Placeholder function for Lens interaction."
    return false
end