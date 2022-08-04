struct Mirror{T} <: AbstractMesh
    mesh::Mesh{T}
end

function interact(mirror::Mirror, beam::Beam, fID)
    normal = orthogonal3d(mirror, fID)
    oray = beam.rays[end]
    npos = oray.pos + oray.len * oray.dir
    append!(beam.rays, [Ray(npos, Float64.(reflection(oray.dir, normal)))])
    return true
end

function reflection(dir, normal)
    return dir - 2 * dot(dir, normal) * normal
end

struct Lens{T} <: AbstractMesh
    mesh::Mesh{T}
    ref_index::Function
end

function interact(lens::Lens, beam::Beam, fID)
    # Check dir. of ray and surface normal
    normal = orthogonal3d(lens, fID)
    if dot(beam.rays[end].dir, normal) < 0
        @debug "Outside lens"
        n1 = 1.0
        n2 = lens.ref_index(beam.λ)
    else
        @debug "Inside lens"
        n1 = lens.ref_index(beam.λ)
        n2 = 1.0
        normal *= -1
    end
    # Calculate new dir. and pos.
    ndir = refraction(beam.rays[end].dir, normal, n1, n2)
    npos = beam.rays[end].pos + beam.rays[end].len * beam.rays[end].dir
    append!(beam.rays, [Ray(npos, ndir)])
    return true
end

function refraction(dir, normal, n1, n2)
    # dir and normal must have unit length!
    n = n1 / n2
    cosθi = -dot(normal, dir)
    sinθt² = n^2 * (1 - cosθi^2)
    # Check for total reflection
    if sinθt² > 1.0
        @debug "Total reflection"
        return reflection(dir, normal)
    end
    cosθt = sqrt(1 - sinθt²)
    return n * dir + (n * cosθi - cosθt) * normal
end
