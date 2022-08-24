interact(entity::AbstractEntity, beam::Beam) = false, missing, missing
interact!(entity::AbstractEntity, beam::Beam) = false

struct Mirror{T} <: AbstractMesh
    mesh::Mesh{T}
end

function interact(mirror::Mirror, beam::Beam)
    normal = beam.rays[end].intersection.n
    oray = beam.rays[end]
    npos = oray.pos + oray.intersection.t * oray.dir
    append!(beam.rays, [Ray(npos, Float64.(reflection3d(oray.dir, normal)))])
    return true
end

struct Prism{T} <: AbstractMesh
    mesh::Mesh{T}
    ref_index::Function
end

function interact(prism::Prism, beam::Beam)
    # Check dir. of ray and surface normal
    normal = beam.rays[end].intersection.n
    if dot(beam.rays[end].dir, normal) < 0
        @debug "Outside prism"
        n1 = 1.0
        n2 = prism.ref_index(beam.λ)
    else
        @debug "Inside prism"
        n1 = prism.ref_index(beam.λ)
        n2 = 1.0
        normal *= -1
    end
    # Calculate new dir. and pos.
    ndir = refraction3d(beam.rays[end].dir, normal, n1, n2)
    npos = beam.rays[end].pos + beam.rays[end].intersection.t * beam.rays[end].dir
    append!(beam.rays, [Ray(npos, ndir)])
    return true
end

struct BallLens{T} <: AbstractSphere
    sphere::Sphere{T}
    ref_index::Function
end

"""
    interact(lens::BallLens, beam::Beam, fID)

Placeholder interaction scheme for a refractive sphere.
"""
function interact(blens::BallLens, beam::Beam)
    dir = beam.rays[end].dir
    npos = beam.rays[end].pos + beam.rays[end].intersection.t * dir
    # Normal equation for a sphere at point p: n = |p-s|
    normal = beam.rays[end].intersection.n
    if dot(dir, normal) < 0
        @debug "Outside lens"
        n1 = 1.0
        n2 = blens.ref_index(beam.λ)
    else
        @debug "Inside lens"
        n1 = 1.0
        n2 = blens.ref_index(beam.λ)
        normal *= -1
    end
    ndir = refraction3d(dir, normal, n1, n2)
    append!(beam.rays, [Ray(npos, ndir)])
    return true
end