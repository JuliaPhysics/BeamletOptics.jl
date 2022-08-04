mutable struct MoellerTrumboreAlgorithm{T}
    kϵ::T
    lϵ::T
end

function MoellerTrumboreAlgorithm(kϵ, lϵ; T=Float64)
    @assert kϵ >= 0 "kϵ must be ≥ 0 ..."
    return MoellerTrumboreAlgorithm{T}(kϵ, lϵ)
end

"""
    MoellerTrumboreAlgorithm(face::Matrix, ray::Ray)

A culling implementation of the **Möller-Trumbore algorithm** for ray-triangle-intersection.
This algorithm evaluates the possible intersection between a `ray` and a `face` that is defined by three vertices.
If no intersection occurs, `Inf` is returned. `kϵ` is the abort threshold for backfacing and non-intersecting triangles.
`lϵ` is the threshold for negative values of `t`.
This algorithm is fast due to multiple breakout conditions.
"""
function (f::MoellerTrumboreAlgorithm)(face, ray::Ray, E1, E2, Pv, Tv, Qv)
    V1 = @view face[1, :]
    V2 = @view face[2, :]
    V3 = @view face[3, :]

    fast_sub3d!(E1, V2, V1)
    fast_sub3d!(E2, V3, V1)
    fast_cross3d!(Pv, ray.dir, E2)
    Det = fast_dot3d(E1, Pv)
    # Check if ray is backfacing or missing face
    if abs(Det) < f.kϵ
        return Inf
    end
    # Compute normalized u and reject if less than 0 or greater than 1
    fast_sub3d!(Tv, ray.pos, V1)
    invDet = 1 / Det
    u = fast_dot3d(Tv, Pv) * invDet
    if (u < 0) || (u > 1)
        return Inf
    end
    # Compute normalized v and reject if less than 0 or greater than 1
    fast_cross3d!(Qv, Tv, E1)
    v = fast_dot3d(ray.dir, Qv) * invDet
    if (v < 0) || (u + v > 1)
        return Inf
    end
    # Compute t (type def. for t to avoid Any)
    t::Float64 = fast_dot3d(E2, Qv) * invDet
    # Return intersection only if "in front of" ray origin
    if t < f.lϵ
        return Inf
    end
    return t
end

const ray_triangle_intersection = MoellerTrumboreAlgorithm(1e-9, 1e-9)

"""
    intersect3d(object::Mesh, ray::Ray)

This function is a generic implementation to check if a ray intersects the object mesh.\\
If true, the **distance** `t` is returned, where the location of intersection is `ray.pos+t*ray.dir`.\\
In addition, the face index of the mesh intersection is returned.
"""
function intersect3d(object::Mesh, ray::Ray)
    numEl = size(object.faces, 1)
    t0::Float64 = Inf
    fID::Int = 0
    # allocate all intermediate vectors once (note that this is NOT THREAD-SAFE)
    T = promote_type(eltype(object.faces), eltype(ray.dir))
    E1 = Vector{T}(undef, 3)
    E2 = Vector{T}(undef, 3)
    Pv = Vector{T}(undef, 3)
    Tv = Vector{T}(undef, 3)
    Qv = Vector{T}(undef, 3)
    for i = 1:numEl
        face = @views object.vertices[object.faces[i, :], :]
        t = ray_triangle_intersection(face, ray, E1, E2, Pv, Tv, Qv)
        # Return closest intersection
        if t < t0
            t0 = t
            fID = i
        end
    end
    return t0, fID
end
