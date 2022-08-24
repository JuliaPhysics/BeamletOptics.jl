mutable struct MoellerTrumboreAlgorithm{T}
    kϵ::T
    lϵ::T
end

function MoellerTrumboreAlgorithm(kϵ, lϵ; T=Float64)
    @assert kϵ >= 0 && lϵ > 0 "kϵ must be ≥ 0 and lϵ > 0..."
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
function (f::MoellerTrumboreAlgorithm)(face, ray::Ray{T}, E1, E2, Pv, Tv, Qv) where T
    V1 = @view face[1, :]
    V2 = @view face[2, :]
    V3 = @view face[3, :]

    fast_sub3d!(E1, V2, V1)
    fast_sub3d!(E2, V3, V1)
    fast_cross3d!(Pv, ray.dir, E2)
    Det = fast_dot3d(E1, Pv)
    # Check if ray is backfacing or missing face
    if abs(Det) < f.kϵ
        # Adjust type of Inf
        return T(Inf) # typemax(T) ?
    end
    # Compute normalized u and reject if less than 0 or greater than 1
    fast_sub3d!(Tv, ray.pos, V1)
    invDet = 1 / Det
    u = fast_dot3d(Tv, Pv) * invDet
    if (u < 0) || (u > 1)
        return T(Inf)
    end
    # Compute normalized v and reject if less than 0 or greater than 1
    fast_cross3d!(Qv, Tv, E1)
    v = fast_dot3d(ray.dir, Qv) * invDet
    if (v < 0) || (u + v > 1)
        return T(Inf)
    end
    # Compute t (type def. for t to avoid Any)
    t::T = fast_dot3d(E2, Qv) * invDet
    # Return intersection only if "in front of" ray origin
    if t < f.lϵ
        return T(Inf)
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
function intersect3d(object::Mesh{M}, ray::Ray{R}) where {M, R}
    numEl = size(object.faces, 1)
    # allocate all intermediate vectors once (note that this is NOT THREAD-SAFE)
    T = promote_type(M, R)
    fID::Int = 0
    t0::T = Inf
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
    if isinf(t0)
        return NoIntersection(T)
    else
        face = @views object.vertices[object.faces[fID, :], :]
        normal = orthogonal3d(object, fID) # fast_cross3d((face[2, :] - face[1, :]), (face[3, :] - face[1, :]))
        return Intersection{T}(t0, normalize3d(T.(normal)), missing)
    end
end

"""
    intersect3d(sphere::AbstractSphere, ray::Ray)

Intersection algorithm for `sphere`ical objects. Returns a face ID of 0 if missed and 1 if hit.\\
This function currently only works for **intersections from outside of the sphere**, but not within the sphere!
Based on [this example](https://gamedev.stackexchange.com/questions/96459/fast-ray-sphere-collision-code).
"""
function intersect3d(sphere::Sphere{S}, ray::Ray{R}) where {S, R}
    T = promote_type(S, R)
    m = ray.pos - sphere.pos
    b = fast_dot3d(m, ray.dir)
    c = fast_dot3d(m, m) - sphere.radius^2
    # Exit if ray origin outside sphere (c > 0) and ray pointing away from sphere (b > 0)
    if c > 0.0 && b > 0.0
        return NoIntersection(T)
    end
    # Exit if negative discriminant (corresponds to ray missing sphere)
    discr = b^2 - c
    if discr < 0.0
        return NoIntersection(T)
    end
    # Ray intersects, compute t
    t = -b - sqrt(discr)
    # If t < 0, ray started inside sphere
    dir = ray.pos .+ t .* ray.dir
    normal = dir .- sphere.pos
    return Intersection{T}(t, normalize3d(T.(normal)), missing)
end
