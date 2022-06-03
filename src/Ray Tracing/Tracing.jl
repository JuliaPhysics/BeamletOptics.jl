"""
    moeller_trumbore_algorithm(face::Matrix, ray::Ray; kϵ=1e-5)

An implementation of the (Möller-Trumbore algorithm)[https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection].
This algorithm evaluates the possible intersection between a ray and a face that is defined by three vertices. If no intersection occurs, Inf is returned.
kϵ is the abort threshold for backfacing and non-intersecting triangles. 
This algorithm is fast due to multiple breakout conditions.
"""
function moeller_trumbore_algorithm(face::Matrix, ray::Ray; kϵ=1e-6)
    V1 = face[1,:]
    V2 = face[2,:]
    V3 = face[3,:]

    E1 = V2-V1
    E2 = V3-V1
    Pv = cross(ray.dir, E2)
    Det = dot(E1, Pv)
    invDet = 1/Det
    # Check if ray is backfacing or missing face
    if abs(Det) < kϵ # (Det < kϵ) && (abs(Det) < kϵ)?
        return Inf
    end
    # Compute normalized u and reject if less than 0 or greater than 1
    Tv = ray.pos-V1
    u = dot(Tv, Pv) * invDet
    if (u < 0) || (u > 1)
        return Inf
    end
    # Compute normalized v and reject if less than 0 or greater than 1
    Qv = cross(Tv, E1)
    v = dot(ray.dir, Qv) * invDet
    if (v < 0) || (u+v > 1)
        return Inf
    end
    # Compute t (type def. for t to avoid Any)
    t::Float64 = dot(E2, Qv) * invDet
    # Return intersection only if "in front of" ray origin
    if t < kϵ
        return Inf
    end
    return t
end

"""
    intersect3d(object::Geometry, ray::Ray)

This function is a generic implementation to check if a ray intersects the object geometry.
If true, the **distance** `t` is returned, where the location of intersection is `ray.pos+t*ray.dir`.
In addition, the face index of the mesh intersection is returned.
"""
function intersect3d(object::Geometry, ray::Ray)
    numEl = size(object.faces)[1]
    t0 = Inf
    fID::Int = 0
    for i = 1:numEl
        face = object.vertices[object.faces[i,:],:]
        t = moeller_trumbore_algorithm(face, ray)
        # Return closest intersection
        if t < t0
            t0 = t
        end
    end
    return t0, fID
end
