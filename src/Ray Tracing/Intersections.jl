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
function (f::MoellerTrumboreAlgorithm)(face, ray::AbstractRay{T}, E1, E2, Pv, Tv, Qv) where T
    V1 = @view face[1, :]
    V2 = @view face[2, :]
    V3 = @view face[3, :]

    fast_sub3d!(E1, V2, V1)
    fast_sub3d!(E2, V3, V1)
    fast_cross3d!(Pv, direction(ray), E2)
    Det = fast_dot3d(E1, Pv)
    # Check if ray is backfacing or missing face
    if abs(Det) < f.kϵ
        # Adjust type of Inf
        return T(Inf) # typemax(T) ?
    end
    # Compute normalized u and reject if less than 0 or greater than 1
    fast_sub3d!(Tv, position(ray), V1)
    invDet = 1 / Det
    u = fast_dot3d(Tv, Pv) * invDet
    if (u < 0) || (u > 1)
        return T(Inf)
    end
    # Compute normalized v and reject if less than 0 or greater than 1
    fast_cross3d!(Qv, Tv, E1)
    v = fast_dot3d(direction(ray), Qv) * invDet
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
If true, the **distance** `t` is returned, where the location of intersection is `position(ray)+t*direction(ray)`.\\
In addition, the face index of the mesh intersection is returned.
"""
function intersect3d(object::AbstractMesh{M}, ray::AbstractRay{R}) where {M, R}
    numEl = size(faces(object), 1)
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
        face = @views vertices(object)[faces(object)[i, :], :]
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
        face = @views vertices(object)[faces(object)[fID, :], :]
        normal = orthogonal3d(object, fID)
        return Intersection{T}(t0, normalize3d(T.(normal)), missing)
    end
end

"""
    intersect3d(sphere::AbstractSphere, ray::Ray)

Intersection algorithm for `sphere`ical objects.
"""
function intersect3d(sphere::Sphere{S}, ray::Ray{R}) where {S, R}
    T = promote_type(S, R)
    intersections = ray_sphere_intersections(sphere, ray)

    # find the physically more relevant intersection (i.e. the one which is in direction of ray)
    if all(i -> i === NoIntersection(T), intersections) || all(i -> i.t < eps(T), intersections)
        # No intersections or both intersections are behind the ray
        return NoIntersection(T)
    elseif intersections[1].t < eps(T)
        return intersections[2]
    elseif intersections[2].t < eps(T)
        return intersections[1]
    else
        return argmin(i -> i.t, intersections)
    end
end

"""
    ray_sphere_intersections(sphere::Sphere{S}, ray::Ray{R})

General intersection function which returns **ALL** intersections of the ray with the sphere.
This is useful for consecutive methods in order to decide which of the found intersections
is the relevant one.
"""
function ray_sphere_intersections(sphere::Sphere{S}, ray::Ray{R}) where {S, R}
    T = promote_type(S, R)
    m = position(ray) - position(sphere)
    b = fast_dot3d(m, direction(ray))
    c = fast_dot3d(m, m) - sphere.radius^2
    # Exit if ray origin outside sphere (c > 0) and ray pointing away from sphere (b > 0)
    if c > 0.0 && b > 0.0
        return NoIntersection(T), NoIntersection(T)
    end
    # Exit if negative discriminant (corresponds to ray missing sphere)
    discr = b^2 - c
    if discr < 0.0
        return NoIntersection(T), NoIntersection(T)
    end

    if discr ≈ zero(T)
        # tangent case, only one intersection. Is this one relevant optically?
        t = -b
        dir = position(ray) .+ t .* direction(ray)
        normal = (dir.- position(sphere)) * sign(sphere.radius)
        return Intersection{T}(t, normalize3d(T.(normal)), missing), NoIntersection(T)
    else
        # Ray intersects twice, compute both t's
        sqrt_discr = sqrt(discr)
        t1, t2 = -b - sqrt_discr, -b + sqrt_discr

        dir1 = position(ray) .+ t1 .* direction(ray)
        normal1 = (dir1 .- position(sphere)) * sign(sphere.radius)

        dir2 = position(ray) .+ t2 .* direction(ray)
        normal2 = (dir2 .- position(sphere)) * sign(sphere.radius)
        return Intersection{T}(t1, normalize3d(T.(normal1)), missing), Intersection{T}(t2, normalize3d(T.(normal2)), missing)
    end
end

"""
    intersect3d(sphere::AbstractSphere, ray::Ray)

Intersection algorithm for `lens` objects.
"""
function intersect3d(lens::SingletLens, ray::Ray{T}) where T
    intersections = intersect3d(lens, lens.front, ray), intersect3d(lens, lens.back, ray)

    return argmin(i -> i.t, intersections)
end

"""
    intersect3d(lens::Lens, s::ConicSurface, ray::Ray)

Intersection implementation for a lens with conic surfaces. This intersects the surface as
equivalent spheres. This has to be implemented by more specific conic surfaces to handle more
complex cases (e.g. aspheres or cylindric lenses).
"""
function intersect3d(lens::SingletLens, s::ConicSurface, ray::Ray{T}) where T
    # get test sphere
    sph = sphere(lens, s)

    # get sphere intersections
    intersections = ray_sphere_intersections(sph, ray)
    if all(i -> i === NoIntersection(T), intersections) || all(i -> i.t < eps(T), intersections)
        # No intersections or both intersections are behind the ray
        return NoIntersection(T)
    else
        # Here the code deviates between the full sphere and a real lens. We are no longer
        # looking for the shortest interaction path but the path which hits an actual physical part of the lens

        # get the position of the intersections
        pos = [position(ray) .+ i.t*direction(ray) for i in intersections]
        _sign = (s === lens.front ? -1 : 1)
        n = _sign*normal(lens)
        bpos = position(lens) .- lens.edge_thickness/2 .* n
        distances = [norm3d(p_i - bpos) for p_i in pos]

        # FIXME: This needs a handling of the chip zone! If the ray hits within mechanical_semi_diameter
        # but within the chip zone a plane intersection should be calculated instead. This
        # should be easy to fix...
        if all(d -> d < mechanical_semi_diameter(s), distances)
            # return the shorter one
            if first(intersections).t < last(intersections).t && first(intersections).t > eps(T)
                return first(intersections)
            else
                return last(intersections)
            end
        else
            for (i, d) in enumerate(distances)
                if d < mechanical_semi_diameter(s) && intersections[i].t > eps(T)
                    return intersections[i]
                end
            end
        end
    end

    # if this is reached, no intersection falls within the mechanical aperture
    return NoIntersection(T)
end

"""
    intersect3d(lens::Lens, s::AsphericSurface, ray::Ray)

Intersection implementation for a lens with an aspheric surface.
"""
function intersect3d(lens::SingletLens, s::AsphericSurface, ray::Ray{T}, iter_max=20, ϵ=1e-9) where T
    @warn "Unfinished! Aspheric surface seems to give incorrect normals at the moment!"

    # test intersection with the lens plane
    intersection = intersect3d(lens, PlanarSurface(mechanical_semi_diameter(s)), ray, s === lens.front ? -1 : 1)
    # any intersection at all!?
    intersection === NoIntersection(T) && return intersection

    # now get the intersection with the sphere as first order approximation
    intersection = intersect3d(lens, s.sphere, ray)

    # now we need to use an interative approach to converge onto the asphere surface.
    # First step: Transform the intersection ray/point to the lens coordinate system
    B_asph = base_transform(orientation(lens))
    r0 = B_asph * (position(ray) .+ intersection.t * direction(ray))
    lpos = B_asph * position(lens)
    r0d = direction(ray)
    c1 = 1/radius_of_curvature(s)
    i = 0
    r_iter, r_last = copy(r0), copy(r0)
    asph_normal = zeros(T, 3)

    # iterate an approximation of the intersection until the change is smaller than `ϵ`
    # or the max number of iterations is reached
    while i < iter_max && (i == 0 || norm3d(r_iter - r_last) > ϵ)
        i += 1

        # store the z-coordinate
        z0 = r_iter[3] - lpos[3]
        r = √((r_iter[1] - lpos[1])^2 + (r_iter[2] - lpos[2])^2)

        # replace z component by sag value at the current point
        z1star = (B_asph*sag_coordinate(lens, s, r))[3]
        # calculate the gradient
        asph_normal[3] = √(1-c1*r*(1+conic_constant(s)))
        # FIXME: Ignores odd coefficients
        aux = (c1+asph_normal[3]+sum(i->s.A_even[i-1]*r^(2i-2), 2:(length(s.A_even)+1)))
        asph_normal[1] = - r_iter[1]*aux
        asph_normal[2] = - r_iter[2]*aux

        # get the translation factor G for the new estimate
        G = asph_normal[3]*(z1star - z0) / (fast_dot3d(r0d, asph_normal))

        # generate a new approximate intersection point
        r_last .= r_iter
        r_iter .= G*r0d .+ r_last
    end

    # check the intersection sign
    if sign(fast_dot3d(r_iter - position(ray), direction(ray))) == -1
        return NoIntersection(T)
    end
    # check if the ray is within the numerical error margin
    if norm3d(r_iter - B_asph*position(ray)) < ϵ
        return NoIntersection(T)
    end
    # check if the ray hits the mechanical aperture
    _sign = (s === lens.front ? -1 : 1)*sign(radius_of_curvature(s))
    edge_pos = B_asph * (position(lens) +_sign*lens.edge_thickness/2*normal(lens))
    if (norm3d(r_iter - edge_pos) > mechanical_semi_diameter(s))
        return NoIntersection(T)
    end

    # prepare normal vector
    r = √((r_iter[1] - lpos[1])^2 + (r_iter[2] - lpos[2])^2)
    asph_normal[3] = √(1-c1*r*(1+conic_constant(s)))
    aux = (c1+asph_normal[3]+mapreduce(i->s.A_even[i-1]*r^(2i-2), +, 2:(length(s.A_even)+1)))
    asph_normal[1] = - r_iter[1]*aux
    asph_normal[2] = - r_iter[2]*aux
    normalize3d!(asph_normal)
    asph_normal .*= -sign(radius_of_curvature(s))

    return Intersection{T}(norm3d(r_iter - position(ray)), asph_normal, missing)
end

function sag_coordinate(l::SingletLens, s::AsphericSurface, r::Real)
    _sign = (s === l.front ? -1 : 1)*sign(radius_of_curvature(s))
    position(l) + _sign*(l.edge_thickness/2 + sag(s, clear_semi_diameter(s) - r))*normal(l)
end

"""
    intersect3d(lens::Lens, s::PlanarSurface, ray::Ray)

Intersection implementation for a lens with planar surfaces.
"""
function intersect3d(lens::SingletLens, s::PlanarSurface, ray::Ray{T}, _sign=(s === lens.front ? -1 : 1)) where T
    n = _sign*normal(lens)
    p0 = position(lens) .+ lens.edge_thickness/2 .* n

    nom = fast_dot3d((p0 .- position(ray)), n)
    denom = fast_dot3d(direction(ray), n)

    if iszero(denom)
        # line and plane are parallel
        if iszero(nom) && norm3d(position(ray) - p0) < mechanical_semi_diameter(s)
            # line is contained in plane
            return Intersection{T}(zero(T), T.(n), missing)
        else
            # line is not contained within plane --> no intersection
            return NoIntersection(T)
        end
    else
        # general case, single point of intersection
        t = nom / denom

        if t > eps(T)
            pos = position(ray) .+ t*direction(ray)
            if norm3d(pos - p0) < mechanical_semi_diameter(s)
                return Intersection{T}(t, T.(n), missing)
            end
        end
    end

    # if this is reached, no intersection falls within the mechanical aperture
    return NoIntersection(T)
end

"""
    intersect3d(lens::Lens, s::CylinderSurface, ray::Ray)

Intersection implementation for a lens with cylinder surfaces.
"""
function intersect3d(lens::SingletLens, s::CylinderSurface, ray::Ray{T}) where T
    _sign = (s === lens.front ? -1 : 1)
    p = position(lens) .+ _sign*mechanical_semi_diameter(s) .* s.d
    # aliases
    m = position(ray) - p
    n = direction(ra)
    d = s.d
    md = fast_dot3d(m, d)
    nd = fast_dot3d(n, d)
    # segment outside either end of the cylinder
    (md < eps(T) && md + nd < eps(T)) && return NoIntersection(T)
    (md > one(T) && md + nd > one(T)) && return NoIntersection(T)
    mn = fast_dot3d(m, n)
    a = one(T) - nd * nd
    k = fast_dot3d(m,m) - s.R^2
    c = k - md^2

    if abs(a) < eps(T)
        # ray is parallel to cylinder axis
        c > eps(T) && return NoIntersection(T)
        # ray intersects, find the intersection type
        if md < eps(T)
            # endcap 1
            t = -mn
        elseif md > one(T)
            # endcap 2
            t = nd - mn
        else
            # inside cylinder
            t = zero(T)
        end

        dir = position(ray) .+ t .* direction(ray)
        normal = (dir .- position(sphere)) * sign(s.R)
        return Intersection{T}(t, normalize3d(T.(normal)), missing)
    end

    b = mn - nd * md
    discr = b^2 - a*c
    if discr < eps(T)
        return NoIntersection(T)
    end

    t = (-b - √(discr)) / a
    if (md + t*nd < eps(T))
        if nd < eps(T)
            return NoIntersection(T)
        end
        t = -md / nd
        if k + 2*t(mn + t) > eps(T)
            return NoIntersection(T)
        end
    elseif md + t*nd > one(T)
        if nd > eps(T)
            return NoIntersection(T)
        end

        t = (one(T) - md) / nd
        if k + one(T) - 2*md + t*(2*(mn-nd) + t) > eps(T)
            return NoIntersection(T)
        end
    end

    dir = position(ray) .+ t .* direction(ray)
    normal = (dir .- fast_dot3d(dir, d) .* d) * sign(s.R)
    return Intersection{T}(t, normalize3d(T.(normal)), missing)
end
