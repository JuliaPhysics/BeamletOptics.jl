struct Ray
    pos::Vector
    dir::Vector
    function Ray(pos, dir)
        @assert norm(dir) != 0 "Illegal vector for direction"
        new(pos, dir/norm(dir))
    end
end

struct Beamlet
    chief::Ray
    divergence::Ray
    waist::Ray

    # Beamlet constructor
    function Beamlet(chief::Ray, λ, w0)
        # All rays are initialized parallel to the x,y-plane
        # Divergence angle in rad
        θ = λ/(π*w0)
        # Divergence ray
        divergence = Ray(chief.pos, rotate3d([0,0,1], θ)*chief.dir)
        # Waist ray
        waist = Ray(chief.pos+orthogonal3d(chief.dir,[0,0,1])*w0, chief.dir)
        new(chief, divergence, waist)
    end
end

mutable struct Geometry
    vertices::Matrix
    faces::Matrix
    pos::Vector
    dir::Vector
    function Geometry(data)
        # Read data from shitty mesh format to matrix format
        numEl = length(data)
        vertices = Matrix{Float64}(undef, numEl*3, 3)
        faces = Matrix{Int}(undef, numEl, 3)
        for i = 1:numEl
            for j = 1:3
                vertices[3(i-1)+j, :] = data[i][j]
                faces[i,j] = 3(i-1)+j
            end
        end
        # Initialize geometry at origin with orientation [1,0,0] scaled to mm
        new(vertices*1e-3, faces, [0,0,0], [1,0,0])
    end
end

struct Mirror
    pos::Vector
    dir::Vector
    function Mirror(pos, dir)
        @assert norm(dir) != 0 "Illegal vector for direction"
        new(pos, dir/norm(dir))
    end
end

"""
    orthogonal3d(target::Vector, reference::Vector)

Returns a vector with unit length that is perpendicular to the target and an additional
reference vector. Vector orientation according to right-hand rule.
"""
function orthogonal3d(target::Vector, reference::Vector)
    n = cross(target, reference)
    n /= norm(n)
    return n
end

"""
    rotate3d(reference::Ray)

Returns the rotation matrix that will rotate a vector around the reference axis at an angle
θ in radians. Vector length is maintained. Rotation in clockwise direction?
"""
function rotate3d(reference::Vector, θ)
    cost = cos(θ)
    sint = sin(θ)
    ux, uy, uz = reference
    R = [
        cost+ux^2*(1-cost) ux*uy*(1-cost)-uz*sint ux*uz*(1-cost)+uy*sint;
        uy*ux*(1-cost)+uz*sint cost+uy^2*(1-cost) uy*uz*(1-cost)-ux*sint;
        uz*ux*(1-cost)-uy*sint uz*uy*(1-cost)+ux*sint cost+uz^2*(1-cost)
    ]
    return R
end

"""
    align3d(start::Vector, target::Vector)

Returns the rotation matrix R that will align the start vector to be parallel to the target vector.   
Based on ['Avoiding Trigonometry'](https://gist.github.com/kevinmoran/b45980723e53edeb8a5a43c49f134724) by Íñigo Quílez. The resulting matrix
was transposed due to column/row major issues. Vector length is maintained. This function is very fast.
"""
function align3d(start::Vector, target::Vector)
    start /= norm(start)
    target /= norm(target)
    rx, ry, rz = cross(target, start)
    # if start and target are already (almost) parallel return unity
    if (abs(rx) < 1e-9) & (abs(ry) < 1e-9) & (abs(rz) < 1e-9)
        return Matrix(1.0I,3,3)
    end
    cosA = dot(start, target)
    k = 1/(1+cosA)
    R = [
        rx^2*k+cosA rx*ry*k+rz rx*rz*k-ry;
        ry*rx*k-rz ry^2*k+cosA ry*rz*k+rx;
        rz*rx*k+ry rz*ry*k-rx rz^2*k+cosA
    ]
    return R
end

"""
    translate3d!(object::Geometry, offset::Vector)

Mutating function that translates the vertices of an object in relation to the offset vector.
In addition, the objection position vector is overwritten to reflect the new "center of gravity".
"""
function translate3d!(object::Geometry, offset::Vector)
    object.pos += offset
    object.vertices .+= offset'
end

"""
    xrotate3d!(object::Geometry, θ)

Mutating function that rotates the object geometry around the **x-axis**. 
The rotation is performed around the "center of gravtiy" axis.
"""
function xrotate3d!(object::Geometry, θ)
    # Calculate rotation matrix
    R = rotate3d([1,0,0], θ)
    # Translate geometry to origin, rotate, retranslate
    object.vertices = (object.vertices .- object.pos')*R .+ object.pos'
end

"""
    yrotate3d!(object::Geometry, θ)

Mutating function that rotates the object geometry around the **y-axis**. 
The rotation is performed around the "center of gravtiy" axis.
"""
function yrotate3d!(object::Geometry, θ)
    # Calculate rotation matrix
    R = rotate3d([0,1,0], θ)
    # Translate geometry to origin, rotate, retranslate
    object.vertices = (object.vertices .- object.pos')*R .+ object.pos'
end

"""
    zrotate3d!(object::Geometry, θ)

Mutating function that rotates the object geometry around the **z-axis**. 
The rotation is performed around the "center of gravtiy" axis.
"""
function zrotate3d!(object::Geometry, θ)
    # Calculate rotation matrix
    R = rotate3d([0,0,1], θ)
    # Translate geometry to origin, rotate, retranslate
    object.vertices = (object.vertices .- object.pos')*R .+ object.pos'
end

"""
    intersect3d(face::Matrix, ray::Ray; kϵ=1e-5)

An implementation of the (Möller-Trumbore algorithm)[https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection].
This algorithm evaluates the possible intersection between a ray and a face that is defined by three vertices. If no intersection occurs, Inf is returned.
kϵ is the abort threshold for backfacing and non-intersecting triangles. 
This algorithm is fast due to multiple breakout conditions.
"""
function intersect3d(face::Matrix, ray::SCDI.Ray; kϵ=1e-6)
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
        return inf 
    end
    return t
end

"""
    intersect3d(object::Geometry, ray::SCDI.Ray)

This function is a generic implementation to check if a ray intersects the object geometry.
If true, the distance **t** is returned, where the location of intersection is ray.pos+t*ray.dir.
"""
function intersect3d(object::Geometry, ray::Ray)
    numEl = size(object.faces)[1]
    t0 = Inf
    for i = 1:numEl
        face = object.vertices[object.faces[i,:],:]
        t = intersect3d(face, ray)
        # Return closest intersection
        if t < t0
            t0 = t
        end
    end
    return t0
end