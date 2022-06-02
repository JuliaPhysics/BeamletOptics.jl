"""
    Ray{T<:Number}

A mutable struct that contains a **position vector** `pos`, a **directional vector** `dir`
and a **length** variable `t` that are used to describe a generic ray as `pos+t*dir`.
The directional vector is required/adjusted to have unit length, i.e. `abs(dir) == 1`.
"""
mutable struct Ray{T<:Number}
    pos::Vector{T}
    dir::Vector{T}
    len::T
    @doc """
        Ray{T}(pos::Vector{T}, dir::Vector{T}) where T

    Parametric type constructor for struct Ray. Takes in a position vector `pos` and directional vector `dir`, which is scaled
    to unit length. The initial ray length `t` is set to `Inf`.
    """
    function Ray{T}(pos::Vector{T}, dir::Vector{T}) where T
        @assert norm(dir) != 0 "Illegal vector for direction"
        new{T}(pos, dir/norm(dir), Inf)
    end
end

"""
    Ray(pos::Vector{T}, dir::Vector{G}) where {T<:Union{Int, Float64}, G<:Union{Int, Float64}}

Concrete implementation of Ray struct for `Float64` data type. Accepts all combinations where `pos` and `dir`
are of type `Float64` and/or the primitive `Int`. Allows seperate implementation for `Float32`, etc.
Promotes input types to `Float64`.
"""
function Ray(pos::Vector{T}, dir::Vector{G}) where {T<:Union{Int, Float64}, G<:Union{Int, Float64}}
    return Ray{Float64}(Float64.(pos), Float64.(dir))
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

"""
    Geometry{T<:Number}

Contains the STL mesh information for an arbitrary object, that is the `vertices` that make up the mesh and
a matrix of `faces`, i.e. the connectivity matrix of the mesh. The data is read in using the `FileIO.jl` and 
`MeshIO.jl` packages. Translations and rotations of the mesh are directly saved in absolute coordinates in the
vertex matrix. For reference, a positional (`pos`) and directional (`dir`) vector are stored.  
"""
mutable struct Geometry{T<:Number}
    vertices::Matrix{T}
    faces::Matrix{Int}
    pos::Vector{T}
    dir::Vector{T}
    @doc """
        Geometry{T}(data; scale=1e-3) where T

    Parametric type constructor for struct Geometry. Takes data of type `GeometryBasics.Mesh` and extracts the 
    vertices and faces. The mesh is initialized at the global origin.
    """
    function Geometry{T}(data; scale=1e-3) where T
        # Read data from shitty mesh format to matrix format
        numEl = length(data)
        vertices = Matrix{T}(undef, numEl*3, 3)
        faces = Matrix{Int}(undef, numEl, 3)
        for i = 1:numEl
            for j = 1:3
                vertices[3(i-1)+j, :] = data[i][j]
                faces[i,j] = 3(i-1)+j
            end
        end
        # Initialize geometry at origin with orientation [1,0,0] scaled to mm
        # Origin and direction are converted to type T
        new{T}(vertices*scale, faces, T.([0,0,0]), T.([1,0,0]))
    end
end

"""
    Geometry(data; scale=1e-3)

Concrete implementation of Geometry struct. Data type of Geometry is variably selected based on
type of vertex data (i.e `Float32`). Mesh data is scaled by factor 1e-3, assuming m scale. 
"""
function Geometry(data; scale=1e-3)
    T = typeof(data[1][1][1])
    return Geometry{T}(data; scale=scale)
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
reference vector. Vector orientation is determined according to right-hand rule.
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
        return Inf
    end
    return t
end

"""
    intersect3d(object::Geometry, ray::SCDI.Ray)

This function is a generic implementation to check if a ray intersects the object geometry.
If true, the **distance** `t` is returned, where the location of intersection is `ray.pos+t*ray.dir`.
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