"""
    AbstractMesh <: AbstractShape

A generic type for an object whose volume can be described by a mesh. Must have a field `mesh` of type `Mesh`. See also `Mesh{T}`.
"""
abstract type AbstractMesh{T} <: AbstractShape{T} end

vertices(mesh::AbstractMesh) = mesh.vertices
vertices!(mesh::AbstractMesh, vertices) = (mesh.vertices .= vertices)

faces(mesh::AbstractMesh) = mesh.faces
faces!(::AbstractMesh, ::Any) = nothing

scale(mesh::AbstractMesh) = mesh.scale
scale!(mesh::AbstractMesh, scale) = (mesh.scale = scale)

"""
    Mesh{T} <: AbstractMesh{T}

Contains the STL mesh information for an arbitrary object, that is the `vertices` that make up the mesh and
a matrix of `faces`, i.e. the connectivity matrix of the mesh. The data is read in using the `FileIO.jl` and
`MeshIO.jl` packages. Translations and rotations of the mesh are directly saved in absolute coordinates in the
vertex matrix. For orientation and translation tracking, a `pos`itional and `dir`ectional matrix are stored.

# Fields
- `vertices`: (m x 3)-matrix that stores the edge points of all triangles
- `faces`: (n x 3)-matrix that stores the connectivity data for all faces
- `dir`: (3 x 3)-matrix that represents the current orientation of the mesh
- `pos`: 3-element vector that is used as the mesh location reference
- `scale`: scalar value that represents the current scale of the original mesh
"""
mutable struct Mesh{T} <: AbstractMesh{T}
    id::UUID
    vertices::Matrix{T}
    faces::Matrix{Int}
    dir::Mat{3, 3, T}
    pos::Point3{T}
    scale::T
end

"""
    Mesh(mesh)

Parametric type constructor for struct Mesh. Takes data of type `GeometryBasics.Mesh` and extracts the
vertices and faces. The mesh is initialized at the global origin. Data type of Mesh is variably selected based on
type of vertex data (i.e `Float32`). Mesh data is scaled by factor 1e-3, assuming m scale.
"""
function Mesh(mesh)
    # Determine mesh data type (i.e. Float32)
    T = typeof(mesh[1][1][1])
    # Read data from shitty mesh format to matrix format
    numEl = length(mesh)
    vertices = Matrix{T}(undef, numEl * 3, 3)
    faces = Matrix{Int}(undef, numEl, 3)
    for i in 1:numEl
        for j in 1:3
            vertices[3(i - 1) + j, :] = mesh[i][j]
            faces[i, j] = 3(i - 1) + j
        end
    end
    # Initialize mesh at origin with orientation [1,0,0] scaled to mm
    # Origin and direction are converted to type T
    scale::T = 1e-3
    return Mesh{T}(uuid4(),
        vertices * scale,
        faces,
        Matrix{T}(I, 3, 3),
        Point3{T}(0),
        scale)
end

"""
    translate3d!(object::AbstractMesh, offset)

Mutating function that translates the vertices of an mesh in relation to the offset vector.
In addition, the mesh position vector is overwritten to reflect the new "center of gravity".
"""
function translate3d!(object::AbstractMesh, offset)
    position!(object, position(object) .+ offset)
    vertices!(object, vertices(object) .+ offset')
    return nothing
end

"""
    rotate3d!(shape::AbstractMesh, axis, θ)

Mutating function that rotates the mesh around the specified rotation `axis` by the angle `θ`.
"""
function rotate3d!(shape::AbstractMesh, axis, θ)
    # Calculate rotation matrix
    R = rotate3d(axis, θ)
    # Translate mesh to origin, rotate, retranslate
    vertices!(shape, (vertices(shape) .- position(shape)') * R' .+ position(shape)')
    orientation!(shape, orientation(shape) * R)
    return nothing
end

function xrotate3d!(shape::AbstractMesh{T}, θ) where {T}
    rotate3d!(shape, @SVector(T[one(T), zero(T), zero(T)]), θ)
end
function yrotate3d!(shape::AbstractMesh{T}, θ) where {T}
    rotate3d!(shape, @SVector(T[zero(T), one(T), zero(T)]), θ)
end
function zrotate3d!(shape::AbstractMesh{T}, θ) where {T}
    rotate3d!(shape, @SVector(T[zero(T), zero(T), one(T)]), θ)
end

"""
    scale3d!(object::AbstractMesh, scale)

Allows rescaling of mesh data around "center of gravity".
"""
function scale3d!(object::AbstractMesh, scale)
    # Translate mesh to origin, scale, return to orig. pos.
    vertices!(object, (vertices(object) .- position(object)') .* scale .+ position(object)')
    scale!(object, scale)
    return nothing
end

"""
    reset_translation3d!(object::AbstractMesh)

Resets all previous translations and returns the mesh back to the global origin.
"""
function reset_translation3d!(object::AbstractMesh{T}) where {T}
    vertices!(object, vertices(object) .- position(object)')
    position!(object, zeros(T, 3))
    return nothing
end

"""
    reset_rotation3d!(object::AbstractMesh)

Resets all previous rotations around the current offset.
"""
function reset_rotation3d!(object::AbstractMesh{T}) where {T}
    R = inv(orientation(object))
    vertices!(object, (vertices(object) .- position(object)') * R .+ position(object)')
    orientation!(object, Matrix{T}(I, 3, 3))
    return nothing
end

"""
    set_new_origin3d!(object::AbstractMesh)

Resets the mesh `dir`ectional matrix and `pos`ition vector to their initial values.\\
**Warning: this operation is non-reversible!**
"""
function set_new_origin3d!(object::AbstractMesh{T}) where {T}
    orientation!(object, Matrix{T}(I, 3, 3))
    position!(object, zeros(T, 3))
    return nothing
end

"""
    normal3d(object::AbstractMesh, fID::Int)

Returns a vector with unit length that is perpendicular to the target `face`` according to
the right-hand rule. The vertices must be listed row-wise within the face matrix.
"""
function normal3d(object::AbstractMesh, fID::Int)
    face = vertices(object)[faces(object)[fID, :], :]
    n = cross((face[2, :] - face[1, :]), (face[3, :] - face[1, :]))
    return normalize(n)
end

"""
    MoellerTrumboreAlgorithm(face::Matrix, ray::Ray)

A culling implementation of the **Möller-Trumbore algorithm** for ray-triangle-intersection.
This algorithm evaluates the possible intersection between a `ray` and a `face` that is defined by three vertices.
If no intersection occurs, `Inf` is returned. `kϵ` is the abort threshold for backfacing and non-intersecting triangles.
`lϵ` is the threshold for negative values of `t`.
This algorithm is fast due to multiple breakout conditions.
"""
function MoellerTrumboreAlgorithm(face, ray::AbstractRay{T}; kϵ = 1e-9, lϵ = 1e-9) where {T}
    V1 = Point3(face[1, 1], face[1, 2], face[1, 3])
    V2 = Point3(face[2, 1], face[2, 2], face[2, 3])
    V3 = Point3(face[3, 1], face[3, 2], face[3, 3])

    E1 = V2 - V1
    E2 = V3 - V1
    Pv = cross(direction(ray), E2)
    Det = dot(E1, Pv)
    # Check if ray is backfacing or missing face
    if abs(Det) < kϵ
        # Adjust type of Inf
        return T(Inf) # typemax(T) ?
    end
    # Compute normalized u and reject if less than 0 or greater than 1
    Tv = position(ray) - V1
    invDet = 1 / Det
    u = dot(Tv, Pv) * invDet
    if (u < 0) || (u > 1)
        return T(Inf)
    end
    # Compute normalized v and reject if less than 0 or greater than 1
    Qv = cross(Tv, E1)
    v = dot(direction(ray), Qv) * invDet
    if (v < 0) || (u + v > 1)
        return T(Inf)
    end
    # Compute t (type def. for t to avoid Any)
    t::T = dot(E2, Qv) * invDet
    # Return intersection only if "in front of" ray origin
    if t < lϵ
        return T(Inf)
    end
    return t
end

"""
    intersect3d(object::Mesh, ray::Ray)

This function is a generic implementation to check if a ray intersects the object mesh.\\
"""
function intersect3d(object::AbstractMesh{M},
        ray::AbstractRay{R}) where {M <: Real, R <: Real}
    numEl = size(faces(object), 1)
    # allocate all intermediate vectors once (note that this is NOT THREAD-SAFE)
    T = promote_type(M, R)
    fID::Int = 0
    t0::T = Inf
    for i in 1:numEl
        face = @views vertices(object)[faces(object)[i, :], :]
        t = MoellerTrumboreAlgorithm(face, ray)
        # Return closest intersection
        if t < t0
            t0 = t
            fID = i
        end
    end
    if isinf(t0)
        return nothing
    else
        face = @views vertices(object)[faces(object)[fID, :], :]
        normal = normal3d(object, fID)
        return Intersection{T}(t0, normalize(T.(normal)), nothing)
    end
end
