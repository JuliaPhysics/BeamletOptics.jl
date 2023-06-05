"""
    AbstractMesh <: AbstractShape

A generic type for an object whose volume can be described by a mesh. Must have a field `mesh` of type `Mesh`. See also `Mesh{T}`.
"""
abstract type AbstractMesh{T} <: AbstractShape{T} end

vertices(mesh::AbstractMesh) = mesh.vertices
vertices!(mesh::AbstractMesh, vertices) = (mesh.vertices .= vertices)

faces(mesh::AbstractMesh) = mesh.faces
faces!(mesh::AbstractMesh, faces) = nothing

orientation(mesh::AbstractMesh) = mesh.dir
orientation!(mesh::AbstractMesh, dir) = (mesh.dir .= dir)

position(mesh::AbstractMesh) = mesh.pos
position!(mesh::AbstractMesh, pos) = (mesh.pos .= pos)

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
    dir::Matrix{T}
    pos::Vector{T}
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
    for i = 1:numEl
        for j = 1:3
            vertices[3(i-1)+j, :] = mesh[i][j]
            faces[i, j] = 3(i - 1) + j
        end
    end
    # Initialize mesh at origin with orientation [1,0,0] scaled to mm
    # Origin and direction are converted to type T
    scale::T = 1e-3
    return Mesh{T}(uuid4(), vertices * scale, faces, Matrix{T}(I, 3, 3), T.([0, 0, 0]), scale)
end

"""
    translate3d!(object::AbstractMesh, offset::Vector)

Mutating function that translates the vertices of an mesh in relation to the offset vector.
In addition, the mesh position vector is overwritten to reflect the new "center of gravity".
"""
function translate3d!(object::AbstractMesh, offset::Vector)
    position!(object, position(object) .+= offset)
    vertices!(object, vertices(object) .+= offset')
    return nothing
end

"""
    xrotate3d!(object::AbstractMesh, θ)

Mutating function that rotates the mesh around the **x-axis**.
The rotation is performed around the "center of gravity" axis.
"""
function xrotate3d!(object::AbstractMesh, θ)
    # Calculate rotation matrix
    R = rotate3d([1, 0, 0], θ)
    # Translate mesh to origin, rotate, retranslate
    vertices!(object, (vertices(object) .- position(object)') * R .+ position(object)')
    orientation!(object, orientation(object) * R)
    return nothing
end

"""
    yrotate3d!(object::AbstractMesh, θ)

Mutating function that rotates the mesh around the **y-axis**.
The rotation is performed around the "center of gravity" axis.
"""
function yrotate3d!(object::AbstractMesh, θ)
    # Calculate rotation matrix
    R = rotate3d([0, 1, 0], θ)
    # Translate mesh to origin, rotate, retranslate
    vertices!(object, (vertices(object) .- position(object)') * R .+ position(object)')
    orientation!(object, orientation(object) * R)
    return nothing
end

"""
    zrotate3d!(object::AbstractMesh, θ)

Mutating function that rotates the mesh around the **z-axis**.
The rotation is performed around the "center of gravity" axis.
"""
function zrotate3d!(object::AbstractMesh, θ)
    # Calculate rotation matrix
    R = rotate3d([0, 0, 1], θ)
    # Translate mesh to origin, rotate, retranslate
    vertices!(object, (vertices(object) .- position(object)') * R .+ position(object)')
    orientation!(object, orientation(object) * R)
    return nothing
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
function reset_translation3d!(object::AbstractMesh{T}) where T
    vertices!(object, vertices(object) .- position(object)')
    position!(object, zeros(T, 3))
    return nothing
end

"""
    reset_rotation3d!(object::AbstractMesh)

Resets all previous rotations around the current offset.
"""
function reset_rotation3d!(object::AbstractMesh{T}) where T
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
function set_new_origin3d!(object::AbstractMesh{T}) where T
    orientation!(object, Matrix{T}(I, 3, 3))
    position!(object, zeros(T, 3))
    return nothing
end

"""
    orthogonal3d(object::AbstractMesh, fID::Int)

Returns a vector with unit length that is perpendicular to the target `face`` according to
the right-hand rule. The vertices must be listed row-wise within the face matrix.
"""
function orthogonal3d(object::AbstractMesh, fID::Int)
    face = vertices(object)[faces(object)[fID, :], :]
    n = fast_cross3d((face[2, :] - face[1, :]), (face[3, :] - face[1, :]))
    return normalize3d(n)
end