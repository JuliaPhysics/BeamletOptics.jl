"""
    Mesh{T<:Number}

Contains the STL mesh information for an arbitrary object, that is the `vertices` that make up the mesh and
a matrix of `faces`, i.e. the connectivity matrix of the mesh. The data is read in using the `FileIO.jl` and
`MeshIO.jl` packages. Translations and rotations of the mesh are directly saved in absolute coordinates in the
vertex matrix. For orientation and translation tracking, a positional (`pos`) and directional (`dir`) matrix are stored.
"""
mutable struct Mesh{T<:Number}
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
    return Mesh{T}(vertices * scale, faces, Matrix{T}(I, 3, 3), T.([0, 0, 0]), scale)
end

"""
    orthogonal3d(object::Mesh, fID::Int)

Returns a vector with unit length that is perpendicular to the target face according to
the right-hand rule. The vertices must be listed row-wise within the face matrix.
"""
function orthogonal3d(object::Mesh, fID::Int)
    face = object.vertices[object.faces[fID, :], :]
    n = cross((face[2, :] - face[1, :]), (face[3, :] - face[1, :]))
    n /= norm(n)
    return n
end

"""
    translate3d!(mesh::Mesh, offset::Vector)

Mutating function that translates the vertices of an mesh in relation to the offset vector.
In addition, the mesh position vector is overwritten to reflect the new "center of gravity".
"""
function translate3d!(mesh::Mesh, offset::Vector)
    mesh.pos += offset
    mesh.vertices .+= offset'
    return nothing
end

"""
    xrotate3d!(mesh::Mesh, θ)

Mutating function that rotates the mesh around the **x-axis**.
The rotation is performed around the "center of gravity" axis.
"""
function xrotate3d!(mesh::Mesh, θ)
    # Calculate rotation matrix
    R = rotate3d([1, 0, 0], θ)
    # Translate mesh to origin, rotate, retranslate
    mesh.vertices .= (mesh.vertices .- mesh.pos') * R .+ mesh.pos'
    mesh.dir *= R
    return nothing
end

"""
    yrotate3d!(mesh::Mesh, θ)

Mutating function that rotates the mesh around the **y-axis**.
The rotation is performed around the "center of gravity" axis.
"""
function yrotate3d!(mesh::Mesh, θ)
    # Calculate rotation matrix
    R = rotate3d([0, 1, 0], θ)
    # Translate mesh to origin, rotate, retranslate
    mesh.vertices .= (mesh.vertices .- mesh.pos') * R .+ mesh.pos'
    mesh.dir *= R
    return nothing
end

"""
    zrotate3d!(mesh::Mesh, θ)

Mutating function that rotates the mesh around the **z-axis**.
The rotation is performed around the "center of gravity" axis.
"""
function zrotate3d!(mesh::Mesh, θ)
    # Calculate rotation matrix
    R = rotate3d([0, 0, 1], θ)
    # Translate mesh to origin, rotate, retranslate
    mesh.vertices .= (mesh.vertices .- mesh.pos') * R .+ mesh.pos'
    mesh.dir *= R
    return nothing
end

"""
    scale3d!(mesh::Mesh, scale)

Allows rescaling of mesh data around "center of gravity".
"""
function scale3d!(mesh::Mesh, scale)
    # Translate mesh to origin, scale, return to orig. pos.
    mesh.vertices .= (mesh.vertices .- mesh.pos') .* scale .+ mesh.pos'
    mesh.scale = scale
    return nothing
end

"""
    reset_translation3d!(mesh::Mesh)

Resets all previous translations and returns the mesh back to the global origin.
"""
function reset_translation3d!(mesh::Mesh)
    mesh.vertices .-= mesh.pos'
    mesh.pos = zeros(eltype(mesh.pos), 3)
    return nothing
end

"""
    reset_rotation3d!(mesh::Mesh)

Resets all previous rotations around the current offset.
"""
function reset_rotation3d!(mesh::Mesh)
    R = inv(mesh.dir)
    mesh.vertices = (mesh.vertices .- mesh.pos') * R .+ mesh.pos'
    mesh.dir = Matrix{eltype(mesh.dir)}(I, 3, 3)
    return nothing
end

abstract type AbstractMesh end

# Enforces that objects have to have the field mesh or implement `mesh`.
mesh(object::AbstractMesh) = object.mesh

translate3d!(object::AbstractMesh, offset::Vector) = translate3d!(mesh(object), offset)
xrotate3d!(object::AbstractMesh, θ) = xrotate3d!(mesh(object), θ)
yrotate3d!(object::AbstractMesh, θ) = yrotate3d!(mesh(object), θ)
zrotate3d!(object::AbstractMesh, θ) = zrotate3d!(mesh(object), θ)
scale3d!(object::AbstractMesh, scale) = scale3d!(mesh(object), scale)
# Mesh intersection
intersect3d(object::AbstractMesh, ray::Ray) = intersect3d(mesh(object), ray)
# Utils
orthogonal3d(object::AbstractMesh, fID::Int) = orthogonal3d(mesh(object), fID)
reset_translation3d!(object::AbstractMesh) = reset_translation3d!(mesh(object))
reset_rotation3d!(object::AbstractMesh) = reset_rotation3d!(mesh(object))