"""
    AbstractMesh <: AbstractShape

A generic type for an shape whose volume can be described by a mesh. Must have a field `mesh` of type `Mesh`. See also `Mesh{T}`.
"""
abstract type AbstractMesh{T} <: AbstractShape{T} end

vertices(mesh::AbstractMesh) = mesh.vertices
vertices!(mesh::AbstractMesh, vertices) = (mesh.vertices .= vertices)

faces(mesh::AbstractMesh) = mesh.faces
faces!(::AbstractMesh, ::Any) = nothing

scale(mesh::AbstractMesh) = mesh.scale
scale!(mesh::AbstractMesh, scale) = (mesh.scale = scale)

"""
    Mesh <: AbstractMesh

Contains the STL mesh information for an arbitrary shape, that is the `vertices` that make up the mesh and
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
    vertices::Matrix{T}
    faces::Matrix{Int}
    dir::SMatrix{3, 3, T, 9}
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
    return Mesh{T}(
        vertices * scale,
        faces,
        Matrix{T}(I, 3, 3),
        Point3{T}(0),
        scale)
end

"""
    translate3d!(mesh::AbstractMesh, offset)

Mutating function that translates the vertices of an mesh in relation to the offset vector.
In addition, the mesh position vector is overwritten to reflect the new "center of gravity".
"""
function translate3d!(mesh::AbstractMesh, offset)
    position!(mesh, position(mesh) .+ offset)
    vertices!(mesh, vertices(mesh) .+ offset')
    return nothing
end

"""
    rotate3d!(mesh, axis, θ)

Mutating function that rotates the `mesh` around the specified rotation `axis` by the angle `θ`.
"""
function rotate3d!(mesh::AbstractMesh, axis, θ)
    # Calculate rotation matrix
    R = rotate3d(axis, θ)
    # Translate mesh to origin, rotate (counter-clockwise), retranslate
    vertices!(mesh, (vertices(mesh) .- position(mesh)') * R' .+ position(mesh)')
    orientation!(mesh, R * orientation(mesh))
    return nothing
end

function xrotate3d!(mesh::AbstractMesh{T}, θ) where {T}
    rotate3d!(mesh, @SArray(T[one(T), zero(T), zero(T)]), θ)
end
function yrotate3d!(mesh::AbstractMesh{T}, θ) where {T}
    rotate3d!(mesh, @SArray(T[zero(T), one(T), zero(T)]), θ)
end
function zrotate3d!(mesh::AbstractMesh{T}, θ) where {T}
    rotate3d!(mesh, @SArray(T[zero(T), zero(T), one(T)]), θ)
end

"""
    align3d!(mesh, target_axis)

Aligns the local `mesh` y-axis onto the `target_axis`.
"""
function align3d!(mesh::AbstractMesh, target_axis)
    # Calculate rotation matrix
    R = align3d(orientation(mesh)[:,2], target_axis)
    # Translate mesh to origin, rotate (counter-clockwise), retranslate
    vertices!(mesh, (vertices(mesh) .- position(mesh)') * R' .+ position(mesh)')
    orientation!(mesh, orientation(mesh) * R)
    return nothing
end

"""
    scale3d!(mesh::AbstractMesh, scale)

Allows rescaling of mesh data around "center of gravity".
"""
function scale3d!(mesh::AbstractMesh, scale)
    # Translate mesh to origin, scale, return to orig. pos.
    vertices!(mesh, (vertices(mesh) .- position(mesh)') .* scale .+ position(mesh)')
    scale!(mesh, scale)
    return nothing
end

"""
    reset_translation3d!(mesh::AbstractMesh)

Resets all previous translations and returns the mesh back to the global origin.
"""
function reset_translation3d!(mesh::AbstractMesh{T}) where {T}
    translate3d!(mesh, -position(mesh))
    return nothing
end

"""
    reset_rotation3d!(mesh::AbstractMesh)

Resets all previous rotations around the current offset.
"""
function reset_rotation3d!(mesh::AbstractMesh{T}) where {T}
    # Calculate rotation reset angle (thx LLMs)
    R = orientation(mesh)
    θ = acos(clamp((tr(R)-1)/2, -1, 1))
    if iszero(θ)
        return nothing
    end
    # Calculate rotation reset axis
    axis = 1/(2*sin(θ)) * [R[3,2]-R[2,3], R[1,3]-R[3,1], R[2,1]-R[1,2]]
    # Reset mesh rotation
    rotate3d!(mesh, axis, -θ)
    # Reset orientation field
    orientation!(mesh, Matrix{T}(I, 3, 3))
    return nothing
end

"""
    set_new_origin3d!(mesh::AbstractMesh)

Resets the mesh `dir`ectional matrix and `pos`ition vector to their initial values.\\
**Warning: this operation is non-reversible!**
"""
function set_new_origin3d!(mesh::AbstractMesh{T}) where {T}
    orientation!(mesh, Matrix{T}(I, 3, 3))
    position!(mesh, zeros(T, 3))
    return nothing
end

"""
    normal3d(mesh::AbstractMesh, fID::Int)

Returns a vector with unit length that is perpendicular to the target `face`` according to
the right-hand rule. The vertices must be listed row-wise within the face matrix.
"""
function normal3d(mesh::AbstractMesh{T}, fID::Int) where{T}
    @views begin
        face = vertices(mesh)[faces(mesh)[fID, :], :]
        n = cross(
            (Point3{T}(face[2, :]) - Point3{T}(face[1, :])),
            (Point3{T}(face[3, :]) - Point3{T}(face[1, :]))
        )
    end
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
    if (u < 0 - kϵ) || (u > 1 + kϵ)
        return T(Inf)
    end
    # Compute normalized v and reject if less than 0 or greater than 1
    Qv = cross(Tv, E1)
    v = dot(direction(ray), Qv) * invDet
    if (v < 0 - kϵ) || (u + v > 1 + kϵ)
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
    intersect3d(mesh::Mesh, ray::Ray)

This function is a generic implementation to check if a ray intersects the shape mesh.\\
"""
function intersect3d(mesh::AbstractMesh{M},
        ray::AbstractRay{R}) where {M <: Real, R <: Real}
    numEl = size(faces(mesh), 1)
    # allocate all intermediate vectors once (note that this is NOT THREAD-SAFE)
    T = promote_type(M, R)
    fID::Int = 0
    t0::T = Inf
    for i in 1:numEl
        face = @views vertices(mesh)[faces(mesh)[i, :], :]
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
        face = @views vertices(mesh)[faces(mesh)[fID, :], :]
        normal = normal3d(mesh, fID)
        return Intersection(t0, normalize(T.(normal)), mesh)
    end
end

## A collection of mesh constructors

"""
    RectangularFlatMesh(width, height)

Creates a 2D rectangular [`Mesh`](@ref) that is centered around the origin and aligned with respect to the **y-axis**.
**Vertex normals are parallel to the positive y-axis.**

# Inputs

- `width`: width along the x-axis in [m]
- `height`: height along the z-axis in [m]
"""
function RectangularFlatMesh(width::W, height::H) where {W<:Real, H<:Real}
    T = promote_type(W, H)
    x = width/2
    z = height/2
    vertices = [
        x 0 z
        x 0 -z
        -x 0 -z
        -x 0 z
    ]
    faces = [
        1 2 4
        2 3 4
    ]
    return Mesh{T}(
        vertices,
        faces,
        Matrix{T}(I, 3, 3),
        zeros(T, 3),
        one(T)
    )
end

"""
    QuadraticFlatMesh(width)

Creates a 2D quadratic [`Mesh`](@ref). Refer to [`RectangularFlatMesh`](@ref) for more information.
"""
QuadraticFlatMesh(width::Real) = RectangularFlatMesh(width, width)

"""
    CircularFlatMesh(radius, n)

Creates a 2D rectangular [`Mesh`](@ref) that is centered around the origin and aligned with respect to the **negative y-axis**.

# Inputs

- `radius`: of the mesh in [m]
- `n`: slice discretization factor (higher equals better resolution)
"""
function CircularFlatMesh(radius::T, n::Int=30) where T<:Real
    # calculate vertices
    xs = [cos(x)*radius for x in LinRange(0, 2pi * (n-1)/n, n)]
    ys = zeros(n+1)
    zs = [sin(x)*radius for x in LinRange(0, 2pi * (n-1)/n, n)]
    # add center vertex
    pushfirst!(xs, 0)
    pushfirst!(zs, 0)
    vertices = Matrix{Float64}(undef, n+1, 3)
    vertices[:, 1] = xs
    vertices[:, 2] = ys
    vertices[:, 3] = zs
    # sort faces
    faces = Matrix{Int64}(undef, n, 3)
    for i = 2:n+1
        faces[i-1, :] = [1, i, i+1]
    end
    # correct last entry
    faces[end] = 2
    return BeamletOptics.Mesh{T}(
        vertices,
        faces,
        Matrix{T}(I, 3, 3),
        zeros(T, 3),
        one(T)
    )
end

"""
    CuboidMesh(x, y, z, θ=π/2)

Constructs the [`Mesh`](@ref) of a rectangular cuboid as per the dimensions specified by `x`, `y` and `z`.
In addition, one side of the mesh can be tilted by an angle `θ` in order to generate the mesh of a rhomb.
The mesh is initialized such that one corner of the cuboid lies at the origin.

# Arguments

- `x, y, z`: dimensions for the cube in [m]
- `θ`: parallel tilt angle
"""
function CuboidMesh(x::X, y::Y, z::Z, θ::Real=π/2) where {X<:Real, Y<:Real, Z<:Real}
    T = promote_type(X, Y, Z)
    dx = cos(θ)*y
    vertices = [
        0 0 0
        x 0 0
        x + dx y 0
        0 + dx y 0
        0 + dx y z
        x + dx y z
        x 0 z
        0 0 z
    ]
    faces = [
        1 3 2
        1 4 3
        3 4 5
        3 5 6
        2 3 6
        2 6 7
        1 8 5
        1 5 4
        6 5 8
        6 8 7
        1 7 8
        1 2 7
    ]
    return Mesh{T}(
        T.(vertices),
        faces,
        Matrix{T}(I, 3, 3),
        zeros(T, 3),
        one(T))
end

"""Refer to [`CuboidMesh`](@ref)."""
function CubeMesh(scale::Real)
    scale = Float64(scale)
    return CuboidMesh(scale, scale, scale)
end
