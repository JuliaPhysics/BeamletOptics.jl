"""
    AbstractGeometry

Dispatch handle for all types that require geometry representation, i.e. translation, rotation, etc. If a subtype inherits from `AbstractGeometry`, 
it is **required to have a composite type field** `geometry::Geometry` in order for this design scheme to work.
"""
abstract type AbstractGeometry end

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
    scale::T
    @doc """
        Geometry{T}(mesh) where T

    Parametric type constructor for struct Geometry. Takes data of type `GeometryBasics.Mesh` and extracts the 
    vertices and faces. The mesh is initialized at the global origin.
    """
    function Geometry{T}(mesh) where T
        # Read data from shitty mesh format to matrix format
        numEl = length(mesh)
        vertices = Matrix{T}(undef, numEl*3, 3)
        faces = Matrix{Int}(undef, numEl, 3)
        for i = 1:numEl
            for j = 1:3
                vertices[3(i-1)+j, :] = mesh[i][j]
                faces[i,j] = 3(i-1)+j
            end
        end
        # Initialize geometry at origin with orientation [1,0,0] scaled to mm
        # Origin and direction are converted to type T
        scale::T = 1e-3
        new{T}(vertices*scale, faces, T.([0,0,0]), T.([1,0,0]), scale)
    end
end

"""
    Geometry(mesh)

Concrete implementation of Geometry struct. Data type of Geometry is variably selected based on
type of vertex data (i.e `Float32`). Mesh data is scaled by factor 1e-3, assuming m scale. 
"""
function Geometry(mesh)
    T = typeof(mesh[1][1][1])
    return Geometry{T}(mesh)
end

"""
    @Geometry

A macro that automatically bequeaths the AbstractGeometry dispatch handle and includes the `x.geometry` field for type composition
with the Geometry type. This is an experimental feature and might be removed in the future. 

```julia
struct MyTyp
    x
end
```

expands to

```julia
struct MyTyp <: AbstractGeometry
    geometry::Geometry
    x
end
```
"""
macro Geometry(structure)
    # Insert inheritance from AbstractGeometry
    structure.args[2] = Expr(:<:, structure.args[2], AbstractGeometry)
    # Insert composite type field geometry::Geometry (requires promotion of Expr)
    f1, f2 = promote(structure.args[3].args, [:(geometry::Geometry)])
    structure.args[3].args = pushfirst!(f1, f2[1])
    # Check if type has internal constructor (and modify)
    try
        if structure.args[3].args[end].head === :function
            # Point to struct and constructor function call
            fcall = structure.args[3].args[end].args[1]
            ncall = structure.args[3].args[end].args[2].args[end]
            # Insert mesh variable at first position in fcall and ncall
            fcall.args = append!([fcall.args[1], :mesh], fcall.args[2:end])
            ncall.args = append!([ncall.args[1], :(Geometry(mesh))], ncall.args[2:end])
        end
    catch
        # do nothing
        @info "@Geometry: no internal constructor was found for " structure.args[2]
    end
    eval(structure)
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
    translate3d!(object::T, offset::Vector) where T<:AbstractGeometry

Wrapper for translate3d!(object::Geometry, offset::Vector). 
"""
function translate3d!(object::T, offset::Vector) where T<:AbstractGeometry
    translate3d!(object.geometry, offset)
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
    xrotate3d!(object::T, θ) where T<:AbstractGeometry

Wrapper for xrotate3d!(object::Geometry, θ).
"""
function xrotate3d!(object::T, θ) where T<:AbstractGeometry
    xrotate3d!(object.geometry, θ)
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
    yrotate3d!(object::T, θ) where T<:AbstractGeometry
    
Wrapper for yrotate3d!(object::Geometry, θ).
"""
function yrotate3d!(object::T, θ) where T<:AbstractGeometry
    yrotate3d!(object.geometry, θ)
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
    zrotate3d!(object::T, θ) where T<:AbstractGeometry

Wrapper for zrotate3d!(object::Geometry, θ).
"""
function zrotate3d!(object::T, θ) where T<:AbstractGeometry
    zrotate3d!(object.geometry, θ)
end

"""
    scale3d!(object::Geometry, scale)

Allows rescaling of mesh data around "center of gravity".
"""
function scale3d!(object::Geometry, scale)
    # Translate geometry to origin, scale, return to orig. pos.
    object.vertices = (object.vertices .- object.pos') .* scale .+ object.pos'
    object.scale = scale
end

"""
    scale3d!(object::T, scale) where T<:AbstractGeometry

Wrapper for scale3d!(object::Geometry, scale).
"""
function scale3d!(object::T, scale) where T<:AbstractGeometry
    scale3d!(object.geometry, scale)
end