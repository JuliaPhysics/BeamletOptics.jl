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
    scale3d!(object::Geometry, scale)

Allows rescaling of mesh data around "center of gravity".
"""
function scale3d!(object::Geometry, scale)
    # Translate geometry to origin, scale, return to orig. pos.
    object.vertices = (object.vertices .- object.pos') .* scale .+ object.pos'
    object.scale = scale
end

"""
    @Geometry

A macro that automatically includes the `geometry` field for type composition with the `Geometry` type.
Additionally, all wrapper functions required to access Geometry manipulation features are created. 
However, these can be overwritten if another form of dispatch is necessary.

```julia
struct MyType
    x
end
```

expands to

```julia
struct MyType
    geometry::Geometry
    x
end
```
"""
macro Geometry(type)
    @assert type.head === :struct "@Geometry only works with structs!"
    # Insert geometry field into struct
    @debug "Composition of type $(type.args[2]) with field geometry::Geometry..." type.args[2]
    pushfirst!(type.args[3].args, :(geometry::Geometry))
    try type.args[3].args[end].head === :function
        # Insert geometry variable into constructor
        @debug "Modified internal constructor of type $(type.args[2])..."
        insert!(type.args[3].args[end].args[1].args, 2, :(geometry::Geometry))
        insert!(type.args[3].args[end].args[2].args[end].args, 2, :(geometry::Geometry))
    catch
        # do nothing
    end
    eval(type)

    # Extract type name
    if typeof(type.args[2]) == Expr
        tname = type.args[2].args[1]
    else
        tname = type.args[2]
    end

    functions = quote
        # Euclidic manipulation operations
        translate3d!(object::eval($(tname)), offset::Vector) = translate3d!(object.geometry, offset)
        xrotate3d!(object::eval($(tname)), θ) = xrotate3d!(object.geometry, θ)
        yrotate3d!(object::eval($(tname)), θ) = yrotate3d!(object.geometry, θ)
        zrotate3d!(object::eval($(tname)), θ) = zrotate3d!(object.geometry, θ)
        scale3d!(object::eval($(tname)), scale) = scale3d!(object.geometry, scale)
        # Mesh intersection
        intersect3d(object::eval($(tname)), ray::Ray) = intersect3d(object.geometry, ray)
        # Utils
        orthogonal3d(object::eval($(tname)), fID::Int) = orthogonal3d(object.geometry, fID)
    end
    @debug "Generating wrapper functions for type $(tname)"
    eval(functions)
end