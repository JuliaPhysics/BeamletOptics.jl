"""
    RetroMesh(scale::Real; T = Float64)

Creates an open tetrahedral [Mesh](@ref) with edges derived from the vertices of a unit cube.
Can be scaled with a `scale` factor. The data type for the vertices and internal computations
can be adjusted using `T` (default: Float64).
"""
function RetroMesh(scale::Real; T = Float64)
    vertices = [0 0 0
        1 0 0
        0 1 0
        0 0 1]
    faces = [1 3 2
        1 4 3
        1 2 4]
    return Mesh{T}(
        vertices .* scale,
        faces,
        Matrix{T}(I, 3, 3),
        T.([0, 0, 0]),
        scale)
end

"""
    Retroreflector

A `Retroreflector` reflects incoming rays back toward their source, independent of the incident angle.
The shape is represented by a tetrahedral [`Mesh`](@ref).

# Fields

- `mesh`: shape of the `Retroreflector`
"""
struct Retroreflector{T} <: AbstractReflectiveOptic{T, Mesh{T}}
    mesh::Mesh{T}
end

shape(rr::Retroreflector) = rr.mesh

"""
    Retroreflector(scale)

Spawns a [`Retroreflector`](@ref).

# Inputs

- `scale`: a scaling factor for the size of the retroreflector, e.g. `1e-3` for 1 mm
"""
Retroreflector(scale::Real) = Retroreflector(RetroMesh(scale))
