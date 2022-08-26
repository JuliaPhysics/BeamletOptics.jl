### Helper functions for testing

mutable struct Cube{T} <: SCDI.AbstractMesh{T}
    mesh::SCDI.Mesh{T}
    dimension::T
end

function Cube(scale::Real; T=Float64)
    vertices = [
        0 0 0
        1 0 0
        1 1 0
        0 1 0
        0 1 1
        1 1 1
        1 0 1
        0 0 1
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
    return Cube{T}(
        SCDI.Mesh{T}(
            vertices .* scale,
            faces,
            Matrix{T}(I, 3, 3),
            T.([0, 0, 0]),
            scale
        ),
        scale
    )
end
