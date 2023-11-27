### Helper functions for testing

function Cube(scale::Real; T = Float64)
    vertices = [0 0 0
        1 0 0
        1 1 0
        0 1 0
        0 1 1
        1 1 1
        1 0 1
        0 0 1]
    faces = [1 3 2
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
        1 2 7]
    return SCDI.Mesh{T}(uuid4(),
        vertices .* scale,
        faces,
        Matrix{T}(I, 3, 3),
        T.([0, 0, 0]),
        scale)
end

struct ReflectiveCube{S <: SCDI.AbstractShape} <: SCDI.AbstractReflectiveOptic
    id::UUID
    shape::S
end

function RetroMesh(scale::Real; T = Float64)
    vertices = [0 0 0
        1 0 0
        0 1 0
        0 0 1]
    faces = [1 3 2
        1 4 3
        1 2 4]
    return SCDI.Mesh{T}(uuid4(),
        vertices .* scale,
        faces,
        Matrix{T}(I, 3, 3),
        T.([0, 0, 0]),
        scale)
end

struct RetroReflector{S <: SCDI.AbstractShape} <: SCDI.AbstractReflectiveOptic
    id::UUID
    shape::S
end

RetroReflector(scale) = RetroReflector(uuid4(), RetroMesh(scale))

function PlanoMirror(scale::T) where {T <: Real}
    sz = 0.5
    vertices = [sz 0 sz
        sz 0 -sz
        -sz 0 -sz
        -sz 0 sz]
    faces = [1 2 4
        2 3 4]
    shape = SCDI.Mesh{T}(uuid4(),
        vertices .* scale,
        faces,
        Matrix{T}(I, 3, 3),
        T.([0, 0, 0]),
        scale)
    return SCDI.Mirror(uuid4(), shape)
end
