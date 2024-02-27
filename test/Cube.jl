### Helper functions for testing

function Cube(scale::NTuple{3, T}, θ::Real=π/2) where T<:Real
    x = scale[1]
    y = scale[2]
    z = scale[3]
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
    return SCDI.Mesh{T}(uuid4(),
        vertices,
        faces,
        Matrix{T}(I, 3, 3),
        zeros(T, 3),
        one(T))
end

function Cube(scale::Real)
    scale = Float64(scale)
    return Cube((scale, scale, scale))
end

function Rhombus(_width::T, _length::Real, _height::Real) where T<:Real
    w = _width/2
    l = _length/2
    h = _height/2
    vertices = [
        -w 0 -h/2
        0 -l -h/2
        w  0 -h/2
        0  l -h/2
        w  0  h/2
        0 -l  h/2
        -w 0  h/2
        0  l  h/2
    ]
    faces = [
        1 3 2
        1 4 3
        3 4 5
        3 5 6
        2 3 6
        2 6 7
        8 4 1
        8 5 4
        6 5 8
        6 8 7
        1 7 8
        1 2 7
    ]
    return SCDI.Mesh{T}(uuid4(),
        vertices,
        faces,
        Matrix{T}(I, 3, 3),
        zeros(T, 3),
        one(T))
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

ellipse(t, a, b, c) = a + b * cos(t) + c * sin(t)

function SCDI.render_beam!(axis, agb::SCDI.AstigmaticGaussianBeamlet; flen = 0.1, color=:red)
    l = length(agb) + flen
    for z = LinRange(0, l, 100)
        ~, ~, ~, ~, B = SCDI.gauss_parameters(beam, z)
        a = B[:, 1]
        b = B[:, 2]
        c = B[:, 3]
        t = LinRange(0, 2pi, 20)
        d = ellipse.(t, Ref(a), Ref(b), Ref(c))
        lines!(axis, d, color=color)
    end
end