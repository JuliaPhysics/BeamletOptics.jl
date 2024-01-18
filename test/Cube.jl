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

#
mutable struct AstigmaticGaussianBeamlet{T} <: SCDI.AbstractBeam{T, SCDI.PolarizedRay{T}}
    id::UUID
    c::SCDI.Beam{T, SCDI.PolarizedRay{T}}       # chief
    wx::SCDI.Beam{T, SCDI.Ray{T}}               # waist x
    wy::SCDI.Beam{T, SCDI.Ray{T}}               # waist y
    dx::SCDI.Beam{T, SCDI.Ray{T}}               # div. x
    dy::SCDI.Beam{T, SCDI.Ray{T}}               # div. y
    parent::SCDI.Nullable{AstigmaticGaussianBeamlet{T}}
    children::Vector{AstigmaticGaussianBeamlet{T}}
end

function AstigmaticGaussianBeamlet(
    c::SCDI.Beam{T, SCDI.PolarizedRay{T}},
    wx::SCDI.Beam{T, SCDI.Ray{T}},
    wy::SCDI.Beam{T, SCDI.Ray{T}},
    dx::SCDI.Beam{T, SCDI.Ray{T}},
    dy::SCDI.Beam{T, SCDI.Ray{T}}
    ) where {T <: Real}
return AstigmaticGaussianBeamlet{T}(uuid4(), c, wx, wy, dx, dy,
    nothing,
    Vector{AstigmaticGaussianBeamlet{T}}())
end

function AstigmaticGaussianBeamlet(
    position,
    direction,
    λ = 1000e-9,
    w0 = 1e-3;
    M2 = 1,
    E0 = [0, 0, 1])
    # Create orthogonal vectors for construction purposes (right-handed)
    s1 = SCDI.normal3d(direction)
    s2 = cross(direction, s1)
    # Divergence angle in rad
    θ = SCDI.divergence_angle(λ, w0, M2)
    # Waist rays
    wx = SCDI.Ray(position + s1 * w0, direction, λ)
    wy = SCDI.Ray(position + s2 * w0, direction, λ)
    # Divergence ray
    div_dir_x = normalize(direction + s1 * tan(θ))
    div_dir_y = normalize(direction + s2 * tan(θ))
    dx = SCDI.Ray(position, div_dir_x, λ)
    dy = SCDI.Ray(position, div_dir_y, λ)
    # Chief ray
    c = SCDI.PolarizedRay(position, direction, λ, E0)
    return AstigmaticGaussianBeamlet(SCDI.Beam(c), SCDI.Beam(wx), SCDI.Beam(wy), SCDI.Beam(dx), SCDI.Beam(dy))
end