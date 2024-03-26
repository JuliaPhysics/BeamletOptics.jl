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

# function parabasal_ray_parameters(agb::SCDI.AstigmaticGaussianBeamlet, p0::AbstractArray, i::Int)
#     pn = agb.c.rays[i].dir
#     h1_real, u1_real = ray_to_plane_projection(p0, pn, agb.dxp.rays[i].pos, agb.dxp.rays[i].dir)
#     h1_imag, u1_imag = ray_to_plane_projection(p0, pn, agb.wxp.rays[i].pos, agb.wxp.rays[i].dir)
#     h1 = h1_real + im * h1_imag
#     u1 = u1_real + im * u1_imag

#     h2_real, u2_real = ray_to_plane_projection(p0, pn, agb.dyp.rays[i].pos, agb.dyp.rays[i].dir)
#     h2_imag, u2_imag = ray_to_plane_projection(p0, pn, agb.wyp.rays[i].pos, agb.wyp.rays[i].dir)
#     h2 = h2_real + im * h2_imag
#     u2 = u2_real + im * u2_imag

#     return h1, u1, h2, u2, p0
# end

# function parabasal_ray_parameters(agb::SCDI.AstigmaticGaussianBeamlet, z::Real)
#     p0, i = SCDI.point_on_beam(agb, z)
#     return parabasal_ray_parameters(agb, p0, i)
# end

# function waist_parameters(agb::SCDI.AstigmaticGaussianBeamlet, z::Real)
#     h1, ~, h2, ~, p0 = parabasal_ray_parameters(agb, z)
#     w1 = real(h1) + imag(h1)
#     w2 = real(h2) + imag(h2)
#     return p0, w1, w2
# end

# lagrange_invariant(h1, u1, h2, u2) = abs(dot(h1, u2) - dot(h2, u1))

# function lagrange_invariant(agb::SCDI.AstigmaticGaussianBeamlet; flen=0.1, n=100)
#     l = length(agb) + flen
#     zs = LinRange(0, l, n)
#     invariant = similar(zs)
#     for (i, z) in enumerate(zs)
#         h1, u1, h2, u2, ~ = parabasal_ray_parameters(beam_2, z)
#         invariant[i] = lagrange_invariant(h1, u1, h2, u2)
#     end
#     return invariant, zs
# end

# """See https://en.wikipedia.org/wiki/Ellipse"""
# ellipse(t, a, b, c) = a + b * cos(t) + c * sin(t)

# function SCDI.render_beam!(axis, agb::SCDI.AstigmaticGaussianBeamlet; flen = 0.1, color=:red, toggle=true)
#     t = LinRange(0, 2pi, 20)
#     l = length(agb) + flen
#     for z = LinRange(0, l, 100)
#         p0, w1, w2 = waist_parameters(agb, z)
#         d = ellipse.(t, Ref(p0), Ref(w1), Ref(w2))
#         lines!(axis, d, color=color)
#     end
# end

# function SCDI.electric_field(agb::SCDI.AstigmaticGaussianBeamlet, r, z)
#     p0, i = SCDI.point_on_beam(agb, z)
#     # test if r is point on orthogonal plane at z
#     if !isorthogonal(SCDI.direction(SCDI.rays(agb.c)[i]), r)
#         return error("r must lie in plane at p0/dir")
#     end
#     # calculate reduced field as per Greynolds/Worku
#     h1, u1, h2, u2, ~ = parabasal_ray_parameters(agb, p0, i)
#     w = norm(cross(h1, r))^2 * dot(u2, r) - norm(cross(h2,r ))^2 * dot(u1, r)
#     w /= 2 * norm(cross(h1, h2))^2
#     k = 2π/SCDI.wavelength(SCDI.rays(agb.c)[i])
#     ψ = 1/sqrt(norm(cross(h1, h2)))^2 * exp(im * k * w)
#     return ψ
# end