# This file contains a collection of convenience constructors for various optical elements.

function RectangularPlanoMirror2D(scale::T) where {T <: Real}
    shape = QuadraticFlatMesh(scale)
    return Mirror(uuid4(), shape)
end

function RoundPlanoMirror(diameter::D, thickness::T) where {D<:Real,T<:Real}
    shape = CylinderSDF(diameter/2, thickness/2)
    return Mirror(uuid4(), shape)
end

function RectangularPlanoMirror(width::W, height::H, thickness::T) where {W<:Real,H<:Real,T<:Real}
    shape = CuboidMesh((width, thickness, height))
    translate3d!(shape, [
        -width/2,
        -thickness/2,
        -height/2,
    ])
    set_new_origin3d!(shape)
    return Mirror(uuid4(), shape)
end

function SquarePlanoMirror(width::W, thickness::T) where {W<:Real,T<:Real}
    return RectangularPlanoMirror(width, width, thickness)
end

struct ReflectiveCube{S <: AbstractShape} <: AbstractReflectiveOptic
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
    return Mesh{T}(uuid4(),
        vertices .* scale,
        faces,
        Matrix{T}(I, 3, 3),
        T.([0, 0, 0]),
        scale)
end

struct RetroReflector{S <: AbstractShape} <: AbstractReflectiveOptic
    id::UUID
    shape::S
end

RetroReflector(scale) = RetroReflector(uuid4(), RetroMesh(scale))

"""
    RectangularPlateBeamSplitter(width, thickness, n=1.5; eps_separation=1e-9)

# FIXME
"""
function RectangularPlateBeamSplitter(width::Real, thickness::Real, n::Function; eps_separation=1e-9)
    coating = SCDI.ThinBeamSplitter(width)
    shape = SCDI.CuboidMesh((width, thickness, width))
    SCDI.translate3d!(shape, [
        -width/2,
        -thickness/2,
        -width/2,
    ])
    SCDI.translate3d!(coating, [0, thickness/2 + eps_separation, 0])
    substrate = SCDI.Prism(uuid4(), shape, n)
    return SCDI.ObjectGroup([coating, substrate])
end

RectangularPlateBeamSplitter(w::Real, t::Real, n::Real; eps_separation=1e-9) = RectangularPlateBeamSplitter(w, t, λ->n; eps_separation)

function RectangularCompensatorPlate(width::W, height::H, thickness::T, n::Function) where {W<:Real,H<:Real,T<:Real}
    shape = CuboidMesh((width, thickness, height))
    translate3d!(shape, [
        -width/2,
        -thickness/2,
        -height/2,
    ])
    set_new_origin3d!(shape)
    return Prism(uuid4(), shape, n)
end

RectangularCompensatorPlate(w::Real, h::Real, t::Real, n::Real) = RectangularCompensatorPlate(w, h, t, λ->n)