# Waveplate Component
# Reuses the generic Coating system under the hood with a JonesCoating model.

"""
    Waveplate(shape, JMat)
    Waveplate(diameter, retardance; fast_axis_angle=0.0)
    Waveplate(width, height, retardance; fast_axis_angle=0.0)
    Waveplate(diameter, thickness, n, retardance; fast_axis_angle=0.0)
    Waveplate(width, height, thickness, n, retardance; fast_axis_angle=0.0)

Ideal retarder/waveplate. The device is represented by a zero-thickness
surface whose fast axis is aligned with the local `x`-axis and slow axis with the local `z`-axis.

Rotate the underlying shape around the y-axis to align the axes with the desired global directions.
"""
function Waveplate(shape::AbstractShape{T}, JMat::GlobalJonesBasis) where {T}
    model = JonesCoating(GlobalJonesBasis{Complex{T}}(JMat))
    return Coating(shape, model)
end

function Waveplate(shape::AbstractShape{T}, retardance::Real) where {T}
    JMat = XZBasis(exp(-im * retardance / 2), 0, 0, exp(im * retardance / 2))
    return Waveplate(shape, JMat)
end

function Waveplate(shape::AbstractShape{T}, f::Function) where {T}
    test_val = f(1e-6)
    if test_val isa Real
        JMat = λ -> XZBasis(exp(-im * f(λ) / 2), 0, 0, exp(im * f(λ) / 2))
        model = JonesCoating(JMat)
    else
        model = JonesCoating(f)
    end
    return Coating(shape, model)
end

function Waveplate(diameter::Real, retardance::Union{Real, Function, GlobalJonesBasis}; fast_axis_angle::Real=0.0)
    shape = CircularFlatMesh(diameter / 2)
    yrotate3d!(shape, -fast_axis_angle)
    return Waveplate(shape, retardance)
end

function Waveplate(width::Real, height::Real, retardance::Union{Real, Function, GlobalJonesBasis}; fast_axis_angle::Real=0.0)
    shape = RectangularFlatMesh(width, height)
    zrotate3d!(shape, π)
    set_new_origin3d!(shape)
    yrotate3d!(shape, -fast_axis_angle)
    return Waveplate(shape, retardance)
end

#=
Thick/Plate Waveplate types and methods
=#

abstract type AbstractPlateWaveplate{T} <: AbstractPlateComponent{T} end

"""
    RectangularPlateWaveplate{T, M} <: AbstractPlateWaveplate{T}

A rectangular waveplate with a bulk glass substrate and a single coated retarder face.
"""
struct RectangularPlateWaveplate{T, M} <: AbstractPlateWaveplate{T}
    substrate::Prism{T, BoxSDF{T}}
    coating::Coating{T, Mesh{T}, M}
end

function RectangularPlateWaveplate(
        width::Real,
        height::Real,
        thickness::Real,
        n::RefractiveIndex,
        retardance::Union{Real, Function};
        fast_axis_angle::Real=0.0
)
    substrate_shape = BoxSDF(width, thickness, height)
    substrate = Prism(substrate_shape, n)
    translate3d!(substrate, [0, thickness/2, 0])
    
    coating_shape = RectangularFlatMesh(width, height)
    zrotate3d!(coating_shape, π)
    set_new_origin3d!(coating_shape)
    yrotate3d!(coating_shape, -fast_axis_angle)
    
    JMat = if retardance isa Function
        λ -> XZBasis(exp(-im * retardance(λ) / 2), 0, 0, exp(im * retardance(λ) / 2))
    else
        XZBasis(exp(-im * retardance / 2), 0, 0, exp(im * retardance / 2))
    end
    coating = Waveplate(coating_shape, JMat)
    
    T = promote_type(typeof(width), typeof(height), typeof(thickness))
    M = typeof(coating.model)
    return RectangularPlateWaveplate{T, M}(substrate, coating)
end

"""
    RoundPlateWaveplate{T, M} <: AbstractPlateWaveplate{T}

A round (cylindrical) waveplate with a bulk glass substrate and a single coated retarder face.
"""
struct RoundPlateWaveplate{T, M} <: AbstractPlateWaveplate{T}
    substrate::Prism{T, PlanoSurfaceSDF{T}}
    coating::Coating{T, Mesh{T}, M}
end

function RoundPlateWaveplate(
        diameter::Real,
        thickness::Real,
        n::RefractiveIndex,
        retardance::Union{Real, Function};
        fast_axis_angle::Real=0.0
)
    substrate_shape = PlanoSurfaceSDF(thickness, diameter)
    substrate = Prism(substrate_shape, n)
    
    coating_shape = CircularFlatMesh(diameter / 2)
    yrotate3d!(coating_shape, -fast_axis_angle)
    
    JMat = if retardance isa Function
        λ -> XZBasis(exp(-im * retardance(λ) / 2), 0, 0, exp(im * retardance(λ) / 2))
    else
        XZBasis(exp(-im * retardance / 2), 0, 0, exp(im * retardance / 2))
    end
    coating = Waveplate(coating_shape, JMat)
    
    T = promote_type(typeof(diameter), typeof(thickness))
    M = typeof(coating.model)
    return RoundPlateWaveplate{T, M}(substrate, coating)
end

# Factory mapping to round/rectangular plate waveplates
function Waveplate(diameter::Real, thickness::Real, n::RefractiveIndex, retardance::Union{Real, Function}; fast_axis_angle::Real=0.0)
    return RoundPlateWaveplate(diameter, thickness, n, retardance; fast_axis_angle=fast_axis_angle)
end

function Waveplate(width::Real, height::Real, thickness::Real, n::RefractiveIndex, retardance::Union{Real, Function}; fast_axis_angle::Real=0.0)
    return RectangularPlateWaveplate(width, height, thickness, n, retardance; fast_axis_angle=fast_axis_angle)
end

#=
Waveplate Factory/Convenience Constructors
=#

"""
    HalfWavePlate(diameter; fast_axis_angle=0.0)
    HalfWavePlate(width, height; fast_axis_angle=0.0)
    HalfWavePlate(diameter, thickness, n; fast_axis_angle=0.0)
    HalfWavePlate(width, height, thickness, n; fast_axis_angle=0.0)

Creates a Half-Wave Plate (retardance of π) of flat or thick, round or rectangular geometry.
"""
function HalfWavePlate(diameter::Real; fast_axis_angle::Real=0.0)
    return Waveplate(diameter, π; fast_axis_angle=fast_axis_angle)
end

function HalfWavePlate(width::Real, height::Real; fast_axis_angle::Real=0.0)
    return Waveplate(width, height, π; fast_axis_angle=fast_axis_angle)
end

function HalfWavePlate(diameter::Real, thickness::Real, n::RefractiveIndex; fast_axis_angle::Real=0.0)
    return Waveplate(diameter, thickness, n, π; fast_axis_angle=fast_axis_angle)
end

function HalfWavePlate(width::Real, height::Real, thickness::Real, n::RefractiveIndex; fast_axis_angle::Real=0.0)
    return Waveplate(width, height, thickness, n, π; fast_axis_angle=fast_axis_angle)
end

"""
    QuarterWavePlate(diameter; fast_axis_angle=0.0)
    QuarterWavePlate(width, height; fast_axis_angle=0.0)
    QuarterWavePlate(diameter, thickness, n; fast_axis_angle=0.0)
    QuarterWavePlate(width, height, thickness, n; fast_axis_angle=0.0)

Creates a Quarter-Wave Plate (retardance of π/2) of flat or thick, round or rectangular geometry.
"""
function QuarterWavePlate(diameter::Real; fast_axis_angle::Real=0.0)
    return Waveplate(diameter, π/2; fast_axis_angle=fast_axis_angle)
end

function QuarterWavePlate(width::Real, height::Real; fast_axis_angle::Real=0.0)
    return Waveplate(width, height, π/2; fast_axis_angle=fast_axis_angle)
end

function QuarterWavePlate(diameter::Real, thickness::Real, n::RefractiveIndex; fast_axis_angle::Real=0.0)
    return Waveplate(diameter, thickness, n, π/2; fast_axis_angle=fast_axis_angle)
end

function QuarterWavePlate(width::Real, height::Real, thickness::Real, n::RefractiveIndex; fast_axis_angle::Real=0.0)
    return Waveplate(width, height, thickness, n, π/2; fast_axis_angle=fast_axis_angle)
end

# Spelling aliases for Waveplates branch compatibility
const HalfWaveplate = HalfWavePlate
const QuarterWaveplate = QuarterWavePlate
