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

abstract type AbstractPlateWaveplate{T, N} <: AbstractRefractiveOptic{T, N} end

"""
    RectangularPlateWaveplate{T, S, N, C} <: AbstractPlateWaveplate{T, N}

A rectangular waveplate with a bulk glass substrate and a retarder coating on its front face.
"""
struct RectangularPlateWaveplate{T, S <: BoxSDF{T}, N <: RefractiveIndex, C <: Tuple} <: AbstractPlateWaveplate{T, N}
    shape::S
    n::N
    coatings::C
    function RectangularPlateWaveplate(
            shape::S, n::N, coatings::C = ()) where {T <: Real, S <: BoxSDF{T}, N <: RefractiveIndex, C <: Tuple}
        test_refractive_index_function(n)
        return new{T, S, N, C}(shape, n, coatings)
    end
end

"""
    RoundPlateWaveplate{T, S, N, C} <: AbstractPlateWaveplate{T, N}

A round (cylindrical) waveplate with a bulk glass substrate and a retarder coating on its front face.
"""
struct RoundPlateWaveplate{T, S <: PlanoSurfaceSDF{T}, N <: RefractiveIndex, C <: Tuple} <: AbstractPlateWaveplate{T, N}
    shape::S
    n::N
    coatings::C
    function RoundPlateWaveplate(
            shape::S, n::N, coatings::C = ()) where {T <: Real, S <: PlanoSurfaceSDF{T}, N <: RefractiveIndex, C <: Tuple}
        test_refractive_index_function(n)
        return new{T, S, N, C}(shape, n, coatings)
    end
end

# Compatibility accessors
substrate(p::AbstractPlateWaveplate) = p
coating(p::AbstractPlateWaveplate) = get_matching_coating(coatings(p), shape(p), [0.0, 0.0, 0.0], [0.0, -1.0, 0.0])

function RectangularPlateWaveplate(
        width::Real,
        height::Real,
        thickness::Real,
        n::RefractiveIndex,
        retardance::Union{Real, Function};
        fast_axis_angle::Real=0.0,
        back_coating=nothing
)
    T = float(promote_type(typeof(width), typeof(height), typeof(thickness)))
    substrate_shape = BoxSDF(T(width), T(thickness), T(height))
    translate3d!(substrate_shape, [T(0), T(thickness/2), T(0)])
    if fast_axis_angle != 0.0
        yrotate3d!(substrate_shape, -T(fast_axis_angle))
    end
    
    JMat = if retardance isa Function
        λ -> XZBasis(exp(-im * retardance(λ) / 2), 0, 0, exp(im * retardance(λ) / 2))
    else
        XZBasis(exp(-im * retardance / 2), 0, 0, exp(im * retardance / 2))
    end
    coat_wp = JonesCoating(JMat)

    coatings_list = if isnothing(back_coating)
        (:front => coat_wp,)
    else
        (:front => coat_wp, :back => _coating_model(back_coating))
    end

    return RectangularPlateWaveplate(substrate_shape, n, coatings_list)
end

function RoundPlateWaveplate(
        diameter::Real,
        thickness::Real,
        n::RefractiveIndex,
        retardance::Union{Real, Function};
        fast_axis_angle::Real=0.0,
        back_coating=nothing
)
    T = float(promote_type(typeof(diameter), typeof(thickness)))
    substrate_shape = PlanoSurfaceSDF(T(thickness), T(diameter))
    if fast_axis_angle != 0.0
        yrotate3d!(substrate_shape, -T(fast_axis_angle))
    end
    
    JMat = if retardance isa Function
        λ -> XZBasis(exp(-im * retardance(λ) / 2), 0, 0, exp(im * retardance(λ) / 2))
    else
        XZBasis(exp(-im * retardance / 2), 0, 0, exp(im * retardance / 2))
    end
    coat_wp = JonesCoating(JMat)

    coatings_list = if isnothing(back_coating)
        (:front => coat_wp,)
    else
        (:front => coat_wp, :back => _coating_model(back_coating))
    end

    return RoundPlateWaveplate(substrate_shape, n, coatings_list)
end

# Factory mapping to round/rectangular plate waveplates
function Waveplate(diameter::Real, thickness::Real, n::RefractiveIndex, retardance::Union{Real, Function}; fast_axis_angle::Real=0.0, back_coating=nothing)
    return RoundPlateWaveplate(diameter, thickness, n, retardance; fast_axis_angle=fast_axis_angle, back_coating=back_coating)
end

function Waveplate(width::Real, height::Real, thickness::Real, n::RefractiveIndex, retardance::Union{Real, Function}; fast_axis_angle::Real=0.0, back_coating=nothing)
    return RectangularPlateWaveplate(width, height, thickness, n, retardance; fast_axis_angle=fast_axis_angle, back_coating=back_coating)
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
