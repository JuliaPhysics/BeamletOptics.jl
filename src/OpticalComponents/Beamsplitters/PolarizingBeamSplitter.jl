# Polarizing Beam Splitter (PBS) Components
# Reuses the generic Coating system under the hood with a SimpleBeamsplitterCoating.

"""
    PolarizingBeamSplitter(shape)
    PolarizingBeamSplitter(width, height)
    PolarizingBeamSplitter(diameter)

Ideal polarizing plate that separates an incoming `PolarizedRay` into
transmitted and reflected beams. The device is represented by a zero-thickness
surface whose orientation sets the splitting axes:

- local `x`-axis → transmitted (Ex) polarization component
- local `z`-axis → reflected (Ez) polarization component

Rotate the underlying shape to align these axes with the desired global
polarization directions. Incoming rays are assumed to strike the plate from the
positive local `y` direction.
"""
function PolarizingBeamSplitter(shape::AbstractShape{T}) where {T}
    J_t = XZBasis(one(T), zero(T), zero(T), zero(T))
    J_r = XZBasis(zero(T), zero(T), zero(T), -one(T))
    model = JonesCoating(J_t, J_r; behavior = Splitting())
    return Coating(shape, model)
end

function PolarizingBeamSplitter(width::Real, height::Real)
    return PolarizingBeamSplitter(RectangularFlatMesh(width, height))
end

function PolarizingBeamSplitter(diameter::Real)
    return PolarizingBeamSplitter(CircularFlatMesh(diameter / 2))
end

#=
Plate Polarizing Beamsplitters
=#

"""
    RectangularPolarizingPlateBeamsplitter{T, M} <: AbstractPlateBeamsplitter{T}

A plate beamsplitter with a rectangular substrate and an ideal polarizing
splitting coating.
"""
struct RectangularPolarizingPlateBeamsplitter{T, M} <: AbstractPlateBeamsplitter{T}
    substrate::Prism{T, BoxSDF{T}}
    coating::Coating{T, Mesh{T}, M}
end

"""
    RectangularPolarizingPlateBeamsplitter(width, height, thickness, n)
"""
function RectangularPolarizingPlateBeamsplitter(
        width::Real,
        height::Real,
        thickness::Real,
        n::RefractiveIndex
)
    T = float(promote_type(typeof(width), typeof(height), typeof(thickness)))
    substrate_shape = BoxSDF(T(width), T(thickness), T(height))
    substrate = Prism(substrate_shape, n)
    translate3d!(substrate, [T(0), T(thickness / 2), T(0)])
    coating = PolarizingBeamSplitter(T(width), T(height))
    zrotate3d!(coating, T(π))
    M = typeof(coating.model)
    return RectangularPolarizingPlateBeamsplitter{T, M}(substrate, coating)
end

"""
    RoundPolarizingPlateBeamsplitter{T, M} <: AbstractPlateBeamsplitter{T}

A plate beamsplitter with a cylindrical substrate and an ideal polarizing
splitting coating.
"""
struct RoundPolarizingPlateBeamsplitter{T, M} <: AbstractPlateBeamsplitter{T}
    substrate::Prism{T, PlanoSurfaceSDF{T}}
    coating::Coating{T, Mesh{T}, M}
end

"""
    RoundPolarizingPlateBeamsplitter(diameter, thickness, n)
"""
function RoundPolarizingPlateBeamsplitter(
        diameter::Real,
        thickness::Real,
        n::RefractiveIndex
)
    T = float(promote_type(typeof(diameter), typeof(thickness)))
    substrate_shape = PlanoSurfaceSDF(T(thickness), T(diameter))
    substrate = Prism(substrate_shape, n)
    coating = PolarizingBeamSplitter(CircularFlatMesh(T(diameter) / 2))
    M = typeof(coating.model)
    return RoundPolarizingPlateBeamsplitter{T, M}(substrate, coating)
end
