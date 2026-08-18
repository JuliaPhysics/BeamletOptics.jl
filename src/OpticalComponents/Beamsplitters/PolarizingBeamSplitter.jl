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
    T = float(promote_type(typeof(width), typeof(height)))
    return PolarizingBeamSplitter(RectangularFlatMesh(T(width), T(height)))
end

function PolarizingBeamSplitter(diameter::Real)
    T = float(typeof(diameter))
    return PolarizingBeamSplitter(CircularFlatMesh(T(diameter) / 2))
end

#=
Plate Polarizing Beamsplitters
=#

"""
    RectangularPolarizingPlateBeamsplitter{T, S, N, C} <: AbstractPlateBeamsplitter{T, N}

A plate beamsplitter with a rectangular substrate and an ideal polarizing splitting coating attached to its front face.
"""
struct RectangularPolarizingPlateBeamsplitter{T, S <: BoxSDF{T}, N <: RefractiveIndex, C <: Tuple} <: AbstractPlateBeamsplitter{T, N}
    shape::S
    n::N
    coatings::C
    function RectangularPolarizingPlateBeamsplitter(
            shape::S, n::N, coatings::C = ()) where {T <: Real, S <: BoxSDF{T}, N <: RefractiveIndex, C <: Tuple}
        test_refractive_index_function(n)
        return new{T, S, N, C}(shape, n, coatings)
    end
end

"""
    RoundPolarizingPlateBeamsplitter{T, S, N, C} <: AbstractPlateBeamsplitter{T, N}

A plate beamsplitter with a cylindrical substrate and an ideal polarizing splitting coating attached to its front face.
"""
struct RoundPolarizingPlateBeamsplitter{T, S <: PlanoSurfaceSDF{T}, N <: RefractiveIndex, C <: Tuple} <: AbstractPlateBeamsplitter{T, N}
    shape::S
    n::N
    coatings::C
    function RoundPolarizingPlateBeamsplitter(
            shape::S, n::N, coatings::C = ()) where {T <: Real, S <: PlanoSurfaceSDF{T}, N <: RefractiveIndex, C <: Tuple}
        test_refractive_index_function(n)
        return new{T, S, N, C}(shape, n, coatings)
    end
end

_attach_coatings(pbs::RectangularPolarizingPlateBeamsplitter, c_tuple; deepcopy_shape::Bool = false) =
    RectangularPolarizingPlateBeamsplitter(deepcopy_shape ? deepcopy(pbs.shape) : pbs.shape, pbs.n, c_tuple)

_attach_coatings(pbs::RoundPolarizingPlateBeamsplitter, c_tuple; deepcopy_shape::Bool = false) =
    RoundPolarizingPlateBeamsplitter(deepcopy_shape ? deepcopy(pbs.shape) : pbs.shape, pbs.n, c_tuple)

"""
    RectangularPolarizingPlateBeamsplitter(width, height, thickness, n; back_coating=nothing)

A plate beamsplitter with a rectangular substrate and an ideal polarizing splitting coating attached to its front face.
"""
function RectangularPolarizingPlateBeamsplitter(
        width::Real,
        height::Real,
        thickness::Real,
        n::RefractiveIndex;
        back_coating=nothing
)
    T = float(promote_type(typeof(width), typeof(height), typeof(thickness)))
    substrate_shape = BoxSDF(T(width), T(thickness), T(height))
    translate3d!(substrate_shape, [T(0), T(thickness / 2), T(0)])

    J_t = XZBasis(one(T), zero(T), zero(T), zero(T))
    J_r = XZBasis(zero(T), zero(T), zero(T), -one(T))
    coat_pbs = JonesCoating(J_t, J_r; behavior = Splitting())

    coatings_list = if isnothing(back_coating)
        (:front => coat_pbs,)
    else
        (:front => coat_pbs, :back => _coating_model(back_coating))
    end

    return RectangularPolarizingPlateBeamsplitter(substrate_shape, n, coatings_list)
end

"""
    RoundPolarizingPlateBeamsplitter(diameter, thickness, n; back_coating=nothing)

A plate beamsplitter with a cylindrical substrate and an ideal polarizing splitting coating attached to its front face.
"""
function RoundPolarizingPlateBeamsplitter(
        diameter::Real,
        thickness::Real,
        n::RefractiveIndex;
        back_coating=nothing
)
    T = float(promote_type(typeof(diameter), typeof(thickness)))
    substrate_shape = PlanoSurfaceSDF(T(thickness), T(diameter))

    J_t = XZBasis(one(T), zero(T), zero(T), zero(T))
    J_r = XZBasis(zero(T), zero(T), zero(T), -one(T))
    coat_pbs = JonesCoating(J_t, J_r; behavior = Splitting())

    coatings_list = if isnothing(back_coating)
        (:front => coat_pbs,)
    else
        (:front => coat_pbs, :back => _coating_model(back_coating))
    end

    return RoundPolarizingPlateBeamsplitter(substrate_shape, n, coatings_list)
end


