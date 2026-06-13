struct PolarizingCubeBeamsplitterShape{T} <: AbstractShape{T} end

"""
    PolarizingCubeBeamsplitter <: AbstractBeamsplitter

A cuboid polarizing beamsplitter where the splitting interaction occurs between two [`RightAnglePrism`](@ref)s.
For more information refer to the [`AbstractPlateBeamsplitter`](@ref) docs.

# Fields

- `front`: the forward facing substrate, represented by a coated prism ([`CoatedRefractive`](@ref) wrapping a [`RightAnglePrism`](@ref))
- `back`: the backward facing substrate, represented by a [`RightAnglePrism`](@ref)
"""
struct PolarizingCubeBeamsplitter{T, F, B} <: AbstractObjectGroup{T}
    front::F
    back::B
end

shape_trait_of(::PolarizingCubeBeamsplitter) = MultiShape()

shape(cbs::PolarizingCubeBeamsplitter) = (cbs.front, cbs.back)

objects(cbs::PolarizingCubeBeamsplitter) = (cbs.front, cbs.back)

refractive_index(cbs::PolarizingCubeBeamsplitter, λ::Real) = refractive_index(cbs.front, λ)

"""
    PolarizingCubeBeamsplitter(leg_length, n)

Creates a [`PolarizingCubeBeamsplitter`](@ref). The cuboid is centered at the origin. The splitter 
coating is orientated at a 45° angle with respect to the y-axis.

# Inputs

- `leg_length`: the x-, y- and z-edge length in [m]
- `n`: the [`RefractiveIndex`](@ref) of the front and back prism
"""
function PolarizingCubeBeamsplitter(
        leg_length::Real,
        n::RefractiveIndex
    )
    front_prism = RightAnglePrism(leg_length, leg_length, n)
    back = RightAnglePrism(leg_length, leg_length, n)
    zrotate3d!(back, deg2rad(180))

    coat_bs = SimpleBeamsplitterCoating(
        1.0, 0.0, # rs, rp
        0.0, 1.0  # ts, tp
    )

    T = typeof(leg_length)
    coatings = (:hypotenuse => coat_bs,)
    front = CoatedRefractive{T, typeof(front_prism), typeof(coatings)}(front_prism, coatings)

    return PolarizingCubeBeamsplitter{T, typeof(front), typeof(back)}(front, back)
end
