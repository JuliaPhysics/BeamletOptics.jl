"""Placeholder type for the  shape of a [`CubeBeamsplitter`](@ref)"""
struct CubeBeamsplitterShape{T} <: AbstractShape{T} end

"""
    CubeBeamsplitter <: AbstractBeamsplitter

A cuboid beamsplitter where the splitting interaction occurs between two [`RightAnglePrism`](@ref)s.
For more information refer to the [`AbstractPlateBeamsplitter`](@ref) docs.

# Fields

- `front`: the forward facing substrate, represented by a coated prism (coated [`RightAnglePrism`](@ref))
- `back`: the backward facing substrate, represented by a [`RightAnglePrism`](@ref)

# Additional information

!!! info "Hints and interaction logic"
    In order to model gap-free beam propagation, the `interact3d` model relies heavily on the [`Hint`](@ref)-API.
    If the `front` or `back` substrate is hit, the `Hint` will ensure that the beam intersects the `coating`.
"""
struct CubeBeamsplitter{T, F, B} <: AbstractObjectGroup{T}
    front::F
    back::B
end

shape_trait_of(::CubeBeamsplitter) = MultiShape()

shape(cbs::CubeBeamsplitter) = (cbs.front, cbs.back)

objects(cbs::CubeBeamsplitter) = (cbs.front, cbs.back)

refractive_index(cbs::CubeBeamsplitter, λ::Real) = refractive_index(cbs.front, λ)

"""
    CubeBeamsplitter(leg_length, n; reflectance=0.5)

Creates a [`CubeBeamsplitter`](@ref). The cuboid is centered at the origin. The splitter 
coating is orientated at a 45° angle with respect to the y-axis.

# Inputs

- `leg_length`: the x-, y- and z-edge length in [m]
- `n`: the [`RefractiveIndex`](@ref) of the front and back prism

# Keywords 

- `reflectance`: defines the splitting ratio in [-], i.e. R = 0 ... 1.0
"""
function CubeBeamsplitter(
        leg_length::Real,
        n::RefractiveIndex;
        reflectance::Real=0.5
    )
    T = float(typeof(leg_length))
    front_prism = RightAnglePrism(T(leg_length), T(leg_length), n)
    back = RightAnglePrism(T(leg_length), T(leg_length), n)
    zrotate3d!(back, deg2rad(180))

    coat_bs = SimpleBeamsplitterCoating(
        sqrt(reflectance), sqrt(reflectance),
        sqrt(1.0 - reflectance), sqrt(1.0 - reflectance)
    )

    front = with_coatings(front_prism, :hypotenuse => coat_bs)

    return CubeBeamsplitter{T, typeof(front), typeof(back)}(front, back)
end