"""
    CubeBeamsplitter{T} <: AbstractObject{T}

A cuboid beamsplitter where the splitting interaction occurs between two [`RightAnglePrism`](@ref)s.
For more information refer to the [`AbstractPlateBeamsplitter`](@ref) docs.

# Fields

- `front`: the forward facing substrate, represented by a [`RightAnglePrism`](@ref)
- `back`: the backward facing substrate, represented by a [`RightAnglePrism`](@ref)
- `coating`: a rectangular [`ThinBeamsplitter`](@ref) that represents the splitting interface
"""
struct CubeBeamsplitter{T} <: AbstractObject{T}
    front::Prism{T, RightAnglePrismSDF{T}}
    back::Prism{T, RightAnglePrismSDF{T}}
    coating::ThinBeamsplitter{T, Mesh{T}}
end

shape_trait_of(::CubeBeamsplitter) = MultiShape()

shape(cbs::CubeBeamsplitter) = (cbs.front, cbs.back, cbs.coating)

Base.position(cbs::CubeBeamsplitter) = position(cbs.coating)
orientation(cbs::CubeBeamsplitter) = orientation(cbs.front)

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
    front = RightAnglePrism(leg_length, leg_length, n)
    back = RightAnglePrism(leg_length, leg_length, n)
    bs = ThinBeamsplitter(√2*leg_length, leg_length; reflectance)
    zrotate3d!(back, deg2rad(180))
    zrotate3d!(bs, deg2rad(180-45))
    set_new_origin3d!(shape(bs))
    return CubeBeamsplitter(front, back, bs)
end