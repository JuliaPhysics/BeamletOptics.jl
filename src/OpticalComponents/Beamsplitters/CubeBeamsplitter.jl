"""Placeholder type for the  shape of a [`CubeBeamsplitter`](@ref)"""
struct CubeBeamsplitterShape{T} <: AbstractShape{T} end

"""
    CubeBeamsplitter <: AbstractBeamsplitter

A cuboid beamsplitter where the splitting interaction occurs between two [`RightAnglePrism`](@ref)s.
For more information refer to the [`AbstractBeamsplitter`](@ref) docs.

# Fields

- `front`: the forward facing substrate, represented by a [`RightAnglePrism`](@ref)
- `back`: the backward facing substrate, represented by a [`RightAnglePrism`](@ref)
- `coating`: the thin layer at the hypotenuse that represents the splitting interface, a rectangular
  [`ThinBeamsplitter`](@ref) or [`CoatingLayer`](@ref), e.g. of a [`PolarizingCubeBeamsplitter`](@ref)

# Additional information

!!! info "Interaction logic"
    The prisms and the coating meet at the hypotenuse. A ray that reaches it hits all three at once, and the
    coating splits it, the transmitted beam passing from one prism into the other, see [`transition`](@ref).
    Moved apart, e.g. via `translate3d!(cbs.front, ...)`, the parts leave an air gap that the trace sees.
"""
struct CubeBeamsplitter{T, C <: AbstractObject{T}} <: AbstractBeamsplitter{T}
    front::Prism{T, RightAnglePrismSDF{T}}
    back::Prism{T, RightAnglePrismSDF{T}}
    coating::C
end

shape_trait_of(::CubeBeamsplitter) = MultiShape()

shape(cbs::CubeBeamsplitter) = (cbs.front, cbs.back, cbs.coating)

refractive_index(cbs::CubeBeamsplitter, λ::Real) = refractive_index(cbs.front, λ)

"""
    CubeBeamsplitter(leg_length, n; reflectance=0.5)

Creates a [`CubeBeamsplitter`](@ref). The cuboid is centered at the origin. The splitter 
coating is orientated at a 45° angle with respect to the y-axis.

# Inputs

- `leg_length`: the x-, y- and z-edge length in [m]
- `n`: the [`AbstractRefractiveIndex`](@ref) of the front and back prism

# Keywords 

- `reflectance`: defines the power reflectance R (dimensionless) of the coating, in the open interval (0, 1)
"""
function CubeBeamsplitter(
        leg_length::Real,
        n::AbstractRefractiveIndex;
        reflectance::Real=0.5
    )
    return cube_beamsplitter(leg_length, n, ThinBeamsplitter(√2 * leg_length, leg_length; reflectance))
end

# A cube of two right-angle prisms of the index `n` with the thin `layer` at the hypotenuse, which is
# built in the x-z plane of the size √2 leg_length × leg_length
function cube_beamsplitter(leg_length::Real, n::AbstractRefractiveIndex, layer::AbstractObject)
    front = RightAnglePrism(leg_length, leg_length, n)
    back = RightAnglePrism(leg_length, leg_length, n)
    zrotate3d!(back, deg2rad(180))
    zrotate3d!(layer, deg2rad(180 - 45))
    set_new_origin3d!(shape(layer))
    return CubeBeamsplitter(front, back, layer)
end

"""
    PolarizingCubeBeamsplitter(leg_length, n; Rs = 1, Rp = 0)

Creates a polarizing [`CubeBeamsplitter`](@ref): the coating at its hypotenuse reflects the
s-component of the field and transmits the p-component, see [`PolarizingCoating`](@ref). Placed as
built, a beam along +y is split in the x-y plane: the field along x (p) is transmitted along +y, the
field along z (s) is reflected along -x. `Rs` and `Rp` are the power reflectances of both components,
e.g. `Rp = 0.01` for a cube with an extinction ratio of 1:100 in transmission.

# Inputs

- `leg_length`: the x-, y- and z-edge length in [m]
- `n`: the [`AbstractRefractiveIndex`](@ref) of the front and back prism
"""
function PolarizingCubeBeamsplitter(leg_length::Real, n::AbstractRefractiveIndex; Rs::Real = 1, Rp::Real = 0)
    layer = CoatingLayer(RectangularFlatMesh(√2 * leg_length, leg_length), PolarizingCoating(Rs, Rp))
    return cube_beamsplitter(leg_length, n, layer)
end
