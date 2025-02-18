"""
    RectangularCompensatorPlate(width, height, thickness, n)

Creates a compensator plate (modeled as a [`Prism`](@ref)) that can be used to remove
parallel beam offsets created by e.g. the [`RectangularPlateBeamsplitter`](@ref).
The compensator is aligned with the positive y-axis. The first surface lies at the origin.

# Inputs

- `width`: compensator width along the x-axis in [m]
- `height`: compensator height along the z-axis in [m]
- `thickness`: compensator thickness along the y-axis in [m]
- `n`: the [`RefractiveIndex`](@ref) of the substrate
"""
function RectangularCompensatorPlate(width::W, height::H, thickness::T, n::RefractiveIndex) where {W<:Real,H<:Real,T<:Real}
    shape = CuboidMesh(width, thickness, height)
    translate3d!(shape, [
        -width/2,
        0,
        -height/2,
    ])
    set_new_origin3d!(shape)
    return Prism(shape, n)
end

RectangularCompensatorPlate(w::Real, h::Real, t::Real, n::Real) = RectangularCompensatorPlate(w, h, t, λ->n)