"""
    SphericalLens(r1, r2, l, d=1inch, n=ConstantRefractiveIndex(1.5))

Creates a spherical [`Lens`](@ref) based on:

- `r1`: front radius
- `r2`: back radius
- `l`: lens thickness
- `d`: lens diameter, default is one inch
- `n`: [`AbstractRefractiveIndex`](@ref) of the lens material

# Notes

!!! info "Radius of curvature (ROC) sign"
    The ROC is defined to be positive if the center is to the right of the surface. Otherwise it is negative.

!!! info "Thin lenses"
    If `l` is set to zero, a [`ThinLens`](@ref) will be created. However, note that the actual lens thickness will be different from zero.
"""
function SphericalLens(r1::Real, r2::Real, l::Real, d::Real = 1inch, n::AbstractRefractiveIndex = ConstantRefractiveIndex(1.5))
    # Test for thin lens
    if iszero(l)
        return ThinLens(r1, r2, d, n)
    end
    # Create lens
    return Lens(
        SphericalSurface(r1, d),
        SphericalSurface(r2, d),
        l,
        n
    )
end


"""
    ThinLens(R1::Real, R2::Real, d::Real, n::AbstractRefractiveIndex)

Directly creates an ideal spherical thin [`Lens`](@ref) with radii of curvature `R1` and `R2` and diameter `d`
and [`AbstractRefractiveIndex`](@ref) `n`.
"""
function ThinLens(R1::Real, R2::Real, d::Real, n::AbstractRefractiveIndex)
    shape = ThinLensSDF(R1, R2, d)
    return Lens(shape, n)
end
