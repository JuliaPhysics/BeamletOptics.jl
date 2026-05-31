abstract type AbstractDoubletRefractiveOptic{
    T,
    F <: AbstractShape{T},
    B <: AbstractShape{T},
    N1 <: RefractiveIndex,
    N2 <: RefractiveIndex
} <: AbstractObjectGroup{T} end

"""
    DoubletLens

Represents a two-component cemented doublet lens with two respective refractive indices `n = n(λ)`.
See also [`SphericalDoubletLens`](@ref).

# Fields

- `front`: front [`Lens`](@ref) component
- `back`: back [`Lens`](@ref) component

# Additional information

!!! warning "Air gap"
    This component type strongly assumes that both lenses are mounted fully flush with respect to each other. 
    Gaps between the components might lead to incorrect results.
"""
struct DoubletLens{T, F<:AbstractShape{T}, B<:AbstractShape{T}, N1<:RefractiveIndex, N2<:RefractiveIndex} <: AbstractDoubletRefractiveOptic{T, F, B, N1, N2}
    front::Lens{T, F, N1}
    back::Lens{T, B, N2}
end

function DoubletLens(front::Lens{T, F, N1}, back::Lens{T, B, N2}) where {T, F, B, N1, N2}
    # Check if they are flush and aligned along the optical axis
    rel_pos = position(back) - position(front)
    opt_axis = orientation(front)[:, 2]
    expected_offset = thickness(front)
    if !isapprox(rel_pos, expected_offset * opt_axis, atol=1e-5)
        @warn "Doublet components are not aligned flush along the optical axis (expected offset $(expected_offset)m along orientation axis, got $(rel_pos)m). This may cause tracing errors at coincident boundaries."
    end
    if !isapprox(orientation(front), orientation(back), atol=1e-5)
        @warn "Doublet components do not have matching orientation. This may cause tracing errors at coincident boundaries."
    end
    return DoubletLens{T, F, B, N1, N2}(front, back)
end

shape_trait_of(::DoubletLens) = MultiShape()

shape(dl::DoubletLens) = (dl.front, dl.back)

objects(dl::DoubletLens) = (dl.front, dl.back)

Base.position(dl::DoubletLens) = position(dl.front)
orientation(dl::DoubletLens) = orientation(dl.front)

thickness(dl::DoubletLens) = thickness(shape(dl.front)) + thickness(shape(dl.back))

"""
    SphericalDoubletLens(r1, r2, r3, l1, l2, d, n1, n2)

Generates a two-component "cemented" doublet lens consisting of two spherical lenses.
For radii sign definition, refer to the [`SphericalLens`](@ref) constructor.

# Arguments

- `r1`: radius of curvature for first surface
- `r2`: radius of curvature for second (cemented) surface
- `r3`: radius of curvature for third surface
- `l1`: first lens thickness
- `l2`: second lens thickness
- `d`: lens diameter
- `n1`: first lens [`RefractiveIndex`](@ref)
- `n2`: second lens [`RefractiveIndex`](@ref)
"""
function SphericalDoubletLens(r1, r2, r3, l1, l2, d, n1, n2)
    # Generate "cemented" front and back spherical lenses
    front = SphericalLens(r1, r2, l1, d, n1)
    back = SphericalLens(r2, r3, l2, d, n2)
    # Move doublet parts into position
    translate3d!(back, [0, thickness(shape(front)), 0])
    return DoubletLens(front, back)
end