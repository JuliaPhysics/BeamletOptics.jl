"""
    DoubletLens

Represents a two-component cemented doublet lens with two respective refractive indices `n = n(λ)`.
See also [`SphericalDoubletLens`](@ref) and the surface-based constructor
`DoubletLens(s1, s2, s3, l1, l2, n1, n2)`.

# Fields

- `front`: front [`Lens`](@ref) component
- `back`: back [`Lens`](@ref) component

# Additional information

!!! warning "Air gap"
    This component type strongly assumes that both lenses are mounted fully flush with respect to each other. 
    Gaps between the components might lead to incorrect results.
"""
struct DoubletLens{T, F<:AbstractShape{T}, B<:AbstractShape{T}, N1<:RefractiveIndex, N2<:RefractiveIndex} <: AbstractObject{T}
    front::Lens{T, F, N1}
    back::Lens{T, B, N2}
end

shape_trait_of(::DoubletLens) = MultiShape()

shape(dl::DoubletLens) = (dl.front, dl.back)

Base.position(dl::DoubletLens) = position(dl.front)
orientation(dl::DoubletLens) = orientation(dl.front)

thickness(dl::DoubletLens) = thickness(shape(dl.front)) + thickness(shape(dl.back))

"""
    DoubletLens(s1, s2, s3, l1, l2, n1, n2)

Generates a two-component "cemented" [`DoubletLens`](@ref) from three surface specifications, e.g.
[`SphericalSurface`](@ref), [`EvenAsphericalSurface`](@ref) or [`CircularFlatSurface`](@ref).
The front element is `Lens(s1, s2, l1, n1)` and the back element is `Lens(s2, s3, l2, n2)`,
refer to the surface-based [`Lens`](@ref) constructor for the construction of the elements.
The front vertex lies at the origin and the back element is translated along +y by `l1`, so that
both elements share the cemented surface `s2`.

# Arguments

- `s1`: first surface
- `s2`: second (cemented) surface
- `s3`: third surface
- `l1`: first lens center thickness in m
- `l2`: second lens center thickness in m
- `n1`: first lens [`RefractiveIndex`](@ref)
- `n2`: second lens [`RefractiveIndex`](@ref)

# Additional information

!!! info "Radius of curvature (ROC) sign definition"
    The ROC is defined to be positive if the center is to the right of the surface, i.e. at +y. Otherwise it is negative.

!!! info "Clear apertures"
    Each surface has its own clear aperture and mechanical diameter, which are leveled per element as for a single [`Lens`](@ref).

!!! warning "Supported surfaces"
    Only rotationally symmetric surfaces are supported, cylindrical surfaces are not. The limitations of the
    surface-based [`Lens`](@ref) constructor apply to each element: a meniscus element, whose center thickness
    does not exceed the sagitta of its convex surface, must consist of spherical surfaces.

!!! warning "Total internal reflection"
    Total internal reflection at the cemented interface between the two elements is not modeled correctly.
"""
function DoubletLens(
        s1::AbstractRotationallySymmetricSurface,
        s2::AbstractRotationallySymmetricSurface,
        s3::AbstractRotationallySymmetricSurface,
        l1::Real,
        l2::Real,
        n1::RefractiveIndex,
        n2::RefractiveIndex)
    # Generate "cemented" front and back lenses that share the surface s2
    front = Lens(s1, s2, l1, n1)
    back = Lens(s2, s3, l2, n2)
    # Move doublet parts into position
    translate3d!(back, [0, thickness(shape(front)), 0])
    return DoubletLens(front, back)
end

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
- `n1`: second lens [`RefractiveIndex`](@ref)
"""
function SphericalDoubletLens(r1, r2, r3, l1, l2, d, n1, n2)
    # Generate "cemented" front and back spherical lenses
    front = SphericalLens(r1, r2, l1, d, n1)
    back = SphericalLens(r2, r3, l2, d, n2)
    # Move doublet parts into position
    translate3d!(back, [0, thickness(shape(front)), 0])
    return DoubletLens(front, back)
end

function interact3d(system::AbstractSystem, dl::DoubletLens, beam::Beam{T, R}, ray::R) where {T <: Real, R <: AbstractRay{T}}
    # Interaction logic: if front is hit, hint to back and vice versa
    if shape(intersection(ray)) === shape(dl.front)
        i = interact3d(system, dl.front, beam, ray)
        hint = Hint(dl, dl.back.shape)
    elseif shape(intersection(ray)) === shape(dl.back)
        i = interact3d(system, dl.back, beam, ray)
        hint = Hint(dl, dl.front.shape)
    end
    return BeamInteraction(hint, i.ray)
end