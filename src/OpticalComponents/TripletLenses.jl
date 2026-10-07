"""
    TripletLens

Represents a three-component cemented triplet lens with three respective refractive indices `n = n(λ)`.
See also [`SphericalTripletLens`](@ref) and the surface-based constructor
`TripletLens(s1, s2, s3, s4, l1, l2, l3, n1, n2, n3)`.

# Fields

- `front`: front [`Lens`](@ref) component
- `middle`: middle [`Lens`](@ref) component
- `back`: back [`Lens`](@ref) component

# Additional information

!!! info "Clear apertures"
    If the surfaces need different clear apertures (e.g. a steep last surface), use the surface-based
    constructor `TripletLens(s1, s2, s3, s4, l1, l2, l3, n1, n2, n3)` with e.g. `SphericalSurface(r, d)`
    for each surface. [`SphericalTripletLens`](@ref) uses one diameter for all surfaces.

!!! warning "Air gap"
    This component type strongly assumes that all three lenses are mounted fully flush with respect to each other.
    Gaps between the components might lead to incorrect results.

!!! warning "Total internal reflection"
    Total internal reflection at a cemented interface between two elements is not modeled correctly.
"""
struct TripletLens{T, F<:AbstractShape{T}, M<:AbstractShape{T}, B<:AbstractShape{T},
                   N1<:RefractiveIndex, N2<:RefractiveIndex, N3<:RefractiveIndex} <: AbstractObject{T}
    front::Lens{T, F, N1}
    middle::Lens{T, M, N2}
    back::Lens{T, B, N3}
end

shape_trait_of(::TripletLens) = MultiShape()

shape(tl::TripletLens) = (tl.front, tl.middle, tl.back)

Base.position(tl::TripletLens) = position(tl.front)
orientation(tl::TripletLens) = orientation(tl.front)

thickness(tl::TripletLens) = thickness(shape(tl.front)) + thickness(shape(tl.middle)) + thickness(shape(tl.back))

"""
    TripletLens(s1, s2, s3, s4, l1, l2, l3, n1, n2, n3)

Generates a three-component "cemented" [`TripletLens`](@ref) from four surface specifications, e.g.
[`SphericalSurface`](@ref), [`EvenAsphericalSurface`](@ref) or [`CircularFlatSurface`](@ref).
The front, middle and back elements are `Lens(s1, s2, l1, n1)`, `Lens(s2, s3, l2, n2)` and
`Lens(s3, s4, l3, n3)`, refer to the surface-based [`Lens`](@ref) constructor for the construction
of the elements. The front vertex lies at the origin and the middle and back elements are translated
along +y by `l1` and `l1 + l2`, so that neighbouring elements share the cemented surfaces `s2` and `s3`.

# Arguments

- `s1`: first surface
- `s2`: second (first cemented) surface
- `s3`: third (second cemented) surface
- `s4`: fourth surface
- `l1`: first lens center thickness in m
- `l2`: second lens center thickness in m
- `l3`: third lens center thickness in m
- `n1`: first lens [`RefractiveIndex`](@ref)
- `n2`: second lens [`RefractiveIndex`](@ref)
- `n3`: third lens [`RefractiveIndex`](@ref)

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
    Total internal reflection at a cemented interface between two elements is not modeled correctly.
"""
function TripletLens(
        s1::AbstractRotationallySymmetricSurface,
        s2::AbstractRotationallySymmetricSurface,
        s3::AbstractRotationallySymmetricSurface,
        s4::AbstractRotationallySymmetricSurface,
        l1::Real,
        l2::Real,
        l3::Real,
        n1::RefractiveIndex,
        n2::RefractiveIndex,
        n3::RefractiveIndex)
    # Generate "cemented" front, middle and back lenses that share the surfaces s2 and s3
    front = Lens(s1, s2, l1, n1)
    middle = Lens(s2, s3, l2, n2)
    back = Lens(s3, s4, l3, n3)
    # Move triplet parts into position
    translate3d!(middle, [0, thickness(shape(front)), 0])
    translate3d!(back, [0, thickness(shape(front)) + thickness(shape(middle)), 0])
    return TripletLens(front, middle, back)
end

"""
    SphericalTripletLens(r1, r2, r3, r4, l1, l2, l3, d, n1, n2, n3)

Generates a three-component "cemented" triplet lens consisting of three spherical lenses.
The middle and back lens are translated along +y so that they sit flush.
For radii sign definition, refer to the [`SphericalLens`](@ref) constructor.

# Arguments

- `r1`: radius of curvature for first surface
- `r2`: radius of curvature for second (first cemented) surface
- `r3`: radius of curvature for third (second cemented) surface
- `r4`: radius of curvature for fourth surface
- `l1`: first lens thickness
- `l2`: second lens thickness
- `l3`: third lens thickness
- `d`: lens diameter
- `n1`: first lens [`RefractiveIndex`](@ref)
- `n2`: second lens [`RefractiveIndex`](@ref)
- `n3`: third lens [`RefractiveIndex`](@ref)

# Additional information

`Inf` gives a plano surface.
"""
function SphericalTripletLens(r1, r2, r3, r4, l1, l2, l3, d, n1, n2, n3)
    # Generate "cemented" front, middle and back spherical lenses
    front = SphericalLens(r1, r2, l1, d, n1)
    middle = SphericalLens(r2, r3, l2, d, n2)
    back = SphericalLens(r3, r4, l3, d, n3)
    # Move triplet parts into position
    translate3d!(middle, [0, thickness(shape(front)), 0])
    translate3d!(back, [0, thickness(shape(front)) + thickness(shape(middle)), 0])
    return TripletLens(front, middle, back)
end

function interact3d(system::AbstractSystem, tl::TripletLens, beam::Beam{T, R}, ray::R) where {T <: Real, R <: AbstractRay{T}}
    # Interaction logic: front/back always hint the middle element. The middle element
    # hints the neighbour that lies ahead along the ray's direction of travel.
    hit = shape(intersection(ray))
    if hit === shape(tl.front)
        i = interact3d(system, tl.front, beam, ray)
        next = tl.middle
    elseif hit === shape(tl.back)
        i = interact3d(system, tl.back, beam, ray)
        next = tl.middle
    elseif hit === shape(tl.middle)
        i = interact3d(system, tl.middle, beam, ray)
        isnothing(i) && return nothing
        s = dot(direction(i.ray), orientation(tl)[:, 2])
        next = s ≥ 0 ? tl.back : tl.front
    else
        error("TripletLens: intersected shape is not part of this lens")
    end
    isnothing(i) && return nothing
    return BeamInteraction(Hint(tl, shape(next)), i.ray)
end
