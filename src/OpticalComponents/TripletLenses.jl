"""
    TripletLens

Represents a three-component cemented triplet lens with three respective refractive indices `n = n(λ)`.
See also [`SphericalTripletLens`](@ref).

# Fields

- `front`: front [`Lens`](@ref) component
- `middle`: middle [`Lens`](@ref) component
- `back`: back [`Lens`](@ref) component

# Additional information

!!! info "Clear apertures"
    If the surfaces need different clear apertures (e.g. a steep last surface), build the elements via
    `Lens(SphericalSurface(r1, d1), SphericalSurface(r2, d2), l, n)` and pass them to `TripletLens`.
    [`SphericalTripletLens`](@ref) uses one diameter for all surfaces.

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
