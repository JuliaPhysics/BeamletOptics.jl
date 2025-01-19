abstract type AbstractDoubletRefractiveOptic{T, F <: AbstractShape{T}, B <: AbstractShape{T}, F1 <: Function, F2 <: Function} <: AbstractRefractiveOptic{T, F, F1} end

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
struct DoubletLens{T, F<:AbstractShape{T}, B<:AbstractShape{T}, F1, F2} <: AbstractDoubletRefractiveOptic{T, F, B, F1, F2}
    front::Lens{T, F, F1}
    back::Lens{T, B, F2}
end

position(dl::DoubletLens) = position(dl.front)
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
- `n1`: first lens refractive index fct.
- `n1`: second lens refractive index fct.
"""
function SphericalDoubletLens(r1, r2, r3, l1, l2, d, n1, n2)
    # Generate "cemented" front and back spherical lenses
    front = SphericalLens(r1, r2, l1, d, n1)
    back = SphericalLens(r2, r3, l2, d, n2)
    # Move doublet parts into position
    translate3d!(back, [0, thickness(shape(front)), 0])
    return DoubletLens(front, back)
end

function intersect3d(dl::DoubletLens, ray::AbstractRay)
    i_f = intersect3d(dl.front.shape, ray)
    i_b = intersect3d(dl.back.shape, ray)
    # Determine which intersections are valid and handle accordingly
    if isnothing(i_f)
        return isnothing(i_b) ? nothing : (object!(i_b, dl); i_b)
    elseif isnothing(i_b)
        return (object!(i_f, dl); i_f)
    else
        closest_intersection = length(i_f) < length(i_b) ? i_f : i_b
        object!(closest_intersection, dl)
        return closest_intersection
    end
end

function interact3d(system::AbstractSystem, dl::DoubletLens, beam::Beam{T, R}, ray::R) where {T <: Real, R <: Ray{T}}
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

function translate3d!(dl::DoubletLens, offset)
    translate3d!(dl.front, offset)
    translate3d!(dl.back, offset)
    return nothing
end

function rotate3d!(dl::DoubletLens, axis, θ)
    # Pivot elements around front lens center point
    rotate3d!(dl.front, axis, θ)
    rotate3d!(dl.back, axis, θ)
    v = position(dl.back) - position(dl.front)
    R = rotate3d(axis, θ)
    v = (R * v) - v
    translate3d!(dl.back, v)
    return nothing
end

function render_object!(ax, dl::DoubletLens)
    render_object!(ax, dl.front)
    render_object!(ax, dl.back)
    return nothing
end
