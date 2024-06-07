struct DoubletLens{T, F<:AbstractShape{T}, B<:AbstractShape{T}} <: AbstractRefractiveOptic{T}
    front::Lens{F}
    back::Lens{B}
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
    # This code is stupid but whatever...
    i_f = intersect3d(dl.front.shape, ray)
    i_b = intersect3d(dl.back.shape, ray)
    # Test if doublet is missed entirely
    if isnothing(i_f) & isnothing(i_b)
        return nothing
    end
    # Test if only front part is hit
    if isnothing(i_b)
        object!(i_f, dl)
        return i_f
    end
    # Test if only back part is hit
    if isnothing(i_f)
        object!(i_b, dl)
        return i_b
    end
    # If both parts are hit, return closest intersection
    if length(i_f) < length(i_b)
        object!(i_f, dl)
        return i_f
    else
        object!(i_b, dl)
        return i_b
    end
    return nothing
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