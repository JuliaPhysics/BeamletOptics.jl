"""
    AbstractLens{T} <: AbstractObject{T}

Abstract super-type for all lens-like objects.

# Implementation

- `normal`: Return a vector normal to the center plane of the lens
"""
abstract type AbstractLens{T} <: AbstractObject{T} end

"""
    normal(l::AbstractLens)

Default implementation for the lens normal. Returns the 3rd vector of the lens orientation.
For a non-tilted lens this is identical to the vector ``\\boldsymbol{e}_z``.
"""
normal(l::AbstractLens) = @view orientation(l)[:,3]

"""
    Lens{T <: AbstractSurface, V <: AbstractSurface, W, F <: Function} <: AbstractLens{W}

Concrete implementation of a lens with two surfaces `front` and `back`.
The lens is characterized by its `normal` vector, its position `pos` and an `edge_thickness`.
The material of this lens currently only affects the refractive index which is passed as
function handle in `ref_index`.
"""
mutable struct SingletLens{T <: AbstractSurface, V <: AbstractSurface, M<:AbstractMatrix{<:AbstractFloat}, W, F <: Function} <: AbstractLens{W}
    front::T
    back::V
    pos::Vector{W}
    dir::M
    edge_thickness::W
    ref_index::F
end

function SingletLens(front::T, back::V, edge_thickness::W=0.0, ref_index=x->1.5) where {T<:AbstractSurface, V<:AbstractSurface, W}
    return SingletLens(
        front,
        back,
        zeros(W, 3),
        Matrix{W}(I, 3, 3),
        edge_thickness,
        ref_index
    )
end

"""
    center_thickness(l::Lens)

Returns the center thickness of the lens.
"""
center_thickness(l::SingletLens) = (sag(l.front, clear_semi_diameter(l.front)) +
                            l.edge_thickness +
                            sag(l.back, clear_semi_diameter(l.back)))


"""
    center_of_curvature(l::Lens, s::ConicSurface)

Returns the center of curvature of the surface `s` of the given lens as 3D-coordinate.
"""
function center_of_curvature(l::SingletLens, s::Union{CylinderSurface,ConicSurface})
    # get the sign which is dependent on the surface type (front or back) and the sign
    # of the radius. Radius can be negative for concave surfaces.
    _sign = (s === l.front ? -1 : 1)*sign(radius_of_curvature(s))
    # Note this assumes that normal is of length one but this is not asserted at the moment
    return l.pos .+ _sign*(l.edge_thickness/2 + sag(s, clear_semi_diameter(s)) - radius_of_curvature(s)) .* normal(l)
end

"""
    sphere(l::Lens, s::ConicSurface)

Returns the equivalent sphere of the conic surface `s` of the lens. This sphere is useful
for intersection testing.
"""
sphere(l::SingletLens, s::Union{CylinderSurface,ConicSurface}) = Sphere(center_of_curvature(l, s), radius_of_curvature(s))
