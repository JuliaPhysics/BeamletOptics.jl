"""
    AbstractLens{T} <: AbstractObject{T}

Abstract super-type for all lens-like objects. Currently this imposes no defined
interface methods.
"""
abstract type AbstractLens{T} <: AbstractObject{T} end

"""
    Lens{T <: AbstractSurface, V <: AbstractSurface, W, F <: Function} <: AbstractLens{W}

Concrete implementation of a lens with two surfaces `front` and `back`.
The lens is characterized by its `normal` vector, its position `pos` and an `edge_thickness`.
The material of this lens currently only affects the refractive index which is passed as
function handle in `ref_index`.
"""
mutable struct Lens{T <: AbstractSurface, V <: AbstractSurface, W, F <: Function} <: AbstractLens{W}
    front::T
    back::V
    normal::Vector{W}
    pos::Vector{W}
    edge_thickness::W
    ref_index::F
end

"""
    center_thickness(l::Lens)

Returns the center thickness of the lens.
"""
center_thickness(l::Lens) = (sag(l.front, clear_semi_diameter(l.front)) +
                            l.edge_thickness +
                            sag(l.back, clear_semi_diameter(l.back)))


"""
    center_of_curvature(l::Lens, s::ConicSurface)

Returns the center of curvature of the surface `s` of the given lens as 3D-coordinate.
"""
function center_of_curvature(l::Lens, s::ConicSurface)
    # get the sign which is dependent on the surface type (front or back) and the sign
    # of the radius. Radius can be negative for concave surfaces.
    _sign = (s === l.front ? -1 : 1)*sign(radius_of_curvature(s))
    # Note this assumes that normal is of length one but this is not asserted at the moment
    return l.pos .+ _sign*(l.edge_thickness/2 + sag(s, clear_semi_diameter(s)) - radius_of_curvature(s)) .* l.normal
end

"""
    sphere(l::Lens, s::ConicSurface)

Returns the equivalent sphere of the conic surface `s` of the lens. This sphere is useful
for intersection testing.
"""
sphere(l::Lens, s::ConicSurface) = Sphere(center_of_curvature(l, s), radius_of_curvature(s))

# AbstractObject Interface defintions --> TODO
