"""
    Prism{T, S <: AbstractShape{T}, N <: AbstractRefractiveIndex} <: AbstractRefractiveOptic{T, N}

Essentially represents the same functionality as [`Lens`](@ref).
Refer to its documentation.
"""
struct Prism{T, S <: AbstractShape{T}, N <: AbstractRefractiveIndex} <: AbstractRefractiveOptic{T, N}
    shape::S
    n::N
    function Prism(shape::S, n::N) where {T<:Real, S<:AbstractShape{T}, N<:AbstractRefractiveIndex}
        return new{T, S, N}(shape, n)
    end
end

thickness(p::Prism) = thickness(shape(p))

"""
    RightAnglePrism(leg_length, height, n)

Creates a right angle symmetric [`Prism`](@ref), see [`RightAnglePrismSDF`](@ref). The legs are the
faces with the normals -x and -y, the hypotenuse has the normal `(1, 1, 0) / √2`. The local origin
is the center of the hypotenuse face, the right angle edge lies at `(-leg_length / 2, -leg_length / 2)`.

As constructed, a beam along +y enters the leg with the normal -y, is totally internally reflected
at the hypotenuse (for `n > √2`) and leaves along -x; after `yrotate3d!(prism, π)` it leaves along +x.
`zrotate3d!(prism, 5π / 4)` turns the hypotenuse normal to -y, so that the hypotenuse faces a beam
along +y and returns it antiparallel.

# Arguments

- `leg_length`: dimension in x- and y-direction in [m]
- `height`: extent along the local z-axis in [m]
- `n`: [`AbstractRefractiveIndex`](@ref) of the prism
"""
function RightAnglePrism(leg_length::Real, height::Real, n::AbstractRefractiveIndex)
    shape = RightAnglePrismSDF(leg_length, height)
    return Prism(shape, n)
end

"""
    Prism(vertices, height, n)

Creates a [`Prism`](@ref) with a convex polygonal cross-section, extruded along the local z-axis.
The local origin is the origin of the given coordinates, i.e. the prism is placed where the
`vertices` are, and rotations act around it. Like the [`RightAnglePrism`](@ref) the prism is
*not aligned* with the y-axis, orient it with the kinematic functions. See [`PolygonPrismSDF`](@ref)
for the geometry. See also [`EquilateralPrism`](@ref) and [`DovePrism`](@ref).

# Arguments

- `vertices`: corners `(x, y)` of the cross-section in the local x-y-plane in [m], as a vector of
  2D points or tuples. Both vertex orders are accepted. The polygon must be strictly convex.
- `height`: extent along the local z-axis in [m], the prism is centered at `z = 0`
- `n`: [`AbstractRefractiveIndex`](@ref) of the prism
"""
function Prism(vertices::AbstractVector, height::Real, n::AbstractRefractiveIndex)
    return Prism(PolygonPrismSDF(vertices, height), n)
end

"""
    EquilateralPrism(side, height, n)

Creates an equilateral triangular [`Prism`](@ref), e.g. for dispersion. The local origin is the
centroid of the triangle. The apex points along the local +y-axis and the base face has the
normal -y. A beam travelling along +x or -x through the faces next to the apex at (approximately)
minimum deviation is deflected towards the base. For a beam along +y, `zrotate3d!(prism, π / 2 - θ)`
with `θ = asin(n sin 30°) - 30°` sets minimum deviation, the beam is deflected towards +x.

# Arguments

- `side`: edge length of the triangular cross-section in [m], must be positive
- `height`: extent along the local z-axis in [m]
- `n`: [`AbstractRefractiveIndex`](@ref) of the prism
"""
function EquilateralPrism(side::Real, height::Real, n::AbstractRefractiveIndex)
    side > 0 || throw(ArgumentError("the side length must be positive"))
    r = side * sqrt(oftype(float(side), 3)) / 6 # inradius
    return Prism([(-side / 2, -r), (side / 2, -r), (zero(r), 2r)], height, n)
end

"""
    DovePrism(length, aperture, height, n)

Creates a [`Prism`](@ref) with an isosceles trapezoid cross-section and 45° end faces, e.g. a
Dove prism that inverts an image without deviating a beam along its axis. The long axis is the
local y-axis, so that a beam travelling along +y enters the 45° end face. The base face with the
normal -x is the one that totally internally reflects. The local origin is the center of the
cross-section's bounding box. A beam along +y leaves the prism parallel to the y-axis, but in
general laterally shifted.

# Arguments

- `length`: length of the long base face in [m], must be larger than `2 * aperture`
- `aperture`: distance between the base face and the parallel short face in [m], must be positive
- `height`: extent along the local z-axis in [m]
- `n`: [`AbstractRefractiveIndex`](@ref) of the prism
"""
function DovePrism(length::Real, aperture::Real, height::Real, n::AbstractRefractiveIndex)
    aperture > 0 || throw(ArgumentError("the aperture must be positive"))
    length > 2aperture || throw(ArgumentError("the length must be larger than twice the aperture"))
    a, l = aperture / 2, length / 2
    return Prism([(-a, -l), (a, -(l - aperture)), (a, l - aperture), (-a, l)], height, n)
end
