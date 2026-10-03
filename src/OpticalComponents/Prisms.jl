"""
    Prism{T, S <: AbstractShape{T}, N <: RefractiveIndex} <: AbstractRefractiveOptic{T, N}

Essentially represents the same functionality as [`Lens`](@ref).
Refer to its documentation.
"""
struct Prism{T, S <: AbstractShape{T}, N <: RefractiveIndex} <: AbstractRefractiveOptic{T, N}
    shape::S
    n::N
    function Prism(shape::S, n::N) where {T<:Real, S<:AbstractShape{T}, N<:RefractiveIndex}
        test_refractive_index_function(n)
        return new{T, S, N}(shape, n)
    end
end

thickness(p::Prism) = thickness(shape(p))

"""
    RightAnglePrism(leg_length, height, n)

Creates a right angle symmetric [`Prism`](@ref). The prism is *not aligned* with the y-axis.

# Inputs

- `leg_length`: dimension in x- and y-direction in [m]
- `height`: in [m]
- `n`: [`RefractiveIndex`](@ref) of the prism, or a number for a constant refractive index
"""
function RightAnglePrism(leg_length::Real, height::Real, n::RefractiveIndex)
    shape = RightAnglePrismSDF(leg_length, height)
    return Prism(shape, n)
end
RightAnglePrism(leg_length::Real, height::Real, n::Real) = RightAnglePrism(leg_length, height, λ -> n)

"""
    Prism(vertices, height, n)

Creates a [`Prism`](@ref) with a convex polygonal cross-section, extruded along the local z-axis.

# Inputs

- `vertices`: corners `(x, y)` of the cross-section in the local x-y-plane in [m], as a vector of
  2D points or tuples. Both vertex orders are accepted. The polygon must be strictly convex.
- `height`: extent along the local z-axis in [m], the prism is centered at `z = 0`
- `n`: [`RefractiveIndex`](@ref) of the prism, or a number for a constant refractive index

The local origin is the origin of the given coordinates, i.e. the prism is placed where the
`vertices` are, and rotations act around it. Like the [`RightAnglePrism`](@ref) the prism is
*not aligned* with the y-axis, orient it with the kinematic functions. See [`PolygonPrismSDF`](@ref)
for the geometry. See also [`EquilateralPrism`](@ref) and [`DovePrism`](@ref).
"""
function Prism(vertices::AbstractVector, height::Real, n::RefractiveIndex)
    return Prism(PolygonPrismSDF(vertices, height), n)
end
Prism(vertices::AbstractVector, height::Real, n::Real) = Prism(vertices, height, λ -> n)

"""
    EquilateralPrism(side, height, n)

Creates an equilateral triangular [`Prism`](@ref), e.g. for dispersion.

# Inputs

- `side`: edge length of the triangular cross-section in [m]
- `height`: extent along the local z-axis in [m]
- `n`: [`RefractiveIndex`](@ref) of the prism, or a number for a constant refractive index

The local origin is the centroid of the triangle. The apex points along the local +y-axis and the
base face has the normal -y. A beam travelling along +x or -x through the faces next to the apex
at (approximately) minimum deviation is deflected towards the base.
"""
function EquilateralPrism(side::Real, height::Real, n::RefractiveIndex)
    r = side * sqrt(3) / 6 # inradius
    return Prism([(-side / 2, -r), (side / 2, -r), (zero(r), 2r)], height, n)
end
EquilateralPrism(side::Real, height::Real, n::Real) = EquilateralPrism(side, height, λ -> n)

"""
    DovePrism(length, aperture, height, n)

Creates a [`Prism`](@ref) with an isosceles trapezoid cross-section and 45° end faces, e.g. a
Dove prism that inverts an image without deviating a beam along its axis.

# Inputs

- `length`: length of the long base face in [m], must be larger than `2 * aperture`
- `aperture`: distance between the base face and the parallel short face in [m]
- `height`: extent along the local z-axis in [m]
- `n`: [`RefractiveIndex`](@ref) of the prism, or a number for a constant refractive index

The long axis is the local y-axis, so that a beam travelling along +y enters the 45° end face.
The base face (normal -x) is the one that totally internally reflects. The local origin is the
center of the cross-section's bounding box.
"""
function DovePrism(length::Real, aperture::Real, height::Real, n::RefractiveIndex)
    length > 2aperture || throw(ArgumentError("the length must be larger than twice the aperture"))
    a, l = aperture / 2, length / 2
    return Prism([(-a, -l), (a, -(l - aperture)), (a, l - aperture), (-a, l)], height, n)
end
DovePrism(length::Real, aperture::Real, height::Real, n::Real) = DovePrism(length, aperture, height, λ -> n)
