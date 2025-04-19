"""
    Spotdetector <: AbstractDetector

Simple 2D screen that stores the intersection point of incoming [`Beam`](@ref)s.
The intersection points are stored in local coordinates of the detector with respect to the screen origin.

# Fields

- `shape`: a 2D [`QuadraticFlatMesh`](@ref) that is **aligned with the negative y-axis**
- `data`: stores the intersection points as `Point2`
- `hw`: half-width of the detector plane

# Additional information

!!! info "Normal vector"
    Check the normal vector orientation of the detector plane if the spot diagram looks mirrored.

!!! warning "Reset behavior"
    Spot diagram data must be manually reset between traces via [`empty!`](@ref)
"""
mutable struct Spotdetector{T} <: AbstractDetector{T, Mesh{T}}
    const shape::Mesh{T}
    data::Vector{Point2{T}}
    hw::T
end

push!(sd::Spotdetector{T}, new::Point2{T}) where T = push!(sd.data, new)

"""
    Spotdetector(width)

Generates a quadratic rectangular 2D [`Spotdetector`](@ref) that is aligned with the **negative y-axis**.
Refer to the type docs for more information.

# Inputs:

- `width`: edge length in [m]
"""
function Spotdetector(width::W) where W<:AbstractFloat
    # Spawn mesh, align with neg. y-axis, empty data field
    shape = QuadraticFlatMesh(width)
    zrotate3d!(shape, π)
    data = Vector{Point2{W}}()
    return Spotdetector(shape, data, width/2)
end

"""Resets the stored spot diagram data"""
empty!(sd::Spotdetector) = empty!(sd.data)

function interact3d(::AbstractSystem, sd::Spotdetector, beam::Beam{T, R}, ray::R) where {T <: Real, R <: AbstractRay{T}}
    # Calculate intersection in global coordinates
    hit_pos = position(ray) + length(ray) * direction(ray)
    # Transform into local detector coordinates
    loc_pos = hit_pos - position(sd)
    x = dot(loc_pos, orientation(sd)[:,1])
    z = dot(loc_pos, orientation(sd)[:,3])
    # Push point into data field
    d = Point2{T}(x, z)
    push!(sd, d)
    return nothing
end