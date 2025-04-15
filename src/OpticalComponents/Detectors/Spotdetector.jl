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
    Spot diagram data must be manually reset between traces via [`reset_detector!`](@ref)
"""
mutable struct Spotdetector{T} <: AbstractDetector{T, Mesh{T}}
    const shape::Mesh{T}
    data::Vector{Point2{T}}
    hw::T
end

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
reset_detector!(sd::Spotdetector{T}) where T = (sd.data = Vector{Point2{T}}())

function interact3d(::AbstractSystem, sd::Spotdetector, beam::Beam{T, R}, ray::R) where {T <: Real, R <: AbstractRay{T}}
    # Calculate intersection in global coordinates
    hit_pos = position(ray) + length(ray) * direction(ray)
    # Transform into local detector coordinates
    loc_pos = hit_pos - position(sd)
    x = dot(loc_pos, orientation(sd)[:,1])
    z = dot(loc_pos, orientation(sd)[:,3])
    # Push point into data field
    d = Point2{T}(x, z)
    push!(sd.data, d)
    return nothing
end

"""
    create_spot_diagram(system, beam, aperture; n_rings, n_rays)

Calculates the spot diagram for concentric rings of **collimated beams** within the specified `aperture` diameter.
Uses retracing of an input `beam`. **System must be aligned onto the global y-axis.**

!!! info
    The input `system` must feature **one** [`Spotdetector`](@ref) at the approx. focal plane. 
    The detector is automatically reset when calling this function.

# Inputs

- `aperture`: maximum aperture diameter in [m]
- `n_rings`: maximum number of concentric rings
- `n_rays`: total number of retracing runs, i.e. "spawned beams"
"""
function create_spot_diagram(system::AbstractSystem, beam::Beam{T}, aperture::Real; n_rings::Int=20, n_rays::Int=1000) where T
    # test if detector present
    _objects = collect(objects(system))
    has_spotdetector = typeof.(_objects) .<: Spotdetector
    if !any(has_spotdetector)
        error("No Spotdetector in system")
    end
    obj_index = findall(has_spotdetector)
    if length(obj_index) > 1
        error("Only one Spotdetector allowed")
    end
    sd = _objects[obj_index[1]]
    # reset detector
    reset_detector!(sd)
    # setup concentric beam rings
    r_max = aperture/2
    rs = LinRange(0, r_max, n_rings)
    # determine length step ds via total circumference
    u_total = 0
    for (i, r) in enumerate(rs)
        i == 0 && continue
        u_total += 2π*r
    end
    ds = u_total/n_rays
    # Retrace input beam
    for r in rs
        n_rays = round(Int, 2π*r / ds)
        xs = [cos(x)*r for x in LinRange(0, 2pi * (n_rays - 1) / n_rays, n_rays)]
        zs = [sin(x)*r for x in LinRange(0, 2pi * (n_rays - 1) / n_rays, n_rays)]
        for i in eachindex(xs)
            y = position(first(beam.rays))[2]
            position!(first(beam.rays), Point3{Float64}(xs[i], y, zs[i]))
            solve_system!(system, beam)
        end 
    end
    return nothing
end