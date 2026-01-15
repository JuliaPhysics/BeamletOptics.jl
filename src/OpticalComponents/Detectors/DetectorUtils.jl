abstract type AbstractCenterAlgorithm end

struct Centroid <: AbstractCenterAlgorithm end
struct MinMax <: AbstractCenterAlgorithm end

"""
    calc_center_point(::AbstractCenterAlgorithm, xs::Vector{T}, zs::Vector{T}, projection_factor::Vector{T}) where T

Calculates a central 2D point based on a distribution of points, specified by `xs` and `zs`.
The algorithm can be selected by specifying a concrete `AbstractCenterAlgorithm`.
In addition, projection factors for each hit point can be passed.

Returns `(x0, z0)`.
"""
calc_center_point(::ACA, ::Any, ::Any, ::Any) where ACA<:AbstractCenterAlgorithm = throw(ErrorException("calc_center_point not available for $ACA"))

function calc_center_point(::Centroid, xs::Vector{T}, zs::Vector{T}, projection_factor::Vector{T}) where T
    w_sum = zero(T)
    wx_sum = zero(T)
    wz_sum = zero(T)
    @inbounds for idx in eachindex(xs, zs, projection_factor)
        w = T(projection_factor[idx])
        w_sum += w
        wx_sum = muladd(w, xs[idx], wx_sum)
        wz_sum = muladd(w, zs[idx], wz_sum)
    end
    if iszero(w_sum)
        x0 = (minimum(xs) + maximum(xs)) / 2
        z0 = (minimum(zs) + maximum(zs)) / 2
    else
        x0 = wx_sum / w_sum
        z0 = wz_sum / w_sum
    end
    return x0, z0
end

function calc_center_point(::MinMax, xs::Vector{T}, zs::Vector{T}, ::Vector{T}) where T
    x0 = (minimum(xs) + maximum(xs)) / 2
    z0 = (minimum(zs) + maximum(zs)) / 2
    return x0, z0
end

"""
    calc_local_pos(pd::Detector; kwargs...)

Calculates the hit position of all registered hits in local detector coordinates.
"""
calc_local_pos(pd::Detector; kwargs...) = calc_local_pos(pd, hits(pd); kwargs...)

function calc_local_pos(::Detector, ::Vector{B}) where B<:AbstractDetectorHit
    throw(ErrorException("calc_local_pos not available for $B"))
end

function calc_local_pos(
        position::AbstractArray,
        local_x::AbstractArray,
        local_z::AbstractArray,
        hits_3D::Vector{Point3{T}}
    ) where T
    # Transform global into local detector coordinates
    hits_2D = Vector{Point2{T}}(undef, length(hits_3D))
    for (i, hit) in enumerate(hits_3D)
        loc_pos = hit - position
        @views x = dot(loc_pos, local_x)
        @views z = dot(loc_pos, local_z)
        hits_2D[i] = Point2{T}(T(x), T(z))
    end
    return hits_2D
end

function calc_local_pos(
        pd::Detector,
        hits::Vector{<:AbstractRayHit{T}}
    ) where T
    pd_pos = position(pd)
    # left handed coord. sys. for correct (x, z) orientation
    local_x = -orientation(pd)[:, 1]
    local_z = orientation(pd)[:, 3]
    hits_3D = hit_point.(hits)
    return calc_local_pos(pd_pos, local_x, local_z, hits_3D)
end

"""
    calc_local_pos(pd::Detector, hits::Vector{GaussianBeamletHit}; crop_factor=1, num_spots=50)

Calculates a projected circle or ellipse of 2D hit spots for each hit of a [`GaussianBeamlet`](@ref)
on the [`Detector`](@ref). The spot coordinates are returned in a left-handed (x, z) coordinate system
where the detector surface normal points towards the incoming beamlets.

# Keyword arguments

- `crop_factor=1`: scales the beam waist radius used to determine the bounding box
- `num_spots=50`: determines the number of 2D hits used to determine the bounding circle/ellipse
"""
function calc_local_pos(
        pd::Detector,
        hits::Vector{GaussianBeamletHit{T}};
        # kwargs
        crop_factor::Real=one(T),
        num_spots::Int=50
    ) where T
    # Calculate waist projections in local (x, z) coordinates
    ts = LinRange(0, 2pi, num_spots)
    p0 = position(pd)
    # Left-handed coord. sys. due to pd mesh rotation
    ex = -orientation(pd)[:,1]
    ey = orientation(pd)[:,2]
    ez = orientation(pd)[:,3]
    pts_2D = Vector{Point2{T}}()
    # Caclulate spots
    for hit in hits
        # Determine waist radius at intersection point
        waist_radius, ~, ~, ~ = gauss_parameters(hit, hit_point(hit))
        waist_radius *= crop_factor
        # build basis vectors
        dir = direction(hit)
        nd1 = normal3d(dir)
        nd2 = cross(dir, nd1)
        origin = hit_point(hit)
        # Calculate 3D waist circle of points
        pts_3D = ellipse(ts, origin, waist_radius*nd1, waist_radius*nd2)
        # Project 3D points onto plane, calculate local (x, z) coords
        for (i, pt) in enumerate(pts_3D)
            dist = line_plane_distance3d(p0, ey, pt, dir)
            pts_3D[i] = pt + dist * dir
        end
        new_pts = calc_local_pos(p0, ex, ez, pts_3D)
        append!(pts_2D, new_pts)
    end
    return pts_2D
end

"""
    calc_local_lims(pd::Detector; kwargs...)

Calculates the limiting values for the flat E-field evaluation grid based on the detector data.
"""
calc_local_lims(pd::Detector; kwargs...) = calc_local_lims(pd, hits(pd); kwargs...)

function calc_local_lims(::Detector, ::Vector{B}) where B<:AbstractDetectorHit
    throw(ErrorException("calc_local_lims not available for $B"))
end

"""
    calc_local_lims(pd::Detector, hits::Vector{<:AbstractRayHit}; crop_factor=1, center=Centroid())

Compute a symmetric [x_min,x_max]×[z_min,z_max] box around the hit positions weighted centroid for ray-based spot diagrams.

• If `center==Centroid()` (the default), uses
    x0 = ∑ wᵢ·xᵢ / ∑ wᵢ,  z0 = ∑ wᵢ·yᵢ / ∑ wᵢ
  with wᵢ = projection_factor.
• If `center==MinMax()`, falls back to the midpoint of [min,max].

Returns `(x_min, x_max, z_min, z_max)`.
"""
function calc_local_lims(
        pd::Detector,
        hits::Vector{<:AbstractRayHit{T}};
        # kwargs
        crop_factor::Real=one(T),
        center::AbstractCenterAlgorithm=Centroid()
    ) where T
    # get all hit‐points in local (x,y)
    local_hits = calc_local_pos(pd)
    xs = getindex.(local_hits, 1)
    zs = getindex.(local_hits, 2)
    # calculate center
    x0, z0 = calc_center_point(center, xs, zs, projection_factor.(hits))
    # half‐widths
    dx = maximum(x->abs(x - x0), xs)
    dy = maximum(y->abs(y - z0), zs)
    hwx, hwy = dx*crop_factor, dy*crop_factor
    # return rect. grid limits
    return x0 - hwx, x0 + hwx, z0 - hwy, z0 + hwy
end

"""
    calc_local_lims(pd::Detector, hits::Vector{GaussianBeamletHit}; crop_factor=1, num_spots=50, kwargs...)

Computes a 2D bounding box around the circular or elliptical waist of [`GaussianBeamlet`](@ref) hits.
Assumes that the beamlet is approximately cylindrical around its optical axis at the point of intersection.

# Keyword arguments

For available keyword args., refer to the corresponding [`calc_local_pos`](@ref) function.
"""
function calc_local_lims(
        pd::Detector,
        ::Vector{GaussianBeamletHit{G}};
        # kwargs
        crop_factor::Real=one(G),
        num_spots::Int=50
    ) where G
    local_hits = calc_local_pos(pd; crop_factor, num_spots)
    xs = getindex.(local_hits, 1)
    zs = getindex.(local_hits, 2)
    # min/max limits
    return minimum(xs), maximum(xs), minimum(zs), maximum(zs)
end
