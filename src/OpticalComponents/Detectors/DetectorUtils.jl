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
    local_x, local_z = local_axes(pd)
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
    ex, ez = local_axes(pd)
    ey = direction(pd)
    pts_2D = Vector{Point2{T}}()
    sizehint!(pts_2D, length(hits) * num_spots)
    
    # Pre-calculate trig
    sin_t = sin.(ts)
    cos_t = cos.(ts)

    # Calculate spots
    for hit in hits
        # Use cached w_max (radius at detector)
        waist_radius = hit.w_max * crop_factor
        
        # build basis vectors
        dir = direction(hit)
        nd1 = normal3d(dir)
        nd2 = cross(dir, nd1)
        origin = hit_point(hit)
        
        # Calculate spots and project
        for k in 1:num_spots
            # 3D waist point
            pt = origin + waist_radius * (cos_t[k] * nd1 + sin_t[k] * nd2)
            # Project onto detector plane along ray direction
            dist = line_plane_distance3d(p0, ey, pt, dir)
            isnothing(dist) && continue
            p_proj = pt + dist * dir
            
            # Local coordinates
            v = p_proj - p0
            push!(pts_2D, Point2{T}(dot(v, ex), dot(v, ez)))
        end
    end
    return pts_2D
end

"""
    calc_local_pos(pd::Detector, hits::Vector{AstigmaticGaussianBeamletHit}; crop_factor=1, num_spots=50)

Calculates a projected ellipse of 2D hit spots for each hit of an [`AstigmaticGaussianBeamlet`](@ref)
on the [`Detector`](@ref).
"""
function calc_local_pos(
        pd::Detector,
        hits::Vector{AstigmaticGaussianBeamletHit{T}};
        # kwargs
        crop_factor::Real=one(T),
        num_spots::Int=50
    ) where T
    # Calculate waist projections in local (x, z) coordinates
    ts = LinRange(0, 2pi, num_spots)
    p0 = position(pd)
    # Left-handed coord. sys. due to pd mesh rotation
    ex, ez = local_axes(pd)
    ey = direction(pd)
    pts_2D = Vector{Point2{T}}()
    sizehint!(pts_2D, length(hits) * num_spots)

    sin_t = sin.(ts)
    cos_t = cos.(ts)

    for hit in hits
        # Use cached parabasal parameters at segment start
        # Propagate to detector (end of segment)
        # chief segment length
        l_seg = hit.l_seg
        h1 = (hit.h1 + l_seg * hit.u1) * crop_factor
        h2 = (hit.h2 + l_seg * hit.u2) * crop_factor
        
        dir = direction(hit)
        origin = hit_point(hit)
        
        for k in 1:num_spots
            pt = origin + cos_t[k] * h1 + sin_t[k] * h2
            # Project onto detector plane along ray direction
            dist = line_plane_distance3d(p0, ey, pt, dir)
            isnothing(dist) && continue
            p_proj = pt + dist * dir
            
            v = p_proj - p0
            push!(pts_2D, Point2{T}(real(dot(v, ex)), real(dot(v, ez))))
        end
    end
    return pts_2D
end

"""
    calc_local_lims(pd::Detector)

Calculates the limiting values for the flat E-field evaluation grid based on the detector data,
with the default `crop_factor` of the hit type.
"""
calc_local_lims(pd::Detector) = calc_local_lims(pd, hits(pd), nothing)

function calc_local_lims(
        ::Detector,
        ::Vector{B},
        crop_factor;
        # kwargs
        center::AbstractCenterAlgorithm=Centroid()
    ) where B
    error(lazy"calc_local_lims not defined for $B")
end

"""
    calc_local_lims(pd::Detector, ::Nothing, crop_factor; center=Centroid())

The window of a [`Detector`](@ref) without hits: the bounding box of its active surface in
local coordinates. `crop_factor` and `center` are ignored, since there are no hits to scale or
center on.
"""
function calc_local_lims(
        pd::Detector{T},
        ::Nothing,
        crop_factor;
        # kwargs
        center::AbstractCenterAlgorithm=Centroid()
    ) where T
    corners = Point3{T}.(eachrow(vertices(shape(pd))))
    local_corners = calc_local_pos(position(pd), local_axes(pd)..., corners)
    xs = getindex.(local_corners, 1)
    zs = getindex.(local_corners, 2)
    return minimum(xs), maximum(xs), minimum(zs), maximum(zs)
end

"""
    calc_local_lims(pd::Detector, hits::Vector{<:AbstractRayHit}, crop_factor; center=Centroid())

Compute a symmetric [x_min,x_max]×[z_min,z_max] box around the hit positions weighted centroid for ray-based spot diagrams.

• If `center==Centroid()` (the default), uses
    x0 = ∑ wᵢ·xᵢ / ∑ wᵢ,  z0 = ∑ wᵢ·yᵢ / ∑ wᵢ
  with wᵢ = projection_factor.
• If `center==MinMax()`, falls back to the midpoint of [min,max].

`crop_factor` scales the half-widths of the box, `nothing` selects the default of `1`.

Returns `(x_min, x_max, z_min, z_max)`.
"""
function calc_local_lims(
        pd::Detector,
        hits::Vector{<:AbstractRayHit{T}},
        crop_factor::Real;
        # kwargs
        center::AbstractCenterAlgorithm=Centroid()
    ) where T
    # get provided hit‐points in local (x,y)
    local_hits = calc_local_pos(pd, hits)
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

calc_local_lims(
    pd::Detector,
    hits::Vector{<:AbstractRayHit{T}},
    ::Nothing;
    # kwargs
    center::AbstractCenterAlgorithm=Centroid()
) where T = calc_local_lims(pd, hits, one(T); center)

"""
    calc_local_lims(pd::Detector, hits::Vector{<:AbstractBeamletHit}, crop_factor; center=Centroid())

Computes a 2D bounding box around the circular or elliptical waist of [`GaussianBeamlet`](@ref) hits.
Assumes that the beamlet is approximately cylindrical around its optical axis at the point of intersection.
The box reaches `crop_factor` beam radii beyond each hit point, divided by the `cos θ` of the
angle of incidence, which stretches the footprint of the beam on the detector. `nothing` selects
the default `crop_factor` of `3`.

The `center` keyword is ignored and only accepted so that all hit types share one signature.
"""
function calc_local_lims(
        pd::Detector,
        hits::Vector{H},
        crop_factor::Real;
        # kwargs
        center::AbstractCenterAlgorithm=Centroid()
    ) where {T, H <: AbstractBeamletHit{T}}

    pos_pd = position(pd)
    ex, ez = local_axes(pd)
    ey = direction(pd)

    xmin = T(Inf)
    xmax = T(-Inf)
    zmin = T(Inf)
    zmax = T(-Inf)
    
    for hit in hits
        # Center in local coords
        p_hit = hit_point(hit) - pos_pd
        x = dot(p_hit, ex)
        z = dot(p_hit, ez)
        
        # Footprint stretched by 1/cos θ
        w = hit.w_max * crop_factor / abs(dot(direction(hit), ey))
        
        xmin = min(xmin, x - w)
        xmax = max(xmax, x + w)
        zmin = min(zmin, z - w)
        zmax = max(zmax, z + w)
    end
    
    # Handle empty or invalid hits
    if xmin > xmax
        return (T(0), T(0), T(0), T(0))
    end
    
    return (xmin, xmax, zmin, zmax)
end

calc_local_lims(
    pd::Detector,
    hits::Vector{<:AbstractBeamletHit},
    ::Nothing;
    # kwargs
    center::AbstractCenterAlgorithm=Centroid()
) = calc_local_lims(pd, hits, 3; center)
