"""
    calc_local_pos(pd::Detector)

Calculates the hit position of all registered hits in local detector coordinates.
"""
calc_local_pos(pd::Detector) = calc_local_pos(pd, hits(pd))

function calc_local_pos(
        pd::Detector,
        hits::Vector{<:AbstractRayHit{T}}
    ) where T
    pd_pos = position(pd)
    local_x = orientation(pd)[:, 1]
    local_z = orientation(pd)[:, 3]
    # Transform global into local detector coordinates
    hits_2D = Vector{Point2{T}}(undef, length(hits))
    for (i, hit) in enumerate(hits)
        _hit = position(hit)
        loc_pos = _hit - pd_pos
        @views x = dot(loc_pos, local_x)
        @views z = dot(loc_pos, local_z)
        hits_2D[i] = Point2{T}(T(x), T(z))
    end
    return hits_2D
end

function calc_local_pos(::Detector, ::Vector{B}) where B<:AbstractBeamletHit
    throw(ErrorException("calc_local_pos not available for $B"))
end

"""
    calc_local_lims(pd::Detector; kwargs...)

Calculates the limiting values for the flat E-field evaluation grid based on the detector data.
"""
calc_local_lims(pd::Detector; kwargs...) = calc_local_lims(pd, hits(pd); kwargs...)

"""
    calc_local_lims(pd::Detector, hits::Vector{<:AbstractRayHit}; crop_factor=1, center=:centroid)

Compute a symmetric [x_min,x_max]×[z_min,z_max] box around the hit positions weighted centroid.

• If `center==:centroid` (the default), uses
    x0 = ∑ wᵢ·xᵢ / ∑ wᵢ,  z0 = ∑ wᵢ·yᵢ / ∑ wᵢ
  with wᵢ = projection_factor.
• If `center==:bbox`, falls back to the midpoint of [min,max].

Returns `(x_min, x_max, z_min, z_max)`.
"""
function calc_local_lims(
        pd::Detector,
        hits::Vector{<:AbstractRayHit{T}};
        # kwargs
        crop_factor::Real=one(T),
        center::Symbol=:centroid
    ) where T
    # get all hit‐points in local (x,y)
    local_hits = calc_local_pos(pd)
    xs = getindex.(local_hits, 1)
    zs = getindex.(local_hits, 2)

    # choose center
    if center == :centroid
        w_sum = sum(projection_factor, hits)
        x0 = sum(x -> (projection_factor(x[1]) * x[2]), zip(hits, xs)) / w_sum
        z0 = sum(x -> (projection_factor(x[1]) * x[2]), zip(hits, zs)) / w_sum
    else
        x0 = (minimum(xs) + maximum(xs)) / 2
        z0 = (minimum(zs) + maximum(zs)) / 2
    end

    # half‐widths
    dx = maximum(x->abs(x - x0), xs)
    dy = maximum(y->abs(y - z0), zs)
    hwx, hwy = dx*crop_factor, dy*crop_factor

    return x0 - hwx, x0 + hwx, z0 - hwy, z0 + hwy
end

function calc_local_lims(::Detector, ::Vector{B}) where B<:AbstractBeamletHit
    throw(ErrorException("calc_local_lims not available for $B"))
end

function calc_local_lims(
        pd::Detector,
        hits::Vector{GaussianBeamletHit{G}};
        # kwargs
        kwargs...
    ) where G
    
end
