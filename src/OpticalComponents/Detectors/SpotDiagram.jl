spot_diagram(d::Detector) = spot_diagram(d, hits(d))

function spot_diagram(d::Detector, hits::Vector{<:AbstractRayHit{T}}) where T
    res = Vector{Point2{T}}(undef, length(hits))
    # Transform global into local detector coordinates
    for (i, h) in enumerate(hits)
        ray = h.ray
        # Global hit pos
        hit_pos = position(ray) + length(ray) * direction(ray)
        loc_pos = hit_pos - position(d)
        x = dot(loc_pos, orientation(d)[:,1])
        z = dot(loc_pos, orientation(d)[:,3])
        res[i] = Point2{T}(x, z)
    end
    return res
end

function spot_diagram(::Detector, hits::Vector{B}) where B<:AbstractBeamletHit
    throw(ErrorException("Spot diagram not available for $B"))
end