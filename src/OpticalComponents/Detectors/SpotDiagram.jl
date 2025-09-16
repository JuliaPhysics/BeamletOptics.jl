spot_diagram(d::Detector) = spot_diagram(d, hits(d))

function spot_diagram(detector::Detector, hits::Vector{<:AbstractRayHit{T}}) where T
    res = Vector{Point2{T}}(undef, length(hits))
    # Transform global into local detector coordinates
    for (i, hit) in enumerate(hits)
        # Global hit pos
        hit_pos = position(hit)
        loc_pos = hit_pos - position(detector)
        x = dot(loc_pos, orientation(detector)[:,1])
        z = dot(loc_pos, orientation(detector)[:,3])
        res[i] = Point2{T}(x, z)
    end
    return res
end

function spot_diagram(::Detector, hits::Vector{B}) where B<:AbstractBeamletHit
    throw(ErrorException("Spot diagram not available for $B"))
end