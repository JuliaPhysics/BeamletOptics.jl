spot_diagram(d::Detector) = spot_diagram(d, hits(d))

function spot_diagram(pd::Detector, ::Vector{<:AbstractRayHit{T}}) where T
    return calc_local_pos(pd)
end

function spot_diagram(::Detector, hits::Vector{B}) where B<:AbstractBeamletHit
    throw(ErrorException("Spot diagram not available for $B"))
end