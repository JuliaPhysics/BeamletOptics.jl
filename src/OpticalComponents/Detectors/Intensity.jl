intensity(d::Detector) = intensity(d, hits(d))

function intensity(d::Detector, hits::Vector{GaussianBeamletHit{T}}) where T
    
end