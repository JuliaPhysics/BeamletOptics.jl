@Geometry struct Mirror end

@Geometry struct Lens
    ref_index::Float64
    function Lens(ref_index)
        @assert ref_index >= 1 "It is assumed that n>=1!"
        new(Float64(ref_index))
    end
end