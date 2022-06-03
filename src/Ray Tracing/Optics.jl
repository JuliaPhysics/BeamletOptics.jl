@Geometry struct Mirror end

function intersect3d(object::Mirror, ray::Ray)
    return intersect3d(object.geometry, ray)
end

@Geometry struct Lens
    ref_index::Float64
    function Lens(ref_index)
        @assert ref_index >= 1 "It is assumed that n>=1!"
        new(Float64(ref_index))
    end
end

function intersect3d(object::Lens, ray::Ray)
    return intersect3d(object.geometry, ray)
end