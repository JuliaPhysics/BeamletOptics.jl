function _check_hole_diameter(hole_diameter, diameter)
    if !(0 < hole_diameter < diameter)
        throw(ArgumentError("hole_diameter must satisfy 0 < hole_diameter < diameter (got $hole_diameter, diameter $diameter)"))
    end
    return nothing
end

function _pierce_substrate(substrate, diameter, y_lo::T, y_hi::T, hole_diameter) where {T}
    hd = T(hole_diameter)
    _check_hole_diameter(hd, T(diameter))
    margin = T(10e-3)
    half_height = (y_hi - y_lo) / 2 + margin
    y_center = (y_hi + y_lo) / 2
    bore = CylinderSDF(hd / 2, half_height)
    translate3d!(bore, [zero(T), y_center, zero(T)])
    return substrate - bore
end
