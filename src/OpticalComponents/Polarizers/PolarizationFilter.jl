# Polarization Filter Component
# Reuses the generic Coating system under the hood with a JonesCoating model.

"""
    PolarizationFilter(shape; cutoff_strength)
    PolarizationFilter(edge_length; cutoff_strength)

Ideal linear polarizer. The device is represented by a zero-thickness
surface whose transmission axis is aligned with the local `x`-axis (transmitting local Ex)
and blocking axis with the local `z`-axis (blocking local Ez).

Rotate the underlying shape around the y-axis to align the axes with the desired global directions.
"""
function PolarizationFilter(shape::AbstractShape{T}; cutoff_strength=eps(T)) where {T}
    # Jones matrix transmitting x and blocking z
    model = JonesCoating(XZBasis(1.0, 0.0, 0.0, 0.0))
    return Coating(shape, model)
end

function PolarizationFilter(edge_length::Real; cutoff_strength=eps())
    T = float(typeof(edge_length))
    shape = QuadraticFlatMesh(T(edge_length))
    zrotate3d!(shape, T(π))
    set_new_origin3d!(shape)
    return PolarizationFilter(shape; cutoff_strength)
end
