# Jones Matrix calculation helpers

function _calculate_global_E0(object::AbstractObject, ray::PolarizedRay,
        out_dir::AbstractArray, J::GlobalJonesBasis)
    in_dir = direction(ray)
    E0 = polarization(ray)
    # Transform Jones matrix according to global object orientation
    R = orientation(object)
    P = R * J * transpose(R)

    Q_in = I - in_dir * transpose(in_dir)
    Q_out = I - out_dir * transpose(out_dir)
    P = Q_out * P * Q_in + out_dir * transpose(in_dir)
    return P * E0
end

include("PolarizationFilter.jl")