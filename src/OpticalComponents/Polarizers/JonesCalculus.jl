abstract type AbstractJonesPolarizer{T, S} <: AbstractObject{T, S} end

function _calculate_global_E0(object::AbstractJonesPolarizer, ray::PolarizedRay, out_dir::AbstractArray, J::GlobalJonesBasis)
    in_dir = direction(ray)
    E0 = polarization(ray)
    if !isparallel3d(in_dir, out_dir)
        # Not sure how good this implementation works
        throw(ErrorException("Thin film polarizer _calculate_global_E0 only works if no direction change occurs."))
    end
    # Transform Jones matrix according to global object orientation
    R = orientation(object)
    P = R * J * transpose(R)
    # Calculate projection of P into ray-E0-transverse plane
    # Note: this in general violates P*in_dir = out_dir
    Q = I - in_dir * transpose(in_dir)
    # technically Q == transpose(Q)
    P = Q * P * transpose(Q)
    return P*E0
end

include("PolarizationFilter.jl")