struct PolarizationFilter{T, S <: AbstractShape{T}} <: AbstractObject{T, S}
    shape::S
end

function PolarizationFilter(size::T) where {T <: Real}
    shape = QuadraticFlatMesh(size)
    return PolarizationFilter(shape)
end

function interact3d(::AbstractSystem,
        polfilter::PolarizationFilter,
        ::Beam{T, R},
        ray::R) where {T <: Real, R <: PolarizedRay{T}}
    npos = position(ray) + length(ray) * direction(ray)
    ndir = direction(ray)
    o = orientation(polfilter)

    # Jones polarization filter matrix
    J_polfilter = @SArray [1 0 0; 0 0 0; 0 0 1]
    J = inv(o) * J_polfilter * o # Update jones matrice according to filter orientation in space

    E0 = _calculate_global_E0(direction(ray), ndir, J, polarization(ray))

    return BeamInteraction{T, R}(nothing,
        PolarizedRay{T}(npos, ndir, nothing, wavelength(ray), refractive_index(ray), E0))
end