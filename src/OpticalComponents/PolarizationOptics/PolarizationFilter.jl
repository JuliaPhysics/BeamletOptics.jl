struct PolarizationFilter{T, S <: AbstractShape{T}} <: AbstractObject{T,S}
    shape::S
    JMat::GlobalJonesBasis{T}
end

function PolarizationFilter(shape::S, J::GlobalJonesBasis{TJ}) where {TS, S<:AbstractShape{TS}, TJ}
    return PolarizationFilter{TS,S}(shape, GlobalJonesBasis{TS}(J))
end

function PolarizationFilter(size::T, JMat=XZBasis(1, 0, 0, 0)) where {T <: Real}
    shape = QuadraticFlatMesh(size)
    return PolarizationFilter(shape, JMat)
end

function interact3d(::AbstractSystem,
        polfilter::PolarizationFilter,
        ::Beam{T, R},
        ray::R) where {T <: Real, R <: PolarizedRay{T}}
    npos = position(ray) + length(ray) * direction(ray)
    ndir = direction(ray)

    E0 = _calculate_global_E0(polfilter, ray, ndir, polfilter.JMat)

    return BeamInteraction{T, R}(nothing,
        PolarizedRay{T}(npos, ndir, nothing, wavelength(ray), refractive_index(ray), E0))
end
