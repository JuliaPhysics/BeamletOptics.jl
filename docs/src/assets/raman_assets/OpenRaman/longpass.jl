"""
    LongpassDichroicMirror{T, N} <: BeamletOptics.AbstractObject{T}

A plane-parallel dichroic mirror that reflects everything below its `cuton` wavelength and
transmits everything above it. Transmission is modelled as a regular refractive plate of
refractive index `n`, i.e. the beam is refracted at both surfaces.

# Fields

- `shape`: a `PlanoSurfaceSDF` representing the substrate
- `n`: the substrate [`BeamletOptics.RefractiveIndex`](@ref)
- `cuton`: the cut-on wavelength in [m]
"""
struct LongpassDichroicMirror{T, N <: BMO.RefractiveIndex} <: BMO.AbstractObject{T}
    shape::BMO.PlanoSurfaceSDF{T}
    n::N
    cuton::T
end

BMO.refractive_index(lpm::LongpassDichroicMirror) = lpm.n
BMO.refractive_index(lpm::LongpassDichroicMirror, λ::Real)::Float64 = lpm.n(λ)
# the substrate is a refracting body, which the transitions at its surfaces take into account
BMO.body_refractive_index(lpm::LongpassDichroicMirror, λ) = BMO.refractive_index(lpm, λ)

"""
    LongpassDichroicMirror(diameter, thickness, n, cuton)

Spawns a round [`LongpassDichroicMirror`](@ref) at the global origin. The reflective face
is normal to the local y-axis.
"""
function LongpassDichroicMirror(diameter::D, thickness::T, n::BMO.RefractiveIndex, cuton::Real) where {D, T}
    TT = promote_type(D, T)
    shape = BMO.PlanoSurfaceSDF(TT(thickness), TT(diameter))
    return LongpassDichroicMirror(shape, n, TT(cuton))
end

BMO.render!(ax::LScene, lpm::LongpassDichroicMirror;
    color=:orange, transparency=true, alpha=0.5, kwargs...) = BMO.render!(ax, BMO.shape(lpm); color, transparency, alpha, kwargs...)

function BMO.interact3d(
    system::BMO.AbstractSystem,
    ldm::LongpassDichroicMirror,
    beam::BMO.Beam{T, R},
    ray::R) where {T <: Real, R <: Ray{T}}
    # test cut-on wavelength
    if BMO.wavelength(ray) < ldm.cuton
        # below cut-on: behave like a plane mirror
        normal = BMO.normal3d(BMO.intersection(ray))
        npos = BMO.position(ray) + BMO.length(ray) * BMO.direction(ray)
        ndir = BMO.reflection3d(BMO.direction(ray), normal)
        return BMO.BeamInteraction{T, R}(
            nothing,
            Ray{T}(npos, ndir, nothing, BMO.wavelength(ray), BMO.refractive_index(ray))
        )
    else
        # above cut-on: behave like a refractive plate. The transition at the hit holds the media on
        # both sides, the substrate and the medium of the system or of a touching body.
        return BMO.cross_boundary(BMO.transition(system, ray), beam, ray)
    end
end
