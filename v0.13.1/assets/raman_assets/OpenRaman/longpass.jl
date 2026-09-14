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
    ::BMO.Beam{T, R},
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
        # above cut-on: behave like a refractive plate
        normal = BMO.normal3d(BMO.intersection(ray))
        lambda = BMO.wavelength(ray)
        if BMO.isentering(ray)
            # Entering ldm
            n1 = BMO.refractive_index(ray)
            n2 = BMO.refractive_index(ldm, lambda)
            # Hint to test ldm again
            hint = BMO.Hint(ldm)
        else
            # Exiting ldm
            n1 = BMO.refractive_index(ldm, lambda)
            n2 = BMO.refractive_index(system, lambda)
            hint = nothing
            # Flip normal for refraction3d
            normal = -normal
        end
        # Calculate new dir. and pos.
        ndir, TIR = BMO.refraction3d(BMO.direction(ray), normal, n1, n2)
        npos = BMO.position(ray) + BMO.length(ray) * BMO.direction(ray)
        # In case of TIR, update hint and n2
        if TIR
            hint = BMO.Hint(ldm)
            n2 = BMO.refractive_index(ldm, lambda)
        end
        return BMO.BeamInteraction{T, R}(
            hint,
            BMO.Ray{T}(npos, ndir, nothing, BMO.wavelength(ray), n2)
        )
    end
end
