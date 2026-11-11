"""
    LongpassDichroicMirror{T, N} <: BeamletOptics.AbstractObject{T}

A plane-parallel dichroic mirror that reflects everything below its `cuton` wavelength and
transmits everything above it. It is a body of the substrate `n` whose front face carries a
longpass coating; the transmitted light is refracted at both surfaces.

# Fields

- `shape`: a `PlanoSurfaceSDF` representing the substrate
- `n`: the substrate [`BeamletOptics.AbstractRefractiveIndex`](@ref)
- `coatings`: the coatings of its faces, see [`surface_model!`](@ref)
"""
struct LongpassDichroicMirror{T, N <: BMO.AbstractRefractiveIndex} <: BMO.AbstractObject{T}
    shape::BMO.PlanoSurfaceSDF{T}
    n::N
    coatings::BMO.FaceModels
end

BMO.refractive_index(lpm::LongpassDichroicMirror) = lpm.n
# the substrate is a body, whose medium the transitions at its surfaces take into account
BMO.medium(lpm::LongpassDichroicMirror) = BMO.IsotropicMedium(lpm.n)
# its faces carry coatings
BMO.face_models(lpm::LongpassDichroicMirror) = lpm.coatings

"""
    LongpassDichroicMirror(diameter, thickness, n, cuton)

Spawns a round [`LongpassDichroicMirror`](@ref) at the global origin. The reflective face
is normal to the local y-axis.
"""
function LongpassDichroicMirror(diameter::D, thickness::T, n::BMO.AbstractRefractiveIndex, cuton::Real) where {D, T}
    TT = promote_type(D, T)
    shape = BMO.PlanoSurfaceSDF(TT(thickness), TT(diameter))
    lpm = LongpassDichroicMirror(shape, n, BMO.FaceModels())
    longpass = SpectralCoating((λ, θ) -> λ < cuton ? 1.0 : 0.0; splitting = Dominant())
    surface_model!(lpm, FrontFace(), longpass)
    return lpm
end

BMO.render!(ax::LScene, lpm::LongpassDichroicMirror;
    color=:orange, transparency=true, alpha=0.5, kwargs...) = BMO.render!(ax, BMO.shape(lpm); color, transparency, alpha, kwargs...)
