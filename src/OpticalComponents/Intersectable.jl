"""
    IntersectableObject

A passive [`AbstractObject`](@ref) which can be hit by a beam. In this case, the object acts like a hard target and blocks the beam path.

# Fields

- `shape`: an [`AbstractShape`](@ref)
- `stop`: Boolean flag indicating if rays can pass through the object, default is `true`
"""
mutable struct IntersectableObject{T, S <: AbstractShape{T}} <: AbstractObject{T, S}
    const shape::S
    stop::Bool
end

stop(io::IntersectableObject) = io.stop

set_new_origin3d!(d::IntersectableObject) = set_new_origin3d!(d.shape)

function interact3d(::AbstractSystem, ::IntersectableObject, ::AbstractBeam, ::AbstractRay)
    @warn "Continue tracing not implemented for this beam type..."
    return nothing
end

function interact3d(::AbstractSystem, io::IntersectableObject, beam::Beam{T, R}, ray::R) where {T <: Real, R <: Ray{T}}
    if stop(io)
        return nothing
    end

    npos = position(ray) + length(ray) * direction(ray)
    ndir = direction(ray)

    return BeamInteraction{T, R}(
        nothing,
        Ray{T}(
            npos,
            ndir, 
            nothing,
            wavelength(ray),
            refractive_index(ray)
        )
    )
end

IntersectableObject(loadpath::String, stop=true) = IntersectableObject(Mesh(load(loadpath)), stop)