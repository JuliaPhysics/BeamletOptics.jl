"""
    BeamsplitterCoating{T} <: AbstractSurfaceModel

A lossless, non-polarizing coating that splits a ray into a transmitted and a reflected beam with the
field factors `reflectance` and `transmittance` (R² + T² = 1). It acts between the media of the
boundary, see [`Transition`](@ref): the transmitted beam is refracted into the medium behind it, e.g.
the substrate that it is coated on, the reflected beam stays in the medium before it. Beyond the
critical angle of these media, the ray is totally reflected and not split.

A [`ThinBeamsplitter`](@ref) carries this coating, see [`surface_model`](@ref).
"""
struct BeamsplitterCoating{T <: Real} <: AbstractSurfaceModel
    reflectance::T
    transmittance::T
end

reflectance(coating::BeamsplitterCoating) = coating.reflectance
transmittance(coating::BeamsplitterCoating) = coating.transmittance

"""
    ThinBeamsplitter <: AbstractBeamsplitter

Represents a 2D beam-splitting device: a thin layer whose [`BeamsplitterCoating`](@ref) acts at the
boundary on which it lies, see [`surface_model`](@ref).

# Fields

- `shape`: 2D [`AbstractShape`](@ref) at which the splitting process occurs (e.g. a 2D-[`Mesh`](@ref))
- `coating`: the [`BeamsplitterCoating`](@ref) with the field factors `reflectance` and `transmittance`

!!! warning
    The stored `reflectance` and `transmittance` are amplitudes. They should be calculated from the
    power reflectance ``R`` such that `reflectance^2 + transmittance^2 = 1`, as checked by `isvalid`.
"""
struct ThinBeamsplitter{T <: Real, S <: AbstractShape{T}} <: AbstractBeamsplitter{T}
    shape::S
    coating::BeamsplitterCoating{T}
end

# The factors in the number type of the `shape`, e.g. of a shape built from dual numbers for
# automatic differentiation w.r.t. its pose, while the reflectance is a `Float64`
function ThinBeamsplitter(shape::AbstractShape{T}, reflectance::Real, transmittance::Real) where {T}
    return ThinBeamsplitter{T, typeof(shape)}(shape, BeamsplitterCoating(T(reflectance), T(transmittance)))
end

surface_model(bs::ThinBeamsplitter) = bs.coating
reflectance(bs::ThinBeamsplitter) = reflectance(surface_model(bs))
transmittance(bs::ThinBeamsplitter) = transmittance(surface_model(bs))

"""
    ThinBeamsplitter(width, height; reflectance=0.5)

Creates a zero-thickness, lossless, non-polarizing 2D rectangular [`ThinBeamsplitter`](@ref) where

- `width`: is the x-dir. edge length in [m]
- `height`: is the z-dir. edge length in [m]
- `reflectance`: kw-arg in the open interval (0, 1) that determines how much light is **reflected**, i.e. 0.7 for a 70:30 splitter

# Additional information

!!! info "Reflectance"
    The input `reflectance` is the power reflectance ``R``. The amplitudes ``\\sqrt{R}`` and ``\\sqrt{1 - R}``
    are stored in the fields `reflectance` and `transmittance`, such that their squares sum to one.

!!! warning "Reflection phase jump"
    The reflection phase jump θᵣ = π applies to a beam that faces the normal of the splitter, see
    [`channels`](@ref) of its [`BeamsplitterCoating`](@ref).
"""
function ThinBeamsplitter(width::Real, height::Real; reflectance::Real = 0.5)
    if reflectance ≥ 1 || reflectance ≈ 0
        error("Splitting ratio ∈ (0, 1)!")
    end
    shape = RectangularFlatMesh(width, height)
    Reflected = sqrt(reflectance)
    Transmitted = sqrt(1 - Reflected^2)
    return ThinBeamsplitter(shape, Reflected, Transmitted)
end

"""
    RoundThinBeamsplitter(diameter; reflectance=0.5)

Creates a zero-thickness, 2D round [`ThinBeamsplitter`](@ref) with the specified `diameter` in [m].
For more information, refer to the [`ThinBeamsplitter`](@ref) constructor.
"""
function RoundThinBeamsplitter(diameter::Real; reflectance::Real = 0.5)
    if reflectance ≥ 1 || reflectance ≈ 0
        error("Splitting ratio ∈ (0, 1)!")
    end
    shape = CircularFlatMesh(diameter / 2)
    R = sqrt(reflectance)
    T = sqrt(1 - R^2)
    return ThinBeamsplitter(shape, R, T)
end

function ThinBeamsplitter(width::Real; reflectance::Real = 0.5)
    ThinBeamsplitter(width, width; reflectance)
end

Base.isvalid(bs::AbstractBeamsplitter) = reflectance(bs)^2 + transmittance(bs)^2 ≈ 1


splitting_of(::BeamsplitterCoating) = AllChannels()

"""
    channels(coating::BeamsplitterCoating, t::Transition, ray, cs)

The channels of the `coating`: the transmitted channel is refracted into the medium behind the
boundary with the field factor `transmittance(coating)`, the reflected one stays in the medium before
it with `reflectance(coating)`. Beyond the critical angle of the transition `t`, the ray is totally
reflected and not split, the channels `cs` of the bare boundary are returned.

The field of a transmitted [`PolarizedRay`](@ref) is scaled for the medium behind the boundary: the
coating transmits a fraction of the power, and the power of a beam is proportional to
`n cos(θ) |E|²`.

# Reflection phase jump

The reflection phase jump is θᵣ = π for a ray that faces the normal of the coating, i.e.
`isentering`, and 0 from the other side, as required by the Stokes relations (r' = -r) for a lossless
coating: only then do both outputs of an interferometer sum to its input. It applies to the field factor of beams
without polarization, e.g. a [`GaussianBeamlet`](@ref), and to the field of a [`PolarizedRay`](@ref),
whose s- and p-components are reflected with `-r` and `r` from the other side, like an internal Fresnel
reflection. The coating is assumed to have an unambiguous normal, i.e. a 2D mesh.
"""
function channels(coating::BeamsplitterCoating, t::Transition, ray::AbstractRay, cs)
    transmitted, reflected = cs
    # Total internal reflection: the ray does not reach the coating
    islit(transmitted) || return cs
    r, τ = reflectance(coating), transmittance(coating)
    s = isentering(ray) ? -1 : 1
    # The transmitted field scaled for the medium behind the boundary
    normal = normal3d(t)
    cos_in = abs(cos_angle3d(direction(ray), normal))
    cos_out = abs(cos_angle3d(direction(transmitted), normal))
    f = τ * sqrt(refractive_index_in(t) * cos_in / (refractive_index_out(t) * cos_out))
    return (transfer(transmitted, ray, SPBasis(f, 0, 0, f), τ),
        transfer(reflected, ray, SPBasis(-s * r, 0, 0, s * r), s * r))
end
