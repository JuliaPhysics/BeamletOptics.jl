"""
    AbstractSide

The side of a boundary into which a [`BoundaryChannel`](@ref) leads: [`Reflected`](@ref) back into the
medium before the boundary, or [`Transmitted`](@ref) into the medium behind it. The refractive index,
the hint and the enclosure of the new ray follow from the side.
"""
abstract type AbstractSide end

"""The side of a [`BoundaryChannel`](@ref) back into the medium before the boundary."""
struct Reflected <: AbstractSide end

"""The side of a [`BoundaryChannel`](@ref) into the medium behind the boundary."""
struct Transmitted <: AbstractSide end

"""
    BoundaryChannel{T, S, M}

A channel of a boundary: a direction into which the ray leaves the boundary on the `side`
([`Reflected`](@ref) or [`Transmitted`](@ref)), and the transfer of its field. A surface model
transforms the channels of the bare boundary, see [`channels`](@ref); the solver builds the new rays
from them.

# Fields

- `side`: the [`AbstractSide`](@ref) of the boundary
- `dir`: the direction of the new ray
- `P`: the 3×3 transfer of the field of a [`PolarizedRay`](@ref) in global coordinates,
  `E_out = P * E_in`, or `nothing` for a ray without a field
- `a`: the complex amplitude factor for beams without polarization, e.g. a
  [`GaussianBeamlet`](@ref); a channel with `a = 0` carries no light

A channel of a bare boundary has `a = 1` on the side into which the ray continues and `a = 0` on the
other side: a beam without polarization is refracted without Fresnel losses.
"""
struct BoundaryChannel{T, S <: AbstractSide, M}
    side::S
    dir::Point3{T}
    P::M
    a::Complex{T}
end

BoundaryChannel(side::AbstractSide, dir::AbstractVector{T}, P, a::Number) where {T} =
    BoundaryChannel{T, typeof(side), typeof(P)}(side, Point3{T}(dir), P, Complex{T}(a))

side(c::BoundaryChannel) = c.side
direction(c::BoundaryChannel) = c.dir
amplitude(c::BoundaryChannel) = c.a

"""Returns a copy of the channel `c` with the field transfer `P` and the amplitude factor `a`."""
transfer(c::BoundaryChannel, P, a::Number) = BoundaryChannel(side(c), direction(c), P, a)

"""The refractive index on the `side` of the boundary of the transition `t`."""
refractive_index(t::Transition, ::Reflected) = refractive_index_in(t)
refractive_index(t::Transition, ::Transmitted) = refractive_index_out(t)

"""The hint into the medium on the `side` of the boundary of the transition `t`."""
side_hint(t::Transition, ::Reflected) = hint_in(t)
side_hint(t::Transition, ::Transmitted) = hint_out(t)

"""
    channels(t::Transition, ray)

Returns the channels of the bare boundary of the transition `t` for the `ray`, as the tuple
`(reflected, transmitted)` of [`BoundaryChannel`](@ref)s. The `ray` is refracted from the medium before into
the medium behind the boundary via [`refraction3d`](@ref); beyond the critical angle it is totally
reflected, and the transmitted channel carries no light. The field of a [`PolarizedRay`](@ref) is
transferred with the [`fresnel_coefficients`](@ref). Between equal indices the ray passes unchanged.
"""
function channels(t::Transition, ray::Ray)
    d = direction(ray)
    n1, n2 = refractive_index_in(t), refractive_index_out(t)
    reflected = reflection3d(d, normal3d(t))
    n1 == n2 && return (BoundaryChannel(Reflected(), reflected, nothing, 0), BoundaryChannel(Transmitted(), d, nothing, 1))
    refracted, TIR = refraction3d(d, normal3d(t), n1, n2)
    TIR && return (BoundaryChannel(Reflected(), reflected, nothing, 1), BoundaryChannel(Transmitted(), reflected, nothing, 0))
    return (BoundaryChannel(Reflected(), reflected, nothing, 0), BoundaryChannel(Transmitted(), refracted, nothing, 1))
end

function channels(t::Transition, ray::PolarizedRay{T}) where {T}
    d = direction(ray)
    normal = normal3d(t)
    n1, n2 = refractive_index_in(t), refractive_index_out(t)
    reflected = reflection3d(d, normal)
    O = zero(SMatrix{3, 3, Complex{T}, 9})
    if n1 == n2
        P = SMatrix{3, 3, Complex{T}, 9}(I)
        return (BoundaryChannel(Reflected(), reflected, O, 0), BoundaryChannel(Transmitted(), d, P, 1))
    end
    # Fresnel coefficients from the cosine, smooth for AD at normal incidence
    rs, rp, ts, tp = _fresnel_coefficients(cos_angle3d(d, -normal), n2 / n1)
    Pr = _calculate_global_E0(d, reflected, normal3d(intersection(ray)), SPBasis(-rs, 0, 0, rp))
    if is_internally_reflected(rp, rs)
        return (BoundaryChannel(Reflected(), reflected, Pr, 1), BoundaryChannel(Transmitted(), reflected, O, 0))
    end
    refracted, _ = refraction3d(d, normal, n1, n2)
    Pt = _calculate_global_E0(d, refracted, normal3d(intersection(ray)), SPBasis(ts, 0, 0, tp))
    return (BoundaryChannel(Reflected(), reflected, Pr, 0), BoundaryChannel(Transmitted(), refracted, Pt, 1))
end

"""
    channels(model::AbstractSurfaceModel, t::Transition, ray, cs)

Returns the channels of the boundary of the transition `t` for the `ray` with the surface `model`,
from the channels `cs` of the boundary without it, see [`BoundaryChannel`](@ref). This is the hook of a
surface model; [`Uncoated`](@ref) returns `cs`.
"""
channels(model::AbstractSurfaceModel, ::Transition, ::AbstractRay, cs) =
    error(lazy"channels not defined for $(typeof(model))")

channels(::Uncoated, ::Transition, ::AbstractRay, cs) = cs

"""
    AbstractSplitting

Which [`BoundaryChannel`](@ref)s of a boundary the solver traces, see [`splitting_of`](@ref):
[`Dominant`](@ref) or [`AllChannels`](@ref).
"""
abstract type AbstractSplitting end

"""The channel with the largest amplitude factor is traced, the beam is not split."""
struct Dominant <: AbstractSplitting end

"""Every channel that carries light is traced as a child beam."""
struct AllChannels <: AbstractSplitting end

"""
    splitting_of(model::AbstractSurfaceModel)

Returns which channels of a boundary with the surface `model` the solver traces: [`Dominant`](@ref)
(default) or [`AllChannels`](@ref).
"""
splitting_of(::AbstractSurfaceModel) = Dominant()

"""
    select_channels(splitting, cs)

Returns the channels of `cs` that are traced: for [`Dominant`](@ref) the channel with the largest
amplitude factor, or `nothing` if none carries light.
"""
function select_channels(::Dominant, (reflected, transmitted)::Tuple{BoundaryChannel, BoundaryChannel})
    c = abs2(amplitude(reflected)) > abs2(amplitude(transmitted)) ? reflected : transmitted
    iszero(amplitude(c)) && return nothing
    return c
end

"""
    realize(t::Transition, beam::Beam, ray, c)

Returns the interaction of the `ray` that leaves the boundary of the transition `t` in the channel `c`:
a [`BeamInteraction`](@ref) with the new ray and the hint into its medium, or `nothing` without a
channel.
"""
realize(::Transition, ::Beam, ::AbstractRay, ::Nothing) = nothing

function realize(t::Transition, ::Beam{T, R}, ray::R, c::BoundaryChannel) where {T <: Real, R <: AbstractRay{T}}
    return BeamInteraction{T, R}(side_hint(t, side(c)), outgoing_ray(t, ray, c))
end

"""Returns the new ray that leaves the boundary of the transition `t` in the channel `c`."""
outgoing_ray(t::Transition, ray::Ray{T}, c::BoundaryChannel) where {T} =
    Ray{T}(hit_point(ray), direction(c), nothing, wavelength(ray), refractive_index(t, side(c)))

function outgoing_ray(t::Transition, ray::PolarizedRay{T}, c::BoundaryChannel) where {T}
    E = _transverse(c.P * polarization(ray), direction(c))
    return PolarizedRay{T}(hit_point(ray), direction(c), nothing, wavelength(ray),
        refractive_index(t, side(c)), E)
end
