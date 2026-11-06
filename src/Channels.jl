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

"""
    transfer(c::BoundaryChannel, P, a)
    transfer(c::BoundaryChannel, ray, J::LocalJonesBasis, a)

Returns a copy of the channel `c` with the field transfer `P` and the amplitude factor `a`. With a
Jones matrix `J` in the s- and p-basis of the `ray` and the channel, e.g. `SPBasis(rs, 0, 0, rp)`, the
field transfer is that of `J` in global coordinates, or `nothing` for a ray without a
field. A surface model builds its channels with it.
"""
transfer(c::BoundaryChannel, P, a::Number) = BoundaryChannel(side(c), direction(c), P, a)
transfer(c::BoundaryChannel, ray::AbstractRay, J::LocalJonesBasis, a::Number) =
    transfer(c, jones_transfer(ray, direction(c), J), a)

# The global field transfer of the Jones matrix `J` from the `ray` into the direction `dir`
jones_transfer(::AbstractRay, dir, ::LocalJonesBasis) = nothing
jones_transfer(ray::PolarizedRay, dir, J::LocalJonesBasis) =
    _calculate_global_E0(direction(ray), dir, normal3d(intersection(ray)), J)

"""The refractive index on the `side` of the boundary of the transition `t`."""
refractive_index(t::Transition, ::Reflected) = refractive_index_in(t)
refractive_index(t::Transition, ::Transmitted) = refractive_index_out(t)

"""The hint into the medium on the `side` of the boundary of the transition `t`."""
side_hint(t::Transition, ::Reflected) = hint_in(t)
side_hint(t::Transition, ::Transmitted) = hint_out(t)

"""
    channels(t::Transition, ray)

Returns the channels of the bare boundary of the transition `t` for the `ray`, as the tuple
`(transmitted, reflected)` of [`BoundaryChannel`](@ref)s, the order of the `children` of a split beam.
The `ray` is refracted from the medium before into the medium behind the boundary via
[`refraction3d`](@ref); beyond the critical angle it is totally reflected, and the transmitted channel
carries no light. The field of a [`PolarizedRay`](@ref) is
transferred with the [`fresnel_coefficients`](@ref). Between equal indices the ray passes unchanged.
"""
function channels(t::Transition, ray::Ray)
    d = direction(ray)
    n1, n2 = refractive_index_in(t), refractive_index_out(t)
    reflected = reflection3d(d, normal3d(t))
    n1 == n2 && return (BoundaryChannel(Transmitted(), d, nothing, 1),
        BoundaryChannel(Reflected(), reflected, nothing, 0))
    refracted, TIR = refraction3d(d, normal3d(t), n1, n2)
    TIR && return (BoundaryChannel(Transmitted(), reflected, nothing, 0),
        BoundaryChannel(Reflected(), reflected, nothing, 1))
    return (BoundaryChannel(Transmitted(), refracted, nothing, 1),
        BoundaryChannel(Reflected(), reflected, nothing, 0))
end

function channels(t::Transition, ray::PolarizedRay{T}) where {T}
    d = direction(ray)
    normal = normal3d(t)
    n1, n2 = refractive_index_in(t), refractive_index_out(t)
    reflected = reflection3d(d, normal)
    O = zero(SMatrix{3, 3, Complex{T}, 9})
    if n1 == n2
        P = SMatrix{3, 3, Complex{T}, 9}(I)
        return (BoundaryChannel(Transmitted(), d, P, 1),
            BoundaryChannel(Reflected(), reflected, O, 0))
    end
    # Fresnel coefficients from the cosine, smooth for AD at normal incidence
    rs, rp, ts, tp = _fresnel_coefficients(cos_angle3d(d, -normal), n2 / n1)
    Pr = _calculate_global_E0(d, reflected, normal3d(intersection(ray)), SPBasis(-rs, 0, 0, rp))
    if is_internally_reflected(rp, rs)
        return (BoundaryChannel(Transmitted(), reflected, O, 0),
            BoundaryChannel(Reflected(), reflected, Pr, 1))
    end
    refracted, _ = refraction3d(d, normal, n1, n2)
    Pt = _calculate_global_E0(d, refracted, normal3d(intersection(ray)), SPBasis(ts, 0, 0, tp))
    return (BoundaryChannel(Transmitted(), refracted, Pt, 1),
        BoundaryChannel(Reflected(), reflected, Pr, 0))
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
amplitude factor, or `nothing` if none carries light; for [`AllChannels`](@ref) all of `cs`, of which
[`realize`](@ref) traces those that carry light.
"""
function select_channels(::Dominant, cs::Tuple)
    c = reduce((a, b) -> abs2(amplitude(a)) ≥ abs2(amplitude(b)) ? a : b, cs)
    return islit(c) ? c : nothing
end

select_channels(::AllChannels, cs::Tuple) = cs

"""Returns `true` if the channel `c` carries light."""
islit(c::BoundaryChannel) = !iszero(amplitude(c))

# The child beams of the tuple `kids` whose channel carries light, `lit`, and whose power is not below
# the power cutoff, as a vector of type `B`
function lit_children(::Type{B}, kids::Tuple, lit::Tuple, root::AbstractBeam) where {B}
    traced = map((kid, l) -> l && !(power_fraction(kid, root) < get_power_cutoff()), kids, lit)
    v = Vector{B}(undef, count(traced))
    k = 0
    for i in eachindex(kids)
        traced[i] || continue
        v[k += 1] = kids[i]
    end
    return v
end

"""
    power_fraction(beam, root)

The fraction of the power of the `root` beam that the child `beam` carries at its start. Beams of
[`Ray`](@ref)s carry no power, their fraction is one.
"""
power_fraction(::Beam{T, <:Ray}, ::AbstractBeam) where {T} = one(T)
power_fraction(beam::Beam{T, <:PolarizedRay}, root::Beam) where {T} =
    ray_flux(first(rays(beam))) / ray_flux(first(rays(root)))
power_fraction(beamlet::GaussianBeamlets, root::GaussianBeamlets) =
    optical_power(beamlet) / optical_power(root)

# Power flux of a polarized ray, up to a constant factor
ray_flux(ray::PolarizedRay) = refractive_index(ray) * sum(abs2, polarization(ray))

"""
    realize(t::Transition, beam::Beam, ray, c)

Returns the interaction of the `ray` that leaves the boundary of the transition `t` in the channel `c`:
a [`BeamInteraction`](@ref) with the new ray and the hint into its medium, or `nothing` without a
channel.
"""
realize(::Transition, ::Beam, ::AbstractRay, ::Nothing) = nothing

function realize(t::Transition, ::Beam{T, R}, ray::R,
        c::BoundaryChannel) where {T <: Real, R <: AbstractRay{T}}
    return BeamInteraction{T, R}(side_hint(t, side(c)), outgoing_ray(t, ray, c))
end

# Several channels: each that carries light starts a child beam, unless its power lies below the power
# cutoff; a single one continues the beam
function realize(t::Transition, beam::Beam{T, R}, ray::R,
        cs::Tuple) where {T <: Real, R <: AbstractRay{T}}
    num_lit = count(islit, cs)
    num_lit == 0 && return nothing
    num_lit == 1 && return realize(t, beam, ray, select_channels(Dominant(), cs))
    kids = map(c -> Beam(outgoing_ray(t, ray, c)), cs)
    children!(beam, lit_children(Beam{T, R}, kids, map(islit, cs), AbstractTrees.getroot(beam)))
    return nothing
end

"""Returns the new ray that leaves the boundary of the transition `t` in the channel `c`."""
outgoing_ray(t::Transition, ray::Ray{T}, c::BoundaryChannel) where {T} =
    Ray{T}(hit_point(ray), direction(c), nothing, wavelength(ray), refractive_index(t, side(c)))

function outgoing_ray(t::Transition, ray::PolarizedRay{T}, c::BoundaryChannel) where {T}
    E = _transverse(c.P * polarization(ray), direction(c))
    return PolarizedRay{T}(hit_point(ray), direction(c), nothing, wavelength(ray),
        refractive_index(t, side(c)), E)
end

"""
    beamlet_interaction(system, splitting, model, beamlet, ray_id)

The interaction of a [`GaussianBeamlet`](@ref) or an [`AstigmaticGaussianBeamlet`](@ref) with the
surface `model` at the segment `ray_id`. For [`Dominant`](@ref), each component ray interacts with the
model at its own hit. For [`AllChannels`](@ref), the channels of the chief ray decide: each channel
that carries light starts a child beamlet, whose component rays leave their own hits on the same side;
a single one continues the beamlet. The field factor of a child is that of the channel of the chief
ray, its waist is the waist of the `beamlet` at the boundary.
"""
function beamlet_interaction(system::AbstractSystem, ::Dominant, model::AbstractSurfaceModel,
        beamlet::GaussianBeamlets, ray_id::Int)
    return beamlet_interaction(beamlet) do b
        ray = rays(b)[ray_id]
        interact3d(system, model, transition(system, ray), b, ray)
    end
end

function beamlet_interaction(system::AbstractSystem, ::AllChannels, model::AbstractSurfaceModel,
        beamlet::GaussianBeamlets, ray_id::Int)
    # New rays and amplitude factors per component and channel
    hits = map(_component_beams(beamlet)) do b
        ray = rays(b)[ray_id]
        t = transition(system, ray)
        cs = channels(model, t, ray, channels(t, ray))
        (map(c -> outgoing_ray(t, ray, c), cs), map(amplitude, cs))
    end
    a = last(first(hits))
    lit = map(!iszero, a)
    num_lit = count(lit)
    num_lit == 0 && return nothing
    num_lit == 1 && return continued_beamlet(system, model, beamlet, ray_id, findfirst(lit))
    kids = ntuple(i -> child_beamlet(beamlet, map(h -> first(h)[i], hits), a[i]), Val(length(a)))
    children!(beamlet, lit_children(typeof(beamlet), kids, lit, AbstractTrees.getroot(beamlet)))
    return nothing
end

# A single channel of the chief ray carries light, e.g. under total internal reflection: each
# component ray continues in this channel at its own hit
function continued_beamlet(system::AbstractSystem, model::AbstractSurfaceModel,
        beamlet::GaussianBeamlets, ray_id::Int, i::Int)
    return beamlet_interaction(beamlet) do b
        ray = rays(b)[ray_id]
        t = transition(system, ray)
        realize(t, b, ray, channels(model, t, ray, channels(t, ray))[i])
    end
end

"""
    child_beamlet(beamlet, starts, a)

Returns the child of the `beamlet` whose component beams start with the rays `starts`, with the
amplitude factor `a` of the channel of its chief ray. The waist of a [`GaussianBeamlet`](@ref) is its waist at
the boundary.
"""
function child_beamlet(gauss::GaussianBeamlet, starts::Tuple, a::Number)
    chief, waist, divergence = map(Beam, starts)
    w0 = gauss_parameters(gauss, length(gauss))[4]
    E0 = a * electric_field(gauss) * (beam_waist(gauss) / w0)
    return GaussianBeamlet(chief, waist, divergence, wavelength(gauss), w0, E0)
end

child_beamlet(::AstigmaticGaussianBeamlet, starts::Tuple, ::Number) =
    AstigmaticGaussianBeamlet(map(Beam, starts)...)
