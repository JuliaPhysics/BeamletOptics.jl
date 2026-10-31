"""
    GaussianBeamlet{T} <: AbstractBeam{T, Ray{T}}

Ray representation of the **stigmatic** Gaussian beam as per J. Arnaud (1985). The beam quality `M2` is fully considered via the divergence angle.
The formalism for the beam parameter calculation is based on the following publications:

**Jacques Arnaud, "Representation of Gaussian beams by complex rays." Appl. Opt. 24, 538-543 (1985)**

and

**Donald DeJager and Mark Noethen, "Gaussian beam parameters that use Coddington-based Y-NU paraprincipal ray tracing," Appl. Opt. 31, 2199-2205 (1992)**

# Fields

- `chief`: a [`Beam`](@ref) of [`Ray`](@ref)s to store the chief ray
- `waist`: a [`Beam`](@ref) of [`Ray`](@ref)s to store the waist ray
- `divergence`: a [`Beam`](@ref) of [`Ray`](@ref)s to store the divergence ray
- `λ`: beam wavelength in [m]
- `w0`: local beam waist radius in [m]
- `E0`: complex reference amplitude in [V/m] of the field, normalized to the power of the beamlet: `P = |E0|² π w0² / (4 Z₀)`
  for any medium. The physical field of a segment in a medium of refractive index `n` is `E0 / √n` (times the
  profile), so that the irradiance `n |E|² / (2 Z₀)` of a beam of fixed power does not depend on the medium, see [`electric_field`](@ref)
- `parent`: reference to the parent beam, if any ([`Nullable`](@ref) to account for the root beam which has no parent)
- `children`: vector of child beams, each child beam represents a branching or bifurcation of the original beam, i.e. beam-splitting

# Additional information

!!! info "Beam parameters"
    Parameters of the beam, e.g. ``w(z)`` or ``R(z)``, can be obtained through the [`gauss_parameters`](@ref) function.

!!! warning "Astigmatism and abberations"
    It is assumed, but not forbidden, that the optical system contains non-flat or non-parabolic beam-surface-interactions that cause the beam to obtain
    astigmatism or higher-order abberations. These can not be represented by the `GaussianBeamlet`.
"""
mutable struct GaussianBeamlet{T} <: AbstractBeam{T, Ray{T}}
    chief::Beam{T, Ray{T}}
    waist::Beam{T, Ray{T}}
    divergence::Beam{T, Ray{T}}
    λ::T
    w0::T
    E0::Complex{T}
    parent::Nullable{GaussianBeamlet{T}}
    children::Vector{GaussianBeamlet{T}}
end

function GaussianBeamlet(chief::Beam{T, Ray{T}},
        waist::Beam{T, Ray{T}},
        div::Beam{T, Ray{T}},
        λ::T,
        w0::T,
        E0::Complex{T}) where {T <: Real}
    return GaussianBeamlet{T}(
        chief,
        waist,
        div,
        λ,
        w0,
        E0,
        nothing,
        Vector{GaussianBeamlet{T}}())
end

"""
    GaussianBeamletInteraction <: AbstractInteraction

This type is used to store the new beamlet section resulting from on optical interaction
between a [`GaussianBeamlet`](@ref) and some [`AbstractObject`](@ref).
Uses the hint of the `chief` beam.

# Fields

- `chief`: [`Beam`](@ref) interaction
- `waist`: [`Beam`](@ref) interaction
- `divergence`: [`Beam`](@ref) interaction
"""
struct GaussianBeamletInteraction{R <: Real} <: AbstractInteraction
    chief::BeamInteraction{R}
    waist::BeamInteraction{R}
    divergence::BeamInteraction{R}
end

hint(i::GaussianBeamletInteraction) = hint(i.chief)
hint!(i::GaussianBeamletInteraction, new_hint::Nullable{Hint}) = hint!(i.chief, new_hint)

wavelength(beam::GaussianBeamlet) = beam.λ
beam_waist(beam::GaussianBeamlet) = beam.w0
electric_field(beam::GaussianBeamlet) = beam.E0
electric_field!(beam::GaussianBeamlet{T}, new) where {T} = (beam.E0 = Complex{T}(new))

"""Return a tuple of the `chief`, `waist` and `divergence` beams of the [`GaussianBeamlet`](@ref)."""
_component_beams(g::GaussianBeamlet) = (g.chief, g.waist, g.divergence)

"""
    parent!(beam::GaussianBeamlet, parent::GaussianBeamlet)

Ensures that the GaussianBeamlet knows about its parent beam. In addition, links the chief beams of child and parent.
Important for correct functioning of [`point_on_beam`](@ref) and [`length`](@ref).
"""
function parent!(child::GaussianBeamlet, parent::GaussianBeamlet)
    child.parent = parent
    child.chief.parent = parent.chief
    return nothing
end

"""
    interact3d(system::AbstractSystem, object::AbstractObject, gauss::GaussianBeamlet{R}, ray_id::Int)

Generic dispatch for the [`interact3d`](@ref) method of a [`GaussianBeamlet`](@ref) with an [`AbstractObject`](@ref).
Unless a more concrete implementation exists, the interaction of the Gaussian is assumed to be the interaction of the
chief, waist and divergence rays with an object.

# Returns

The `interact3d` method for the [`GaussianBeamlet`](@ref) must return a `GaussianBeamletInteraction`.
"""
function interact3d(system::AbstractSystem,
        object::AbstractObject,
        gauss::GaussianBeamlet{R},
        ray_id::Int) where {R}
    i_c = interact3d(system, object, gauss.chief, rays(gauss.chief)[ray_id])
    i_w = interact3d(system, object, gauss.waist, rays(gauss.waist)[ray_id])
    i_d = interact3d(system, object, gauss.divergence, rays(gauss.divergence)[ray_id])
    if any(isnothing, (i_c, i_w, i_d))
        return nothing
    end
    return GaussianBeamletInteraction{R}(i_c, i_w, i_d)
end

function Base.push!(gauss::GaussianBeamlet{T},
        interaction::GaussianBeamletInteraction{T}) where {T}
    push!(gauss.chief, interaction.chief)
    push!(gauss.waist, interaction.waist)
    push!(gauss.divergence, interaction.divergence)
    return nothing
end

"""
    _beams_hits_same_shape(gauss, id)

Tests if all rays at section `id` of `gauss` hit the same object shape.
Returns `true` or `false`.
"""
@inline function _beams_hits_same_shape(gauss::GaussianBeamlet, id::Int)::Bool
    c = intersection(rays(gauss.chief)[id])
    w = intersection(rays(gauss.waist)[id])
    d = intersection(rays(gauss.divergence)[id])
    are_nothing = isnothing.((c, w, d))
    if any(are_nothing)
        return all(are_nothing)
    end
    return shape(c) === shape(w) === shape(d)
end

"""
    GaussianBeamlet(position, direction, λ, w0; kwargs...)

Constructs a Gaussian beamlet at its waist with the specified beam parameters.

# Arguments

The following inputs and arguments can be used to configure the beamlet:

## Inputs

- `position`: origin of the beamlet
- `direction`: direction of the beamlet
- `λ`: wavelength of the beamlet in [m]. Default value is 1000 nm.
- `w0`: beam waist (radius) in [m]. Default value is 1 mm.

## Keyword Arguments

- `M2`: beam quality factor. Default is 1
- `P0`: beam total power in [W]. Default is 1 mW
- `z0`: beam waist offset in [m]. Default is 0 m
- `support`: [`Nullable`](@ref) support vector for the construction of the waist and div rays. Must be
  orthogonal to `direction`, has any non-zero length and is normalized internally. Default is `nothing`,
  i.e. a vector orthogonal to `direction` is chosen.

# Additional information

!!! tip "Waist offset"
    The `z0` keyword arg. can be used in order to spawn a beam where the waist is not located at the
    specified `position`, but rather at an offset `z0` in [m] along the chief ray axis.

!!! info "Support vector"
    In order to calculate the basis vectors required for the beamlet construction, a random orthogonal vector is chosen.
    If results fluctuate due to the randomness of this vector, make sure to specify a fixed orthogonal `support` vector.
    The waist and divergence rays start from the chief ray along this vector. A `support` that is not
    orthogonal to the `direction` (the dot product of the unit vectors exceeds
    `get_orthogonality_threshold()`, default `1e-10`) throws an `ArgumentError`, since it would
    give a beamlet with a smaller waist radius and a shifted waist.
"""
function GaussianBeamlet(
        position::AbstractArray{P},
        direction::AbstractArray{D},
        λ::L = get_default_wavelength(),
        w0::Real = get_default_waist();
        M2::Real = 1,
        P0::Real = get_default_power(),
        z0::Real = 0,
        support::Nullable{AbstractArray} = nothing
) where {P <: Real, D <: Real, L <: Real}
    # T = promote_type(P, D, L)
    dir = normalize(direction)
    # Create orthogonal vector for construction purposes (right-handed)
    if isnothing(support)
        s1 = normal3d(dir)
    else
        s1 = support
    end
    s1 = normalize(s1)
    # Unit vectors: test independent of support length
    if !isorthogonal3d(dir, s1; atol = Config.get_orthogonality_threshold())
        throw(ArgumentError("the support vector must be orthogonal to the direction (dot product of the unit vectors: $(dot(dir, s1)))"))
    end
    # Divergence angle in rad
    tanθ = tan(divergence_angle(λ, w0, M2))
    # Waist ray
    wst = Ray(position + s1 * w0, dir, λ)
    # Divergence ray (Δz determines waist pos. at z0 along optical axis)
    div_dir = normalize(dir + s1 * tanθ)
    Δz = -z0 * tanθ
    div = Ray(position + s1 * Δz, div_dir, λ)
    # Chief ray
    chf = Ray(position, dir, λ)
    # Calculate E0 based on P0, assume zero initial phase offset
    I0 = 2 * P0 / (π * w0^2)
    E0 = electric_field(I0)
    return GaussianBeamlet(
        Beam(chf),
        Beam(wst),
        Beam(div),
        λ,
        w0,
        E0
    )
end

"""
    gauss_parameters(gauss::GaussianBeamlet, z; hint::Union{Nothing, Tuple{Int, Vector{<:Real}}}=nothing)

Calculate the local waist radius and Gouy phase of an unastigmatic Gaussian beamlet at a specific cartesian distance `z` based on the method of J. Arnaud (1985) and D. DeJager (1992).

# Arguments

- `gauss`: the GaussianBeamlet object for which parameters are to be calculated.
- `z`: the position along the beam at which to calculate the parameters.
- `hint`: an optional hint parameter for the relevant point/index of the appropriate beam segment. If not provided, the function will automatically select the ray.

# Returns

- `w`: local radius
- `R`: curvature, i.e. 1/r where r is the radius of curvature
- `ψ`: Gouy phase (note that -atan definition is used), collected along the beam from the start
  of the root beam: continuous at every surface and from a parent beam to its children, i.e. the
  Gouy phase that [`electric_field`](@ref) contains. The same holds for an
  [`AstigmaticGaussianBeamlet`](@ref).
- `w0`: local beam waist radius
"""
gauss_parameters(gauss::GaussianBeamlet, z::Real; hint = point_on_beam(gauss, z)) =
    gauss_parameters(gauss, z, _gouy_before(gauss, last(hint)); hint)

"""
    gauss_parameters(gauss::GaussianBeamlet, z, ψ0; hint)

As [`gauss_parameters`](@ref)`(gauss, z)`, but the returned Gouy phase is the one relative to the
waist of the local segment plus `ψ0`. With `ψ0 = 0` it is the local Gouy phase, which jumps at every
surface, e.g. at a lens; `gauss_parameters(gauss, z)` passes the Gouy phase that the beam has
collected in front of the segment.
"""
function gauss_parameters(
        gauss::GaussianBeamlet,
        z::Real,
        ψ0::Real;
        hint = point_on_beam(gauss, z)
)
    p0, index = hint
    chief = rays(gauss.chief)[index]
    w, R, ψ, w0 = gauss_parameters(GaussianSegment(gauss, index), dot(p0 - position(chief), direction(chief)))
    return (w, R, ψ + ψ0, w0)
end

"""
    RayCrossing{T}

Offset `y0 + l * v` from the chief ray at distance `l` and direction `a` of an auxiliary ray.
"""
struct RayCrossing{T}
    y0::Point3{T}
    v::Point3{T}
    a::Point3{T}
end

function RayCrossing{T}(p::AbstractVector, d::AbstractVector, ray::AbstractRay) where {T}
    a = direction(ray)
    c = dot(d, a)
    return RayCrossing{T}(position(ray) + dot(p - position(ray), d) / c * a - p, a / c - d, a)
end

"""
    GaussianSegment{T}

The waist and divergence rays of one segment of a [`GaussianBeamlet`](@ref), as two `RayCrossing`s
relative to its chief ray. Both rays are straight within the segment, hence the crossings hold exactly.
`n` is the refractive index of the segment and `λ` the vacuum wavelength.
[`gauss_parameters`](@ref)`(segment, l)` evaluates it, a [`GaussianBeamletHit`](@ref) stores it.
"""
struct GaussianSegment{T}
    waist::RayCrossing{T}
    divergence::RayCrossing{T}
    n::T
    λ::T
end

function GaussianSegment(gauss::GaussianBeamlet{T}, index::Integer) where {T}
    chief = rays(gauss.chief)[index]
    p, d = position(chief), direction(chief)
    waist, divergence = map(b -> RayCrossing{T}(p, d, rays(b)[index]), (gauss.waist, gauss.divergence))
    return GaussianSegment{T}(waist, divergence, refractive_index(chief), wavelength(gauss))
end

# Height and slope (cotangent of the angle between offset and direction) of an aux ray at distance `l`
function height_slope(c::RayCrossing, l)
    y = c.y0 + l * c.v
    h = norm(y)
    ŷ = y / h
    return h, dot(ŷ, c.a) / norm(cross(ŷ, c.a))
end

"""
    gauss_parameters(segment::GaussianSegment, l)

Returns `(w, R, ψ, w0)` of the [`GaussianSegment`](@ref) at the distance `l` behind the start of
its chief ray, as [`gauss_parameters`](@ref)`(gauss, z)` does, but with the Gouy phase `ψ` relative
to the local waist of the segment.
"""
function gauss_parameters(s::GaussianSegment, l::Real)
    # Height and slope of aux rays
    (y_w, m_w), (y_d, m_d) = map(c -> height_slope(c, l), (s.waist, s.divergence))
    # Beam parameters as per Arnaud (1985) and DeJager (1992)
    n = s.n
    H = abs(n * (y_w * m_d - y_d * m_w))
    # Test optical invariant
    λ = s.λ
    if !isapprox(H, λ / π, atol = Config.get_invariant_threshold())
        H = λ / π
    end
    E_kt = y_d * m_d + y_w * m_w
    F_kt = sqrt(m_d^2 + m_w^2)
    w = sqrt(y_d^2 + y_w^2)
    R = E_kt / w^2
    z = E_kt / F_kt^2
    ψ = -atan(1, √(1 / (R * z) - 1))
    w0 = H / (n * F_kt)
    # Catch NaNs and correct Gouy phase sign based on curvature sign
    isnan(R) && (R = zero(R))
    isnan(ψ) && (ψ = zero(ψ))
    isnan(w0) && (w0 = w)
    R < 0 && (ψ = -ψ)
    return (w, R, ψ, w0)
end

"""
    waist_parameters(gb, z)

Compute the position and cross-section axes at distance `z` along the beam.
For a stigmatic `GaussianBeamlet`, the cross-section is circular, and the axes are
determined by the orientation of the auxiliary waist ray.
Returns the vectors `(p0, w1, w2)` and scalar waist radii `(w0, w0)`.
"""
function waist_parameters(gb::GaussianBeamlet, z::Real)
    p0, i = point_on_beam(gb, z)
    w, _, _, _ = gauss_parameters(gb, z, hint = (p0, i))
    chief = rays(gb.chief)[i]
    waist = rays(gb.waist)[i]
    dir = direction(chief)
    # Basis vector from aux waist ray
    vec = position(waist) - p0
    # Project to transverse plane to be sure
    vec = vec - dot(vec, dir) * dir
    if norm(vec) < 1e-12
        s1 = normal3d(dir)
    else
        s1 = normalize(vec)
    end
    s2 = cross(dir, s1)
    # Get waist size for unification
    _, _, _, w0 = gauss_parameters(gb, z)
    return p0, w * s1, w * s2, w0, w0
end

function waist_parameters(gauss::GaussianBeamlet{G}, zs::AbstractArray) where {G}
    n = length(zs)
    p0 = Vector{Point3{G}}(undef, n)
    w1 = Vector{Point3{G}}(undef, n)
    w2 = Vector{Point3{G}}(undef, n)
    w01 = Vector{G}(undef, n)
    w02 = Vector{G}(undef, n)
    @inbounds for i in 1:n
        p0[i], w1[i], w2[i], w01[i], w02[i] = waist_parameters(gauss, zs[i])
    end
    return (p0, w1, w2, w01, w02)
end

"""
    gauss_parameters(gauss::GaussianBeamlet, zs::AbstractArray)

Return the parameters of the `GaussianBeamlet` along the specified positions in `zs`.
"""
function gauss_parameters(gauss::GaussianBeamlet{G}, zs::AbstractArray) where {G}
    n = length(zs)
    w = Vector{G}(undef, n)
    R = Vector{G}(undef, n)
    ψ = Vector{G}(undef, n)
    w0 = Vector{G}(undef, n)
    @inbounds for i in 1:n
        w[i], R[i], ψ[i], w0[i] = gauss_parameters(gauss, zs[i])
    end
    return (w, R, ψ, w0)
end

"""
    electric_field(gauss::GaussianBeamlet, r, z)

Calculates the electric field phasor [V/m] of the [`GaussianBeamlet`](@ref) at the radial and longitudinal positions `r` and `z`.
It is the physical field, like that of an [`AstigmaticGaussianBeamlet`](@ref): with the irradiance `n |E|² / (2 Z₀)`
of a medium of refractive index `n`, a beam of fixed power has the amplitude `∝ 1/√n`, hence the field of a
segment in a medium is the reference amplitude `E0` of the beamlet divided by `√n`.
The phase contains the [`optical_path_length`](@ref) from the start of the beam (parents included) up to `z`: within a medium
of refractive index `n` it advances with the wavenumber `n k₀`, along the beam and across it. The Gouy phase is collected
along the beam as well, i.e. it is continuous at every surface and from a parent beam to its children, and equals the `ψ`
of [`gauss_parameters`](@ref).

!!! warning
    Note that `z` and `r` must be specified as cartesian distances. Using the optical path length for `z` can lead to false results.
"""
function electric_field(gauss::GaussianBeamlet, r, z; hint = point_on_beam(gauss, z))
    point, index = hint
    return segment_field(gauss, r, z, point, index, length_before(gauss, index),
        path_excess_before(gauss, index), _gouy_before(gauss, index))
end

# Gouy phase relative to the local waist at the start and at the end of segment `index`
local_gouy_start(gauss::GaussianBeamlet{T}, index::Integer) where {T} =
    gauss_parameters(gauss, zero(T), zero(T); hint = (position(rays(gauss.chief)[index]), index))[3]
function local_gouy_end(gauss::GaussianBeamlet{T}, index::Integer) where {T}
    ray = rays(gauss.chief)[index]
    return gauss_parameters(gauss, zero(T), zero(T);
        hint = (hit_point(ray), index))[3]
end

"""
    _gouy_before(gauss::GaussianBeamlet, index)

Gouy phase collected in front of segment `index`, parents included.
"""
function _gouy_before(gauss::GaussianBeamlet{T}, index::Integer) where {T}
    p = gauss.parent
    ψ0 = zero(T)
    if !isnothing(p)
        last_id = length(rays(p.chief))
        ψ0 = _gouy_before(p, last_id) + local_gouy_end(p, last_id) - local_gouy_start(gauss, 1)
    end
    for i in 2:index
        ψ0 += local_gouy_end(gauss, i - 1) - local_gouy_start(gauss, i)
    end
    return ψ0
end

# Field of segment `index` at distance `z` and radius `r` from the length, path and Gouy offsets in front of it
function segment_field(gauss::GaussianBeamlet, r, z, point, index, l0, Δl0, ψ0)
    w, R, ψ, w0 = gauss_parameters(gauss, z, ψ0; hint = (point, index))
    k0 = wavenumber(wavelength(gauss))
    n = refractive_index(gauss, index)
    # Local field strength from E0 * w0 = const, reduced by √n
    E0 = electric_field(gauss) * (beam_waist(gauss) / w0) / sqrt(n)
    # Optical path up to the point
    opl = z + Δl0 + (n - 1) * (z - l0)
    return E0 * w0 / w * exp(-r^2 / w^2) * cis(k0 * opl + ψ + n * k0 * r^2 * R / 2)
end

"""
    optical_power(gauss::GaussianBeamlet)

Power of the beamlet in [W] at its start. The reference amplitude `E0` is normalized to the power, hence
the power does not depend on the medium of the beamlet.
"""
function optical_power(gauss::GaussianBeamlet)
    intensity(electric_field(gauss)) / 2 * π * beam_waist(gauss)^2
end

"""
    isparaxial(system, gb::GaussianBeamlet, threshold=π/4)

Tests the angle between the waist and divergence beams and refractive surfaces.
A target threshold of π/4 or 45° is assumed before abberations become dominant.
"""
isparaxial(system::AbstractSystem, gb::GaussianBeamlet, threshold::Real = π / 4) = isparaxial(
    system,
    gb.waist,
    threshold) & isparaxial(
    system,
    gb.divergence,
    threshold)

"""
    istilted(system::System, gb::GaussianBeamlet)

Tests if refractive elements are tilted with respect to the beamlet optical axis, i.e. introduce simple astigmatism.
"""
istilted(system::AbstractSystem, gb::GaussianBeamlet) = !isparaxial(system, gb.chief, 0)

"""
    rayleigh_range(g::GaussianBeamlet; M2=1)

Returns the Rayleigh range for the **first** beam section of the [`GaussianBeamlet`](@ref) `g`.
Note: `M2` is not stored in `g` during construction and must be specified by the user.
"""
function rayleigh_range(g::GaussianBeamlet; M2 = 1)
    λ = wavelength(g)
    w0 = beam_waist(g)
    return rayleigh_range(λ, w0, M2)
end
