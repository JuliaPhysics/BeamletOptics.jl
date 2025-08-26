"""
    GaussianBeamlet{T} <: AbstractBeam{T, Ray{T}}

Ray representation of the **stigmatic** Gaussian beam as per J. Arnaud (1985). The beam quality `M2` is fully considered via the divergence angle.
The formalism for the beam parameter calculation is based on the following publications:

**Jacques Arnaud, "Representation of Gaussian beams by complex rays," Appl. Opt. 24, 538-543 (1985)**

and

**Donald DeJager and Mark Noethen, "Gaussian beam parameters that use Coddington-based Y-NU paraprincipal ray tracing," Appl. Opt. 31, 2199-2205 (1992)**

# Fields

- `chief`: a [`Beam`](@ref) of [`Ray`](@ref)s to store the chief ray
- `waist`: a [`Beam`](@ref) of [`Ray`](@ref)s to store the waist ray
- `divergence`: a [`Beam`](@ref) of [`Ray`](@ref)s to store the divergence ray
- `λ`: beam wavelength in [m]
- `w0`: local beam waist radius in [m]
- `E0`: complex field value in [V/m]
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
wavelength!(beam::GaussianBeamlet, new) = (beam.λ = new)
beam_waist(beam::GaussianBeamlet) = beam.w0
electric_field(beam::GaussianBeamlet) = beam.E0
electric_field!(beam::GaussianBeamlet{T}, new) where T = (beam.E0 = Complex{T}(new))
refractive_index(beam::GaussianBeamlet, id::Int) = refractive_index(rays(beam.chief)[id])
function refractive_index!(beam::GaussianBeamlet, id::Int, n_new::Real)
    refractive_index!(rays(beam.chief)[id], n_new)
    refractive_index!(rays(beam.waist)[id], n_new)
    refractive_index!(rays(beam.divergence)[id], n_new)
    return nothing
end

Base.length(gauss::GaussianBeamlet) = length(gauss.chief)
optical_path_length(gauss::GaussianBeamlet) = optical_path_length(gauss.chief)

isentering(beam::GaussianBeamlet, id::Int) = isentering(rays(beam.chief)[id])

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

function Base.replace!(gauss::GaussianBeamlet{T},
        interaction::GaussianBeamletInteraction{T},
        index::Int) where {T}
    replace!(gauss.chief, interaction.chief, index)
    replace!(gauss.waist, interaction.waist, index)
    replace!(gauss.divergence, interaction.divergence, index)
    return nothing
end

function _modify_beam_head!(old::GaussianBeamlet{T},
        new::GaussianBeamlet{T}) where {T <: Real}
    _modify_beam_head!(old.chief, new.chief)
    _modify_beam_head!(old.waist, new.waist)
    _modify_beam_head!(old.divergence, new.divergence)
    wavelength!(old, wavelength(new))
    electric_field!(old, electric_field(new))
end

_last_beam_intersection(gauss::GaussianBeamlet) = intersection(last(rays(gauss.chief)))

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
- `support`: [`Nullable`](@ref) support vector for the construction of the waist and div rays

# Additional information

!!! tip "Waist offset"
    The `z0` keyword arg. can be used in order to spawn a beam where the waist is not located at the
    specified `position`, but rather at an offset `z0` in [m] along the chief ray axis.

!!! info "Support vector"
    In order to calculate the basis vectors required for the beamlet construction, a random orthogonal vector is chosen.
    If results fluctuate due to the randomness of this vector, make sure to specify a fixed orthogonal `support` vector.
"""
function GaussianBeamlet(
    position::AbstractArray{P},
    direction::AbstractArray{D}, 
    λ::L = 1e-6,
    w0::Real = 1e-3;
    M2::Real = 1,
    P0::Real = 1e-3,
    z0::Real = 0,
    support::Nullable{AbstractArray} = nothing
    ) where {P<:Real, D<:Real, L<:Real}
    # T = promote_type(P, D, L)
    dir = normalize(direction)
    # Create orthogonal vector for construction purposes (right-handed)
    if isnothing(support)
        s1 = normal3d(dir)
    else
        # FIXME check for orthogonality (see Pol. astigm. beamlet MR)
        s1 = support
    end
    s1 = normalize(s1)
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

global gb_constructor_warn = true

#FIXME remove deprecated constructor until end of 2025
function GaussianBeamlet(chief::Ray{T}, λ=1e-6, w0=1e-3; M2=1, P0=1e-3, support = [0,0,1]) where T
    if gb_constructor_warn
        @warn "The GaussianBeamlet(::Ray, ...) constructor will be deprecated in the future"
        global gb_constructor_warn = false
    end
    s1 = normal3d(direction(chief), support)
    return GaussianBeamlet(
        position(chief),
        direction(chief),
        λ,
        w0;
        M2,
        P0,
        support = s1
    )
end

point_on_beam(gauss::GaussianBeamlet, t::Real) = point_on_beam(gauss.chief, t)

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
- `ψ`: Gouy phase (note that -atan definition is used)
- `w0`: local beam waist radius
"""
function gauss_parameters(
        gauss::GaussianBeamlet,
        z::Real;
        hint = point_on_beam(gauss, z)
    )
    p0, index = hint
    chief = gauss.chief.rays[index]
    div = gauss.divergence.rays[index]
    waist = gauss.waist.rays[index]
    #=
    Divergence ray height and slope (same for waist ray)
    - find divergence ray "height" and "slope" at intersection point y0 with target plane at p0 of chief ray
    - ray height "y_d" is length between p0 and y0
    - ray slope "m_d" is angle between vector p0 -> y0 and divergence ray direction -> gives unambiguous angle for signed ray slope calculation
    - fails if y_d is zero -> catch R=Inf, ψ=0 and H=λ/π
    =#
    intersection_len = line_plane_distance3d(p0,
        direction(chief),
        position(div),
        direction(div))
    y0 = position(div) + intersection_len * direction(div) - p0
    y_d = norm(y0)
    y0 /= y_d
    m_d = tan(π / 2 - angle3d(y0, direction(div)))
    # Waist ray height and slope
    intersection_len = line_plane_distance3d(p0,
        direction(chief),
        position(waist),
        direction(waist))
    y0 = position(waist) + intersection_len * direction(waist) - p0
    y_w = norm(y0)
    y0 /= y_w
    m_w = tan(π / 2 - angle3d(y0, direction(waist)))
    # Beam parameters as per Arnaud (1985) and DeJager (1992)
    n = refractive_index(chief)
    H = abs(n * (y_w * m_d - y_d * m_w))
    # Test optical invariant
    λ = wavelength(gauss)
    if !isapprox(H, λ / π, atol = 1e-6)
        H = λ / π
        # println("H not fulfilled at z=$z")
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
    return w, R, ψ, w0
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
This function also considers phase changes due to changes in the [`optical_path_length`](@ref) of the beamlet.

!!! warning
    Note that `z` and `r` must be specified as cartesian distances. Using the optical path length for `z` can lead to false results.
"""
function electric_field(gauss::GaussianBeamlet, r, z)
    point, index = point_on_beam(gauss, z)
    w, R, ψ, w0 = gauss_parameters(gauss, z, hint = (point, index))
    k = wave_number(wavelength(gauss))
    # Calculate new local field strength based on E0*w0 = const.
    E0 = electric_field(gauss) * (beam_waist(gauss) / w0)
    # Calculate phase change due to optical path length
    Δl = optical_path_length(gauss) - length(gauss)
    # Note: geometrical length changes considered in `electric_field` call below 
    ref_ϕ = Δl / wavelength(gauss) * 2π
    return electric_field(r, z, E0, w0, w, k, ψ, R) * exp(im*ref_ϕ)
end

optical_power(gauss::GaussianBeamlet) = intensity(electric_field(gauss)) / 2 * π * beam_waist(gauss)^2

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

function isparentbeam(beam::GaussianBeamlet, ray::AbstractRay)
    c = isparentbeam(beam.chief, ray)
    w = isparentbeam(beam.waist, ray)
    d = isparentbeam(beam.divergence, ray)
    return any((c, w, d))
end

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
