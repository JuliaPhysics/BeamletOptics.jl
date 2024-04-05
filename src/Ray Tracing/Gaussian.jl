"""
    GaussianBeamlet{T} <: AbstractBeam{T, Ray{T}}

Ray representation of the **stigmatic** Gaussian beam as per J. Arnaud (1985). The beam quality `M2` is fully considered via the divergence angle. 
Formalism for beam parameter calculation based on publications:

**Jacques Arnaud, "Representation of Gaussian beams by complex rays," Appl. Opt. 24, 538-543 (1985)**

and

**Donald DeJager and Mark Noethen, "Gaussian beam parameters that use Coddington-based Y-NU paraprincipal ray tracing," Appl. Opt. 31, 2199-2205 (1992)**

# Fields

- `id`: beam ID (uuid4)
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
    It is assumed, but not forbidden, that the optical system contains non-symmetric optical elements that cause the beam to obtain
    astigmatism or higher-order abberations. These can not be represented by the `GaussianBeamlet`.
    Refer to **FIXME** for more information.
"""
mutable struct GaussianBeamlet{T} <: AbstractBeam{T, Ray{T}}
    id::UUID
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
    return GaussianBeamlet{T}(uuid4(),
        chief,
        waist,
        div,
        λ,
        w0,
        E0,
        nothing,
        Vector{GaussianBeamlet{T}}())
end

struct GaussianBeamletInteraction{R <: Real} <: AbstractInteraction
    chief::BeamInteraction{R}
    waist::BeamInteraction{R}
    divergence::BeamInteraction{R}
end

hint(interaction::GaussianBeamletInteraction) = hint(interaction.chief)

wavelength(beam::GaussianBeamlet) = beam.λ
beam_waist(beam::GaussianBeamlet) = beam.w0
beam_amplitude(beam::GaussianBeamlet) = beam.E0

Base.length(gauss::GaussianBeamlet; opl::Bool=false) = length(gauss.chief; opl)

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
end

_last_beam_intersection(gauss::GaussianBeamlet) = intersection(last(rays(gauss.chief)))

"""
    GaussianBeamlet(chief::Ray{T}, λ=1000e-9, w0=1e-3; M2=1, P0=1e-3, support=[0, 0, 1])

Construct a Gaussian beamlet at its waist with a specified beam diameter.

# Arguments

- `chief`: Chief ray defining the origin of the Gaussian beamlet.
- `λ`: Wavelength of the beamlet in [m]. Default value is 1000 nm.
- `w0`: Beam waist (radius) in [m]. Default value is 1 mm.

# Keyword Arguments

- `M2`: beam quality factor
- `P0`: beam total power in [W]
- `support`: Support vector that can be adjusted for beamlet construction.
"""
function GaussianBeamlet(chief::Ray{T},
        λ = 1000e-9,
        w0 = 1e-3;
        M2 = 1,
        P0 = 1e-3,
        support = [0, 0, 1]) where {T}
    # Create orthogonal vector for construction purposes
    bob_the_builder = normal3d(direction(chief), support)
    # Divergence angle in rad
    θ = divergence_angle(λ, w0, M2)
    # Waist ray
    ξ = Ray(position(chief) + bob_the_builder * w0, direction(chief), λ)
    # Divergence ray
    div_dir = normalize(direction(chief) + bob_the_builder * tan(θ))
    η = Ray(position(chief), div_dir, λ)
    # Ensure that lambda is correct
    wavelength!(chief, λ)
    # Calculate E0 based on P0, assume zero initial phase offset
    I0 = 2 * P0 / (π * w0^2)
    E0 = electric_field(I0)
    return GaussianBeamlet(Beam(chief), Beam(ξ), Beam(η), λ, w0, E0)
end

point_on_beam(gauss::GaussianBeamlet, t::Real) = point_on_beam(gauss.chief, t)

"""
    paraprincipal_ray_parameters(c::AbstractRay, w::AbstractRay, d::AbstractRay, p0)

Calculates the heights and slopes of the waist and divergence ray as per the publication:

**DeJager, Donald, and Mark Noethen. "Gaussian beam parameters that use Coddington-based Y-Nu paraprincipal ray tracing." Applied Optics 31.13 (1992): 2199-2205.**

To calculate a signed slope in 3D, the vector between p0 and the orthogonal point on the div./waist ray is used as a reference, yielding an angle around 90°.
Larger than 90° is treated as a positive slope and vice versa. This method fails if the div. or waist ray height is zero, which occurs at beam waists!
In this, checks on the derived values are necessary. See also [`gauss_parameters`](@ref)

# Arguments
- `c`: [`AbstractRay`](@ref) representing the chief ray
- `w`: [`AbstractRay`](@ref) representing the waist ray
- `d`: [`AbstractRay`](@ref) representing the divergence ray
- `p0`: point on chief ray for which the ray heights and slopes are calculated

# Returns
- `h_w`: waist ray height
- `s_w`: waist ray slope
- `h_d`: div. ray height
_ `s_d`: div. ray slope
"""
function paraprincipal_ray_parameters(c::AbstractRay, w::AbstractRay, d::AbstractRay, p0)
    # waist ray height and slope at point p0 along chief ray
    l = line_plane_distance3d(p0, direction(c), position(w), direction(w))
    y = position(w) + l * direction(w) - p0
    h_w = norm(y)
    y /= h_w
    s_w = tan(π / 2 - angle3d(y, direction(w)))
    # divergence ray height and slope at point p0 along chief ray
    l = line_plane_distance3d(p0, direction(c), position(d), direction(d))
    y = position(d) + l * direction(d) - p0
    h_d = norm(y)
    y /= h_d
    s_d = tan(π / 2 - angle3d(y, direction(d)))
    return h_w, s_w, h_d, s_d
end

"""
    gauss_parameters(c::AbstractRay, w::AbstractRay, d::AbstractRay, p0)

Calculates the parameters of the stigmatic Gaussian beam as described by D. DeJager (1992).

# Arguments
- `c`: [`AbstractRay`](@ref) representing the chief ray
- `w`: [`AbstractRay`](@ref) representing the waist ray
- `d`: [`AbstractRay`](@ref) representing the divergence ray
- `p0`: point on chief ray for which the ray heights and slopes are calculated

# Returns
- `w`: local radius
- `R`: curvature, i.e. 1/r where r is the radius of curvature
- `ψ`: Gouy phase (note that -atan definition is used)
- `w0`: local beam waist radius
- `Hn`: optical invariant H divided by ref. index n
"""
function gauss_parameters(c::AbstractRay, w::AbstractRay, d::AbstractRay, p0)
    # heights and slopes, optical invariant H/n
    h_w, s_w, h_d, s_d = paraprincipal_ray_parameters(c, w, d, p0) 
    Hn = abs(h_w * s_d - h_d * s_w)
    # beam parameters
    E_kt = h_d * s_d + h_w * s_w
    F_kt = sqrt(s_d^2 + s_w^2)
    w = sqrt(h_d^2 + h_w^2)
    R = E_kt / w^2
    z = E_kt / F_kt^2
    ψ = -atan(1, √(1 / (R * z) - 1))
    w0 = Hn / F_kt
    # Catch NaNs at beam waist and correct Gouy phase sign based on curvature sign
    isnan(R) ? R = zero(R) : nothing
    isnan(ψ) ? ψ = zero(ψ) : nothing
    isnan(w0) ? w0 = w : nothing
    R < 0 ? ψ = -ψ : nothing
    return w, R, ψ, w0, Hn
end

"""
    gauss_parameters(gauss::GaussianBeamlet, z; hint=nothing)

Calculate the local waist radius and Gouy phase of an unastigmatic Gaussian beamlet at a specific distance `z` based on the method of J. Arnaud (1985) and D. DeJager (1992).

# Arguments

- `gauss`: the GaussianBeamlet for which parameters are to be calculated.
- `z`: the position along the beam at which to calculate the parameters.
- `hint`: an optional hint parameter `Union{Nothing, Tuple{Int, Vector{<:Real}}}` for the relevant point/index of the appropriate beam segment. If not provided, the function will automatically select the ray.

# Returns

- `w`: local radius
- `R`: curvature, i.e. 1/r where r is the radius of curvature
- `ψ`: Gouy phase (note that -atan definition is used)
- `w0`: local beam waist radius
"""
function gauss_parameters(gauss::GaussianBeamlet,
        z::Real;
        hint = nothing)
    # If hint not provided, find values automatically
    if isnothing(hint)
        p0, index = point_on_beam(gauss, z)
    else
        p0, index = hint
    end
    # Select relevant beam section
    chief = gauss.chief.rays[index]
    div = gauss.divergence.rays[index]
    waist = gauss.waist.rays[index]
    # Calculate beam parameters
    λ = wavelength(gauss)
    w, R, ψ, w0, Hn = gauss_parameters(chief, waist, div, p0)
    # Test optical invariant
    if !isapprox(Hn, λ / π, atol = 1e-6)
        println("H/n not fulfilled at z=$z")
    end
    return w, R, ψ, w0
end

"""
    gauss_parameters(gauss::GaussianBeamlet, zs::AbstractArray)

Return the parameters of the `GaussianBeamlet` along the specified positions in `zs`.
"""
function gauss_parameters(gauss::GaussianBeamlet{G}, zs::AbstractArray) where {G}
    w = Vector{G}(undef, length(zs))
    R = Vector{G}(undef, length(zs))
    ψ = Vector{G}(undef, length(zs))
    w0 = Vector{G}(undef, length(zs))
    for (i, z) in enumerate(zs)
        w[i], R[i], ψ[i], w0[i] = gauss_parameters(gauss, z)
    end
    return w, R, ψ, w0
end

"""
    electric_field(gauss::GaussianBeamlet, r, z)

Calculates the electric field phasor [V/m] of `gauss` at the radial and longitudinal positions `r` and `z`.
"""
function electric_field(gauss::GaussianBeamlet, r, z)
    point, index = point_on_beam(gauss, z)
    w, R, ψ, w0 = gauss_parameters(gauss, z, hint = (point, index))
    k = wave_number(wavelength(gauss))
    # Calculate new local field strength based on E0*w0 = const.
    E0 = beam_amplitude(gauss) * (beam_waist(gauss) / w0)
    return electric_field(r, z, E0, w0, w, k, ψ, R)
end

"""
    isparaxial(system, gb::GaussianBeamlet, threshold=π/4)

Tests the angle between the waist and divergence beams and refractive surfaces.
A target threshold of π/4 or 45° is assumed before abberations become dominant.
"""
isparaxial(system::AbstractSystem, gb::GaussianBeamlet, threshold::Real = π / 4) = isparaxial(system,
    gb.waist,
    threshold) & isparaxial(system,
    gb.divergence,
    threshold)

"""
    istilted(system::System, gb::GaussianBeamlet)

Tests if refractive elements are tilted with respect to the beamlet optical axis, i.e. introduce simple astigmatism.
"""
istilted(system::AbstractSystem, gb::GaussianBeamlet) = !isparaxial(system, gb.chief, 0)

function isparentbeam(beam::GaussianBeamlet, ray_id)
    c = isparentbeam(beam.chief, ray_id)
    w = isparentbeam(beam.waist, ray_id)
    d = isparentbeam(beam.divergence, ray_id)
    return any((c, w, d))
end
