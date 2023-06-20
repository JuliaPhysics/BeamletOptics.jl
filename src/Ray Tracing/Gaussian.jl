"""
    GaussianBeamlet

Ray representation of the **unastigmatic** Gaussian beam as per J. Arnaud (1985).
The beam quality `M2` is fully considered via the divergence angle.
"""
mutable struct GaussianBeamlet{T} <: AbstractEntity
    id::UUID
    chief::Beam{T}
    waist::Beam{T}
    divergence::Beam{T}
    λ::T
    w0::T
    M2::T
    E0::Complex{T}
end

wavelength(beam::GaussianBeamlet) = beam.λ
beam_waist(beam::GaussianBeamlet) = beam.w0
beam_amplitude(beam::GaussianBeamlet) = beam.E0
beam_quality(beam::GaussianBeamlet) = beam.M2

"""
    GaussianBeamlet(chief::Ray{T}, λ=1064, w0=1; support=[0,0,1], M2=1)

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
function GaussianBeamlet(chief::Ray{T}, λ=1000e-9, w0=1e-3; M2=1, P0=1e-3, support=[0, 0, 1]) where {T}
    # Create orthogonal vector for construction purposes
    bob_the_builder = normal3d(direction(chief), support)
    # Divergence angle in rad
    θ = divergence_angle(λ, w0, M2)
    # Waist ray
    ξ = Ray(position(chief) + bob_the_builder * w0, direction(chief), λ)
    # Divergence ray
    div_dir = direction(chief) + bob_the_builder * tan(θ)
    normalize3d!(div_dir)
    η = Ray(position(chief), div_dir, λ)
    # Ensure that lambda is correct
    chief.parameters.λ = λ
    # Calculate E0 based on P0, assume zero initial phase offset
    I0 = 2*P0/(π*w0^2)
    E0 = electric_field(I0)
    return GaussianBeamlet{T}(uuid4(), Beam(chief), Beam(ξ), Beam(η), λ, w0, M2, E0)
end

point_on_beam(gauss::GaussianBeamlet, t::Real) = point_on_beam(gauss.chief, t)
point_on_beam!(point::AbstractVector, gauss::GaussianBeamlet, t::Real) = point_on_beam!(point, gauss.chief, t)

"""
    gauss_parameters(gauss::GaussianBeamlet, z; hint::Union{Nothing, Tuple{Int, Vector{<:Real}}}=nothing)

Calculate the local waist radius and Gouy phase of an unastigmatic Gaussian beamlet at a specific distance `z` based on the method of J. Arnaud (1985) and D. DeJager (1992).

# Arguments
- `gauss`: the GaussianBeamlet object for which parameters are to be calculated.
- `z`: the position along the beam at which to calculate the parameters.
- `hint`: an optional hint parameter for the relevant point/index of the appropriate beam segmnent. If not provided, the function will automatically select the ray.

# Returns
- `w`: local radius
- `R`: curvature, i.e. 1/r where r is the radius of curvature
- `ψ`: Gouy phase (note that -atan definition is used)
- `w0`: local beam waist radius
"""
function gauss_parameters(gauss::GaussianBeamlet, z::Real, y0::AbstractVector{Float64}=zeros(3); hint::Nullable{Tuple{Vector{<:Real},Int}}=nothing)
    if isnothing(hint)
        p0, index = point_on_beam(gauss, z)
    else
        p0, index = hint
    end
    chief = gauss.chief.rays[index]
    n = refractive_index(chief)
    λ = wavelength(gauss)
    #=
    Divergence ray height and slope (same for waist ray)
    - find divergence ray "height" and "slope" at intersection point y0 with target plane at p0 of chief ray
    - ray height "y_d" is length between p0 and y0
    - ray slope "m_d" is angle between vector p0 -> y0 and divergence ray direction -> gives unambiguous angle for signed ray slope calculation
    - fails if y_d is zero -> catch R=Inf, ψ=0 and H=λ/π
    =#
    div = gauss.divergence.rays[index]
    intersection_len = line_plane_distance3d(p0, direction(chief), position(div), direction(div))
    @. y0 = $position(div) + intersection_len * $direction(div) - p0    
    y_d = norm3d(y0)
    y0 ./= y_d
    m_d = tan(π / 2 - angle3d(y0, direction(div)))
    # Waist ray height and slope
    waist = gauss.waist.rays[index]
    intersection_len = line_plane_distance3d(p0, direction(chief), position(waist), direction(waist))
    @. y0 = $position(waist) + intersection_len * $direction(waist) - p0    
    y_w = norm3d(y0)
    y0 ./= y_w
    m_w = tan(π / 2 - angle3d(y0, direction(waist)))
    # Beam parameters as per Arnaud (1985) and DeJager (1992)
    H = abs(n * (y_w * m_d - y_d * m_w))
    # Test optical invariant
    if !isapprox(H, λ / π, atol=1e-6)
        H = λ / π
        # println("H not fulfilled at z=$z")
    end
    E_kt = y_d * m_d + y_w * m_w
    F_kt = sqrt(m_d^2 + m_w^2)
    w = sqrt(y_d^2 + y_w^2)
    R = E_kt / w^2
    z = E_kt / F_kt^2
    ψ = -atan(1, √(1/(R*z) - 1))
    w0 = H / (n * F_kt)
    # Catch NaNs and correct Gouy phase sign based on curvature sign
    isnan(R) ? R = zero(R) : nothing
    isnan(ψ) ? ψ = zero(ψ) : nothing
    isnan(w0) ? w0 = w : nothing
    R < 0 ? ψ = -ψ : nothing
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
    electric_field(gauss::GaussianBeamlet, r, z, g_b=zeros(3), p_b=zeros(3))

Calculates the electric field phasor [V/m] of `gauss` at the radial and longitudinal positions `r` and `z`.
Optionally, buffer vectors `g_b` and `p_b` can be passed.  
"""
function electric_field(gauss::GaussianBeamlet, r, z, g_b=zeros(3), p_b=zeros(3))
    point, index = point_on_beam!(p_b, gauss, z)
    w, R, ψ, w0 = gauss_parameters(gauss, z, g_b, hint=(point, index))
    k = wave_number(wavelength(gauss))
    # Calculate new local field strength based on E0*w0 = const.
    E0 = beam_amplitude(gauss) * (beam_waist(gauss)/w0)
    return electric_field(r, z, E0, w0, w, k, ψ, R)
end

"""
    isparaxial(system, gb::GaussianBeamlet, threshold=π/4)

Tests the angle between the waist and divergence beams and refractive surfaces. 
A target threshold of π/4 or 45° is assumed before abberations become dominant.
"""
isparaxial(system::AbstractSystem, gb::GaussianBeamlet, threshold::Real=π/4) = isparaxial(system, gb.waist, threshold) & isparaxial(system, gb.divergence, threshold)

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
