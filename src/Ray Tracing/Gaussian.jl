"""
    GaussianBeamlet

Ray representation of the unastigmatic Gaussian beam as per J. Arnaud (1985).
The beam quality `M2` is fully considered via the divergence angle.
"""
mutable struct GaussianBeamlet{T} <: AbstractEntity
    id::UUID
    chief::Beam{T}
    waist::Beam{T}
    divergence::Beam{T}
    w0::T
    M2::T
    zR::T
end

beam_waist(beam::GaussianBeamlet) = beam.w0
beam_quality(beam::GaussianBeamlet) = beam.M2
rayleigh_range(beam::GaussianBeamlet) = beam.zR

"""
    GaussianBeamlet(chief::Ray{T}, λ=1064, w0=1; support=[0,0,1], M2=1)

Construct a Gaussian beamlet at its waist with a specified beam diameter.

# Arguments
- `chief`: Chief ray defining the origin of the Gaussian beamlet.
- `λ`: Wavelength of the beamlet in nanometers. Default value is 1064.
- `w0`: Beam waist (radius) in millimeters. Default value is 1.

# Keyword Arguments
- `support`: Support vector that can be adjusted for beamlet construction.
"""
function GaussianBeamlet(chief::Ray{T}, λ=1064, w0=1; support=[0,0,1], M2=1) where T
    # Create orthogonal vector for construction purposes
    bob_the_builder = orthogonal3d(direction(chief), support)
    # Divergence angle in rad
    λ *= 1e-9
    w0 *= 1e-3
    θ = divergence_angle(λ, w0, M2)
    zR = rayleigh_range(λ, w0, M2)
    # Waist ray
    ξ = Ray(position(chief) + bob_the_builder * w0, direction(chief), λ)
    # Divergence ray
    div_dir = direction(chief) + bob_the_builder * tan(θ)
    normalize3d!(div_dir)
    η = Ray(position(chief), div_dir, λ)
    # Ensure that lambda is correct
    chief.parameters.λ = λ 
    return GaussianBeamlet{T}(uuid4(), Beam(chief), Beam(ξ), Beam(η), w0, M2, zR)
end

point_on_beam(gauss::GaussianBeamlet, t::Real) = point_on_beam(gauss.chief, t)

"""
    gauss_parameters(gauss::GaussianBeamlet, z; hint::Union{Nothing, Tuple{Int, Vector{<:Real}}}=nothing)

Calculate the local waist radius and Gouy phase of a Gaussian beamlet at a specific distance `z` based on the method of J. Arnaud (1985) and D. DeJager (1992).

# Arguments
- `gauss`: the GaussianBeamlet object for which parameters are to be calculated.
- `z`: the position along the beam at which to calculate the parameters.
- `hint`: an optional hint parameter for the index of the ray to consider. If not provided, the function will automatically select the appropriate ray.

# Returns
- `w`: waist radius
- `R`: radius of curvature
- `ψ`: Gouy phase (note that -atan definition is used)
"""
function gauss_parameters(gauss::GaussianBeamlet, z::Real; hint::Nullable{Tuple{Vector{<:Real},Int}}=nothing)
    if isnothing(hint)
        p0, index = point_on_beam(gauss, z)
    else
        p0, index = hint
    end
    chief = gauss.chief.rays[index]
    #=
    Divergence ray height and slope
    - find divergence ray "height" and "slope" at intersection point y0 with target plane at p0 of chief ray
    - ray height "y_d" is length between p0 and y0
    - ray slope "m_d" is angle between vector p0 -> y0 and divergence ray direction -> gives unambiguous angle for signed ray slope calculation
    - fails if y_d is zero -> catch R=Inf and ψ=0
    =#
    div = gauss.divergence.rays[index]
    intersection = intersect3d(p0, direction(chief), div)
    y0 = position(div) + length(intersection) * direction(div)
    dy = y0 - p0
    y_d = norm3d(dy)
    normalize3d!(dy)
    m_d = tan(π/2 - angle3d(dy, direction(div)))
    # Waist ray height and slope
    waist = gauss.waist.rays[index]
    intersection = intersect3d(p0, direction(chief), waist)
    y0 = position(waist) + length(intersection) * direction(waist)
    dy = y0 - p0
    y_w = norm3d(dy)
    normalize3d!(dy)
    m_w = tan(π/2 - angle3d(dy, direction(waist)))
    # Beam parameters as per Arnaud (1985) and DeJager (1992)
    temp = y_d * m_d + y_w * m_w # E_kt
    w = sqrt(y_d^2 + y_w^2)
    R = w^2 / (temp)
    z = temp / (m_d^2 + m_w^2)
    ψ = -atan(1, √(R/z - 1))
    # Correct Gouy phase sign based on curvature sign and catch Inf/NaN values for R and ψ
    if R < 0
        ψ = -ψ
    end
    if isinf(R)
        R = NaN
    end
    if isnan(ψ)
        ψ = 0
    end
    return w, R, ψ
end

"""
    gauss_parameters(gauss::GaussianBeamlet, zs::AbstractArray)

Return the parameters of the `GaussianBeamlet` along the specified positions in `zs`.
"""
function gauss_parameters(gauss::GaussianBeamlet{G}, zs::AbstractArray) where G
    w = Vector{G}(undef, length(zs))
    R = Vector{G}(undef, length(zs))
    ψ = Vector{G}(undef, length(zs))
    for (i, z) in enumerate(zs)
        w[i], R[i], ψ[i] = gauss_parameters(gauss, z)
    end
    return w, R, ψ
end

function electric_field(beam::GaussianBeamlet, r, z)
    point, index = point_on_beam(beam, z)
    ray_parameters = parameters(beam.chief.rays[index])
    E0 = ray_parameters.I
    w0 = beam_waist(beam)
    w, R, ψ = gauss_parameters(beam, z, hint=(point, index))
    k = wave_number(ray_parameters.λ)
    return electric_field(r, z, E0, w0, w, k, ψ, R)
end

function isparentbeam(beam::GaussianBeamlet, ray_id)
    c = isparentbeam(beam.chief, ray_id)
    w = isparentbeam(beam.waist, ray_id)
    d = isparentbeam(beam.divergence, ray_id)
    return any([c,w,d])
end
