"""
    reflection3d(dir, normal)

Calculates the reflection between an input vector `dir` and surface `normal` vector in R³.
Vectors `dir` and `normal` must have **unit length**!
"""
function reflection3d(dir, normal)
    return dir - 2 * dot(dir, normal) * normal
end

"""
    refraction3d(dir, normal, n1, n2)

Calculates the refraction between an input vector `dir` and surface `normal` vector in R³.
`n1` is the "outside" refractive index and `n2` is the "inside" refractive index.
The function returns the new direction of propagation and a boolean flag to indicate if internal refraction has occured.

Vectors `dir` and `normal` must have **unit length**!

# Total internal reflection

If the critical angle for n1, n2 and the incident angle is reached, the ray is reflected internally instead!

# Arguments

- `dir`: direction vector of incoming ray
- `normal`: surface normal at point of intersection
- `n1`: index of ref. before refraction
- `n2`: index of ref. after refraction
"""
function refraction3d(dir::AbstractArray, normal::AbstractArray, n1::Real, n2::Real)
    # dir and normal must have unit length!
    isapprox(norm(dir), 1) || throw(ArgumentError("dir must have  unit length"))
    isapprox(norm(normal), 1) || throw(ArgumentError(lazy"norm must have  unit length: $(norm(normal))"))
    n = n1 / n2
    cosθi = -dot(normal, dir)
    sinθt² = n^2 * (1 - cosθi^2)
    # Check for total reflection
    if sinθt² > 1.0
        return (reflection3d(dir, normal), true)
    end
    cosθt = sqrt(1 - sinθt²)
    return (@. n * dir + (n * cosθi - cosθt) * normal, false)
end

"""
    lensmakers_eq(R1, R2, n)

Calculates the thin lens focal length based on the radius of curvature `R1`/`R2` and the lens refractive index `n`.
If center of sphere is on left then R < 0. If center of sphere is on right then R > 0.
"""
lensmakers_eq(R1, R2, n) = 1 / ((n - 1) * (1 / R1 - 1 / R2))



rayleigh_range(λ, w0, M2, n=1) = π * n * w0^2 / λ / M2

beam_waist(z, w0, zr) = w0 * sqrt(1 + (z / zr)^2)

gouy_phase(z, zr) = -atan(z / zr) # some definitions with +

wavefront_curvature(z, zr) = z / (z^2 + zr^2) # 1/r

divergence_angle(λ, w0, M2) = M2 * λ / (π * w0)

wave_number(λ) = 2π / λ

"""
    electric_field(r, z, E0, w0, w, k, ψ, R) -> ComplexF64

Computes the analytical complex electric field distribution of a stigmatic TEM₀₀ Gaussian beam which is described by:

```math
E(r,z) = {E_0}\\frac{{{w_0}}}{{w(z)}}\\exp\\left( { - \\frac{{{r^2}}}{{w{{(z)}^2}}}} \\right)\\exp\\left(i\\left[ {kz + \\psi + \\frac{{k{r^2}}}{2 R(z)}} \\right] \\right)
```

# Arguments

- `r`: radial distance from beam origin
- `z`: axial distance from beam origin
- `E0`: peak electric field amplitude
- `w0`: waist radius
- `w`: local beam radius
- `k`: wave number, equal to `2π/λ`
- `ψ`: Gouy phase shift (defined as ``-\\text{atan}\\left(\\frac{z}{z_r}\\right)`` !)
- `R`: wavefront curvature, i.e. 1/r (radius of curvature)
"""
electric_field(r::Real, z::Real, E0, w0, w, k, ψ, R) = E0 * w0 / w * exp(-r^2 / w^2) * exp(im * (k * z + ψ + (k * r^2 * R) / 2))

function electric_field(r::Real, z::Real, E0, w0, λ, M2=1)
    zr = rayleigh_range(λ, w0, M2)
    w = beam_waist(z, w0, zr)
    k = wave_number(λ)
    ψ = gouy_phase(z, zr)
    R = wavefront_curvature(z, zr)
    return electric_field(r, z, E0, w0, w, k, ψ, R)
end

"""
    electric_field(I::Real, Z=Z_vacuum, ϕ=0)

Calculates the E-field phasor in [V/m] for a given intensity `I` and phase ϕ. Vacuum wave impedance is assumed.
"""
electric_field(I::Real, Z=Z_vacuum, ϕ=0) = sqrt(2 * I * Z) * exp(im * ϕ)

"Calculates the intensity in [W/m²] for a given complex electric field phasor `E`. Vacuum wave impedance is assumed."
intensity(E::Number, Z=Z_vacuum) = abs2(E) / (2 * Z)

"""
fresnel_coefficients(θ, n)

Calculates the complex Fresnel coefficients for reflection and transmission based on the incident angle `θ` in [rad]
and the refractive index ratio `n = n₂ / n₁`. Returns rₛ, rₚ, tₛ and tₚ.

# Signs

!!! info
    The signs of rₛ, rₚ are based on the definition by Fowles (1975, 2nd Ed. p. 44) and Peatross (2015, 2023 Ed. p. 78)
"""
function fresnel_coefficients(θ::T, n::Number) where T
    C = Complex{T}
    cost = cos(θ)
    n2s2 = sqrt(C(n^2 - sin(θ)^2))
    # Calculate coefficients for reflection/transmission
    rs = (cost - n2s2) / (cost + n2s2)
    rp = (-n^2 * cost + n2s2) / (n^2 * cost + n2s2)
    ts = rs + 1
    tp = 2*n*cost / (n^2 * cost + n2s2)
    return rs, rp, ts, tp
end

function fresnel_coefficients(θ::AbstractArray{T}, n::Number) where T
    rs = Vector{Complex{T}}(undef, length(θ))
    rp = Vector{Complex{T}}(undef, length(θ))
    ts = Vector{Complex{T}}(undef, length(θ))
    tp = Vector{Complex{T}}(undef, length(θ))
    for (i, theta) in enumerate(θ)
        rs[i], rp[i], ts[i], tp[i] = fresnel_coefficients(theta, n)
    end
    return rs, rp, ts, tp
end

is_internally_reflected(rp::Number, rs::Number) = isapprox(abs2(rs), 1, atol=1e-6) && isapprox(abs2(rp), 1, atol=1e-6)

"""
    sag(r::Real, l::Real)

Calculates the sag of a cut circle with radius `r` and chord length `l`
"""
sag(r::Real, l::Real) = r - sqrt(r^2 - 0.25 * l^2)

function check_sag(r, d)
    if abs(2r) < d
        throw(ArgumentError("Radius of curvature (r = $(r)) must be ≥ than half the diameter (d = $(d)) or an illegal shape results!"))
    else
        return nothing
    end
end

"""
    numerical_aperture(θ, n=1)

Returns the `NA` for a opening half-angle `θ` and scalar ref. index `n`.
For more information refer to [this website](https://www.rp-photonics.com/numerical_aperture.html).
"""
numerical_aperture(θ::Real, n::Real=1.0) = n*sin(θ)

"""
    visibility(opt_pwr)

Calculates e.g. the interferometric contrast from a series of optical power measurements.
For more information go [here](https://en.wikipedia.org/wiki/Interferometric_visibility).
"""
function visibility(opt_pwr)
    I_max = maximum(opt_pwr)
    I_min = minimum(opt_pwr)
    return (I_max - I_min) / (I_max + I_min)
end 