"""
    AbstractMedium

Abstract supertype for optical propagation media governing bulk properties (refractive index, dispersion, attenuation, gain, thermo-optic and piezo-optic responses).

# Core Interface

Subtypes of `AbstractMedium` should implement:
- `refractive_index(medium, λ; T=293.15, P=101325.0, kwargs...)`: Real part of the refractive index ``n(\\lambda, T, P)`` for ray tracing and Snell refraction.
- `complex_refractive_index(medium, λ; T=293.15, P=101325.0, kwargs...)`: Full complex refractive index ``\\tilde{n} = n + i\\kappa`` for thin-film, absorption, and gain calculations.
- `extinction_coefficient(medium, λ; kwargs...)`: Imaginary part ``\\kappa`` representing absorption (``\\kappa > 0``) or gain (``\\kappa < 0``).
- `absorption_coefficient(medium, λ; kwargs...)`: Linear attenuation/gain coefficient ``\\alpha = \\frac{4\\pi \\kappa}{\\lambda}`` [1/m].
"""
abstract type AbstractMedium end

"""
    Ambient <: AbstractMedium

Represents ambient space (e.g. Air or Vacuum), dynamically resolving to the enclosing system's ambient medium.
"""
struct Ambient <: AbstractMedium end

refractive_index(::Ambient, λ::Real = 0.0; kwargs...) = 1.0
complex_refractive_index(::Ambient, λ::Real = 0.0; kwargs...) = ComplexF64(1.0, 0.0)
extinction_coefficient(::Ambient, λ::Real = 0.0; kwargs...) = 0.0
absorption_coefficient(::Ambient, λ::Real = 0.0; kwargs...) = 0.0

"""
    IsotropicMedium{D, A} <: AbstractMedium

An isotropic optical medium characterized by a dispersion formula (real or complex) and bulk attenuation/gain.

# Fields
- `name::Symbol`: Material name (e.g. `:N_BK7`, `:FusedSilica`, `:Water`, `:Gold`, `:Constant`)
- `dispersion::D`: Dispersion function `(λ; T, P, kwargs...) -> n` or `(λ) -> n` (Real or Complex), or a constant Number (`Real` or `Complex`)
- `attenuation::A`: Bulk attenuation function `(λ; T, P, kwargs...) -> α` \\[1/m\\] (or `nothing`)
"""
struct IsotropicMedium{D, A} <: AbstractMedium
    name::Symbol
    dispersion::D
    attenuation::A
end

IsotropicMedium(name::Symbol, dispersion) = IsotropicMedium(name, dispersion, nothing)
IsotropicMedium(dispersion::Number) = IsotropicMedium(dispersion isa Complex ? :ComplexConstant : :Constant, dispersion, nothing)
IsotropicMedium(dispersion) = IsotropicMedium(:Dispersive, dispersion, nothing)

# Evaluates dispersion callable or constant with zero overhead and full inlining
@inline _eval_dispersion(disp::Number, λ::Real) = disp
@inline _eval_dispersion(disp::Number, λ::Real, kwargs) = disp
@inline _eval_dispersion(disp, λ::Real) = disp(λ)
@inline function _eval_dispersion(disp, λ::Real, kwargs)
    if isempty(kwargs)
        return disp(λ)
    else
        return disp(λ; kwargs...)
    end
end

"""
    complex_refractive_index(m::IsotropicMedium, λ=0.0; kwargs...) -> ComplexF64

Returns the full complex refractive index ``\\tilde{n} = n + i\\kappa``.
"""
@inline function complex_refractive_index(m::IsotropicMedium, λ::Real = 0.0; kwargs...)::ComplexF64
    val = _eval_dispersion(m.dispersion, λ, kwargs)
    if val isa Complex
        return ComplexF64(val)
    else
        n_real = Float64(val)
        if m.attenuation !== nothing && λ > 0
            α = _eval_dispersion(m.attenuation, λ, kwargs)
            κ = Float64(α) * λ / (4 * π)
            return ComplexF64(n_real, κ)
        else
            return ComplexF64(n_real, 0.0)
        end
    end
end

"""
    refractive_index(m::IsotropicMedium, λ=0.0; kwargs...) -> Float64

Returns the real part of the refractive index ``n = \\text{Re}(\\tilde{n})``.
"""
@inline function refractive_index(m::IsotropicMedium, λ::Real = 0.0; kwargs...)::Float64
    val = _eval_dispersion(m.dispersion, λ, kwargs)
    return Float64(real(val))
end

"""
    extinction_coefficient(m::AbstractMedium, λ=0.0; kwargs...) -> Float64

Returns the extinction coefficient ``\\kappa = \\text{Im}(\\tilde{n})`` (positive for absorption, negative for gain).
"""
@inline function extinction_coefficient(m::AbstractMedium, λ::Real = 0.0; kwargs...)::Float64
    return Float64(imag(complex_refractive_index(m, λ; kwargs...)))
end

"""
    absorption_coefficient(m::AbstractMedium, λ=0.0; kwargs...) -> Float64

Returns the linear attenuation coefficient ``\\alpha = \\frac{4\\pi \\kappa}{\\lambda}`` [1/m].
"""
function absorption_coefficient(m::AbstractMedium, λ::Real = 0.0; kwargs...)::Float64
    if λ <= 0
        return 0.0
    end
    κ = extinction_coefficient(m, λ; kwargs...)
    return (4 * π * κ) / λ
end

"""
    AbstractAnisotropicMedium <: AbstractMedium

Abstract supertype for anisotropic and birefringent optical media where refractive index and dielectric response depend on direction and polarization.
"""
abstract type AbstractAnisotropicMedium <: AbstractMedium end

"""
    UniaxialMedium{No, Ne, C, A} <: AbstractAnisotropicMedium

A uniaxial birefringent crystal characterized by ordinary index ``n_o``, extraordinary index ``n_e``,
and optic axis vector ``\\hat{\\mathbf{c}}``.

# Fields
- `name::Symbol`: Crystal material name (e.g. `:Calcite`, `:BBO`, `:Quartz`, `:LiNbO3`)
- `dispersion_o::No`: Ordinary dispersion formula or constant ``n_o(\\lambda)``
- `dispersion_e::Ne`: Extraordinary dispersion formula or constant ``n_e(\\lambda)``
- `optic_axis::C`: Optic axis unit vector ``\\hat{\\mathbf{c}} \\in \\mathbb{R}^3``
- `attenuation::A`: Bulk attenuation formula (or `nothing`)
"""
struct UniaxialMedium{No, Ne, C <: Point3, A} <: AbstractAnisotropicMedium
    name::Symbol
    dispersion_o::No
    dispersion_e::Ne
    optic_axis::C
    attenuation::A
end

function UniaxialMedium(
    name::Symbol,
    no,
    ne,
    optic_axis::AbstractArray = Point3(0.0, 0.0, 1.0);
    attenuation = nothing
)
    c_axis = normalize(Point3{Float64}(optic_axis))
    return UniaxialMedium(name, no, ne, c_axis, attenuation)
end

function UniaxialMedium(
    no::Union{Number, Function, DiscreteRefractiveIndex, SellmeierEquation},
    ne,
    optic_axis::AbstractArray = Point3(0.0, 0.0, 1.0);
    attenuation = nothing
)
    return UniaxialMedium(:Uniaxial, no, ne, optic_axis; attenuation)
end

"""
    BiaxialMedium{Nx, Ny, Nz, A} <: AbstractAnisotropicMedium

A biaxial anisotropic optical crystal characterized by three principal refractive indices ``n_x, n_y, n_z``.
"""
struct BiaxialMedium{Nx, Ny, Nz, A} <: AbstractAnisotropicMedium
    name::Symbol
    dispersion_x::Nx
    dispersion_y::Ny
    dispersion_z::Nz
    attenuation::A
end

function BiaxialMedium(
    name::Symbol,
    nx,
    ny,
    nz;
    attenuation = nothing
)
    return BiaxialMedium(name, nx, ny, nz, attenuation)
end

function BiaxialMedium(
    nx,
    ny,
    nz;
    attenuation = nothing
)
    return BiaxialMedium(:Biaxial, nx, ny, nz; attenuation)
end

optic_axis(m::UniaxialMedium) = m.optic_axis
is_uniaxial(::UniaxialMedium) = true
is_uniaxial(::AbstractMedium) = false
is_biaxial(::BiaxialMedium) = true
is_biaxial(::AbstractMedium) = false

"""
    refractive_index_o(m::UniaxialMedium, λ=0.0; kwargs...) -> Float64

Returns the ordinary refractive index ``n_o(\\lambda)``.
"""
@inline function refractive_index_o(m::UniaxialMedium, λ::Real = 0.0; kwargs...)::Float64
    val = _eval_dispersion(m.dispersion_o, λ, kwargs)
    return Float64(real(val))
end

"""
    refractive_index_e(m::UniaxialMedium, λ=0.0; kwargs...) -> Float64

Returns the principal extraordinary refractive index ``n_e(\\lambda)``.
"""
@inline function refractive_index_e(m::UniaxialMedium, λ::Real = 0.0; kwargs...)::Float64
    val = _eval_dispersion(m.dispersion_e, λ, kwargs)
    return Float64(real(val))
end

"""
    birefringence(m::UniaxialMedium, λ=0.0; kwargs...) -> Float64

Returns the crystal birefringence ``\\Delta n = n_e(\\lambda) - n_o(\\lambda)``.
Positive for positive uniaxial crystals (e.g. Quartz), negative for negative uniaxial crystals (e.g. Calcite).
"""
@inline function birefringence(m::UniaxialMedium, λ::Real = 0.0; kwargs...)::Float64
    return refractive_index_e(m, λ; kwargs...) - refractive_index_o(m, λ; kwargs...)
end

@inline function refractive_index(m::UniaxialMedium, λ::Real = 0.0; kwargs...)::Float64
    return refractive_index_o(m, λ; kwargs...)
end

"""
    refractive_index(m::UniaxialMedium, λ::Real, k_dir::AbstractArray; kwargs...) -> Float64

Calculates the angle-dependent extraordinary refractive index ``n_e(\\theta)`` for propagation wavevector direction `k_dir`:
```math
\\frac{1}{n_e(\\theta)^2} = \\frac{\\cos^2\\theta}{n_o^2} + \\frac{\\sin^2\\theta}{n_e^2}
```
where ``\\cos\\theta = |\\hat{\\mathbf{k}} \\cdot \\hat{\\mathbf{c}}|``.
"""
function refractive_index(m::UniaxialMedium, λ::Real, k_dir::AbstractArray; kwargs...)::Float64
    no = refractive_index_o(m, λ; kwargs...)
    ne = refractive_index_e(m, λ; kwargs...)
    k_unit = normalize(Point3{Float64}(k_dir))
    c_unit = optic_axis(m)
    cosθ = abs(dot(k_unit, c_unit))
    sinθ2 = max(0.0, 1.0 - cosθ^2)
    denom = ne^2 * cosθ^2 + no^2 * sinθ2
    if denom <= 0
        return no
    end
    return (no * ne) / sqrt(denom)
end

"""
    dielectric_tensor(medium, λ=0.0; kwargs...) -> SMatrix{3, 3, Float64, 9}

Returns the 3×3 relative permittivity tensor ``\\boldsymbol{\\varepsilon}_r(\\lambda)``.
"""
function dielectric_tensor(m::Ambient, λ::Real = 0.0; kwargs...)
    return SMatrix{3, 3, Float64, 9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
end

function dielectric_tensor(m::IsotropicMedium, λ::Real = 0.0; kwargs...)
    n2 = refractive_index(m, λ; kwargs...)^2
    return SMatrix{3, 3, Float64, 9}(n2, 0.0, 0.0, 0.0, n2, 0.0, 0.0, 0.0, n2)
end

function dielectric_tensor(m::UniaxialMedium, λ::Real = 0.0; kwargs...)
    no2 = refractive_index_o(m, λ; kwargs...)^2
    ne2 = refractive_index_e(m, λ; kwargs...)^2
    c = optic_axis(m)
    # ε_r = ε_o * I + (ε_e - ε_o) * (c ⊗ c)
    Δε = ne2 - no2
    return SMatrix{3, 3, Float64, 9}(
        no2 + Δε * c[1] * c[1], Δε * c[2] * c[1],       Δε * c[3] * c[1],
        Δε * c[1] * c[2],       no2 + Δε * c[2] * c[2], Δε * c[3] * c[2],
        Δε * c[1] * c[3],       Δε * c[2] * c[3],       no2 + Δε * c[3] * c[3]
    )
end

function dielectric_tensor(m::BiaxialMedium, λ::Real = 0.0; kwargs...)
    nx2 = Float64(real(_eval_dispersion(m.dispersion_x, λ, kwargs)))^2
    ny2 = Float64(real(_eval_dispersion(m.dispersion_y, λ, kwargs)))^2
    nz2 = Float64(real(_eval_dispersion(m.dispersion_z, λ, kwargs)))^2
    return SMatrix{3, 3, Float64, 9}(nx2, 0.0, 0.0, 0.0, ny2, 0.0, 0.0, 0.0, nz2)
end

"""
    medium_from(x) -> AbstractMedium

Converts a medium, constant refractive index, dispersion function, or optical object into an `AbstractMedium`.
"""
medium_from(m::AbstractMedium) = m
medium_from(n::Number) = IsotropicMedium(n)
medium_from(::AbstractObject) = Ambient()
medium_from(f) = IsotropicMedium(f)

