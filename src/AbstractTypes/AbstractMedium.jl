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
- `attenuation::A`: Bulk attenuation function `(λ; T, P, kwargs...) -> α` [1/m] (or `nothing`)
"""
struct IsotropicMedium{D, A} <: AbstractMedium
    name::Symbol
    dispersion::D
    attenuation::A
end

IsotropicMedium(name::Symbol, dispersion) = IsotropicMedium(name, dispersion, nothing)
IsotropicMedium(dispersion::Number) = IsotropicMedium(dispersion isa Complex ? :ComplexConstant : :Constant, dispersion, nothing)
IsotropicMedium(dispersion) = IsotropicMedium(:Dispersive, dispersion, nothing)

# Evaluates dispersion callable or constant, handling keyword args gracefully
@inline function _eval_dispersion(disp::Number, λ::Real; kwargs...)
    return disp
end

@inline function _eval_dispersion(disp::Function, λ::Real; T::Real = 293.15, P::Real = 101325.0, kwargs...)
    # Try calling with environmental kwargs, fallback to 1-arg if not accepted
    try
        return disp(λ; T = T, P = P, kwargs...)
    catch e
        if e isa MethodError
            return disp(λ)
        else
            rethrow(e)
        end
    end
end

@inline function _eval_dispersion(disp, λ::Real; kwargs...)
    # Generic callable functor (SellmeierEquation, DiscreteRefractiveIndex, etc.)
    return disp(λ)
end

"""
    complex_refractive_index(m::IsotropicMedium, λ=0.0; T=293.15, P=101325.0, kwargs...) -> ComplexF64

Returns the full complex refractive index ``\\tilde{n} = n + i\\kappa``.
"""
function complex_refractive_index(m::IsotropicMedium, λ::Real = 0.0; T::Real = 293.15, P::Real = 101325.0, kwargs...)::ComplexF64
    val = _eval_dispersion(m.dispersion, λ; T = T, P = P, kwargs...)
    if val isa Complex
        return ComplexF64(val)
    else
        n_real = Float64(val)
        if m.attenuation !== nothing && λ > 0
            α = _eval_dispersion(m.attenuation, λ; T = T, P = P, kwargs...)
            κ = Float64(α) * λ / (4 * π)
            return ComplexF64(n_real, κ)
        else
            return ComplexF64(n_real, 0.0)
        end
    end
end

"""
    refractive_index(m::IsotropicMedium, λ=0.0; T=293.15, P=101325.0, kwargs...) -> Float64

Returns the real part of the refractive index ``n = \\text{Re}(\\tilde{n})``.
"""
function refractive_index(m::IsotropicMedium, λ::Real = 0.0; T::Real = 293.15, P::Real = 101325.0, kwargs...)::Float64
    val = _eval_dispersion(m.dispersion, λ; T = T, P = P, kwargs...)
    return Float64(real(val))
end

"""
    extinction_coefficient(m::AbstractMedium, λ=0.0; kwargs...) -> Float64

Returns the extinction coefficient ``\\kappa = \\text{Im}(\\tilde{n})`` (positive for absorption, negative for gain).
"""
function extinction_coefficient(m::AbstractMedium, λ::Real = 0.0; kwargs...)::Float64
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
    medium_from(x) -> AbstractMedium

Converts a medium, constant refractive index, or dispersion function into an `AbstractMedium`.
"""
medium_from(m::AbstractMedium) = m
medium_from(n::Number) = IsotropicMedium(n)
medium_from(f) = IsotropicMedium(f)
