"""
    DiscreteRefractiveIndex{T}

Represents a incomplete set of dispersion data where for each exact wavelength one refractive index value is stored in the `data` field.
Can be called like a function `n = n(λ)`. Does not interpolate between data points.
Refer to [`RefractiveIndex`](@ref) for more information.
"""
struct DiscreteRefractiveIndex{T}
    data::Dict{T, T}
end

"""
    DiscreteRefractiveIndex(λs, n)

Creates a [`DiscreteRefractiveIndex`](@ref) dictionary where each wavelength in `λs` is mapped onto an exact exact refractive index in `ns`.

# Inputs

- `λs`: array of wavelengths
- `ns`: array of refractive indices
"""
function DiscreteRefractiveIndex(λs::AbstractArray{L}, ns::AbstractArray{N}) where {L, N}
    if length(λs) != length(ns)
        throw(ArgumentError("Number of wavelengths must match number of ref. indices"))
    end
    T = promote_type(L, N)
    d = Dict(zip(λs, ns))
    return DiscreteRefractiveIndex{T}(d)
end

(dri::DiscreteRefractiveIndex)(λ) = dri.data[λ]

"[`DiscreteRefractiveIndex`](@ref) passes test by default"
test_refractive_index_function(::DiscreteRefractiveIndex) = nothing

"""
    test_refractive_index_function(input)

Tests if `input` is callable with a single `Real` argument for the wavelength `λ` and
returns a single `Number` value (real or complex) for the refractive index `n`.
"""
function test_refractive_index_function(input)
    # Test function compat for the following types of λ
    Ts = [Int, Float32, Float64]
    try
        for T in Ts
            # Test if input accepts types
            answer = input(one(T))
            # Test if input returns single number value (real or complex)
            if !isa(answer, Number)
                error()
            end
        end
    catch
        error_msg = "Ref. index must be callable with a single Real argument and return a single Number result."
        throw(ArgumentError(error_msg))
    end
    return nothing
end

"""
    SellmeierEquation

A parametric type representing the **six-coefficient Sellmeier equation** for a transparent dielectric material.
The Sellmeier equation models the wavelength dependence of the refractive index `n(λ)` in the material's transparency window.
This type can provide `n(λ)` via a functor call, e.g.:

```julia
NBK7 = SellmeierEquation(...)
n_532 = NBK7(532e-9)
```

!!! info
    When initializing this type, the data should be provided with the classic **μm-based coefficients**.
    However, when calling the function, use SI-units, e.g. `532e-9` for 532 nm. 

# Fields

- `B1`, `B2`, `B3` : dimensionless Sellmeier coefficients.
- `C1`, `C2`, `C3` : squared resonance wavelengths in μm².
"""
struct SellmeierEquation{T<:Real}
    B1::T
    B2::T
    B3::T
    C1::T
    C2::T
    C3::T
end

"Returns the ref. index `n(λ)` for the six coefficient Sellmeier equation."
function (SE::SellmeierEquation)(λ)
    λ *= 1e6 # m to μm
    n² = 1 + (SE.B1*λ^2)/(λ^2 - SE.C1) +
             (SE.B2*λ^2)/(λ^2 - SE.C2) +
             (SE.B3*λ^2)/(λ^2 - SE.C3)
    return sqrt(n²)
end

"""
    AbsorbingMedium(n_real, α)

A callable medium wrapper representing a material with real base refractive index `n_real`
and bulk linear absorption coefficient `α` in [1/m] (Lambert-Beer: ``I(L) = I_0 e^{-\\alpha L}``).
`n_real` and `α` can be constant numbers or functions of wavelength `λ -> n_real(λ)`, `λ -> α(λ)`.

When called with wavelength `λ`, returns the complex refractive index:
```math
\\tilde{n}(\\lambda) = n_{\\mathrm{real}}(\\lambda) + i \\kappa(\\lambda) = n_{\\mathrm{real}}(\\lambda) + i \\frac{\\alpha(\\lambda) \\lambda}{4\\pi}
```

# Sign convention
Optics standard: ``\\tilde{n} = n + i\\kappa`` with ``\\kappa > 0`` for absorbing media (damping factor ``e^{-\\alpha L} = e^{-4\\pi \\kappa L / \\lambda}``).
"""
struct AbsorbingMedium{N, A}
    n_real::N
    alpha::A
end

function (m::AbsorbingMedium)(λ::Real)
    n_r = m.n_real isa Function ? real(m.n_real(λ)) : real(m.n_real)
    α_val = m.alpha isa Function ? m.alpha(λ) : m.alpha
    κ = (α_val * λ) / (4π)
    return complex(n_r, κ)
end

test_refractive_index_function(::AbsorbingMedium) = nothing

"""
    GainMedium(n_real, g)

A callable medium wrapper representing an active laser gain material with real base refractive index `n_real`
and small-signal gain coefficient `g` in [1/m] (where `g > 0` amplifies power as ``I(L) = I_0 e^{+g L}``).
`n_real` and `g` can be constant numbers or functions of wavelength `λ -> n_real(λ)`, `λ -> g(λ)`.

When called with wavelength `λ`, returns the complex refractive index:
```math
\\tilde{n}(\\lambda) = n_{\\mathrm{real}}(\\lambda) - i \\frac{g(\\lambda) \\lambda}{4\\pi}
```
"""
struct GainMedium{N, G}
    n_real::N
    gain::G
end

function (m::GainMedium)(λ::Real)
    n_r = m.n_real isa Function ? real(m.n_real(λ)) : real(m.n_real)
    g_val = m.gain isa Function ? m.gain(λ) : m.gain
    κ = -(g_val * λ) / (4π)
    return complex(n_r, κ)
end

test_refractive_index_function(::GainMedium) = nothing

"""
    RefractiveIndex

Union type that represents valid means to pass a refractive index `n` to e.g. [`AbstractObject`](@ref)s.
The core assumption is that:

1. the refractive index is callable with a **single** `Number` argument `λ` to represent the wavelength in [m]
2. the return value is a **single** `Number` value (real or complex) for the refractive index

Refer to e.g. [`DiscreteRefractiveIndex`](@ref), [`SellmeierEquation`](@ref), [`AbsorbingMedium`](@ref), and [`GainMedium`](@ref).
"""
const RefractiveIndex = Union{Function, DiscreteRefractiveIndex, SellmeierEquation, AbsorbingMedium, GainMedium}

"""
    real_refractive_index(n)
    real_refractive_index(obj, λ)

Returns the real part of the refractive index ``\\text{Re}(\\tilde{n})``.
"""
@inline real_refractive_index(n::Number) = real(n)
@inline real_refractive_index(obj, λ::Real) = real(refractive_index(obj, λ))

"""
    extinction_coefficient(n)
    extinction_coefficient(obj, λ)

Returns the imaginary part of the refractive index (extinction coefficient) ``\\kappa = \\text{Im}(\\tilde{n})``.
"""
@inline extinction_coefficient(n::Number) = imag(n)
@inline extinction_coefficient(obj, λ::Real) = imag(refractive_index(obj, λ))

"""
    absorption_coefficient(n, λ)
    absorption_coefficient(obj, λ)

Returns the bulk linear absorption coefficient ``\\alpha = \\frac{4\\pi \\kappa}{\\lambda}`` in [1/m].
Positive for absorbing media, negative for gain media.
"""
@inline absorption_coefficient(n::Number, λ::Real) = (4π * imag(n)) / λ
@inline absorption_coefficient(obj, λ::Real) = absorption_coefficient(refractive_index(obj, λ), λ)

"""
    bulk_attenuation_factor(n, λ, L)

Computes the power and field attenuation factors `(att_power, att_field)` for light propagating
over physical path length `L` in a medium with refractive index `n` at wavelength `λ`.
- For non-absorbing media (``\\text{Im}(n) = 0``), returns `(1.0, 1.0)`.
- For absorbing media (``\\kappa = \\text{Im}(n) > 0``), ``\\text{att\\_power} = e^{-\\alpha L}`` and ``\\text{att\\_field} = e^{-\\frac{\\alpha}{2} L}``.
- For gain media (``\\kappa < 0``), ``\\text{att\\_power} = e^{+g L}`` and ``\\text{att\\_field} = e^{+\\frac{g}{2} L}``.
"""
@inline function bulk_attenuation_factor(n::Number, λ::Real, L::Real)
    κ = imag(n)
    if iszero(κ)
        T = typeof(real(n) * L)
        return (one(T), one(T))
    end
    arg_field = -2π * κ * L / λ
    att_field = exp(arg_field)
    att_power = exp(2 * arg_field)
    return (att_power, att_field)
end