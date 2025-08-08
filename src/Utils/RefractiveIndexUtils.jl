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
returns a single `Real` value for the refractive index `n`.
"""
function test_refractive_index_function(input)
    # Test function compat for the following types of λ
    Ts = [Int, Float32, Float64]
    try
        for T in Ts
            # Test if input accepts types
            answer = input(one(T))
            # Test if input returns single real value
            if !isa(answer, Real)
                error()
            end
        end
    catch
        error_msg = "Ref. index must be callable with a single Real argument and return a single real result."
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
    RefractiveIndex

Union type that represents valid means to pass a refractive index `n` to e.g. [`AbstractObject`](@ref)s.
The core assumption is that:

1. the refractive index is callable with a **single** `Number` argument `λ` to represent the wavelength in [m]
2. the return value is a **single** `Number` value for the refractive index

Refer to e.g. [`DiscreteRefractiveIndex`](@ref).
"""
const RefractiveIndex = Union{Function, DiscreteRefractiveIndex, SellmeierEquation}