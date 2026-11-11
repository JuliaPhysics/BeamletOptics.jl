"""
    AbstractRefractiveIndex

A generic representation of the refractive index `n(λ)` of a material, passed as `n` to e.g.
[`Lens`](@ref), [`Prism`](@ref) and the beamsplitters. See [`ConstantRefractiveIndex`](@ref),
[`DiscreteRefractiveIndex`](@ref) and [`SellmeierEquation`](@ref).

# Implementation reqs.

Subtypes of `AbstractRefractiveIndex` must implement the following:

## Functions

- `(n::AbstractRefractiveIndex)(λ)`: returns the refractive index as a `Real` for the wavelength
  `λ` in \\[m\\], for any `Real` type of `λ` within the domain of the model
- `extinction_coefficient(n, λ)`: optional, the extinction coefficient κ of an absorbing material,
  the imaginary part of its complex index `n(λ) + iκ`; 0 by default, see
  [`AbsorbingRefractiveIndex`](@ref)

# Additional information

!!! info "Domain"
    A model may be defined for some wavelengths only, e.g. a [`DiscreteRefractiveIndex`](@ref)
    only at its stored wavelengths.
"""
abstract type AbstractRefractiveIndex end

(n::AbstractRefractiveIndex)(λ) = error(lazy"refractive index n(λ) not defined for $(typeof(n))")

"""
    extinction_coefficient(n::AbstractRefractiveIndex, λ)

Returns the extinction coefficient κ of the material of the refractive index `n` at the wavelength `λ`
in \\[m\\]: the imaginary part of its complex index `n(λ) + iκ`, with the field `∝ exp(i(kz - ωt))`
of BMO. A field decays by `exp(-2πκz/λ)`. Zero for a transparent material (default).

!!! info
    The thin layers of a [`ThinFilmCoating`](@ref) absorb with κ; the volume of a body does not
    attenuate the light.
"""
extinction_coefficient(n::AbstractRefractiveIndex, λ::Real) = zero(float(typeof(n(λ))))

"""
    AbsorbingRefractiveIndex(n, κ)

The refractive index of an absorbing material: the real index `n`, an
[`AbstractRefractiveIndex`](@ref), and the extinction coefficient `κ`, a number or an
`AbstractRefractiveIndex` that returns κ for a wavelength, e.g. a [`DiscreteRefractiveIndex`](@ref)
of measured data. See [`extinction_coefficient`](@ref).

```julia
silver = AbsorbingRefractiveIndex(ConstantRefractiveIndex(0.135), 3.985)   # at 633 nm
```
"""
struct AbsorbingRefractiveIndex{
    N <: AbstractRefractiveIndex, K <: Union{Real, AbstractRefractiveIndex}} <: AbstractRefractiveIndex
    n::N
    κ::K
end

(ari::AbsorbingRefractiveIndex)(λ) = ari.n(λ)

extinction_coefficient(ari::AbsorbingRefractiveIndex, λ::Real) = value_at(ari.κ, λ)

# A constant or a function of the wavelength `λ`
value_at(x::Real, ::Real) = x
value_at(f::AbstractRefractiveIndex, λ::Real) = f(λ)

"""Returns the complex refractive index `n(λ) + iκ` of the material of the refractive index `n`."""
complex_index(n::AbstractRefractiveIndex, λ::Real) = complex(n(λ), extinction_coefficient(n, λ))

"""
    ConstantRefractiveIndex{T}

Represents a refractive index `n` that is the same for all wavelengths, i.e. without dispersion.
Can be called like a function `n = n(λ)`.
"""
struct ConstantRefractiveIndex{T <: Real} <: AbstractRefractiveIndex
    n::T
end

(cri::ConstantRefractiveIndex)(λ) = cri.n

"""
    DiscreteRefractiveIndex{T}

Represents a incomplete set of dispersion data where for each exact wavelength one refractive index value is stored in the `data` field.
Can be called like a function `n = n(λ)`. Does not interpolate between data points.
Refer to [`AbstractRefractiveIndex`](@ref) for more information.
"""
struct DiscreteRefractiveIndex{T} <: AbstractRefractiveIndex
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
struct SellmeierEquation{T<:Real} <: AbstractRefractiveIndex
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
