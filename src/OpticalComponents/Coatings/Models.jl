# Coatings models, traits, and structs

"""
    AbstractCoating{T} <: AbstractObject{T}

Base abstract type for all optical coating models and interface properties.
"""
abstract type AbstractCoating{T} <: AbstractObject{T} end

is_thin_interface(::AbstractCoating) = true

# Coating behavior traits
"""
    CoatingBehavior

Abstract trait representing the physical interaction behavior of a coating boundary.
The core physics engine dispatches ray and beamlet propagation logic based on this trait.
Subtypes are `Transmissive`, `Reflective`, and `Splitting`.
"""
abstract type CoatingBehavior end

"""
    Transmissive <: CoatingBehavior

Trait indicating a coating boundary primarily intended to be traversed by light.
"""
struct Transmissive <: CoatingBehavior end

"""
    Reflective <: CoatingBehavior

Trait indicating a coating boundary intended to reflect light back into the incident medium.
"""
struct Reflective <: CoatingBehavior end

"""
    Splitting <: CoatingBehavior

Trait indicating a beam-splitting coating boundary that explicitly splits each incoming ray into both a transmitted and a reflected child ray.
"""
struct Splitting <: CoatingBehavior end

# Coating models
"""
    coating_behavior(coating, ray)
    coating_behavior(coating)

Query the physical interaction behavior trait of a coating (`Transmissive`, `Reflective`, or `Splitting`).
The two-argument version allows custom coating models to dynamically determine their behavior based on the incident ray parameters (e.g. angle of incidence, wavelength).
Falls back to `coating_behavior(coating)` for static behaviors.
"""
coating_behavior(c, ray) = coating_behavior(c)

"""
    Uncoated

Represents an uncoated dielectric boundary. Behaves as `Transmissive` by default.
"""
struct Uncoated end
coating_behavior(::Uncoated) = Transmissive()

"""
    SimpleARCoating(R::Float64 = 0.0)

A simplified Anti-Reflective (AR) coating model defined by a constant power reflectance `R`.
Behaves as `Transmissive`.
"""
struct SimpleARCoating
    R::Float64
end
SimpleARCoating() = SimpleARCoating(0.0)
coating_behavior(::SimpleARCoating) = Transmissive()

"""
    SimpleHRCoating(R::Float64 = 1.0)

A simplified Highly Reflective (HR) coating model defined by a constant power reflectance `R`.
Behaves as `Reflective`.
"""
struct SimpleHRCoating
    R::Float64
end
SimpleHRCoating() = SimpleHRCoating(1.0)
coating_behavior(::SimpleHRCoating) = Reflective()

"""
    SimpleBeamsplitterCoating(rs, rp, ts, tp)

A simple beam-splitting coating defined by fixed complex amplitude reflection (`rs`, `rp`) and transmission (`ts`, `tp`) coefficients.
Behaves as `Splitting`.

!!! note "Note on Transmission Coefficients"
    The provided `ts` and `tp` values act effectively as the square roots of the power transmissivity,
    rather than pure electric field amplitude transmission coefficients. The implementation implicitly
    scales the resulting Jones matrix by the boundary impedance factor to ensure that ray intensity
    strictly maps to optical power.
"""
struct SimpleBeamsplitterCoating
    rs::ComplexF64
    rp::ComplexF64
    ts::ComplexF64
    tp::ComplexF64
end
coating_behavior(::SimpleBeamsplitterCoating) = Splitting()

"""
    JonesCoating(jones_trans, jones_refl = XZBasis(0,0,0,0); behavior = Transmissive())

A fully generalized coating defined by its transmitted and reflected Jones matrices.
The Jones matrices can be constants (e.g. `SPBasis`) or functions taking wavelength `λ` and returning a `JonesMatrix`.
"""
struct JonesCoating{JT, JR}
    jones_trans::JT
    jones_refl::JR
    behavior::CoatingBehavior
end
JonesCoating(trans) = JonesCoating(trans, XZBasis(0, 0, 0, 0), Transmissive())
JonesCoating(trans, refl; behavior = Transmissive()) = JonesCoating(trans, refl, behavior)
coating_behavior(c::JonesCoating) = c.behavior

"""
    ThinFilmCoating(ns::Vector, ds::Vector{<:Real}; behavior = Transmissive())
    ThinFilmCoating(n, d::Real; behavior = Transmissive())

A physical thin-film interference coating model defined by a stack of refractive indices `ns` and physical thicknesses `ds` (in mm).
Indices can be static real numbers or dispersion functions `f(λ) -> n`.
Calculates exact complex Fresnel coefficients using the characteristic transfer matrix method.
"""
struct ThinFilmCoating{N <: Tuple, D <: Tuple, B <: CoatingBehavior}
    ns::N
    ds::D
    behavior::B
end
function ThinFilmCoating(
        n::Union{Number, Function}, d::Real; behavior::CoatingBehavior = Transmissive())
    return ThinFilmCoating{typeof((n,)), typeof((Float64(d),)), typeof(behavior)}((n,), (Float64(d),), behavior)
end
function ThinFilmCoating(ns::Union{AbstractVector, Tuple}, ds::Union{AbstractVector, Tuple};
        behavior::CoatingBehavior = Transmissive())
    ns_t = Tuple(ns)
    ds_t = Tuple(Float64.(ds))
    return ThinFilmCoating{typeof(ns_t), typeof(ds_t), typeof(behavior)}(ns_t, ds_t, behavior)
end
coating_behavior(c::ThinFilmCoating) = c.behavior

# fresnel_coefficients interface
fresnel_coefficients(::Uncoated, θi, λ, n1, n2) = fresnel_coefficients(θi, n2 / n1)

function get_jones_matrix(::Uncoated, θi, λ, n1, n2, is_reflected; from_front::Bool = true)
    rs, rp, ts, tp = fresnel_coefficients(θi, n2 / n1)
    if is_reflected
        return SPBasis(-rs, 0, 0, rp)
    else
        return SPBasis(ts, 0, 0, tp)
    end
end

function fresnel_coefficients(c::SimpleARCoating, θi, λ, n1, n2)
    R = c.R
    T = 1.0 - R
    rs = -sqrt(R)
    rp = sqrt(R)
    ts = sqrt(T)
    tp = sqrt(T)
    return rs, rp, ts, tp
end

function get_jones_matrix(
        c::SimpleARCoating, θi, λ, n1, n2, is_reflected; from_front::Bool = true)
    rs, rp, ts, tp = fresnel_coefficients(c, θi, λ, n1, n2)
    if is_reflected
        return SPBasis(-rs, 0, 0, rp)
    else
        return SPBasis(ts, 0, 0, tp)
    end
end

function fresnel_coefficients(c::SimpleHRCoating, θi, λ, n1, n2)
    R = c.R
    T = 1.0 - R
    rs = -sqrt(R)
    rp = sqrt(R)
    ts = sqrt(T)
    tp = sqrt(T)
    return rs, rp, ts, tp
end

function get_jones_matrix(
        c::SimpleHRCoating, θi, λ, n1, n2, is_reflected; from_front::Bool = true)
    rs, rp, ts, tp = fresnel_coefficients(c, θi, λ, n1, n2)
    if is_reflected
        return SPBasis(-rs, 0, 0, rp)
    else
        return SPBasis(ts, 0, 0, tp)
    end
end

fresnel_coefficients(c::SimpleBeamsplitterCoating, θi, λ, n1, n2) = (c.rs, c.rp, c.ts, c.tp)

function get_jones_matrix(c::SimpleBeamsplitterCoating, θi, λ, n1, n2, is_reflected; from_front::Bool=true)
    if is_reflected
        rs = from_front ? c.rs : -c.rs
        rp = from_front ? c.rp : -c.rp
        return SPBasis(-rs, 0, 0, rp)
    else
        # Scale transmission coefficients to conserve power across index-mismatched boundaries
        sinθt² = (n1 / n2)^2 * sin(θi)^2
        if sinθt² >= 1.0
            return SPBasis(0.0, 0.0, 0.0, 0.0)
        else
            cosθt = sqrt(1.0 - sinθt²)
            scale = sqrt((n1 * cos(θi)) / (n2 * cosθt))
            return SPBasis(c.ts * scale, 0, 0, c.tp * scale)
        end
    end
end

@inline _eval_jones(j::GlobalJonesBasis, ::Real) = j
@inline _eval_jones(j::LocalJonesBasis, ::Real) = j
@inline _eval_jones(j::Function, λ::Real) = j(λ)

function get_jones_matrix(
        c::JonesCoating, θi, λ, n1, n2, is_reflected; from_front::Bool = true)
    if is_reflected
        return _eval_jones(c.jones_refl, λ)
    else
        J = _eval_jones(c.jones_trans, λ)
        if n1 ≈ n2
            return J
        else
            sinθt² = (n1 / n2)^2 * sin(θi)^2
            if sinθt² >= 1.0
                return typeof(J)(zero(static_data(J)))
            else
                cosθt = sqrt(1.0 - sinθt²)
                scale = sqrt((n1 * cos(θi)) / (n2 * cosθt))
                return scale * J
            end
        end
    end
end

function fresnel_coefficients(c::ThinFilmCoating, θi::Real, λ::Real,
        n1::Number, n2::Number; from_front::Bool = true)
    sinθi = sin(θi)
    cosθi = cos(θi)

    cosθt = sqrt(Complex(1.0 - (n1 / n2)^2 * sinθi^2))

    η1s = n1 * cosθi
    η2s = n2 * cosθt

    η1p = n1 / cosθi
    η2p = n2 / cosθt

    # Initialize characteristic transfer matrices to Identity
    # Using scalar variables rather than matrices to prevent heap allocation in tight raymarching loops
    m11_s = complex(1.0, 0.0)
    m12_s = complex(0.0, 0.0)
    m21_s = complex(0.0, 0.0)
    m22_s = complex(1.0, 0.0)

    m11_p = complex(1.0, 0.0)
    m12_p = complex(0.0, 0.0)
    m21_p = complex(0.0, 0.0)
    m22_p = complex(1.0, 0.0)

    N_layers = length(c.ns)
    layer_indices = from_front ? (1:N_layers) : (N_layers:-1:1)

    for j in layer_indices
        nj_val = ComplexF64(c.ns[j] isa Function ? c.ns[j](λ) : c.ns[j])
        dj = c.ds[j]

        cosθj = sqrt(1.0 - (n1 / nj_val)^2 * sinθi^2)
        δj = (2π / λ) * nj_val * dj * cosθj

        ηjs = nj_val * cosθj
        ηjp = nj_val / cosθj

        cos_δj = cos(δj)
        sin_δj = sin(δj)

        # Mjs components
        Mjs11 = cos_δj
        Mjs12 = im * sin_δj / ηjs
        Mjs21 = im * ηjs * sin_δj
        Mjs22 = cos_δj

        # Ms = Ms * Mjs
        new_m11_s = m11_s * Mjs11 + m12_s * Mjs21
        new_m12_s = m11_s * Mjs12 + m12_s * Mjs22
        new_m21_s = m21_s * Mjs11 + m22_s * Mjs21
        new_m22_s = m21_s * Mjs12 + m22_s * Mjs22

        m11_s, m12_s, m21_s, m22_s = new_m11_s, new_m12_s, new_m21_s, new_m22_s

        # Mjp components
        Mjp11 = cos_δj
        Mjp12 = im * sin_δj / ηjp
        Mjp21 = im * ηjp * sin_δj
        Mjp22 = cos_δj

        # Mp = Mp * Mjp
        new_m11_p = m11_p * Mjp11 + m12_p * Mjp21
        new_m12_p = m11_p * Mjp12 + m12_p * Mjp22
        new_m21_p = m21_p * Mjp11 + m22_p * Mjp21
        new_m22_p = m21_p * Mjp12 + m22_p * Mjp22

        m11_p, m12_p, m21_p, m22_p = new_m11_p, new_m12_p, new_m21_p, new_m22_p
    end

    Bs = m11_s + η2s * m12_s
    Cs = m21_s + η2s * m22_s
    rs = (η1s * Bs - Cs) / (η1s * Bs + Cs)
    ts = 2.0 * η1s / (η1s * Bs + Cs)

    Bp = m11_p + η2p * m12_p
    Cp = m21_p + η2p * m22_p
    rp = (η1p * Bp - Cp) / (η1p * Bp + Cp)
    tp = (2.0 * η1p / (η1p * Bp + Cp)) * (cosθi / cosθt)

    return rs, rp, ts, tp
end

function get_jones_matrix(
        c::ThinFilmCoating, θi, λ, n1, n2, is_reflected; from_front::Bool = true)
    rs, rp, ts, tp = fresnel_coefficients(c, θi, λ, n1, n2; from_front = from_front)
    if is_reflected
        return SPBasis(-rs, 0, 0, rp)
    else
        return SPBasis(ts, 0, 0, tp)
    end
end
