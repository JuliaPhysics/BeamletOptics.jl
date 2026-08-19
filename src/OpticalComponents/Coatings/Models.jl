# Coatings models, traits, and structs

"""
    AbstractCoating{T} <: AbstractObject{T}

Base abstract type for standalone coating *objects* — thin coated boundaries that exist as
independent optical components in a system. This is distinct from coating *models* (such as
`SimpleARCoating`, `ThinFilmCoating`, etc.) which define the physical interaction properties
but are not optical objects themselves.
"""
abstract type AbstractCoating{T} <: AbstractObject{T} end

is_thin_interface(::AbstractCoating) = true

# Coating behavior traits
"""
    CoatingBehavior

Abstract trait representing the physical interaction behavior of a coating boundary.
The core physics engine dispatches ray and beamlet propagation logic based on this trait.
Subtypes are `Transmissive`, `Reflective`, `Splitting`, and `Absorptive`.
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

"""
    Absorptive <: CoatingBehavior

Trait indicating a coating boundary that absorbs incident light. Rays hitting an absorptive
coating are terminated (the interaction returns `nothing`).
"""
struct Absorptive <: CoatingBehavior end

# Coating models
"""
    coating_behavior(coating, ray)
    coating_behavior(coating)

Query the physical interaction behavior trait of a coating (`Transmissive`, `Reflective`, or `Splitting`).
The two-argument version allows custom coating models to dynamically determine their behavior based on the incident ray parameters (e.g. angle of incidence, wavelength).
Falls back to `coating_behavior(coating)` for static behaviors.
"""
coating_behavior(c, ray) = coating_behavior(c)

# Surface and Coating models, traits, and structs

"""
    AbstractSurfaceModel

Base abstract type for all optical surface physics models (coatings, diffusers, diffractive surfaces, absorbing apertures, etc.).
"""
abstract type AbstractSurfaceModel end

"""
    AbstractCoatingModel <: AbstractSurfaceModel

Base abstract type for thin-film, Jones, and dielectric coating models.
"""
abstract type AbstractCoatingModel <: AbstractSurfaceModel end

"""
    Uncoated

Represents an uncoated dielectric boundary. Behaves as `Transmissive` by default.
"""
struct Uncoated <: AbstractCoatingModel end
coating_behavior(::Uncoated) = Transmissive()

"""
    is_coated(coating_model)

Predicate indicating whether a coating model represents an active coating layer rather than an uncoated boundary (`Uncoated`).
"""
is_coated(::Uncoated) = false
is_coated(::AbstractCoatingModel) = true
is_coated(::Any) = true

@inline get_coating_behavior(c) = coating_behavior(c)

"""
    SimpleARCoating(R::Real = 0.0)

A simplified Anti-Reflective (AR) coating model defined by a constant power reflectance `R`.
Behaves as `Transmissive`.

!!! warning "Angle Independence"
    This model returns constant coefficients regardless of angle of incidence, wavelength,
    or refractive index ratio. For angle-dependent AR behavior, use [`ThinFilmCoating`](@ref).
"""
struct SimpleARCoating{T <: Real} <: AbstractCoatingModel
    R::T
end
SimpleARCoating() = SimpleARCoating(0.0)
coating_behavior(::SimpleARCoating) = Transmissive()

"""
    SimpleHRCoating(R::Real = 1.0)

A simplified Highly Reflective (HR) coating model defined by a constant power reflectance `R`.
Behaves as `Reflective`.

!!! warning "Angle Independence"
    This model returns constant coefficients regardless of angle of incidence, wavelength,
    or refractive index ratio. For angle-dependent HR behavior, use [`ThinFilmCoating`](@ref).
"""
struct SimpleHRCoating{T <: Real} <: AbstractCoatingModel
    R::T
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
struct SimpleBeamsplitterCoating{C <: Complex} <: AbstractCoatingModel
    rs::C
    rp::C
    ts::C
    tp::C
end
function SimpleBeamsplitterCoating(rs::Number, rp::Number, ts::Number, tp::Number)
    T = promote_type(typeof(complex(rs)), typeof(complex(rp)), typeof(complex(ts)), typeof(complex(tp)))
    return SimpleBeamsplitterCoating{T}(T(rs), T(rp), T(ts), T(tp))
end
coating_behavior(::SimpleBeamsplitterCoating) = Splitting()

"""
    JonesCoating(jones_trans, jones_refl = XZBasis(0,0,0,0); behavior = Transmissive())

A fully generalized coating defined by its transmitted and reflected Jones matrices.
The Jones matrices can be constants (e.g. `SPBasis`, `XZBasis`) or functions taking wavelength `λ` and returning a `JonesMatrix`.

!!! note "Transmission Scaling & Power Conservation"
    When transmitting across an interface with different refractive indices (`n1 != n2`),
    the transmitted Jones matrix is automatically scaled by the optical impedance factor
    `sqrt((n1 * cos(θi)) / (n2 * cos(θt)))`. This ensures that the ray's electric field magnitude
    `|E|²` directly maps to conserved optical power (Poynting flux).
"""
struct JonesCoating{JT, JR} <: AbstractCoatingModel
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

A physical thin-film interference coating model defined by a stack of refractive indices `ns` and physical thicknesses `ds` (in same unit as wavelength, typically m).
Indices can be static real numbers or dispersion functions `f(λ) -> n`.
Calculates exact complex Fresnel coefficients using the characteristic transfer matrix method.

!!! tip "Empty Layer Stack"
    Passing empty vectors `ThinFilmCoating(Float64[], Float64[])` creates a bare interface
    that yields the standard Fresnel coefficients for the surrounding media. This is useful
    for modeling splitting at an uncoated index boundary.
"""
struct ThinFilmCoating{N <: Tuple, D <: Tuple, B <: CoatingBehavior} <: AbstractCoatingModel
    ns::N
    ds::D
    behavior::B
end
function ThinFilmCoating(
        n::Union{Number, Function}, d::Real; behavior::CoatingBehavior = Transmissive())
    return ThinFilmCoating((n,), (d,), behavior)
end
function ThinFilmCoating(ns::Union{AbstractVector, Tuple}, ds::Union{AbstractVector, Tuple};
        behavior::CoatingBehavior = Transmissive())
    ns_t = Tuple(ns)
    ds_t = Tuple(ds)
    return ThinFilmCoating(ns_t, ds_t, behavior)
end
coating_behavior(c::ThinFilmCoating) = c.behavior

# fresnel_coefficients interface
fresnel_coefficients(::Uncoated, θi, λ, n1, n2) = fresnel_coefficients(θi, n2 / n1)

function get_jones_matrix(::Uncoated, θi, λ, n1, n2, is_reflected; from_front::Bool = true, local_p = nothing)
    rs, rp, ts, tp = fresnel_coefficients(θi, n2 / n1)
    if is_reflected
        return SPBasis(-rs, 0, 0, rp)
    else
        return SPBasis(ts, 0, 0, tp)
    end
end

"""
    unpolarized_reflectance(J)

Compute the average power reflectance from a Jones reflection matrix for unpolarized light.

!!! note "Absorption"
    For non-absorbing coatings, `T = 1 - R` holds exactly. For absorbing coatings (complex
    refractive index layers in a `ThinFilmCoating`), the absorbed fraction is accounted for such
    that `T + R < 1`. Both `Ray` and `PolarizedRay` derive `T` via `coating_transmittance`
    and `R` via `coating_reflectance` so that absorbed power `A = 1 - R - T` is properly deducted.
"""
@inline unpolarized_reflectance(J) = 0.5 * (abs2(J[1, 1]) + abs2(J[1, 2]) + abs2(J[2, 1]) + abs2(J[2, 2]))

"""
    unpolarized_transmittance(J, n1=1.0, n2=1.0, θi=0.0)

Compute the average power transmittance from a Jones transmission matrix `J` for unpolarized light,
accounting for the medium admittance ratio `Re(n2*cosθt) / Re(n1*cosθi)`.
"""
@inline function unpolarized_transmittance(J, n1::Number = 1.0, n2::Number = 1.0, θi::Real = 0.0)
    sinθt2 = (n1 / n2)^2 * sin(θi)^2
    cosθt = sqrt(Complex(1.0 - sinθt2))
    cosθi = cos(θi)
    denom = real(n1 * cosθi)
    factor = denom > 0 ? max(0.0, real(n2 * cosθt) / denom) : 1.0
    return factor * 0.5 * (abs2(J[1, 1]) + abs2(J[1, 2]) + abs2(J[2, 1]) + abs2(J[2, 2]))
end

"""
    coating_reflectance(coating, θi, λ, n1, n2; from_front=true, local_p=nothing)

Compute the unpolarized power reflectance R ∈ [0, 1] of a coating model at angle of incidence `θi`.
"""
function coating_reflectance(c, θi, λ, n1, n2; from_front::Bool = true, local_p = nothing)
    J_r = get_jones_matrix(c, θi, λ, n1, n2, true; from_front = from_front, local_p = local_p)
    return clamp(unpolarized_reflectance(J_r), 0.0, 1.0)
end

"""
    coating_transmittance(coating, θi, λ, n1, n2; from_front=true, local_p=nothing)

Compute the unpolarized power transmittance T ∈ [0, 1] of a coating model at angle of incidence `θi`.
"""
function coating_transmittance(c, θi, λ, n1, n2; from_front::Bool = true, local_p = nothing)
    J_t = get_jones_matrix(c, θi, λ, n1, n2, false; from_front = from_front, local_p = local_p)
    return clamp(unpolarized_transmittance(J_t, n1, n2, θi), 0.0, 1.0)
end

# Shared implementation for simplified constant-reflectance coating models
const SimpleReflectanceCoating = Union{SimpleARCoating, SimpleHRCoating}

function fresnel_coefficients(c::SimpleReflectanceCoating, θi, λ, n1, n2)
    R = c.R
    T = 1.0 - R
    rs = -sqrt(R)
    rp = sqrt(R)
    ts = sqrt(T)
    tp = sqrt(T)
    return rs, rp, ts, tp
end

function get_jones_matrix(
        c::SimpleReflectanceCoating, θi, λ, n1, n2, is_reflected; from_front::Bool = true, local_p = nothing)
    rs, rp, ts, tp = fresnel_coefficients(c, θi, λ, n1, n2)
    if is_reflected
        return SPBasis(-rs, 0, 0, rp)
    else
        return SPBasis(ts, 0, 0, tp)
    end
end

coating_reflectance(c::SimpleReflectanceCoating, θi, λ, n1, n2; from_front::Bool = true, local_p = nothing) = c.R
coating_transmittance(c::SimpleReflectanceCoating, θi, λ, n1, n2; from_front::Bool = true, local_p = nothing) = 1.0 - c.R

fresnel_coefficients(c::SimpleBeamsplitterCoating, θi, λ, n1, n2) = (c.rs, c.rp, c.ts, c.tp)

function get_jones_matrix(c::SimpleBeamsplitterCoating, θi, λ, n1, n2, is_reflected; from_front::Bool = true, local_p = nothing)
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
        c::JonesCoating, θi, λ, n1, n2, is_reflected; from_front::Bool = true, local_p = nothing)
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

"""
    coating_properties(coating, λ)
    coating_properties(coating, λ, local_p)

Query coating optical properties at wavelength `λ` and local spatial position `local_p`.
Falls back to `coating_properties(coating, λ)` for spatially uniform coatings.
"""
coating_properties(c::AbstractCoatingModel, λ::Real) = ()

@inline _eval_refractive_index(n::Function, λ::Real) = complex(n(λ))
@inline _eval_refractive_index(n, λ::Real) = complex(n)

function coating_properties(c::ThinFilmCoating, λ::Real)
    n_vals = map(n -> _eval_refractive_index(n, λ), c.ns)
    return n_vals, c.ds
end

coating_properties(c::AbstractCoatingModel, λ::Real, local_p::AbstractVector) = coating_properties(c, λ)
coating_properties(c::AbstractCoatingModel, λ::Real, ::Nothing) = coating_properties(c, λ)

function _fresnel_coefficients_matrix(n_vals, d_vals, θi::Real, λ::Real,
        n1::Number, n2::Number; from_front::Bool = true)
    sinθi = sin(θi)
    cosθi = cos(θi)

    cosθt = sqrt(Complex(1.0 - (n1 / n2)^2 * sinθi^2))

    η1s = n1 * cosθi
    η2s = n2 * cosθt

    η1p = n1 / cosθi
    η2p = n2 / cosθt

    # Initialize characteristic transfer matrices to Identity
    m11_s = complex(1.0, 0.0)
    m12_s = complex(0.0, 0.0)
    m21_s = complex(0.0, 0.0)
    m22_s = complex(1.0, 0.0)

    m11_p = complex(1.0, 0.0)
    m12_p = complex(0.0, 0.0)
    m21_p = complex(0.0, 0.0)
    m22_p = complex(1.0, 0.0)

    N_layers = length(n_vals)
    layer_indices = from_front ? (1:N_layers) : (N_layers:-1:1)

    for j in layer_indices
        nj_val = n_vals[j]
        dj = d_vals[j]

        cosθj = sqrt(Complex(1.0 - (n1 / nj_val)^2 * sinθi^2))
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

function fresnel_coefficients(c::ThinFilmCoating, θi::Real, λ::Real,
        n1::Number, n2::Number; from_front::Bool = true)
    ns, ds = coating_properties(c, λ)
    return _fresnel_coefficients_matrix(ns, ds, θi, λ, n1, n2; from_front = from_front)
end

function get_jones_matrix(
        c::ThinFilmCoating, θi, λ, n1, n2, is_reflected; from_front::Bool = true, local_p = nothing)
    rs, rp, ts, tp = fresnel_coefficients(c, θi, λ, n1, n2; from_front = from_front)
    if is_reflected
        return SPBasis(-rs, 0, 0, rp)
    else
        return SPBasis(ts, 0, 0, tp)
    end
end

"""
    GradedThinFilmCoating(thickness_func, base_coating::ThinFilmCoating)

A spatially graded thin-film coating where layer thicknesses are scaled spatially by
`thickness_func(local_p)` at hit point `local_p`.
"""
struct GradedThinFilmCoating{F, C <: ThinFilmCoating} <: AbstractCoatingModel
    thickness_func::F
    base_coating::C
end

coating_behavior(g::GradedThinFilmCoating) = coating_behavior(g.base_coating)

function coating_properties(g::GradedThinFilmCoating, λ::Real, local_p::AbstractVector)
    factor = g.thickness_func(local_p)
    ns, ds = coating_properties(g.base_coating, λ)
    return ns, ds .* factor
end

function get_jones_matrix(
        g::GradedThinFilmCoating, θi, λ, n1, n2, is_reflected; from_front::Bool = true, local_p = nothing)
    factor = isnothing(local_p) ? 1.0 : g.thickness_func(local_p)
    ns, ds = coating_properties(g.base_coating, λ)
    ds_graded = ds .* factor
    rs, rp, ts, tp = _fresnel_coefficients_matrix(ns, ds_graded, θi, λ, n1, n2; from_front = from_front)
    if is_reflected
        return SPBasis(-rs, 0, 0, rp)
    else
        return SPBasis(ts, 0, 0, tp)
    end
end

"""
    CompositeSurfaceModel(models::AbstractSurfaceModel...)

Composes multiple optical surface physics models sequentially.
Reflectance and transmittance matrices are multiplied across sub-models.
"""
struct CompositeSurfaceModel{T <: Tuple{Vararg{AbstractSurfaceModel}}} <: AbstractSurfaceModel
    models::T
end

CompositeSurfaceModel(models::AbstractSurfaceModel...) = CompositeSurfaceModel(models)

@inline _behavior_precedence(::Absorptive) = 4
@inline _behavior_precedence(::Splitting) = 3
@inline _behavior_precedence(::Reflective) = 2
@inline _behavior_precedence(::Transmissive) = 1

@inline _dominant_behavior(b1::CoatingBehavior, b2::CoatingBehavior) =
    _behavior_precedence(b1) >= _behavior_precedence(b2) ? b1 : b2

function coating_behavior(comp::CompositeSurfaceModel)
    return isempty(comp.models) ? Transmissive() : reduce(_dominant_behavior, map(coating_behavior, comp.models))
end

function get_jones_matrix(
        comp::CompositeSurfaceModel, θi, λ, n1, n2, is_reflected; from_front::Bool = true, local_p = nothing)
    J_tot = get_jones_matrix(comp.models[1], θi, λ, n1, n2, is_reflected; from_front = from_front, local_p = local_p)
    for i in 2:length(comp.models)
        J_i = get_jones_matrix(comp.models[i], θi, λ, n1, n2, is_reflected; from_front = from_front, local_p = local_p)
        J_tot = J_i * J_tot
    end
    return J_tot
end

# Base.show methods for coating models
Base.show(io::IO, ::MIME"text/plain", ::Uncoated) = print(io, "Uncoated")
Base.show(io::IO, ::MIME"text/plain", c::SimpleARCoating) = print(io, "SimpleARCoating(R = ", c.R, ")")
Base.show(io::IO, ::MIME"text/plain", c::SimpleHRCoating) = print(io, "SimpleHRCoating(R = ", c.R, ")")
Base.show(io::IO, ::MIME"text/plain", c::SimpleBeamsplitterCoating) = print(io, "SimpleBeamsplitterCoating(rs = ", c.rs, ", rp = ", c.rp, ", ts = ", c.ts, ", tp = ", c.tp, ")")
Base.show(io::IO, ::MIME"text/plain", c::JonesCoating) = print(io, "JonesCoating(behavior = ", c.behavior, ")")
Base.show(io::IO, ::MIME"text/plain", c::ThinFilmCoating) = print(io, "ThinFilmCoating(", length(c.ns), " layers, behavior = ", c.behavior, ")")
Base.show(io::IO, ::MIME"text/plain", c::GradedThinFilmCoating) = print(io, "GradedThinFilmCoating(base = ", c.base_coating, ")")
Base.show(io::IO, ::MIME"text/plain", c::CompositeSurfaceModel) = print(io, "CompositeSurfaceModel(", length(c.models), " models)")

