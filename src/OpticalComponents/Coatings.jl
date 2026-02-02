"""
    AbstractCoating

Abstract type for all optical coatings.
Subtypes must implement `coating_coefficients(coating, n_in, n_out, λ, θ)`.
"""
abstract type AbstractCoating end

"""
    Uncoated <: AbstractCoating

Represents an uncoated surface. Uses standard Fresnel equations.
"""
struct Uncoated <: AbstractCoating end

"""
    get_coating_at(coating, optic, point, normal) -> AbstractCoating

Returns the coating at a specific point on the surface.
For standard uniform coatings, this simply returns the coating object itself.
For spatially varying coatings, it returns the specific sub-coating at that point.

# Arguments
- `coating`: the coating object
- `optic`: the optical component being interacted with
- `point`: the intersection point in **global coordinates** (will be transformed if needed by selector)
- `normal`: the surface normal in **global coordinates**
"""
function get_coating_at(coating::AbstractCoating, optic, point, normal)
    return coating
end

"""
    SpatiallyVaryingCoating{F} <: AbstractCoating

A coating that varies across the surface.
The field `selector` is a function `f(point_local, normal_local) -> AbstractCoating`.
The wrapper automatically transforms global coordinates to the optic's local SDF frame before calling the selector.

# Example
```julia
# AR on left half (x<0), Mirror on right half (x>=0)
c = SpatiallyVaryingCoating((p, n) -> p[1] < 0 ? ar_coating : mirror_coating)
```
"""
struct SpatiallyVaryingCoating{F} <: AbstractCoating
    selector::F
end

function get_coating_at(c::SpatiallyVaryingCoating, optic, point, normal)
    # Transform to local coordinates
    # We assume optic has a `shape` field that is an AbstractSDF.
    # _world_to_sdf is available from SDFs module (which is included in Components.jl)
    p_local = _world_to_sdf(optic.shape, point)

    # We also need to rotate the normal.
    # _world_to_sdf does `T * (p - pos)`. The rotation part is T.
    # So `n_local = T * n_global`.
    sdf_shape = optic.shape
    T = transposed_orientation(sdf_shape)
    n_local = T * normal

    return c.selector(p_local, n_local)
end

"""
    coating_coefficients(coating, n_in, n_out, λ, θ) -> (rs, rp, ts, tp)

Computes the complex reflection and transmission coefficients for s- and p-polarization.

# Arguments
- `coating`: The coating object.
- `n_in`: Refractive index of the incident medium.
- `n_out`: Refractive index of the substrate/exit medium.
- `λ`: Wavelength in [m].
- `θ`: Angle of incidence in [rad].
"""
function coating_coefficients(::Uncoated, n_in, n_out, λ, θ)
    return fresnel_coefficients(θ, n_out / n_in)
end

"""
    SimpleCoating{T, R} <: AbstractCoating

An idealized coating defined by explicit transmission and reflection functions or values.
Useful for modeling generic AR or HR coatings without physical layer details.
"""
struct SimpleCoating{T, R} <: AbstractCoating
    transmission::T # Function (λ, θ) -> T or constant
    reflection::R   # Function (λ, θ) -> R or constant
end

function SimpleCoating(R::Real)
    T = 1.0 - R
    return SimpleCoating((λ, θ) -> T, (λ, θ) -> R)
end

function coating_coefficients(c::SimpleCoating, n_in, n_out, λ, θ)
    # Get power coefficients (Intensity)
    R_power = c.reflection isa Function ? c.reflection(λ, θ) : c.reflection
    T_power = c.transmission isa Function ? c.transmission(λ, θ) : c.transmission

    # We need to return amplitude coefficients (rs, rp, ts, tp).
    # Since we only have power/intensity info, we assume no phase shift (0 phase)
    # and equal split for s/p for simplicity, unless we want to make this more complex later.
    # Note: R = |r|^2, T = |t|^2 * (n_out * cos(theta_t)) / (n_in * cos(theta_i))

    # Fundamental check: T + R <= 1 (allowing for absorption)
    if R_power + T_power > 1.0 + 1e-9
        @warn "SimpleCoating: T + R > 1 ($T_power + $R_power)"
    end

    # Amplitude magnitudes
    # For now, simplistic approximation: r = sqrt(R), t_adjusted...
    # Actually, often users just want "50/50 splitter" or "HR mirror".
    # Let's keep phase 0 for reflection and transmission for now.

    r_mag = sqrt(R_power)

    # Transmission amplitude t is trickier due to impedance mismatch factors in T equation.
    # T = (n2 cos t) / (n1 cos i) * |t|^2
    # But usually "Transmission T" implies the power transmission.
    # OpticalComponents.jl often expects Fresnel coeffs which are AMPlITUDE coeffs.
    # However, for energy tracking `interact3d`, we usually care about the power mostly for splitting?
    # Actually `PolarizedRay` tracks E-field, so we MUST return Amplitude coefficients.

    # If we assume this SimpleCoating is on a wrapper interface where we "overwrite" the physics:
    # We can just return effective r and t that would yield this R and T.

    # Calculate angle of refraction for impedance factor
    sin_theta_t = (n_in / n_out) * sin(θ)
    if abs(sin_theta_t) > 1
        # TIR condition, if this were a real interface.
        # But SimpleCoating might force transmission?
        # Let's fallback to standard physics if inconsistent.
        cos_theta_t = 0.0 # pure imaginary in reality
    else
        cos_theta_t = sqrt(1 - sin_theta_t^2)
    end

    # Impedance factor ratio for T -> t conversion
    # T_power = (n_out * cos_theta_t) / (n_in * cos(θ)) * |t|^2
    # -> |t| = sqrt( T_power * (n_in * cos(θ)) / (n_out * cos_theta_t) )

    if isapprox(cos_theta_t, 0, atol = 1e-9) || isapprox(cos(θ), 0, atol = 1e-9)
        # Avoid div by zero. In TIR or grazing incidence.
        # If TIR, T should be 0.
        t_mag = 0.0
    else
        impedance_ratio = (n_in * cos(θ)) / (n_out * cos_theta_t)
        t_mag = sqrt(T_power * impedance_ratio)
    end

    # Restore "Mirror-like" phase behavior (Diag[-1, 1]) for high reflection
    # This matches expectations for metallic mirrors.
    rs = -r_mag + 0im
    rp = r_mag + 0im

    # Transmission
    ts = t_mag + 0im
    tp = t_mag + 0im

    return (rs, rp, ts, tp)
end

"""
    ThinFilmLayer{M, T}

A single layer in a multilayer stack.
"""
struct ThinFilmLayer{M, T}
    material::M     # RefractiveIndex (function or number)
    thickness::T    # Physical thickness in [m]
end

"""
    MultilayerCoating{L} <: AbstractCoating

A physical coating consisting of a stack of thin film layers.
Calculated using the Transfer Matrix Method (TMM).
The layers are ordered from the incident medium towards the substrate.

# Example
```julia
# Quarter-wave layer of MgF2 (n=1.38) at 550nm
d = 550e-9 / (4 * 1.38)
layer = ThinFilmLayer(1.38, d)
coating = MultilayerCoating([layer])
```
"""
struct MultilayerCoating{L <: AbstractVector{<:ThinFilmLayer}} <: AbstractCoating
    layers::L
    design_λ::Float64 # Optional reference wavelength, could be used for scaling
end

function MultilayerCoating(layers::Vector{<:ThinFilmLayer})
    return MultilayerCoating(layers, 550e-9) # Default visual
end

"""
    _layer_matrix(layer, n, λ, θ, pol)

Calculates the characteristic matrix for a single thin film layer.
"""
function _layer_matrix(layer::ThinFilmLayer, n, λ, θ, pol::Symbol)
    # Phase thickness
    δ = 2π / λ * n * layer.thickness * cos(θ)

    # Admittance
    η = if pol === :s
        n * cos(θ)
    else # p-polarization
        n / cos(θ)
    end

    # Characteristic matrix
    return @SMatrix [cos(δ) (im / η)*sin(δ);
                     im*η*sin(δ) cos(δ)]
end

function coating_coefficients(c::MultilayerCoating, n_in, n_out, λ, θ)
    # Transfer Matrix Method (Abeles Matrix formulation)
    # We follow the convention where the matrix relates E and H fields.

    # beta is the transverse propagation constant (invariant across layers)
    # beta = n * sin(theta)
    beta = n_in * sin(θ)

    # Initialize system matrix (Identity)
    # M = [A B; C D]
    M_s = @SMatrix [1.0+0.0im 0.0+0.0im; 0.0+0.0im 1.0+0.0im]
    M_p = @SMatrix [1.0+0.0im 0.0+0.0im; 0.0+0.0im 1.0+0.0im]

    # Vacuum impedance (normalized to 1 for optical admittance space usually,
    # but let's stick to standard admittance eta = n (TE) or n (TM) variations).
    # Common TMM notation uses optical admittance Y.
    # For s-pol (TE): Y = n * cos(theta)
    # For p-pol (TM): Y = n / cos(theta)

    for layer in c.layers
        n = layer.material isa Function ? layer.material(λ) : layer.material

        # Calculate angle in layer using Snell's invariant beta
        # n sin(theta) = beta
        # cos(theta) = sqrt(1 - sin^2(theta)) = sqrt(1 - (beta/n)^2)
        # Note: can be complex if TIR occurs (evanescent wave)

        cos_theta = sqrt(Complex(1 - (beta / n)^2))

        # Phase thickness delta = (2 pi / lambda) * n * d * cos(theta)
        delta = (2 * π / λ) * n * layer.thickness * cos_theta

        # Optical admittances for this layer
        eta_s = n * cos_theta
        eta_p = n / cos_theta

        # Layer characteristic matrix
        # m = [cos(delta)       (i/eta)*sin(delta)]
        #     [i*eta*sin(delta)  cos(delta)       ]

        i_sin_delta = im * sin(delta)
        cos_delta = cos(delta)

        # Update s-pol matrix
        m_s = @SMatrix [cos_delta (1 / eta_s)*i_sin_delta;
                        eta_s*i_sin_delta cos_delta]

        # Update p-pol matrix
        m_p = @SMatrix [cos_delta (1 / eta_p)*i_sin_delta;
                        eta_p*i_sin_delta cos_delta]

        M_s = M_s * m_s
        M_p = M_p * m_p
    end

    # Calculate final reflection and transmission coefficients
    # System relates input (in) to output (out) fields.
    # We typically match admittance at boundaries.
    # Y0 (input), Ysub (substrate/output)

    # Incident medium angle
    cos_theta_in = cos(θ)
    eta_in_s = n_in * cos_theta_in
    eta_in_p = n_in / cos_theta_in

    # Substrate medium angle
    # n_out sin(theta_out) = beta
    cos_theta_out = sqrt(Complex(1 - (beta / n_out)^2))
    eta_out_s = n_out * cos_theta_out
    eta_out_p = n_out / cos_theta_out

    # Coefficient calculation from the system matrix M = [A B; C D]
    # r = (Y0*M11 + Y0*Ysub*M12 - M21 - Ysub*M22) / (Y0*M11 + Y0*Ysub*M12 + M21 + Ysub*M22)
    # t = 2*Y0 / (Y0*M11 + Y0*Ysub*M12 + M21 + Ysub*M22)

    # S-polarization
    A, B = M_s[1, 1], M_s[1, 2]
    C, D = M_s[2, 1], M_s[2, 2]

    denom_s = eta_in_s * A + eta_in_s * eta_out_s * B + C + eta_out_s * D
    rs = (eta_in_s * A + eta_in_s * eta_out_s * B - C - eta_out_s * D) / denom_s
    ts = 2 * eta_in_s / denom_s

    # P-polarization
    A, B = M_p[1, 1], M_p[1, 2]
    C, D = M_p[2, 1], M_p[2, 2]

    denom_p = eta_in_p * A + eta_in_p * eta_out_p * B + C + eta_out_p * D
    rp = (eta_in_p * A + eta_in_p * eta_out_p * B - C - eta_out_p * D) / denom_p # Note: Sign convention check might be needed here to match Fresnel
    tp = 2 * eta_in_p / denom_p

    # Sign convention adjustment to match BeamletOptics Fresnel:
    # Fresnel usually r_p = (n2 cos1 - n1 cos2) / (n2 cos1 + n1 cos2) vs (n1 cos2 - n2 cos1)
    # Our TMM derivation follows Macleod or standard EM texts.
    # Often p-pol reflection definition varies by -1 factor.
    # Let's trust this standard derivation for now. Users dealing with multi-layers usually work with Power (R, T) anyway where sign doesn't matter.
    # But for interferometry it matters.
    # For Uncoated interface (no layers), M=Identity.
    # r = (Y0 - Ysub)/(Y0 + Ysub).
    # s-pol: (n1 cos1 - n2 cos2) / (...) -> Matches standard Fresnel
    # p-pol: (n1/c1 - n2/c2) / (...) = (n1 c2 - n2 c1) / (...) -> This is -rs (if normal incidence).
    # Standard Fresenl rp = (n2 c1 - n1 c2) / ...
    # So our TMM p-pol r is likely -1 * StandardFresnel p-pol.
    # We should probably flip rp to match the codebase convention if we want `Uncoated` layers to match `Uncoated` type.

    return rs, -rp, ts, tp
end
