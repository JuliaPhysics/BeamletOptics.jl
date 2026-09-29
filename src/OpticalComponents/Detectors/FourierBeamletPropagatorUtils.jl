# ==============================================================================
# Zero-Allocation Trapezoidal Numerical Quadrature
# ==============================================================================

"""
    _trapz_integrate(coords::Tuple{AbstractVector}, I::AbstractVector{T}) where {T}
    _trapz_integrate(coords::Tuple{AbstractVector, AbstractVector}, I::AbstractMatrix{T}) where {T}
    _trapz_integrate(coords::Tuple{AbstractVector, AbstractVector, AbstractVector}, I::AbstractArray{T, 3}) where {T}
    _trapz_integrate(coords::Tuple{AbstractVector, AbstractVector, AbstractVector, AbstractVector}, I::AbstractArray{T, 4}) where {T}
    _trapz_integrate(coords::AbstractVector, I::AbstractArray)

Evaluates the numerical integral of array `I` over the coordinate grid `coords`
(1D, 2D, 3D, or 4D).

# Arguments
- `coords`: Tuple of coordinate vectors along each dimension.
- `I`: Array of function/intensity values matching grid dimensions.

# Returns
- Scalar integrated value of type `T` (or corresponding float type).
"""
@inline function _trapz_integrate(coords::Tuple{AbstractVector}, I::AbstractVector{T}) where {T}
    x = coords[1]
    Nx = length(x)
    Nx == length(I) || throw(DimensionMismatch("Grid and intensity dimension mismatch: length(x) = $Nx vs length(I) = $(length(I))"))
    Nx < 2 && return zero(T)
    s = zero(T)
    @inbounds for i in 1:(Nx - 1)
        s += (x[i + 1] - x[i]) * (I[i] + I[i + 1])
    end
    return 0.5 * s
end

@inline function _trapz_integrate(coords::Tuple{AbstractVector, AbstractVector}, I::AbstractMatrix{T}) where {T}
    x, y = coords
    Nx = length(x)
    Ny = length(y)
    (Nx, Ny) == size(I) || throw(DimensionMismatch("Grid and intensity dimension mismatch: ($Nx, $Ny) vs $(size(I))"))
    (Nx < 2 || Ny < 2) && return zero(T)
    s = zero(T)
    @inbounds for j in 1:(Ny - 1)
        dy = y[j + 1] - y[j]
        s_inner = zero(T)
        for i in 1:(Nx - 1)
            s_inner += (x[i + 1] - x[i]) * (I[i, j] + I[i + 1, j] + I[i, j + 1] + I[i + 1, j + 1])
        end
        s += dy * s_inner
    end
    return 0.25 * s
end

@inline function _trapz_integrate(coords::Tuple{AbstractVector, AbstractVector, AbstractVector}, I::AbstractArray{T, 3}) where {T}
    x, y, z = coords
    Nx = length(x)
    Ny = length(y)
    Nz = length(z)
    (Nx, Ny, Nz) == size(I) || throw(DimensionMismatch("Grid and intensity dimension mismatch: ($Nx, $Ny, $Nz) vs $(size(I))"))
    (Nx < 2 || Ny < 2 || Nz < 2) && return zero(T)
    s = zero(T)
    @inbounds for k in 1:(Nz - 1)
        dz = z[k + 1] - z[k]
        s_k = zero(T)
        for j in 1:(Ny - 1)
            dy = y[j + 1] - y[j]
            s_j = zero(T)
            for i in 1:(Nx - 1)
                c_sum = (I[i, j, k] + I[i + 1, j, k] + I[i, j + 1, k] + I[i + 1, j + 1, k] +
                         I[i, j, k + 1] + I[i + 1, j, k + 1] + I[i, j + 1, k + 1] + I[i + 1, j + 1, k + 1])
                s_j += (x[i + 1] - x[i]) * c_sum
            end
            s_k += dy * s_j
        end
        s += dz * s_k
    end
    return 0.125 * s
end

@inline function _trapz_integrate(coords::Tuple{AbstractVector, AbstractVector, AbstractVector, AbstractVector}, I::AbstractArray{T, 4}) where {T}
    x, y, z, t = coords
    Nx = length(x)
    Ny = length(y)
    Nz = length(z)
    Nt = length(t)
    (Nx, Ny, Nz, Nt) == size(I) || throw(DimensionMismatch("Grid and intensity dimension mismatch: ($Nx, $Ny, $Nz, $Nt) vs $(size(I))"))
    (Nx < 2 || Ny < 2 || Nz < 2 || Nt < 2) && return zero(T)
    s = zero(T)
    @inbounds for l in 1:(Nt - 1)
        dt = t[l + 1] - t[l]
        s_l = zero(T)
        for k in 1:(Nz - 1)
            dz = z[k + 1] - z[k]
            s_k = zero(T)
            for j in 1:(Ny - 1)
                dy = y[j + 1] - y[j]
                s_j = zero(T)
                for i in 1:(Nx - 1)
                    c_sum = (I[i, j, k, l]     + I[i + 1, j, k, l]     + I[i, j + 1, k, l]     + I[i + 1, j + 1, k, l] +
                             I[i, j, k + 1, l]   + I[i + 1, j, k + 1, l]   + I[i, j + 1, k + 1, l]   + I[i + 1, j + 1, k + 1, l] +
                             I[i, j, k, l + 1]   + I[i + 1, j, k, l + 1]   + I[i, j + 1, k, l + 1]   + I[i + 1, j + 1, k, l + 1] +
                             I[i, j, k + 1, l + 1] + I[i + 1, j, k + 1, l + 1] + I[i, j + 1, k + 1, l + 1] + I[i + 1, j + 1, k + 1, l + 1])
                    s_j += (x[i + 1] - x[i]) * c_sum
                end
                s_k += dy * s_j
            end
            s_l += dz * s_k
        end
        s += dt * s_l
    end
    return 0.0625 * s
end

@inline function _trapz_integrate(coords::AbstractVector, I::AbstractArray)
    return _trapz_integrate(Tuple(coords), I)
end

# ==============================================================================
# Numerical Validation and Mathematical Primitives
# ==============================================================================

"""
    safe_finite(x::T, default::T=zero(T)) where {T<:Real}

Validates and returns a finite scalar value, falling back to a default if non-finite.

# Arguments
- `x::T`: Scalar input value of type `T`.
- `default::T`: Fallback value of type `T` used if `x` is non-finite. Defaults to `zero(T)`.

# Returns
- Value of type `T`.
"""
@inline function safe_finite(x::T, default::T=zero(T)) where {T<:Real}
    return isfinite(x) ? x : default
end

@inline function safe_finite(x::Real, default::Real)
    T = promote_type(typeof(x), typeof(default))
    return safe_finite(T(x), T(default))
end

"""
    safe_pos(x::T, fallback::T) where {T<:Real}

Validates and returns a strictly positive finite scalar value, falling back to a default otherwise.

# Arguments
- `x::T`: Scalar input value of type `T`.
- `fallback::T`: Fallback value of type `T` used if `x` is non-finite or non-positive.

# Returns
- Strictly positive value of type `T`.
"""
@inline function safe_pos(x::T, fallback::T) where {T<:Real}
    return (isfinite(x) && x > zero(T)) ? x : fallback
end

@inline function safe_pos(x::Real, fallback::Real)
    T = promote_type(typeof(x), typeof(fallback))
    return safe_pos(T(x), T(fallback))
end

"""
    safe_inv(x::T, default::T=zero(T)) where {T<:Real}

Safely inverts a scalar value, returning a default value if the input or result is zero or non-finite.

# Arguments
- `x::T`: Scalar input value of type `T`.
- `default::T`: Fallback value of type `T` used if `x` cannot be inverted. Defaults to `zero(T)`.

# Returns
- Reciprocal scalar value of type `T`.
"""
@inline function safe_inv(x::T, default::T=zero(T)) where {T<:Real}
    if isfinite(x) && !iszero(x)
        ix = inv(x)
        return isfinite(ix) ? T(ix) : default
    end
    return default
end

@inline function safe_inv(x::Real, default::Real)
    T = promote_type(typeof(x), typeof(default))
    return safe_inv(T(x), T(default))
end

"""
    _beamlet_defocus_distance(C::T, w::T, wl::T) where {T<:Real}

Calculates the paraxial defocus distance from wavefront curvature, beam waist, and wavelength.

# Arguments
- `C::T`: Wavefront curvature scalar in 1/m of type `T`.
- `w::T`: Beam waist radius in meters of type `T`.
- `wl::T`: Optical wavelength in meters of type `T`.

# Returns
- Paraxial defocus distance in meters of type `T`.
"""
@inline function _beamlet_defocus_distance(C::T, w::T, wl::T)::T where {T<:Real}
    if !isfinite(C) || !isfinite(w) || w <= zero(T) || !isfinite(wl) || wl <= zero(T)
        return zero(T)
    end
    zr_inv = wl / (T(π) * w^2)
    denom = C^2 + zr_inv^2
    res = denom > zero(T) ? (C / denom) : zero(T)
    return isfinite(res) ? res : zero(T)
end

@inline function _beamlet_defocus_distance(C::Real, w::Real, wl::Real)
    T = promote_type(typeof(float(C)), typeof(float(w)), typeof(float(wl)))
    return _beamlet_defocus_distance(T(C), T(w), T(wl))
end

# ==============================================================================
# Layer 4: Observable Extraction & Detector Readout
# ==============================================================================

"""
    _calculate_pixel_area(grid::SpatialGrid{T, D}) where {T, D}

Calculates the spatial integration pixel area or length element from grid coordinate ranges.

Parameters:
- `grid`: Spatial detector grid of type `SpatialGrid{T, D}`.

Returns:
- Integration cell measure (length in 1D, area in 2D or higher) of type `T`.

Contracts & Invariants:
- Uses standard range stepping.
- Returned pixel area is strictly positive.
"""
@inline function _calculate_pixel_area(grid::SpatialGrid{T, D}) where {T, D}
    if D == 1
        return abs(step(grid.ranges[1]))
    else
        return abs(step(grid.ranges[1])) * abs(step(grid.ranges[2]))
    end
end

"""
    DetectorObservables{T <: Real, A <: AbstractArray}

Container holding physical optical observables extracted from the synthesized field.

Parameters:
- `intensity`: Spatial optical intensity distribution array of type A in W/m^2.
- `phase`: Spatial wavefront phase distribution array of type A in radians.
- `stokes`: Four-element tuple of Stokes parameter arrays (S0, S1, S2, S3) of type A in W/m^2.
- `power`: Total integrated optical power across the detector aperture of type T in Watts.

Invariants:
- Optical intensity values satisfy non-negativity (I >= 0).
- Stokes polarization parameters satisfy consistency (S0^2 >= S1^2 + S2^2 + S3^2).
- Total integrated optical power satisfies non-negativity (power >= 0).
- Array dimensions match the underlying spatial field grid dimensions.
"""
mutable struct DetectorObservables{T <: Real, A <: AbstractArray}
    intensity::A
    phase::A
    stokes::NTuple{4, A}
    power::T
end

"""
    DetectorObservables(intensity::AbstractArray{T}, phase::AbstractArray{T}, stokes::NTuple{4, AbstractArray{T}}, power::Real = zero(T)) where {T <: Real}
    DetectorObservables(dims::Dims{D}, ::Type{T} = Float64) where {D, T <: Real}
    DetectorObservables(dims::Int...; T::Type = Float64)

Constructs an observable container holding intensity, phase, Stokes parameters, and integrated optical power.

Parameters:
- `intensity`: Spatial intensity array in W/m^2.
- `phase`: Spatial phase angle array in radians.
- `stokes`: 4-tuple of Stokes polarization parameter arrays (S0, S1, S2, S3) in W/m^2.
- `power`: Total integrated optical power in Watts.
- `dims`: Array dimensions for preallocating uninitialized observable arrays.
- `T`: Floating-point scalar precision type, defaults to `Float64`.

Returns:
- New `DetectorObservables{T, A}` container.
"""
@inline DetectorObservables(intensity::AbstractArray{T}, phase::AbstractArray{T}, stokes::NTuple{4, AbstractArray{T}}, power::Real) where {T <: Real} =
    DetectorObservables{T, typeof(intensity)}(intensity, phase, stokes, T(power))

@inline DetectorObservables(intensity::AbstractArray{T}, phase::AbstractArray{T}, stokes::NTuple{4, AbstractArray{T}}) where {T <: Real} =
    DetectorObservables{T, typeof(intensity)}(intensity, phase, stokes, zero(T))

@inline function DetectorObservables(
    dims::Dims{D},
    ::Type{T_val}
) where {D, T_val <: Real}
    intensity = Array{T_val, D}(undef, dims)
    phase = Array{T_val, D}(undef, dims)
    stokes = (
        Array{T_val, D}(undef, dims),
        Array{T_val, D}(undef, dims),
        Array{T_val, D}(undef, dims),
        Array{T_val, D}(undef, dims)
    )
    return DetectorObservables{T_val, Array{T_val, D}}(intensity, phase, stokes, zero(T_val))
end

@inline DetectorObservables(dims::Dims{D}; T::Type{T_val} = Float64) where {D, T_val <: Real} =
    DetectorObservables(dims, T_val)

@inline DetectorObservables(dims::Int...; T::Type{T_val} = Float64) where {T_val <: Real} =
    DetectorObservables(dims, T_val)

"""
    calculate_observables!(obs::DetectorObservables{T}, d::FourierBeamletPropagator, field::AbstractArray{SVector{3, Complex{T}}, D}; pixel_area::Real = 1.0, Z::Real = Z_vacuum) where {T <: Real, D}

In-place mutation of preallocated observable arrays without heap allocations.

Parameters:
- `obs`: Target observable container of type `DetectorObservables{T}` with preallocated destination arrays.
- `d`: Propagator instance of type `FourierBeamletPropagator`.
- `field`: Complex 3D electric field array of type `AbstractArray{SVector{3, Complex{T}}, D}` in V/m.
- `pixel_area`: Spatial integration cell area of type `Real` in m^2, strictly positive, defaults to 1.0.
- `Z`: Characteristic wave impedance of type `Real` in Ohms, strictly positive, defaults to `Z_vacuum`.

Returns:
- Mutated `DetectorObservables{T}` container with populated observables and total integrated power.

Contracts & Invariants:
- Zero heap allocations during evaluation (@allocated == 0).
- Input observable arrays are mutated directly in-place.
- In-place array dimensions must match the input field dimensions.
- Intensity, Stokes S0 parameter, and integrated power satisfy non-negativity.
- Stokes parameters satisfy the polarization consistency relation.
"""
function calculate_observables!(
    obs::DetectorObservables{T}, d::FourierBeamletPropagator,
    field::AbstractArray{SVector{3, Complex{T}}, D};
    pixel_area::Real = 1.0, Z::Real = Z_vacuum
) where {T <: Real, D}
    if size(obs.intensity) != size(field) ||
       size(obs.phase) != size(field) ||
       size(obs.stokes[1]) != size(field) ||
       size(obs.stokes[2]) != size(field) ||
       size(obs.stokes[3]) != size(field) ||
       size(obs.stokes[4]) != size(field)
        throw(DimensionMismatch("DetectorObservables array dimensions $(size(obs.intensity)) do not match field dimensions $(size(field))."))
    end

    inv_2Z = inv(2 * T(Z))
    pix_T = T(pixel_area)

    sx = zero(T)
    sy = zero(T)
    sz = zero(T)
    @inbounds for i in eachindex(field)
        Ef = field[i]
        sx += abs2(Ef[1])
        sy += abs2(Ef[2])
        sz += abs2(Ef[3])
    end
    dom_idx = (sx >= sy && sx >= sz) ? 1 : (sy >= sz ? 2 : 3)

    ex, ez = _detector_local_axes(d, T)
    en = normalize(cross(ex, ez))

    I_arr = obs.intensity
    phi_arr = obs.phase
    s0_arr = obs.stokes[1]
    s1_arr = obs.stokes[2]
    s2_arr = obs.stokes[3]
    s3_arr = obs.stokes[4]

    P_sum = zero(T)
    @inbounds for i in eachindex(field)
        Ef = field[i]
        Eh = dot(Ef, ex)
        Ev = dot(Ef, ez)
        En = dot(Ef, en)

        Eh_sq = abs2(Eh)
        Ev_sq = abs2(Ev)
        En_sq = abs2(En)

        E_trans_sq = (Eh_sq + Ev_sq) * inv_2Z
        I_val = max(zero(T), E_trans_sq + En_sq * inv_2Z)
        I_arr[i] = I_val
        P_sum += I_val

        phi_arr[i] = angle(Ef[dom_idx])

        s0_val = max(zero(T), E_trans_sq)
        s1_val = (Eh_sq - Ev_sq) * inv_2Z

        Eh_conj_Ev = conj(Eh) * Ev
        s2_val = (2 * real(Eh_conj_Ev)) * inv_2Z
        s3_val = (2 * imag(Eh_conj_Ev)) * inv_2Z

        s_norm_sq = s1_val^2 + s2_val^2 + s3_val^2
        s_norm = sqrt(s_norm_sq)
        s0_val = max(s0_val, s_norm)
        while s0_val^2 < s_norm_sq
            s0_val = nextfloat(s0_val)
        end

        s0_arr[i] = s0_val
        s1_arr[i] = s1_val
        s2_arr[i] = s2_val
        s3_arr[i] = s3_val
    end

    total_power = max(zero(T), P_sum * pix_T)
    obs.power = total_power
    return obs
end

"""
    calculate_observables!(obs::DetectorObservables{T}, d::FourierBeamletPropagator, field::AbstractArray{SVector{3, Complex{T}}, D}, grid::SpatialGrid{T, D}; Z::Real = Z_vacuum) where {T <: Real, D}

In-place mutation of preallocated observable arrays using grid cell spacing to determine integration cell area.

Parameters:
- `obs`: Target observable container of type DetectorObservables{T} with preallocated destination arrays.
- `d`: Propagator instance of type FourierBeamletPropagator.
- `field`: Complex 3D electric field array of type AbstractArray{SVector{3, Complex{T}}, D} in V/m.
- `grid`: Spatial detector grid of type SpatialGrid{T, D} defining coordinate sampling intervals.
- `Z`: Characteristic wave impedance of type Real in Ohms, strictly positive (Z > 0), defaults to Z_vacuum.

Output:
- Mutated DetectorObservables{T} container with populated observables and total integrated power.

Contracts & Invariants:
- Zero heap allocations during evaluation (@allocated == 0).
- Spatial cell integration area is determined by grid coordinate ranges.
- In-place array dimensions must match the input field dimensions.
"""
function calculate_observables!(
    obs::DetectorObservables{T}, d::FourierBeamletPropagator,
    field::AbstractArray{SVector{3, Complex{T}}, D}, grid::SpatialGrid{T, D};
    Z::Real = Z_vacuum
) where {T <: Real, D}
    return calculate_observables!(obs, d, field; pixel_area = _calculate_pixel_area(grid), Z = Z)
end

"""
    calculate_observables(d::FourierBeamletPropagator, field::AbstractArray{SVector{3, Complex{T}}, D}; pixel_area::Real = 1.0, Z::Real = Z_vacuum) where {T <: Real, D}

Extraction of physical optical observables from a synthesized complex vector electric field.

Parameters:
- `d`: Propagator instance of type `FourierBeamletPropagator`.
- `field`: Complex 3D electric field array of type `AbstractArray{SVector{3, Complex{T}}, D}` in V/m.
- `pixel_area`: Spatial integration cell area of type `Real` in m^2, strictly positive, defaults to 1.0.
- `Z`: Characteristic wave impedance of type `Real` in Ohms, strictly positive, defaults to `Z_vacuum`.

Returns:
- Newly allocated `DetectorObservables{T}` container holding intensity, phase, Stokes parameters, and total power.

Contracts & Invariants:
- Dimensions of all output observable arrays match the input field dimensions.
- Energy conservation holds within Parseval limits.
- Intensity, Stokes S0 parameter, and integrated power satisfy non-negativity.
- Stokes parameters satisfy the polarization consistency relation.
"""
function calculate_observables(
    d::FourierBeamletPropagator, field::AbstractArray{SVector{3, Complex{T}}, D};
    pixel_area::Real = 1.0, Z::Real = Z_vacuum
) where {T <: Real, D}
    obs = DetectorObservables(size(field), T)
    return calculate_observables!(obs, d, field; pixel_area = pixel_area, Z = Z)
end

"""
    calculate_observables(d::FourierBeamletPropagator, field::AbstractArray{SVector{3, Complex{T}}, D}, grid::SpatialGrid{T, D}; Z::Real = Z_vacuum) where {T <: Real, D}

Extraction of physical optical observables using detector spatial grid cell spacing to determine integration cell area.

Parameters:
- `d`: Propagator instance of type FourierBeamletPropagator.
- `field`: Complex 3D electric field array of type AbstractArray{SVector{3, Complex{T}}, D} in V/m.
- `grid`: Spatial detector grid of type SpatialGrid{T, D} defining coordinate sampling intervals.
- `Z`: Characteristic wave impedance of type Real in Ohms, strictly positive (Z > 0), defaults to Z_vacuum.

Output:
- Newly allocated DetectorObservables{T} container holding intensity, phase, Stokes parameters, and total power.

Contracts & Invariants:
- Spatial cell integration area is determined by grid coordinate ranges.
- Dimensions of all output observable arrays match the input field dimensions.
- Intensity, Stokes S0 parameter, and integrated power satisfy non-negativity.
"""
function calculate_observables(
    d::FourierBeamletPropagator, field::AbstractArray{SVector{3, Complex{T}}, D},
    grid::SpatialGrid{T, D}; Z::Real = Z_vacuum
) where {T <: Real, D}
    return calculate_observables(d, field; pixel_area = _calculate_pixel_area(grid), Z = Z)
end

"""
    calculate_observables(d::FourierBeamletPropagator, grid::SpatialGrid{T, D}; Z::Real = Z_vacuum) where {T <: Real, D}

Extraction of physical optical observables directly from propagator and spatial grid.

Parameters:
- `d`: Propagator instance of type FourierBeamletPropagator.
- `grid`: Spatial detector grid of type SpatialGrid{T, D} defining coordinate sampling intervals.
- `Z`: Characteristic wave impedance of type Real in Ohms, strictly positive (Z > 0), defaults to Z_vacuum.

Output:
- Newly allocated DetectorObservables{T} container holding intensity, phase, Stokes parameters, and total power.

Contracts & Invariants:
- Synthesizes the electric field on the spatial grid via pure NUFFT.
- Dimensions of all output observable arrays match the input grid dimensions.
- Intensity, Stokes S0 parameter, and integrated power satisfy non-negativity.
"""
function calculate_observables(
    d::FourierBeamletPropagator, grid::SpatialGrid{T, D};
    Z::Real = Z_vacuum
) where {T <: Real, D}
    field = synthesize_field(d, grid)
    return calculate_observables(d, field, grid; Z = Z)
end

"""
    calculate_observables(d::FourierBeamletPropagator, ranges::NTuple{D, AbstractRange}; Z::Real = Z_vacuum) where {D}

Extraction of physical optical observables directly from propagator and coordinate range tuples.

Parameters:
- `d`: Propagator instance of type FourierBeamletPropagator.
- `ranges`: Coordinate range tuples of type NTuple{D, AbstractRange} defining coordinate sampling intervals.
- `Z`: Characteristic wave impedance of type Real in Ohms, strictly positive (Z > 0), defaults to Z_vacuum.

Output:
- Newly allocated DetectorObservables container holding intensity, phase, Stokes parameters, and total power.

Contracts & Invariants:
- Synthesizes the electric field on the spatial grid via pure NUFFT.
- Dimensions of all output observable arrays match the range lengths.
- Intensity, Stokes S0 parameter, and integrated power satisfy non-negativity.
"""
function calculate_observables(
    d::FourierBeamletPropagator, ranges::NTuple{D, AbstractRange};
    Z::Real = Z_vacuum
) where {D}
    return calculate_observables(d, SpatialGrid(ranges); Z = Z)
end

"""
    calculate_observables!(obs::DetectorObservables{T}, d::FourierBeamletPropagator, grid::SpatialGrid{T, D}; Z::Real = Z_vacuum) where {T <: Real, D}

In-place evaluation of physical optical observables directly into preallocated arrays.

Parameters:
- `obs`: Target observable container of type DetectorObservables{T} with preallocated destination arrays.
- `d`: Propagator instance of type FourierBeamletPropagator.
- `grid`: Spatial detector grid of type SpatialGrid{T, D} defining coordinate sampling intervals.
- `Z`: Characteristic wave impedance of type Real in Ohms, strictly positive (Z > 0), defaults to Z_vacuum.

Output:
- Mutated DetectorObservables{T} container with populated observables and total integrated power.

Contracts & Invariants:
- Input observable arrays are mutated directly in-place.
- In-place array dimensions must match the input grid dimensions.
- Intensity, Stokes S0 parameter, and integrated power satisfy non-negativity.
"""
function calculate_observables!(
    obs::DetectorObservables{T}, d::FourierBeamletPropagator,
    grid::SpatialGrid{T, D}; Z::Real = Z_vacuum
) where {T <: Real, D}
    field = synthesize_field(d, grid)
    return calculate_observables!(obs, d, field, grid; Z = Z)
end

"""
    calculate_observables!(obs::DetectorObservables{T}, d::FourierBeamletPropagator, ranges::NTuple{D, AbstractRange}; Z::Real = Z_vacuum) where {T <: Real, D}

In-place evaluation of physical optical observables directly into preallocated arrays from coordinate range tuples.

Parameters:
- `obs`: Target observable container of type DetectorObservables{T} with preallocated destination arrays.
- `d`: Propagator instance of type FourierBeamletPropagator.
- `ranges`: Coordinate range tuples of type NTuple{D, AbstractRange} defining coordinate sampling intervals.
- `Z`: Characteristic wave impedance of type Real in Ohms, strictly positive (Z > 0), defaults to Z_vacuum.

Output:
- Mutated DetectorObservables{T} container with populated observables and total integrated power.

Contracts & Invariants:
- Input observable arrays are mutated directly in-place.
- In-place array dimensions must match the range lengths.
- Intensity, Stokes S0 parameter, and integrated power satisfy non-negativity.
"""
function calculate_observables!(
    obs::DetectorObservables{T}, d::FourierBeamletPropagator,
    ranges::NTuple{D, AbstractRange}; Z::Real = Z_vacuum
) where {T <: Real, D}
    return calculate_observables!(obs, d, SpatialGrid{T, D}(ranges); Z = Z)
end

# Default fallback overloads when propagator is omitted
"""
    calculate_observables(field::AbstractArray{SVector{3, Complex{T}}, D}; pixel_area::Real = 1.0, Z::Real = Z_vacuum) where {T <: Real, D}
    calculate_observables(field::AbstractArray{SVector{3, Complex{T}}, D}, grid::SpatialGrid{T, D}; Z::Real = Z_vacuum) where {T <: Real, D}
    calculate_observables!(obs::DetectorObservables{T}, field::AbstractArray{SVector{3, Complex{T}}, D}; pixel_area::Real = 1.0, Z::Real = Z_vacuum) where {T <: Real, D}
    calculate_observables!(obs::DetectorObservables{T}, field::AbstractArray{SVector{3, Complex{T}}, D}, grid::SpatialGrid{T, D}; Z::Real = Z_vacuum) where {T <: Real, D}

Evaluates optical observables directly from a synthesized electric field array when propagator context is omitted.

Parameters:
- `obs`: Target observable container of type `DetectorObservables{T}`.
- `field`: Complex 3D electric field array in V/m.
- `grid`: Spatial detector grid of type `SpatialGrid{T, D}`.
- `pixel_area`: Spatial integration cell area in m^2, defaults to 1.0.
- `Z`: Characteristic impedance in Ohms, defaults to vacuum impedance `Z_vacuum`.

Returns:
- `DetectorObservables{T}` containing computed physical observables.
"""
calculate_observables(field::AbstractArray{SVector{3, Complex{T}}, D}; pixel_area::Real = 1.0, Z::Real = Z_vacuum) where {T <: Real, D} =
    calculate_observables(FourierBeamletPropagator(one(T)), field; pixel_area = pixel_area, Z = Z)

calculate_observables(field::AbstractArray{SVector{3, Complex{T}}, D}, grid::SpatialGrid{T, D}; Z::Real = Z_vacuum) where {T <: Real, D} =
    calculate_observables(FourierBeamletPropagator(one(T)), field, grid; Z = Z)

calculate_observables!(obs::DetectorObservables{T}, field::AbstractArray{SVector{3, Complex{T}}, D}; pixel_area::Real = 1.0, Z::Real = Z_vacuum) where {T <: Real, D} =
    calculate_observables!(obs, FourierBeamletPropagator(one(T)), field; pixel_area = pixel_area, Z = Z)

calculate_observables!(obs::DetectorObservables{T}, field::AbstractArray{SVector{3, Complex{T}}, D}, grid::SpatialGrid{T, D}; Z::Real = Z_vacuum) where {T <: Real, D} =
    calculate_observables!(obs, FourierBeamletPropagator(one(T)), field, grid; Z = Z)

export DetectorObservables, calculate_observables, calculate_observables!

# ==============================================================================
# Fourier Beamlet Propagator Hit Accessors and Primitives
# ==============================================================================

@inline function optical_path_length(hit::GaussianBeamletHit{T}) where {T}
    chief_beam = hit.gauss.chief
    p_parent = chief_beam.parent
    opl = if !isnothing(p_parent)
        optical_path_length(p_parent)::T
    else
        p_gauss = hit.gauss.parent
        if !isnothing(p_gauss)
            optical_path_length(p_gauss)::T
        else
            zero(T)
        end
    end
    chief_rays = rays(chief_beam)
    target_id = clamp(hit.id, 0, length(chief_rays))
    for i in 1:target_id
        ray = chief_rays[i]
        isnothing(intersection(ray)) && break
        opl += optical_path_length(ray)::T
    end
    return opl
end
optical_path_length(hit::AstigmaticGaussianBeamletHit) = hit.l0 + optical_path_length(hit.agb.c.rays[hit.id]) + hit.Δl

wavelength(hit::AbstractRayHit) = wavelength(hit.ray)
wavelength(hit::GaussianBeamletHit) = wavelength(hit.gauss)
wavelength(hit::AstigmaticGaussianBeamletHit) = wavelength(hit.agb)
@inline wavelength(hit::AbstractDetectorHit) = wavenumber(hit) > 0 ? (2π / wavenumber(hit)) : 1.0

polarization(hit::RayHit{T}) where {T} = polarization(hit.ray)
function polarization(hit::GaussianBeamletHit{T}) where {T}
    chief_ray = hit.gauss.chief.rays[hit.id]
    if chief_ray isa PolarizedRay
        return polarization(chief_ray)
    else
        dir = direction(hit)
        n_perp = transverse_polarization(dir)
        amp = hit.gauss.E0
        return Point3{Complex{T}}(Complex{T}(n_perp[1] * amp), Complex{T}(n_perp[2] * amp), Complex{T}(n_perp[3] * amp))
    end
end
function polarization(hit::AstigmaticGaussianBeamletHit{T}) where {T}
    chief_ray = hit.agb.c.rays[hit.id]
    if chief_ray isa PolarizedRay
        return polarization(chief_ray)
    else
        dir = direction(hit)
        n_perp = transverse_polarization(dir)
        amp = hit.E_ref_amp
        return Point3{Complex{T}}(Complex{T}(n_perp[1] * amp), Complex{T}(n_perp[2] * amp), Complex{T}(n_perp[3] * amp))
    end
end

@inline function transverse_polarization(dir::AbstractVector{T}) where T
    v = normalize(SVector{3,T}(dir[1], dir[2], dir[3]))
    if v[3] > -T(0.99999)
        inv_1_pz = one(T) / (one(T) + v[3])
        ux = one(T) - v[1]^2 * inv_1_pz
        uy = -v[1] * v[2] * inv_1_pz
        uz = -v[1]
        return Point3{T}(ux, uy, uz)
    else
        inv_1_mz = one(T) / (one(T) - v[3])
        ux = one(T) - v[1]^2 * inv_1_mz
        uy = -v[1] * v[2] * inv_1_mz
        uz = v[1]
        return Point3{T}(ux, uy, uz)
    end
end
@inline transverse_polarization(dir::Point3{T}) where T = transverse_polarization(SVector{3,T}(dir[1], dir[2], dir[3]))

