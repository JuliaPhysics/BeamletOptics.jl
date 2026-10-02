"""
    GaussianBeamletDecomposition(pos, dir, λ, w0; n_grid=20)

Decomposes a macroscopic Gaussian beam of waist `w0` into an `n_grid × n_grid`
array of microscopic `AstigmaticGaussianBeamlet`s. This allows the accurate
tracing of large Gaussian beams through highly aberrative or non-paraxial
optical systems, perfectly preserving higher-order phase aberrations via the
coherent superposition of the sub-beamlets.

# Arguments
- `pos`: Waist center position of the macroscopic beam.
- `dir`: Propagation direction.
- `λ`: Wavelength.
- `w0`: Macroscopic beam waist.
- `n_grid`: Number of beamlets along one axis (default `20` yields `400` total beamlets).
- `overlap`: Scaling factor for the sub-waist relative to grid spacing (default `1.2` ensures smooth overlap).
- `basis`: Optional tuple `(ex, ey)` of the macroscopic sampling grid axes, i.e. the 3D directions corresponding to the `x` and `y` grid axes. `ex` must not be zero or parallel to `dir`; its component normal to `dir` becomes the local x-axis of the group [`orientation`](@ref).
- `randomize_axes`: If `true`, the internal principal axes of each individual beamlet are randomly rotated. This averages out numerical grid-alignment biases and is essential for preserving rotational symmetry in focused spots.
- `rng`: Random number generator to use for `randomize_axes`.
- `P0`: Total power of the macroscopic Gaussian beam in [W].
- `E0`: Optional Jones vector defining the polarization and initial phase of the beam. If `nothing`, defaults to linear polarization along the first grid axis.
- `threshold`: Relative amplitude below which beamlets are not spawned (default `1e-4`).
"""
function GaussianBeamletDecomposition(
        pos::AbstractArray{P},
        dir::AbstractArray{D1},
        λ::L,
        w0::W;
        n_grid::Int = 20,
        overlap::Float64 = 1.2,
        basis::Union{Nothing, Tuple{AbstractVector, AbstractVector}} = nothing,
        randomize_axes::Bool = false,
        rng = Random.GLOBAL_RNG,
        P0::Real = get_default_power(),
        E0 = nothing,
        threshold::Float64 = 1e-4
) where {P <: Real, D1 <: Real, W <: Real, L <: Real}
    T = promote_type(P, D1, W, L)
    dir_n = normalize(dir)
    # Use provided basis or fallback to automatic orthogonal basis
    e1 = isnothing(basis) ? normal3d(dir_n) : normalize(basis[1])
    e2 = isnothing(basis) ? normalize(cross(dir_n, e1)) : normalize(basis[2])

    # Determine macroscopic peak amplitude
    if isnothing(E0)
        I0 = (2 * P0) / (π * w0^2)
        E_phasor = T(electric_field(I0))
        pol = e1
    else
        E_phasor = one(T)
        pol = E0
    end

    # We tile over a 4*w0 x 4*w0 area to capture the macroscopic Gaussian tails
    D = 4 * w0
    Δ = D / n_grid

    # The sub-waist should be slightly larger than grid spacing to ensure smooth overlap
    w0s = Δ * overlap

    beams = Vector{AstigmaticGaussianBeamlet{T}}()
    xs = LinRange(-D / 2 + Δ / 2, D / 2 - Δ / 2, n_grid)
    ys = LinRange(-D / 2 + Δ / 2, D / 2 - Δ / 2, n_grid)

    for x in xs
        for y in ys
            offset = x * e1 + y * e2
            # Macroscopic Gaussian envelope amplitude weighting
            amp = exp(-(x^2 + y^2) / w0^2)
            # Threshold to avoid tracing zero-amplitude beamlets
            if amp > threshold
                # Normalization factor for coherent superposition
                norm_factor = (Δ^2) / (π * w0s^2)
                # Standard linear polarization along e1
                # Support vector for the beamlet axes
                if randomize_axes
                    base_s = normal3d(dir_n)
                    ortho_s = cross(dir_n, base_s)
                    phi = 2π * rand(rng)
                    local_support = base_s * cos(phi) + ortho_s * sin(phi)
                else
                    local_support = nothing
                end
                b = AstigmaticGaussianBeamlet(
                    pos + offset, dir_n, λ, w0s; E0 = (amp * norm_factor * E_phasor) .* pol, support = local_support)
                push!(beams, b)
            end
        end
    end
    # group orientation from the first grid axis, orthogonalized w.r.t. dir (sampling unchanged)
    e1_o = isnothing(basis) ? e1 : sampling_basis(dir_n, basis[1], T)
    return AstigmaticBeamGroup(beams, pos, _group_orientation(dir_n, e1_o, T))
end

@inline _wrap_phase(Δ::Real) = mod2pi(Δ + π) - π

"""
    _eikonal_direction(Δφ, i, j, nx, ny, dx, dy, k, e1, e2, dir_n, ∇φ0 = (0, 0))

Local propagation direction at grid index `(i, j)` from the Eikonal equation
`∇φ = k sin θ`. `Δφ(a, b)` returns the phase difference from grid index `a` to the
neighboring index `b` (tuples); central differences, one-sided at the grid edges, give the
transverse direction cosines `sin θx`, `sin θy` along `e1`, `e2`. `∇φ0` is a phase
gradient along `e1`, `e2` that is known analytically and added to the measured one (e.g.
that of a reference sphere, see [`_sphere_phase`](@ref)). Returns `nothing` if the
gradient is too steep to represent a non-grazing propagation direction
(`sin²θx + sin²θy ≥ 1`) or is `NaN`. Shared by both [`WavefrontBeamletDecomposition`](@ref) methods, which differ only
in how they measure `Δφ` (scalar phase or vector overlap).
"""
function _eikonal_direction(Δφ, i, j, nx, ny, dx, dy, k, e1, e2, dir_n, ∇φ0 = (0, 0))
    ip, im_ = min(i + 1, nx), max(i - 1, 1)
    jp, jm = min(j + 1, ny), max(j - 1, 1)
    dφ_dx = (Δφ((i, j), (ip, j)) + Δφ((im_, j), (i, j))) / ((ip - im_) * dx) + ∇φ0[1]
    dφ_dy = (Δφ((i, j), (i, jp)) + Δφ((i, jm), (i, j))) / ((jp - jm) * dy) + ∇φ0[2]

    sin_θx = dφ_dx / k
    sin_θy = dφ_dy / k
    # grazing directions (sin²θ = 1) are rejected too: they have no component along
    # `dir_n` to carry the field's normal component (see `_transverse_field`)
    if isnan(sin_θx) || isnan(sin_θy) || (sin_θx^2 + sin_θy^2 >= 1.0)
        return nothing
    end
    cos_θz = sqrt(1.0 - sin_θx^2 - sin_θy^2)
    return normalize(sin_θx * e1 + sin_θy * e2 + cos_θz * dir_n)
end

# Largest non-NaN value, so that masked (NaN) samples do not hide the rest of the field.
_finite_maximum(A) = maximum(a for a in A if !isnan(a); init = zero(eltype(A)))

"""
    _transverse_field(E1, E2, e1, e2, n, d)

3D field vector with the tangential components `E1`, `E2` along `e1`, `e2` of a plane with
normal `n`, completed by the normal component that makes it transverse to the propagation
direction `d` (`E·d = 0`). Keeps the given tangential field of an oblique beam instead of
projecting it, which would scale it by `cos²θ`.
"""
function _transverse_field(E1, E2, e1, e2, n, d)
    En = -(dot(d, e1) * E1 + dot(d, e2) * E2) / dot(d, n)
    return E1 * e1 + E2 * e2 + En * n
end

"""
    _sphere_phase(x, y, R, k)
    _sphere_phase_gradient(x, y, R, k)

Phase `sign(R) k (√(ρ² + R²) − |R|)` of a reference sphere with radius `R` about the axis
through `x = y = 0` (`ρ² = x² + y²`), written without cancellation for `ρ ≪ |R|`, and its
gradient along `x`, `y`. Both are zero for `R = Inf`.
"""
_sphere_phase(x, y, R, k) = sign(R) * k * (x^2 + y^2) / (sqrt(x^2 + y^2 + R^2) + abs(R))
_sphere_phase_gradient(x, y, R, k) = (sign(R) * k / sqrt(x^2 + y^2 + R^2)) .* (x, y)

function _check_reference_radius(R)
    !isnan(R) && !iszero(R) ||
        throw(ArgumentError("R must be nonzero; use R = Inf for no reference sphere"))
    return nothing
end

# Waist radius and waist position along the beam, relative to a plane on which a Gaussian
# beam has the radius `w` and the wavefront curvature `c = 1/R`. From the complex beam
# parameter 1/q = c − i λ/(π w²): the plane lies Re(q) behind the waist.
function _waist(w, c, λ)
    q = 1 / complex(c, -λ / (π * w^2))
    return sqrt(λ * imag(q) / π), -real(q)
end

"""
    WavefrontBeamletDecomposition(x, y, amplitude, phase, dir, λ; threshold=1e-4)

Decomposes an arbitrary complex scalar field (defined by a spatial `amplitude` and `phase`
distribution on a 2D grid `x` and `y`) into a collection of `AstigmaticGaussianBeamlet`s.

This function uses the Eikonal approximation to map the local phase gradient into a
local propagation direction for each beamlet. It allows users to import arbitrary,
aberrated, or custom beam profiles (e.g. from a camera or a wavefront sensor) and
propagate them through a `BeamletOptics` system.

The grid is placed in the plane normal to `dir` through the global origin, i.e. the group
`center` is `(0, 0, 0)`. Use [`translate3d!`](@ref) to move the decomposed field to its
actual position.

# Arguments
- `x`, `y`: Vectors defining the 1D spatial coordinates of the 2D field grid.
- `amplitude`: 2D array of field amplitudes `(length(x) × length(y))`.
- `phase`: 2D array of phase values in radians.
- `dir`: The macroscopic reference propagation direction (e.g., `[0, 1, 0]`).
- `λ`: Wavelength.
- `threshold`: Relative amplitude threshold below which beamlets are not spawned (saves computation).
- `overlap`: Scaling factor for the sub-waist relative to grid spacing (default `1.2` ensures smooth overlap).
- `basis`: Optional tuple `(ex, ey)` of the macroscopic sampling grid axes, i.e. the 3D directions corresponding to the `x` and `y` input axes. `ex` must not be zero or parallel to `dir`; its component normal to `dir` becomes the local x-axis of the group [`orientation`](@ref).
- `randomize_axes`: If `true`, the internal principal axes of each individual beamlet are randomly rotated. This averages out numerical grid-alignment biases and helps preserve rotational symmetry in focused patterns.
- `rng`: Random number generator to use for `randomize_axes`.
- `E0`: Optional reference polarization vector or Jones vector. If `nothing`, defaults to linear polarization along the first grid axis.
"""
function WavefrontBeamletDecomposition(
        x::AbstractVector{P1},
        y::AbstractVector{P2},
        amplitude::AbstractMatrix{A},
        phase::AbstractMatrix{Ph},
        dir::AbstractArray{D},
        λ::L;
        threshold::Float64 = 1e-4,
        overlap::Float64 = 1.2,
        basis::Union{Nothing, Tuple{AbstractVector, AbstractVector}} = nothing,
        randomize_axes::Bool = false,
        rng = Random.GLOBAL_RNG,
        E0 = nothing
) where {P1 <: Real, P2 <: Real, A <: Real, Ph <: Real, D <: Real, L <: Real}
    T = promote_type(P1, P2, A, Ph, D, L)

    nx, ny = length(x), length(y)
    if size(amplitude) != (nx, ny) || size(phase) != (nx, ny)
        throw(DimensionMismatch("Amplitude and phase arrays must match the dimensions of x and y."))
    end

    dir_n = normalize(dir)
    # Use provided basis or fallback to automatic orthogonal basis
    e1_v = isnothing(basis) ? normal3d(dir_n) : normalize(basis[1])
    e2_v = isnothing(basis) ? normalize(cross(dir_n, e1_v)) : normalize(basis[2])

    dx = nx > 1 ? (x[2] - x[1]) : 1.0
    dy = ny > 1 ? (y[2] - y[1]) : 1.0

    # Per-axis sub-waists for correct overlap on non-square grids
    w0s_x = T(dx * overlap)
    w0s_y = T(dy * overlap)
    k = 2π / λ

    beams = Vector{AstigmaticGaussianBeamlet{T}}()
    sizehint!(beams, ceil(Int, nx * ny * 0.1)) # Conservative estimate
    beams_lock = ReentrantLock()
    max_amp = _finite_maximum(amplitude)
    max_amp > 0 || throw(ArgumentError("cannot decompose a field without nonzero samples"))
    Δφ(a, b) = _wrap_phase(phase[b...] - phase[a...])

    Threads.@threads for i in 1:nx
        for j in 1:ny
            amp = amplitude[i, j]
            if isnan(amp) || amp < max_amp * threshold
                continue
            end
            ph = phase[i, j]
            if isnan(ph)
                continue
            end

            # Local propagation direction from the Eikonal equation (shared helper, see
            # `_eikonal_direction`)
            local_dir = _eikonal_direction(Δφ, i, j, nx, ny, dx, dy, k, e1_v, e2_v, dir_n)
            if isnothing(local_dir)
                @warn lazy"Phase gradient too steep or NaN at ($i, $j); skipping."
                continue
            end

            # Position
            pos = x[i] * e1_v + y[j] * e2_v

            # Complex amplitude (E0 vector)
            # The E0 vector MUST be orthogonal to local_dir.
            # We project the macroscopic polarization (E0 or e1_v) onto the plane orthogonal to local_dir:
            base_pol = isnothing(E0) ? e1_v : E0
            # `local_dir` is real, so the (possibly complex) parallel coefficient is the
            # bilinear `dot(local_dir, base_pol)`, not `dot(base_pol, local_dir)` (its
            # conjugate for complex base_pol): `dot` conjugates its first argument.
            pol_axis = base_pol .- dot(local_dir, base_pol) .* local_dir
            if norm(pol_axis) < 1e-6
                pol_axis = e2_v .- dot(local_dir, e2_v) .* local_dir
            end

            # Normalization factor for power conservation:
            # Each sub-beamlet represents a cell of area dx*dy. The footprint of a beamlet
            # tilted by θ against `dir` is stretched by 1/cos θ on the sampling plane, so
            # the Gaussian overlap integral there is π*w0_x*w0_y / cos θ and we scale by
            # 1/S = (dx*dy) cos θ / (π*w0_x*w0_y).
            norm_factor = (dx * dy) * dot(local_dir, dir_n) / (π * w0s_x * w0s_y)
            E0_complex = normalize(pol_axis) * (amp * exp(im * ph) * norm_factor)

            # Support vector for the beamlet axes
            if randomize_axes
                # Randomly rotate the basis around the local direction
                base_s = normal3d(local_dir)
                ortho_s = cross(local_dir, base_s)
                phi = 2π * rand(rng)
                local_support = base_s * cos(phi) + ortho_s * sin(phi)
            else
                # Align the principal axes with the sampling basis so that
                # w0s_x/w0s_y belong to the e1_v/e2_v grid axes as intended.
                # Without this the constructor falls back to normal3d(local_dir),
                # which is unrelated to `basis` and transposes the two waists.
                s1 = cross(local_dir, e2_v)
                local_support = norm(s1) < 1e-6 ? nothing : normalize(s1)
            end

            b = AstigmaticGaussianBeamlet(pos, local_dir, λ, w0s_x, w0s_y; E0 = E0_complex, support = local_support)
            lock(beams_lock) do
                push!(beams, b)
            end
        end
    end

    # grid coordinates `x`, `y` are given relative to the global origin
    e1_o = isnothing(basis) ? e1_v : sampling_basis(dir_n, basis[1], T)
    return AstigmaticBeamGroup(beams, zeros(T, 3), _group_orientation(dir_n, e1_o, T))
end

"""
    WavefrontBeamletDecomposition(x, y, Eu, Ev, dir, λ; threshold=1e-4)

Vector-field counterpart of [`WavefrontBeamletDecomposition`](@ref)`(x, y, amplitude,
phase, dir, λ)`: decomposes a transverse complex vector field, given as its two
components `Eu`, `Ev` along the sampling axes `e1`, `e2` (`basis`, default `e1 =
normal3d(dir)`, `e2 = dir × e1`), into a collection of `AstigmaticGaussianBeamlet`s.

Unlike the scalar method (one amplitude/phase pair plus a single macroscopic
polarization `E0`), this method carries an independent complex amplitude per polarization
component and grid point, i.e. an arbitrary, spatially varying polarization state
(including elliptical and depolarizing-looking patterns coming from separate `Eu`, `Ev`
measurements). This is the minimal generic addition needed for
`BeamletOptics.WavefrontBeamletDecomposition(::OpticsBase.PlaneField)`
(`ext/BeamletOpticsOpticsBaseExt.jl`), which hands in the plane's tangential `E` field
directly; it is not extension-specific and works standalone.

# Arguments

- `x`, `y`: 1D spatial coordinates of the field grid, along `e1`, `e2`, in \\[m\\].
- `Eu`, `Ev`: `length(x) × length(y)` complex matrices, the field components along `e1`,
  `e2` in \\[V/m\\] (peak amplitude). Both must have the full spatial phase.
- `dir`: macroscopic reference propagation direction (e.g. `[0, 1, 0]`).
- `λ`: wavelength in \\[m\\].
- `threshold`: relative amplitude (of `hypot.(abs.(Eu), abs.(Ev))`) below which no beamlet
  is spawned.
- `overlap`, `basis`, `randomize_axes`, `rng`: as in the scalar method.
- `R`: radius in \\[m\\] of a reference sphere that has been removed from `Eu`, `Ev`:
  the physical field is `E .* cis(sign(R) k (√(x² + y² + R²) − |R|))` with `k = 2π/λ`
  (`R > 0`: diverging from the point at the distance `R` behind `x = y = 0`, `R < 0`:
  converging to the point `|R|` in front of it). Default `Inf`: `Eu`, `Ev` are the
  physical field. The sphere enters analytically, so the grid only has to resolve the
  phase of `Eu`, `Ev`, not that of the sphere: every beamlet travels along the ray of
  the sphere through its sample (plus the gradient of the sampled phase) and starts with
  the wavefront curvature of the sphere.

# Local propagation direction and polarization

The local propagation direction at each grid point is found from the Eikonal equation,
with the phase difference between neighboring samples taken from the complex vector
overlap `E(a)ᴴ E(b)`. This is continuous where the polarization changes and assumes that
both components share one local wavevector, i.e. that the field is not the superposition
of two waves with genuinely different directions at that point. The phase step between
neighbors must stay below `π`: a field tilted by `θ` against `dir` needs
`dx < λ / (2 sin θ)`, while the beamlets (waist `overlap * dx`) must stay much larger
than `λ`, so keep the sampling plane close to normal to the beam.

The beamlet's `E0` keeps the given tangential components `Eu`, `Ev` and gets the normal
component that makes it transverse to the local direction (`E·d = 0`). Its amplitude is
scaled by `cos θ = d·dir`, because the footprint of a tilted beamlet on the sampling plane
is stretched by `1/cos θ`. The sum of the beamlets then reproduces the sampled tangential
field for oblique beams as well.

Throws an `ArgumentError` for a field that is zero everywhere and for `R = 0`.

The grid is placed in the plane normal to `dir` through the global origin; use
[`translate3d!`](@ref) to move it to its actual position.
"""
function WavefrontBeamletDecomposition(
        x::AbstractVector{P1},
        y::AbstractVector{P2},
        Eu::AbstractMatrix{<:Complex},
        Ev::AbstractMatrix{<:Complex},
        dir::AbstractArray{D},
        λ::L;
        threshold::Float64 = 1e-4,
        overlap::Float64 = 1.2,
        basis::Union{Nothing, Tuple{AbstractVector, AbstractVector}} = nothing,
        randomize_axes::Bool = false,
        rng = Random.GLOBAL_RNG,
        R::Real = Inf
) where {P1 <: Real, P2 <: Real, D <: Real, L <: Real}
    T = promote_type(P1, P2, real(eltype(Eu)), real(eltype(Ev)), D, L)

    nx, ny = length(x), length(y)
    if size(Eu) != (nx, ny) || size(Ev) != (nx, ny)
        throw(DimensionMismatch("Eu and Ev arrays must match the dimensions of x and y."))
    end
    _check_reference_radius(R)

    dir_n = normalize(dir)
    e1_v = isnothing(basis) ? normal3d(dir_n) : normalize(basis[1])
    e2_v = isnothing(basis) ? normalize(cross(dir_n, e1_v)) : normalize(basis[2])

    dx = nx > 1 ? (x[2] - x[1]) : 1.0
    dy = ny > 1 ? (y[2] - y[1]) : 1.0

    w0s_x = T(dx * overlap)
    w0s_y = T(dy * overlap)
    k = 2π / λ

    amp = sqrt.(abs2.(Eu) .+ abs2.(Ev))
    # Phase difference between neighbors from the complex vector overlap E(a)ᴴ E(b): it
    # does not depend on which component dominates, so it stays continuous where the
    # polarization changes. The reference sphere is not part of it: its phase and phase
    # gradient enter analytically below, so the grid need not resolve them.
    Δφ(a, b) = angle(conj(Eu[a...]) * Eu[b...] + conj(Ev[a...]) * Ev[b...])

    beams = Vector{AstigmaticGaussianBeamlet{T}}()
    sizehint!(beams, ceil(Int, nx * ny * 0.1))
    beams_lock = ReentrantLock()
    max_amp = _finite_maximum(amp)
    max_amp > 0 || throw(ArgumentError("cannot decompose a field without nonzero samples"))

    Threads.@threads for i in 1:nx
        for j in 1:ny
            a = amp[i, j]
            if isnan(a) || a < max_amp * threshold
                continue
            end
            local_dir = _eikonal_direction(Δφ, i, j, nx, ny, dx, dy, k, e1_v, e2_v, dir_n,
                _sphere_phase_gradient(x[i], y[j], R, k))
            if isnothing(local_dir)
                @warn lazy"Phase gradient too steep or NaN at ($i, $j); skipping."
                continue
            end

            pos = x[i] * e1_v + y[j] * e2_v

            # Normalization factor for power conservation, including the cos θ of the
            # stretched footprint of a tilted beamlet (see the scalar method)
            norm_factor = (dx * dy) * dot(local_dir, dir_n) / (π * w0s_x * w0s_y)
            E0_complex = _transverse_field(Eu[i, j], Ev[i, j], e1_v, e2_v, dir_n, local_dir) .*
                         (norm_factor * cis(_sphere_phase(x[i], y[j], R, k)))

            if randomize_axes
                base_s = normal3d(local_dir)
                ortho_s = cross(local_dir, base_s)
                phi = 2π * rand(rng)
                local_support = base_s * cos(phi) + ortho_s * sin(phi)
            else
                s1 = cross(local_dir, e2_v)
                local_support = norm(s1) < 1e-6 ? nothing : normalize(s1)
            end

            # The beamlet starts with the radii w0s and the wavefront curvature of the
            # reference sphere at its sample (a waist at the sample for R = Inf)
            c_ref = sign(R) / sqrt(x[i]^2 + y[j]^2 + R^2)
            w0_x, z0_x = _waist(w0s_x, c_ref, λ)
            w0_y, z0_y = _waist(w0s_y, c_ref, λ)
            b = AstigmaticGaussianBeamlet(pos, local_dir, λ, w0_x, w0_y; E0 = E0_complex,
                support = local_support, z0_x, z0_y)
            lock(beams_lock) do
                push!(beams, b)
            end
        end
    end

    e1_o = isnothing(basis) ? e1_v : sampling_basis(dir_n, basis[1], T)
    return AstigmaticBeamGroup(beams, zeros(T, 3), _group_orientation(dir_n, e1_o, T))
end

"""
    GaussianModeDecomposition(x, y, Eu, Ev, dir, λ; basis = nothing, origin = zeros(3), R = Inf)

Fits a single [`AstigmaticGaussianBeamlet`](@ref) to a transverse complex vector field,
given by its components `Eu`, `Ev` along the sampling axes `e1`, `e2` of a plane normal to
`dir`. This is the lowest-order Gaussian mode decomposition: the counterpart of
[`WavefrontBeamletDecomposition`](@ref) for fields that are close to one Gaussian mode and
too small to be tiled into many beamlets, e.g. the output of a single-mode fiber
(beamlets must be much larger than λ).

# Arguments

- `x`, `y`: 1D sample coordinates along `e1`, `e2` in \\[m\\], relative to `origin`.
- `Eu`, `Ev`: `length(x) × length(y)` complex field components along `e1`, `e2` in
  \\[V/m\\] (peak amplitude), with the full spatial phase and travelling along `dir`.
- `dir`: propagation direction, the normal of the sampling plane.
- `λ`: vacuum wavelength in \\[m\\]. The plane must lie in vacuum or air (`n = 1`), where
  BMO starts beams.
- `basis`: optional tuple `(ex, ey)` of the sampling axes `e1`, `e2`, as in
  [`WavefrontBeamletDecomposition`](@ref). Default: `e1 = normal3d(dir)`, `e2 = dir × e1`.
- `origin`: global position of the sample `x = y = 0` in \\[m\\].
- `R`: radius in \\[m\\] of a reference sphere that has been removed from `Eu`, `Ev`, as
  in [`WavefrontBeamletDecomposition`](@ref)`(x, y, Eu, Ev, dir, λ)`. Default `Inf`:
  `Eu`, `Ev` are the physical field. The sphere enters the fit analytically (tilt at the
  centroid, curvature `1/R`), so the grid only has to resolve the phase of `Eu`, `Ev`.

# Fit

1. Polarization: the dominant eigenvector of the coherence matrix `Σ E Eᴴ`; the field
   is projected onto it.
2. Position and waist: centroid and principal axes of the second intensity moments;
   the beam radius along a principal axis is twice the rms width.
3. Direction: the mean transverse wave vector, from the intensity-weighted phase step
   between neighboring samples.
4. Curvature: the quadratic phase from the moments `⟨rᵢ ∂ⱼφ⟩` of the phase gradient
   (phase steps between neighbors), taken along the principal axes of the intensity. Twist between curvature and intensity axes
   (general astigmatism) is not represented.
5. Amplitude and phase: the projection of the field onto the fitted mode.

The beamlet starts at the centroid on the plane, with waist radii and waist positions
along its principal axes that reproduce the fitted radii and curvatures there. The
widths are measured in the plane, so the tilt of the fitted direction against `dir` must
be small.

The beamlet carries only the projected part of the field: `optical_power(beamlet)`
divided by the power of the input is the fraction of the power in the fitted Gaussian
mode, i.e. the quality of the approximation.
"""
function GaussianModeDecomposition(
        x::AbstractVector{<:Real},
        y::AbstractVector{<:Real},
        Eu::AbstractMatrix{<:Number},
        Ev::AbstractMatrix{<:Number},
        dir::AbstractVector{<:Real},
        λ::Real;
        basis::Union{Nothing, Tuple{AbstractVector, AbstractVector}} = nothing,
        origin::AbstractVector{<:Real} = zeros(3),
        R::Real = Inf
)
    nx, ny = length(x), length(y)
    size(Eu) == size(Ev) == (nx, ny) ||
        throw(DimensionMismatch("Eu and Ev arrays must match the dimensions of x and y."))
    nx > 2 && ny > 2 || throw(ArgumentError("the field needs at least 3 × 3 samples"))
    _check_reference_radius(R)

    dir_n = normalize(dir)
    e1 = isnothing(basis) ? normal3d(dir_n) : normalize(basis[1])
    e2 = isnothing(basis) ? normalize(cross(dir_n, e1)) : normalize(basis[2])
    dx, dy = x[2] - x[1], y[2] - y[1]
    k = 2π / λ

    # 1. dominant polarization (Jones vector in e1, e2) and the projected scalar field
    # Σ E Eᴴ: element (1, 2) is Σ Eu conj(Ev) = dot(Ev, Eu), since `dot` conjugates its
    # first argument
    Jm = [sum(abs2, Eu) dot(Ev, Eu); dot(Eu, Ev) sum(abs2, Ev)]
    p = eigen(Hermitian(Jm)).vectors[:, 2]
    ψ = conj(p[1]) .* Eu .+ conj(p[2]) .* Ev
    # Reference sphere: its paraxial part k ρ²/(2R) about x = y = 0, which the grid need
    # not resolve, enters the fit analytically below (tilt at the centroid, curvature 1/R,
    # phase at the centroid). Only the smooth rest is applied to the samples.
    ψ = ψ .* cis.(_sphere_phase.(x, y', R, k) .- k .* (x .^ 2 .+ y' .^ 2) ./ (2R))

    # 2. centroid and second moments of the intensity
    I_ = abs2.(ψ)
    P = sum(I_)
    P > 0 || throw(ArgumentError("cannot fit a Gaussian mode to a zero field"))
    xc = sum(I_ .* x) / P
    yc = sum(I_ .* y') / P
    X = x .- xc
    Y = y' .- yc
    M = [sum(I_ .* X .^ 2) sum(I_ .* X .* Y); sum(I_ .* X .* Y) sum(I_ .* Y .^ 2)] ./ P

    # 3. mean transverse wave vector from the phase step between neighbors, which is
    # exact for a linear phase (finite differences of ψ would give sin(kx Δ)/Δ)
    kx = angle(sum(conj.(ψ[1:(end - 1), :]) .* ψ[2:end, :])) / dx
    ky = angle(sum(conj.(ψ[:, 1:(end - 1)]) .* ψ[:, 2:end])) / dy
    # wave vector of the physical field: the reference sphere adds its gradient at the
    # centroid
    Kx, Ky = kx + k * xc / R, ky + k * yc / R
    Kx^2 + Ky^2 < k^2 || throw(ArgumentError("the mean phase gradient exceeds the wavenumber"))

    # 4. quadratic phase (k/2) rᵀ C r of the untilted field: ⟨rᵢ ∂ⱼφ⟩ = k (M C)ᵢⱼ. The
    # phase gradient is the phase step between neighbors, taken at their midpoint, which
    # is exact for a quadratic phase (finite differences of ψ would be biased by
    # sin(Δφ)/Δφ where the wavefront is strongly curved).
    ψt = ψ .* cis.(-(kx .* X .+ ky .* Y))
    Sx = conj.(ψt[1:(end - 1), :]) .* ψt[2:end, :]      # neighbors along x
    Sy = conj.(ψt[:, 1:(end - 1)]) .* ψt[:, 2:end]      # neighbors along y
    ∂x, ∂y = angle.(Sx) ./ dx, angle.(Sy) ./ dy
    Wx, Wy = abs.(Sx), abs.(Sy)                         # ≈ intensity at the midpoints
    Xx, Yy = (X[1:(end - 1)] .+ X[2:end]) ./ 2, (Y[1:(end - 1)] .+ Y[2:end]) ./ 2
    moment(W, R, D) = sum(W .* R .* D) / sum(W)
    G = [moment(Wx, Xx, ∂x) moment(Wy, X, ∂y)
         moment(Wx, Y, ∂x) moment(Wy, Yy', ∂y)]
    C = M \ G ./ k
    C = (C + C') ./ 2

    # principal axes of the intensity; beam radius w = 2 σ, curvature 1/R along each axis
    σ2, A = eigen(Symmetric(M))
    ws = 2 .* sqrt.(σ2)
    curv = [A[:, i]' * C * A[:, i] for i in 1:2]
    # waists of the physical field: the reference sphere adds 1/R to both curvatures
    w0_1, z0_1 = _waist(ws[1], curv[1] + 1 / R, λ)
    w0_2, z0_2 = _waist(ws[2], curv[2] + 1 / R, λ)

    # 5. projection onto the fitted mode (unit amplitude and zero phase at the centroid),
    # taken on the samples; the reference sphere adds its phase at the centroid
    ξ1 = A[1, 1] .* X .+ A[2, 1] .* Y
    ξ2 = A[1, 2] .* X .+ A[2, 2] .* Y
    g = exp.(-ξ1 .^ 2 ./ ws[1]^2 .- ξ2 .^ 2 ./ ws[2]^2) .*
        cis.(k .* (curv[1] .* ξ1 .^ 2 .+ curv[2] .* ξ2 .^ 2) ./ 2 .+ kx .* X .+ ky .* Y)
    a = dot(g, ψ) / sum(abs2, g) * cis(k * (xc^2 + yc^2) / (2R))

    d = normalize(Kx / k * e1 + Ky / k * e2 + sqrt(1 - (Kx^2 + Ky^2) / k^2) * dir_n)
    s = A[1, 1] * e1 + A[2, 1] * e2
    support = normalize(s - dot(s, d) * d)
    position = origin + xc * e1 + yc * e2
    E0 = _transverse_field(a * p[1], a * p[2], e1, e2, dir_n, d)
    return AstigmaticGaussianBeamlet(position, d, λ, w0_1, w0_2;
        E0, support, z0_x = z0_1, z0_y = z0_2)
end
