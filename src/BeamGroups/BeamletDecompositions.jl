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
- `n_grid`: Number of beamlets along one axis (default `20` yields `400` total beamlets). The grid
  spans `4 w0`; the sub-waist `overlap ⋅ 4 w0 / n_grid` of the beamlets must be smaller than `w0`,
  i.e. `n_grid > 4 overlap`, else an `ArgumentError` is thrown. The superposition of the beamlets
  is the Gaussian beam of the waist `w0` and the power `P0`: the sampled envelope is deconvolved
  with the profile of a beamlet.
- `overlap`: Scaling factor for the sub-waist relative to grid spacing (default `1.2` ensures smooth overlap).
- `basis`: Optional tuple `(ex, ey)` of the macroscopic sampling grid axes, i.e. the 3D directions corresponding to the `x` and `y` grid axes. `ex` must not be zero or parallel to `dir`; its component normal to `dir` becomes the local x-axis of the group [`orientation`](@ref).
- `randomize_axes`: If `true`, the internal principal axes of each individual beamlet are randomly rotated. This averages out numerical grid-alignment biases and is essential for preserving rotational symmetry in focused spots.
- `rng`: Random number generator to use for `randomize_axes`.
- `P0`: Total power of the macroscopic Gaussian beam in [W].
- `E0`: Optional Jones vector defining the polarization and initial phase of the beam; its component transverse to `dir` is used, its magnitude does not matter, the power is `P0`. If `nothing`, defaults to linear polarization along the first grid axis.
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
    dir_n = normalize(T.(dir))
    e1, e2 = grid_axes(dir_n, basis)

    # Macroscopic peak amplitude, from P0 also for a given polarization E0
    I0 = (2 * P0) / (π * w0^2)
    E_phasor = T(electric_field(I0))
    pol = isnothing(E0) ? e1 : jones_field(E0, dir_n, one(T))

    # Tile 4 w0 × 4 w0 to capture the tails
    D = 4 * w0
    Δ = D / n_grid

    # Sub-waist slightly larger than the grid spacing
    w0s = Δ * overlap
    w0s < w0 || throw(ArgumentError("the sub-waist $(w0s) of the beamlets must be smaller than the waist w0 = $(w0), increase n_grid or decrease overlap"))
    # Deconvolve the envelope by the beamlet profile
    w_env = sqrt(w0^2 - w0s^2)
    E_phasor *= w0^2 / w_env^2

    beams = Vector{AstigmaticGaussianBeamlet{T}}()
    xs = LinRange(-D / 2 + Δ / 2, D / 2 - Δ / 2, n_grid)
    ys = LinRange(-D / 2 + Δ / 2, D / 2 - Δ / 2, n_grid)

    for x in xs
        for y in ys
            offset = x * e1 + y * e2
            # Deconvolved envelope
            amp = exp(-(x^2 + y^2) / w_env^2)
            # Skip beamlets below the threshold
            if amp > threshold
                # Normalize for coherent superposition
                norm_factor = (Δ^2) / (π * w0s^2)
                local_support = randomize_axes ? random_support(dir_n, rng) : nothing
                b = AstigmaticGaussianBeamlet(
                    pos + offset, dir_n, λ, w0s; E0 = (amp * norm_factor * E_phasor) .* pol, support = local_support)
                push!(beams, b)
            end
        end
    end
    return AstigmaticBeamGroup(beams, pos, _group_orientation(dir_n, sampling_basis(dir_n, basis, T), T))
end

"""
    eikonal_direction(Δφ, i, j, nx, ny, dx, dy, k, e1, e2, dir_n, ∇φ0 = (0, 0))

Local propagation direction at the grid index `(i, j)` from the phase steps `Δφ(a, b)` between
neighbors plus the known gradient `∇φ0`, or `nothing` if it is grazing, evanescent or `NaN`.
"""
function eikonal_direction(Δφ::Function, i::Integer, j::Integer, nx::Integer, ny::Integer, dx, dy, k, e1, e2, dir_n, ∇φ0 = (0, 0))
    ip, im_ = min(i + 1, nx), max(i - 1, 1)
    jp, jm = min(j + 1, ny), max(j - 1, 1)
    dφ_dx = (Δφ((i, j), (ip, j)) + Δφ((im_, j), (i, j))) / ((ip - im_) * dx) + ∇φ0[1]
    dφ_dy = (Δφ((i, j), (i, jp)) + Δφ((i, jm), (i, j))) / ((jp - jm) * dy) + ∇φ0[2]

    sin_θx = dφ_dx / k
    sin_θy = dφ_dy / k
    # Reject grazing directions too
    if isnan(sin_θx) || isnan(sin_θy) || (sin_θx^2 + sin_θy^2 >= 1.0)
        return nothing
    end
    cos_θz = sqrt(1.0 - sin_θx^2 - sin_θy^2)
    return normalize(sin_θx * e1 + sin_θy * e2 + cos_θz * dir_n)
end

"""
    transverse_field(E1, E2, e1, e2, n, d)

The field with the tangential components `E1`, `E2` along `e1`, `e2` of a plane with the normal
`n`, completed by the normal component that makes it transverse to the direction `d`.
"""
function transverse_field(E1::Number, E2, e1, e2, n, d)
    En = -(dot(d, e1) * E1 + dot(d, e2) * E2) / dot(d, n)
    return E1 * e1 + E2 * e2 + En * n
end

"""
    sphere_phase(x, y, R, k)

Phase `sign(R) k (√(ρ² + R²) − |R|)` of a reference sphere with the radius `R` about the axis
through `x = y = 0`, without cancellation for `ρ ≪ |R|`. Zero for `R = Inf`.
"""
sphere_phase(x::Real, y, R, k) = sign(R) * k * (x^2 + y^2) / (sqrt(x^2 + y^2 + R^2) + abs(R))

function check_reference_radius(R::Real)
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
    wavefront_beamlets(x, y, amp, pos, dir, λ, T, Δφ, ∇φ0, skip, field, waist; threshold, overlap, basis, randomize_axes, rng)

The beamlet group of a sampled field with the amplitude `amp` on the grid `x`, `y` at `pos`.
The calling method supplies the phase step `Δφ(a, b)`, the known gradient `∇φ0(i, j, k)`, an
extra `skip(i, j)`, the field `field(i, j, d, c, axes, k)` scaled by `c` and the waist
`waist(w0s, i, j) -> (w0, z0)` of a sub-waist `w0s`.
"""
function wavefront_beamlets(
        x::AbstractVector,
        y::AbstractVector,
        amp::AbstractMatrix,
        pos::AbstractArray,
        dir::AbstractArray,
        λ::Real,
        ::Type{T},
        Δφ,
        ∇φ0,
        skip,
        field,
        waist;
        threshold::Float64,
        overlap::Float64,
        basis::Union{Nothing, Tuple{AbstractVector, AbstractVector}},
        randomize_axes::Bool,
        rng
) where {T <: Real}
    nx, ny = length(x), length(y)
    dir_n = normalize(T.(dir))
    e1_v, e2_v = grid_axes(dir_n, basis)

    dx = nx > 1 ? (x[2] - x[1]) : 1.0
    dy = ny > 1 ? (y[2] - y[1]) : 1.0

    # Per-axis sub-waists for non-square grids
    w0s_x = T(dx * overlap)
    w0s_y = T(dy * overlap)
    k = wavenumber(λ)

    # Largest non-NaN amplitude
    max_amp = maximum(a for a in amp if !isnan(a); init = zero(eltype(amp)))
    max_amp > 0 || throw(ArgumentError("cannot decompose a field without nonzero samples"))

    beams = Vector{AstigmaticGaussianBeamlet{T}}()
    sizehint!(beams, ceil(Int, nx * ny * 0.1))
    beams_lock = ReentrantLock()

    Threads.@threads for i in 1:nx
        for j in 1:ny
            a = amp[i, j]
            if isnan(a) || a < max_amp * threshold || skip(i, j)
                continue
            end
            local_dir = eikonal_direction(Δφ, i, j, nx, ny, dx, dy, k, e1_v, e2_v, dir_n, ∇φ0(i, j, k))
            if isnothing(local_dir)
                @warn lazy"Phase gradient too steep or NaN at ($i, $j); skipping."
                continue
            end

            p = x[i] * e1_v + y[j] * e2_v
            # Normalize power, cos θ of the stretched footprint included
            norm_factor = (dx * dy) * dot(local_dir, dir_n) / (π * w0s_x * w0s_y)
            E0 = field(i, j, local_dir, norm_factor, (e1_v, e2_v, dir_n), k)

            if randomize_axes
                local_support = random_support(local_dir, rng)
            else
                # Align the principal axes with the grid axes
                s1 = cross(local_dir, e2_v)
                # Support parallel to e2
                local_support = norm(s1) < 1e-6 ? nothing : normalize(s1)
            end

            w0_x, z0_x = waist(w0s_x, i, j)
            w0_y, z0_y = waist(w0s_y, i, j)
            b = AstigmaticGaussianBeamlet(p, local_dir, λ, w0_x, w0_y; E0, support = local_support, z0_x, z0_y)
            lock(beams_lock) do
                push!(beams, b)
            end
        end
    end

    # Place the grid at pos
    group = AstigmaticBeamGroup(beams, zeros(T, 3), _group_orientation(dir_n, sampling_basis(dir_n, basis, T), T))
    translate3d!(group, pos)
    return group
end

"""
    WavefrontBeamletDecomposition(x, y, amplitude, phase, pos, dir, λ; threshold=1e-4)

Decomposes an arbitrary complex scalar field (defined by a spatial `amplitude` and `phase`
distribution on a 2D grid `x` and `y`) into a collection of `AstigmaticGaussianBeamlet`s.

This function uses the Eikonal approximation to map the local phase gradient into a
local propagation direction for each beamlet. It allows users to import arbitrary,
aberrated, or custom beam profiles (e.g. from a camera or a wavefront sensor) and
propagate them through a `BeamletOptics` system.

The grid lies in the plane normal to `dir` through `pos`, the position of the sample
`x = y = 0` and the center of the group.

# Arguments
- `x`, `y`: Vectors defining the 1D spatial coordinates of the 2D field grid, relative to `pos`.
- `amplitude`: 2D array of field amplitudes `(length(x) × length(y))`.
- `phase`: 2D array of phase values in radians.
- `pos`: global position of the sample `x = y = 0` in \\[m\\].
- `dir`: The macroscopic reference propagation direction (e.g., `[0, 1, 0]`).
- `λ`: Wavelength.
- `threshold`: Relative amplitude threshold below which beamlets are not spawned (saves computation).
- `overlap`: Scaling factor for the sub-waist relative to grid spacing (default `1.2` ensures smooth overlap).
- `basis`: Optional tuple `(ex, ey)` of the macroscopic sampling grid axes, i.e. the 3D directions corresponding to the `x` and `y` input axes. `ex` must not be zero or parallel to `dir`; its component normal to `dir` becomes the local x-axis of the group [`orientation`](@ref).
- `randomize_axes`: If `true`, the internal principal axes of each individual beamlet are randomly rotated. This averages out numerical grid-alignment biases and helps preserve rotational symmetry in focused patterns.
- `rng`: Random number generator to use for `randomize_axes`.
- `E0`: Optional reference polarization vector or Jones vector. If `nothing`, defaults to linear polarization along the first grid axis.

!!! warning "Sampling"
    The phase step between neighboring samples must stay below `π`: a field tilted by `θ`
    against `dir` needs `dx < λ / (2 sin θ)`, while the beamlets (waist `overlap * dx`) must
    stay much larger than `λ`. Keep the sampling plane close to normal to the beam.
"""
function WavefrontBeamletDecomposition(
        x::AbstractVector{P1},
        y::AbstractVector{P2},
        amplitude::AbstractMatrix{A},
        phase::AbstractMatrix{Ph},
        pos::AbstractArray{Ps},
        dir::AbstractArray{D},
        λ::L;
        threshold::Float64 = 1e-4,
        overlap::Float64 = 1.2,
        basis::Union{Nothing, Tuple{AbstractVector, AbstractVector}} = nothing,
        randomize_axes::Bool = false,
        rng = Random.GLOBAL_RNG,
        E0 = nothing
) where {P1 <: Real, P2 <: Real, A <: Real, Ph <: Real, Ps <: Real, D <: Real, L <: Real}
    T = promote_type(P1, P2, A, Ph, Ps, D, L)

    nx, ny = length(x), length(y)
    if size(amplitude) != (nx, ny) || size(phase) != (nx, ny)
        throw(DimensionMismatch("Amplitude and phase arrays must match the dimensions of x and y."))
    end

    Δφ(a, b) = rem2pi(float(phase[b...] - phase[a...]), RoundNearest)
    function field(i, j, d, c, axes, k)
        e1, e2 = axes
        # Project the polarization normal to d, `dot` conjugates its first argument
        base_pol = isnothing(E0) ? e1 : E0
        pol_axis = base_pol .- dot(d, base_pol) .* d
        # Polarization axis parallel to d
        if norm(pol_axis) < 1e-6
            pol_axis = e2 .- dot(d, e2) .* d
        end
        return normalize(pol_axis) * (amplitude[i, j] * exp(im * phase[i, j]) * c)
    end
    return wavefront_beamlets(x, y, amplitude, pos, dir, λ, T, Δφ,
        (i, j, k) -> (0, 0), (i, j) -> isnan(phase[i, j]), field, (w0s, i, j) -> (w0s, 0);
        threshold, overlap, basis, randomize_axes, rng)
end

"""
    WavefrontBeamletDecomposition(x, y, Eu, Ev, pos, dir, λ; threshold=1e-4)

Vector-field counterpart of [`WavefrontBeamletDecomposition`](@ref)`(x, y, amplitude,
phase, pos, dir, λ)`: decomposes a transverse complex vector field, given as its two
components `Eu`, `Ev` along the sampling axes `e1`, `e2` (`basis`, default `e1 =
normal3d(dir)`, `e2 = dir × e1`), into a collection of `AstigmaticGaussianBeamlet`s.

Unlike the scalar method (one amplitude/phase pair plus a single macroscopic
polarization `E0`), this method carries an independent complex amplitude per polarization
component and grid point, i.e. an arbitrary, spatially varying polarization state
(including elliptical and depolarizing-looking patterns coming from separate `Eu`, `Ev`
measurements).

# Arguments

- `x`, `y`: 1D spatial coordinates of the field grid, along `e1`, `e2`, relative to `pos`,
  in \\[m\\].
- `Eu`, `Ev`: `length(x) × length(y)` complex matrices, the field components along `e1`,
  `e2` in \\[V/m\\] (peak amplitude). Both must have the full spatial phase.
- `pos`: global position of the sample `x = y = 0` in \\[m\\], the center of the group.
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

Throws an `ArgumentError` for a field that is zero everywhere and for `R = 0`.

!!! warning "Sampling"
    The phase step between neighboring samples must stay below `π`: a field tilted by `θ`
    against `dir` needs `dx < λ / (2 sin θ)`, while the beamlets (waist `overlap * dx`) must
    stay much larger than `λ`. Keep the sampling plane close to normal to the beam. Both
    components must share one local wavevector at each sample.
"""
function WavefrontBeamletDecomposition(
        x::AbstractVector{P1},
        y::AbstractVector{P2},
        Eu::AbstractMatrix{<:Complex},
        Ev::AbstractMatrix{<:Complex},
        pos::AbstractArray{Ps},
        dir::AbstractArray{D},
        λ::L;
        threshold::Float64 = 1e-4,
        overlap::Float64 = 1.2,
        basis::Union{Nothing, Tuple{AbstractVector, AbstractVector}} = nothing,
        randomize_axes::Bool = false,
        rng = Random.GLOBAL_RNG,
        R::Real = Inf
) where {P1 <: Real, P2 <: Real, Ps <: Real, D <: Real, L <: Real}
    T = promote_type(P1, P2, real(eltype(Eu)), real(eltype(Ev)), Ps, D, L)

    nx, ny = length(x), length(y)
    if size(Eu) != (nx, ny) || size(Ev) != (nx, ny)
        throw(DimensionMismatch("Eu and Ev arrays must match the dimensions of x and y."))
    end
    check_reference_radius(R)

    amp = sqrt.(abs2.(Eu) .+ abs2.(Ev))
    # Phase step from the vector overlap E(a)ᴴ E(b), without the reference sphere
    Δφ(a, b) = angle(conj(Eu[a...]) * Eu[b...] + conj(Ev[a...]) * Ev[b...])
    # Gradient of the reference sphere
    ∇φ0(i, j, k) = (sign(R) * k / sqrt(x[i]^2 + y[j]^2 + R^2)) .* (x[i], y[j])
    field(i, j, d, c, axes, k) = transverse_field(Eu[i, j], Ev[i, j], axes..., d) .*
                                 (c * cis(sphere_phase(x[i], y[j], R, k)))
    # Curvature of the reference sphere at the sample
    waist(w0s, i, j) = _waist(w0s, sign(R) / sqrt(x[i]^2 + y[j]^2 + R^2), λ)
    return wavefront_beamlets(x, y, amp, pos, dir, λ, T, Δφ, ∇φ0, (i, j) -> false, field, waist;
        threshold, overlap, basis, randomize_axes, rng)
end

"""
    GaussianModeDecomposition(x, y, Eu, Ev, pos, dir, λ; basis = nothing, R = Inf)

Fits a single [`AstigmaticGaussianBeamlet`](@ref) to a transverse complex vector field,
given by its components `Eu`, `Ev` along the sampling axes `e1`, `e2` of a plane normal to
`dir`. This is the lowest-order Gaussian mode decomposition: the counterpart of
[`WavefrontBeamletDecomposition`](@ref) for fields that are close to one Gaussian mode and
too small to be tiled into many beamlets, e.g. the output of a single-mode fiber
(beamlets must be much larger than λ).

# Arguments

- `x`, `y`: 1D sample coordinates along `e1`, `e2` in \\[m\\], relative to `pos`.
- `Eu`, `Ev`: `length(x) × length(y)` complex field components along `e1`, `e2` in
  \\[V/m\\] (peak amplitude), with the full spatial phase and travelling along `dir`.
- `pos`: global position of the sample `x = y = 0` in \\[m\\].
- `dir`: propagation direction, the normal of the sampling plane.
- `λ`: vacuum wavelength in \\[m\\]. The plane must lie in vacuum or air (`n = 1`), where
  BMO starts beams.
- `basis`: optional tuple `(ex, ey)` of the sampling axes `e1`, `e2`, as in
  [`WavefrontBeamletDecomposition`](@ref). Default: `e1 = normal3d(dir)`, `e2 = dir × e1`.
- `R`: radius in \\[m\\] of a reference sphere that has been removed from `Eu`, `Ev`, as
  in [`WavefrontBeamletDecomposition`](@ref)`(x, y, Eu, Ev, pos, dir, λ)`. Default `Inf`:
  `Eu`, `Ev` are the physical field. The sphere enters the fit analytically (tilt at the
  centroid, curvature `1/R`), so the grid only has to resolve the phase of `Eu`, `Ev`.

The beamlet starts at the centroid on the plane, with waist radii and waist positions
along its principal axes that reproduce the fitted radii and curvatures there.

The beamlet carries only the projected part of the field: `optical_power(beamlet)`
divided by the power of the input is the fraction of the power in the fitted Gaussian
mode, i.e. the quality of the approximation.

!!! warning "Limits"
    The widths are measured in the plane, so the tilt of the fitted direction against `dir`
    must be small. A twist between the axes of curvature and intensity (general astigmatism)
    is not represented.
"""
function GaussianModeDecomposition(
        x::AbstractVector{<:Real},
        y::AbstractVector{<:Real},
        Eu::AbstractMatrix{<:Number},
        Ev::AbstractMatrix{<:Number},
        pos::AbstractVector{<:Real},
        dir::AbstractVector{<:Real},
        λ::Real;
        basis::Union{Nothing, Tuple{AbstractVector, AbstractVector}} = nothing,
        R::Real = Inf
)
    nx, ny = length(x), length(y)
    size(Eu) == size(Ev) == (nx, ny) ||
        throw(DimensionMismatch("Eu and Ev arrays must match the dimensions of x and y."))
    nx > 2 && ny > 2 || throw(ArgumentError("the field needs at least 3 × 3 samples"))
    check_reference_radius(R)

    dir_n = normalize(dir)
    e1, e2 = grid_axes(dir_n, basis)
    dx, dy = x[2] - x[1], y[2] - y[1]
    k = wavenumber(λ)

    # 1. Dominant polarization, `dot` conjugates its first argument
    Jm = [sum(abs2, Eu) dot(Ev, Eu); dot(Eu, Ev) sum(abs2, Ev)]
    p = eigen(Hermitian(Jm)).vectors[:, 2]
    ψ = conj(p[1]) .* Eu .+ conj(p[2]) .* Ev
    # Reference sphere without its paraxial part
    ψ = ψ .* cis.(sphere_phase.(x, y', R, k) .- k .* (x .^ 2 .+ y' .^ 2) ./ (2R))

    # 2. Centroid and second moments of the intensity
    I_ = abs2.(ψ)
    P = sum(I_)
    P > 0 || throw(ArgumentError("cannot fit a Gaussian mode to a zero field"))
    xc = sum(I_ .* x) / P
    yc = sum(I_ .* y') / P
    X = x .- xc
    Y = y' .- yc
    M = [sum(I_ .* X .^ 2) sum(I_ .* X .* Y); sum(I_ .* X .* Y) sum(I_ .* Y .^ 2)] ./ P

    # 3. Mean transverse wave vector from the phase steps between neighbors
    kx = angle(sum(conj.(ψ[1:(end - 1), :]) .* ψ[2:end, :])) / dx
    ky = angle(sum(conj.(ψ[:, 1:(end - 1)]) .* ψ[:, 2:end])) / dy
    # Add the gradient of the reference sphere at the centroid
    Kx, Ky = kx + k * xc / R, ky + k * yc / R
    Kx^2 + Ky^2 < k^2 || throw(ArgumentError("the mean phase gradient exceeds the wavenumber"))

    # 4. Quadratic phase (k/2) rᵀ C r of the untilted field from ⟨rᵢ ∂ⱼφ⟩ = k (M C)ᵢⱼ
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

    # Principal axes of the intensity, beam radius w = 2 σ
    σ2, A = eigen(Symmetric(M))
    ws = 2 .* sqrt.(σ2)
    curv = [A[:, i]' * C * A[:, i] for i in 1:2]
    # Add the curvature of the reference sphere
    w0_1, z0_1 = _waist(ws[1], curv[1] + 1 / R, λ)
    w0_2, z0_2 = _waist(ws[2], curv[2] + 1 / R, λ)

    # 5. Projection onto the fitted mode, plus the phase of the reference sphere at the centroid
    ξ1 = A[1, 1] .* X .+ A[2, 1] .* Y
    ξ2 = A[1, 2] .* X .+ A[2, 2] .* Y
    g = exp.(-ξ1 .^ 2 ./ ws[1]^2 .- ξ2 .^ 2 ./ ws[2]^2) .*
        cis.(k .* (curv[1] .* ξ1 .^ 2 .+ curv[2] .* ξ2 .^ 2) ./ 2 .+ kx .* X .+ ky .* Y)
    a = dot(g, ψ) / sum(abs2, g) * cis(k * (xc^2 + yc^2) / (2R))

    d = normalize(Kx / k * e1 + Ky / k * e2 + sqrt(1 - (Kx^2 + Ky^2) / k^2) * dir_n)
    s = A[1, 1] * e1 + A[2, 1] * e2
    support = normalize(s - dot(s, d) * d)
    E0 = transverse_field(a * p[1], a * p[2], e1, e2, dir_n, d)
    return AstigmaticGaussianBeamlet(pos + xc * e1 + yc * e2, d, λ, w0_1, w0_2;
        E0, support, z0_x = z0_1, z0_y = z0_2)
end
