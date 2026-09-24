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
    e1_o = isnothing(basis) ? e1 : _sampling_basis(dir_n, basis[1], T)
    return AstigmaticBeamGroup(beams, pos, _group_orientation(dir_n, e1_o, T))
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
    max_amp = maximum(amplitude)

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

            # Compute local phase gradients using central differences
            # Use mod2pi wrap to avoid slow angle(exp(im*x))
            _wrap = Δ -> mod2pi(Δ + π) - π
            
            if i > 1 && i < nx
                dφ_dx = (_wrap(phase[i + 1, j] - ph) +
                         _wrap(ph - phase[i - 1, j])) / (2dx)
            elseif i == 1
                dφ_dx = _wrap(phase[i + 1, j] - ph) / dx
            else
                dφ_dx = _wrap(ph - phase[i - 1, j]) / dx
            end

            if j > 1 && j < ny
                dφ_dy = (_wrap(phase[i, j + 1] - ph) +
                         _wrap(ph - phase[i, j - 1])) / (2dy)
            elseif j == 1
                dφ_dy = _wrap(phase[i, j + 1] - ph) / dy
            else
                dφ_dy = _wrap(ph - phase[i, j - 1]) / dy
            end

            # Convert phase gradient to angular deviation (Eikonal equation: ∇φ = k * sin(θ))
            sin_θx = dφ_dx / k
            sin_θy = dφ_dy / k

            # Ensure valid angles and handle NaNs
            if isnan(sin_θx) || isnan(sin_θy) || (sin_θx^2 + sin_θy^2 > 1.0)
                @warn lazy"Phase gradient too steep or NaN at ($i, $j); skipping."
                continue
            end

            cos_θz = sqrt(1.0 - sin_θx^2 - sin_θy^2)

            # Construct the local direction vector
            local_dir = sin_θx * e1_v + sin_θy * e2_v + cos_θz * dir_n
            local_dir = normalize(local_dir)

            # Position
            pos = x[i] * e1_v + y[j] * e2_v

            # Complex amplitude (E0 vector)
            # The E0 vector MUST be orthogonal to local_dir.
            # We project the macroscopic polarization (E0 or e1_v) onto the plane orthogonal to local_dir:
            base_pol = isnothing(E0) ? e1_v : E0
            pol_axis = base_pol .- dot(base_pol, local_dir) .* local_dir
            if norm(pol_axis) < 1e-6
                pol_axis = e2_v .- dot(e2_v, local_dir) .* local_dir
            end

            # Normalization factor for power conservation:
            # Each sub-beamlet represents a cell of area dx*dy. The Gaussian overlap
            # integral is π*w0_x*w0_y, so we scale by 1/S = (dx*dy) / (π*w0_x*w0_y).
            norm_factor = (dx * dy) / (π * w0s_x * w0s_y)
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
    e1_o = isnothing(basis) ? e1_v : _sampling_basis(dir_n, basis[1], T)
    return AstigmaticBeamGroup(beams, zeros(T, 3), _group_orientation(dir_n, e1_o, T))
end
