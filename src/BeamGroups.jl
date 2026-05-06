"""
    PointSource <: AbstractBeamGroup

Represents a cone of [`Beam`](@ref)s being emitted from a single point in space.

# Fields

- `beams`: a vector of all [`Beam`](@ref)s originating from the source
- `NA`: the [`numerical_aperture`](@ref) of the point source spread angle

# Functions

- `numerical_aperture`: returns the NA of the source
"""
struct PointSource{T, R <: AbstractRay{T}} <: AbstractBeamGroup{T, R}
    beams::Vector{Beam{T, R}}
    NA::T
end

numerical_aperture(ps::PointSource) = ps.NA

"""
    PointSource(pos, dir, θ, λ; num_rings, num_rays)

Spawns a point source of [`Beam`](@ref)s at the specified `pos`ition and `dir`ection.
The point source is modelled as a collection of concentric beam fans centered around the center beam.
The amount of beam rings between the center ray and half-spread-angle `θ` can be specified via `num_rings`.

!!! info
    Note that for correct sampling, the number of rays should be atleast 20x the number of rings.

# Arguments

The following inputs and arguments can be used to configure the [`PointSource`](@ref):

## Inputs

- `pos`: center beam starting position
- `dir`: center beam starting direction
- `θ`: half spread angle in rad
- `λ = 1e-6`: wavelength in [m], default val. is 1000 nm

## Keyword Arguments

- `num_rings`: number of concentric beam rings, default is 10
- `num_rays`: total number of rays in the source, default is 100x num_rings

!!! warning
    The orthogonal basis vectors for the beam generation are generated randomly.
"""
function PointSource(
        pos::AbstractArray{P},
        dir::AbstractArray{D},
        θ::H,
        λ::L = 1e-6;
        num_rings::Int = 10,
        num_rays::Int = 100 * num_rings
) where {P <: Real, D <: Real, H <: Real, L <: Real}
    T = promote_type(P, D, H, L)
    if num_rays < num_rings * 20
        throw(ErrorException("No. of rays should be atleast 20x no. of rings (passed: $num_rays, req: $(num_rings*20))"))
    end
    if θ ≥ pi
        throw(ErrorException("Point source opening half-angle θ must be ≤ π"))
    end
    # define basis vectors
    dir = normalize(dir)
    b1 = normal3d(dir) # random seed
    b2 = normal3d(dir, b1)
    θ_NA = LinRange(0, θ, num_rings)
    # define buffer
    beams = Vector{Beam{T, Ray{T}}}()
    push!(beams, Beam(Ray(pos, dir, λ)))
    num_rays -= 1
    # calculate total accumulated circumference of all rings
    ndirs = [rotate3d(b2, step(θ_NA) * i) * dir for i in eachindex(θ_NA[2:end])]
    circm = norm.(ndirs .- dot.(ndirs, Ref(dir)) .* Ref(dir)) .* Ref(2π)
    total = sum(circm)
    ds = total / num_rays
    # calculate number of rays per ring
    n_rays = round.(Int, circm / ds)
    # correct n_rays to match num_rays
    n_rays[end] += (num_rays - sum(n_rays))
    for (i, ndir) in enumerate(ndirs)
        numEl = n_rays[i]
        if iszero(numEl)
            continue
        end
        dphi = 2π / numEl
        RotMat = rotate3d(dir, dphi)
        cdir = ndir
        for _ in 1:numEl
            push!(beams, Beam(pos, cdir, λ))
            # rotate vector (not-thread safe!)
            cdir = RotMat * cdir
        end
    end
    NA = numerical_aperture(θ)
    return PointSource(beams, NA)
end

"""
    CollimatedSource <: AbstractBeamGroup

Represents a parallel bundle of [`Beam`](@ref)s being emitted from a disk in space.

# Fields

- `beams`: a vector of all [`Beam`](@ref)s originating from the source
- `diameter`: the diameter of the outermost beam ring

# Functions

- `diameter`: returns the diameter of the source
"""
struct CollimatedSource{T, R <: AbstractRay{T}} <: AbstractBeamGroup{T, R}
    beams::Vector{Beam{T, R}}
    diameter::T
end

diameter(cs::CollimatedSource) = cs.diameter

"""
    CollimatedSource(pos, dir, diameter, λ; num_rings, num_rays)

Spawns a bundle of collimated [`Beam`](@ref)s at the specified `pos`ition and `dir`ection.
The source is modelled as a ring of concentric beam rings around the center beam.
The amount of beam rings between the center ray and outer `diameter` can be specified via `num_rings`.

!!! info
    Note that for correct sampling, the number of rays should be atleast 20x the number of rings.

# Arguments

The following inputs and arguments can be used to configure the [`CollimatedSource`](@ref):

## Inputs

- `pos`: center beam starting position
- `dir`: center beam starting direction
- `diameter`: outer beam bundle diameter in [m]
- `λ = 1e-6`: wavelength in [m], default val. is 1000 nm

## Keyword Arguments

- `num_rings`: number of concentric beam rings, default is 10
- `num_rays`: total number of rays in the source, default is 100x num_rings

!!! warning
    The orthogonal basis vectors for the beam generation are generated randomly.
"""
function CollimatedSource(
        pos::AbstractArray{P},
        dir::AbstractArray{D1},
        diameter::D2,
        λ::L = 1e-6;
        num_rings::Int = 10,
        num_rays::Int = 100 * num_rings
) where {P <: Real, D1 <: Real, D2 <: Real, L <: Real}
    T = promote_type(P, D1, D2, L)
    if num_rays < num_rings * 20
        throw(ErrorException("No. of rays should be atleast 20x no. of rings (passed: $num_rays, req: $(num_rings*20))"))
    end
    # define buffer
    beams = Vector{Beam{T, Ray{T}}}()
    push!(beams, Beam(Ray(pos, dir, λ)))
    num_rays -= 1
    # setup concentric beam ring radii
    b1 = normal3d(dir) # random seed
    r_max = diameter / 2
    radii = LinRange(0, r_max, num_rings)[2:end]
    # calculate total accumulated circumference of all rings
    circm = radii * 2π
    total = sum(circm)
    ds = total / num_rays
    # calculate number of rays per ring
    n_rays = round.(Int, circm / ds)
    # correct n_rays to match num_rays
    n_rays[end] += (num_rays - sum(n_rays))
    # Generate beam rings
    for (i, r) in enumerate(radii)
        numEl = n_rays[i]
        if iszero(numEl)
            continue
        end
        dphi = 2π / numEl
        RotMat = rotate3d(dir, dphi)
        helper = b1 * r
        for _ in 1:numEl
            push!(beams, Beam(pos + helper, dir, λ))
            helper = RotMat * helper
        end
    end
    return CollimatedSource(beams, T(diameter))
end

"""
    UniformDiscSource(pos, dir, diameter, λ; num_rays=1_000)

Generates a ray fan with *equal area per ray* across a circular pupil
using the deterministic sunflower (Fibonacci) pattern.

!!! note
    This is merely a [`CollimatedSource`](@ref) constructor which uses Fibonacci sampling
    instead of a linear grid.

# Arguments

The following inputs and arguments can be used to configure the underlying [`CollimatedSource`](@ref):

## Inputs

- `pos`: center beam starting position
- `dir`: center beam starting direction
- `diameter`: outer beam bundle diameter in [m]
- `λ = 1e-6`: wavelength in [m]

## Keyword Arguments

- `num_rays=1000`: total number of rays in the source
"""
function UniformDiscSource(
        pos::AbstractArray{P},
        dir::AbstractArray{D1},
        diameter::D2,
        λ::L = 1e-6;
        # kwargs
        num_rays::Int = 1_000
) where {P <: Real, D1 <: Real, D2 <: Real, L <: Real}
    T = promote_type(P, D1, D2, L)
    R = diameter / 2
    φ0 = 2π / (1 + √5)         # golden angle
    beams = Vector{Beam{T, Ray{T}}}(undef, num_rays)
    # orthogonal basis in the pupil plane
    e1 = normal3d(dir)
    e2 = normalize(cross(dir, e1))
    for k in 0:(num_rays - 1)
        ρ = √((k + 0.5) / num_rays)     # equal-area radius
        φ = k * φ0
        r = R * ρ
        x = r * cos(φ) * e1 + r * sin(φ) * e2
        beams[k + 1] = Beam(pos + x, dir, λ)
    end
    return CollimatedSource(beams, T(diameter))
end

"""
    AstigmaticBeamGroup{T, R} <: AbstractBeamGroup{T, R}

A generic container for groups of [`AstigmaticGaussianBeamlet`](@ref)s.
"""
struct AstigmaticBeamGroup{T, R <: AbstractRay{T}} <: AbstractBeamGroup{T, R}
    beams::Vector{AstigmaticGaussianBeamlet{T}}
end

function AstigmaticBeamGroup(beams::Vector{AstigmaticGaussianBeamlet{T}}) where {T}
    return AstigmaticBeamGroup{T, PolarizedRay{T}}(beams)
end

"""
    CollimatedGaussianBeamletSource(pos, dir, D, λ, w0s; n_grid=20)


Generates an `n_grid × n_grid` square array of `AstigmaticGaussianBeamlet`s, all pointing
in the same `dir`ection and having the same uniform amplitude. This is used to model
macroscopic plane waves or flat-top beams passing through hard apertures (e.g., for
Fraunhofer or Fresnel diffraction).

# Arguments
- `pos`: Center position of the source plane.
- `dir`: Propagation direction.
- `D`: Total width/height of the square aperture.
- `λ`: Wavelength.
- `w0s`: Sub-waist of each individual beamlet. For smooth overlap, `w0s ≈ D / n_grid`.
- `n_grid`: Number of beamlets along one axis (default `20` yields `400` total beamlets).
- `basis`: Orientation of the macroscopic sampling grid. Determines the 3D directions corresponding to the `x` and `y` grid axes.
- `randomize_axes`: If `true`, the internal principal axes of each individual beamlet are randomly rotated. This averages out numerical grid-alignment biases and is essential for preserving rotational symmetry in focused spots (e.g. Airy disks).
- `rng`: Random number generator to use for `randomize_axes`.
"""
function CollimatedGaussianBeamletSource(
        pos::AbstractArray{P},
        dir::AbstractArray{D1},
        D::D2,
        λ::L,
        w0s::W;
        n_grid::Int = 20,
        basis::Union{Nothing, Tuple{AbstractVector, AbstractVector}} = nothing,
        randomize_axes::Bool = false,
        rng = Random.GLOBAL_RNG
) where {P <: Real, D1 <: Real, D2 <: Real, L <: Real, W <: Real}
    T = promote_type(P, D1, D2, L, W)
    dir_n = normalize(dir)
    # Use provided basis or fallback to automatic orthogonal basis
    e1 = isnothing(basis) ? normal3d(dir_n) : normalize(basis[1])
    e2 = isnothing(basis) ? normalize(cross(dir_n, e1)) : normalize(basis[2])

    beams = Vector{AstigmaticGaussianBeamlet{T}}()
    Δ = D / n_grid

    xs = LinRange(-D / 2 + Δ / 2, D / 2 - Δ / 2, n_grid)
    ys = LinRange(-D / 2 + Δ / 2, D / 2 - Δ / 2, n_grid)

    for x in xs
        for y in ys
            offset = x * e1 + y * e2
            # Support vector for the beamlet axes
            if randomize_axes
                base_s = normal3d(dir_n)
                ortho_s = cross(dir_n, base_s)
                phi = 2π * rand(rng)
                local_support = base_s * cos(phi) + ortho_s * sin(phi)
            else
                local_support = nothing
            end
            b = AstigmaticGaussianBeamlet(pos + offset, dir_n, λ, w0s; support = local_support)
            push!(beams, b)
        end
    end
    return AstigmaticBeamGroup(beams)
end

"""
    GaussianBeamletDecomposition(pos, dir, w0, λ; n_grid=20)

Decomposes a macroscopic Gaussian beam of waist `w0` into an `n_grid × n_grid`
array of microscopic `AstigmaticGaussianBeamlet`s. This allows the accurate
tracing of large Gaussian beams through highly aberrative or non-paraxial
optical systems, perfectly preserving higher-order phase aberrations via the
coherent superposition of the sub-beamlets.

# Arguments
- `pos`: Waist center position of the macroscopic beam.
- `dir`: Propagation direction.
- `w0`: Macroscopic beam waist.
- `λ`: Wavelength.
- `n_grid`: Number of beamlets along one axis (default `20` yields `400` total beamlets).
- `overlap`: Scaling factor for the sub-waist relative to grid spacing (default `1.2` ensures smooth overlap).
- `basis`: Orientation of the macroscopic sampling grid. Determines the 3D directions corresponding to the `x` and `y` grid axes.
- `randomize_axes`: If `true`, the internal principal axes of each individual beamlet are randomly rotated. This averages out numerical grid-alignment biases and is essential for preserving rotational symmetry in focused spots.
- `rng`: Random number generator to use for `randomize_axes`.
- `P0`: Total power of the macroscopic Gaussian beam in [W].
- `E0`: Optional Jones vector defining the polarization and initial phase of the beam. If `nothing`, defaults to linear polarization along the first grid axis.
- `threshold`: Relative amplitude below which beamlets are not spawned (default `1e-4`).
"""
function GaussianBeamletDecomposition(
        pos::AbstractArray{P},
        dir::AbstractArray{D1},
        w0::W,
        λ::L;
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
    return AstigmaticBeamGroup(beams)
end

"""
    SphericalGaussianBeamletSource(pos, dir, θ, λ; num_rings=10, num_rays=100*num_rings)

Decomposes a macroscopic spherical wave (or a highly divergent/focused beam) into a cone
of `AstigmaticGaussianBeamlet`s originating from a single point `pos`.

The angular spread of the beamlets is bounded by the half-angle `θ`. To ensure a smooth
far-field interference pattern without speckle, the sub-waist `w0s` of each beamlet is
automatically calculated such that their far-field divergence perfectly overlaps with
adjacent beamlets in the grid.

# Arguments
- `pos`: Origin point of the spherical wave.
- `dir`: Central propagation direction.
- `θ`: Half spread angle in radians.
- `λ`: Wavelength.
- `num_rings`: Number of concentric angular rings.
- `num_rays`: Total number of beamlets to generate.
- `overlap`: Scaling factor for the sub-waist divergence (default `1.2` ensures smooth overlap).
- `basis`: Optional reference vector (e.g. `[1,0,0]`) to define the starting azimuthal angle for the source rings.
- `randomize_axes`: If `true`, the internal principal axes of each individual beamlet are randomly rotated. This averages out numerical artifacts in the far-field focus.
- `rng`: Random number generator to use for `randomize_axes`.
- `P0`: Total power of the spherical source in [W].
- `E0`: Optional Jones vector defining the polarization and phase of the spherical wave.
"""
function SphericalGaussianBeamletSource(
        pos::AbstractArray{P},
        dir::AbstractArray{D},
        θ::H,
        λ::L;
        num_rings::Int = 10,
        num_rays::Int = 100 * num_rings,
        overlap::Float64 = 1.2,
        basis::Union{Nothing, AbstractVector} = nothing,
        randomize_axes::Bool = false,
        rng = Random.GLOBAL_RNG,
        P0::Real = get_default_power(),
        E0 = nothing
) where {P <: Real, D <: Real, H <: Real, L <: Real}
    T = promote_type(P, D, H, L)
    if num_rays < num_rings * 20
        throw(ErrorException("No. of beamlets should be atleast 20x no. of rings (passed: $num_rays, req: $(num_rings*20))"))
    end
    if θ ≥ pi
        throw(ErrorException("Spherical source opening half-angle θ must be ≤ π"))
    end

    # Calculate optimal sub-waist for smooth far-field overlap
    # Angular spacing between rings:
    Δθ = θ / num_rings
    # We want beamlet divergence θ_div = λ / (π * w0s) to be approx overlap * Δθ
    w0s = λ / (π * overlap * Δθ)

    dir_n = normalize(dir)
    b1 = isnothing(basis) ? normal3d(dir_n) : normalize(basis - dot(basis, dir_n) * dir_n)
    b2 = normalize(cross(dir_n, b1))
    θ_NA = LinRange(0, θ, num_rings)

    # Initial ray count for power distribution
    total_rays = num_rays
    P_sub = P0 / total_rays

    beams = Vector{AstigmaticGaussianBeamlet{T}}()

    # Central beamlet
    push!(beams, AstigmaticGaussianBeamlet(pos, dir_n, λ, w0s; P0 = P_sub, E0 = E0))
    num_rays -= 1

    # Calculate circumference weights
    ndirs = [rotate3d(b2, step(θ_NA) * i) * dir_n for i in eachindex(θ_NA[2:end])]
    circm = norm.(ndirs .- dot.(ndirs, Ref(dir_n)) .* Ref(dir_n)) .* Ref(2π)
    total = sum(circm)
    ds = total / num_rays

    n_rays = round.(Int, circm / ds)
    n_rays[end] += (num_rays - sum(n_rays))

    for (i, ndir) in enumerate(ndirs)
        numEl = n_rays[i]
        if iszero(numEl)
            continue
        end
        dphi = 2π / numEl
        RotMat = rotate3d(dir_n, dphi)
        cdir = ndir
        for _ in 1:numEl
            # Polarization is set via the E0 kwarg. We could try to project it properly,
            # but default linear works for paraxial-ish cones.
            # Support vector for the beamlet axes
            if randomize_axes
                base_s = normal3d(cdir)
                ortho_s = cross(cdir, base_s)
                phi = 2π * rand(rng)
                local_support = base_s * cos(phi) + ortho_s * sin(phi)
            else
                local_support = nothing
            end
            push!(beams, AstigmaticGaussianBeamlet(pos, cdir, λ, w0s; P0 = P_sub, E0 = E0, support = local_support))
            cdir = RotMat * cdir
        end
    end

    return AstigmaticBeamGroup(beams)
end

"""
    WavefrontBeamletDecomposition(x, y, amplitude, phase, dir, λ; threshold=1e-4)

Decomposes an arbitrary complex scalar field (defined by a spatial `amplitude` and `phase`
distribution on a 2D grid `x` and `y`) into a collection of `AstigmaticGaussianBeamlet`s.

This function uses the Eikonal approximation to map the local phase gradient into a
local propagation direction for each beamlet. It allows users to import arbitrary,
aberrated, or custom beam profiles (e.g. from a camera or a wavefront sensor) and
propagate them through a `BeamletOptics` system.

# Arguments
- `x`, `y`: Vectors defining the 1D spatial coordinates of the 2D field grid.
- `amplitude`: 2D array of field amplitudes `(length(x) × length(y))`.
- `phase`: 2D array of phase values in radians.
- `dir`: The macroscopic reference propagation direction (e.g., `[0, 1, 0]`).
- `λ`: Wavelength.
- `threshold`: Relative amplitude threshold below which beamlets are not spawned (saves computation).
- `overlap`: Scaling factor for the sub-waist relative to grid spacing (default `1.2` ensures smooth overlap).
- `basis`: Orientation of the macroscopic sampling grid. Determines the 3D directions corresponding to the `x` and `y` input axes.
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
    max_amp = maximum(amplitude)

    for i in 1:nx
        for j in 1:ny
            amp = amplitude[i, j]
            ph = phase[i, j]
            if isnan(amp) || isnan(ph) || amp < max_amp * threshold
                continue
            end

            # Compute local phase gradients using central differences (fallback to forward/backward at edges)
            if i > 1 && i < nx
                dφ_dx = (angle(exp(im * (phase[i + 1, j] - phase[i, j]))) +
                         angle(exp(im * (phase[i, j] - phase[i - 1, j])))) / (2dx)
            elseif i == 1
                dφ_dx = angle(exp(im * (phase[i + 1, j] - phase[i, j]))) / dx
            else
                dφ_dx = angle(exp(im * (phase[i, j] - phase[i - 1, j]))) / dx
            end

            if j > 1 && j < ny
                dφ_dy = (angle(exp(im * (phase[i, j + 1] - phase[i, j]))) +
                         angle(exp(im * (phase[i, j] - phase[i, j - 1])))) / (2dy)
            elseif j == 1
                dφ_dy = angle(exp(im * (phase[i, j + 1] - phase[i, j]))) / dy
            else
                dφ_dy = angle(exp(im * (phase[i, j] - phase[i, j - 1]))) / dy
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
                local_support = nothing
            end

            b = AstigmaticGaussianBeamlet(pos, local_dir, λ, w0s_x, w0s_y; E0 = E0_complex, support = local_support)
            push!(beams, b)
        end
    end

    return AstigmaticBeamGroup(beams)
end
