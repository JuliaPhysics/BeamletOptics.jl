"""
    PointSource <: AbstractBeamGroup

Represents a cone of [`Beam`](@ref)s being emitted from a single point in space.

# Fields

- `beams`: a vector of all [`Beam`](@ref)s originating from the source
- `NA`: the [`numerical_aperture`](@ref) of the point source spread angle

# Functions

- `numerical_aperture`: returns the NA of the source
"""
struct PointSource{T, R<:AbstractRay{T}} <: AbstractBeamGroup{T,R}
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
        num_rings::Int=10,
        num_rays::Int=100*num_rings,
    ) where {P<:Real, D<:Real, H<:Real, L<:Real}
    T = promote_type(P, D, H, L)
    if num_rays < num_rings*20
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
    ndirs = [rotate3d(b2, step(θ_NA)*i) * dir for i in eachindex(θ_NA[2:end])]
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
        for _ = 1:numEl
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
struct CollimatedSource{T, R<:AbstractRay{T}} <: AbstractBeamGroup{T,R}
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
        num_rings::Int=10,
        num_rays::Int=100*num_rings,
    ) where {P<:Real, D1<:Real, D2<:Real, L<:Real}
    T = promote_type(P, D1, D2, L)
    if num_rays < num_rings*20
        throw(ErrorException("No. of rays should be atleast 20x no. of rings (passed: $num_rays, req: $(num_rings*20))"))
    end
    # define buffer
    beams = Vector{Beam{T, Ray{T}}}()
    push!(beams, Beam(Ray(pos, dir, λ)))
    num_rays -= 1
    # setup concentric beam ring radii
    b1 = normal3d(dir) # random seed
    r_max = diameter/2
    radii = LinRange(0, r_max, num_rings)[2:end]
    # calculate total accumulated circumference of all rings
    circm = radii*2π
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
        for _ = 1:numEl
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
        λ::L=1e-6;
        # kwargs
        num_rays::Int=1_000
    ) where {P<:Real, D1<:Real, D2<:Real, L<:Real}
    T  = promote_type(P, D1, D2, L)
    R  = diameter / 2
    φ0 = 2π / (1 + √5)         # golden angle
    beams = Vector{Beam{T,Ray{T}}}(undef, num_rays)
    # orthogonal basis in the pupil plane
    e1 = normal3d(dir)
    e2 = normalize(cross(dir, e1))
    for k in 0:num_rays-1
        ρ  = √((k + 0.5)/num_rays)     # equal-area radius
        φ  = k * φ0
        r  = R * ρ
        x  =  r * cos(φ) * e1 + r * sin(φ) * e2
        beams[k+1] = Beam(pos + x, dir, λ)
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
"""
function CollimatedGaussianBeamletSource(
        pos::AbstractArray{P},
        dir::AbstractArray{D1},
        D::D2,
        λ::L,
        w0s::W;
        n_grid::Int=20
    ) where {P<:Real, D1<:Real, D2<:Real, L<:Real, W<:Real}
    T = promote_type(P, D1, D2, L, W)
    dir_n = normalize(dir)
    # Create orthogonal basis for the square grid
    e1 = normal3d(dir_n)
    e2 = normalize(cross(dir_n, e1))
    
    beams = Vector{AstigmaticGaussianBeamlet{T}}()
    Δ = D / n_grid
    
    xs = LinRange(-D / 2 + Δ / 2, D / 2 - Δ / 2, n_grid)
    ys = LinRange(-D / 2 + Δ / 2, D / 2 - Δ / 2, n_grid)
    
    for x in xs
        for y in ys
            offset = x * e1 + y * e2
            b = AstigmaticGaussianBeamlet(pos + offset, dir_n, λ, w0s)
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
"""
function GaussianBeamletDecomposition(
        pos::AbstractArray{P},
        dir::AbstractArray{D1},
        w0::W,
        λ::L;
        n_grid::Int=20,
        overlap::Float64=1.2
    ) where {P<:Real, D1<:Real, W<:Real, L<:Real}
    T = promote_type(P, D1, W, L)
    dir_n = normalize(dir)
    e1 = normal3d(dir_n)
    e2 = normalize(cross(dir_n, e1))
    
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
            if amp > 1e-4
                # Normalization factor for coherent superposition
                norm_factor = (Δ^2) / (π * w0s^2)
                # Standard linear polarization along e1
                b = AstigmaticGaussianBeamlet(pos + offset, dir_n, λ, w0s; E0 = amp * e1 * norm_factor)
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
"""
function SphericalGaussianBeamletSource(
        pos::AbstractArray{P},
        dir::AbstractArray{D},
        θ::H,
        λ::L;
        num_rings::Int=10,
        num_rays::Int=100*num_rings,
        overlap::Float64=1.2
    ) where {P<:Real, D<:Real, H<:Real, L<:Real}
    T = promote_type(P, D, H, L)
    if num_rays < num_rings*20
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
    b1 = normal3d(dir_n)
    b2 = normal3d(dir_n, b1)
    θ_NA = LinRange(0, θ, num_rings)
    
    beams = Vector{AstigmaticGaussianBeamlet{T}}()
    
    # Central beamlet
    push!(beams, AstigmaticGaussianBeamlet(pos, dir_n, λ, w0s))
    num_rays -= 1
    
    # Calculate circumference weights
    ndirs = [rotate3d(b2, step(θ_NA)*i) * dir_n for i in eachindex(θ_NA[2:end])]
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
        for _ = 1:numEl
            # Polarization is set via the E0 kwarg. We could try to project it properly, 
            # but default linear works for paraxial-ish cones.
            push!(beams, AstigmaticGaussianBeamlet(pos, cdir, λ, w0s))
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
- `basis`: Optional `Tuple(e1, e2)` of orthogonal vectors defining the local coordinate system of the input grid. If `nothing`, an automatic basis is generated.
"""
function WavefrontBeamletDecomposition(
    x::AbstractVector{P1},
    y::AbstractVector{P2},
    amplitude::AbstractMatrix{A},
    phase::AbstractMatrix{Ph},
    dir::AbstractArray{D},
    λ::L;
    threshold::Float64=1e-4,
    overlap::Float64=1.2,
    basis::Union{Nothing, Tuple{AbstractVector, AbstractVector}}=nothing
) where {P1<:Real, P2<:Real, A<:Real, Ph<:Real, D<:Real, L<:Real}
    
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
    
    # The waist is chosen to perfectly overlap with adjacent grid points
    w0s = T(max(dx, dy) * overlap)
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
            # Compute local phase gradients using phase-unwrapped differences
            if i > 1 && i < nx
                dφ_dx = (angle(exp(im * (phase[i+1, j] - phase[i, j]))) + angle(exp(im * (phase[i, j] - phase[i-1, j])))) / (2dx)
            elseif i == 1
                dφ_dx = angle(exp(im * (phase[i+1, j] - phase[i, j]))) / dx
            else
                dφ_dx = angle(exp(im * (phase[i, j] - phase[i-1, j]))) / dx
            end

            if j > 1 && j < ny
                dφ_dy = (angle(exp(im * (phase[i, j+1] - phase[i, j]))) + angle(exp(im * (phase[i, j] - phase[i, j-1])))) / (2dy)
            elseif j == 1
                dφ_dy = angle(exp(im * (phase[i, j+1] - phase[i, j]))) / dy
            else
                dφ_dy = angle(exp(im * (phase[i, j] - phase[i, j-1]))) / dy
            end
            
            # Convert phase gradient to angular deviation (Eikonal equation: ∇φ = k * sin(θ))
            sin_θx = dφ_dx / k
            sin_θy = dφ_dy / k
            
            # Ensure valid angles and handle NaNs
            if isnan(sin_θx) || isnan(sin_θy) || (sin_θx^2 + sin_θy^2 > 1.0)
                @warn "Phase gradient too steep or NaN at ($i, $j); skipping."
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
            # We project the macroscopic polarization (e1_v) onto the plane orthogonal to local_dir:
            pol_axis = e1_v .- dot(e1_v, local_dir) .* local_dir
            if norm(pol_axis) < 1e-6
                pol_axis = e2_v .- dot(e2_v, local_dir) .* local_dir
            end
            
            # Normalization factor for power conservation:
            # When summing Gaussians exp(-r^2/w0^2) on a grid with spacing dx, dy,
            # the coherent sum is approximately S = π * w0^2 / (dx * dy).
            # To ensure the reconstructed field has amplitude 'amp', we must scale by 1/S.
            norm_factor = (dx * dy) / (π * w0s^2)
            E0_complex = normalize(pol_axis) * (amp * exp(im * ph) * norm_factor)
            
            b = AstigmaticGaussianBeamlet(pos, local_dir, λ, w0s; E0 = E0_complex)
            push!(beams, b)
        end
    end
    
    return AstigmaticBeamGroup(beams)
end
