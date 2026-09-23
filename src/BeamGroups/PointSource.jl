"""
    PointSource <: AbstractBeamGroup

Represents a cone of [`Beam`](@ref)s being emitted from a single point in space.

# Fields

- `beams`: a vector of all [`Beam`](@ref)s originating from the source
- `NA`: the [`numerical_aperture`](@ref) of the point source spread angle
- `center`: source position, pivot for rotations
- `orientation`: right-handed orthonormal matrix, columns are the sampling reference vector, the central source direction and their cross product, see [`AbstractBeamGroup`](@ref)

# Functions

- `numerical_aperture`: returns the NA of the source
"""
mutable struct PointSource{T, R <: AbstractRay{T}} <: AbstractBeamGroup{T, R}
    beams::Vector{Beam{T, R}}
    NA::T
    center::Point3{T}
    orientation::SMatrix{3, 3, T, 9}
end

"""
    PointSource(beams, NA, pos, dir::AbstractVector)
    PointSource(beams, NA, pos, orientation::AbstractMatrix)

Wraps existing `beams` into a [`PointSource`](@ref) with the numerical aperture `NA` and the
source position `pos`. The group orientation is either derived from the central direction
`dir` (the sampling reference vector is then picked deterministically via [`normal3d`](@ref)),
or passed explicitly as a right-handed orthonormal 3x3 `orientation` matrix whose second column
is the central direction. An invalid `orientation` throws an `ArgumentError`.
"""
function PointSource(beams::Vector{Beam{T, R}}, NA, pos, dir::AbstractVector) where {T, R <: AbstractRay{T}}
    d = normalize(dir)
    M = _group_orientation(d, _sampling_basis(d, nothing, T), T)
    _check_kinematic_members(beams)
    return PointSource{T, R}(beams, T(NA), Point3{T}(pos), M)
end

function PointSource(beams::Vector{Beam{T, R}}, NA, pos, orientation::AbstractMatrix) where {T, R <: AbstractRay{T}}
    _check_kinematic_members(beams)
    return PointSource{T, R}(beams, T(NA), Point3{T}(pos), _check_orientation(orientation, T))
end
numerical_aperture(ps::PointSource) = ps.NA

"""
    PointSource(pos, dir, θ, λ; num_rings, num_rays, basis)

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
- `θ`: half spread angle in rad, must be `< π`
- `λ = 1e-6`: wavelength in [m], default val. is 1000 nm

## Keyword Arguments

- `num_rings`: number of concentric beam rings, default is 10
- `num_rays`: total number of rays in the source, default is 100x num_rings
- `basis`: Optional reference vector (e.g. `[1,0,0]`) to define the starting azimuthal angle for the source rings.

!!! info "Reproducible sampling"
    If no `basis` is passed, the orthogonal basis vectors are derived from `dir` deterministically,
    so two sources sharing the same `dir`, `θ`, `num_rings` and `num_rays` sample exactly the same
    ray directions. Pass a `basis` to rotate the azimuthal sampling of a source about its own axis,
    e.g. to interleave several otherwise identical sources.
"""
function PointSource(
        pos::AbstractArray{P},
        dir::AbstractArray{D},
        θ::H,
        λ::L = 1e-6;
        num_rings::Int = 10,
        num_rays::Int = 100 * num_rings,
        basis::Union{Nothing, AbstractVector} = nothing
) where {P <: Real, D <: Real, H <: Real, L <: Real}
    T = promote_type(P, D, H, L)
    if num_rays < num_rings * 20
        throw(ErrorException("No. of rays should be atleast 20x no. of rings (passed: $num_rays, req: $(num_rings*20))"))
    end
    if θ ≥ pi
        throw(ErrorException("Point source opening half-angle θ must be < π"))
    end
    # define basis vectors
    dir = normalize(dir)
    b1 = _sampling_basis(dir, basis, T)
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
    return PointSource(beams, NA, pos, _group_orientation(dir, b1, T))
end

"""
    UniformPointSource(pos, dir, θ, λ; num_rays=1_000, basis)

Generates a cone of [`Beam`](@ref)s emitted from `pos` with *equal solid angle per ray*
across the spherical cap `0 ≤ ϑ ≤ θ` around `dir`, using the deterministic sunflower
(Fibonacci) pattern. The polar angle `ϑₖ` of the `k`-th ray (`k = 0 … N-1`) follows from
`cos ϑₖ = 1 - (k + ½)/N ⋅ (1 - cos θ)`, its azimuth is `k` times the golden angle `π(3 - √5)`.

!!! note
    This is merely a [`PointSource`](@ref) constructor which uses Fibonacci sampling
    instead of concentric rings. Unlike [`PointSource`](@ref), there is no dedicated center
    beam along `dir`, i.e. all rays are equally weighted samples of the cap.

# Arguments

The following inputs and arguments can be used to configure the underlying [`PointSource`](@ref):

## Inputs

- `pos`: starting position of all beams
- `dir`: central source direction, i.e. the cone axis
- `θ`: half spread angle in rad, must be `< π`
- `λ = 1e-6`: wavelength in [m], default val. is 1000 nm

## Keyword Arguments

- `num_rays=1000`: total number of rays in the source, must be `≥ 1`
- `basis`: Optional reference vector (e.g. `[1,0,0]`) to define the starting azimuthal angle of the sunflower pattern.
  Must not be zero or parallel to `dir`.

!!! info "Reproducible sampling"
    If no `basis` is passed, the orthogonal basis vectors are derived from `dir` deterministically,
    so two sources sharing the same `dir`, `θ` and `num_rays` sample exactly the same ray
    directions. Pass a `basis` to rotate the sunflower pattern about its own axis, e.g. to
    interleave several otherwise identical sources.
"""
function UniformPointSource(
        pos::AbstractArray{P},
        dir::AbstractArray{D},
        θ::H,
        λ::L = 1e-6;
        num_rays::Int = 1_000,
        basis::Union{Nothing, AbstractVector} = nothing
) where {P <: Real, D <: Real, H <: Real, L <: Real}
    T = promote_type(P, D, H, L)
    if θ ≥ pi
        throw(ErrorException("Point source opening half-angle θ must be < π"))
    end
    if num_rays < 1
        throw(ErrorException("No. of rays must be at least 1 (passed: $num_rays)"))
    end
    beams = Vector{Beam{T, Ray{T}}}(undef, num_rays)
    dir = normalize(dir)
    # orthogonal basis normal to the cone axis
    e1 = _sampling_basis(dir, basis, T)
    e2 = normal3d(dir, e1)
    one_minus_cosθ = 1 - cos(θ)
    for k in 0:(num_rays - 1)
        cosϑ = 1 - (k + 0.5) / num_rays * one_minus_cosθ    # equal solid angle
        sinϑ = sqrt(max(0, (1 - cosϑ) * (1 + cosϑ)))
        φ = k * _GOLDEN_ANGLE
        cdir = cosϑ * dir + sinϑ * (cos(φ) * e1 + sin(φ) * e2)
        beams[k + 1] = Beam(pos, cdir, λ)
    end
    NA = numerical_aperture(θ)
    return PointSource(beams, NA, pos, _group_orientation(dir, e1, T))
end