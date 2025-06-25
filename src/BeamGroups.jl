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
    UniformDiscSource(pos, dir, diameter, λ; N=1_000)

Generates `N` rays with *equal area per ray* across a circular pupil
using the deterministic sunflower (Fibonacci) pattern.

!!! note
    This is merely a [`CollimatedSource`](@ref) constructor which uses Fibonacci sampling
    instead of a linear grid.
"""
function UniformDiscSource(pos, dir, D, λ=1e-6; N=1_000)
    T  = promote_type(eltype(pos), eltype(dir), typeof(D), typeof(λ))
    R  = D / 2
    φ0 = 2π / (1 + √5)         # golden angle
    beams = Vector{Beam{T,Ray{T}}}(undef, N)
    # orthogonal basis in the pupil plane
    e1 = normal3d(dir)
    e2 = normalize(cross(dir, e1))
    for k in 0:N-1
        ρ  = √((k + 0.5)/N)     # equal-area radius
        φ  = k * φ0
        r  = R * ρ
        x  =  r * cos(φ) * e1 + r * sin(φ) * e2
        beams[k+1] = Beam(pos + x, dir, λ)
    end
    return CollimatedSource(beams, T(D))
end
