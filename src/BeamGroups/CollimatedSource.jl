"""
    CollimatedSource <: AbstractBeamGroup

Represents a parallel bundle of [`Beam`](@ref)s being emitted from a disk in space.

# Fields

- `beams`: a vector of all [`Beam`](@ref)s originating from the source
- `diameter`: the diameter of the outermost beam ring
- `center`: source position, pivot for rotations
- `orientation`: right-handed orthonormal matrix, columns are the sampling reference vector, the central source direction and their cross product, see [`AbstractBeamGroup`](@ref)

# Functions

- `diameter`: returns the diameter of the source
"""
mutable struct CollimatedSource{T, R <: AbstractRay{T}} <: AbstractBeamGroup{T, R}
    beams::Vector{Beam{T, R}}
    diameter::T
    center::Point3{T}
    orientation::SMatrix{3, 3, T, 9}
end

"""
    CollimatedSource(beams, diameter, pos, dir::AbstractVector)
    CollimatedSource(beams, diameter, pos, orientation::AbstractMatrix)

Wraps existing `beams` into a [`CollimatedSource`](@ref) with the given `diameter` and the
source position `pos`. The group orientation is either derived from the central direction
`dir` (the sampling reference vector is then picked deterministically via [`normal3d`](@ref)),
or passed explicitly as a right-handed orthonormal 3x3 `orientation` matrix whose second column
is the central direction. An invalid `orientation` throws an `ArgumentError`.
"""
function CollimatedSource(beams::Vector{Beam{T, R}}, diameter, pos, dir::AbstractVector) where {T, R <: AbstractRay{T}}
    d = normalize(dir)
    M = _group_orientation(d, _sampling_basis(d, nothing, T), T)
    _check_kinematic_members(beams)
    return CollimatedSource{T, R}(beams, T(diameter), Point3{T}(pos), M)
end

function CollimatedSource(beams::Vector{Beam{T, R}}, diameter, pos, orientation::AbstractMatrix) where {T, R <: AbstractRay{T}}
    _check_kinematic_members(beams)
    return CollimatedSource{T, R}(beams, T(diameter), Point3{T}(pos), _check_orientation(orientation, T))
end
diameter(cs::CollimatedSource) = cs.diameter

"""
    CollimatedSource(pos, dir, diameter, λ; num_rings, num_rays, basis)

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
- `basis`: Optional reference vector (e.g. `[1,0,0]`) to define the starting azimuthal angle for the beam rings.

!!! info "Reproducible sampling"
    If no `basis` is passed, the orthogonal basis vectors spanning the pupil plane are derived
    from `dir` deterministically, so two sources sharing the same `dir`, `diameter`,
    `num_rings` and `num_rays` sample exactly the same ray positions. Pass a `basis` to rotate
    the azimuthal sampling of a source about its own axis, e.g. to interleave several
    otherwise identical sources.
"""
function CollimatedSource(
        pos::AbstractArray{P},
        dir::AbstractArray{D1},
        diameter::D2,
        λ::L = 1e-6;
        num_rings::Int = 10,
        num_rays::Int = 100 * num_rings,
        basis::Union{Nothing, AbstractVector} = nothing
) where {P <: Real, D1 <: Real, D2 <: Real, L <: Real}
    T = promote_type(P, D1, D2, L)
    if num_rays < num_rings * 20
        throw(ErrorException("No. of rays should be atleast 20x no. of rings (passed: $num_rays, req: $(num_rings*20))"))
    end
    # ensure normalization
    dir = normalize(dir)
    # define buffer
    beams = Vector{Beam{T, Ray{T}}}()
    push!(beams, Beam(Ray(pos, dir, λ)))
    num_rays -= 1
    # setup concentric beam ring radii
    b1 = _sampling_basis(dir, basis, T)
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
    return CollimatedSource(beams, diameter, pos, _group_orientation(dir, b1, T))
end

"""
    UniformDiscSource(pos, dir, diameter, λ; num_rays=1_000, basis)

Generates a ray fan with *equal area per ray* across a circular pupil
using the deterministic sunflower (Fibonacci) pattern.

The radius of the `k`-th ray (`k = 0 … N-1`) is `ρₖ = diameter/2 ⋅ √((k + ½)/N)`, its azimuth is
`k` times the golden angle `π(3 - √5)`.

!!! note
    This is merely a [`CollimatedSource`](@ref) constructor which uses Fibonacci sampling
    instead of concentric rings. Unlike [`CollimatedSource`](@ref), there is no dedicated
    center beam at `pos`, i.e. all rays are equally weighted samples of the pupil.

# Arguments

The following inputs and arguments can be used to configure the underlying [`CollimatedSource`](@ref):

## Inputs

- `pos`: center of the pupil disc, i.e. the source position
- `dir`: starting direction of all beams
- `diameter`: outer beam bundle diameter in [m]
- `λ = 1e-6`: wavelength in [m]

## Keyword Arguments

- `num_rays=1000`: total number of rays in the source
- `basis`: Optional reference vector (e.g. `[1,0,0]`) to define the starting azimuthal angle of the sunflower pattern.

!!! info "Reproducible sampling"
    If no `basis` is passed, the orthogonal basis vectors spanning the pupil plane are derived
    from `dir` deterministically, so two sources sharing the same `dir`, `diameter` and
    `num_rays` sample exactly the same ray positions. Pass a `basis` to rotate the sunflower
    pattern about its own axis, e.g. to interleave several otherwise identical sources.
"""
function UniformDiscSource(
        pos::AbstractArray{P},
        dir::AbstractArray{D1},
        diameter::D2,
        λ::L = 1e-6;
        # kwargs
        num_rays::Int = 1_000,
        basis::Union{Nothing, AbstractVector} = nothing
) where {P <: Real, D1 <: Real, D2 <: Real, L <: Real}
    T = promote_type(P, D1, D2, L)
    R = diameter / 2
    beams = Vector{Beam{T, Ray{T}}}(undef, num_rays)
    dir = normalize(dir)
    # orthogonal basis in the pupil plane
    e1 = _sampling_basis(dir, basis, T)
    e2 = normal3d(dir, e1)
    for k in 0:(num_rays - 1)
        ρ = √((k + 0.5) / num_rays)     # equal-area radius
        φ = k * _GOLDEN_ANGLE
        r = R * ρ
        x = r * cos(φ) * e1 + r * sin(φ) * e2
        beams[k + 1] = Beam(pos + x, dir, λ)
    end
    return CollimatedSource(beams, diameter, pos, _group_orientation(dir, e1, T))
end
