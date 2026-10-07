"""
    CollimatedSource <: AbstractBeamGroup

Represents a parallel bundle of [`Beam`](@ref)s being emitted from a disk in space.

# Fields

- `beams`: a vector of all [`Beam`](@ref)s originating from the source. The first one is the center beam,
  which starts at the source position along the central source direction, unless the source wraps given beams.
- `diameter`: the diameter of the outermost beam ring
- `center`: source position, pivot for rotations
- `orientation`: right-handed orthonormal matrix, columns are the sampling reference vector, the central source direction and their cross product, see [`AbstractBeamGroup`](@ref)
- `sampling`: how the rays were sampled (rings, sunflower, line, or given beams), see [`set_num_rays!`](@ref)

# Functions

- `diameter`: returns the diameter of the source
- [`set_num_rays!`](@ref): regenerates the rays with another number of rays
"""
mutable struct CollimatedSource{T, R <: AbstractRay{T}} <: AbstractBeamGroup{T, R}
    beams::Vector{Beam{T, R}}
    diameter::T
    center::Point3{T}
    orientation::SMatrix{3, 3, T, 9}
    sampling::AbstractSampling
end

"""
    CollimatedSource(beams, diameter, pos, dir::AbstractVector)
    CollimatedSource(beams, diameter, pos, orientation::AbstractMatrix)

Wraps existing `beams` into a [`CollimatedSource`](@ref) with the given `diameter` and the
source position `pos`. The group orientation is either derived from the central direction
`dir` (the sampling reference vector is then picked deterministically via [`normal3d`](@ref)),
or passed explicitly as a right-handed orthonormal 3x3 `orientation` matrix whose second column
is the central direction. An invalid `orientation` throws an `ArgumentError`. The rays of such a
source can not be regenerated with [`set_num_rays!`](@ref).
"""
function CollimatedSource(beams::Vector{Beam{T, R}}, diameter, pos, dir::AbstractVector) where {T, R <: AbstractRay{T}}
    d = normalize(dir)
    return CollimatedSource(beams, diameter, pos, _group_orientation(d, sampling_basis(d, nothing, T), T))
end

CollimatedSource(beams::Vector{<:Beam}, diameter, pos, orientation::AbstractMatrix) =
    CollimatedSource(beams, diameter, pos, orientation, NoSampling())

function CollimatedSource(beams::Vector{Beam{T, R}}, diameter, pos, orientation::AbstractMatrix,
        sampling::AbstractSampling) where {T, R <: AbstractRay{T}}
    _check_kinematic_members(beams)
    return CollimatedSource{T, R}(beams, T(diameter), Point3{T}(pos), _check_orientation(orientation, T), sampling)
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
    # ensure normalization
    dir = normalize(dir)
    b1 = sampling_basis(dir, basis, T)
    sampling = DiscRings(num_rings)
    beams = source_beams(sampling, pos, dir, b1, diameter, λ, num_rays, T)
    return CollimatedSource(beams, diameter, pos, _group_orientation(dir, b1, T), sampling)
end


"""
    UniformDiscSource(pos, dir, diameter, λ; num_rays=1_000, basis)

Generates a ray fan with *equal area per ray* across a circular pupil
using the deterministic sunflower (Fibonacci) pattern.

The spiral starts with the center ray at `pos` (`k = 0`). The radius of the `k`-th ray
(`k = 1 … N-1`) is `ρₖ = diameter/2 ⋅ √((k + ½)/N)`, its azimuth is `k` times the golden angle
`π(3 - √5)`.

!!! note
    This is merely a [`CollimatedSource`](@ref) constructor which uses Fibonacci sampling
    instead of concentric rings. All rays are equally weighted samples of the pupil: ray `k`
    represents the part of the disc between the radii `diameter/2 ⋅ √(k/N)` and
    `diameter/2 ⋅ √((k + 1)/N)`. The center ray represents the innermost of these parts, a disc,
    by its center.

# Arguments

The following inputs and arguments can be used to configure the underlying [`CollimatedSource`](@ref):

## Inputs

- `pos`: center of the pupil disc, i.e. the source position
- `dir`: starting direction of all beams
- `diameter`: outer beam bundle diameter in [m]
- `λ = 1e-6`: wavelength in [m]

## Keyword Arguments

- `num_rays=1000`: total number of rays in the source, must be `≥ 1`
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
    dir = normalize(dir)
    e1 = sampling_basis(dir, basis, T)
    beams = source_beams(DiscSunflower(), pos, dir, e1, diameter, λ, num_rays, T)
    return CollimatedSource(beams, diameter, pos, _group_orientation(dir, e1, T), DiscSunflower())
end

"""
    UniformLineSource(pos, dir, width, λ; num_rays=101, basis)

Generates a plane sheet of parallel [`Beam`](@ref)s, i.e. a two-dimensional collimated source:
the rays start equidistantly on a line of the length `width` centered at `pos` and travel along
`dir`. The line lies along `basis`, such that all rays lie in the plane spanned by `dir` and
`basis`. This is useful for layout sketches and for tracing a single (e.g. meridional or
sagittal) section of a system.

The first ray is the center ray at `pos`. It is followed by `h = (N - 1) ÷ 2` rays on each side
of it, at `pos + j/h ⋅ width/2 ⋅ basis` for `j = -h … h` without `j = 0`, i.e. from one marginal
ray at the end of the line to the other. A single ray starts at `pos`.

!!! info "Use an odd number of rays"
    A pattern that is symmetric to its center ray has an odd number of rays. An even `num_rays`
    yields the pattern of `num_rays - 1` rays and a second ray at `pos`.

!!! note
    This is merely a [`CollimatedSource`](@ref) constructor which samples one diameter of the
    disc instead of concentric rings. `width` is stored as the `diameter` of the source.

!!! warning
    The rays sample a line, not an area: they are not equal-area samples of a circular pupil.
    For a point spread function or an [`intensity`](@ref) use [`UniformDiscSource`](@ref).

# Arguments

The following inputs and arguments can be used to configure the underlying [`CollimatedSource`](@ref):

## Inputs

- `pos`: center of the line, i.e. the source position
- `dir`: starting direction of all beams
- `width`: length of the line in [m], i.e. the distance between the marginal rays
- `λ = 1e-6`: wavelength in [m]

## Keyword Arguments

- `num_rays=101`: total number of rays in the source, must be `≥ 1`
- `basis`: direction of the line (e.g. `[1,0,0]`), projected into the plane normal to `dir`.
  Must not be zero or parallel to `dir`.

!!! info "Orientation of the line"
    If no `basis` is passed, the direction of the line is derived from `dir` deterministically
    via [`normal3d`](@ref), which is rarely the section of interest: pass a `basis` to choose
    the plane of the rays. The line can also be turned about `dir` afterwards with
    [`rotate3d!`](@ref).
"""
function UniformLineSource(
        pos::AbstractArray{P},
        dir::AbstractArray{D1},
        width::D2,
        λ::L = 1e-6;
        # kwargs
        num_rays::Int = 101,
        basis::Union{Nothing, AbstractVector} = nothing
) where {P <: Real, D1 <: Real, D2 <: Real, L <: Real}
    T = promote_type(P, D1, D2, L)
    dir = normalize(dir)
    e1 = sampling_basis(dir, basis, T)
    beams = source_beams(DiscLine(), pos, dir, e1, width, λ, num_rays, T)
    return CollimatedSource(beams, width, pos, _group_orientation(dir, e1, T), DiscLine())
end

