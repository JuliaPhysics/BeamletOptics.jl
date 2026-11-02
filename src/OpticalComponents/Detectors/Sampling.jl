"""
    AbstractDetectorSampling

Abstract supertype of the samplings of a [`Detector`](@ref): a sampling decides where on the
detector the readouts [`electric_field`](@ref), [`intensity`](@ref) and [`optical_power`](@ref)
are evaluated, in which window, and with which weights the power is integrated.

# Implementation reqs.

## Functions

- `sample_points(pd::Detector, sampling::S) -> DetectorSamples`: the points, weights and coordinates
  of the sampling on `pd`, see [`sample_points`](@ref BeamletOptics.sample_points). The window
  is the business of the sampling: it may look at the hits of `pd`, e.g. as
  [`RegularGrid`](@ref) does, or take a fixed one, as [`RectilinearGrid`](@ref) does. All
  readouts then work for the new sampling.

Concrete subtypes include:

- [`RegularGrid`](@ref) — `num_points×num_points` points in a window around the hits
- [`RectilinearGrid`](@ref) — the grid `xs × zs` of one's own, with any spacing
"""
abstract type AbstractDetectorSampling end

"""
    DetectorSamples(coordinates, points, weights)

The points at which the readouts of a [`Detector`](@ref) are evaluated, as returned by
[`sample_points`](@ref BeamletOptics.sample_points):

- `points`: an array of local detector coordinates `(x, z)` in \\[m\\], each a vector
  (e.g. a `Point2`). The readouts return their values in the shape of `points`.
- `weights`: the detector area in \\[m²\\] each point stands for, in the shape of `points`, as
  an array or a lazy broadcast (`Broadcast.broadcasted`), which is only evaluated by
  [`optical_power`](@ref). It is `Σ wᵢ Iᵢ` with the irradiance `Iᵢ` at the points, so the
  weights are the quadrature rule of the sampling.
- `coordinates`: a tuple that the readouts return in front of the values, e.g. `(xs, zs)` for
  a grid.

The local coordinates `(x, z)` form a left-handed system with the detector normal, see
[`electric_field`](@ref)`(pd, sampling)`.
"""
struct DetectorSamples{C <: Tuple, P <: AbstractArray, W}
    coordinates::C
    points::P
    # An array or a lazy broadcast
    weights::W
end

"""
    sample_points(pd::Detector, sampling::AbstractDetectorSampling) -> DetectorSamples

The points at which the readouts of the [`Detector`](@ref) `pd` are evaluated for the sampling
`sampling`, with their weights and the coordinates that the readouts return, see
[`DetectorSamples`](@ref BeamletOptics.DetectorSamples). Every
[`AbstractDetectorSampling`](@ref BeamletOptics.AbstractDetectorSampling) implements a method;
the fallback throws an error.
"""
sample_points(::Detector, sampling::AbstractDetectorSampling) =
    error(lazy"sample_points not defined for $(typeof(sampling))")

"""
    RectilinearGrid(xs::AbstractVector, zs::AbstractVector)

A [`Detector`](@ref) sampling on the grid `xs × zs` of local detector coordinates in \\[m\\],
with any spacing: the readouts return `(xs, zs, values)` with `values[i, j]` at
`(xs[i], zs[j])`. [`optical_power`](@ref) integrates with the trapezoidal rule along both axes.
The window is `xs × zs` as given, independent of the hits. Both axes must be non-empty and
ascending, otherwise an `ArgumentError` is thrown: a descending axis would integrate to a
negative power.
"""
struct RectilinearGrid{X <: AbstractVector{<:Real}, Z <: AbstractVector{<:Real}} <:
       AbstractDetectorSampling
    xs::X
    zs::Z
    function RectilinearGrid(xs::X, zs::Z) where {X <: AbstractVector{<:Real},
            Z <: AbstractVector{<:Real}}
        (isempty(xs) || isempty(zs)) &&
            throw(ArgumentError("RectilinearGrid needs at least one point per axis"))
        (issorted(xs) && issorted(zs)) ||
            throw(ArgumentError("RectilinearGrid needs ascending xs and zs"))
        return new{X, Z}(xs, zs)
    end
end

function sample_points(::Detector, g::RectilinearGrid)
    points = [Point2(x, z) for x in g.xs, z in g.zs]
    # Lazy outer product of weights
    weights = Broadcast.instantiate(Broadcast.broadcasted(*,
        trapezoid_weights(g.xs), transpose(trapezoid_weights(g.zs))))
    return DetectorSamples((g.xs, g.zs), points, weights)
end

# Weights of the trapezoidal rule on the nodes `xs`, zero for a single node
function trapezoid_weights(xs::AbstractVector)
    w = zeros(float(eltype(xs)), length(xs))
    for i in firstindex(xs):(lastindex(xs) - 1)
        h = (xs[i + 1] - xs[i]) / 2
        w[i - firstindex(xs) + 1] += h
        w[i - firstindex(xs) + 2] += h
    end
    return w
end

"""
    RegularGrid(; num_points = 100, center = Centroid(), crop_factor, x_min, x_max, z_min,
        z_max, x0_shift = 0, z0_shift = 0)

The default sampling of a [`Detector`](@ref): `num_points×num_points` equidistant points in a
window around its hits, in local detector coordinates. The readouts return `(xs, zs, values)`
with `values[i, j]` at `(xs[i], zs[j])`; [`optical_power`](@ref) integrates with the
trapezoidal rule. The window follows from the hits when the readout is evaluated, see
[`calc_local_lims`](@ref). A detector that no beam hit has the window of its active surface.

# Keyword arguments

- `num_points::Int = 100`: number of points per axis.
- `center::AbstractCenterAlgorithm = Centroid()`: for ray hits, how the window is centred:
  `Centroid()` uses the projection-weighted centroid, `MinMax()` the mid-point of the bounding
  box.
- `crop_factor::Real`: scales the width of the window around the hits; values >1 expand, <1
  shrink. The default is `1` for ray hits (the box around the hit points) and `3` for beamlet
  hits (three beam radii around each hit, which holds the power of a Gaussian beam except for
  a part of `10⁻⁸`).
- `x_min, x_max, z_min, z_max`: bounds of the window in \\[m\\]. A given bound replaces its side
  of the window around the hits, the others stay automatic. If a lower bound ends up above its
  upper one, the readout throws an `ArgumentError`.
- `x0_shift = 0, z0_shift = 0`: shift the whole window along `x` or `z` in \\[m\\].
"""
struct RegularGrid{C <: AbstractCenterAlgorithm, F <: Nullable{Real}, XL <: Nullable{Real},
    XH <: Nullable{Real}, ZL <: Nullable{Real}, ZH <: Nullable{Real}, S <: Real} <:
       AbstractDetectorSampling
    num_points::Int
    crop_factor::F
    center::C
    x_min::XL
    x_max::XH
    z_min::ZL
    z_max::ZH
    x0_shift::S
    z0_shift::S
end

function RegularGrid(;
        # grid
        num_points::Integer = 100,
        center::AbstractCenterAlgorithm = Centroid(),
        # window
        crop_factor::Nullable{Real} = nothing,
        x_min::Nullable{Real} = nothing,
        x_max::Nullable{Real} = nothing,
        z_min::Nullable{Real} = nothing,
        z_max::Nullable{Real} = nothing,
        # shift
        x0_shift::Real = 0,
        z0_shift::Real = 0
    )
    num_points ≥ 1 || throw(ArgumentError(
        "RegularGrid needs at least one point per axis, got num_points = $num_points"))
    return RegularGrid(Int(num_points), crop_factor, center, x_min, x_max, z_min, z_max,
        promote(x0_shift, z0_shift)...)
end

function sample_points(pd::Detector, g::RegularGrid)
    x_lo, x_hi, z_lo, z_hi = calc_local_lims(pd, hits(pd), g.crop_factor;
        center = g.center)
    x_lo, x_hi = something(g.x_min, x_lo), something(g.x_max, x_hi)
    z_lo, z_hi = something(g.z_min, z_lo), something(g.z_max, z_hi)
    x_lo ≤ x_hi || throw(ArgumentError("RegularGrid: x_min = $x_lo lies above x_max = $x_hi"))
    z_lo ≤ z_hi || throw(ArgumentError("RegularGrid: z_min = $z_lo lies above z_max = $z_hi"))
    xs = LinRange(x_lo, x_hi, g.num_points) .+ g.x0_shift
    zs = LinRange(z_lo, z_hi, g.num_points) .+ g.z0_shift
    return sample_points(pd, RectilinearGrid(xs, zs))
end
