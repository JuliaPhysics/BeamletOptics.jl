"""
    surface_part(shape::AbstractShape, point)

The part of `shape` that owns the boundary at `point`: `shape` itself, or for an
[`AbstractCompositeSDF`](@ref) the (recursively resolved) operand whose surface is closest to `point`.
This distinguishes e.g. the optical faces of a lens from its mechanical rim, which are all
reported by the same [`Intersection`](@ref) `shape`.
"""
surface_part(shape::AbstractShape, point) = shape

function surface_part(shape::AbstractCompositeSDF, point)
    ops = operands(shape)
    return surface_part(ops[argmin(abs(sdf(op, point)) for op in ops)], point)
end

# the traversed surface of a ray: its object and the boundary part of the object that it hits
function _hit_surface(ray::AbstractRay)
    isect = intersection(ray)
    isnothing(isect) && return nothing
    s = shape(isect)
    isnothing(s) && return (object(isect), nothing)
    return (object(isect), surface_part(s, position(ray) + length(isect) * direction(ray)))
end

# the chief beam, i.e. the scalar-ray path, of a beam or beamlet
_chief_beam(b::Beam) = b
_chief_beam(b::GaussianBeamlet) = b.chief
_chief_beam(b::AstigmaticGaussianBeamlet) = b.c

"""
    hit_sequence(beam::AbstractBeam)

The surfaces that the chief rays of the `beam` and of all its child beams hit, in depth-first
order, as tuples of the hit object and the boundary part of the object (see [`surface_part`](@ref))
with `nothing` for a ray that leaves the system without an intersection. Two beams of a
non-vignetted bundle that were traced through the same system have the same hit sequence, see
[`vignetted`](@ref). For a [`GaussianBeamlet`](@ref) or [`AstigmaticGaussianBeamlet`](@ref) the
`chief` ray is traversed.
"""
function hit_sequence(beam::AbstractBeam)
    seq = Any[]
    for b in PreOrderDFS(beam), ray in rays(_chief_beam(b))
        push!(seq, _hit_surface(ray))
    end
    return seq
end

_same_sequence(a, b) = length(a) == length(b) && all(((x, y),) -> x === y, zip(a, b))

"""
    vignetted(group::AbstractBeamGroup, reference::AbstractBeam) -> Vector{Int}

Indices of the beams of the traced `group` that are vignetted with respect to the traced
`reference` beam, e.g. the axial or chief ray of a bundle. A beam is vignetted if the sequence of
surfaces that it hits (object and boundary part of the object, including its child beams, and a
ray that leaves the system without intersection) differs from the one of the `reference`, e.g. because it misses a lens or is stopped
by an aperture. The criterion derives from the traced geometry only, no order of the objects in the
[`System`](@ref) is assumed. A beam that is not vignetted by an object but changes its path
otherwise, e.g. by total internal reflection at a steep lens rim, is also reported. The `group` and the `reference` can consist of [`Beam`](@ref)s,
[`GaussianBeamlet`](@ref)s or [`AstigmaticGaussianBeamlet`](@ref)s, in which case the chief rays are
compared.

# Examples

```julia
src = UniformDiscSource([0, 0, 0], [0, 1, 0], 20e-3)
axis = Beam([0, 0, 0], [0, 1, 0])
solve_system!(system, src)
solve_system!(system, axis)
vignetted(src, axis)
```
"""
function vignetted(group::AbstractBeamGroup, reference::AbstractBeam)
    ref = hit_sequence(reference)
    return findall(b -> !_same_sequence(hit_sequence(b), ref), beams(group))
end

"""
    clear_aperture(system, pos, dir; λ = 1e-6, rings = 6, azimuths = 24, rtol = 1e-3,
        d_start = 1e-3, d_max = 1.0, kwargs...) -> Float64

The largest diameter in [m] of a collimated bundle of rays that passes the `system` without
vignetting. The bundle is centered on the axis through `pos` along `dir` and its diameter is
measured in the plane through `pos` normal to `dir`.

A ray is vignetted if it hits another sequence of objects than the axial ray (also including the
objects hit by child beams of beam splitting) or leaves the system early, see [`vignetted`](@ref). The
solver finds the sequence by tracing, so no order of the objects in the `system` is assumed. The
largest clear diameter is found by doubling the diameter `d_start` until the bundle is vignetted
and bisecting between the last clear and the first vignetted diameter until they differ by less
than the relative tolerance `rtol`. The returned diameter is the last diameter that was traced
without vignetting. If the bundle is not vignetted at the diameter `d_max`, `d_max` is returned,
and `0.0` if already the smallest traced bundle is vignetted. Throws an `ArgumentError` if the axial ray
does not hit any object.

# Keyword Arguments

- `λ = 1e-6`: wavelength in [m]
- `rings = 6`: number of rings of the bundle with equidistant radii; the outermost ring has the
  radius of the tested diameter
- `azimuths = 24`: number of rays per ring, equidistant in the azimuth with a staggered start angle
  from ring to ring
- `rtol = 1e-3`: relative tolerance of the diameter
- `d_start = 1e-3`, `d_max = 1.0`: first diameter and upper limit of the search in [m]
- `kwargs...`: passed on to [`solve_system!`](@ref), e.g. `depth_max`

!!! warning "Limitations"
    The result is sampling-based: vignetting that only affects the bundle between the sampled radii
    and azimuths (e.g. a narrow obstruction) is not seen, increase `rings` and `azimuths` to
    check. The search assumes that the clear region is a disc around the axis, i.e. that if a
    bundle is clear then every smaller bundle is as well; a central obstruction violates this. The
    input is collimated: for a diverging or converging field use the chief ray of a
    [`PointSource`](@ref) and [`vignetted`](@ref). The rays of the bundle are traced
    as [`Beam`](@ref)s of scalar [`Ray`](@ref)s, no Gaussian beam clipping is considered.

The probe traces do not change the state of the `system`: the hits stored in its
[`Detector`](@ref)s (also inside object groups) are saved before and restored afterwards, also if
an error is thrown.

# Examples

```julia
D = clear_aperture(system, [0, -0.1, 0], [0, 1, 0]; λ = 633e-9)
```
"""
function clear_aperture(system::AbstractSystem, pos::AbstractVector, dir::AbstractVector;
        λ::Real = 1e-6, rings::Int = 6, azimuths::Int = 24, rtol::Real = 1e-3,
        d_start::Real = 1e-3, d_max::Real = 1.0, kwargs...)
    rings >= 1 || throw(ArgumentError("rings must be at least 1"))
    azimuths >= 3 || throw(ArgumentError("azimuths must be at least 3"))
    d_start > 0 || throw(ArgumentError("d_start must be positive"))
    d_max >= d_start || throw(ArgumentError("d_max must not be smaller than d_start"))
    T = float(promote_type(eltype(pos), eltype(dir), typeof(λ), typeof(d_start), typeof(d_max)))
    d = normalize(T.(dir))
    e1 = sampling_basis(d, nothing, T)
    e2 = normal3d(d, e1)

    # the probe traces must not leave data in stateful components, e.g. detector hits
    detectors = Detector[o for o in objects(system) if o isa Detector]
    saved = map(hits, detectors)
    try
        foreach(empty!, detectors)
        axis = Beam(pos, d, λ)
        solve_system!(system, axis; kwargs...)
        ref = hit_sequence(axis)
        all(isnothing, ref) && throw(ArgumentError("the axial ray does not hit any object of the system"))

        function isclear(D)
            bundle = Beam{T, Ray{T}}[]
            for k in 1:rings, j in 1:azimuths
                φ = 2π * (j - 1 + (isodd(k) ? 0 : 0.5)) / azimuths
                r = D / 2 * k / rings
                push!(bundle, Beam(Ray(pos + r * (cos(φ) * e1 + sin(φ) * e2), d, λ)))
            end
            src = CollimatedSource(bundle, D, pos, d)
            solve_system!(system, src; progress = false, kwargs...)
            return all(b -> _same_sequence(hit_sequence(b), ref), beams(src))
        end

        # bracket the limit, then bisect
        lo, hi = zero(T), T(d_start)
        while isclear(hi)
            lo = hi
            hi >= d_max && return T(d_max)
            hi = min(2hi, T(d_max))
        end
        for _ in 1:64
            hi - lo <= rtol * hi && break
            mid = (lo + hi) / 2
            isclear(mid) ? (lo = mid) : (hi = mid)
        end
        return lo
    finally
        foreach(hits!, detectors, saved)
    end
end
