"""
    hit_sequence(beam::Beam)

The objects that the rays of the `beam` and of all its child beams hit, in depth-first order, with
`nothing` for a ray that leaves the system without an intersection. Two beams of a non-vignetted
bundle that were traced through the same system have the same hit sequence, see [`vignetted`](@ref).
"""
function hit_sequence(beam::Beam)
    seq = Any[]
    for b in PreOrderDFS(beam), ray in rays(b)
        isect = intersection(ray)
        push!(seq, isnothing(isect) ? nothing : object(isect))
    end
    return seq
end

_same_sequence(a, b) = length(a) == length(b) && all(((x, y),) -> x === y, zip(a, b))

"""
    vignetted(group::AbstractBeamGroup, reference::Beam) -> Vector{Int}

Indices of the beams of the traced `group` that are vignetted with respect to the traced
`reference` beam, e.g. the axial or chief ray of a bundle. A beam is vignetted if the sequence of
objects that it hits (including its child beams, and a ray that leaves the system without
intersection) differs from the one of the `reference`, e.g. because it misses a lens or is stopped
by an aperture. The criterion derives from the traced geometry only, no order of the objects in the
[`System`](@ref) is assumed. A beam that is not vignetted by an object but changes its path
otherwise, e.g. by total internal reflection at a steep lens rim, is also reported.

# Examples

```julia
src = UniformDiscSource([0, 0, 0], [0, 1, 0], 20e-3)
axis = Beam([0, 0, 0], [0, 1, 0])
solve_system!(system, src)
solve_system!(system, axis)
vignetted(src, axis)
```
"""
function vignetted(group::AbstractBeamGroup, reference::Beam)
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
end
