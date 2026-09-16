# =============================================================================
# Polarization field-vector curve helpers for `render!(...; show_polarization)`
#
# Ported from the standalone prototype in `scripts/polarization_field_trace.jl`
# (see that file's header for the phase/time convention). The curve is a
# static t=0 snapshot, Re{E_perp * exp(i*k_vis*s)}, evaluated along the
# accumulated optical path length s -- it is a qualitative visualization, not
# a physically accurate instantaneous field animation.
# =============================================================================

"""
    _transverse(E, d)

Projects a complex field vector `E` onto the plane orthogonal to the real,
unit ray direction `d`: `E - (d⋅E)·d`. Uses the non-conjugating product
(`sum(d .* E)`, not `LinearAlgebra.dot`, which would conjugate `E`).
"""
_transverse(E, d) = E .- sum(d .* E) .* d

"""
    _plotted_length(ray, flen)

Plotted length of `ray`: `flen` if the ray has no intersection, otherwise its
geometric length.
"""
_plotted_length(ray, flen) = isnothing(BMO.intersection(ray)) ? flen : length(ray)

"""
    _field_segment!(pts, p, d, E⊥, S, n, L, amp, k_vis, λ_vis, ppl)

Samples the field-vector curve `x(t) = p + t·d + amp(t)·Re{E⊥·exp(i·k_vis·(S + n·t))}`
for `t ∈ [0, L]` and pushes the resulting `Point3f`s onto `pts`, followed by a
`Point3f(NaN, NaN, NaN)` segment separator. `amp` is a function `t -> Real`.
Returns the accumulated optical path length `S + n*L`. Zero-length segments
(`L <= eps(float(L))`) push nothing and return `S` unchanged.
"""
function _field_segment!(pts, p, d, E⊥, S, n, L, amp, k_vis, λ_vis, ppl)
    if L <= eps(float(L))
        return S
    end
    N = max(2, ceil(Int, ppl * n * L / λ_vis) + 1)
    for j in 0:(N - 1)
        t = j * L / (N - 1)
        u = amp(t) .* real.(E⊥ .* cis(k_vis * (S + n * t)))
        push!(pts, Point3f(p .+ t .* d .+ u))
    end
    push!(pts, Point3f(NaN, NaN, NaN))
    return S + n * L
end

"""
    _render_field_curve!(axis, pts; color, linewidth)

Draws the NaN-separated point vector `pts` as a single `lines!` call. Does
nothing if `pts` is empty.
"""
function _render_field_curve!(axis, pts; color, linewidth)
    if isempty(pts)
        return nothing
    end
    lines!(axis, getindex.(pts, 1), getindex.(pts, 2), getindex.(pts, 3); color, linewidth)
    return nothing
end

"""
    _check_polarized(x)

Throws an `ArgumentError` unless `x` is a `BMO.PolarizedRay`, or a `BMO.Beam`
whose ray type is a `BMO.PolarizedRay`.
"""
function _check_polarized(x)
    ok = x isa BMO.Beam ? x isa BMO.Beam{<:Any, <:BMO.PolarizedRay} : x isa BMO.PolarizedRay
    if !ok
        throw(ArgumentError("show_polarization = true requires PolarizedRays, got $(typeof(x))"))
    end
    return nothing
end

"""
    _polarization_points(ray::BMO.PolarizedRay; flen, λ_vis = nothing, amplitude = nothing, ppl = 32)

Field-vector curve sample points for a single `PolarizedRay`, see
[`_field_segment!`](@ref).
"""
function _polarization_points(ray::BMO.PolarizedRay; flen, λ_vis = nothing, amplitude = nothing, ppl = 32)
    L = _plotted_length(ray, flen)
    if L <= 0
        return Point3f[]
    end
    λ_vis = something(λ_vis, L / 20)
    σ = something(amplitude, λ_vis / 4)
    k_vis = 2π / λ_vis

    d = BMO.direction(ray)
    E⊥ = _transverse(BMO.polarization(ray), d)
    Emax = max(norm(E⊥), eps())

    pts = Point3f[]
    _field_segment!(pts, BMO.position(ray), d, E⊥, 0.0, BMO.refractive_index(ray), L,
        t -> σ / Emax, k_vis, λ_vis, ppl)
    return pts
end

"""
    _polarization_points(beam::BMO.Beam; flen, λ_vis = nothing, amplitude = nothing, ppl = 32)

Field-vector curve sample points for a whole `Beam` tree (chief + child
beams), with phase continuity across segments.
"""
function _polarization_points(beam::BMO.Beam; flen, λ_vis = nothing, amplitude = nothing, ppl = 32)
    # Pass 1: Emax and total plotted length, over the same segments as pass 2.
    Emax = 0.0
    L_tot = 0.0
    for b in PreOrderDFS(beam)
        for r in BMO.rays(b)
            Emax = max(Emax, norm(_transverse(BMO.polarization(r), BMO.direction(r))))
            is_last = isnothing(BMO.intersection(r))
            L = is_last ? flen : length(r)
            L_tot += L
            is_last && break
        end
    end
    if L_tot <= 0
        return Point3f[]
    end
    Emax = max(Emax, eps())
    λ_vis = something(λ_vis, L_tot / 20)
    σ = something(amplitude, λ_vis / 4)
    k_vis = 2π / λ_vis

    # Pass 2: build the curve, following prototype `field_trace` exactly.
    pts = Point3f[]
    for b in PreOrderDFS(beam)
        p_beam = b.parent
        S = isnothing(p_beam) ? 0.0 : BMO.optical_path_length(p_beam)
        for r in BMO.rays(b)
            p = BMO.position(r)
            d = BMO.direction(r)
            n = BMO.refractive_index(r)
            E⊥ = _transverse(BMO.polarization(r), d)

            is_last = isnothing(BMO.intersection(r))
            L = is_last ? flen : length(r)

            S = _field_segment!(pts, p, d, E⊥, S, n, L, t -> σ / Emax, k_vis, λ_vis, ppl)
            is_last && break
        end
    end
    return pts
end

"""
    _polarization_points(agb::BMO.AstigmaticGaussianBeamlet; flen, λ_vis = nothing, scale = 1.0, ppl = 32)

Field-vector curve sample points along the chief ray of an
`AstigmaticGaussianBeamlet` tree. The curve amplitude at each sample is the
local mean 1/e² beam radius, `scale * (‖b‖ + ‖c‖)/2 / Emax`, where `(_, b, c)
= waist_parameters(child, z)`. Gouy phase and phase-front curvature are
ignored; the curve is a qualitative visualization only.
"""
function _polarization_points(agb::BMO.AstigmaticGaussianBeamlet; flen, λ_vis = nothing, scale = 1.0, ppl = 32)
    # Pass 1: Emax (over chief rays) and total plotted length.
    Emax = 0.0
    L_tot = 0.0
    for child in PreOrderDFS(agb)
        for ray in BMO.rays(child.c)
            Emax = max(Emax, norm(_transverse(BMO.polarization(ray), BMO.direction(ray))))
            is_last = isnothing(BMO.intersection(ray))
            L = is_last ? flen : length(ray)
            L_tot += L
            is_last && break
        end
    end
    if L_tot <= 0
        return Point3f[]
    end
    Emax = max(Emax, eps())
    λ_vis = something(λ_vis, L_tot / 20)
    k_vis = 2π / λ_vis

    pts = Point3f[]
    for child in PreOrderDFS(agb)
        parent_agb = child.parent
        l = isnothing(parent_agb) ? 0.0 : length(parent_agb)
        S = isnothing(parent_agb) ? 0.0 : BMO.optical_path_length(parent_agb)

        for ray in BMO.rays(child.c)
            p = BMO.position(ray)
            d = BMO.direction(ray)
            n = BMO.refractive_index(ray)
            E⊥ = _transverse(BMO.polarization(ray), d)

            is_last = isnothing(BMO.intersection(ray))
            L = is_last ? flen : length(ray)

            l0 = l
            amp = function (t)
                (_, b, c) = BMO.waist_parameters(child, l0 + t)
                return scale * (norm(b) + norm(c)) / 2 / Emax
            end

            S = _field_segment!(pts, p, d, E⊥, S, n, L, amp, k_vis, λ_vis, ppl)
            l += L
            is_last && break
        end
    end
    return pts
end
