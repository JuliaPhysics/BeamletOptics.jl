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

const _MAX_PERIODS = 2000

"""
    _resolve_λ_vis(λ_vis, L_plot)

Resolve the visualization wavelength: `λ_vis` if given, else `L_plot / 20`.
Clamped from below to `L_plot / _MAX_PERIODS` so that a `pol_λ` mistakenly
set to the physical ray wavelength cannot generate millions of points.
"""
function _resolve_λ_vis(λ_vis, L_plot)
    λ = something(λ_vis, L_plot / 20)
    λ_min = L_plot / _MAX_PERIODS
    if λ < λ_min
        @warn "pol_λ is too fine for the plotted path; clamping. pol_λ is a \
               visualization wavelength, not the physical ray wavelength." λ λ_min maxlog=1
        λ = λ_min
    end
    return λ
end

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
    λ_vis = _resolve_λ_vis(λ_vis, L)
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
    _polarization_points(beam::BMO.Beam{<:Any, <:BMO.PolarizedRay}; flen, λ_vis = nothing, amplitude = nothing, ppl = 32)

Field-vector curve sample points for a whole `Beam` tree (chief + child
beams), with phase continuity across segments.
"""
function _polarization_points(beam::BMO.Beam{<:Any, <:BMO.PolarizedRay}; flen, λ_vis = nothing, amplitude = nothing, ppl = 32)
    # Pass 1: Emax and the max plotted length over any branch, over the same segments as pass 2.
    Emax = 0.0
    L_max = 0.0
    for b in PreOrderDFS(beam)
        for r in BMO.rays(b)
            Emax = max(Emax, norm(_transverse(BMO.polarization(r), BMO.direction(r))))
            isnothing(BMO.intersection(r)) && break
        end
        L_max = max(L_max, length(b) + flen)
    end
    if L_max <= 0
        return Point3f[]
    end
    Emax = max(Emax, eps())
    λ_vis = _resolve_λ_vis(λ_vis, L_max)
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
    _polarization_points(agb::BMO.AstigmaticGaussianBeamlet; flen, λ_vis = nothing, scale = 1.0,
                         focus_exponent = 1.0, gain_max = 10.0, ppl = 32)

Field-vector curve sample points along the chief ray of an
`AstigmaticGaussianBeamlet` tree. The curve amplitude follows the on-axis field
amplitude of the beamlet, clamped to a maximum gain relative to the reference
amplitude,

    a(z) = scale * r_ref * min((A_ref / A(z))^(focus_exponent / 2), gain_max) / Emax,

with `A = ‖b‖·‖c‖` the product of the 1/e² semi-axes from
`(_, b, c) = waist_parameters(child, z)`, and `r_ref`, `A_ref` the mean radius
and semi-axis product at the start of the root beamlet. For `focus_exponent = 1`
this is the physical scaling `E ∝ √(w0x·w0y / (wx·wy))`, so the curve is raised
where the beam is compressed (focus) and flattened where it expands, saturating
at `gain_max` so the curve stays on-screen through a tight focus;
`focus_exponent = 0` gives a constant amplitude. Gouy phase and phase-front
curvature are ignored; the curve is a qualitative visualization only.
"""
function _polarization_points(agb::BMO.AstigmaticGaussianBeamlet; flen, λ_vis = nothing, scale = 1.0,
        focus_exponent = 1.0, gain_max = 10.0, ppl = 32)
    # Pass 1: Emax (over chief rays) and the max plotted length over any branch.
    Emax = 0.0
    L_max = 0.0
    for child in PreOrderDFS(agb)
        for ray in BMO.rays(child.c)
            Emax = max(Emax, norm(_transverse(BMO.polarization(ray), BMO.direction(ray))))
            isnothing(BMO.intersection(ray)) && break
        end
        L_max = max(L_max, length(child) + flen)
    end
    if L_max <= 0
        return Point3f[]
    end
    Emax = max(Emax, eps())
    λ_vis = _resolve_λ_vis(λ_vis, L_max)
    k_vis = 2π / λ_vis

    # Reference beam size at the start of the root beamlet
    (_, b_ref, c_ref) = BMO.waist_parameters(agb, 0.0)
    r_ref = (norm(b_ref) + norm(c_ref)) / 2
    A_ref = norm(b_ref) * norm(c_ref)

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
                A = norm(b) * norm(c)
                gain = if A_ref <= 0
                    one(A_ref)          # degenerate source: no focus scaling
                elseif A <= 0
                    gain_max            # exact caustic: saturate, never divide by zero
                else
                    (A_ref / A)^(focus_exponent / 2)
                end
                return scale * r_ref * min(gain, gain_max) / Emax
            end

            S = _field_segment!(pts, p, d, E⊥, S, n, L, amp, k_vis, λ_vis, ppl)
            l += L
            is_last && break
        end
    end
    return pts
end
