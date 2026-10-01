# Tests of BeamletOpticsOpticsBaseExt: BMO detector hits -> OpticsBase.PlaneField and
# back via WavefrontBeamletDecomposition. Evaluated with `using OpticsBase` loaded (see
# runtests.jl), which triggers the package extension.
module TestOpticsBaseExt

using BeamletOptics
using OpticsBase
using Test
using LinearAlgebra
using Logging: with_logger, NullLogger

const BMO = BeamletOptics
const Z0 = OpticsBase.VACUUM_IMPEDANCE

const λ = 1.064e-6
const w0 = 0.5e-3
const P0 = 1e-3
const zdet = 0.2

# Detector at distance `z` along +y, optionally rotated by `tilt` about BMO's z axis.
function detector_at(z; tilt = 0.0)
    det = BMO.Detector(0.05)
    tilt == 0 || BMO.zrotate3d!(det, tilt)
    BMO.translate3d!(det, [0, z, 0])
    return det
end

function trace!(det, beams...)
    system = BMO.System([det])
    for b in beams
        BMO.solve_system!(system, b)
    end
    return det
end

agb_at(; kw...) = BMO.AstigmaticGaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0; P0,
    support = [1.0, 0, 0], kw...)

@testset "PlaneField: collimated Gaussian beam at (near) its waist" begin
    det = trace!(detector_at(zdet), agb_at())
    f = PlaneField(det; size = (96, 96), spacing = (w0 / 10, w0 / 10), progress = false)

    @test f isa OpticsBase.PlaneField
    @test f.n == 1
    @test isinf(f.R)

    # Total power reproduces the source power
    @test power(f) ≈ P0 rtol = 1e-6

    # All power travels forward; essentially nothing backward
    @test power(backward(f)) ≈ 0 atol = 1e-9 * P0
    @test power(forward(f)) ≈ power(f) rtol = 1e-9

    # H is consistent with E for a beam close to normal incidence: it matches the local
    # plane wave rule built from the *total* E to good accuracy (the two constructions
    # only differ once individual beamlet directions differ noticeably from `n`, see the
    # tilted testset below).
    Hu_para = -(f.n / Z0) .* f.E[:, :, 2]
    Hv_para = (f.n / Z0) .* f.E[:, :, 1]
    scale = maximum(abs, f.H)
    @test maximum(abs.(f.H[:, :, 1] .- Hu_para)) <= 1e-6 * scale
    @test maximum(abs.(f.H[:, :, 2] .- Hv_para)) <= 1e-6 * scale

    # The E profile matches the analytic Gaussian at this z (BMO's own closed-form
    # paraxial beamlet field, so this is essentially exact).
    zR = π * w0^2 / λ
    w_z = w0 * sqrt(1 + (zdet / zR)^2)
    Epeak2 = 4 * Z0 * P0 / (π * w_z^2)
    us = collect(OpticsBase.coordinates(f, 1))
    vs = collect(OpticsBase.coordinates(f, 2))
    n0 = length(us) ÷ 2 + 1
    I(i, j) = abs2(f.E[i, j, 1]) + abs2(f.E[i, j, 2])
    @test I(n0, n0) ≈ Epeak2 rtol = 1e-5
    for k in (5, 15, 25)
        r2 = us[n0 + k]^2 + vs[n0]^2
        @test I(n0 + k, n0) ≈ Epeak2 * exp(-2 * r2 / w_z^2) rtol = 1e-3
    end
end

@testset "PlaneField: stigmatic GaussianBeamlet (scalar, placed along u)" begin
    gb = BMO.GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0; P0, support = [1.0, 0, 0])
    det = trace!(detector_at(zdet), gb)
    f = PlaneField(det; size = (64, 64), spacing = (w0 / 8, w0 / 8), progress = false)
    @test power(f) ≈ P0 rtol = 1e-6
    # scalar convention: everything along u, nothing along v
    @test maximum(abs, f.E[:, :, 2]) == 0
end

@testset "PlaneField: phase is not conjugated (exp(−iωt), Gouy phase −atan)" begin
    # Far from the waist, so Gouy phase and curvature are large enough to tell a
    # conjugated field apart. The sampled field must equal BMO's own beamlet field.
    z = 5 * π * w0^2 / λ
    gb = BMO.GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0; P0, support = [1.0, 0, 0])
    det = trace!(detector_at(z), gb)
    f = PlaneField(det; size = (65, 65), spacing = (w0, w0), progress = false)
    c = 33
    r = 4w0
    @test f.E[c, c, 1] ≈ BMO.electric_field(gb, 0.0, z) rtol = 1e-6
    @test f.E[c + 4, c, 1] ≈ BMO.electric_field(gb, r, z) rtol = 1e-6
    # analytic exp(−iωt) Gaussian: phase kz − atan(z/zR) on axis
    k, zR = 2π / λ, π * w0^2 / λ
    @test angle(f.E[c, c, 1] * cis(-k * z)) ≈ -atan(z / zR) atol = 1e-6

    # A complex E0 of an astigmatic beamlet rotates the field by its phase once
    E_real = [1.0, 0, 0] * BMO.electric_field(2P0 / (π * w0^2))
    sample(E0) = PlaneField(
        trace!(detector_at(zdet), agb_at(; E0)); size = (16, 16), spacing = (w0 / 4, w0 / 4),
        progress = false)
    f1, f2 = sample(E_real), sample(cis(1.0) * E_real)
    @test f2.E ≈ cis(1.0) .* f1.E
    @test f2.H ≈ cis(1.0) .* f1.H
end

@testset "PlaneField: tilted beam" begin
    tilt = deg2rad(10)
    det = trace!(detector_at(zdet; tilt), agb_at())
    f = PlaneField(det; size = (96, 96), spacing = (w0 / 8, w0 / 8), progress = false)

    # Energy conservation: the whole beam crosses the (effectively infinite) tilted plane,
    # so the net Poynting flux is P0. The flux density drops by cos(tilt), the footprint
    # grows by 1/cos(tilt). This checks that E and H are the physical fields, without the
    # √cos projection factor of BMO's own detector sum.
    @test power(f) ≈ P0 rtol = 1e-4

    # H built per-hit from `d_i x E_i` (this extension's construction) differs measurably
    # from the paraxial single-direction construction `H = (n/Z0) n x E_total` once the
    # beamlet direction is not ~parallel to the plane normal `n`; the two would coincide
    # only at normal incidence.
    Hu_para = -(f.n / Z0) .* f.E[:, :, 2]
    Hv_para = (f.n / Z0) .* f.E[:, :, 1]
    scale = maximum(abs, f.H)
    rel_diff = max(maximum(abs.(f.H[:, :, 1] .- Hu_para)),
        maximum(abs.(f.H[:, :, 2] .- Hv_para))) / scale
    @test rel_diff > 1e-3
end

@testset "PlaneField: default frame and origin" begin
    det = trace!(detector_at(zdet), agb_at())
    f = PlaneField(det; size = (8, 8), spacing = (w0 / 8, w0 / 8), progress = false)
    R = BMO.orientation(det)
    @test f.origin ≈ BMO.position(det)
    @test f.axes[:, 3] ≈ -R[:, 2]              # n points downstream
    @test f.axes[:, 1] ≈ -R[:, 1]               # u = BMO local x
    @test f.axes[:, 2] ≈ cross(f.axes[:, 3], f.axes[:, 1])
    @test LinearAlgebra.det(f.axes) ≈ 1.0 atol=1e-12
end

@testset "PlaneField: errors" begin
    @test_throws ArgumentError PlaneField(detector_at(zdet); size = (8, 8), spacing = (1e-5, 1e-5))

    det = trace!(detector_at(zdet), BMO.Beam(BMO.Ray([0.0, 0, 0], [0.0, 1.0, 0.0], λ)))
    @test_throws ArgumentError PlaneField(det; size = (8, 8), spacing = (1e-5, 1e-5))

    det2 = trace!(detector_at(zdet), agb_at(), BMO.AstigmaticGaussianBeamlet(
        [1e-3, 0, 0], [0.0, 1, 0], 2λ, w0; P0, support = [1.0, 0, 0]))
    @test_throws ArgumentError PlaneField(det2; size = (8, 8), spacing = (1e-5, 1e-5))
end

@testset "WavefrontBeamletDecomposition(::PlaneField): round trip" begin
    det = trace!(detector_at(zdet), agb_at())
    f = PlaneField(det; size = (128, 128), spacing = (w0 / 16, w0 / 16), progress = false)

    group = BMO.WavefrontBeamletDecomposition(f)
    @test group isa BMO.AstigmaticBeamGroup
    @test length(BMO.beams(group)) > 0

    # Trace the decomposed field and, as the reference, the original beam to the same
    # plane further downstream
    z2 = zdet + 0.05
    sample(det) = PlaneField(det; size = (128, 128), spacing = (w0 / 16, w0 / 16),
        progress = false)
    f2 = sample(trace!(detector_at(z2), group))
    ref = sample(trace!(detector_at(z2), agb_at()))

    # Power is reproduced up to the beamlet-grid discretization error of the
    # decomposition (about 0.5% at this sampling)
    @test power(f2) ≈ power(ref) rtol = 0.01

    # Complex field, including the absolute phase: normalized correlation with the
    # directly traced beam
    c = dot(ref.E, f2.E) / (norm(ref.E) * norm(f2.E))
    @test abs(c) > 0.999
    @test abs(angle(c)) < 0.02

    @test_throws ArgumentError BMO.WavefrontBeamletDecomposition(f; basis = ([1.0,0,0],[0.0,0,1.0]))
end

@testset "GaussianModeDecomposition(::PlaneField): analytic astigmatic mode" begin
    # Plane at the origin with the default detector frame (n = +y); detector one
    # wavelength behind it, so the propagation phase is 2π and diffraction negligible.
    det = detector_at(λ)
    axes = Base.get_extension(BMO, :BeamletOpticsOpticsBaseExt)._default_plane_axes(det)
    N, Δ = 128, 3e-6
    k = 2π / λ
    w1, w2, R1, R2 = 40e-6, 60e-6, 0.1, -0.2        # radii and curvatures along ξ1, ξ2
    θ = deg2rad(30)                                  # principal axes rotated against (u, v)
    uc, vc, tilt = 30e-6, -20e-6, 2e-3               # centroid and tilt along u
    jones = normalize([1.0, 0.5im])
    amp = 100 * cis(0.7)
    cs = ((0:(N - 1)) .- N ÷ 2) .* Δ
    ψ = [begin
             X, Y = x - uc, y - vc
             ξ1, ξ2 = cos(θ) * X + sin(θ) * Y, -sin(θ) * X + cos(θ) * Y
             amp * exp(-ξ1^2 / w1^2 - ξ2^2 / w2^2) *
             cis(k * (ξ1^2 / R1 + ξ2^2 / R2) / 2 + k * tilt * X)
         end
         for x in cs, y in cs]
    E = cat(jones[1] .* ψ, jones[2] .* ψ; dims = 3)
    f = OpticsBase.PlaneField(E, (Δ, Δ), zeros(3), axes, λ)

    beamlet = BMO.GaussianModeDecomposition(f)
    @test BMO.optical_power(beamlet) / power(f) ≈ 1 rtol = 1e-3

    g = PlaneField(trace!(det, beamlet); size = (N, N), spacing = (Δ, Δ), progress = false)
    c = dot(f.E, g.E) / (norm(f.E) * norm(g.E))
    @test abs(c) > 0.9999
    @test abs(angle(c)) < 1e-3
    @test power(g) ≈ power(f) rtol = 1e-3
end

@testset "GaussianModeDecomposition(::PlaneField): propagates like the original beam" begin
    # A beamlet far from its waist (curved wavefront) is sampled, fitted and traced on;
    # the result must match the original beamlet traced directly to the same plane.
    z1, z2 = 0.5, 0.8
    sample(det) = PlaneField(det; size = (96, 96), spacing = (w0 / 12, w0 / 12),
        progress = false)
    f = sample(trace!(detector_at(z1), agb_at()))
    beamlet = BMO.GaussianModeDecomposition(f)
    @test BMO.optical_power(beamlet) / power(f) ≈ 1 rtol = 1e-4

    fit = sample(trace!(detector_at(z2), beamlet))
    ref = sample(trace!(detector_at(z2), agb_at()))
    c = dot(ref.E, fit.E) / (norm(ref.E) * norm(fit.E))
    @test abs(c) > 0.9999
    @test abs(angle(c)) < 1e-3
    @test power(fit) ≈ power(ref) rtol = 1e-4

    # only planes in air: BMO starts beams in vacuum
    fn = OpticsBase.PlaneField(f.E, f.H, f.spacing, f.origin, f.axes, f.λ; n = 1.5)
    @test_throws ArgumentError BMO.GaussianModeDecomposition(fn)
end

@testset "Oblique fields keep their tangential components" begin
    # A plane wave tilted by θ = 10° against the plane normal, sampled finely enough to
    # resolve its phase ramp. Every beamlet must reproduce the sampled tangential field
    # and be transverse to its own direction; a projection onto the transverse plane
    # would scale the tangential field by cos²θ.
    θ = deg2rad(10)
    k = 2π / λ
    e1, e2, n = [1.0, 0, 0], [0.0, 0, -1], [0.0, 1, 0]
    xs = collect(((0:15) .- 8) .* 0.5e-6)
    Eu = [complex(cis(k * sin(θ) * x)) for x in xs, y in xs]
    Ev = 0.3im .* Eu
    overlap = 1.2
    group = BMO.WavefrontBeamletDecomposition(xs, xs, Eu, Ev, n, λ; basis = (e1, e2),
        overlap)
    nf = 1 / (π * overlap^2)                     # amplitude normalization per beamlet
    for b in BMO.beams(group)
        r = BMO.first_ray(b)
        E0 = BMO.polarization(b)
        i = argmin(abs.(xs .- dot(BMO.position(r), e1)))
        j = argmin(abs.(xs .- dot(BMO.position(r), e2)))
        @test dot(e1, E0) ≈ Eu[i, j] * nf rtol = 1e-9
        @test dot(e2, E0) ≈ Ev[i, j] * nf rtol = 1e-9
        @test abs(dot(BMO.direction(r), E0)) < 1e-12 * norm(E0)
        @test dot(BMO.direction(r), n) ≈ cos(θ) rtol = 1e-3
    end

    # the Gaussian mode fit keeps the tangential field at the centroid as well
    # (wide enough window: a truncated Gaussian biases the moments)
    xg = collect(((0:63) .- 32) .* 0.5e-6)
    ψ = [exp(-(x^2 + y^2) / (2e-6)^2) * cis(k * sin(θ) * x) for x in xg, y in xg]
    c = 33                                       # index of x = 0
    f = OpticsBase.PlaneField(cat(ψ, 0.3im .* ψ; dims = 3), (0.5e-6, 0.5e-6), zeros(3),
        hcat(e1, e2, n), λ)
    m = BMO.GaussianModeDecomposition(f)
    E0 = BMO.polarization(m)
    @test dot(e1, E0) ≈ ψ[c, c] rtol = 1e-3
    @test dot(e2, E0) ≈ 0.3im * ψ[c, c] rtol = 1e-3
    @test abs(dot(BMO.direction(BMO.first_ray(m)), E0)) < 1e-12 * norm(E0)
end

@testset "PlaneField: stigmatic beamlet on a tilted plane is a transverse wave" begin
    # The scalar beamlet is polarized along u made transverse to the beam, here in the
    # plane of incidence (TM). For a transverse wave the tangential admittance is then
    # Hv / Eu = Y / cos θ (with E along u itself it would be Y cos θ).
    tilt = deg2rad(10)
    gb = BMO.GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0; P0, support = [1.0, 0, 0])
    f = PlaneField(trace!(detector_at(zdet; tilt), gb); size = (17, 17),
        spacing = (w0 / 4, w0 / 4), progress = false)
    c = 9
    @test f.H[c, c, 2] / f.E[c, c, 1] ≈ (1 / Z0) / cos(tilt) rtol = 1e-9
    @test power(f) ≈ P0 rtol = 2e-2              # coarse sampling of the footprint
end

@testset "WavefrontBeamletDecomposition: phase gradient independent of the dominant component" begin
    # Eu = A e^{iφ}, Ev = i B e^{iφ}: one flat wavefront whose polarization changes from u
    # to v across x = 0, where |Eu| = |Ev|. Every beamlet must travel along `dir`.
    xs = collect(range(-200e-6, 200e-6; length = 41))
    A = [exp(-((x + 50e-6)^2 + y^2) / (80e-6)^2) for x in xs, y in xs]
    B = [exp(-((x - 50e-6)^2 + y^2) / (80e-6)^2) for x in xs, y in xs]
    φ = cis(0.3)
    group = BMO.WavefrontBeamletDecomposition(xs, xs, complex.(A) .* φ, im .* B .* φ,
        [0.0, 1, 0], λ)
    dirs = [BMO.direction(BMO.first_ray(b)) for b in BMO.beams(group)]
    @test maximum(d -> norm(d - [0, 1, 0]), dirs) < 1e-9
end

@testset "Decompositions reject fields without forward-travelling light" begin
    f = PlaneField(trace!(detector_at(zdet), agb_at()); size = (16, 16),
        spacing = (w0 / 4, w0 / 4), progress = false)
    backward_only = OpticsBase.PlaneField(f.E, -f.H, f.spacing, f.origin, f.axes, f.λ)
    @test_throws ArgumentError BMO.WavefrontBeamletDecomposition(backward_only)
    @test_throws ArgumentError BMO.GaussianModeDecomposition(backward_only)
    xs = collect(range(-1e-4, 1e-4; length = 8))
    @test_throws ArgumentError BMO.WavefrontBeamletDecomposition(xs, xs, zeros(8, 8),
        zeros(8, 8), [0.0, 1, 0], λ)
    # BMO starts beams in vacuum: a field in a medium is rejected by both decompositions
    in_glass = OpticsBase.PlaneField(f.E, f.H, f.spacing, f.origin, f.axes, f.λ; n = 1.5)
    @test_throws ArgumentError BMO.WavefrontBeamletDecomposition(in_glass)
    @test_throws ArgumentError BMO.GaussianModeDecomposition(in_glass)
end

@testset "Decompositions: masked samples and grazing directions" begin
    # one NaN sample must not hide the rest of the field
    xs = collect(((0:15) .- 8) .* (w0 / 4))
    Eu = [complex(exp(-(x^2 + y^2) / w0^2)) for x in xs, y in xs]
    Eu[3, 3] = NaN
    group = with_logger(NullLogger()) do
        BMO.WavefrontBeamletDecomposition(xs, xs, Eu, zero(Eu), [0.0, 1, 0], λ)
    end
    @test length(BMO.beams(group)) > 0
    @test all(b -> all(isfinite, BMO.polarization(b)), BMO.beams(group))
    amplitude = abs.(Eu)
    phase = zeros(16, 16)
    group = with_logger(NullLogger()) do
        BMO.WavefrontBeamletDecomposition(xs, xs, amplitude, phase, [0.0, 1, 0], λ)
    end
    @test length(BMO.beams(group)) > 0

    # a phase gradient of exactly k is grazing: no direction along `dir` to carry the
    # normal field component, so it is rejected
    e1, e2, n = [1.0, 0, 0], [0.0, 0, -1], [0.0, 1, 0]
    @test BMO._eikonal_direction((a, b) -> 1.0, 2, 2, 3, 3, 1.0, 1.0, 1.0, e1, e2, n) ===
          nothing
    # sin θ = 0.5 along both e1 and e2
    @test BMO._eikonal_direction((a, b) -> 0.5, 2, 2, 3, 3, 1.0, 1.0, 1.0, e1, e2, n) ≈
          0.5 * e1 + 0.5 * e2 + sqrt(0.5) * n
end

end # module
