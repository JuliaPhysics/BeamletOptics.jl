# Tests of BeamletOpticsOpticsBaseExt: BMO detector hits -> OpticsBase.PlaneField and
# back via WavefrontBeamletDecomposition. Evaluated with `using OpticsBase` loaded (see
# runtests.jl), which triggers the package extension.
module TestOpticsBaseExt

using BeamletOptics
using OpticsBase
using Test
using LinearAlgebra

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

end # module
