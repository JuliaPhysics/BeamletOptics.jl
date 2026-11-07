module TestPolarizingComponents

using BeamletOptics
using Test
using LinearAlgebra
using AbstractTrees

const BMO = BeamletOptics
const mm = 1e-3
const λ0 = 633e-9
const n_glass = 1.5
const fresnel_t = 1 - ((n_glass - 1) / (n_glass + 1))^2

flux(r) = BMO.refractive_index(r) * norm(BMO.polarization(r))^2
polarized(E; p = [0, 0, 0.0], d = [0, 1.0, 0]) = Beam(PolarizedRay(p, d, λ0, complex.(E)))

# The power of the leaves of the beam tree by their final direction, rounded
function exits(beam)
    out = Dict{Vector{Float64}, Float64}()
    for leaf in Leaves(beam)
        r = last(rays(leaf))
        k = round.(collect(BMO.direction(r)), digits = 6) .+ 0.0
        out[k] = get(out, k, 0.0) + flux(r)
    end
    return out
end

# Two exit dictionaries agree in their directions and, approximately, in their powers
same_exits(a, b) = keys(a) == keys(b) && all(isapprox(a[k], b[k]; atol = 1e-12) for k in keys(a))

function pbs_at(p; kwargs...)
    pbs = PolarizingCubeBeamsplitter(20mm, ConstantRefractiveIndex(n_glass); kwargs...)
    translate3d!(pbs, p)
    return pbs
end

@testset "Polarizing cube beamsplitter" begin
    system = System([pbs_at([0, 50mm, 0])])
    # p (x) is transmitted, s (z) reflected, each through two uncoated faces of the cube
    for (E, dir) in (([1, 0, 0], [0, 1, 0]), ([0, 0, 1], [-1, 0, 0]))
        beam = polarized(E)
        solve_system!(system, beam)
        @test same_exits(exits(beam), Dict(Float64.(dir) => fresnel_t^2))
    end
    # Light at 45° splits in half
    beam = polarized([1, 0, 1] / √2)
    solve_system!(system, beam)
    @test same_exits(exits(beam), Dict([0.0, 1, 0] => fresnel_t^2 / 2, [-1.0, 0, 0] => fresnel_t^2 / 2))
    # A finite extinction ratio leaks p into the reflected port
    beam = polarized([1, 0, 0])
    solve_system!(System([pbs_at([0, 50mm, 0]; Rp = 0.01)]), beam)
    @test exits(beam)[[-1.0, 0, 0]] ≈ 0.01 * fresnel_t^2
    # Placement does not matter: turned about the beam, the cube reflects what is s for it
    pbs = pbs_at([0, 50mm, 0])
    yrotate3d!(pbs, π / 2)
    beam = polarized([1, 0, 0])
    solve_system!(System([pbs]), beam)
    @test only(keys(exits(beam)))[3] ≈ 1
    @test_throws ArgumentError PolarizingCoating(1.5, 0)
end

@testset "Wave plates" begin
    # A half-wave plate turned by θ turns linear x by 2θ
    for θ in (0, π / 8, π / 6, π / 4)
        hwp = HalfWavePlate(25mm)
        yrotate3d!(hwp, θ)
        translate3d!(hwp, [0, 20mm, 0])
        @test fast_axis(hwp) ≈ [cos(θ), 0, -sin(θ)]
        beam = polarized([1, 0, 0])
        solve_system!(System([hwp]), beam)
        E = BMO.polarization(last(rays(beam)))
        @test abs.(E) ≈ abs.([cos(2θ), 0, sin(2θ)]) atol = 1e-12
    end
    # A quarter-wave plate at 45° makes circular light, two make a half-wave plate
    plates = [QuarterWavePlate(25mm) for _ in 1:2]
    for (k, qwp) in enumerate(plates)
        yrotate3d!(qwp, π / 4)
        translate3d!(qwp, [0, 10mm * k, 0])
    end
    beam = polarized([1, 0, 0])
    solve_system!(System(plates[1:1]), beam)
    E = BMO.polarization(last(rays(beam)))
    @test abs(E[1]) ≈ abs(E[3]) ≈ 1 / √2
    @test abs(angle(E[3] / E[1])) ≈ π / 2
    beam = polarized([1, 0, 0])
    solve_system!(System(plates), beam)
    @test abs.(BMO.polarization(last(rays(beam)))) ≈ [0, 0, 1] atol = 1e-12
    # A zero-order plate is a half-wave plate at half its design wavelength
    qwp = QuarterWavePlate(25mm; λ0 = 2λ0)
    yrotate3d!(qwp, π / 8)
    translate3d!(qwp, [0, 20mm, 0])
    beam = polarized([1, 0, 0])
    solve_system!(System([qwp]), beam)
    @test abs.(BMO.polarization(last(rays(beam)))) ≈ [1, 0, 1] / √2 atol = 1e-12
    @test_throws ArgumentError HalfWavePlate(0)
end

@testset "An optical isolator returns the light to the side port" begin
    # PBS, quarter-wave plate at 45° and a mirror: the light passes the plate twice, a half wave,
    # and returns along z, which the PBS reflects away from the source
    pbs = pbs_at([0, 50mm, 0])
    qwp = QuarterWavePlate(25mm)
    yrotate3d!(qwp, π / 4)
    translate3d!(qwp, [0, 80mm, 0])
    mirror = SquarePlanoMirror2D(30mm)
    zrotate3d!(mirror, π)
    translate3d!(mirror, [0, 100mm, 0])
    beam = polarized([1, 0, 0])
    solve_system!(System([pbs, qwp, mirror]), beam)
    out = exits(beam)
    @test collect(keys(out)) == [[1.0, 0, 0]]
    @test only(values(out)) ≈ fresnel_t^4
end

end # module
