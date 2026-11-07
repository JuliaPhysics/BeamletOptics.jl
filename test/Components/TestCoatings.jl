module TestCoatings

using BeamletOptics
using Test
using LinearAlgebra
using AbstractTrees

const BMO = BeamletOptics
const mm = 1e-3
const λ0 = 550e-9
const ns = 1.52

# Quarter-wave layers of a high and a low index
const H = ConstantRefractiveIndex(2.35)
const L = ConstantRefractiveIndex(1.38)

# A 5 mm glass plate at y = 50 mm, and the position of its front face on the y-axis
function glass_plate(n = ns)
    plate = RectangularCompensatorPlate(40mm, 40mm, 5mm, ConstantRefractiveIndex(n))
    translate3d!(plate, [0, 50mm, 0])
    beam = Beam([0, 0, 0], [0, 1, 0], λ0)
    solve_system!(System([plate]), beam)
    return plate, position(rays(beam)[2])
end

# A coating layer flush on the front face at `p`, its normal away from the plate (towards -y)
function front_coating(model, p)
    layer = CoatingLayer(BMO.RectangularFlatMesh(30mm, 30mm), model)
    zrotate3d!(layer, π)
    translate3d!(layer, p)
    return layer
end

@testset "Thin films against analytic values" begin
    # The bare boundary: the Fresnel coefficients
    for cos0 in (1.0, 0.8, 0.3)
        amplitudes, _ = BMO.thin_film_coefficients((), 1.0, ns, cos0, λ0)
        @test all(isapprox.(amplitudes, BMO._fresnel_coefficients(cos0, ns); atol = 1e-15))
    end
    # Quarter-wave anti-reflection layer, Hecht eq. 9.106
    _, (Rs, _, Ts, _) = BMO.thin_film_coefficients((QuarterWaveLayer(L, λ0),), 1.0, ns, 1.0, λ0)
    @test Rs ≈ ((ns - 1.38^2) / (ns + 1.38^2))^2
    @test Rs + Ts ≈ 1
    # Quarter-wave mirror (HL)ᴺH, Macleod eq. 5.2
    for N in (1, 4, 8)
        hr = (QuarterWaveLayer(H, λ0), Iterators.flatten(
            (QuarterWaveLayer(L, λ0), QuarterWaveLayer(H, λ0)) for _ in 1:N)...)
        _, (R, _, T, _) = BMO.thin_film_coefficients(hr, 1.0, ns, 1.0, λ0)
        y = (2.35 / 1.38)^(2N) * 2.35^2 / ns
        @test R ≈ ((1 - y) / (1 + y))^2
        @test R + T ≈ 1
    end
    # Lossless at oblique incidence, s and p
    stack = (QuarterWaveLayer(H, λ0), QuarterWaveLayer(L, λ0))
    _, (Rs, Rp, Ts, Tp) = BMO.thin_film_coefficients(stack, 1.0, ns, cos(deg2rad(50)), 600e-9)
    @test Rs + Ts ≈ 1
    @test Rp + Tp ≈ 1
    # Reciprocity: the same transmittance from both sides
    cos_s = sqrt(1 - (sin(deg2rad(30)) / ns)^2)
    _, (_, _, T1, _) = BMO.thin_film_coefficients(stack, 1.0, ns, cos(deg2rad(30)), λ0)
    _, (_, _, T2, _) = BMO.thin_film_coefficients(reverse(stack), ns, 1.0, cos_s, λ0)
    @test T1 ≈ T2
    # A thick silver layer reflects like bulk silver at normal incidence
    silver = AbsorbingRefractiveIndex(ConstantRefractiveIndex(0.135), 3.985)
    _, (R, _, T, _) = BMO.thin_film_coefficients((ThinFilmLayer(silver, 1e-6),), 1.0, ns, 1.0, 633e-9)
    @test R ≈ ((0.135 - 1)^2 + 3.985^2) / ((0.135 + 1)^2 + 3.985^2)
    @test T < 1e-20
    @test extinction_coefficient(silver, 633e-9) == 3.985
    @test extinction_coefficient(L, λ0) == 0
    @test_throws ArgumentError ThinFilmLayer(L, -1e-9)
end

@testset "Coatings on a glass plate" begin
    plate, p = glass_plate()
    E = [1.0 + 0im, 0, 0]
    # Power behind the front face, relative to the input
    function power_inside(model)
        beam = Beam(PolarizedRay([0, 0, 0], [0, 1, 0], λ0, E))
        solve_system!(System([plate, front_coating(model, p)]), beam)
        r = rays(beam)[2]
        return BMO.refractive_index(r) * norm(BMO.polarization(r))^2
    end
    @test power_inside(AntireflectionCoating()) ≈ 1
    @test power_inside(ThinFilmCoating(QuarterWaveLayer(L, λ0); splitting = BMO.Dominant())) ≈
          1 - ((ns - 1.38^2) / (ns + 1.38^2))^2
    # The bare plate loses the Fresnel reflection
    beam = Beam(PolarizedRay([0, 0, 0], [0, 1, 0], λ0, E))
    solve_system!(System([plate]), beam)
    r = rays(beam)[2]
    @test BMO.refractive_index(r) * norm(BMO.polarization(r))^2 ≈ 1 - ((ns - 1) / (ns + 1))^2

    # A mirror coating reflects back
    beam = Beam([0, 0, 0], [0, 1, 0], λ0)
    solve_system!(System([plate, front_coating(MirrorCoating(), p)]), beam)
    @test BMO.direction(last(rays(beam))) ≈ [0, -1, 0]

    # Fresnel ghosts split at the front face: the reflected child carries R of the power
    beam = Beam(PolarizedRay([0, 0, 0], [0, 1, 0], λ0, E))
    solve_system!(System([plate, front_coating(FresnelGhosts(), p)]), beam; depth_max = 2)
    @test length(children(beam)) == 2
    @test BMO.power_fraction(children(beam)[2], beam) ≈ ((ns - 1) / (ns + 1))^2
end

@testset "Beamsplitter coatings take the power reflectance" begin
    bs = BeamsplitterCoating(0.7)
    @test BMO.reflectance(bs)^2 ≈ 0.7
    @test BMO.transmittance(bs)^2 ≈ 0.3
    @test BMO.reflectance(BeamsplitterCoating()) ≈ √0.5
    @test BMO.transmittance(MirrorCoating()) == 0
    @test_throws ArgumentError BeamsplitterCoating(1.2)
end

@testset "Spectral coatings" begin
    table = BMO.LinearTable([500e-9, 600e-9, 700e-9], [0.0, 0.5, 1.0])
    @test table(550e-9) ≈ 0.25
    @test table(400e-9) == 0
    @test table(800e-9, 0.3) == 1
    @test_throws ArgumentError BMO.LinearTable([2.0, 1.0], [0.0, 1.0])

    # A longpass dichroic at 45°: short wavelengths reflect, long ones pass
    dichroic = CoatingLayer(BMO.RectangularFlatMesh(25mm, 25mm),
        SpectralCoating((λ, θ) -> λ < 800e-9 ? 1.0 : 0.0))
    zrotate3d!(dichroic, deg2rad(45))
    translate3d!(dichroic, [0, 50mm, 0])
    system = System([dichroic])
    for (λ, dir) in ((532e-9, [1, 0, 0]), (900e-9, [0, 1, 0]))
        beam = Beam([0, 0, 0], [0, 1, 0], λ)
        solve_system!(system, beam)
        @test abs.(BMO.direction(last(rays(beam)))) ≈ dir
    end
    # The angle in air of a ray in glass: a filter that reflects above 30° in air reflects a ray at
    # 25° in glass (sin = 0.634 in air) but passes one at 15° (sin = 0.388)
    plate = RectangularCompensatorPlate(30mm, 30mm, 5mm, ConstantRefractiveIndex(1.5))
    translate3d!(plate, [0, 50mm, 0])
    surface_model!(plate, FlatFace([0, 1, 0]), SpectralCoating((λ, θ) -> θ > deg2rad(30) ? 1.0 : 0.0))
    for (θ_glass, reflected) in ((deg2rad(25), true), (deg2rad(15), false))
        θ_air = asin(1.5 * sin(θ_glass))
        beam = Beam([0, 45mm, 0], [sin(θ_air), cos(θ_air), 0], λ0)
        solve_system!(System([plate]), beam)
        @test (BMO.direction(rays(beam)[3])[2] < 0) == reflected
    end
    # A tabulated 50:50 point splits
    split = CoatingLayer(BMO.RectangularFlatMesh(25mm, 25mm),
        SpectralCoating([400e-9, 700e-9], [0.5, 0.5]))
    zrotate3d!(split, deg2rad(45))
    translate3d!(split, [0, 50mm, 0])
    gauss = GaussianBeamlet([0, 0, 0], [0, 1, 0], λ0, 1mm)
    solve_system!(System([split]), gauss)
    @test sum(optical_power, children(gauss)) ≈ optical_power(gauss)
end

@testset "A thin-film beamsplitter in a Mach-Zehnder conserves power" begin
    # A free-standing quarter-wave film of high index splits; its phases follow from the film
    film() = CoatingLayer(BMO.RectangularFlatMesh(BMO.inch, BMO.inch),
        ThinFilmCoating(QuarterWaveLayer(H, λ0)))
    m1, m2 = SquarePlanoMirror2D(BMO.inch), SquarePlanoMirror2D(BMO.inch)
    b1, b2 = film(), film()
    translate3d!(b2, [2BMO.inch, 2BMO.inch, 0])
    translate3d!(m1, [0, 2BMO.inch, 0])
    translate3d!(m2, [2BMO.inch, 0, 0])
    zrotate3d!(b1, deg2rad(225))
    zrotate3d!(b2, deg2rad(45))
    zrotate3d!(m1, deg2rad(225))
    zrotate3d!(m2, deg2rad(45))
    system = System([m1, m2, b1, b2])
    for E in ([0, 0, 1.0 + 0im], [1.0 + 0im, 0, 0])
        beam = Beam(PolarizedRay([0, -0.1, 0], [0, 1, 0], λ0, E))
        solve_system!(system, beam)
        ports = Dict{Vector{Float64}, Vector{ComplexF64}}()
        for leaf in Leaves(beam)
            r = first(rays(leaf))
            k = round.(collect(BMO.direction(r)), digits = 6)
            ports[k] = get(ports, k, zeros(ComplexF64, 3)) + BMO.polarization(r)
        end
        @test length(ports) == 2
        @test sum(norm(E)^2 for E in values(ports)) ≈ 1
    end
end

end # module
