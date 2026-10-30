module TestBeamletFieldInMedia

using BeamletOptics
using LinearAlgebra
using Test

const BMO = BeamletOptics
const mm = 1e-3

#=
Field and power of the two beamlet models in a medium. Power convention: the irradiance is
S = n |E|² / (2 Z₀), hence a beam of fixed power has |E| ∝ 1/√n, whichever beamlet model carries it.
=#

const λ = 1e-6
const w0 = 0.5mm
const P0 = 1e-3
const n_glass = 1.5

source(::Type{GaussianBeamlet}) = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0; P0)
source(::Type{AstigmaticGaussianBeamlet}) =
    AstigmaticGaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0; P0)

chief_beam(g::GaussianBeamlet) = g.chief
chief_beam(agb::AstigmaticGaussianBeamlet) = agb.c

# radius of the beam on the axis at the distance z (circular beams)
radius(g::GaussianBeamlet, z) = first(BMO.gauss_parameters(g, z))
radius(agb::AstigmaticGaussianBeamlet, z) = first(BMO.gauss_parameters(agb, z))

on_axis(::GaussianBeamlet) = 0.0
on_axis(::AstigmaticGaussianBeamlet) = [0.0, 0, 0]

# the peak field [V/m] of a circular beam of the power `P` in a medium of the index `n` and radius `w`
peak_field(P, n, w) = sqrt(4 * BMO.Z_vacuum * P / (n * π * w^2))

# a detector at the end of segment `id` of `b`, hit by that segment only: it measures the beam inside
# the medium of the segment without a surface of the detector in the way
function probe(b, id)
    pd = Detector(40mm)
    translate3d!(pd, BMO.hit_point(BMO.rays(chief_beam(b))[id]))
    push!(pd, b isa GaussianBeamlet ? BMO.GaussianBeamletHit(b, id) : BMO.AstigmaticGaussianBeamletHit(b, id))
    return pd
end

power(pd) = optical_power(pd, RegularGrid(num_points = 301, crop_factor = 4))

function plate_beamsplitter()
    pbs = RectangularPlateBeamsplitter(30mm, 30mm, 5mm, _ -> n_glass)
    translate3d!(pbs, [0, 50mm, 0])
    return pbs
end

function glass_plate()
    plate = Lens(RectangularFlatSurface(30mm), RectangularFlatSurface(30mm), 5mm, _ -> n_glass)
    translate3d!(plate, [0, 50mm, 0])
    return plate
end

@testset "Beamlet fields in media" begin
    models = (GaussianBeamlet, AstigmaticGaussianBeamlet)

    @testset "Plate beamsplitter, transmitted arm in the substrate" begin
        # the splitter divides the power equally, the arm in the substrate is segment 1 of the child
        fields = Dict()
        for M in models
            b = source(M)
            solve_system!(System([plate_beamsplitter()]), b)
            c = b.children[1]
            @test BMO.refractive_index(c, 1) == n_glass
            z = 52.5mm
            E = abs(electric_field(c, on_axis(c), z))
            fields[M] = E
            # |E| ∝ 1/√n at the power P0 / 2
            @test E ≈ peak_field(P0 / 2, n_glass, radius(c, z)) rtol = 1e-3
            # the power of the beamlet is that of the arm, in the substrate as well
            @test optical_power(c) ≈ P0 / 2 rtol = 1e-6
            if c isa AstigmaticGaussianBeamlet
                @test BMO.intensity(c, on_axis(c), z) ≈ n_glass * E^2 / (2 * BMO.Z_vacuum)
                @test BMO.intensity(c, on_axis(c), z) ≈ 2 * (P0 / 2) / (π * radius(c, z)^2) rtol = 1e-3
            end
            # a detector in the substrate reads the same power
            @test power(probe(c, 1)) ≈ P0 / 2 rtol = 1e-3
        end
        @test fields[GaussianBeamlet] ≈ fields[AstigmaticGaussianBeamlet] rtol = 1e-3
    end

    @testset "Plate beamsplitter, tilted" begin
        # the beam is elliptical in the substrate, which the stigmatic model does not describe:
        # the power is compared only
        for M in models
            b = source(M)
            pbs = plate_beamsplitter()
            zrotate3d!(pbs, deg2rad(30))
            solve_system!(System([pbs]), b)
            c = b.children[1]
            @test power(probe(c, 1)) ≈ P0 / 2 rtol = 1e-3
        end
    end

    @testset "Thick glass plate" begin
        # Uncoated surfaces take no Fresnel loss for a GaussianBeamlet, whereas the astigmatic model
        # takes the loss of both surfaces for its polarized chief ray. Per surface the field is
        # that of the Fresnel coefficient, and the power inside the glass is the transmitted one.
        T1 = 4 * n_glass / (n_glass + 1)^2
        loss = Dict(GaussianBeamlet => 1.0, AstigmaticGaussianBeamlet => T1)
        fields = Dict()
        for M in models
            b = source(M)
            pd = Detector(40mm)
            translate3d!(pd, [0, 80mm, 0])
            solve_system!(System([glass_plate(), pd]), b)
            # inside the glass
            z = 52.5mm
            E = abs(electric_field(b, on_axis(b), z))
            fields[M] = E / sqrt(loss[M])
            @test E ≈ peak_field(P0 * loss[M], n_glass, radius(b, z)) rtol = 1e-3
            @test power(probe(b, 2)) ≈ P0 * loss[M] rtol = 1e-3
            # in air behind the plate
            @test optical_power(pd) ≈ P0 * loss[M]^2 rtol = 1e-3
            # the field in air is that of the same beam without a medium
            E_air = norm(electric_field(M == GaussianBeamlet ? BMO.GaussianBeamletHit(b, 3) :
                                        BMO.AstigmaticGaussianBeamletHit(b, 3), [0.0, 80mm, 0]))
            @test E_air ≈ peak_field(P0 * loss[M]^2, 1.0, radius(b, 80mm)) rtol = 1e-3
        end
        @test fields[GaussianBeamlet] ≈ fields[AstigmaticGaussianBeamlet] rtol = 1e-3
    end

    @testset "Detector in glass and in air" begin
        # the same beam of power P0, once in air and once in glass: the field is smaller by √n, the
        # detector reads the power in both
        for M in models
            b = source(M)
            solve_system!(System([glass_plate()]), b)
            hit_air = M == GaussianBeamlet ? BMO.GaussianBeamletHit(b, 1) : BMO.AstigmaticGaussianBeamletHit(b, 1)
            hit_glass = M == GaussianBeamlet ? BMO.GaussianBeamletHit(b, 2) : BMO.AstigmaticGaussianBeamletHit(b, 2)
            @test power(probe(b, 1)) ≈ P0 rtol = 1e-3
            T = M == GaussianBeamlet ? 1.0 : 4 * n_glass / (n_glass + 1)^2
            @test power(probe(b, 2)) ≈ T * P0 rtol = 1e-3
            # the field of the hits at the surface point (the beam radius is the same there)
            p = [0.0, 50mm, 0]
            E_air = norm(electric_field(hit_air, p))
            E_glass = norm(electric_field(hit_glass, p))
            @test E_glass / E_air ≈ sqrt(T / n_glass) rtol = 1e-3
        end
    end
end

end
