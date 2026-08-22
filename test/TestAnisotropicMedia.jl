module TestAnisotropicMedia

using BeamletOptics
using Test
using LinearAlgebra
using StaticArrays
using GeometryBasics

const BMO = BeamletOptics

@testset "Anisotropic and Birefringent Media" begin
    @testset "UniaxialMedium Construction and Properties" begin
        # Calcite: negative uniaxial (no ≈ 1.658, ne ≈ 1.486 at ~589nm)
        calcite = UniaxialMedium(:Calcite, 1.658, 1.486, Point3(0.0, 0.0, 1.0))
        @test is_uniaxial(calcite)
        @test !is_biaxial(calcite)
        @test optic_axis(calcite) ≈ Point3(0.0, 0.0, 1.0)
        @test refractive_index_o(calcite) ≈ 1.658
        @test refractive_index_e(calcite) ≈ 1.486
        @test birefringence(calcite) ≈ (1.486 - 1.658) atol=1e-12
        @test birefringence(calcite) < 0 # Negative uniaxial crystal

        # Quartz: positive uniaxial (no ≈ 1.544, ne ≈ 1.553)
        quartz = UniaxialMedium(:Quartz, 1.544, 1.553, [1.0, 0.0, 0.0])
        @test optic_axis(quartz) ≈ Point3(1.0, 0.0, 0.0)
        @test birefringence(quartz) > 0 # Positive uniaxial crystal
    end

    @testset "Direction-Dependent Refractive Index (Index Ellipsoid)" begin
        no = 1.658
        ne = 1.486
        c_axis = Point3(0.0, 0.0, 1.0)
        crystal = UniaxialMedium(:Crystal, no, ne, c_axis)

        # Propagation parallel to optic axis (θ = 0) -> sees ordinary index no
        k_parallel = Point3(0.0, 0.0, 1.0)
        @test refractive_index(crystal, 1e-6, k_parallel) ≈ no atol=1e-12

        # Propagation perpendicular to optic axis (θ = 90°) -> sees principal extraordinary index ne
        k_perp = Point3(1.0, 0.0, 0.0)
        @test refractive_index(crystal, 1e-6, k_perp) ≈ ne atol=1e-12

        # Propagation at 45 degrees
        k_45 = normalize(Point3(1.0, 0.0, 1.0))
        cosθ = cos(π / 4)
        sinθ = sin(π / 4)
        expected_n45 = (no * ne) / sqrt(ne^2 * cosθ^2 + no^2 * sinθ^2)
        @test refractive_index(crystal, 1e-6, k_45) ≈ expected_n45 atol=1e-12
    end

    @testset "Dielectric Tensor Calculations" begin
        # Isotropic medium
        iso = IsotropicMedium(1.5)
        eps_iso = dielectric_tensor(iso)
        @test eps_iso ≈ SMatrix{3, 3, Float64, 9}(2.25, 0, 0, 0, 2.25, 0, 0, 0, 2.25)

        # Ambient
        eps_amb = dielectric_tensor(Ambient())
        @test eps_amb ≈ SMatrix{3, 3, Float64, 9}(1, 0, 0, 0, 1, 0, 0, 0, 1)

        # Uniaxial crystal along z-axis
        no = 1.5
        ne = 1.7
        crystal_z = UniaxialMedium(:CrystalZ, no, ne, Point3(0.0, 0.0, 1.0))
        eps_z = dielectric_tensor(crystal_z)
        @test eps_z ≈ SMatrix{3, 3, Float64, 9}(no^2, 0, 0, 0, no^2, 0, 0, 0, ne^2) atol=1e-12

        # Biaxial crystal
        biax = BiaxialMedium(:Mica, 1.55, 1.58, 1.60)
        @test is_biaxial(biax)
        @test !is_uniaxial(biax)
        eps_biax = dielectric_tensor(biax)
        @test eps_biax ≈ SMatrix{3, 3, Float64, 9}(1.55^2, 0, 0, 0, 1.58^2, 0, 0, 0, 1.60^2) atol=1e-12
    end

    @testset "Fresnel Power Coefficients & Admittance Factor" begin
        n1 = 1.0
        n2 = 1.5
        n_ratio = n2 / n1
        θi = deg2rad(30)

        # Snell transmission angle
        sin_θt = sin(θi) / n_ratio
        θt = asin(sin_θt)

        adm = admittance_factor(θi, n_ratio)
        adm_expected = admittance_factor(θi, θt, n1, n2)
        @test adm ≈ adm_expected atol=1e-12

        Rs, Rp, Ts, Tp = fresnel_power_coefficients(θi, n_ratio)
        @test Rs + Ts ≈ 1.0 atol=1e-12
        @test Rp + Tp ≈ 1.0 atol=1e-12
        @test Rs >= 0 && Ts >= 0
        @test Rp >= 0 && Tp >= 0
    end
end

end
