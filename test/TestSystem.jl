module TestSystem

using BeamletOptics
using Test
# using LinearAlgebra
# using GeometryBasics

const BMO = BeamletOptics

@testset "System" begin
    @testset "Testing implementation" begin
        struct SystemTestBeam{T} <: BMO.AbstractBeam{T, Ray{T}} end
        struct SystemTestObject{T, S} <: BMO.AbstractObject{T} end
        o1 = SystemTestObject{Real, BMO.AbstractShape{Real}}()
        o2 = SystemTestObject{Real, BMO.AbstractShape{Real}}()
        system = System(o1)
        beam = SystemTestBeam{Real}()
        # Test missing implementation warnings
        @test_logs (:warn, "Tracing for $(typeof(beam)) not implemented") BMO.trace_system!(
            system,
            beam)
    end

    # Setup circular multipass cell with flat mirrors
    n_mirrors = 101
    radius = 1
    L = 6 * radius / n_mirrors
    Δθ = 360 / (n_mirrors + 1)
    mirrors = [SquarePlanoMirror2D(L) for _ in 1:n_mirrors]
    θ = 1 * Δθ
    for m in mirrors
        point = radius * [cos(deg2rad(θ)), sin(deg2rad(θ)), 0]
        zrotate3d!(m, deg2rad(θ))
        translate3d!(m, point)
        θ += Δθ
    end
    zrotate3d!.(mirrors, deg2rad(90))

    # Initial ray orientation and position
    dir = [-1, 0, 0]
    Rot = BMO.rotate3d([0, 0, 1], deg2rad(Δθ * 1))
    dir = Vector(Rot * dir)
    origin = [radius, 0, 0] + -1 * dir

    @testset "Testing tracing subroutines" begin
        system = System(mirrors)
        ray = Ray(origin, dir)
        first_obj = mirrors[(n_mirrors + 1) ÷ 2 + 2]
        false_obj = mirrors[(n_mirrors + 1) ÷ 2 + 2 + 1]
        # trace_all
        @test BMO.trace_all(system, ray)[2] === first_obj
        # trace_one
        @test BMO.trace_one(system, ray, BMO.Hint(first_obj))[2] === first_obj
        @test BMO.trace_one(system, ray, BMO.Hint(false_obj))[2] === first_obj
        # tracing step
        obj = BMO.tracing_step!(system, ray, nothing)
        @test obj === first_obj
        @test BMO.intersection(ray) isa Intersection
    end

    @testset "Testing system tracing" begin
        system = System(mirrors)
        first_ray = Ray(origin, dir)
        beam = Beam(first_ray)
        # Test trace_system!
        nmax = 10
        BMO.trace_system!(system, beam, r_max = nmax)
        @test length(BMO.rays(beam)) == nmax
        BMO.trace_system!(system, beam, r_max = 1000000)
        @test length(BMO.rays(beam)) == n_mirrors + 1
        first_ray_dir = BMO.direction(first_ray)
        last_ray_dir = BMO.direction(last(BMO.rays(beam)))
        @test 180 - rad2deg(BMO.angle3d(first_ray_dir, last_ray_dir)) ≈ 2 * Δθ
        @test !isnothing(BMO.intersection(first_ray))
    end

    @testset "Testing StaticSystem tracing" begin
        # same testset as before
        system = StaticSystem(mirrors)
        first_ray = Ray(origin, dir)
        beam = Beam(first_ray)
        # Test trace_system!
        nmax = 10
        BMO.trace_system!(system, beam, r_max = nmax)
        @test length(BMO.rays(beam)) == nmax
        BMO.trace_system!(system, beam, r_max = 1000000)
        @test length(BMO.rays(beam)) == n_mirrors + 1
        first_ray_dir = BMO.direction(first_ray)
        last_ray_dir = BMO.direction(last(BMO.rays(beam)))
        @test 180 - rad2deg(BMO.angle3d(first_ray_dir, last_ray_dir)) ≈ 2 * Δθ
        @test !isnothing(BMO.intersection(first_ray))
    end

    @testset "Dynamic Coincident Boundaries (MultiIntersection)" begin
        lens1 = SphericalLens(Inf, -50e-3, 10e-3, 25.4e-3, 1.5)
        lens2 = SphericalLens(-50e-3, Inf, 10e-3, 25.4e-3, 1.6)
        translate3d!(lens2, [0.0, thickness(lens1), 0.0])
        
        sys = System([lens1, lens2])
        ray = Ray([0.0, -10e-3, 0.0], [0.0, 1.0, 0.0])
        beam = Beam(ray)
        
        solve_system!(sys, beam)
        @test length(rays(beam)) >= 3
    end
end

end # MODULE