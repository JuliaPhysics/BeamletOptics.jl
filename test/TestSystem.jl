module TestSystem

using BeamletOptics
using Test
# using LinearAlgebra
# using GeometryBasics

const BMO = BeamletOptics

@testset "System" begin
    @testset "Testing implementation" begin
        struct SystemTestBeam{T} <: BMO.AbstractBeam{T, Ray{T}} end
        struct SystemTestObject{T, S} <: BMO.AbstractObject{T, S} end
        o1 = SystemTestObject{Real, BMO.AbstractShape{Real}}()
        o2 = SystemTestObject{Real, BMO.AbstractShape{Real}}()
        system = System(o1)
        beam = SystemTestBeam{Real}()
        # Test missing implementation warnings
        @test_logs (:warn, "Tracing for $(typeof(beam)) not implemented") BMO.trace_system!(
            system,
            beam)
        @test_logs (:warn, "Retracing for $(typeof(beam)) not implemented") BMO.retrace_system!(
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
        @test BMO.object(BMO.trace_all(system, ray)) === first_obj
        # trace_one
        @test BMO.object(BMO.trace_one(
            system, ray, BMO.Hint(first_obj))) === first_obj
        @test BMO.object(BMO.trace_one(
            system, ray, BMO.Hint(false_obj))) === first_obj
        # tracing step
        BMO.tracing_step!(system, ray, nothing)
        @test BMO.object(BMO.intersection(ray)) === first_obj
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
        @test BMO.object(BMO.intersection(first_ray)) ===
              mirrors[(n_mirrors + 1) ÷ 2 + 2]
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
        @test BMO.object(BMO.intersection(first_ray)) ===
              mirrors[(n_mirrors + 1) ÷ 2 + 2]
    end

    @testset "Testing system retracing" begin
        system = System(mirrors)
        first_ray = Ray(origin, dir)
        beam = Beam(first_ray)
        t1 = @timed BMO.trace_system!(system, beam, r_max = 1000000)
        t2 = @timed BMO.retrace_system!(system, beam) # for precompilation
        t2 = @timed BMO.retrace_system!(system, beam)
        if t1.time < t2.time
            @warn "Retracing took longer than tracing, something might be bugged...\n   Tracing: $(t1.time) s\n   Retracing: $(t2.time) s"
        end
    end
end

end # MODULE